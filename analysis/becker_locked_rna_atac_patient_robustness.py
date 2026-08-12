#!/usr/bin/env python3
"""Patient-aware sensitivity analyses for locked Becker RNA-ATAC coupling."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT
    / "results"
    / "becker_rna_atac_concordance_locked"
    / "becker_rna_atac_paired_scores.tsv"
)
OUT_DIR = ROOT / "results" / "becker_rna_atac_concordance_locked"
SEED = 20260710
N_PERMUTATIONS = 100_000
N_BOOTSTRAPS = 10_000

LOCKED_ROUTE_PAIRS = [
    (
        "locked_route__wnt_tss",
        "rna_epi__ca_route_signature",
        "atac_tss__wnt_stem",
    ),
    (
        "locked_route__wnt_tss_minus_housekeeping",
        "rna_epi__ca_route_signature",
        "atac_contrast__wnt_stem_minus_housekeeping",
    ),
    (
        "locked_route__wnt_tcf_ascl2_axis",
        "rna_epi__ca_route_signature",
        "atac_tf__wnt_tcf_ascl2_axis",
    ),
    (
        "locked_route__wnt_tcf_ascl2_axis_minus_housekeeping",
        "rna_epi__ca_route_signature",
        "atac_tf_contrast__wnt_tcf_ascl2_axis_minus_housekeeping",
    ),
]

COVARIATES = [
    "is_polyp",
    "rna_epi__proliferation_control",
    "log10_tss_total_counts",
]


def bh(values: pd.Series) -> np.ndarray:
    result = np.full(len(values), np.nan)
    valid = values.notna().to_numpy()
    if valid.any():
        result[valid] = multipletests(values.to_numpy()[valid], method="fdr_bh")[1]
    return result


def analysis_frame(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    keep = ["patient_id", "disease_stage_group"] + columns
    out = frame.loc[
        frame["disease_stage_group"].isin(["normal_unaffected", "polyp"]), keep
    ].copy()
    return out.replace([np.inf, -np.inf], np.nan).dropna()


def model_rows(frame: pd.DataFrame, patient_fixed_effects: bool) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for analysis_id, outcome, predictor in LOCKED_ROUTE_PAIRS:
        data = analysis_frame(frame, [outcome, predictor] + COVARIATES)
        design = data[[predictor] + COVARIATES].astype(float).reset_index(drop=True)
        if patient_fixed_effects:
            patient_terms = pd.get_dummies(
                data["patient_id"].reset_index(drop=True),
                prefix="patient",
                drop_first=True,
                dtype=float,
            )
            design = pd.concat([design, patient_terms], axis=1)
        design = sm.add_constant(design, has_constant="add")
        outcome_values = data[outcome].astype(float).reset_index(drop=True)
        groups = data["patient_id"].reset_index(drop=True)
        fit = sm.OLS(outcome_values, design).fit(
            cov_type="cluster",
            cov_kwds={
                "groups": groups,
                "use_correction": True,
                "df_correction": True,
            },
            use_t=True,
        )
        rows.append(
            {
                "analysis_id": analysis_id,
                "model": (
                    "patient_fixed_effects_with_patient_clustered_se"
                    if patient_fixed_effects
                    else "patient_clustered_se"
                ),
                "outcome": outcome,
                "predictor": predictor,
                "covariates": "+".join(COVARIATES),
                "n_samples": len(data),
                "n_patients": data["patient_id"].nunique(),
                "coef": float(fit.params[predictor]),
                "se_patient_cluster": float(fit.bse[predictor]),
                "ci_low": float(fit.conf_int().loc[predictor, 0]),
                "ci_high": float(fit.conf_int().loc[predictor, 1]),
                "p_value": float(fit.pvalues[predictor]),
                "r_squared": float(fit.rsquared),
            }
        )
    output = pd.DataFrame(rows)
    output["q_value"] = bh(output["p_value"])
    return output


def permutation_p_value(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> float:
    x_rank = stats.rankdata(x)
    y_rank = stats.rankdata(y)
    x_rank = (x_rank - x_rank.mean()) / x_rank.std(ddof=0)
    y_rank = (y_rank - y_rank.mean()) / y_rank.std(ddof=0)
    observed = float(np.mean(x_rank * y_rank))
    exceedances = 0
    for _ in range(N_PERMUTATIONS):
        exceedances += abs(float(np.mean(x_rank * rng.permutation(y_rank)))) >= abs(observed)
    return (exceedances + 1) / (N_PERMUTATIONS + 1)


def patient_median_rows(frame: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows: list[dict[str, object]] = []
    for analysis_id, outcome, predictor in LOCKED_ROUTE_PAIRS:
        data = analysis_frame(frame, [outcome, predictor])
        medians = data.groupby("patient_id", as_index=False)[[outcome, predictor]].median()
        rho, asymptotic_p = stats.spearmanr(medians[outcome], medians[predictor])
        rows.append(
            {
                "analysis_id": analysis_id,
                "unit": "patient_median",
                "outcome": outcome,
                "predictor": predictor,
                "n_samples_before_aggregation": len(data),
                "n_patients": len(medians),
                "spearman_rho": float(rho),
                "p_value_asymptotic": float(asymptotic_p),
                "p_value_patient_permutation": permutation_p_value(
                    medians[outcome].to_numpy(),
                    medians[predictor].to_numpy(),
                    rng,
                ),
                "n_permutations": N_PERMUTATIONS,
                "seed": SEED,
            }
        )
    output = pd.DataFrame(rows)
    output["q_value_patient_permutation"] = bh(output["p_value_patient_permutation"])
    return output


def cluster_bootstrap_rows(frame: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows: list[dict[str, object]] = []
    for analysis_id, outcome, predictor in LOCKED_ROUTE_PAIRS:
        data = analysis_frame(frame, [outcome, predictor])
        patients = data["patient_id"].unique()
        blocks = {patient: part for patient, part in data.groupby("patient_id", sort=False)}
        observed = float(stats.spearmanr(data[outcome], data[predictor]).statistic)
        estimates: list[float] = []
        for _ in range(N_BOOTSTRAPS):
            sampled = rng.choice(patients, size=len(patients), replace=True)
            bootstrap = pd.concat([blocks[patient] for patient in sampled], ignore_index=True)
            estimate = stats.spearmanr(bootstrap[outcome], bootstrap[predictor]).statistic
            if np.isfinite(estimate):
                estimates.append(float(estimate))
        estimates_array = np.asarray(estimates)
        rows.append(
            {
                "analysis_id": analysis_id,
                "unit": "patient_cluster_bootstrap",
                "outcome": outcome,
                "predictor": predictor,
                "n_samples": len(data),
                "n_patients": len(patients),
                "observed_sample_spearman_rho": observed,
                "bootstrap_median_rho": float(np.median(estimates_array)),
                "bootstrap_ci_low": float(np.quantile(estimates_array, 0.025)),
                "bootstrap_ci_high": float(np.quantile(estimates_array, 0.975)),
                "n_bootstraps_requested": N_BOOTSTRAPS,
                "n_bootstraps_valid": len(estimates_array),
                "seed": SEED,
            }
        )
    return pd.DataFrame(rows)


def write_summary(
    cluster_models: pd.DataFrame,
    fixed_models: pd.DataFrame,
    medians: pd.DataFrame,
    bootstraps: pd.DataFrame,
) -> None:
    key = "locked_route__wnt_tss"
    cluster = cluster_models.loc[cluster_models["analysis_id"] == key].iloc[0]
    fixed = fixed_models.loc[fixed_models["analysis_id"] == key].iloc[0]
    median = medians.loc[medians["analysis_id"] == key].iloc[0]
    bootstrap = bootstraps.loc[bootstraps["analysis_id"] == key].iloc[0]
    text = (
        "Patient-aware Becker locked RNA-ATAC robustness\n"
        "================================================\n"
        f"Primary analysis: {key}\n"
        f"Samples/patients: {int(cluster['n_samples'])}/{int(cluster['n_patients'])}\n"
        f"Patient-clustered adjusted coefficient: {cluster['coef']:.6f} "
        f"(95% CI {cluster['ci_low']:.6f} to {cluster['ci_high']:.6f}; "
        f"P={cluster['p_value']:.6g}; BH q={cluster['q_value']:.6g})\n"
        f"Patient-fixed-effect coefficient with clustered SE: {fixed['coef']:.6f} "
        f"(95% CI {fixed['ci_low']:.6f} to {fixed['ci_high']:.6f}; "
        f"P={fixed['p_value']:.6g}; BH q={fixed['q_value']:.6g})\n"
        f"Patient-median Spearman rho: {median['spearman_rho']:.6f} "
        f"(patient-label permutation P={median['p_value_patient_permutation']:.6g}; "
        f"BH q={median['q_value_patient_permutation']:.6g})\n"
        f"Patient-cluster bootstrap sample-level rho: {bootstrap['observed_sample_spearman_rho']:.6f} "
        f"(95% percentile CI {bootstrap['bootstrap_ci_low']:.6f} to "
        f"{bootstrap['bootstrap_ci_high']:.6f}; {int(bootstrap['n_bootstraps_valid'])} valid replicates)\n"
    )
    (OUT_DIR / "becker_locked_rna_atac_patient_robustness_summary.txt").write_text(
        text,
        encoding="utf-8",
    )


def main() -> None:
    frame = pd.read_csv(INPUT, sep="\t")
    cluster_models = model_rows(frame, patient_fixed_effects=False)
    fixed_models = model_rows(frame, patient_fixed_effects=True)
    medians = patient_median_rows(frame)
    bootstraps = cluster_bootstrap_rows(frame)

    cluster_models.to_csv(
        OUT_DIR / "becker_locked_rna_atac_patient_cluster_models.tsv",
        sep="\t",
        index=False,
    )
    fixed_models.to_csv(
        OUT_DIR / "becker_locked_rna_atac_patient_fixed_effect_models.tsv",
        sep="\t",
        index=False,
    )
    medians.to_csv(
        OUT_DIR / "becker_locked_rna_atac_patient_median_correlations.tsv",
        sep="\t",
        index=False,
    )
    bootstraps.to_csv(
        OUT_DIR / "becker_locked_rna_atac_patient_cluster_bootstrap.tsv",
        sep="\t",
        index=False,
    )
    write_summary(cluster_models, fixed_models, medians, bootstraps)


if __name__ == "__main__":
    main()
