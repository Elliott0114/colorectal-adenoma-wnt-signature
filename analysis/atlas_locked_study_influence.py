#!/usr/bin/env python3
"""Study-support and leave-one-study-out checks for locked CRC Atlas models."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "route_signature_locked" / "atlas_locked_signature_donor_scores.tsv"
OUT_DIR = ROOT / "results" / "route_signature_locked"

REFERENCE = "normal_epithelial"
STATES = [
    REFERENCE,
    "polyp_epithelial",
    "polyp_cancer",
    "tumor_epithelial",
    "tumor_cancer",
    "metastasis_epithelial",
    "metastasis_cancer",
]
OUTCOMES = ["score__ca_route_signature", "score__wnt_stem"]


def prepare() -> pd.DataFrame:
    frame = pd.read_csv(INPUT, sep="\t")
    frame = frame.loc[frame["n_cells_sampled"] >= 20].copy()
    frame["log10_n_cells_sampled"] = np.log10(frame["n_cells_sampled"].clip(lower=1))
    return frame


def design_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    carrier = pd.Categorical(frame["carrier_group"], categories=STATES)
    carrier_design = pd.get_dummies(carrier, prefix="carrier_group", drop_first=True, dtype=float)
    observed_studies = sorted(frame["study_id"].unique())
    study = pd.Categorical(frame["study_id"], categories=observed_studies)
    study_design = pd.get_dummies(study, prefix="study_id", drop_first=True, dtype=float)
    design = pd.concat(
        [carrier_design.reset_index(drop=True), study_design.reset_index(drop=True)],
        axis=1,
    )
    design["score__proliferation_control"] = frame["score__proliferation_control"].to_numpy(float)
    design["log10_n_cells_sampled"] = frame["log10_n_cells_sampled"].to_numpy(float)
    varying = design.nunique(dropna=False) > 1
    return sm.add_constant(design.loc[:, varying], has_constant="add")


def fit_global(frame: pd.DataFrame, outcome: str, omitted_study: str) -> list[dict[str, object]]:
    columns = [
        outcome,
        "donor_id",
        "carrier_group",
        "study_id",
        "score__proliferation_control",
        "log10_n_cells_sampled",
    ]
    data = frame[columns].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if omitted_study != "__NONE__":
        data = data.loc[data["study_id"] != omitted_study].copy()
    x = design_matrix(data)
    fit = sm.OLS(data[outcome].to_numpy(float), x).fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data["donor_id"].to_numpy(),
            "use_correction": True,
            "df_correction": True,
        },
        use_t=True,
    )
    rows: list[dict[str, object]] = []
    for state in STATES[1:]:
        term = f"carrier_group_{state}"
        target = data.loc[data["carrier_group"] == state]
        estimable = term in x.columns and len(target) > 0
        rows.append(
            {
                "outcome": outcome,
                "omitted_study": omitted_study,
                "state": state,
                "term": term,
                "estimable": estimable,
                "n_observations": len(data),
                "n_donors": data["donor_id"].nunique(),
                "n_target_observations": len(target),
                "n_target_donors": target["donor_id"].nunique(),
                "n_reference_observations": int((data["carrier_group"] == REFERENCE).sum()),
                "n_reference_donors": data.loc[
                    data["carrier_group"] == REFERENCE, "donor_id"
                ].nunique(),
                "coef": float(fit.params[term]) if estimable else np.nan,
                "se_cluster": float(fit.bse[term]) if estimable else np.nan,
                "ci_low": float(fit.conf_int().loc[term, 0]) if estimable else np.nan,
                "ci_high": float(fit.conf_int().loc[term, 1]) if estimable else np.nan,
                "p_value": float(fit.pvalues[term]) if estimable else np.nan,
            }
        )
    return rows


def leave_one_study_out(frame: pd.DataFrame) -> pd.DataFrame:
    studies = sorted(frame["study_id"].unique())
    rows: list[dict[str, object]] = []
    for outcome in OUTCOMES:
        rows.extend(fit_global(frame, outcome, "__NONE__"))
        for study in studies:
            rows.extend(fit_global(frame, outcome, study))
    return pd.DataFrame(rows)


def influence_summary(results: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (outcome, state), part in results.groupby(["outcome", "state"], sort=False):
        full = part.loc[part["omitted_study"] == "__NONE__"].iloc[0]
        loo = part.loc[(part["omitted_study"] != "__NONE__") & part["estimable"]].copy()
        worst_low = loo.loc[loo["coef"].idxmin()]
        worst_high = loo.loc[loo["coef"].idxmax()]
        chen = loo.loc[loo["omitted_study"] == "Chen_2021_Cell"]
        rows.append(
            {
                "outcome": outcome,
                "state": state,
                "full_coef": full["coef"],
                "full_ci_low": full["ci_low"],
                "full_ci_high": full["ci_high"],
                "full_p_value": full["p_value"],
                "n_eligible_omissions": len(loo),
                "loo_min_coef": loo["coef"].min(),
                "loo_max_coef": loo["coef"].max(),
                "loo_positive_fraction": (loo["coef"] > 0).mean(),
                "loo_ci_excludes_zero_positive_fraction": (loo["ci_low"] > 0).mean(),
                "minimum_target_donors_after_omission": loo["n_target_donors"].min(),
                "lowest_coef_omitted_study": worst_low["omitted_study"],
                "highest_coef_omitted_study": worst_high["omitted_study"],
                "omit_chen_2021_coef": chen.iloc[0]["coef"] if len(chen) else np.nan,
                "omit_chen_2021_ci_low": chen.iloc[0]["ci_low"] if len(chen) else np.nan,
                "omit_chen_2021_ci_high": chen.iloc[0]["ci_high"] if len(chen) else np.nan,
                "omit_chen_2021_p_value": chen.iloc[0]["p_value"] if len(chen) else np.nan,
                "omit_chen_2021_target_donors": (
                    chen.iloc[0]["n_target_donors"] if len(chen) else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def study_support(frame: pd.DataFrame) -> pd.DataFrame:
    observations = (
        frame.groupby(["study_id", "carrier_group"], observed=True)
        .size()
        .rename("n_observations")
    )
    donors = (
        frame.groupby(["study_id", "carrier_group"], observed=True)["donor_id"]
        .nunique()
        .rename("n_donors")
    )
    cells = (
        frame.groupby(["study_id", "carrier_group"], observed=True)["n_cells_sampled"]
        .sum()
        .rename("n_cells_sampled_total")
    )
    return pd.concat([observations, donors, cells], axis=1).reset_index()


def within_study_contrasts(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for study, study_frame in frame.groupby("study_id", sort=True):
        for state in STATES[1:]:
            data = study_frame.loc[
                study_frame["carrier_group"].isin([REFERENCE, state])
            ].copy()
            target_donors = data.loc[data["carrier_group"] == state, "donor_id"].nunique()
            reference_donors = data.loc[
                data["carrier_group"] == REFERENCE, "donor_id"
            ].nunique()
            if target_donors < 3 or reference_donors < 3 or data["donor_id"].nunique() < 8:
                continue
            data["is_target"] = (data["carrier_group"] == state).astype(float)
            for outcome in OUTCOMES:
                columns = [
                    outcome,
                    "is_target",
                    "score__proliferation_control",
                    "log10_n_cells_sampled",
                    "donor_id",
                ]
                model_data = data[columns].replace([np.inf, -np.inf], np.nan).dropna()
                x = sm.add_constant(
                    model_data[
                        ["is_target", "score__proliferation_control", "log10_n_cells_sampled"]
                    ].astype(float),
                    has_constant="add",
                )
                fit = sm.OLS(model_data[outcome].astype(float), x).fit(
                    cov_type="cluster",
                    cov_kwds={
                        "groups": model_data["donor_id"],
                        "use_correction": True,
                        "df_correction": True,
                    },
                    use_t=True,
                )
                rows.append(
                    {
                        "study_id": study,
                        "outcome": outcome,
                        "state": state,
                        "n_observations": len(model_data),
                        "n_donors": model_data["donor_id"].nunique(),
                        "n_target_donors": target_donors,
                        "n_reference_donors": reference_donors,
                        "coef": float(fit.params["is_target"]),
                        "se_cluster": float(fit.bse["is_target"]),
                        "ci_low": float(fit.conf_int().loc["is_target", 0]),
                        "ci_high": float(fit.conf_int().loc["is_target", 1]),
                        "p_value": float(fit.pvalues["is_target"]),
                    }
                )
    return pd.DataFrame(rows)


def write_summary(influence: pd.DataFrame, support: pd.DataFrame) -> None:
    route = influence.loc[influence["outcome"] == "score__ca_route_signature"]
    lines = [
        "CRC Atlas locked-score study influence summary",
        "==============================================",
        f"Studies represented: {support['study_id'].nunique()}",
        "",
    ]
    for _, row in route.iterrows():
        lines.append(
            f"{row['state']}: full beta={row['full_coef']:.6f}; "
            f"LOO range={row['loo_min_coef']:.6f} to {row['loo_max_coef']:.6f}; "
            f"positive in {row['loo_positive_fraction']:.1%} of omissions; "
            f"minimum target donors={int(row['minimum_target_donors_after_omission'])}; "
            f"omit Chen 2021 beta={row['omit_chen_2021_coef']:.6f} "
            f"(95% CI {row['omit_chen_2021_ci_low']:.6f} to "
            f"{row['omit_chen_2021_ci_high']:.6f}; "
            f"target donors={int(row['omit_chen_2021_target_donors'])})."
        )
    (OUT_DIR / "atlas_locked_study_influence_summary.txt").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    frame = prepare()
    support = study_support(frame)
    loo = leave_one_study_out(frame)
    influence = influence_summary(loo)
    within = within_study_contrasts(frame)

    support.to_csv(OUT_DIR / "atlas_locked_state_study_support.tsv", sep="\t", index=False)
    loo.to_csv(OUT_DIR / "atlas_locked_leave_one_study_out.tsv", sep="\t", index=False)
    influence.to_csv(OUT_DIR / "atlas_locked_study_influence.tsv", sep="\t", index=False)
    within.to_csv(OUT_DIR / "atlas_locked_within_study_contrasts.tsv", sep="\t", index=False)
    write_summary(influence, support)


if __name__ == "__main__":
    main()
