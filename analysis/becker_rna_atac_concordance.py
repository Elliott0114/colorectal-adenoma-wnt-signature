#!/usr/bin/env python3
"""Becker sample-level RNA-ATAC concordance analysis.

This analysis connects the Becker snRNA route scores to matched multiome
sample-level promoter-accessibility scores. It is intentionally sample-level:
the goal is to test whether the epithelial WNT/adenoma-memory RNA signal is
coupled to the ATAC regulatory layer without claiming motif activity or direct
TF binding.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "becker_rna_atac_concordance"

RNA_ROUTE_PATH = ROOT / "results" / "becker_route" / "becker_sample_module_scores.tsv"
RNA_SIGNATURE_PATH = ROOT / "results" / "route_signature" / "becker_conventional_route_signature_scores.tsv"
TSS_PATH = ROOT / "results" / "becker_multiome_tss" / "becker_multiome_tss_module_scores.tsv"
TF_AXIS_PATH = ROOT / "results" / "becker_multiome_tf_axis" / "becker_multiome_tf_axis_scores.tsv"


PRIMARY_CORRELATIONS = [
    (
        "rna_epi_wnt_stem__atac_wnt_stem",
        "rna_epi__wnt_stem",
        "atac_tss__wnt_stem",
        "Curated epithelial RNA WNT/stem versus WNT/stem TSS accessibility",
    ),
    (
        "rna_epi_wnt_core_ihc__atac_wnt_core_ihc",
        "rna_epi__wnt_core_ihc",
        "atac_tss__wnt_core_ihc",
        "FFPE-bridge RNA WNT core versus WNT core TSS accessibility",
    ),
    (
        "rna_epi_ca_route__atac_wnt_stem",
        "rna_epi__ca_route_signature",
        "atac_tss__wnt_stem",
        "Data-driven conventional adenoma RNA route versus WNT/stem TSS accessibility",
    ),
    (
        "rna_epi_ca_route__atac_wnt_minus_housekeeping",
        "rna_epi__ca_route_signature",
        "atac_contrast__wnt_stem_minus_housekeeping",
        "Data-driven conventional adenoma RNA route versus housekeeping-adjusted WNT/stem accessibility",
    ),
    (
        "rna_epi_wnt_stem__atac_wnt_tcf_ascl2_axis",
        "rna_epi__wnt_stem",
        "atac_tf__wnt_tcf_ascl2_axis",
        "Curated epithelial RNA WNT/stem versus WNT/TCF/ASCL2 promoter-axis accessibility",
    ),
    (
        "rna_epi_ca_route__atac_wnt_tcf_ascl2_axis",
        "rna_epi__ca_route_signature",
        "atac_tf__wnt_tcf_ascl2_axis",
        "Data-driven conventional adenoma RNA route versus WNT/TCF/ASCL2 promoter-axis accessibility",
    ),
    (
        "rna_epi_ca_route__atac_wnt_tcf_ascl2_minus_housekeeping",
        "rna_epi__ca_route_signature",
        "atac_tf_contrast__wnt_tcf_ascl2_axis_minus_housekeeping",
        "Data-driven conventional adenoma RNA route versus housekeeping-adjusted WNT/TCF/ASCL2 promoter-axis accessibility",
    ),
    (
        "rna_epi_antigen_ifn__atac_antigen_ifn",
        "rna_epi__antigen_presentation_ifn",
        "atac_tss__antigen_presentation_ifn",
        "RNA antigen-presentation/IFN versus antigen-presentation/IFN TSS accessibility",
    ),
    (
        "rna_epi_proliferation__atac_proliferation",
        "rna_epi__proliferation_control",
        "atac_tss__proliferation_control",
        "RNA proliferation-control versus proliferation-control TSS accessibility",
    ),
]


PRIMARY_MODELS = [
    (
        "rna_epi_wnt_stem_from_atac_wnt_stem",
        "rna_epi__wnt_stem",
        "atac_tss__wnt_stem",
        ["is_polyp", "rna_epi__proliferation_control", "log10_tss_total_counts"],
    ),
    (
        "rna_epi_wnt_stem_from_atac_wnt_minus_housekeeping",
        "rna_epi__wnt_stem",
        "atac_contrast__wnt_stem_minus_housekeeping",
        ["is_polyp", "rna_epi__proliferation_control", "log10_tss_total_counts"],
    ),
    (
        "rna_epi_ca_route_from_atac_wnt_stem",
        "rna_epi__ca_route_signature",
        "atac_tss__wnt_stem",
        ["is_polyp", "rna_epi__proliferation_control", "log10_tss_total_counts"],
    ),
    (
        "rna_epi_ca_route_from_atac_wnt_tcf_ascl2_axis",
        "rna_epi__ca_route_signature",
        "atac_tf__wnt_tcf_ascl2_axis",
        ["is_polyp", "rna_epi__proliferation_control", "log10_tss_total_counts"],
    ),
    (
        "rna_epi_ca_route_from_atac_wnt_tcf_ascl2_minus_housekeeping",
        "rna_epi__ca_route_signature",
        "atac_tf_contrast__wnt_tcf_ascl2_axis_minus_housekeeping",
        ["is_polyp", "rna_epi__proliferation_control", "log10_tss_total_counts"],
    ),
]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def bh_adjust(frame: pd.DataFrame, p_col: str = "p_value", q_col: str = "q_value") -> pd.DataFrame:
    out = frame.copy()
    if out.empty or p_col not in out.columns:
        out[q_col] = np.nan
        return out
    mask = out[p_col].notna()
    out[q_col] = np.nan
    if mask.any():
        out.loc[mask, q_col] = multipletests(out.loc[mask, p_col], method="fdr_bh")[1]
    return out


def zscore(series: pd.Series) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce")
    sd = vals.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return vals * np.nan
    return (vals - vals.mean()) / sd


def standardize_inputs() -> pd.DataFrame:
    route = read_tsv(RNA_ROUTE_PATH)
    signature = read_tsv(RNA_SIGNATURE_PATH)
    tss = read_tsv(TSS_PATH)
    tf = read_tsv(TF_AXIS_PATH)

    route_keep = route[
        [
            "geo_accession",
            "sample_title",
            "sample_name_from_title",
            "patient_id",
            "disease_stage_group",
            "familial_adenomatous_polyposis",
            "sex",
            "score_epi__wnt_stem",
            "score_epi__wnt_core_ihc",
            "score_epi__antigen_presentation_ifn",
            "score_epi__proliferation_control",
            "score_all__wnt_stem",
            "score_all__wnt_core_ihc",
            "score_all__antigen_presentation_ifn",
            "score_all__proliferation_control",
            "log10_n_nuclei",
        ]
    ].rename(
        columns={
            "geo_accession": "scrna_geo_accession",
            "sample_name_from_title": "sample_core",
            "score_epi__wnt_stem": "rna_epi__wnt_stem",
            "score_epi__wnt_core_ihc": "rna_epi__wnt_core_ihc",
            "score_epi__antigen_presentation_ifn": "rna_epi__antigen_presentation_ifn",
            "score_epi__proliferation_control": "rna_epi__proliferation_control",
            "score_all__wnt_stem": "rna_all__wnt_stem",
            "score_all__wnt_core_ihc": "rna_all__wnt_core_ihc",
            "score_all__antigen_presentation_ifn": "rna_all__antigen_presentation_ifn",
            "score_all__proliferation_control": "rna_all__proliferation_control",
        }
    )

    sig_keep = signature[
        [
            "geo_accession",
            "score_epi__ca_route_signature",
            "score_all__ca_route_signature",
            "score_epi__ca_route_n_up_present",
            "score_epi__ca_route_n_down_present",
        ]
    ].rename(
        columns={
            "geo_accession": "scrna_geo_accession",
            "score_epi__ca_route_signature": "rna_epi__ca_route_signature",
            "score_all__ca_route_signature": "rna_all__ca_route_signature",
            "score_epi__ca_route_n_up_present": "rna_epi__ca_route_n_up_present",
            "score_epi__ca_route_n_down_present": "rna_epi__ca_route_n_down_present",
        }
    )

    tss_keep = tss[
        [
            "multiome_geo_accession",
            "scrna_geo_accession",
            "sample_core",
            "score__wnt_stem",
            "score__wnt_core_ihc",
            "score__antigen_presentation_ifn",
            "score__proliferation_control",
            "score__housekeeping_control",
            "contrast__wnt_stem_minus_housekeeping",
            "contrast__wnt_core_ihc_minus_housekeeping",
            "log10_tss_total_counts",
        ]
    ].rename(
        columns={
            "score__wnt_stem": "atac_tss__wnt_stem",
            "score__wnt_core_ihc": "atac_tss__wnt_core_ihc",
            "score__antigen_presentation_ifn": "atac_tss__antigen_presentation_ifn",
            "score__proliferation_control": "atac_tss__proliferation_control",
            "score__housekeeping_control": "atac_tss__housekeeping_control",
            "contrast__wnt_stem_minus_housekeeping": "atac_contrast__wnt_stem_minus_housekeeping",
            "contrast__wnt_core_ihc_minus_housekeeping": "atac_contrast__wnt_core_ihc_minus_housekeeping",
        }
    )

    tf_keep = tf[
        [
            "multiome_geo_accession",
            "sample_core",
            "score__wnt_tcf_ascl2_axis",
            "score__intestinal_hnf4_cdx_axis",
            "score__ifn_irf_stat_axis",
            "score__serrated_metaplasia_tf_axis",
            "score__housekeeping_control",
            "contrast__wnt_tcf_ascl2_axis_minus_housekeeping",
            "contrast__ifn_irf_stat_axis_minus_housekeeping",
            "contrast__serrated_metaplasia_tf_axis_minus_housekeeping",
        ]
    ].rename(
        columns={
            "score__wnt_tcf_ascl2_axis": "atac_tf__wnt_tcf_ascl2_axis",
            "score__intestinal_hnf4_cdx_axis": "atac_tf__intestinal_hnf4_cdx_axis",
            "score__ifn_irf_stat_axis": "atac_tf__ifn_irf_stat_axis",
            "score__serrated_metaplasia_tf_axis": "atac_tf__serrated_metaplasia_tf_axis",
            "score__housekeeping_control": "atac_tf__housekeeping_control",
            "contrast__wnt_tcf_ascl2_axis_minus_housekeeping": "atac_tf_contrast__wnt_tcf_ascl2_axis_minus_housekeeping",
            "contrast__ifn_irf_stat_axis_minus_housekeeping": "atac_tf_contrast__ifn_irf_stat_axis_minus_housekeeping",
            "contrast__serrated_metaplasia_tf_axis_minus_housekeeping": "atac_tf_contrast__serrated_metaplasia_tf_axis_minus_housekeeping",
        }
    )

    frame = (
        tss_keep.merge(tf_keep, on=["multiome_geo_accession", "sample_core"], how="left")
        .merge(route_keep, on=["scrna_geo_accession", "sample_core"], how="left", suffixes=("", "_rna"))
        .merge(sig_keep, on="scrna_geo_accession", how="left")
    )
    frame["is_polyp"] = (frame["disease_stage_group"] == "polyp").astype(int)
    frame["is_crc"] = (frame["disease_stage_group"] == "crc").astype(int)
    frame["is_fap"] = (frame["familial_adenomatous_polyposis"] == "Y").astype(int)
    frame["analysis_set"] = np.where(frame["disease_stage_group"].isin(["normal_unaffected", "polyp"]), "normal_polyp", "all")

    score_cols = [c for c in frame.columns if c.startswith(("rna_", "atac_")) and frame[c].dtype.kind in "ifc"]
    for col in score_cols:
        frame[f"paired_z__{col}"] = zscore(frame[col])
    return frame


def subset_frame(frame: pd.DataFrame, subset: str) -> pd.DataFrame:
    if subset == "all_samples":
        return frame.copy()
    if subset == "normal_polyp":
        return frame[frame["disease_stage_group"].isin(["normal_unaffected", "polyp"])].copy()
    if subset == "fap_normal_polyp":
        return frame[(frame["disease_stage_group"].isin(["normal_unaffected", "polyp"])) & (frame["is_fap"] == 1)].copy()
    raise ValueError(f"Unknown subset: {subset}")


def correlation_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for subset in ["all_samples", "normal_polyp", "fap_normal_polyp"]:
        part = subset_frame(frame, subset)
        for analysis_id, x_col, y_col, description in PRIMARY_CORRELATIONS:
            valid = part[[x_col, y_col]].dropna()
            if len(valid) < 5 or valid[x_col].nunique() < 3 or valid[y_col].nunique() < 3:
                rho = pearson = p_spearman = p_pearson = np.nan
            else:
                rho, p_spearman = stats.spearmanr(valid[x_col], valid[y_col])
                pearson, p_pearson = stats.pearsonr(valid[x_col], valid[y_col])
            rows.append(
                {
                    "subset": subset,
                    "analysis_id": analysis_id,
                    "description": description,
                    "rna_feature": x_col,
                    "atac_feature": y_col,
                    "n": len(valid),
                    "spearman_rho": rho,
                    "p_spearman": p_spearman,
                    "pearson_r": pearson,
                    "p_pearson": p_pearson,
                }
            )
    out = pd.DataFrame(rows)
    out = bh_adjust(out, "p_spearman", "q_spearman")
    out = bh_adjust(out, "p_pearson", "q_pearson")
    return out


def ols_model(data: pd.DataFrame, outcome: str, predictor: str, covariates: list[str]) -> dict[str, object]:
    cols = [outcome, predictor] + covariates
    valid = data[cols].replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid) <= len(covariates) + 4:
        return {
            "n": len(valid),
            "predictor_coef": np.nan,
            "predictor_se_hc3": np.nan,
            "predictor_p_value": np.nan,
            "model_r2": np.nan,
        }
    y = valid[outcome].astype(float)
    x = valid[[predictor] + covariates].astype(float)
    x = sm.add_constant(x, has_constant="add")
    model = sm.OLS(y, x).fit(cov_type="HC3")
    return {
        "n": len(valid),
        "predictor_coef": model.params.get(predictor, np.nan),
        "predictor_se_hc3": model.bse.get(predictor, np.nan),
        "predictor_p_value": model.pvalues.get(predictor, np.nan),
        "model_r2": model.rsquared,
    }


def adjusted_models(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for subset in ["normal_polyp", "fap_normal_polyp"]:
        part = subset_frame(frame, subset)
        for model_id, outcome, predictor, covariates in PRIMARY_MODELS:
            result = ols_model(part, outcome, predictor, covariates)
            rows.append(
                {
                    "subset": subset,
                    "model_id": model_id,
                    "outcome": outcome,
                    "predictor": predictor,
                    "covariates": "+".join(covariates),
                    **result,
                }
            )
    out = pd.DataFrame(rows)
    out = bh_adjust(out, "predictor_p_value", "predictor_q_value")
    return out


def residualize(values: pd.Series, covars: pd.DataFrame) -> pd.Series:
    valid = pd.concat([values, covars], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid) <= covars.shape[1] + 3:
        return pd.Series(index=values.index, dtype=float)
    y = valid.iloc[:, 0].astype(float)
    x = sm.add_constant(valid.iloc[:, 1:].astype(float), has_constant="add")
    fit = sm.OLS(y, x).fit()
    out = pd.Series(index=values.index, dtype=float)
    out.loc[valid.index] = fit.resid
    return out


def residual_correlations(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for subset in ["normal_polyp", "fap_normal_polyp"]:
        part = subset_frame(frame, subset)
        covariates = part[["is_polyp", "rna_epi__proliferation_control", "log10_tss_total_counts"]]
        for analysis_id, x_col, y_col, description in PRIMARY_CORRELATIONS:
            x_resid = residualize(part[x_col], covariates)
            y_resid = residualize(part[y_col], covariates)
            valid = pd.concat([x_resid.rename("x"), y_resid.rename("y")], axis=1).dropna()
            if len(valid) < 5 or valid["x"].nunique() < 3 or valid["y"].nunique() < 3:
                rho = p_value = np.nan
            else:
                rho, p_value = stats.spearmanr(valid["x"], valid["y"])
            rows.append(
                {
                    "subset": subset,
                    "analysis_id": analysis_id,
                    "description": description,
                    "rna_feature": x_col,
                    "atac_feature": y_col,
                    "n": len(valid),
                    "residual_spearman_rho": rho,
                    "p_value": p_value,
                    "covariates": "is_polyp+rna_epi__proliferation_control+log10_tss_total_counts",
                }
            )
    out = pd.DataFrame(rows)
    out = bh_adjust(out, "p_value", "q_value")
    return out


def write_summary(
    path: Path,
    frame: pd.DataFrame,
    corrs: pd.DataFrame,
    models: pd.DataFrame,
    residuals: pd.DataFrame,
) -> None:
    def focus_corr(subset: str, analysis_id: str) -> pd.Series:
        match = corrs[(corrs["subset"] == subset) & (corrs["analysis_id"] == analysis_id)]
        if match.empty:
            return pd.Series(dtype=object)
        return match.iloc[0]

    def focus_model(subset: str, model_id: str) -> pd.Series:
        match = models[(models["subset"] == subset) & (models["model_id"] == model_id)]
        if match.empty:
            return pd.Series(dtype=object)
        return match.iloc[0]

    lines = [
        "Becker RNA-ATAC sample-level concordance analysis",
        "=================================================",
        "",
        "Purpose: connect Becker snRNA route-expression scores to matched multiome promoter-accessibility scores at sample level.",
        "Scope: this is RNA/accessibility concordance and does not claim motif deviation, enhancer activity, or direct TF binding.",
        "",
        f"Matched multiome/snRNA samples: {len(frame)}",
        f"Normal/unaffected + polyp main-analysis samples: {frame['disease_stage_group'].isin(['normal_unaffected', 'polyp']).sum()}",
        f"FAP normal/unaffected + polyp sensitivity samples: {((frame['disease_stage_group'].isin(['normal_unaffected', 'polyp'])) & (frame['is_fap'] == 1)).sum()}",
        f"Disease-stage groups: {frame['disease_stage_group'].value_counts().to_dict()}",
        "",
        "Primary RNA-ATAC correlations:",
    ]

    for analysis_id in [
        "rna_epi_wnt_stem__atac_wnt_stem",
        "rna_epi_ca_route__atac_wnt_stem",
        "rna_epi_ca_route__atac_wnt_tcf_ascl2_axis",
        "rna_epi_ca_route__atac_wnt_tcf_ascl2_minus_housekeeping",
    ]:
        row = focus_corr("normal_polyp", analysis_id)
        if not row.empty:
            lines.append(
                f"- normal/polyp {analysis_id}: n = {int(row['n'])}, Spearman rho = {row['spearman_rho']:.3f}, p = {row['p_spearman']:.3g}, BH q = {row['q_spearman']:.3g}."
            )

    lines.extend(["", "Primary adjusted concordance models:"])
    for model_id in [
        "rna_epi_wnt_stem_from_atac_wnt_stem",
        "rna_epi_ca_route_from_atac_wnt_stem",
        "rna_epi_ca_route_from_atac_wnt_tcf_ascl2_axis",
        "rna_epi_ca_route_from_atac_wnt_tcf_ascl2_minus_housekeeping",
    ]:
        row = focus_model("normal_polyp", model_id)
        if not row.empty:
            lines.append(
                f"- normal/polyp {model_id}: n = {int(row['n'])}, predictor coef = {row['predictor_coef']:.3f}, HC3 p = {row['predictor_p_value']:.3g}, BH q = {row['predictor_q_value']:.3g}."
            )

    best_resid = residuals[
        (residuals["subset"] == "normal_polyp")
        & (residuals["analysis_id"].isin(["rna_epi_ca_route__atac_wnt_stem", "rna_epi_ca_route__atac_wnt_tcf_ascl2_axis"]))
    ]
    if not best_resid.empty:
        lines.extend(["", "Residualized sensitivity after lesion group, RNA proliferation, and TSS-depth adjustment:"])
        for _, row in best_resid.iterrows():
            lines.append(
                f"- {row['analysis_id']}: n = {int(row['n'])}, residual Spearman rho = {row['residual_spearman_rho']:.3f}, p = {row['p_value']:.3g}, BH q = {row['q_value']:.3g}."
            )

    lines.extend(
        [
            "",
            "Interpretation boundary:",
            "- Positive concordance supports a coupled expression/accessibility layer for the conventional adenoma-WNT route.",
            "- It should be written as sample-level multiome concordance, not chromVAR, motif binding, enhancer looping, or causality.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paired = standardize_inputs()
    corrs = correlation_table(paired)
    models = adjusted_models(paired)
    residuals = residual_correlations(paired)

    paired.to_csv(OUT_DIR / "becker_rna_atac_paired_scores.tsv", sep="\t", index=False)
    corrs.to_csv(OUT_DIR / "becker_rna_atac_correlations.tsv", sep="\t", index=False)
    models.to_csv(OUT_DIR / "becker_rna_atac_adjusted_models.tsv", sep="\t", index=False)
    residuals.to_csv(OUT_DIR / "becker_rna_atac_residual_correlations.tsv", sep="\t", index=False)
    write_summary(OUT_DIR / "becker_rna_atac_concordance_summary.txt", paired, corrs, models, residuals)


if __name__ == "__main__":
    main()
