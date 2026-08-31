#!/usr/bin/env python3
"""Validate frozen consensus modules in five external cohorts and FFPE tissue."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.meta_analysis import combine_effects


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
RESULT_ROOT = ROOT / "results" / "state_aware_program_v1"
SOURCE_ROOT = RESULT_ROOT / "functional_architecture_v1"
OUT_ROOT = Path(os.environ.get("STATE_AWARE_MODULE_RUN_ROOT", str(SOURCE_ROOT))).resolve()
WGCNA_ROOT = SOURCE_ROOT / "consensus_wgcna"
OUT = OUT_ROOT / "module_external_validation"
MODULE_PATH = SOURCE_ROOT / "consensus_modules.tsv"
SUMMARY_PATH = Path(
    os.environ.get(
        "STATE_AWARE_MODULE_SELECTION_PATH",
        str(WGCNA_ROOT / "module_internal_gate_summary.tsv"),
    )
).resolve()
ROUTE_COLUMN = os.environ.get(
    "STATE_AWARE_MODULE_ROUTE_COLUMN", "internal_gate_pass"
)
ADDENDUM_PATH = (
    ANALYSIS
    / "contracts"
    / "state_aware_functional_architecture_downstream_addendum_v1_2026-08-30.md"
)

sys.path.insert(0, str(ANALYSIS))
import external_sporadic_adenoma_validation as external  # noqa: E402
import gse117606_paired_route_validation as ffpe  # noqa: E402
import translation_reduced_panel_v2_0 as reduced  # noqa: E402
import conventional_route_signature_transfer as recurrence  # noqa: E402


SEED = 20260830
EXPECTED_ADDENDUM_HASH = (
    "32343afe117d09007066fbe01f8fbe7cf4a11ee628dbe6a121f0f814968da3bb"
)
PAIRED_COHORTS = {"GSE8671", "GSE72820"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def truth(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return values.astype(str).str.lower().eq("true")


def module_external_gene_data(
    genes: list[str],
) -> tuple[dict[str, tuple[pd.DataFrame, pd.DataFrame]], int]:
    """Load five cohorts, replacing the compact-panel GPL570 map with full annotation."""
    wanted = set(genes) | set(external.PROLIFERATION_CONTROL)
    data = reduced.external_gene_data(pd.DataFrame({"gene": sorted(wanted)}))
    _, probe_to_symbols = recurrence.read_gpl_mapping(
        recurrence.GPL570_ANNOT, wanted
    )
    mapping = pd.DataFrame(
        [
            {"feature_id": probe, "gene": gene}
            for probe, symbols in probe_to_symbols.items()
            for gene in sorted(symbols & wanted)
        ]
    ).drop_duplicates()
    metadata, expression, _ = external.parse_series_matrix("GSE8671")
    metadata = external.add_design_fields("GSE8671", metadata)
    expression = np.log2(expression.clip(lower=0) + 1)
    gene_expression, _ = external.choose_features_without_labels(
        expression, mapping, wanted, "GSE8671"
    )
    module_gene_coverage = len(set(genes) & set(gene_expression.columns)) / len(genes)
    if module_gene_coverage < 0.80:
        raise RuntimeError(
            "Full GPL570 annotation covered fewer than 80% of routed module genes"
        )
    data["GSE8671"] = (metadata, gene_expression)
    return data, len(mapping)


def score_module(
    metadata: pd.DataFrame,
    expression: pd.DataFrame,
    genes: list[str],
    direction_sign: int,
    cohort: str,
) -> tuple[pd.DataFrame | None, dict[str, object]]:
    present = [gene for gene in genes if gene in expression.columns]
    if present:
        standard_deviation = expression[present].astype(float).std(axis=0, ddof=1)
        present = [
            gene
            for gene in present
            if np.isfinite(standard_deviation[gene]) and standard_deviation[gene] > 0
        ]
    coverage = len(present) / len(genes)
    coverage_row: dict[str, object] = {
        "cohort": cohort,
        "module_size": len(genes),
        "n_measurable": len(present),
        "coverage_fraction": coverage,
        "coverage_gate_pass": len(present) >= 15 and coverage >= 0.50,
        "measurable_genes": ";".join(present),
        "missing_or_zero_variance_genes": ";".join(sorted(set(genes) - set(present))),
    }
    if not coverage_row["coverage_gate_pass"]:
        return None, coverage_row

    selected = expression[present].astype(float)
    z = selected.sub(selected.mean(axis=0), axis=1).div(
        selected.std(axis=0, ddof=1).replace(0, np.nan), axis=1
    )
    scores = metadata.copy().set_index("sample_id")
    scores["module_score_raw"] = z.mean(axis=1)
    scores["module_score_oriented"] = direction_sign * scores["module_score_raw"]
    scores = scores.reset_index()
    scores.insert(0, "cohort", cohort)
    scores["patient_cluster_id"] = cohort + "::" + scores["patient_id"].astype(str)
    scores["n_module_genes_measured"] = len(present)
    return scores, coverage_row


def patient_clustered_effect(scores: pd.DataFrame) -> dict[str, object]:
    data = scores.loc[
        scores["tissue_group"].isin(["normal", "adenoma"]),
        ["patient_cluster_id", "tissue_group", "module_score_oriented"],
    ].dropna().copy()
    score_sd = data["module_score_oriented"].std(ddof=1)
    if len(data) < 6 or not np.isfinite(score_sd) or score_sd <= 0:
        raise RuntimeError("External module score has insufficient non-zero observations")
    data["score_sd"] = (
        data["module_score_oriented"] - data["module_score_oriented"].mean()
    ) / score_sd
    data["is_adenoma"] = data["tissue_group"].eq("adenoma").astype(float)
    design = sm.add_constant(data[["is_adenoma"]], has_constant="add")
    fit = sm.OLS(data["score_sd"], design).fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data["patient_cluster_id"],
            "use_correction": True,
            "df_correction": True,
        },
        use_t=True,
    )
    interval = fit.conf_int().loc["is_adenoma"]
    return {
        "n_samples": len(data),
        "n_patient_clusters": data["patient_cluster_id"].nunique(),
        "oriented_adenoma_effect_sd": float(fit.params["is_adenoma"]),
        "standard_error": float(fit.bse["is_adenoma"]),
        "ci_low": float(interval.iloc[0]),
        "ci_high": float(interval.iloc[1]),
        "p_value": float(fit.pvalues["is_adenoma"]),
    }


def random_effects(frame: pd.DataFrame) -> dict[str, float]:
    yi = frame["oriented_adenoma_effect_sd"].to_numpy(float)
    sei = frame["standard_error"].to_numpy(float)
    keep = np.isfinite(yi) & np.isfinite(sei) & (sei > 0)
    yi = yi[keep]
    sei = sei[keep]
    if len(yi) < 2:
        return {
            "k": len(yi),
            "pooled_effect": np.nan,
            "pooled_se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p_value": np.nan,
            "tau2_paule_mandel": np.nan,
            "i2": np.nan,
        }
    fit = combine_effects(yi, sei**2, method_re="iterated", use_t=True)
    estimate = float(fit.mean_effect_re)
    standard_error = float(fit.sd_eff_w_re)
    critical = float(stats.t.ppf(0.975, df=len(yi) - 1))
    p_value = float(
        2 * stats.t.sf(abs(estimate / standard_error), df=len(yi) - 1)
    )
    return {
        "k": len(yi),
        "pooled_effect": estimate,
        "pooled_se": standard_error,
        "ci_low": estimate - critical * standard_error,
        "ci_high": estimate + critical * standard_error,
        "p_value": p_value,
        "tau2_paule_mandel": float(max(fit.tau2, 0)),
        "i2": float(max(fit.i2, 0)),
    }


def paired_ffpe_effect(scores: pd.DataFrame) -> dict[str, object]:
    subset = scores.loc[scores["tissue_group"].isin(["normal", "adenoma"])].copy()
    wide = subset.pivot_table(
        index="patient_id",
        columns="tissue_group",
        values="module_score_oriented",
        aggfunc="median",
    )
    pairs = wide[["adenoma", "normal"]].dropna()
    difference = pairs["adenoma"] - pairs["normal"]
    test = (
        stats.wilcoxon(difference, zero_method="wilcox", alternative="two-sided")
        if len(difference) and (difference != 0).any()
        else None
    )
    return {
        "n_pairs": len(difference),
        "median_oriented_paired_difference": float(difference.median()),
        "mean_oriented_paired_difference": float(difference.mean()),
        "positive_pair_fraction": float((difference > 0).mean()),
        "p_paired_wilcoxon": float(test.pvalue) if test is not None else np.nan,
    }


def main() -> None:
    if sha256(ADDENDUM_PATH) != EXPECTED_ADDENDUM_HASH:
        raise RuntimeError("The frozen downstream validation addendum changed")
    modules = pd.read_csv(MODULE_PATH, sep="\t")
    summary = pd.read_csv(SUMMARY_PATH, sep="\t")
    if ROUTE_COLUMN not in summary.columns:
        raise RuntimeError(f"Routing column is missing: {ROUTE_COLUMN}")
    eligible_modules = summary.loc[truth(summary[ROUTE_COLUMN]), "module"].tolist()
    OUT.mkdir(parents=True, exist_ok=True)

    if not eligible_modules:
        summary.assign(
            external_coverage_cohorts=0,
            external_gate_pass=False,
            external_gate_reason=f"no module passed routing column {ROUTE_COLUMN}",
        ).to_csv(OUT / "module_external_validation_audit.tsv", sep="\t", index=False)
        (OUT / "module_external_validation_manifest.json").write_text(
            json.dumps(
                {
                    "analysis": "state_aware_module_external_validation_v1",
                    "created_utc": pd.Timestamp.utcnow().isoformat(),
                    "route_eligible_modules": [],
                    "route_column": ROUTE_COLUMN,
                    "external_scoring_performed": False,
                    "input_sha256": {
                        "consensus_modules": sha256(MODULE_PATH),
                        "module_selection": sha256(SUMMARY_PATH),
                        "downstream_addendum": sha256(ADDENDUM_PATH),
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        print("No internally eligible module; external module scoring was not run")
        return

    module_gene_sets = {
        module: modules.loc[modules["module"].eq(module), "gene"].astype(str).tolist()
        for module in eligible_modules
    }
    direction_sign = {
        row.module: 1 if row.heldout_direction == "Up" else -1
        for row in summary.loc[summary["module"].isin(eligible_modules)].itertuples()
    }
    all_genes = sorted(set().union(*map(set, module_gene_sets.values())))
    gene_data, full_gpl570_mapping_rows = module_external_gene_data(all_genes)

    score_parts: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, object]] = []
    effect_rows: list[dict[str, object]] = []
    gene_presence_rows: list[dict[str, object]] = []
    for cohort, (metadata, expression) in gene_data.items():
        for gene in all_genes:
            present = gene in expression.columns and expression[gene].std(ddof=1) > 0
            gene_presence_rows.append(
                {"gene": gene, "platform": cohort, "measurable": bool(present)}
            )
        for module in eligible_modules:
            scores, coverage = score_module(
                metadata,
                expression,
                module_gene_sets[module],
                direction_sign[module],
                cohort,
            )
            coverage["module"] = module
            coverage_rows.append(coverage)
            if scores is None:
                continue
            scores["module"] = module
            score_parts.append(scores)
            effect_rows.append(
                {
                    "module": module,
                    "cohort": cohort,
                    "paired_design": cohort in PAIRED_COHORTS,
                    **patient_clustered_effect(scores),
                }
            )

    effects = pd.DataFrame(effect_rows)
    meta_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    for module in eligible_modules:
        local = effects.loc[effects["module"].eq(module)].copy()
        pooled = random_effects(local)
        local_loo = []
        if len(local) >= 4:
            for cohort in sorted(local["cohort"].unique()):
                value = random_effects(local.loc[local["cohort"].ne(cohort)])
                row = {"module": module, "excluded_cohort": cohort, **value}
                loo_rows.append(row)
                local_loo.append(row)
        all_loo_positive = bool(local_loo) and all(
            row["pooled_effect"] > 0 for row in local_loo
        )
        meta_rows.append(
            {
                "module": module,
                "n_coverage_eligible_cohorts": len(local),
                **pooled,
                "all_leave_one_cohort_out_positive": all_loo_positive,
                "external_gate_pass": len(local) >= 3
                and pooled["ci_low"] > 0
                and all_loo_positive,
            }
        )

    ffpe_metadata = ffpe.parse_metadata()
    ffpe_expression_raw = ffpe.parse_expression()
    if list(ffpe_expression_raw.index) != ffpe_metadata["sample_id"].tolist():
        raise RuntimeError("GSE117606 expression and metadata order mismatch")
    ffpe_mapping = ffpe.map_symbols_to_features(set(all_genes), ffpe_expression_raw)
    ffpe_expression = ffpe.select_gene_expression(ffpe_expression_raw, ffpe_mapping)
    ffpe_rows: list[dict[str, object]] = []
    for gene in all_genes:
        present = gene in ffpe_expression.columns and ffpe_expression[gene].std(ddof=1) > 0
        gene_presence_rows.append(
            {"gene": gene, "platform": "GSE117606", "measurable": bool(present)}
        )
    for module in eligible_modules:
        scores, coverage = score_module(
            ffpe_metadata,
            ffpe_expression,
            module_gene_sets[module],
            direction_sign[module],
            "GSE117606",
        )
        coverage["module"] = module
        coverage_rows.append(coverage)
        if scores is None:
            ffpe_rows.append(
                {
                    "module": module,
                    "coverage_gate_pass": False,
                    "n_pairs": 0,
                }
            )
            continue
        scores["module"] = module
        score_parts.append(scores)
        ffpe_rows.append(
            {
                "module": module,
                "coverage_gate_pass": True,
                **paired_ffpe_effect(scores),
            }
        )

    meta_table = pd.DataFrame(meta_rows)
    validation = summary.merge(meta_table, on="module", how="left")
    validation["external_gate_pass"] = validation["external_gate_pass"].fillna(False)
    validation["routing_status"] = np.where(
        truth(validation[ROUTE_COLUMN]) & validation["external_gate_pass"],
        "pending_technical_orthogonal_and_perturbation_gates",
        np.where(
            truth(validation[ROUTE_COLUMN]),
            "audit_external_gate_failed",
            "audit_analysis_route_failed",
        ),
    )

    pd.concat(score_parts, ignore_index=True, sort=False).to_csv(
        OUT / "module_sample_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(coverage_rows).to_csv(
        OUT / "module_platform_coverage.tsv", sep="\t", index=False
    )
    effects.to_csv(OUT / "module_external_cohort_effects.tsv", sep="\t", index=False)
    meta_table.to_csv(
        OUT / "module_external_random_effects_paule_mandel_audit.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(loo_rows).to_csv(
        OUT / "module_external_leave_one_cohort_out.tsv", sep="\t", index=False
    )
    pd.DataFrame(ffpe_rows).to_csv(
        OUT / "module_ffpe_paired_effects.tsv", sep="\t", index=False
    )
    pd.DataFrame(gene_presence_rows).to_csv(
        OUT / "module_gene_platform_measurability.tsv", sep="\t", index=False
    )
    validation.to_csv(
        OUT / "module_external_validation_paule_mandel_audit.tsv",
        sep="\t",
        index=False,
    )

    manifest = {
        "analysis": "state_aware_module_external_validation_v1",
        "created_utc": pd.Timestamp.utcnow().isoformat(),
        "random_seed": SEED,
        "route_eligible_modules": eligible_modules,
        "route_column": ROUTE_COLUMN,
        "full_gpl570_mapping_rows": full_gpl570_mapping_rows,
        "modules_passing_external_gate": meta_table.loc[
            meta_table["external_gate_pass"], "module"
        ].tolist(),
        "module_score": "direction-oriented mean of all measurable cohort-z-scored genes",
        "coverage_gate": "at least 50% and at least 15 genes",
        "meta_analysis": (
            "Paule-Mandel random effects retained as audit only; "
            "the frozen REML-Knapp-Hartung gate is applied by "
            "state_aware_module_meta_analysis_v1.R"
        ),
        "input_sha256": {
            "consensus_modules": sha256(MODULE_PATH),
            "module_selection": sha256(SUMMARY_PATH),
            "downstream_addendum": sha256(ADDENDUM_PATH),
            "gpl570_full_annotation": sha256(recurrence.GPL570_ANNOT),
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "statsmodels": platform.python_implementation()
            + " "
            + __import__("statsmodels").__version__,
        },
    }
    (OUT / "module_external_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"External module validation complete: {len(eligible_modules)} route eligible; "
        f"{int(meta_table['external_gate_pass'].sum())} passed the external gate"
    )


if __name__ == "__main__":
    main()
