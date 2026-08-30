#!/usr/bin/env python3
"""Validate the frozen state-shared programme and eight-gene readout externally."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
RESULT_ROOT = ROOT / "results" / "state_aware_program_v1"
OUT = RESULT_ROOT / "external_validation"
COMMON_PATH = RESULT_ROOT / "common_effects" / "cross_state_common_effects.tsv.gz"
PANEL_PATH = RESULT_ROOT / "panel_derivation" / "compact_state_shared_panel_frozen.tsv"
ADDENDUM_PATH = (
    ANALYSIS / "contracts" / "state_shared_external_validation_addendum_v1_2026-08-29.md"
)
sys.path.insert(0, str(ANALYSIS))

import external_sporadic_adenoma_validation as external  # noqa: E402
import gse117606_paired_route_validation as ffpe  # noqa: E402
import translation_reduced_panel_v2_0 as reduced  # noqa: E402


EXPECTED_HASHES = {
    "common_effects": "a1ac4b7b67ac279782e04e386971d7463e169cbdbc24d7da4c314e34a4f3e946",
    "compact_panel": "c5997e572342a72da8441df312fba4e3461cacfa5e30d0d8590a3f23ae3d96f0",
    "external_addendum": "ea0ac55c402ba70fe840b721ae878bf8a4d2028593bdf3322a894b86bd005b43",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write(frame: pd.DataFrame, filename: str, compress: bool = False) -> None:
    frame.to_csv(
        OUT / filename,
        sep="\t",
        index=False,
        compression="gzip" if compress else None,
    )


def load_signatures() -> tuple[pd.DataFrame, pd.DataFrame]:
    common = pd.read_csv(COMMON_PATH, sep="\t")
    strict_flag = common["strict_state_shared"]
    if strict_flag.dtype != bool:
        strict_flag = strict_flag.astype(str).str.lower().eq("true")
    full = common.loc[strict_flag, ["gene", "shared_direction"]].copy()
    full["arm"] = full["shared_direction"]
    full["route_weight"] = np.where(full["arm"].eq("up"), 1.0, -1.0)
    full["signature_id"] = "state_shared_1843"
    panel = pd.read_csv(PANEL_PATH, sep="\t").copy()
    panel["signature_id"] = "compact_8"
    if len(full) != 1843 or full["arm"].value_counts().to_dict() != {"down": 959, "up": 884}:
        raise RuntimeError("Frozen full programme has unexpected dimensions")
    if len(panel) != 8 or panel["arm"].value_counts().to_dict() != {"up": 4, "down": 4}:
        raise RuntimeError("Frozen compact readout has unexpected dimensions")
    if panel["validation_outcomes_used"].astype(str).str.lower().ne("false").any():
        raise RuntimeError("Compact readout records use of validation outcomes")
    return full, panel


def score_expression(
    meta: pd.DataFrame,
    expression: pd.DataFrame,
    signature: pd.DataFrame,
    cohort: str,
    require_complete: bool,
) -> tuple[pd.DataFrame, dict[str, object]]:
    up_expected = signature.loc[signature["route_weight"].eq(1), "gene"].astype(str).tolist()
    down_expected = signature.loc[signature["route_weight"].eq(-1), "gene"].astype(str).tolist()
    up_present = [gene for gene in up_expected if gene in expression.columns]
    down_present = [gene for gene in down_expected if gene in expression.columns]
    structural_missing = sorted(
        set(up_expected + down_expected).difference(expression.columns)
    )
    standard_deviation = expression[up_present + down_present].astype(float).std(
        axis=0,
        ddof=1,
    )
    zero_variance = sorted(
        standard_deviation.index[
            ~np.isfinite(standard_deviation) | standard_deviation.eq(0)
        ].tolist()
    )
    if require_complete and (structural_missing or zero_variance):
        raise RuntimeError(
            f"{cohort}: compact readout unmeasurable genes: "
            f"missing={structural_missing}; zero_variance={zero_variance}"
        )
    up = [gene for gene in up_present if gene not in zero_variance]
    down = [gene for gene in down_present if gene not in zero_variance]
    if not up or not down:
        raise RuntimeError(f"{cohort}: at least one programme arm is unmeasurable")
    selected = expression[up + down].astype(float)
    z = selected.sub(selected.mean(axis=0), axis=1).div(
        selected.std(axis=0, ddof=1).replace(0, np.nan),
        axis=1,
    )
    output = meta.copy().set_index("sample_id")
    output["route_score_k50"] = z[up].mean(axis=1) - z[down].mean(axis=1)
    output = output.reset_index()
    output.insert(0, "cohort", cohort)
    output["patient_cluster_id"] = cohort + "::" + output["patient_id"].astype(str)
    coverage = {
        "cohort": cohort,
        "signature_id": signature["signature_id"].iloc[0],
        "expected_up": len(up_expected),
        "measured_up": len(up),
        "up_coverage": len(up) / len(up_expected),
        "expected_down": len(down_expected),
        "measured_down": len(down),
        "down_coverage": len(down) / len(down_expected),
        "structurally_missing_genes": ";".join(structural_missing),
        "zero_variance_genes": ";".join(zero_variance),
    }
    return output, coverage


def external_layers(
    full: pd.DataFrame,
    panel: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    gene_data = reduced.external_gene_data(full[["gene"]])
    base_controls = pd.read_csv(
        ROOT / "results" / "external_sporadic_adenoma_validation" / "sample_scores.tsv",
        sep="\t",
    )[
        [
            "cohort",
            "sample_id",
            "score__proliferation_control",
            "n_proliferation_control_present",
        ]
    ]
    paired_cohorts = {"GSE8671", "GSE72820"}
    score_parts: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, object]] = []
    test_rows: list[dict[str, object]] = []
    pooled_rows: list[dict[str, object]] = []
    adjusted_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []

    for signature, require_complete in ((full, False), (panel, True)):
        signature_id = signature["signature_id"].iloc[0]
        signature_scores: list[pd.DataFrame] = []
        for cohort, (meta, expression) in gene_data.items():
            scores, coverage = score_expression(
                meta,
                expression,
                signature,
                cohort,
                require_complete,
            )
            scores["signature_id"] = signature_id
            signature_scores.append(scores)
            coverage_rows.append(coverage)
            comparison = external.cohort_comparison(
                scores,
                cohort,
                50,
                "adenoma",
                "normal",
                cohort in paired_cohorts,
                "adenoma_vs_normal",
            )
            comparison["expected_up"] = coverage["expected_up"]
            comparison["measured_up"] = coverage["measured_up"]
            comparison["expected_down"] = coverage["expected_down"]
            comparison["measured_down"] = coverage["measured_down"]
            comparison["score_definition"] = "mean(all measurable up) - mean(all measurable down)"
            comparison.pop("signature_size_per_direction", None)
            comparison["signature_id"] = signature_id
            test_rows.append(comparison)
        scores = pd.concat(signature_scores, ignore_index=True, sort=False)
        scores = scores.merge(
            base_controls,
            on=["cohort", "sample_id"],
            how="left",
            validate="one_to_one",
        )
        score_parts.append(scores)
        pooled = external.one_stage_model(scores, 50)
        pooled["expected_up"] = int(signature["route_weight"].eq(1).sum())
        pooled["expected_down"] = int(signature["route_weight"].eq(-1).sum())
        pooled["score_definition"] = "cohort-specific all-measurable programme score"
        pooled.pop("signature_size_per_direction", None)
        pooled["signature_id"] = signature_id
        pooled_rows.append(pooled)
        adjusted = external.proliferation_adjusted_model(scores)
        adjusted["expected_up"] = int(signature["route_weight"].eq(1).sum())
        adjusted["expected_down"] = int(signature["route_weight"].eq(-1).sum())
        adjusted["score_definition"] = "cohort-specific all-measurable programme score"
        adjusted.pop("signature_size_per_direction", None)
        adjusted["signature_id"] = signature_id
        adjusted_rows.append(adjusted)
        for excluded in sorted(scores["cohort"].unique()):
            current = external.proliferation_adjusted_model(
                scores,
                excluded_cohort=excluded,
            )
            current["expected_up"] = int(signature["route_weight"].eq(1).sum())
            current["expected_down"] = int(signature["route_weight"].eq(-1).sum())
            current["score_definition"] = "cohort-specific all-measurable programme score"
            current.pop("signature_size_per_direction", None)
            current["signature_id"] = signature_id
            loo_rows.append(current)

    all_scores = pd.concat(score_parts, ignore_index=True, sort=False)
    full_scores = all_scores.loc[
        all_scores["signature_id"].eq("state_shared_1843"),
        ["cohort", "sample_id", "route_score_k50"],
    ].rename(columns={"route_score_k50": "full_score"})
    compact_scores = all_scores.loc[
        all_scores["signature_id"].eq("compact_8"),
        ["cohort", "sample_id", "route_score_k50"],
    ].rename(columns={"route_score_k50": "compact_score"})
    fidelity_data = full_scores.merge(
        compact_scores,
        on=["cohort", "sample_id"],
        validate="one_to_one",
    )
    fidelity_rows = []
    for cohort, frame in fidelity_data.groupby("cohort", sort=False):
        test = stats.spearmanr(frame["compact_score"], frame["full_score"])
        fidelity_rows.append(
            {
                "cohort": cohort,
                "n_samples": len(frame),
                "spearman_compact_vs_full": float(test.statistic),
                "p_value": float(test.pvalue),
            }
        )
    return {
        "scores": all_scores,
        "coverage": pd.DataFrame(coverage_rows),
        "tests": pd.DataFrame(test_rows),
        "pooled": pd.DataFrame(pooled_rows),
        "adjusted": pd.DataFrame(adjusted_rows),
        "loo": pd.DataFrame(loo_rows),
        "fidelity": pd.DataFrame(fidelity_rows),
    }


def ffpe_layers(
    full: pd.DataFrame,
    panel: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    meta = ffpe.parse_metadata()
    expression = ffpe.parse_expression()
    if list(expression.index) != meta["sample_id"].tolist():
        raise RuntimeError("GSE117606 expression and metadata order mismatch")
    wanted = set(full["gene"]) | set(external.PROLIFERATION_CONTROL)
    mapping = ffpe.map_symbols_to_features(wanted, expression)
    gene_expression = ffpe.select_gene_expression(expression, mapping)
    score_parts: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, object]] = []
    test_rows: list[dict[str, object]] = []
    for signature, require_complete in ((full, False), (panel, True)):
        scores, coverage = score_expression(
            meta,
            gene_expression,
            signature,
            "GSE117606",
            require_complete,
        )
        signature_id = signature["signature_id"].iloc[0]
        scores["signature_id"] = signature_id
        score_parts.append(scores)
        coverage_rows.append(coverage)
        comparison = external.cohort_comparison(
            scores,
            "GSE117606",
            50,
            "adenoma",
            "normal",
            True,
            "conventional_adenoma_vs_adjacent_mucosa",
        )
        comparison["expected_up"] = coverage["expected_up"]
        comparison["measured_up"] = coverage["measured_up"]
        comparison["expected_down"] = coverage["expected_down"]
        comparison["measured_down"] = coverage["measured_down"]
        comparison["score_definition"] = "mean(all measurable up) - mean(all measurable down)"
        comparison.pop("signature_size_per_direction", None)
        comparison["signature_id"] = signature_id
        test_rows.append(comparison)
    scores = pd.concat(score_parts, ignore_index=True, sort=False)
    full_scores = scores.loc[
        scores["signature_id"].eq("state_shared_1843"),
        ["sample_id", "route_score_k50"],
    ].rename(columns={"route_score_k50": "full_score"})
    compact_scores = scores.loc[
        scores["signature_id"].eq("compact_8"),
        ["sample_id", "route_score_k50"],
    ].rename(columns={"route_score_k50": "compact_score"})
    fidelity_data = full_scores.merge(compact_scores, on="sample_id", validate="one_to_one")
    fidelity_test = stats.spearmanr(
        fidelity_data["compact_score"],
        fidelity_data["full_score"],
    )
    fidelity = pd.DataFrame(
        [
            {
                "cohort": "GSE117606",
                "n_samples": len(fidelity_data),
                "spearman_compact_vs_full": float(fidelity_test.statistic),
                "p_value": float(fidelity_test.pvalue),
            }
        ]
    )
    candidates = panel[["gene", "route_weight"]].rename(
        columns={"route_weight": "expected_direction"}
    )
    candidates["direction_stability"] = panel["donor_heldout_selection_frequency"].to_numpy()
    gene_tests = ffpe.paired_gene_tests(meta, gene_expression, candidates)
    return {
        "scores": scores,
        "coverage": pd.DataFrame(coverage_rows),
        "tests": pd.DataFrame(test_rows),
        "fidelity": fidelity,
        "gene_tests": gene_tests,
    }


def perturbation_layers(
    full: pd.DataFrame,
    panel: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    score_parts: list[pd.DataFrame] = []
    apc_parts: list[pd.DataFrame] = []
    tcf_parts: list[pd.DataFrame] = []
    for signature in (full, panel):
        signature_id = signature["signature_id"].iloc[0]
        scores, apc, tcf = reduced.perturbation_analysis(signature)
        scores["signature_id"] = signature_id
        apc["signature_id"] = signature_id
        tcf["signature_id"] = signature_id
        score_parts.append(scores)
        apc_parts.append(apc)
        tcf_parts.append(tcf)
    return {
        "scores": pd.concat(score_parts, ignore_index=True, sort=False),
        "apc": pd.concat(apc_parts, ignore_index=True, sort=False),
        "tcf": pd.concat(tcf_parts, ignore_index=True, sort=False),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    actual_hashes = {
        "common_effects": sha256(COMMON_PATH),
        "compact_panel": sha256(PANEL_PATH),
        "external_addendum": sha256(ADDENDUM_PATH),
    }
    if actual_hashes != EXPECTED_HASHES:
        raise RuntimeError("A frozen external-validation input changed")
    full, panel = load_signatures()

    external_result = external_layers(full, panel)
    ffpe_result = ffpe_layers(full, panel)
    perturbation_result = perturbation_layers(full, panel)

    external_scores = external_result["scores"].rename(
        columns={"route_score_k50": "programme_score"}
    )
    write(external_scores, "external_sample_scores.tsv.gz", compress=True)
    write(external_result["coverage"], "external_gene_coverage.tsv")
    write(external_result["tests"], "external_cohort_tests.tsv")
    write(external_result["pooled"], "external_pooled_models.tsv")
    write(external_result["adjusted"], "external_proliferation_adjusted_models.tsv")
    write(external_result["loo"], "external_adjusted_leave_one_cohort_out.tsv")
    write(external_result["fidelity"], "external_compact_full_fidelity.tsv")
    ffpe_scores = ffpe_result["scores"].rename(
        columns={"route_score_k50": "programme_score"}
    )
    write(ffpe_scores, "ffpe_sample_scores.tsv.gz", compress=True)
    write(ffpe_result["coverage"], "ffpe_gene_coverage.tsv")
    write(ffpe_result["tests"], "ffpe_paired_tests.tsv")
    write(ffpe_result["fidelity"], "ffpe_compact_full_fidelity.tsv")
    write(ffpe_result["gene_tests"], "ffpe_compact_gene_tests.tsv")
    write(perturbation_result["scores"], "perturbation_sample_scores.tsv")
    write(perturbation_result["apc"], "apc_organoid_effects.tsv")
    write(perturbation_result["tcf"], "tcf7l2_clone_effects.tsv")

    manifest = {
        "analysis": "validate_state_shared_external_layers_v1",
        "analysis_date": "2026-08-29",
        "frozen_input_hashes": actual_hashes,
        "full_programme_genes": len(full),
        "compact_readout_genes": len(panel),
        "gene_reselection": False,
        "weight_fitting": False,
        "cutpoint_optimization": False,
        "score": "within-dataset gene z scores; mean(up) minus mean(down)",
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    (OUT / "external_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"output_dir\t{OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
