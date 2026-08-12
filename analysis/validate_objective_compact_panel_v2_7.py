#!/usr/bin/env python3
"""Validate the frozen objective compact panel without gene reselection.

The 12-gene panel is loaded from the discovery-only R selection output before
any validation expression matrix is opened.  The previous biology-guided
10-gene panel is carried only as an internal comparator and reproducibility
control; the 100-gene score is not a candidate in this analysis.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "analysis"
OUT_DIR = ROOT / "results" / "objective_compact_panel_v2_7"
sys.path.insert(0, str(ANALYSIS_DIR))

import conventional_route_signature_transfer as transfer  # noqa: E402
import discovery_locked_route_signature as discovery  # noqa: E402
import external_sporadic_adenoma_validation as external  # noqa: E402
import gse117606_paired_route_validation as ffpe  # noqa: E402
import translation_reduced_panel_v2_0 as reduced  # noqa: E402


OBJECTIVE_PATH = OUT_DIR / "objective_compact_panel_frozen.tsv"
CORE_PATH = (
    ROOT
    / "results"
    / "data_adaptive_panel_pilot_v2_6"
    / "stable_error_controlled_core.tsv"
)
COMPARATOR_PATH = (
    ROOT
    / "results"
    / "programme_transparency_v2_5"
    / "compact_panel_definition_corrected.tsv"
)
PANEL_ORDER = ["objective_12", "biology_guided_10"]
PANEL_LABELS = {
    "objective_12": "Objective Kneedle panel (12)",
    "biology_guided_10": "Previous biology-guided panel (10)",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, name: str, compress: bool = False) -> None:
    frame.to_csv(
        OUT_DIR / name,
        sep="\t",
        index=False,
        compression="gzip" if compress else None,
    )


def load_fixed_panels() -> tuple[pd.DataFrame, pd.DataFrame]:
    objective = pd.read_csv(OBJECTIVE_PATH, sep="\t")
    if len(objective) != 12 or objective["validation_outcomes_used"].astype(str).str.lower().ne("false").any():
        raise RuntimeError("Objective panel is not the expected frozen discovery-only 12-gene panel")
    objective = objective[
        ["gene", "arm", "route_weight", "pair_step", "direction_stability"]
    ].copy()
    objective["panel_id"] = "objective_12"
    objective["panel_label"] = PANEL_LABELS["objective_12"]

    comparator_raw = pd.read_csv(COMPARATOR_PATH, sep="\t")
    comparator = comparator_raw[
        ["gene", "route_weight", "panel_order_within_arm", "direction_stability"]
    ].copy()
    comparator["arm"] = np.where(comparator["route_weight"].eq(1), "up", "down")
    comparator = comparator.rename(columns={"panel_order_within_arm": "pair_step"})
    comparator["panel_id"] = "biology_guided_10"
    comparator["panel_label"] = PANEL_LABELS["biology_guided_10"]

    panels = pd.concat([objective, comparator], ignore_index=True, sort=False)
    core = pd.read_csv(CORE_PATH, sep="\t")
    core_panel = core[["gene", "arm"]].copy()
    core_panel["route_weight"] = np.where(core_panel["arm"].eq("up"), 1.0, -1.0)
    if len(core_panel) != 287:
        raise RuntimeError("Expected the fixed 287-gene discovery core")
    return panels, core_panel


def spearman(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    keep = x.notna() & y.notna()
    if keep.sum() < 3:
        return int(keep.sum()), np.nan, np.nan
    test = stats.spearmanr(x.loc[keep], y.loc[keep])
    return int(keep.sum()), float(test.statistic), float(test.pvalue)


def chen_validation(
    panels: pd.DataFrame, core_panel: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score_parts: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []
    correlation_rows: list[dict[str, object]] = []
    for dataset in ["discovery", "validation"]:
        meta, expression = discovery.load_chen_pseudobulk(dataset)
        core_score = transfer.route_score_from_expression(
            expression,
            core_panel.loc[core_panel["arm"].eq("up"), "gene"].tolist(),
            core_panel.loc[core_panel["arm"].eq("down"), "gene"].tolist(),
        )["score__ca_route_signature"]
        core_frame = meta.copy()
        core_frame["panel_id"] = "full_core_287"
        core_frame["panel_label"] = "Full error-controlled core (287)"
        core_frame["panel_score"] = core_score.to_numpy()
        score_parts.append(core_frame)
        core_metrics = reduced.heldout_metrics(meta, expression, core_panel)
        metric_rows.append(
            {
                "dataset": dataset,
                "panel_id": "full_core_287",
                "panel_label": "Full error-controlled core (287)",
                "n_genes": 287,
                **core_metrics,
            }
        )

        for panel_id in PANEL_ORDER:
            panel = panels.loc[panels["panel_id"].eq(panel_id)].copy()
            result = transfer.route_score_from_expression(
                expression,
                panel.loc[panel["arm"].eq("up"), "gene"].tolist(),
                panel.loc[panel["arm"].eq("down"), "gene"].tolist(),
            )
            frame = meta.copy()
            frame["panel_id"] = panel_id
            frame["panel_label"] = PANEL_LABELS[panel_id]
            frame["panel_score"] = result["score__ca_route_signature"].to_numpy()
            score_parts.append(frame)
            metrics = reduced.heldout_metrics(meta, expression, panel)
            metric_rows.append(
                {
                    "dataset": dataset,
                    "panel_id": panel_id,
                    "panel_label": PANEL_LABELS[panel_id],
                    "n_genes": len(panel),
                    **metrics,
                }
            )
            n, rho, p_value = spearman(frame["panel_score"], core_score)
            correlation_rows.append(
                {
                    "dataset": dataset,
                    "panel_id": panel_id,
                    "panel_label": PANEL_LABELS[panel_id],
                    "n": n,
                    "spearman_rho_vs_full_287_core": rho,
                    "p_value": p_value,
                }
            )
    return (
        pd.concat(score_parts, ignore_index=True, sort=False),
        pd.DataFrame(metric_rows),
        pd.DataFrame(correlation_rows),
    )


def external_validation(
    panels: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    union = panels[["gene"]].drop_duplicates()
    gene_data = reduced.external_gene_data(union)
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
    score_parts: list[pd.DataFrame] = []
    test_rows: list[dict[str, object]] = []
    pooled_rows: list[dict[str, object]] = []
    adjusted_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    paired = {"GSE8671", "GSE72820"}

    for panel_id in PANEL_ORDER:
        panel = panels.loc[panels["panel_id"].eq(panel_id)].copy()
        scores = pd.concat(
            [
                reduced.score_gene_data(meta, expression, panel, cohort)
                for cohort, (meta, expression) in gene_data.items()
            ],
            ignore_index=True,
            sort=False,
        )
        scores = scores.merge(
            base_controls,
            on=["cohort", "sample_id"],
            how="left",
            validate="one_to_one",
        )
        scores["route_score_k50"] = scores["route_score_k5"]
        scores["panel_id"] = panel_id
        scores["panel_label"] = PANEL_LABELS[panel_id]
        score_parts.append(scores)
        for cohort, frame in scores.groupby("cohort", sort=False):
            row = external.cohort_comparison(
                frame,
                cohort,
                50,
                "adenoma",
                "normal",
                cohort in paired,
                "adenoma_vs_normal",
            )
            row.update({"panel_id": panel_id, "panel_label": PANEL_LABELS[panel_id]})
            test_rows.append(row)
        pooled = external.one_stage_model(scores, 50)
        pooled.update({"panel_id": panel_id, "panel_label": PANEL_LABELS[panel_id]})
        pooled_rows.append(pooled)
        adjusted = external.proliferation_adjusted_model(scores)
        adjusted.update({"panel_id": panel_id, "panel_label": PANEL_LABELS[panel_id]})
        adjusted_rows.append(adjusted)
        for excluded in sorted(scores["cohort"].unique()):
            row = external.proliferation_adjusted_model(scores, excluded_cohort=excluded)
            row.update({"panel_id": panel_id, "panel_label": PANEL_LABELS[panel_id]})
            loo_rows.append(row)
    return (
        pd.concat(score_parts, ignore_index=True, sort=False),
        pd.DataFrame(test_rows),
        pd.DataFrame(pooled_rows),
        pd.DataFrame(adjusted_rows),
        pd.DataFrame(loo_rows),
    )


def ffpe_validation(
    panels: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score_parts: list[pd.DataFrame] = []
    test_rows: list[pd.DataFrame] = []
    gene_rows: list[pd.DataFrame] = []
    for panel_id in PANEL_ORDER:
        panel = panels.loc[panels["panel_id"].eq(panel_id)].copy()
        scores, test, _, gene_expression = reduced.ffpe_analysis(panel)
        scores["panel_id"] = panel_id
        scores["panel_label"] = PANEL_LABELS[panel_id]
        score_parts.append(scores)
        test["panel_id"] = panel_id
        test["panel_label"] = PANEL_LABELS[panel_id]
        test_rows.append(test)

        candidates = panel[["gene", "route_weight", "direction_stability"]].rename(
            columns={"route_weight": "expected_direction"}
        )
        gene_test = ffpe.paired_gene_tests(ffpe.parse_metadata(), gene_expression, candidates)
        gene_test["panel_id"] = panel_id
        gene_test["panel_label"] = PANEL_LABELS[panel_id]
        gene_rows.append(gene_test)
    return (
        pd.concat(score_parts, ignore_index=True, sort=False),
        pd.concat(test_rows, ignore_index=True, sort=False),
        pd.concat(gene_rows, ignore_index=True, sort=False),
    )


def perturbation_validation(panels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score_parts: list[pd.DataFrame] = []
    apc_parts: list[pd.DataFrame] = []
    tcf_parts: list[pd.DataFrame] = []
    for panel_id in PANEL_ORDER:
        panel = panels.loc[panels["panel_id"].eq(panel_id)].copy()
        scores, apc, tcf = reduced.perturbation_analysis(panel)
        scores["panel_id"] = panel_id
        scores["panel_label"] = PANEL_LABELS[panel_id]
        score_parts.append(scores)
        apc["panel_id"] = panel_id
        apc["panel_label"] = PANEL_LABELS[panel_id]
        apc_parts.append(apc)
        tcf["panel_id"] = panel_id
        tcf["panel_label"] = PANEL_LABELS[panel_id]
        tcf_parts.append(tcf)
    return (
        pd.concat(score_parts, ignore_index=True, sort=False),
        pd.concat(apc_parts, ignore_index=True, sort=False),
        pd.concat(tcf_parts, ignore_index=True, sort=False),
    )


def qa_checks(
    panels: pd.DataFrame,
    chen_metrics: pd.DataFrame,
    external_pooled: pd.DataFrame,
    ffpe_tests: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add(check: str, passed: bool, observed: object, expected: object) -> None:
        rows.append(
            {"check": check, "passed": bool(passed), "observed": observed, "expected": expected}
        )

    add("objective_panel_frozen_at_12_genes", (panels["panel_id"] == "objective_12").sum() == 12, (panels["panel_id"] == "objective_12").sum(), 12)
    add("objective_panel_balanced", panels.query("panel_id == 'objective_12'").groupby("arm")["gene"].size().nunique() == 1, panels.query("panel_id == 'objective_12'").groupby("arm")["gene"].size().to_dict(), "6 up / 6 down")

    existing_chen = pd.read_csv(
        ROOT / "results" / "translation_reduced_panel_v2_0" / "chen_reduced_performance.tsv",
        sep="\t",
    )
    observed_chen = float(
        chen_metrics.query("dataset == 'validation' and panel_id == 'biology_guided_10'")["auc"].iloc[0]
    )
    expected_chen = float(existing_chen.query("dataset == 'validation'")["auc"].iloc[0])
    add("biology_guided_10_chen_auc_reproduced", abs(observed_chen - expected_chen) < 1e-12, observed_chen, expected_chen)

    observed_external = float(
        external_pooled.query("panel_id == 'biology_guided_10'")["adenoma_coef_sd"].iloc[0]
    )
    expected_external = float(
        pd.read_csv(
            ROOT / "results" / "translation_reduced_panel_v2_0" / "external_reduced_pooled_model.tsv",
            sep="\t",
        )["adenoma_coef_sd"].iloc[0]
    )
    add("biology_guided_10_external_effect_reproduced", abs(observed_external - expected_external) < 1e-12, observed_external, expected_external)

    observed_ffpe = float(
        ffpe_tests.query("panel_id == 'biology_guided_10'")["median_paired_difference"].iloc[0]
    )
    expected_ffpe = float(
        pd.read_csv(
            ROOT / "results" / "translation_reduced_panel_v2_0" / "ffpe_reduced_paired_test.tsv",
            sep="\t",
        )["median_paired_difference"].iloc[0]
    )
    add("biology_guided_10_ffpe_delta_reproduced", abs(observed_ffpe - expected_ffpe) < 1e-12, observed_ffpe, expected_ffpe)
    return pd.DataFrame(rows)


def main() -> None:
    # The frozen-panel checksum is recorded before any validation loader runs.
    frozen_checksum = sha256(OBJECTIVE_PATH)
    panels, core_panel = load_fixed_panels()
    write_tsv(panels, "validation_panel_definitions.tsv")

    chen_scores, chen_metrics, chen_correlations = chen_validation(panels, core_panel)
    external_scores, external_tests, external_pooled, external_adjusted, external_loo = external_validation(panels)
    ffpe_scores, ffpe_tests, ffpe_gene_tests = ffpe_validation(panels)
    perturbation_scores, apc_effects, tcf_effects = perturbation_validation(panels)

    qa = qa_checks(panels, chen_metrics, external_pooled, ffpe_tests)

    write_tsv(chen_scores, "validation_chen_panel_scores.tsv.gz", compress=True)
    write_tsv(chen_metrics, "validation_chen_panel_metrics.tsv")
    write_tsv(chen_correlations, "validation_chen_core_fidelity.tsv")
    write_tsv(external_scores, "validation_external_panel_scores.tsv.gz", compress=True)
    write_tsv(external_tests, "validation_external_cohort_tests.tsv")
    write_tsv(external_pooled, "validation_external_pooled_models.tsv")
    write_tsv(external_adjusted, "validation_external_proliferation_adjusted_models.tsv")
    write_tsv(external_loo, "validation_external_adjusted_leave_one_cohort_out.tsv")
    write_tsv(ffpe_scores, "validation_ffpe_panel_scores.tsv.gz", compress=True)
    write_tsv(ffpe_tests, "validation_ffpe_paired_tests.tsv")
    write_tsv(ffpe_gene_tests, "validation_ffpe_gene_tests.tsv")
    write_tsv(perturbation_scores, "validation_perturbation_sample_scores.tsv")
    write_tsv(apc_effects, "validation_apc_organoid_effects.tsv")
    write_tsv(tcf_effects, "validation_tcf7l2_clone_effects.tsv")
    write_tsv(qa, "validation_pipeline_qa.tsv")

    manifest = {
        "analysis": "validation of the frozen objective compact panel",
        "objective_panel_sha256_before_validation_read": frozen_checksum,
        "panel_frozen_before_validation": True,
        "gene_reselection": False,
        "weight_fitting": False,
        "cutpoint_optimization": False,
        "validation_layers": [
            "Chen held-out",
            "five independent sporadic-adenoma cohorts",
            "GSE117606 paired FFPE",
            "APC-knockout organoids",
            "TCF7L2-knockout CRC clones",
        ],
        "previous_10_gene_panel_role": "internal comparator and reproducibility control only",
        "one_hundred_gene_layer_included": False,
    }
    (OUT_DIR / "validation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    if not qa["passed"].all():
        failed = qa.loc[~qa["passed"], "check"].tolist()
        raise RuntimeError("Validation QA failed: " + ", ".join(failed))


if __name__ == "__main__":
    main()
