#!/usr/bin/env python3
"""Validate the frozen 11-gene strict-majority stability consensus.

This is a method-refinement sensitivity analysis.  The gene set is derived
solely from donor-held-out discovery paths and is checksummed before validation
loaders run.  No validation-driven gene replacement or weighting is allowed.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
OUT = ROOT / "results" / "stability_consensus_panel_v2_8"
PANEL_PATH = OUT / "stability_consensus_panel_frozen.tsv"
CORE_PATH = (
    ROOT
    / "results"
    / "data_adaptive_panel_pilot_v2_6"
    / "stable_error_controlled_core.tsv"
)
sys.path.insert(0, str(ANALYSIS))

import validate_objective_compact_panel_v2_7 as validation  # noqa: E402


PANEL_ID = "stability_consensus_11"
PANEL_LABEL = "Strict-majority stability consensus (11)"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write(frame: pd.DataFrame, filename: str, compress: bool = False) -> None:
    frame.to_csv(
        OUT / filename,
        sep="\t",
        index=False,
        compression="gzip" if compress else None,
    )


def fixed_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    panel = pd.read_csv(PANEL_PATH, sep="\t")
    if len(panel) != 11:
        raise RuntimeError(f"Expected 11 consensus genes, found {len(panel)}")
    if panel["validation_outcomes_used"].astype(str).str.lower().ne("false").any():
        raise RuntimeError("Consensus panel records validation-outcome use")
    panel = panel.copy()
    panel["panel_id"] = PANEL_ID
    panel["panel_label"] = PANEL_LABEL
    core = pd.read_csv(CORE_PATH, sep="\t")[["gene", "arm"]]
    core["route_weight"] = np.where(core["arm"].eq("up"), 1.0, -1.0)
    return panel, core


def qa_checks(
    panel: pd.DataFrame,
    chen: pd.DataFrame,
    external_adjusted: pd.DataFrame,
    ffpe: pd.DataFrame,
    apc: pd.DataFrame,
    tcf: pd.DataFrame,
) -> pd.DataFrame:
    heldout = chen.loc[chen["dataset"].eq("validation")].iloc[0]
    external = external_adjusted.iloc[0]
    ffpe_row = ffpe.iloc[0]
    apc_non_null = apc.loc[apc["expected_direction"].ne(0)]
    tcf_ko = tcf.loc[tcf["genotype"].eq("KO")]
    rows = [
        {
            "check": "strict_majority_panel_has_data_derived_size",
            "passed": len(panel) == 11 and panel["fixed_gene_count_used"].astype(str).str.lower().eq("false").all(),
            "observed": len(panel),
        },
        {
            "check": "both_arms_retained",
            "passed": set(panel["arm"]) == {"up", "down"},
            "observed": panel["arm"].value_counts().to_dict(),
        },
        {
            "check": "heldout_auc_at_least_0_90",
            "passed": float(heldout["auc"]) >= 0.90,
            "observed": float(heldout["auc"]),
        },
        {
            "check": "external_adjusted_ci_excludes_zero",
            "passed": float(external["ci_low"]) > 0,
            "observed": float(external["ci_low"]),
        },
        {
            "check": "ffpe_positive_fraction_at_least_0_90",
            "passed": float(ffpe_row["paired_positive_fraction"]) >= 0.90,
            "observed": float(ffpe_row["paired_positive_fraction"]),
        },
        {
            "check": "apc_non_null_directions_supported",
            "passed": bool(apc_non_null["direction_matches_expected"].fillna(False).all()),
            "observed": f"{int(apc_non_null['direction_matches_expected'].sum())}/{len(apc_non_null)}",
        },
        {
            "check": "tcf7l2_ko_directions_supported",
            "passed": bool(tcf_ko["direction_matches_expected"].fillna(False).all()),
            "observed": f"{int(tcf_ko['direction_matches_expected'].sum())}/{len(tcf_ko)}",
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    checksum_before_validation = sha256(PANEL_PATH)
    panel, core = fixed_inputs()
    validation.PANEL_ORDER = [PANEL_ID]
    validation.PANEL_LABELS = {PANEL_ID: PANEL_LABEL}

    chen_scores, chen_metrics, chen_fidelity = validation.chen_validation(panel, core)
    external_scores, external_tests, external_pooled, external_adjusted, external_loo = (
        validation.external_validation(panel)
    )
    ffpe_scores, ffpe_tests, ffpe_genes = validation.ffpe_validation(panel)
    perturb_scores, apc, tcf = validation.perturbation_validation(panel)
    qa = qa_checks(panel, chen_metrics, external_adjusted, ffpe_tests, apc, tcf)

    write(panel, "validation_panel_definition.tsv")
    write(chen_scores, "validation_chen_scores.tsv.gz", compress=True)
    write(chen_metrics, "validation_chen_metrics.tsv")
    write(chen_fidelity, "validation_chen_fidelity_to_287.tsv")
    write(external_scores, "validation_external_scores.tsv.gz", compress=True)
    write(external_tests, "validation_external_cohort_tests.tsv")
    write(external_pooled, "validation_external_pooled_model.tsv")
    write(external_adjusted, "validation_external_adjusted_model.tsv")
    write(external_loo, "validation_external_adjusted_leave_one_cohort_out.tsv")
    write(ffpe_scores, "validation_ffpe_scores.tsv.gz", compress=True)
    write(ffpe_tests, "validation_ffpe_paired_test.tsv")
    write(ffpe_genes, "validation_ffpe_gene_tests.tsv")
    write(perturb_scores, "validation_perturbation_scores.tsv")
    write(apc, "validation_apc_organoid_effects.tsv")
    write(tcf, "validation_tcf7l2_clone_effects.tsv")
    write(qa, "validation_qa.tsv")

    manifest = {
        "analysis": "validation of frozen strict-majority stability consensus",
        "panel_sha256_before_validation_read": checksum_before_validation,
        "panel_frozen_before_validation": True,
        "gene_count_fixed": False,
        "gene_reselection_after_validation": False,
        "weight_fitting": False,
        "cutpoint_optimization": False,
        "method_refinement_status": "sensitivity analysis prompted by the negative random benchmark of the initial Kneedle-path panel",
        "all_primary_qa_checks_passed": bool(qa["passed"].all()),
    }
    (OUT / "validation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
