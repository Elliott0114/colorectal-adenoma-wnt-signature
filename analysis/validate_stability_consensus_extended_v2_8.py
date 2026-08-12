#!/usr/bin/env python3
"""Run all extended layers for the frozen 11-gene stability consensus."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
OUT = ROOT / "results" / "stability_consensus_panel_v2_8" / "extended_validation"
PANEL_PATH = (
    ROOT
    / "results"
    / "stability_consensus_panel_v2_8"
    / "stability_consensus_panel_frozen.tsv"
)
sys.path.insert(0, str(ANALYSIS))

import validate_objective_panel_extended_layers_v2_7 as extended  # noqa: E402


def write(frame: pd.DataFrame, filename: str) -> None:
    frame.to_csv(OUT / filename, sep="\t", index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checksum = extended.sha256(PANEL_PATH)
    panel = pd.read_csv(PANEL_PATH, sep="\t")
    if len(panel) != 11:
        raise RuntimeError("Expected the frozen 11-gene stability consensus")
    panel["signature_direction"] = np.where(
        panel["arm"].eq("up"), "adenoma_up", "adenoma_down"
    )
    extended.OUT = OUT

    becker = extended.becker_layers(panel)
    atlas = extended.atlas_layers(panel)
    perturbation_spatial = extended.perturbation_and_spatial_layers(panel)
    loo = extended.leave_one_gene_out(panel)
    qa = extended.qa_summary(becker, atlas, perturbation_spatial, loo)

    write(loo, "stability_consensus_leave_one_gene_out.tsv")
    write(qa, "extended_validation_qa.tsv")
    write(
        pd.DataFrame(
            [
                {"layer": "becker", **becker},
                {"layer": "crc_atlas", **atlas},
                {"layer": "perturbation_spatial", **perturbation_spatial},
            ]
        ),
        "extended_validation_summary.tsv",
    )
    manifest = {
        "analysis": "extended validation of frozen strict-majority stability consensus",
        "panel_sha256_before_validation_read": checksum,
        "panel_frozen_before_validation": True,
        "gene_reselection": False,
        "weight_fitting": False,
        "cutpoint_optimization": False,
        "all_strict_qa_checks_passed": bool(qa["passed"].all()),
    }
    (OUT / "extended_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
