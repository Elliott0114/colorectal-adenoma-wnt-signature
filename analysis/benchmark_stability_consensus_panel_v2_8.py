#!/usr/bin/env python3
"""Random-panel benchmark for the frozen 6-up/5-down stability consensus."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
OUT = ROOT / "results" / "stability_consensus_panel_v2_8" / "random_benchmark"
PANEL_PATH = (
    ROOT
    / "results"
    / "stability_consensus_panel_v2_8"
    / "stability_consensus_panel_frozen.tsv"
)
CANDIDATE_PATH = (
    ROOT
    / "results"
    / "objective_compact_panel_v2_7"
    / "portable_protein_coding_candidate_universe.tsv"
)
CORE_PATH = (
    ROOT
    / "results"
    / "data_adaptive_panel_pilot_v2_6"
    / "stable_error_controlled_core.tsv"
)
sys.path.insert(0, str(ANALYSIS))

import benchmark_objective_panel_against_random_v2_7 as benchmark  # noqa: E402


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checksum = benchmark.sha256(PANEL_PATH)
    panel = pd.read_csv(PANEL_PATH, sep="\t")
    candidate = pd.read_csv(CANDIDATE_PATH, sep="\t")
    candidate["route_weight"] = np.where(candidate["arm"].eq("up"), 1, -1)
    core = pd.read_csv(CORE_PATH, sep="\t")
    if panel["arm"].value_counts().to_dict() != {"up": 6, "down": 5}:
        raise RuntimeError("Expected a frozen six-up/five-down consensus panel")

    benchmark.OUT = OUT
    random_weights, genes, _ = benchmark.balanced_panels(
        candidate, n_up=6, n_down=5
    )
    observed_weights = benchmark.objective_weights(genes, panel)
    all_weights = np.column_stack([observed_weights, random_weights])

    chen_auc, chen_fidelity = benchmark.chen_metrics(all_weights, genes, core)
    external_auc, external_min_auc = benchmark.external_metrics(
        all_weights, genes, candidate
    )
    ffpe_delta, ffpe_positive = benchmark.ffpe_metrics(
        all_weights, genes, candidate
    )
    arrays = {
        "chen_heldout_auc": chen_auc,
        "chen_heldout_fidelity_to_287": chen_fidelity,
        "external_pair_weighted_auc": external_auc,
        "external_minimum_cohort_auc": external_min_auc,
        "ffpe_median_paired_delta": ffpe_delta,
        "ffpe_positive_pair_fraction": ffpe_positive,
    }
    random_table = pd.DataFrame(
        {"random_panel_id": np.arange(1, benchmark.N_RANDOM + 1)}
    )
    rows = []
    for metric, values in arrays.items():
        observed = float(values[0])
        random = values[1:]
        random_table[metric] = random
        rows.append(
            {
                "metric": metric,
                "consensus_11_observed": observed,
                "random_median": float(np.median(random)),
                "random_q025": float(np.quantile(random, 0.025)),
                "random_q975": float(np.quantile(random, 0.975)),
                "consensus_percentile_among_random": float(np.mean(random <= observed)),
                "empirical_p_random_at_least_observed": benchmark.empirical_p(
                    random, observed
                ),
            }
        )
    joint = (
        (chen_auc[1:] >= chen_auc[0])
        & (chen_fidelity[1:] >= chen_fidelity[0])
        & (external_auc[1:] >= external_auc[0])
        & (ffpe_positive[1:] >= ffpe_positive[0])
    )
    joint_count = int(joint.sum())
    joint_p = float((joint_count + 1) / (benchmark.N_RANDOM + 1))
    rows.append(
        {
            "metric": "joint_all_four_primary_metrics",
            "consensus_11_observed": 1.0,
            "random_median": np.nan,
            "random_q025": np.nan,
            "random_q975": np.nan,
            "consensus_percentile_among_random": 1 - joint_count / benchmark.N_RANDOM,
            "empirical_p_random_at_least_observed": joint_p,
        }
    )
    random_table.to_csv(
        OUT / "random_panel_validation_metrics.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(rows).to_csv(
        OUT / "random_panel_benchmark_summary.tsv", sep="\t", index=False
    )
    manifest = {
        "analysis": "random benchmark of strict-majority stability consensus",
        "panel_sha256_before_validation_read": checksum,
        "n_random_panels": benchmark.N_RANDOM,
        "random_seed": benchmark.SEED,
        "random_panel_design": "6 up + 5 down without replacement from the same fixed portable candidate universe",
        "joint_random_panels_meeting_or_exceeding_consensus": joint_count,
        "joint_empirical_p": joint_p,
        "gene_reselection": False,
    }
    (OUT / "random_panel_benchmark_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
