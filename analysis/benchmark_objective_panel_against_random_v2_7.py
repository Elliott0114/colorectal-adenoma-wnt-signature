#!/usr/bin/env python3
"""Benchmark the frozen 12-gene panel against balanced random panels.

Random panels contain six up-arm and six down-arm genes drawn without
replacement from the same 62-gene portable, protein-coding candidate universe.
The benchmark is validation-only and cannot alter the frozen panel.
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
ANALYSIS = ROOT / "analysis"
OUT = ROOT / "results" / "objective_compact_panel_v2_7" / "random_benchmark"
PANEL_PATH = (
    ROOT
    / "results"
    / "objective_compact_panel_v2_7"
    / "objective_compact_panel_frozen.tsv"
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
SEED = 20260810
N_RANDOM = 10_000
sys.path.insert(0, str(ANALYSIS))

import conventional_route_signature_transfer as transfer  # noqa: E402
import discovery_locked_route_signature as discovery  # noqa: E402
import translation_reduced_panel_v2_0 as reduced  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def zscore_expression(expression: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    frame = expression.loc[:, genes].astype(float)
    return frame.sub(frame.mean(axis=0), axis=1).div(
        frame.std(axis=0, ddof=1).replace(0, np.nan), axis=1
    )


def balanced_panels(
    candidate: pd.DataFrame, n_up: int = 6, n_down: int = 6
) -> tuple[np.ndarray, list[str], list[str]]:
    up = sorted(candidate.loc[candidate["arm"].eq("up"), "gene"].tolist())
    down = sorted(candidate.loc[candidate["arm"].eq("down"), "gene"].tolist())
    rng = np.random.default_rng(SEED)
    seen: set[tuple[tuple[str, ...], tuple[str, ...]]] = set()
    panels: list[tuple[tuple[str, ...], tuple[str, ...]]] = []
    while len(panels) < N_RANDOM:
        selected_up = tuple(sorted(rng.choice(up, size=n_up, replace=False).tolist()))
        selected_down = tuple(sorted(rng.choice(down, size=n_down, replace=False).tolist()))
        key = (selected_up, selected_down)
        if key not in seen:
            seen.add(key)
            panels.append(key)
    genes = up + down
    gene_pos = {gene: pos for pos, gene in enumerate(genes)}
    weights = np.zeros((len(genes), N_RANDOM), dtype=np.float64)
    for column, (selected_up, selected_down) in enumerate(panels):
        weights[[gene_pos[gene] for gene in selected_up], column] = 1 / n_up
        weights[[gene_pos[gene] for gene in selected_down], column] = -1 / n_down
    panel_rows = []
    for panel_id, (selected_up, selected_down) in enumerate(panels, start=1):
        panel_rows.append(
            {
                "random_panel_id": panel_id,
                "up_genes": ",".join(selected_up),
                "down_genes": ",".join(selected_down),
            }
        )
    pd.DataFrame(panel_rows).to_csv(
        OUT / "random_panel_gene_sets.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    return weights, genes, up


def objective_weights(genes: list[str], panel: pd.DataFrame) -> np.ndarray:
    weights = np.zeros((len(genes), 1), dtype=float)
    gene_pos = {gene: pos for pos, gene in enumerate(genes)}
    arm_sizes = panel.groupby("arm")["gene"].size().to_dict()
    for row in panel.itertuples(index=False):
        weights[gene_pos[row.gene], 0] = float(row.route_weight) / arm_sizes[row.arm]
    return weights


def auc_vector(scores: np.ndarray, positive: np.ndarray, negative: np.ndarray) -> np.ndarray:
    pos = scores[positive, :]
    neg = scores[negative, :]
    pairwise = pos[:, None, :] - neg[None, :, :]
    return ((pairwise > 0).sum(axis=(0, 1)) + 0.5 * (pairwise == 0).sum(axis=(0, 1))) / (
        len(positive) * len(negative)
    )


def spearman_vector(scores: np.ndarray, target: np.ndarray) -> np.ndarray:
    target_rank = stats.rankdata(target).astype(float)
    target_rank = (target_rank - target_rank.mean()) / target_rank.std(ddof=0)
    score_rank = np.apply_along_axis(stats.rankdata, 0, scores).astype(float)
    score_rank = score_rank - score_rank.mean(axis=0, keepdims=True)
    score_rank = score_rank / score_rank.std(axis=0, ddof=0, keepdims=True)
    return np.mean(score_rank * target_rank[:, None], axis=0)


def chen_metrics(
    all_weights: np.ndarray, genes: list[str], core: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray]:
    meta, expression = discovery.load_chen_pseudobulk("validation")
    z = zscore_expression(expression, genes)
    scores = z.to_numpy(float) @ all_weights
    positive = np.flatnonzero(meta["route_group"].eq("conventional_adenoma").to_numpy())
    negative = np.flatnonzero(meta["route_group"].eq("normal").to_numpy())
    auc = auc_vector(scores, positive, negative)
    core_score = transfer.route_score_from_expression(
        expression,
        core.loc[core["arm"].eq("up"), "gene"].tolist(),
        core.loc[core["arm"].eq("down"), "gene"].tolist(),
    )["score__ca_route_signature"].to_numpy(float)
    fidelity = spearman_vector(scores, core_score)
    return auc, fidelity


def external_metrics(
    all_weights: np.ndarray, genes: list[str], candidate: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray]:
    external_data = reduced.external_gene_data(candidate)
    cohort_aucs = []
    weights = []
    for cohort, (meta, expression) in external_data.items():
        z = zscore_expression(expression, genes).reindex(meta["sample_id"])
        scores = z.to_numpy(float) @ all_weights
        positive = np.flatnonzero(meta["tissue_group"].eq("adenoma").to_numpy())
        negative = np.flatnonzero(meta["tissue_group"].eq("normal").to_numpy())
        cohort_aucs.append(auc_vector(scores, positive, negative))
        weights.append(len(positive) * len(negative))
    auc_matrix = np.vstack(cohort_aucs)
    weighted_auc = np.average(auc_matrix, axis=0, weights=np.asarray(weights))
    minimum_auc = np.min(auc_matrix, axis=0)
    return weighted_auc, minimum_auc


def ffpe_metrics(
    all_weights: np.ndarray, genes: list[str], candidate: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray]:
    scores, _, _, expression = reduced.ffpe_analysis(candidate)
    meta = scores[["sample_id", "patient_id", "tissue_group"]].copy()
    z = zscore_expression(expression, genes).reindex(meta["sample_id"])
    panel_scores = z.to_numpy(float) @ all_weights
    sample_pos = {sample: pos for pos, sample in enumerate(meta["sample_id"])}
    paired = (
        meta.loc[meta["tissue_group"].isin(["normal", "adenoma"])]
        .pivot(index="patient_id", columns="tissue_group", values="sample_id")
        .dropna()
    )
    adenoma_rows = np.array([sample_pos[sample] for sample in paired["adenoma"]])
    normal_rows = np.array([sample_pos[sample] for sample in paired["normal"]])
    delta = panel_scores[adenoma_rows, :] - panel_scores[normal_rows, :]
    return np.median(delta, axis=0), np.mean(delta > 0, axis=0)


def empirical_p(random: np.ndarray, observed: float) -> float:
    return float((1 + np.sum(random >= observed)) / (len(random) + 1))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checksum_before_validation = sha256(PANEL_PATH)
    panel = pd.read_csv(PANEL_PATH, sep="\t")
    candidate = pd.read_csv(CANDIDATE_PATH, sep="\t")
    candidate["route_weight"] = np.where(candidate["arm"].eq("up"), 1, -1)
    core = pd.read_csv(CORE_PATH, sep="\t")
    if len(panel) != 12 or len(candidate) != 62 or len(core) != 287:
        raise RuntimeError("Unexpected frozen panel, candidate universe, or core size")

    random_weights, genes, _ = balanced_panels(candidate)
    observed_weights = objective_weights(genes, panel)
    all_weights = np.column_stack([observed_weights, random_weights])

    chen_auc, chen_fidelity = chen_metrics(all_weights, genes, core)
    external_auc, external_min_auc = external_metrics(all_weights, genes, candidate)
    ffpe_delta, ffpe_positive = ffpe_metrics(all_weights, genes, candidate)

    metric_arrays = {
        "chen_heldout_auc": chen_auc,
        "chen_heldout_fidelity_to_287": chen_fidelity,
        "external_pair_weighted_auc": external_auc,
        "external_minimum_cohort_auc": external_min_auc,
        "ffpe_median_paired_delta": ffpe_delta,
        "ffpe_positive_pair_fraction": ffpe_positive,
    }
    rows = []
    random_table = pd.DataFrame({"random_panel_id": np.arange(1, N_RANDOM + 1)})
    for metric, values in metric_arrays.items():
        observed = float(values[0])
        random = values[1:]
        random_table[metric] = random
        rows.append(
            {
                "metric": metric,
                "objective_12_observed": observed,
                "random_median": float(np.median(random)),
                "random_q025": float(np.quantile(random, 0.025)),
                "random_q975": float(np.quantile(random, 0.975)),
                "objective_percentile_among_random": float(np.mean(random <= observed)),
                "empirical_p_random_at_least_observed": empirical_p(random, observed),
            }
        )

    joint = (
        (chen_auc[1:] >= chen_auc[0])
        & (chen_fidelity[1:] >= chen_fidelity[0])
        & (external_auc[1:] >= external_auc[0])
        & (ffpe_positive[1:] >= ffpe_positive[0])
    )
    joint_count = int(joint.sum())
    joint_p = float((joint_count + 1) / (N_RANDOM + 1))
    rows.append(
        {
            "metric": "joint_all_four_primary_metrics",
            "objective_12_observed": 1.0,
            "random_median": np.nan,
            "random_q025": np.nan,
            "random_q975": np.nan,
            "objective_percentile_among_random": 1 - joint_count / N_RANDOM,
            "empirical_p_random_at_least_observed": joint_p,
        }
    )
    random_table.to_csv(
        OUT / "random_panel_validation_metrics.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "random_panel_benchmark_summary.tsv", sep="\t", index=False)

    qa = pd.DataFrame(
        [
            {
                "check": "objective_panel_checksum_recorded_before_validation",
                "passed": bool(checksum_before_validation),
                "observed": checksum_before_validation,
            },
            {
                "check": "ten_thousand_unique_balanced_random_panels",
                "passed": len(random_table) == N_RANDOM,
                "observed": len(random_table),
            },
            {
                "check": "objective_exceeds_random_median_on_all_primary_metrics",
                "passed": all(
                    metric_arrays[key][0] > np.median(metric_arrays[key][1:])
                    for key in [
                        "chen_heldout_auc",
                        "chen_heldout_fidelity_to_287",
                        "external_pair_weighted_auc",
                        "ffpe_positive_pair_fraction",
                    ]
                ),
                "observed": "four primary metrics",
            },
        ]
    )
    qa.to_csv(OUT / "random_panel_benchmark_qa.tsv", sep="\t", index=False)
    manifest = {
        "analysis": "balanced random-panel validation benchmark",
        "panel_sha256_before_validation_read": checksum_before_validation,
        "n_random_panels": N_RANDOM,
        "random_seed": SEED,
        "random_panel_design": "6 up + 6 down without replacement from the fixed 30-up/32-down portable protein-coding universe",
        "panel_reselection": False,
        "joint_random_panels_meeting_or_exceeding_objective": joint_count,
        "joint_empirical_p": joint_p,
    }
    (OUT / "random_panel_benchmark_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
