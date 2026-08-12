#!/usr/bin/env python3
"""Run frozen, validation-only GenKI virtual knockouts for the 12-gene panel.

The virtual knockout output is an unsigned network-impact ranking. The script
therefore tests only pre-existing gene-set coupling and coherence; it never
uses the ranking to select, replace, or reweight signature genes.
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import platform
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
import sklearn
import torch
import torch_geometric
from scipy.spatial.distance import cdist
from scipy.stats import hypergeom, mannwhitneyu, spearmanr
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Data

import GenKI
from GenKI import DataLoader, VGAE_trainer, build_adata, get_distance


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "virtual_knockout_validation_v2_9"
INPUT_DIR = OUT_DIR / "input"
MODEL_DIR = OUT_DIR / "models"
GRN_DIR = OUT_DIR / "grn_cache"
INPUT_H5AD = INPUT_DIR / "chen_validation_conventional_adenoma_balanced_raw.h5ad"
CONTRACT_PATH = INPUT_DIR / "frozen_validation_contract.json"

TCF_EFFECTS = (
    ROOT
    / "results"
    / "objective_compact_panel_v2_7"
    / "validation_tcf7l2_clone_effects.tsv"
)
APC_EFFECTS = (
    ROOT
    / "results"
    / "objective_compact_panel_v2_7"
    / "validation_apc_organoid_effects.tsv"
)
EXTRA_EFFECTS = (
    ROOT
    / "results"
    / "objective_compact_panel_v2_7"
    / "extended_validation"
    / "perturbation_spatial"
    / "perturbation_effect_summary.tsv"
)

EDGE_CUTOFF = 85
PC_COMPONENTS = 5
N_EPOCHS = 100
LATENT_DIMENSIONS = 2
LEARNING_RATE = 7e-4
WEIGHT_DECAY = 9e-4
BETA = 1e-4
MODEL_SEEDS = (20260810, 20260811)
DISTANCE_METRICS = ("KL", "EMD")
PRIMARY_DISTANCE = "KL"
MATCHED_NULL_REPLICATES = 10_000
MATCH_POOL_SIZE = 50
TOP_FRACTION = 0.05


class SeededVGAETrainer(VGAE_trainer):
    """Restore the requested model seed after GenKI's fixed link split.

    GenKI 0.2.1 calls ``seed_everything(42)`` while constructing its fixed
    train/validation/test split, after the requested trainer seed is set. That
    also resets model initialization. Restoring the requested seed here keeps
    the package's fixed split while making initialization-seed sensitivity
    real and auditable.
    """

    def _transform_data(self, *args, **kwargs):
        super()._transform_data(*args, **kwargs)
        if self.seed is not None:
            np.random.seed(self.seed)
            torch.manual_seed(self.seed)


def stable_seed(*parts: str) -> int:
    digest = hashlib.sha256("::".join(parts).encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little")


def bh_adjust(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    result = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return result
    order = valid.sort_values().index
    sorted_p = valid.loc[order].to_numpy(dtype=float)
    n = len(sorted_p)
    adjusted = sorted_p * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[order] = np.minimum(adjusted, 1.0)
    return result


def make_ko_data(wt_data: Data, target_index: int) -> Data:
    edge_index = wt_data.edge_index
    keep_edge = (edge_index[0] != target_index) & (edge_index[1] != target_index)
    x_ko = wt_data.x.clone()
    x_ko[target_index, :] = 0.0
    return Data(
        x=x_ko,
        edge_index=edge_index[:, keep_edge],
        y=wt_data.y,
    )


def impact_percentiles(distance: np.ndarray, target_index: int) -> np.ndarray:
    distance = np.asarray(distance, dtype=float)
    result = np.full(distance.shape, np.nan, dtype=float)
    valid = np.isfinite(distance)
    valid[target_index] = False
    values = distance[valid]
    if len(values) < 2:
        raise ValueError("Fewer than two finite non-target distances")
    ranks = scipy.stats.rankdata(-values, method="average")
    result[valid] = 1.0 - (ranks - 1.0) / (len(values) - 1.0)
    return result


def top_jaccard(a: np.ndarray, b: np.ndarray, fraction: float = TOP_FRACTION) -> float:
    valid = np.isfinite(a) & np.isfinite(b)
    n_top = max(1, math.ceil(valid.sum() * fraction))
    indices = np.flatnonzero(valid)
    top_a = set(indices[np.argsort(a[valid])[-n_top:]])
    top_b = set(indices[np.argsort(b[valid])[-n_top:]])
    return len(top_a & top_b) / len(top_a | top_b)


def matched_null_test(
    impacts: pd.Series,
    members: list[str],
    matching_features: pd.DataFrame,
    *,
    seed: int,
    n_replicates: int = MATCHED_NULL_REPLICATES,
) -> tuple[dict[str, float | int], np.ndarray]:
    members = [gene for gene in members if gene in impacts.index and pd.notna(impacts[gene])]
    universe = impacts.dropna().index
    member_set = set(members)
    candidates = [gene for gene in universe if gene not in member_set]
    if not members or len(candidates) < MATCH_POOL_SIZE:
        raise ValueError("Insufficient members or matched-null candidates")

    feature_frame = matching_features.loc[list(universe)].replace([np.inf, -np.inf], np.nan)
    if feature_frame.isna().any().any():
        raise ValueError("Non-finite matching features")
    scaler = StandardScaler().fit(feature_frame)
    scaled = pd.DataFrame(
        scaler.transform(feature_frame), index=feature_frame.index, columns=feature_frame.columns
    )
    distances = cdist(scaled.loc[members], scaled.loc[candidates], metric="euclidean")
    k = min(MATCH_POOL_SIZE, len(candidates))
    nearest_positions = np.argpartition(distances, kth=k - 1, axis=1)[:, :k]
    candidate_array = np.asarray(candidates, dtype=object)
    nearest_genes = candidate_array[nearest_positions]
    nearest_impact = impacts.loc[nearest_genes.ravel()].to_numpy().reshape(nearest_genes.shape)

    observed = float(impacts.loc[members].mean())
    rng = np.random.default_rng(seed)
    null = np.empty(n_replicates, dtype=np.float64)
    chunk_size = 500
    row_index = np.arange(len(members))[None, :]
    for start in range(0, n_replicates, chunk_size):
        stop = min(start + chunk_size, n_replicates)
        draws = rng.integers(0, k, size=(stop - start, len(members)))
        sampled = nearest_impact[row_index, draws]
        null[start:stop] = sampled.mean(axis=1)

    p_empirical = (1.0 + float(np.count_nonzero(null >= observed))) / (n_replicates + 1.0)
    null_sd = float(np.std(null, ddof=1))
    matched_z = (observed - float(np.mean(null))) / null_sd if null_sd > 0 else np.nan

    background = impacts.loc[candidates].to_numpy(dtype=float)
    member_values = impacts.loc[members].to_numpy(dtype=float)
    mw = mannwhitneyu(member_values, background, alternative="greater")
    auc = float(mw.statistic / (len(member_values) * len(background)))

    n_universe = len(universe)
    n_top = max(1, math.ceil(n_universe * TOP_FRACTION))
    top_genes = set(impacts.loc[universe].nlargest(n_top).index)
    n_top_members = len(member_set & top_genes)
    p_hypergeom = float(
        hypergeom.sf(n_top_members - 1, n_universe, len(members), n_top)
    )

    summary = {
        "n_universe": n_universe,
        "n_members": len(members),
        "observed_mean_impact_percentile": observed,
        "matched_null_mean": float(np.mean(null)),
        "matched_null_sd": null_sd,
        "matched_null_ci_low": float(np.quantile(null, 0.025)),
        "matched_null_ci_high": float(np.quantile(null, 0.975)),
        "matched_z": float(matched_z),
        "p_matched_empirical_one_sided": float(p_empirical),
        "auc_members_vs_background": auc,
        "p_mannwhitney_one_sided": float(mw.pvalue),
        "top_fraction": TOP_FRACTION,
        "n_top_genes": n_top,
        "n_members_in_top": n_top_members,
        "p_hypergeom_top_enrichment": p_hypergeom,
        "match_pool_size": k,
        "matched_null_replicates": n_replicates,
    }
    return summary, null


def aggregate_endpoint(
    tests: pd.DataFrame,
    null_store: dict[tuple[str, str, str], np.ndarray],
    *,
    distance_metric: str,
    endpoint: str,
    targets: list[str],
    gene_set: str,
) -> dict[str, object]:
    mask = (
        tests["distance_metric"].eq(distance_metric)
        & tests["target"].isin(targets)
        & tests["gene_set"].eq(gene_set)
    )
    frame = tests.loc[mask].copy()
    missing = sorted(set(targets) - set(frame["target"]))
    if missing:
        raise ValueError(f"Missing target-level tests for {endpoint}: {missing}")
    frame = frame.set_index("target").loc[targets].reset_index()
    null_matrix = np.vstack(
        [null_store[(distance_metric, target, gene_set)] for target in targets]
    )
    null = null_matrix.mean(axis=0)
    observed = float(frame["observed_mean_impact_percentile"].mean())
    p_empirical = (1.0 + float(np.count_nonzero(null >= observed))) / (len(null) + 1.0)
    return {
        "distance_metric": distance_metric,
        "endpoint": endpoint,
        "gene_set": gene_set,
        "targets": ",".join(targets),
        "n_targets": len(targets),
        "mean_gene_set_size": float(frame["n_members"].mean()),
        "observed_mean_impact_percentile": observed,
        "matched_null_mean": float(np.mean(null)),
        "matched_null_ci_low": float(np.quantile(null, 0.025)),
        "matched_null_ci_high": float(np.quantile(null, 0.975)),
        "matched_z": float((observed - np.mean(null)) / np.std(null, ddof=1)),
        "p_matched_empirical_one_sided": float(p_empirical),
        "target_level_positive_fraction": float(
            (frame["observed_mean_impact_percentile"] > frame["matched_null_mean"]).mean()
        ),
    }


def empirical_direction_calibration() -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    tcf = pd.read_csv(TCF_EFFECTS, sep="\t")
    tcf = tcf[tcf["panel_id"].eq("objective_12")].copy()
    rows.append(
        {
            "evidence": "TCF7L2 genetic perturbation",
            "dataset": "GSE135328",
            "comparison": "TCF7L2 KO/Het versus matched WT",
            "species_model": "human CRC cell lines",
            "n_units": int(len(tcf)),
            "mean_route_score_change": float(tcf["difference_vs_WT"].mean()),
            "median_route_score_change": float(tcf["difference_vs_WT"].median()),
            "expected_direction": "decrease",
            "direction_supported": bool((tcf["difference_vs_WT"] < 0).all()),
            "source_table": str(TCF_EFFECTS.relative_to(ROOT)),
        }
    )

    apc = pd.read_csv(APC_EFFECTS, sep="\t")
    apc = apc[(apc["panel_id"].eq("objective_12")) & (apc["feature"].eq("route_score"))]
    for comparison in ("APC_vs_WT_without_Wnt", "WT_withdrawal", "genotype_by_Wnt_interaction"):
        row = apc[apc["comparison"].eq(comparison)].iloc[0]
        rows.append(
            {
                "evidence": "APC/WNT organoid perturbation",
                "dataset": row["dataset"],
                "comparison": comparison,
                "species_model": "human isogenic colonic organoids",
                "n_units": int(row["n_units"]),
                "mean_route_score_change": float(row["mean_difference"]),
                "median_route_score_change": float(row["median_difference"]),
                "expected_direction": "increase" if row["expected_direction"] > 0 else "decrease",
                "direction_supported": bool(row["direction_matches_expected"]),
                "source_table": str(APC_EFFECTS.relative_to(ROOT)),
            }
        )

    extra = pd.read_csv(EXTRA_EFFECTS, sep="\t")
    extra = extra[extra["feature"].eq("route_score")]
    selected = {
        ("GSE130822", "ascl2_ko_vs_resting_wt"): "ASCL2 genetic perturbation",
        ("GSE171910", "conditional_wnt_silencing"): "conditional WNT silencing",
        ("GSE67186", "apc_restoration_shApc"): "APC restoration",
        ("GSE67186", "apc_restoration_shApc_Kras"): "APC restoration with KRAS",
    }
    for (dataset, comparison), evidence in selected.items():
        row = extra[(extra["dataset"].eq(dataset)) & (extra["comparison"].eq(comparison))].iloc[0]
        rows.append(
            {
                "evidence": evidence,
                "dataset": dataset,
                "comparison": comparison,
                "species_model": f"{row['species']} {row['model_system']}",
                "n_units": int(row["n_units"]),
                "mean_route_score_change": float(row["mean_difference"]),
                "median_route_score_change": float(row["median_difference"]),
                "expected_direction": "increase" if row["expected_direction"] > 0 else "decrease",
                "direction_supported": bool(
                    np.sign(row["mean_difference"]) == np.sign(row["expected_direction"])
                ),
                "source_table": str(EXTRA_EFFECTS.relative_to(ROOT)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    start_time = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    GRN_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )

    n_threads = min(16, os.cpu_count() or 1)
    torch.set_num_threads(n_threads)
    torch.set_num_interop_threads(min(4, n_threads))
    torch.use_deterministic_algorithms(True)

    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    panel_genes = [str(g).upper() for g in contract["frozen_sets"]["panel_genes"]]
    upstream_targets = [
        str(g).upper() for g in contract["prespecified_targets"]["upstream_context"]
    ]
    all_targets = [
        str(g).upper() for g in contract["prespecified_targets"]["all_unique"]
    ]

    raw = ad.read_h5ad(INPUT_H5AD)
    raw.var_names = raw.var_names.astype(str).str.upper()
    if raw.var_names.has_duplicates:
        raise ValueError("Prepared GenKI input has duplicate gene symbols")
    processed = build_adata(
        raw,
        log_normalize=True,
        scale_data=True,
        as_sparse=True,
        uppercase=True,
    )

    loader = DataLoader(
        processed,
        target_gene=[all_targets[0]],
        GRN_file_dir=str(GRN_DIR),
        rebuild_GRN=True,
        cutoff=EDGE_CUTOFF,
        svd_solver="shared",
        random_state=0,
        n_oversamples=20,
        timeit=True,
    )
    wt_data = loader.load_data()
    genes = np.asarray(wt_data.y, dtype=object)
    gene_to_index = {gene: i for i, gene in enumerate(genes)}
    missing_targets = sorted(set(all_targets) - set(gene_to_index))
    if missing_targets:
        raise ValueError(f"Targets absent from GenKI universe: {missing_targets}")

    edge_index = loader.edge_index.numpy()
    degree = np.bincount(edge_index[0], minlength=len(genes)).astype(int)
    # Sparse reductions can differ at the last floating-point bits when their
    # internal summation order changes. Weighted degree is reported only as a
    # diagnostic, so round it before export to make the release byte-stable.
    weighted_degree = np.round(
        np.asarray(np.abs(loader.net).sum(axis=1)).ravel(), decimals=10
    )
    gene_meta = processed.var.loc[genes].copy()
    gene_meta["gene"] = genes
    gene_meta["network_degree"] = degree
    gene_meta["network_weighted_degree"] = weighted_degree
    gene_meta["log1p_mean_raw_count"] = np.log1p(
        pd.to_numeric(gene_meta["mean_raw_count"], errors="coerce").to_numpy()
    )
    gene_meta["log1p_network_degree"] = np.log1p(degree)
    gene_meta.to_csv(OUT_DIR / "genki_network_gene_metrics.tsv", sep="\t", index=False)

    run_rows: list[dict[str, object]] = []
    for target in all_targets:
        idx = gene_to_index[target]
        run_rows.append(
            {
                "target": target,
                "target_role": (
                    "upstream_context_and_panel_member"
                    if target in upstream_targets and target in panel_genes
                    else "upstream_context"
                    if target in upstream_targets
                    else "panel_member"
                ),
                "panel_arm": gene_meta.loc[target, "panel_arm"],
                "detection_fraction": float(gene_meta.loc[target, "detection_fraction"]),
                "mean_raw_count": float(gene_meta.loc[target, "mean_raw_count"]),
                "network_degree": int(degree[idx]),
                "network_weighted_degree": float(weighted_degree[idx]),
            }
        )
    pd.DataFrame(run_rows).to_csv(
        OUT_DIR / "prespecified_knockout_targets.tsv", sep="\t", index=False
    )

    distance_rows: list[pd.DataFrame] = []
    training_rows: list[dict[str, object]] = []
    for seed in MODEL_SEEDS:
        logging.info("Training GenKI seed %s", seed)
        trainer = SeededVGAETrainer(
            wt_data,
            out_channels=LATENT_DIMENSIONS,
            epochs=N_EPOCHS,
            lr=LEARNING_RATE,
            weight_decay=WEIGHT_DECAY,
            beta=BETA,
            verbose=False,
            seed=seed,
        )
        trainer.train()
        epoch, loss, auc, ap = trainer.final_metrics
        training_rows.append(
            {
                "seed": seed,
                "epochs": epoch,
                "final_loss": float(loss),
                "test_edge_auc": float(auc),
                "test_edge_average_precision": float(ap),
            }
        )
        torch.save(
            {
                "state_dict": trainer.model.state_dict(),
                "seed": seed,
                "metrics": trainer.final_metrics,
                "genki_version": GenKI.__version__,
                "seed_fix": "requested seed restored after package fixed link split",
            },
            MODEL_DIR / f"genki_seed_{seed}.pt",
        )
        z_wt_mean, z_wt_var = trainer.get_latent_vars(wt_data)
        for target in all_targets:
            logging.info("Scoring virtual knockout seed=%s target=%s", seed, target)
            target_index = gene_to_index[target]
            ko_data = make_ko_data(wt_data, target_index)
            z_ko_mean, z_ko_var = trainer.get_latent_vars(ko_data)
            for distance_metric in DISTANCE_METRICS:
                distance = get_distance(
                    z_ko_mean,
                    z_ko_var,
                    z_wt_mean,
                    z_wt_var,
                    by=distance_metric,
                )
                impact = impact_percentiles(distance, target_index)
                distance_rows.append(
                    pd.DataFrame(
                        {
                            "seed": seed,
                            "target": target,
                            "distance_metric": distance_metric,
                            "gene": genes,
                            "distance": distance,
                            "impact_percentile": impact,
                            "is_knockout_target": genes == target,
                        }
                    )
                )
            del ko_data

    training = pd.DataFrame(training_rows)
    training.to_csv(OUT_DIR / "genki_training_metrics.tsv", sep="\t", index=False)
    distances = pd.concat(distance_rows, ignore_index=True)
    distances.to_csv(
        OUT_DIR / "genki_gene_impact_by_seed.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    consensus = (
        distances.groupby(["target", "distance_metric", "gene"], observed=True)
        .agg(
            mean_distance=("distance", "mean"),
            mean_impact_percentile=("impact_percentile", "mean"),
            min_impact_percentile=("impact_percentile", "min"),
            max_impact_percentile=("impact_percentile", "max"),
            is_knockout_target=("is_knockout_target", "first"),
        )
        .reset_index()
    )
    consensus.to_csv(
        OUT_DIR / "genki_gene_impact_consensus.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    stability_rows = []
    for target in all_targets:
        for metric in DISTANCE_METRICS:
            frame = distances[
                distances["target"].eq(target) & distances["distance_metric"].eq(metric)
            ].pivot(index="gene", columns="seed", values="impact_percentile")
            a = frame[MODEL_SEEDS[0]].to_numpy(dtype=float)
            b = frame[MODEL_SEEDS[1]].to_numpy(dtype=float)
            valid = np.isfinite(a) & np.isfinite(b)
            rho, p_value = spearmanr(a[valid], b[valid])
            stability_rows.append(
                {
                    "target": target,
                    "distance_metric": metric,
                    "n_ranked_genes": int(valid.sum()),
                    "spearman_rho_between_seeds": float(rho),
                    "spearman_p_value": float(p_value),
                    "top5pct_jaccard_between_seeds": top_jaccard(a, b),
                }
            )
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(OUT_DIR / "genki_seed_stability.tsv", sep="\t", index=False)

    metric_concordance_rows = []
    for target in all_targets:
        frame = consensus[consensus["target"].eq(target)].pivot(
            index="gene", columns="distance_metric", values="mean_impact_percentile"
        )
        valid = frame["KL"].notna() & frame["EMD"].notna()
        rho, p_value = spearmanr(frame.loc[valid, "KL"], frame.loc[valid, "EMD"])
        metric_concordance_rows.append(
            {
                "target": target,
                "n_ranked_genes": int(valid.sum()),
                "spearman_rho_kl_vs_emd": float(rho),
                "spearman_p_value": float(p_value),
            }
        )
    pd.DataFrame(metric_concordance_rows).to_csv(
        OUT_DIR / "genki_distance_metric_sensitivity.tsv", sep="\t", index=False
    )

    core_genes = gene_meta.index[gene_meta["fixed_287_core"].astype(bool)].tolist()
    matching_features = gene_meta.set_index("gene")[[
        "log1p_mean_raw_count",
        "detection_fraction",
        "log1p_network_degree",
    ]]
    test_rows: list[dict[str, object]] = []
    null_store: dict[tuple[str, str, str], np.ndarray] = {}
    for metric in DISTANCE_METRICS:
        for target in all_targets:
            impact = (
                consensus[
                    consensus["target"].eq(target)
                    & consensus["distance_metric"].eq(metric)
                ]
                .set_index("gene")["mean_impact_percentile"]
                .astype(float)
            )
            gene_sets = {
                "measurable_fixed_287_core": [g for g in core_genes if g != target],
                "leave_target_out_fixed_12_panel": [
                    g for g in panel_genes if g != target
                ],
            }
            for gene_set, members in gene_sets.items():
                summary, null = matched_null_test(
                    impact,
                    members,
                    matching_features,
                    seed=stable_seed(metric, target, gene_set),
                )
                summary.update(
                    {
                        "distance_metric": metric,
                        "target": target,
                        "target_role": (
                            "upstream_context_and_panel_member"
                            if target in upstream_targets and target in panel_genes
                            else "upstream_context"
                            if target in upstream_targets
                            else "panel_member"
                        ),
                        "gene_set": gene_set,
                    }
                )
                test_rows.append(summary)
                null_store[(metric, target, gene_set)] = null

    tests = pd.DataFrame(test_rows)
    tests["q_matched_all_target_level"] = tests.groupby("distance_metric", observed=True)[
        "p_matched_empirical_one_sided"
    ].transform(bh_adjust)
    upstream_family = (
        tests["distance_metric"].eq(PRIMARY_DISTANCE)
        & tests["target"].isin(upstream_targets)
    )
    tests.loc[upstream_family, "q_matched_primary_upstream_family"] = bh_adjust(
        tests.loc[upstream_family, "p_matched_empirical_one_sided"]
    )
    panel_family = (
        tests["distance_metric"].eq(PRIMARY_DISTANCE)
        & tests["target"].isin(panel_genes)
        & tests["gene_set"].eq("leave_target_out_fixed_12_panel")
    )
    tests.loc[panel_family, "q_matched_panel_connectivity_family"] = bh_adjust(
        tests.loc[panel_family, "p_matched_empirical_one_sided"]
    )
    tests.to_csv(OUT_DIR / "genki_fixed_gene_set_tests.tsv", sep="\t", index=False)

    aggregate_rows = []
    for metric in DISTANCE_METRICS:
        aggregate_rows.extend(
            [
                aggregate_endpoint(
                    tests,
                    null_store,
                    distance_metric=metric,
                    endpoint="upstream_context_to_core",
                    targets=upstream_targets,
                    gene_set="measurable_fixed_287_core",
                ),
                aggregate_endpoint(
                    tests,
                    null_store,
                    distance_metric=metric,
                    endpoint="upstream_context_to_panel",
                    targets=upstream_targets,
                    gene_set="leave_target_out_fixed_12_panel",
                ),
                aggregate_endpoint(
                    tests,
                    null_store,
                    distance_metric=metric,
                    endpoint="within_panel_virtual_knockout_coherence",
                    targets=panel_genes,
                    gene_set="leave_target_out_fixed_12_panel",
                ),
            ]
        )
    aggregate = pd.DataFrame(aggregate_rows)
    primary_aggregate = aggregate["distance_metric"].eq(PRIMARY_DISTANCE)
    aggregate.loc[primary_aggregate, "q_matched_primary_endpoints"] = bh_adjust(
        aggregate.loc[primary_aggregate, "p_matched_empirical_one_sided"]
    )
    aggregate.to_csv(OUT_DIR / "genki_aggregate_validation_endpoints.tsv", sep="\t", index=False)

    panel_matrix = consensus[
        consensus["distance_metric"].eq(PRIMARY_DISTANCE)
        & consensus["target"].isin(panel_genes)
        & consensus["gene"].isin(panel_genes)
    ].pivot(index="target", columns="gene", values="mean_impact_percentile")
    panel_matrix = panel_matrix.reindex(index=panel_genes, columns=panel_genes)
    panel_matrix.to_csv(OUT_DIR / "genki_panel_knockout_impact_matrix.tsv", sep="\t")

    calibration = empirical_direction_calibration()
    calibration.to_csv(
        OUT_DIR / "empirical_direction_calibration.tsv", sep="\t", index=False
    )

    qa_rows = [
        {
            "check": "frozen_contract_validation_only",
            "passed": contract["analysis_role"] == "validation_not_discovery",
            "detail": contract["analysis_role"],
        },
        {
            "check": "all_prespecified_targets_run",
            "passed": set(all_targets) == set(consensus["target"]),
            "detail": f"{consensus['target'].nunique()}/{len(all_targets)}",
        },
        {
            "check": "two_distinct_model_seeds",
            "passed": distances["seed"].nunique() == 2,
            "detail": ",".join(map(str, sorted(distances["seed"].unique()))),
        },
        {
            "check": "all_panel_genes_measurable",
            "passed": set(panel_genes).issubset(gene_meta.index),
            "detail": f"{len(set(panel_genes) & set(gene_meta.index))}/{len(panel_genes)}",
        },
        {
            "check": "finite_non_target_primary_impacts",
            "passed": bool(
                distances.loc[
                    distances["distance_metric"].eq(PRIMARY_DISTANCE)
                    & ~distances["is_knockout_target"],
                    "impact_percentile",
                ].notna().all()
            ),
            "detail": "KL impact percentiles",
        },
        {
            "check": "matched_null_replicates_complete",
            "passed": all(len(v) == MATCHED_NULL_REPLICATES for v in null_store.values()),
            "detail": f"{len(null_store)} tests x {MATCHED_NULL_REPLICATES}",
        },
        {
            "check": "empirical_direction_calibration_not_genki_signed",
            "passed": bool(calibration["direction_supported"].all()),
            "detail": f"{int(calibration['direction_supported'].sum())}/{len(calibration)} existing contrasts",
        },
    ]
    qa = pd.DataFrame(qa_rows)
    qa.to_csv(OUT_DIR / "genki_validation_qa.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "frozen 12-gene panel virtual-knockout validation",
        "analysis_role": "validation_not_discovery",
        "input_h5ad": str(INPUT_H5AD.relative_to(ROOT)),
        "n_cells": int(raw.n_obs),
        "n_donors": int(raw.obs["donor_id"].nunique()),
        "n_genes": int(raw.n_vars),
        "n_measurable_core_genes": len(core_genes),
        "panel_genes": panel_genes,
        "upstream_targets": upstream_targets,
        "all_knockout_targets": all_targets,
        "model_seeds": list(MODEL_SEEDS),
        "runtime_seed_fix": True,
        "epochs": N_EPOCHS,
        "edge_cutoff_percentile": EDGE_CUTOFF,
        "network_edges": int(wt_data.edge_index.shape[1]),
        "distance_metrics": list(DISTANCE_METRICS),
        "primary_distance": PRIMARY_DISTANCE,
        "matched_null_replicates": MATCHED_NULL_REPLICATES,
        "elapsed_minutes": (time.time() - start_time) / 60.0,
        "software": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "GenKI": GenKI.__version__,
            "torch": torch.__version__,
            "torch_geometric": torch_geometric.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
            "anndata": ad.__version__ if hasattr(ad, "__version__") else "0.12.19",
        },
        "interpretation_limit": (
            "GenKI distances are unsigned network-embedding perturbation magnitudes; "
            "direction is supplied only by existing empirical perturbation contrasts."
        ),
        "qa_all_passed": bool(qa["passed"].all()),
    }
    (OUT_DIR / "genki_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
