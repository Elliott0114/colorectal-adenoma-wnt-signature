#!/usr/bin/env python3
"""Build and transfer a discovery-only, donor-blocked adenoma route signature.

The legacy route signature used both Chen discovery and validation effects during
gene selection.  This analysis makes the validation contract explicit: genes are
selected from the discovery cohort only, using a donor-cluster bootstrap, and the
locked gene set is then evaluated without re-tuning in every downstream cohort.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from conventional_route_signature_transfer import (
    BECKER_DIR,
    CHEN_DATASETS,
    becker_tests_and_models,
    bh_adjust,
    chen_signature_scores,
    compare_groups,
    is_excluded_gene,
    score_atlas,
    score_becker,
    score_gse39582,
)


ROOT = Path(__file__).resolve().parents[1]
LEGACY_DIR = ROOT / "results" / "route_signature"
OUT_DIR = ROOT / "results" / "route_signature_locked"

BOOTSTRAP_REPLICATES = 1_000
RANDOM_SEED = 20260710
MIN_DIRECTION_STABILITY = 0.90
MIN_MEAN_EXPRESSION = 0.001
N_UP = 50
N_DOWN = 50


def load_chen_pseudobulk(dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load cached Chen specimen pseudobulk tables with aligned row order."""
    meta = pd.read_csv(
        LEGACY_DIR / f"chen_{dataset}_specimen_pseudobulk_meta.tsv",
        sep="\t",
    )
    expr = pd.read_csv(
        LEGACY_DIR / f"chen_{dataset}_specimen_pseudobulk_expression.tsv.gz",
        sep="\t",
    )
    if len(meta) != len(expr):
        raise ValueError(f"Chen {dataset}: metadata/expression row mismatch")
    if expr.columns.duplicated().any():
        raise ValueError(f"Chen {dataset}: duplicated gene columns")
    return meta, expr


def donor_route_expression(
    meta: pd.DataFrame,
    expr: pd.DataFrame,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Collapse specimens to one median expression vector per donor and route."""
    keep = meta["route_group"].isin(["normal", "conventional_adenoma"])
    route_meta = meta.loc[keep].reset_index(drop=True)
    route_expr = expr.loc[keep].reset_index(drop=True)
    keys = (
        route_meta[["donor_id", "route_group"]]
        .drop_duplicates()
        .sort_values(["donor_id", "route_group"])
        .reset_index(drop=True)
    )
    values = route_expr.to_numpy(dtype=np.float32)
    rows = []
    for key in keys.itertuples(index=False):
        idx = route_meta["donor_id"].eq(key.donor_id) & route_meta["route_group"].eq(key.route_group)
        rows.append(np.nanmedian(values[idx.to_numpy()], axis=0))
    return keys, np.vstack(rows).astype(np.float32)


def donor_cluster_bootstrap_coefficients(
    keys: pd.DataFrame,
    n_bootstrap: int,
    seed: int,
) -> np.ndarray:
    """Create route-contrast coefficients after resampling whole donor blocks."""
    donors = keys["donor_id"].astype(str).unique()
    row_donors = keys["donor_id"].astype(str).to_numpy()
    is_adenoma = keys["route_group"].eq("conventional_adenoma").to_numpy()
    rng = np.random.default_rng(seed)
    coefficients = np.zeros((n_bootstrap, len(keys)), dtype=np.float32)

    for bootstrap_id in range(n_bootstrap):
        sampled = rng.choice(donors, size=len(donors), replace=True)
        donor_ids, counts = np.unique(sampled, return_counts=True)
        count_by_donor = dict(zip(donor_ids, counts, strict=True))
        weights = np.array([count_by_donor.get(donor, 0) for donor in row_donors], dtype=np.float32)
        adenoma_weights = weights * is_adenoma
        normal_weights = weights * ~is_adenoma
        if adenoma_weights.sum() == 0 or normal_weights.sum() == 0:
            raise RuntimeError("A donor bootstrap replicate lost one route group")
        coefficients[bootstrap_id] = (
            adenoma_weights / adenoma_weights.sum()
            - normal_weights / normal_weights.sum()
        )
    return coefficients


def build_discovery_locked_signature(
    discovery_meta: pd.DataFrame,
    discovery_expr: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select stable genes without reading or using validation outcomes."""
    keys, expression = donor_route_expression(discovery_meta, discovery_expr)
    is_adenoma = keys["route_group"].eq("conventional_adenoma").to_numpy()
    full_effect = expression[is_adenoma].mean(axis=0) - expression[~is_adenoma].mean(axis=0)
    coefficients = donor_cluster_bootstrap_coefficients(
        keys,
        BOOTSTRAP_REPLICATES,
        RANDOM_SEED,
    )
    bootstrap_effects = coefficients @ expression
    positive_probability = (bootstrap_effects > 0).mean(axis=0)
    negative_probability = (bootstrap_effects < 0).mean(axis=0)
    direction_stability = np.where(full_effect >= 0, positive_probability, negative_probability)
    ci_low, ci_high = np.quantile(bootstrap_effects, [0.025, 0.975], axis=0)
    sign_p = 2 * np.minimum(positive_probability, negative_probability)
    sign_p = np.maximum(sign_p, 1 / (BOOTSTRAP_REPLICATES + 1))

    audit = pd.DataFrame(
        {
            "gene": discovery_expr.columns.astype(str),
            "discovery_effect_adenoma_minus_normal": full_effect,
            "bootstrap_effect_median": np.median(bootstrap_effects, axis=0),
            "bootstrap_ci_low": ci_low,
            "bootstrap_ci_high": ci_high,
            "direction_stability": direction_stability,
            "bootstrap_sign_p": sign_p,
            "mean_expression_discovery": expression.mean(axis=0),
        }
    )
    audit["bootstrap_sign_q"] = bh_adjust(audit["bootstrap_sign_p"])
    audit["excluded_gene"] = audit["gene"].map(is_excluded_gene)
    audit["ci_excludes_zero"] = (audit["bootstrap_ci_low"] > 0) | (audit["bootstrap_ci_high"] < 0)
    audit["eligible"] = (
        ~audit["excluded_gene"]
        & (audit["mean_expression_discovery"] > MIN_MEAN_EXPRESSION)
        & (audit["direction_stability"] >= MIN_DIRECTION_STABILITY)
        & audit["ci_excludes_zero"]
    )
    audit["selection_magnitude"] = audit["discovery_effect_adenoma_minus_normal"].abs()

    up = audit[
        audit["eligible"] & (audit["discovery_effect_adenoma_minus_normal"] > 0)
    ].nlargest(N_UP, "selection_magnitude").copy()
    down = audit[
        audit["eligible"] & (audit["discovery_effect_adenoma_minus_normal"] < 0)
    ].nlargest(N_DOWN, "selection_magnitude").copy()
    if len(up) < N_UP or len(down) < N_DOWN:
        raise RuntimeError(
            f"Insufficient stable genes: requested {N_UP}/{N_DOWN}, found {len(up)}/{len(down)}"
        )

    up["signature_direction"] = "adenoma_up"
    down["signature_direction"] = "adenoma_down"
    signature = pd.concat([up, down], ignore_index=True)
    signature["rank_within_direction"] = signature.groupby("signature_direction")[
        "selection_magnitude"
    ].rank(ascending=False, method="first")
    signature["selection_cohort"] = "Chen discovery only"
    signature["validation_used_for_selection"] = False
    signature = signature.sort_values(["signature_direction", "rank_within_direction", "gene"])
    audit["selected"] = audit["gene"].isin(signature["gene"])
    return signature, audit.sort_values(
        ["selected", "selection_magnitude", "gene"],
        ascending=[False, False, True],
    )


def focused_row(
    tests: pd.DataFrame,
    dataset: str,
    comparison: str = "conventional_vs_normal",
) -> dict[str, object]:
    row = tests[
        tests["dataset"].eq(dataset)
        & tests["comparison"].eq(comparison)
        & tests["score"].eq("ca_route_signature")
    ]
    if len(row) != 1:
        raise ValueError(f"Expected one focused Chen row for {dataset}, found {len(row)}")
    return row.iloc[0].to_dict()


def chen_dataset_tests(scores: pd.DataFrame) -> pd.DataFrame:
    """Keep discovery and held-out validation tests statistically separate."""
    score_columns = [
        "score__ca_route_signature",
        "score__ca_route_up",
        "score__ca_route_down",
        "score__wnt_stem",
        "score__proliferation_control",
    ]
    comparisons = [
        ("conventional_vs_normal", "conventional_adenoma", "normal"),
        ("serrated_vs_normal", "serrated", "normal"),
        ("conventional_vs_serrated", "conventional_adenoma", "serrated"),
    ]
    rows = []
    for dataset, frame in list(scores.groupby("dataset", observed=True)) + [("combined", scores)]:
        tests = compare_groups(
            frame[frame["route_group"].isin(["normal", "conventional_adenoma", "serrated"])],
            "route_group",
            score_columns,
            comparisons,
        )
        tests.insert(0, "dataset", dataset)
        rows.append(tests)
    return pd.concat(rows, ignore_index=True)


def chen_discrimination_metrics(scores: pd.DataFrame) -> pd.DataFrame:
    """Report label-free locked-score discrimination without optimizing a cut-point."""
    rows = []
    for dataset, frame in scores.groupby("dataset", observed=True):
        adenoma = frame.loc[
            frame["route_group"].eq("conventional_adenoma"),
            "score__ca_route_signature",
        ].dropna()
        normal = frame.loc[
            frame["route_group"].eq("normal"),
            "score__ca_route_signature",
        ].dropna()
        statistic = mannwhitneyu(adenoma, normal, alternative="two-sided")
        rows.append(
            {
                "dataset": dataset,
                "unit": "specimen",
                "n_adenoma": len(adenoma),
                "n_normal": len(normal),
                "auc_adenoma_vs_normal": float(statistic.statistic / (len(adenoma) * len(normal))),
                "p_mannwhitney": float(statistic.pvalue),
                "cutpoint_optimized": False,
            }
        )
    return pd.DataFrame(rows)


def write_summary(
    signature: pd.DataFrame,
    audit: pd.DataFrame,
    chen_tests: pd.DataFrame,
    discrimination: pd.DataFrame,
    overlap: pd.DataFrame,
) -> None:
    discovery = focused_row(chen_tests, "discovery")
    validation = focused_row(chen_tests, "validation")
    with (OUT_DIR / "discovery_locked_route_signature_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write("Discovery-only donor-blocked conventional adenoma route signature\n")
        handle.write("=" * 76 + "\n\n")
        handle.write("Selection contract\n")
        handle.write("- Gene selection used Chen discovery only.\n")
        handle.write("- Chen validation outcomes were not read by the selection function.\n")
        handle.write(f"- Donor-cluster bootstrap replicates: {BOOTSTRAP_REPLICATES}.\n")
        handle.write(f"- Direction-stability threshold: {MIN_DIRECTION_STABILITY:.2f}.\n")
        handle.write(f"- Locked genes: {len(signature)} ({N_UP} up, {N_DOWN} down).\n")
        handle.write(f"- Eligible stable discovery genes: {int(audit['eligible'].sum())}.\n")
        handle.write(f"- Overlap with legacy two-cohort-selected signature: {int(overlap['in_both'].sum())}/100.\n\n")
        handle.write("Focused Chen performance\n")
        for label, row in [("discovery", discovery), ("validation (held out)", validation)]:
            dataset = "validation" if label.startswith("validation") else "discovery"
            auc = discrimination.loc[discrimination["dataset"].eq(dataset), "auc_adenoma_vs_normal"].iloc[0]
            handle.write(
                f"- {label}: n={int(row['n_a'])} adenoma/{int(row['n_b'])} normal; "
                f"median delta={row['delta_a_minus_b']:.4f}; "
                f"AUC={auc:.3f}; P={row['p_mannwhitney']:.6g}; q={row['q_within_comparison']:.6g}.\n"
            )
        handle.write("\nTop adenoma-up genes\n")
        handle.write(", ".join(signature.loc[signature["signature_direction"].eq("adenoma_up"), "gene"].head(20)))
        handle.write("\n\nTop adenoma-down genes\n")
        handle.write(", ".join(signature.loc[signature["signature_direction"].eq("adenoma_down"), "gene"].head(20)))
        handle.write("\n\nInterpretation boundary\n")
        handle.write(
            "This is a discovery-locked transfer score, not a fitted recurrence predictor. "
            "Downstream cohorts are used for transportability or endpoint validation without gene re-selection.\n"
        )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    discovery_meta, discovery_expr = load_chen_pseudobulk("discovery")
    signature, audit = build_discovery_locked_signature(discovery_meta, discovery_expr)

    # Validation data are loaded only after the signature has been fully locked.
    validation_meta, validation_expr = load_chen_pseudobulk("validation")
    meta_by_dataset = {"discovery": discovery_meta, "validation": validation_meta}
    expr_by_dataset = {"discovery": discovery_expr, "validation": validation_expr}
    chen_scores, _, chen_paired = chen_signature_scores(
        meta_by_dataset,
        expr_by_dataset,
        signature,
    )
    chen_tests = chen_dataset_tests(chen_scores)
    discrimination = chen_discrimination_metrics(chen_scores)

    legacy = pd.read_csv(LEGACY_DIR / "conventional_adenoma_route_signature_genes.tsv", sep="\t")
    overlap = pd.DataFrame({"gene": sorted(set(signature["gene"]) | set(legacy["gene"]))})
    overlap["in_locked"] = overlap["gene"].isin(signature["gene"])
    overlap["in_legacy"] = overlap["gene"].isin(legacy["gene"])
    overlap["in_both"] = overlap["in_locked"] & overlap["in_legacy"]

    signature.to_csv(OUT_DIR / "discovery_locked_signature_genes.tsv", sep="\t", index=False)
    audit.to_csv(OUT_DIR / "discovery_gene_stability_audit.tsv", sep="\t", index=False)
    overlap.to_csv(OUT_DIR / "locked_vs_legacy_signature_overlap.tsv", sep="\t", index=False)
    chen_scores.to_csv(OUT_DIR / "chen_locked_signature_scores.tsv", sep="\t", index=False)
    chen_tests.to_csv(OUT_DIR / "chen_locked_signature_tests.tsv", sep="\t", index=False)
    chen_paired.to_csv(OUT_DIR / "chen_locked_signature_paired_tests.tsv", sep="\t", index=False)
    discrimination.to_csv(OUT_DIR / "chen_locked_signature_discrimination.tsv", sep="\t", index=False)

    becker_scores, becker_availability = score_becker(signature)
    becker_tests, becker_models = becker_tests_and_models(becker_scores)
    becker_scores.to_csv(OUT_DIR / "becker_locked_signature_scores.tsv", sep="\t", index=False)
    becker_availability.to_csv(OUT_DIR / "becker_locked_signature_gene_availability.tsv", sep="\t", index=False)
    becker_tests.to_csv(OUT_DIR / "becker_locked_signature_tests.tsv", sep="\t", index=False)
    becker_models.to_csv(OUT_DIR / "becker_locked_signature_adjusted_models.tsv", sep="\t", index=False)

    atlas_donors, atlas_tests, atlas_models = score_atlas(signature)
    atlas_donors.to_csv(OUT_DIR / "atlas_locked_signature_donor_scores.tsv", sep="\t", index=False)
    atlas_tests.to_csv(OUT_DIR / "atlas_locked_signature_tests.tsv", sep="\t", index=False)
    atlas_models.to_csv(OUT_DIR / "atlas_locked_signature_adjusted_models.tsv", sep="\t", index=False)

    gse_scores, gse_mapping, gse_selected, gse_results = score_gse39582(signature)
    gse_scores.to_csv(OUT_DIR / "gse39582_locked_signature_scores.tsv", sep="\t", index=False)
    gse_mapping.to_csv(OUT_DIR / "gse39582_locked_signature_probe_mapping_candidates.tsv", sep="\t", index=False)
    gse_selected.to_csv(OUT_DIR / "gse39582_locked_signature_selected_probes.tsv", sep="\t", index=False)
    gse_results.to_csv(OUT_DIR / "gse39582_locked_signature_clinical_results.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "discovery_locked_route_signature",
        "selection_dataset": "Chen discovery epithelial specimen pseudobulk",
        "held_out_dataset": "Chen validation epithelial specimen pseudobulk",
        "validation_used_for_selection": False,
        "bootstrap_unit": "donor block",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "random_seed": RANDOM_SEED,
        "direction_stability_threshold": MIN_DIRECTION_STABILITY,
        "signature_size": {"adenoma_up": N_UP, "adenoma_down": N_DOWN},
        "downstream_transfer": ["Chen validation", "Becker snRNA", "CRC Atlas", "GSE39582"],
        "source_cache": str(LEGACY_DIR.relative_to(ROOT)),
        "becker_source_present": BECKER_DIR.exists(),
        "chen_source_paths": {key: str(path.relative_to(ROOT)) for key, path in CHEN_DATASETS.items()},
    }
    (OUT_DIR / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_summary(signature, audit, chen_tests, discrimination, overlap)


if __name__ == "__main__":
    main()
