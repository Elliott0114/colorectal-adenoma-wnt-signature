#!/usr/bin/env python3
"""Evaluate the frozen compact readout with label-independent sample ranks.

Analysis date: 2026-08-30
Random seed: 20260830
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "analysis"
PARENT = ROOT / "results" / "state_aware_program_v1"
OUT_DIR = ROOT / "results" / "state_shared_revision_v2" / "external_rank"
COMMON = PARENT / "common_effects" / "cross_state_common_effects.tsv.gz"
PANEL = PARENT / "panel_derivation" / "compact_state_shared_panel_frozen.tsv"
PORTABLE = (
    PARENT / "panel_derivation" / "portable_state_shared_candidate_universe.tsv"
)
REFERENCE_SCORES = PARENT / "external_validation" / "external_sample_scores.tsv.gz"
CONTRACT = (
    ANALYSIS_DIR
    / "contracts"
    / "state_shared_revision_validation_v2_2026-08-30.md"
)
SEED = 20260830
N_RANDOM = 2000

sys.path.insert(0, str(ANALYSIS_DIR))
import translation_reduced_panel_v2_0 as reduced  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def row_percentile_ranks(expression: pd.DataFrame) -> pd.DataFrame:
    values = expression.to_numpy(dtype=float)
    ranks = stats.rankdata(values, axis=1, method="average") / values.shape[1]
    return pd.DataFrame(ranks, index=expression.index, columns=expression.columns)


def score_from_ranks(
    ranks: pd.DataFrame, up_genes: list[str], down_genes: list[str]
) -> np.ndarray:
    up = [gene for gene in up_genes if gene in ranks.columns]
    down = [gene for gene in down_genes if gene in ranks.columns]
    if not up or not down:
        raise ValueError("Both programme arms must be measurable")
    return ranks[up].mean(axis=1).to_numpy() - ranks[down].mean(axis=1).to_numpy()


def patient_values(frame: pd.DataFrame, group_column: str) -> pd.DataFrame:
    return (
        frame.groupby(["patient_cluster_id", group_column], observed=True)[
            "compact_rank_score"
        ]
        .mean()
        .reset_index()
    )


def contrast(
    frame: pd.DataFrame,
    group_column: str,
    group_a: str,
    group_b: str,
    name: str,
) -> dict[str, object] | None:
    values = patient_values(frame, group_column)
    values = values.loc[values[group_column].isin([group_a, group_b])]
    a = values.loc[values[group_column].eq(group_a), "compact_rank_score"].to_numpy()
    b = values.loc[values[group_column].eq(group_b), "compact_rank_score"].to_numpy()
    if len(a) < 3 or len(b) < 3:
        return None
    test = stats.ttest_ind(a, b, equal_var=False)
    variance = np.var(a, ddof=1) / len(a) + np.var(b, ddof=1) / len(b)
    standard_error = np.sqrt(variance)
    numerator = variance**2
    denominator = (
        (np.var(a, ddof=1) / len(a)) ** 2 / (len(a) - 1)
        + (np.var(b, ddof=1) / len(b)) ** 2 / (len(b) - 1)
    )
    degrees_of_freedom = numerator / denominator
    critical = stats.t.ppf(0.975, degrees_of_freedom)
    difference = float(np.mean(a) - np.mean(b))
    mann_whitney = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "contrast": name,
        "group_a": group_a,
        "group_b": group_b,
        "n_a": len(a),
        "n_b": len(b),
        "mean_difference": difference,
        "standard_error": standard_error,
        "ci_low": difference - critical * standard_error,
        "ci_high": difference + critical * standard_error,
        "welch_t": float(test.statistic),
        "degrees_of_freedom": degrees_of_freedom,
        "welch_p_value": float(test.pvalue),
        "mann_whitney_u": float(mann_whitney.statistic),
        "mann_whitney_p_value": float(mann_whitney.pvalue),
    }


def main() -> None:
    np.random.seed(SEED)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in [COMMON, PANEL, PORTABLE, REFERENCE_SCORES, CONTRACT]:
        if not path.exists():
            raise FileNotFoundError(path)

    common = pd.read_csv(COMMON, sep="\t")
    strict = common.loc[common["strict_state_shared"].astype(bool)].copy()
    strict["arm"] = strict["shared_direction"]
    strict["route_weight"] = np.where(strict["arm"].eq("up"), 1.0, -1.0)
    panel = pd.read_csv(PANEL, sep="\t")
    portable = pd.read_csv(PORTABLE, sep="\t")
    up_genes = sorted(strict.loc[strict["arm"].eq("up"), "gene"])
    down_genes = sorted(strict.loc[strict["arm"].eq("down"), "gene"])
    compact_up = panel.loc[panel["arm"].eq("up"), "gene"].tolist()
    compact_down = panel.loc[panel["arm"].eq("down"), "gene"].tolist()
    if len(strict) != 1843 or len(compact_up) != 4 or len(compact_down) != 4:
        raise ValueError("Frozen programme inputs have unexpected dimensions")

    gene_data = reduced.external_gene_data(strict[["gene"]])
    reference = pd.read_csv(REFERENCE_SCORES, sep="\t")
    reference = reference.loc[
        reference["signature_id"].eq("state_shared_1843"),
        ["cohort", "sample_id", "programme_score"],
    ].rename(columns={"programme_score": "reference_full_score"})

    scored_parts: list[pd.DataFrame] = []
    rank_matrices: dict[str, pd.DataFrame] = {}
    fidelity_rows: list[dict[str, object]] = []
    contrast_rows: list[dict[str, object]] = []
    for cohort, (meta, expression) in gene_data.items():
        ranks = row_percentile_ranks(expression)
        rank_matrices[cohort] = ranks
        scored = meta.copy().set_index("sample_id")
        scored["full_rank_score"] = score_from_ranks(ranks, up_genes, down_genes)
        scored["compact_rank_score"] = score_from_ranks(
            ranks, compact_up, compact_down
        )
        scored = scored.reset_index()
        scored.insert(0, "cohort", cohort)
        scored["patient_cluster_id"] = cohort + "::" + scored["patient_id"].astype(str)
        scored = scored.merge(
            reference.loc[reference["cohort"].eq(cohort)],
            on=["cohort", "sample_id"],
            how="left",
            validate="one_to_one",
        )
        if scored["reference_full_score"].isna().any():
            raise ValueError(f"Reference score alignment failed for {cohort}")
        correlation = stats.spearmanr(
            scored["compact_rank_score"], scored["reference_full_score"]
        )
        fidelity_rows.append(
            {
                "cohort": cohort,
                "n_samples": len(scored),
                "n_measurable_programme_genes": ranks.shape[1],
                "spearman_compact_rank_vs_reference_full": float(
                    correlation.statistic
                ),
                "p_value": float(correlation.pvalue),
            }
        )
        main_contrast = contrast(
            scored, "tissue_group", "adenoma", "normal", "adenoma_vs_normal"
        )
        if main_contrast is not None:
            main_contrast["cohort"] = cohort
            main_contrast["grouping_variable"] = "tissue_group"
            contrast_rows.append(main_contrast)
        if cohort == "GSE40362":
            for definition in [
                ("hyperplastic", "normal", "hyperplastic_vs_normal"),
                ("adenoma", "hyperplastic", "adenoma_vs_hyperplastic"),
            ]:
                result = contrast(scored, "tissue_group", *definition)
                if result is not None:
                    result["cohort"] = cohort
                    result["grouping_variable"] = "tissue_group"
                    contrast_rows.append(result)
        if cohort == "GSE41657":
            for definition in [
                ("low_grade", "normal", "low_grade_adenoma_vs_normal"),
                ("high_grade", "normal", "high_grade_adenoma_vs_normal"),
                ("high_grade", "low_grade", "high_vs_low_grade_adenoma"),
                ("crc", "normal", "crc_vs_normal"),
                ("crc", "high_grade", "crc_vs_high_grade_adenoma"),
            ]:
                result = contrast(scored, "grade_group", *definition)
                if result is not None:
                    result["cohort"] = cohort
                    result["grouping_variable"] = "grade_group"
                    contrast_rows.append(result)
        scored_parts.append(scored)

    scores = pd.concat(scored_parts, ignore_index=True, sort=False)
    fidelity = pd.DataFrame(fidelity_rows)
    contrasts = pd.DataFrame(contrast_rows)
    contrasts["welch_p_value_BH"] = contrasts.groupby("grouping_variable")[
        "welch_p_value"
    ].transform(lambda values: stats.false_discovery_control(values, method="bh"))

    portable_up = portable.loc[
        portable["arm"].eq("up") & portable["objective_selection_eligible"].astype(bool),
        "gene",
    ].tolist()
    portable_down = portable.loc[
        portable["arm"].eq("down") & portable["objective_selection_eligible"].astype(bool),
        "gene",
    ].tolist()
    rng = np.random.default_rng(SEED)
    random_records: list[dict[str, object]] = []
    panel_records: list[dict[str, object]] = []
    for iteration in range(1, N_RANDOM + 1):
        selected_up = rng.choice(portable_up, 4, replace=False).tolist()
        selected_down = rng.choice(portable_down, 4, replace=False).tolist()
        cohort_correlations = []
        for cohort, ranks in rank_matrices.items():
            random_score = score_from_ranks(ranks, selected_up, selected_down)
            full = scores.loc[scores["cohort"].eq(cohort), "reference_full_score"]
            cohort_correlations.append(
                float(stats.spearmanr(random_score, full).statistic)
            )
        random_records.append(
            {
                "random_panel_id": iteration,
                "median_cohort_spearman": float(np.median(cohort_correlations)),
                "minimum_cohort_spearman": float(np.min(cohort_correlations)),
                "maximum_cohort_spearman": float(np.max(cohort_correlations)),
            }
        )
        panel_records.extend(
            {
                "random_panel_id": iteration,
                "arm": arm,
                "gene": gene,
            }
            for arm, genes in (("up", selected_up), ("down", selected_down))
            for gene in genes
        )

    random_benchmark = pd.DataFrame(random_records)
    random_membership = pd.DataFrame(panel_records)
    observed_median = float(
        fidelity["spearman_compact_rank_vs_reference_full"].median()
    )
    random_q95 = float(random_benchmark["median_cohort_spearman"].quantile(0.95))
    benchmark_summary = pd.DataFrame(
        [
            {
                "observed_median_cohort_spearman": observed_median,
                "random_median": random_benchmark[
                    "median_cohort_spearman"
                ].median(),
                "random_q95": random_q95,
                "random_maximum": random_benchmark[
                    "median_cohort_spearman"
                ].max(),
                "empirical_upper_tail_p": (
                    1
                    + random_benchmark["median_cohort_spearman"]
                    .ge(observed_median)
                    .sum()
                )
                / (N_RANDOM + 1),
            }
        ]
    )

    main_route = contrasts.loc[contrasts["contrast"].eq("adenoma_vs_normal")]
    gates = pd.DataFrame(
        {
            "gate": [
                "all_external_cohort_rank_effects_positive",
                "median_external_fidelity_exceeds_random_q95",
            ],
            "passed": [
                bool(main_route["mean_difference"].gt(0).all()),
                bool(observed_median > random_q95),
            ],
        }
    )

    scores.to_csv(
        OUT_DIR / "external_single_sample_rank_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    fidelity.to_csv(OUT_DIR / "external_rank_fidelity.tsv", sep="\t", index=False)
    contrasts.to_csv(
        OUT_DIR / "external_rank_group_contrasts.tsv", sep="\t", index=False
    )
    random_benchmark.to_csv(
        OUT_DIR / "random_eight_gene_benchmark.tsv", sep="\t", index=False
    )
    random_membership.to_csv(
        OUT_DIR / "random_eight_gene_panel_membership.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    benchmark_summary.to_csv(
        OUT_DIR / "random_eight_gene_benchmark_summary.tsv", sep="\t", index=False
    )
    gates.to_csv(OUT_DIR / "quality_gates.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "state_shared_revision_external_rank_v2",
        "created_utc": pd.Timestamp.utcnow().isoformat(),
        "random_seed": SEED,
        "input_sha256": {
            "common_effects": sha256(COMMON),
            "compact_panel": sha256(PANEL),
            "portable_universe": sha256(PORTABLE),
            "reference_scores": sha256(REFERENCE_SCORES),
            "contract": sha256(CONTRACT),
        },
        "score_definition": (
            "within-sample percentile ranks among all measurable frozen programme "
            "genes; mean rank of compact up genes minus compact down genes"
        ),
        "label_dependent_standardisation": False,
        "random_panels": N_RANDOM,
        "quality_gates": dict(zip(gates["gate"], gates["passed"], strict=True)),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }
    with (OUT_DIR / "analysis_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)

    print(
        f"External rank analysis completed: median compact fidelity={observed_median:.3f}; "
        f"random 95th percentile={random_q95:.3f}."
    )


if __name__ == "__main__":
    main()
