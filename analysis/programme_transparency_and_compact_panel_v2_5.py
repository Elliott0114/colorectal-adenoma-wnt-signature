#!/usr/bin/env python3
"""Make the 100-gene derivation and post hoc 10-gene nomination auditable.

This analysis does not replace or reselect the discovery-locked 100-gene
programme.  It adds (i) an explicit discovery selection flow, (ii) rank-based
Hallmark interpretation of the complete discovery universe, and (iii)
descriptive fidelity and robustness analyses for one fixed, biology-guided
10-gene candidate panel.  The random-panel benchmark is deliberately confined
to balanced, high-ranked and platform-complete panels and is not an independent
validation analysis.
"""

from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "analysis"
OUT_DIR = ROOT / "results" / "programme_transparency_v2_5"
RESOURCE_DIR = ROOT / "data_sources" / "reference_gene_sets_v2_5"
sys.path.insert(0, str(ANALYSIS_DIR))

import conventional_route_signature_transfer as transfer  # noqa: E402
import discovery_locked_route_signature as discovery  # noqa: E402
import external_sporadic_adenoma_validation as external  # noqa: E402
import gse117606_paired_route_validation as ffpe  # noqa: E402
import translation_reduced_panel_v2_0 as reduced_v20  # noqa: E402


SEED = 20260710
N_RANDOM_PANELS = 10_000
TOP_RANK_UNIVERSE = 20
MIN_HALLMARK_OVERLAP = 10

UP_GENES = ["OLFM4", "ASCL2", "RNF43", "NKD1", "AXIN2"]
DOWN_GENES = ["FABP1", "CA2", "PCK1", "LGALS4", "AQP8"]
PANEL_GENES = UP_GENES + DOWN_GENES

HALLMARK_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
    "2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt"
)
HALLMARK_SHA256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596"
HALLMARK_PATH = RESOURCE_DIR / "h.all.v2026.1.Hs.symbols.gmt"

BIOLOGICAL_ROLE = {
    "OLFM4": "intestinal stem/progenitor marker and tissue-level protein anchor",
    "ASCL2": "WNT-responsive intestinal stem-cell transcription factor",
    "RNF43": "WNT target and receptor-level negative-feedback regulator",
    "NKD1": "WNT-inducible intracellular negative-feedback regulator",
    "AXIN2": "canonical beta-catenin transcriptional readout and feedback regulator",
    "FABP1": "mature absorptive epithelial lipid-handling component",
    "CA2": "mature absorptive/ion-transport component and tissue-level protein anchor",
    "PCK1": "differentiated epithelial metabolic component",
    "LGALS4": "mature intestinal epithelial differentiation/adhesion component",
    "AQP8": "mature absorptive epithelial water-transport component",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, filename: str) -> None:
    frame.to_csv(OUT_DIR / filename, sep="\t", index=False)


def ensure_hallmark_resource() -> None:
    RESOURCE_DIR.mkdir(parents=True, exist_ok=True)
    if not HALLMARK_PATH.exists():
        opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
        with opener.open(HALLMARK_URL, timeout=60) as response:
            HALLMARK_PATH.write_bytes(response.read())
    observed = sha256(HALLMARK_PATH)
    if observed != HALLMARK_SHA256:
        raise RuntimeError(
            f"Hallmark resource checksum mismatch: expected {HALLMARK_SHA256}, observed {observed}"
        )
    manifest = pd.DataFrame(
        [
            {
                "resource": "MSigDB Hallmark",
                "version": "2026.1.Hs",
                "source_url": HALLMARK_URL,
                "local_path": str(HALLMARK_PATH.relative_to(ROOT)),
                "sha256": observed,
                "use": "rank-based interpretation of the full discovery gene universe",
                "redistribution_note": "obtain directly from MSigDB subject to its terms",
            }
        ]
    )
    manifest.to_csv(RESOURCE_DIR / "resource_manifest.tsv", sep="\t", index=False)


def discovery_selection_flow(audit: pd.DataFrame) -> pd.DataFrame:
    masks = [
        ("Assayed epithelial features", "all discovery pseudobulk features", pd.Series(True, index=audit.index)),
        (
            "Expressed, non-technical features",
            "not prespecified technical/exclusion genes; discovery mean expression > 0.001",
            (~audit["excluded_gene"]) & (audit["mean_expression_discovery"] > 0.001),
        ),
        (
            "Directionally stable features",
            "previous criteria plus donor-bootstrap directional stability >= 0.90",
            (~audit["excluded_gene"])
            & (audit["mean_expression_discovery"] > 0.001)
            & (audit["direction_stability"] >= 0.90),
        ),
        (
            "Eligible stable features",
            "previous criteria plus 95% donor-bootstrap interval excluding zero",
            audit["eligible"],
        ),
        (
            "Discovery-locked reference programme",
            "50 largest absolute discovery effects in each direction",
            audit["selected"],
        ),
    ]
    rows = []
    for order, (stage, criterion, mask) in enumerate(masks, 1):
        subset = audit.loc[mask]
        rows.append(
            {
                "stage_order": order,
                "stage": stage,
                "criterion": criterion,
                "n_features": int(len(subset)),
                "n_adenoma_up": int((subset["discovery_effect_adenoma_minus_normal"] > 0).sum()),
                "n_adenoma_down": int((subset["discovery_effect_adenoma_minus_normal"] < 0).sum()),
                "validation_outcomes_used": False,
            }
        )
    return pd.DataFrame(rows)


def discovery_design_audit() -> pd.DataFrame:
    meta, expr = discovery.load_chen_pseudobulk("discovery")
    keys, _ = discovery.donor_route_expression(meta, expr)
    donors = keys["donor_id"].astype(str)
    paired = keys.groupby("donor_id")["route_group"].nunique().eq(2)
    return pd.DataFrame(
        [
            {
                "discovery_specimens": int(len(meta)),
                "unique_donors": int(donors.nunique()),
                "adenoma_donors": int(keys.loc[keys["route_group"].eq("conventional_adenoma"), "donor_id"].nunique()),
                "normal_donors": int(keys.loc[keys["route_group"].eq("normal"), "donor_id"].nunique()),
                "paired_donors": int(paired.sum()),
                "donor_route_medians": int(len(keys)),
                "bootstrap_replicates": discovery.BOOTSTRAP_REPLICATES,
                "bootstrap_unit": "whole donor",
                "seed": discovery.RANDOM_SEED,
            }
        ]
    )


def read_gmt(path: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                gene_sets[parts[0]] = set(parts[2:])
    return gene_sets


def rank_based_hallmark_enrichment(audit: pd.DataFrame) -> pd.DataFrame:
    universe = audit.loc[
        (~audit["excluded_gene"]) & (audit["mean_expression_discovery"] > 0.001),
        ["gene", "discovery_effect_adenoma_minus_normal"],
    ].drop_duplicates("gene")
    effect = universe.set_index("gene")["discovery_effect_adenoma_minus_normal"]
    background = set(effect.index)
    rows = []
    for term, members in read_gmt(HALLMARK_PATH).items():
        present = sorted(members & background)
        absent = sorted(background - set(present))
        if len(present) < MIN_HALLMARK_OVERLAP:
            continue
        result = stats.mannwhitneyu(
            effect.loc[present], effect.loc[absent], alternative="two-sided", method="asymptotic"
        )
        auc = float(result.statistic / (len(present) * len(absent)))
        rows.append(
            {
                "hallmark": term,
                "n_members_in_discovery_universe": len(present),
                "discovery_universe_n": len(effect),
                "rank_auc": auc,
                "rank_biserial": 2 * auc - 1,
                "enriched_direction": "adenoma_up" if auc > 0.5 else "adenoma_down",
                "p_mann_whitney": float(result.pvalue),
                "member_genes_in_universe": ";".join(present),
            }
        )
    output = pd.DataFrame(rows)
    output["q_bh"] = multipletests(output["p_mann_whitney"], method="fdr_bh")[1]
    output["neg_log10_q"] = -np.log10(output["q_bh"].clip(lower=np.finfo(float).tiny))
    output["display_label"] = (
        output["hallmark"].str.removeprefix("HALLMARK_").str.replace("_", " ", regex=False).str.title()
    )
    return output.sort_values(["q_bh", "rank_biserial"], ascending=[True, False])


def corrected_panel_definition(locked: pd.DataFrame) -> pd.DataFrame:
    panel = locked.loc[locked["gene"].isin(PANEL_GENES)].copy()
    if set(panel["gene"]) != set(PANEL_GENES):
        raise RuntimeError("One or more compact-panel genes are absent from the locked programme")
    panel = panel.rename(columns={"rank_within_direction": "discovery_rank_within_direction"})
    panel["panel_order_within_arm"] = panel["gene"].map(
        {gene: rank for rank, gene in enumerate(UP_GENES, 1)}
        | {gene: rank for rank, gene in enumerate(DOWN_GENES, 1)}
    )
    panel["panel_arm"] = np.where(
        panel["signature_direction"].eq("adenoma_up"),
        "WNT_stem_progenitor_up",
        "mature_differentiation_down",
    )
    panel["route_weight"] = np.where(panel["signature_direction"].eq("adenoma_up"), 1.0, -1.0)
    panel["biological_role"] = panel["gene"].map(BIOLOGICAL_ROLE)
    panel["nomination_status"] = "post_hoc_biology_guided_compact_candidate"
    panel["nomination_constraints"] = (
        "locked-programme member; discovery rank <=20; interpretable arm role; "
        "complete principal-platform coverage; balanced 5+5 architecture"
    )
    panel["weights_fitted_to_outcome"] = False
    panel["cutpoint_optimised"] = False
    panel["gene_swapping_after_fixing_for_evaluation"] = False
    return panel.sort_values(
        ["route_weight", "panel_order_within_arm"], ascending=[False, True]
    )


def scoring_panel(panel: pd.DataFrame) -> pd.DataFrame:
    """Return the compatibility columns expected by existing score functions."""
    output = panel.copy()
    output["rank_within_direction"] = output["panel_order_within_arm"]
    return output


def load_platform_complete_data(
    locked: pd.DataFrame,
) -> tuple[
    dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    all_panel = locked.copy()
    all_panel["panel_arm"] = all_panel["signature_direction"]
    all_panel["route_weight"] = np.where(all_panel["signature_direction"].eq("adenoma_up"), 1, -1)
    external_data = reduced_v20.external_gene_data(all_panel)

    ffpe_meta = ffpe.parse_metadata()
    ffpe_raw = ffpe.parse_expression()
    if list(ffpe_raw.index) != ffpe_meta["sample_id"].tolist():
        raise RuntimeError("GSE117606 expression and metadata order mismatch")
    wanted = set(locked["gene"]) | set(external.PROLIFERATION_CONTROL)
    ffpe_mapping = ffpe.map_symbols_to_features(wanted, ffpe_raw)
    ffpe_expression = ffpe.select_gene_expression(ffpe_raw, ffpe_mapping)

    chen_meta, chen_expression = discovery.load_chen_pseudobulk("validation")
    return external_data, ffpe_meta, ffpe_expression, chen_meta, chen_expression


def platform_coverage_audit(
    locked: pd.DataFrame,
    external_data: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    ffpe_expression: pd.DataFrame,
    chen_expression: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    datasets = {
        "Chen held-out": chen_expression,
        **{cohort: expression for cohort, (_, expression) in external_data.items()},
        "GSE117606 FFPE": ffpe_expression,
    }
    for row in locked.itertuples(index=False):
        record = {
            "gene": row.gene,
            "signature_direction": row.signature_direction,
            "discovery_rank_within_direction": int(row.rank_within_direction),
        }
        for dataset, expression in datasets.items():
            record[f"present__{dataset}"] = row.gene in expression.columns
        present_columns = [key for key in record if key.startswith("present__")]
        record["n_principal_transcriptomic_resources_present"] = int(
            sum(bool(record[key]) for key in present_columns)
        )
        record["n_principal_transcriptomic_resources"] = len(present_columns)
        record["complete_principal_platform_coverage"] = all(
            bool(record[key]) for key in present_columns
        )
        rows.append(record)
    return pd.DataFrame(rows)


def candidate_audit(
    locked: pd.DataFrame,
    coverage: pd.DataFrame,
) -> pd.DataFrame:
    top = locked.loc[locked["rank_within_direction"].le(TOP_RANK_UNIVERSE)].copy()
    top = top.rename(columns={"rank_within_direction": "discovery_rank_within_direction"})
    top = top.merge(
        coverage[
            [
                "gene",
                "n_principal_transcriptomic_resources_present",
                "n_principal_transcriptomic_resources",
                "complete_principal_platform_coverage",
            ]
        ],
        on="gene",
        how="left",
        validate="one_to_one",
    )
    top["in_fixed_10_gene_candidate"] = top["gene"].isin(PANEL_GENES)
    top["biological_role_if_nominated"] = top["gene"].map(BIOLOGICAL_ROLE).fillna("")
    top["reporting_note"] = np.where(
        top["in_fixed_10_gene_candidate"],
        "included in the current post hoc biology-guided candidate",
        "retained in the 100-gene reference; not included in this compact candidate",
    )
    return top.sort_values(["signature_direction", "discovery_rank_within_direction"])


def zscore(expression: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    values = expression[genes].astype(float)
    return values.sub(values.mean(axis=0), axis=1).div(values.std(axis=0, ddof=1), axis=1)


def panel_score_matrix(
    expression: pd.DataFrame,
    genes: list[str],
    up_indices: np.ndarray,
    down_indices: np.ndarray,
) -> np.ndarray:
    z = zscore(expression, genes).to_numpy(dtype=float)
    weights = np.zeros((len(genes), len(up_indices)), dtype=float)
    column_index = np.arange(len(up_indices))[:, None]
    weights[up_indices, column_index] = 1.0 / up_indices.shape[1]
    weights[down_indices, column_index] = -1.0 / down_indices.shape[1]
    return z @ weights


def auc_columns(scores: np.ndarray, is_adenoma: np.ndarray, is_normal: np.ndarray) -> np.ndarray:
    output = np.empty(scores.shape[1], dtype=float)
    n_a = int(is_adenoma.sum())
    n_n = int(is_normal.sum())
    for column in range(scores.shape[1]):
        values = np.concatenate([scores[is_normal, column], scores[is_adenoma, column]])
        ranks = stats.rankdata(values, method="average")
        rank_sum_adenoma = ranks[n_n:].sum()
        u = rank_sum_adenoma - n_a * (n_a + 1) / 2
        output[column] = u / (n_a * n_n)
    return output


def external_fixed_effect_columns(scores: np.ndarray, meta: pd.DataFrame) -> np.ndarray:
    keep = meta["tissue_group"].isin(["normal", "adenoma"]).to_numpy()
    data = meta.loc[keep].reset_index(drop=True)
    values = scores[keep].copy()
    standardised = np.empty_like(values)
    for cohort in data["cohort"].unique():
        index = data["cohort"].eq(cohort).to_numpy()
        subset = values[index]
        standardised[index] = (subset - subset.mean(axis=0)) / subset.std(axis=0, ddof=1)
    x = data["tissue_group"].eq("adenoma").astype(float).to_numpy()
    x_residual = x.copy()
    for cohort in data["cohort"].unique():
        index = data["cohort"].eq(cohort).to_numpy()
        x_residual[index] -= x[index].mean()
    return (x_residual[:, None] * standardised).sum(axis=0) / np.square(x_residual).sum()


def ffpe_positive_fraction_columns(scores: np.ndarray, meta: pd.DataFrame) -> np.ndarray:
    differences = []
    for _, frame in meta.loc[meta["tissue_group"].isin(["normal", "adenoma"])].groupby(
        "patient_id", sort=False
    ):
        normal = frame.index[frame["tissue_group"].eq("normal")].to_numpy()
        adenoma = frame.index[frame["tissue_group"].eq("adenoma")].to_numpy()
        if len(normal) and len(adenoma):
            differences.append(np.median(scores[adenoma], axis=0) - np.median(scores[normal], axis=0))
    delta = np.vstack(differences)
    return (delta > 0).mean(axis=0)


def generate_random_panels(up_genes: list[str], down_genes: list[str]) -> tuple[list[tuple[str, ...]], np.ndarray, np.ndarray]:
    rng = np.random.default_rng(SEED)
    current = tuple(sorted(UP_GENES) + sorted(DOWN_GENES))
    panels: set[tuple[str, ...]] = set()
    while len(panels) < N_RANDOM_PANELS:
        up = tuple(sorted(rng.choice(up_genes, size=5, replace=False).tolist()))
        down = tuple(sorted(rng.choice(down_genes, size=5, replace=False).tolist()))
        panel = up + down
        if panel != current:
            panels.add(panel)
    ordered = sorted(panels)
    gene_order = up_genes + down_genes
    index = {gene: idx for idx, gene in enumerate(gene_order)}
    up_indices = np.array([[index[gene] for gene in panel[:5]] for panel in ordered], dtype=int)
    down_indices = np.array([[index[gene] for gene in panel[5:]] for panel in ordered], dtype=int)
    return ordered, up_indices, down_indices


def random_panel_benchmark(
    locked: pd.DataFrame,
    coverage: pd.DataFrame,
    external_data: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    ffpe_meta: pd.DataFrame,
    ffpe_expression: pd.DataFrame,
    chen_meta: pd.DataFrame,
    chen_expression: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = locked.merge(
        coverage[["gene", "complete_principal_platform_coverage"]],
        on="gene",
        validate="one_to_one",
    )
    eligible = eligible.loc[
        eligible["rank_within_direction"].le(TOP_RANK_UNIVERSE)
        & eligible["complete_principal_platform_coverage"]
    ].copy()
    up_genes = eligible.loc[eligible["signature_direction"].eq("adenoma_up")].sort_values(
        "rank_within_direction"
    )["gene"].tolist()
    down_genes = eligible.loc[eligible["signature_direction"].eq("adenoma_down")].sort_values(
        "rank_within_direction"
    )["gene"].tolist()
    if not set(PANEL_GENES).issubset(set(up_genes + down_genes)):
        raise RuntimeError("The current compact panel is outside the random-panel candidate universe")
    if len(up_genes) < 5 or len(down_genes) < 5:
        raise RuntimeError("Too few platform-complete genes for the balanced random-panel benchmark")

    panels, up_indices, down_indices = generate_random_panels(up_genes, down_genes)
    gene_order = up_genes + down_genes

    chen_matrix = panel_score_matrix(chen_expression, gene_order, up_indices, down_indices)
    chen_auc = auc_columns(
        chen_matrix,
        chen_meta["route_group"].eq("conventional_adenoma").to_numpy(),
        chen_meta["route_group"].eq("normal").to_numpy(),
    )

    external_meta = []
    external_matrices = []
    for cohort, (meta, expression) in external_data.items():
        frame = meta.copy().reset_index(drop=True)
        frame.insert(0, "cohort", cohort)
        external_meta.append(frame)
        external_matrices.append(panel_score_matrix(expression, gene_order, up_indices, down_indices))
    external_meta_frame = pd.concat(external_meta, ignore_index=True, sort=False)
    external_matrix = np.vstack(external_matrices)
    external_effect = external_fixed_effect_columns(external_matrix, external_meta_frame)

    ffpe_matrix = panel_score_matrix(ffpe_expression, gene_order, up_indices, down_indices)
    ffpe_meta_indexed = ffpe_meta.reset_index(drop=True)
    ffpe_positive = ffpe_positive_fraction_columns(ffpe_matrix, ffpe_meta_indexed)

    random_results = pd.DataFrame(
        {
            "panel_id": [f"random_{index:05d}" for index in range(1, len(panels) + 1)],
            "up_genes": [";".join(panel[:5]) for panel in panels],
            "down_genes": [";".join(panel[5:]) for panel in panels],
            "heldout_auc": chen_auc,
            "five_cohort_fixed_effect_sd": external_effect,
            "ffpe_positive_pair_fraction": ffpe_positive,
        }
    )

    current_up = np.array([[gene_order.index(gene) for gene in UP_GENES]], dtype=int)
    current_down = np.array([[gene_order.index(gene) for gene in DOWN_GENES]], dtype=int)
    current_chen = panel_score_matrix(chen_expression, gene_order, current_up, current_down)
    current_external = np.vstack(
        [
            panel_score_matrix(expression, gene_order, current_up, current_down)
            for _, expression in external_data.values()
        ]
    )
    current_ffpe = panel_score_matrix(ffpe_expression, gene_order, current_up, current_down)
    observed = {
        "heldout_auc": float(
            auc_columns(
                current_chen,
                chen_meta["route_group"].eq("conventional_adenoma").to_numpy(),
                chen_meta["route_group"].eq("normal").to_numpy(),
            )[0]
        ),
        "five_cohort_fixed_effect_sd": float(
            external_fixed_effect_columns(current_external, external_meta_frame)[0]
        ),
        "ffpe_positive_pair_fraction": float(
            ffpe_positive_fraction_columns(current_ffpe, ffpe_meta_indexed)[0]
        ),
    }
    summary_rows = []
    for metric, value in observed.items():
        random_values = random_results[metric]
        summary_rows.append(
            {
                "metric": metric,
                "observed_10_gene_panel": value,
                "random_median": float(random_values.median()),
                "random_q025": float(random_values.quantile(0.025)),
                "random_q975": float(random_values.quantile(0.975)),
                "empirical_percentile": float((random_values < value).mean()),
                "one_sided_empirical_tail_ge_observed": float(
                    ((random_values >= value).sum() + 1) / (len(random_values) + 1)
                ),
                "benchmark_panels": len(random_values),
                "benchmark_scope": (
                    f"balanced 5+5 panels sampled from platform-complete top-{TOP_RANK_UNIVERSE} "
                    "genes per locked arm"
                ),
            }
        )
    jointly_ge = (
        (random_results["heldout_auc"] >= observed["heldout_auc"])
        & (
            random_results["five_cohort_fixed_effect_sd"]
            >= observed["five_cohort_fixed_effect_sd"]
        )
        & (
            random_results["ffpe_positive_pair_fraction"]
            >= observed["ffpe_positive_pair_fraction"]
        )
    )
    summary_rows.append(
        {
            "metric": "jointly_meets_or_exceeds_all_three_observed_metrics",
            "observed_10_gene_panel": 1.0,
            "random_median": np.nan,
            "random_q025": np.nan,
            "random_q975": np.nan,
            "empirical_percentile": np.nan,
            "one_sided_empirical_tail_ge_observed": float(
                (jointly_ge.sum() + 1) / (len(jointly_ge) + 1)
            ),
            "benchmark_panels": len(jointly_ge),
            "benchmark_scope": (
                f"balanced 5+5 panels sampled from platform-complete top-{TOP_RANK_UNIVERSE} "
                "genes per locked arm"
            ),
        }
    )
    return random_results, pd.DataFrame(summary_rows)


def compact_panel_external_robustness() -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.read_csv(
        ROOT / "results" / "translation_reduced_panel_v2_0" / "external_reduced_sample_scores.tsv",
        sep="\t",
    )
    full = pd.read_csv(
        ROOT / "results" / "external_sporadic_adenoma_validation" / "sample_scores.tsv",
        sep="\t",
    )[["cohort", "sample_id", "score__proliferation_control"]]
    if "score__proliferation_control" in scores.columns:
        scores = scores.drop(columns="score__proliferation_control")
    scores = scores.merge(full, on=["cohort", "sample_id"], validate="one_to_one")
    adjusted_scores = scores.copy()
    adjusted_scores["route_score_k50"] = adjusted_scores["route_score_k5"]
    adjusted = pd.DataFrame([external.proliferation_adjusted_model(adjusted_scores)])
    adjusted["programme"] = "fixed_10_gene_candidate"

    rows = []
    for cohort in ["__NONE__", *scores["cohort"].unique().tolist()]:
        unadjusted = external.one_stage_model(scores, 5, excluded_cohort=cohort)
        proliferation = external.proliferation_adjusted_model(
            adjusted_scores, excluded_cohort=cohort
        )
        rows.append(
            {
                "excluded_cohort": cohort,
                "n_cohorts": unadjusted["n_cohorts"],
                "unadjusted_effect_sd": unadjusted["adenoma_coef_sd"],
                "unadjusted_ci_low": unadjusted["ci_low"],
                "unadjusted_ci_high": unadjusted["ci_high"],
                "unadjusted_p_value": unadjusted["p_value"],
                "proliferation_adjusted_effect_sd": proliferation["adenoma_coef_sd"],
                "proliferation_adjusted_ci_low": proliferation["ci_low"],
                "proliferation_adjusted_ci_high": proliferation["ci_high"],
                "proliferation_adjusted_p_value": proliferation["p_value"],
            }
        )
    return adjusted, pd.DataFrame(rows)


def gene_evidence_matrix(panel: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    ffpe_genes = pd.read_csv(
        ROOT / "results" / "gse117606_paired_route_validation" / "candidate_gene_paired_tests.tsv",
        sep="\t",
    )
    ffpe_genes = ffpe_genes.loc[ffpe_genes["gene"].isin(PANEL_GENES)].copy()
    recruitment = pd.read_csv(
        ROOT
        / "results"
        / "expanded_public_adenoma_validation"
        / "candidate_recruitment_cluster_summary.tsv",
        sep="\t",
    )
    recruitment = recruitment.loc[recruitment["gene"].isin(PANEL_GENES)].copy()
    proteins = pd.read_csv(
        ROOT
        / "results"
        / "public_adenoma_protein_triangulation"
        / "candidate_public_protein_evidence_matrix.tsv",
        sep="\t",
    )
    proteins = proteins.loc[proteins["gene"].isin(PANEL_GENES)].copy()
    apc = pd.read_csv(
        ROOT / "results" / "perturbation_validation_locked_route" / "gse125472_gene_effects.tsv",
        sep="\t",
    )
    apc = apc.loc[
        apc["gene"].isin(PANEL_GENES)
        & apc["comparison"].isin(["APC_vs_WT_with_Wnt", "APC_vs_WT_without_Wnt"])
    ].copy()
    apc["direction_match"] = np.sign(apc["mean_difference"]) == apc["expected_direction"]
    apc_summary = (
        apc.groupby("gene", as_index=False)
        .agg(
            apc_organoid_contrasts_evaluable=("direction_match", "size"),
            apc_organoid_direction_matches=("direction_match", "sum"),
        )
    )

    evidence = panel[
        [
            "gene",
            "panel_arm",
            "route_weight",
            "discovery_rank_within_direction",
            "discovery_effect_adenoma_minus_normal",
            "bootstrap_ci_low",
            "bootstrap_ci_high",
            "direction_stability",
            "biological_role",
        ]
    ].merge(
        coverage[
            [
                "gene",
                "n_principal_transcriptomic_resources_present",
                "n_principal_transcriptomic_resources",
                "complete_principal_platform_coverage",
            ]
        ],
        on="gene",
        validate="one_to_one",
    )
    evidence = evidence.merge(
        ffpe_genes[
            [
                "gene",
                "n_pairs",
                "median_paired_delta",
                "paired_positive_fraction",
                "gse117606_direction_match",
                "q_paired_wilcoxon_bh_candidates",
            ]
        ].rename(
            columns={
                "n_pairs": "ffpe_n_pairs",
                "median_paired_delta": "ffpe_median_paired_delta",
                "paired_positive_fraction": "ffpe_pair_direction_fraction_raw_above_zero",
                "gse117606_direction_match": "ffpe_direction_matches_expected",
                "q_paired_wilcoxon_bh_candidates": "ffpe_q_bh_candidates",
            }
        ),
        on="gene",
        how="left",
        validate="one_to_one",
    )
    evidence["ffpe_fraction_in_expected_direction"] = np.where(
        evidence["route_weight"].eq(1),
        evidence["ffpe_pair_direction_fraction_raw_above_zero"],
        1 - evidence["ffpe_pair_direction_fraction_raw_above_zero"],
    )
    evidence = evidence.merge(
        recruitment[
            [
                "gene",
                "n_recruitment_clusters_evaluable",
                "n_cluster_direction_matches",
                "cluster_direction_match_fraction",
            ]
        ],
        on="gene",
        how="left",
        validate="one_to_one",
    )
    evidence = evidence.merge(apc_summary, on="gene", how="left", validate="one_to_one")
    evidence = evidence.merge(
        proteins[
            [
                "gene",
                "protein_gene_present",
                "age_sex_adjusted_log2_effect",
                "age_sex_adjusted_q_bh_candidates",
                "pxd002137_adjusted_direction_match",
                "pxd017269_detection_fraction",
                "public_protein_evidence_tier",
            ]
        ],
        on="gene",
        how="left",
        validate="one_to_one",
    )
    evidence["protein_evidence_available"] = evidence["protein_gene_present"].fillna(False)
    return evidence.sort_values(["route_weight", "discovery_rank_within_direction"], ascending=[False, True])


def programme_size_context() -> pd.DataFrame:
    nested = pd.read_csv(
        ROOT
        / "results"
        / "external_sporadic_adenoma_validation"
        / "one_stage_patient_cluster_models.tsv",
        sep="\t",
    )
    nested = nested.loc[nested["excluded_cohort"].eq("__NONE__")].copy()
    reference = float(
        nested.loc[nested["signature_size_per_direction"].eq(50), "adenoma_coef_sd"].iloc[0]
    )
    nested["retained_fraction_vs_50_per_arm"] = nested["adenoma_coef_sd"] / reference
    nested["interpretation"] = np.where(
        nested["signature_size_per_direction"].eq(50),
        "fixed balanced reference breadth; not selected by validation optimisation",
        "prespecified nested sensitivity score",
    )
    return nested


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ensure_hallmark_resource()

    audit = pd.read_csv(
        ROOT / "results" / "route_signature_locked" / "discovery_gene_stability_audit.tsv",
        sep="\t",
    )
    locked = pd.read_csv(
        ROOT / "results" / "route_signature_locked" / "discovery_locked_signature_genes.tsv",
        sep="\t",
    )
    locked["rank_within_direction"] = pd.to_numeric(locked["rank_within_direction"])

    flow = discovery_selection_flow(audit)
    design = discovery_design_audit()
    hallmark = rank_based_hallmark_enrichment(audit)
    panel = corrected_panel_definition(locked)
    external_data, ffpe_meta, ffpe_expression, chen_meta, chen_expression = (
        load_platform_complete_data(locked)
    )
    coverage = platform_coverage_audit(
        locked, external_data, ffpe_expression, chen_expression
    )
    candidates = candidate_audit(locked, coverage)
    random_results, random_summary = random_panel_benchmark(
        locked,
        coverage,
        external_data,
        ffpe_meta,
        ffpe_expression,
        chen_meta,
        chen_expression,
    )
    adjusted, leave_cohort_out = compact_panel_external_robustness()
    evidence = gene_evidence_matrix(panel, coverage)
    nested = programme_size_context()

    write_tsv(flow, "reference_programme_selection_flow.tsv")
    write_tsv(design, "reference_programme_discovery_design.tsv")
    write_tsv(hallmark, "discovery_rank_based_hallmark_enrichment.tsv")
    write_tsv(panel, "compact_panel_definition_corrected.tsv")
    write_tsv(coverage, "locked_programme_principal_platform_coverage_by_gene.tsv")
    write_tsv(candidates, "compact_panel_top20_candidate_audit.tsv")
    write_tsv(random_results, "matched_random_balanced_panel_benchmark.tsv.gz")
    write_tsv(random_summary, "matched_random_balanced_panel_summary.tsv")
    write_tsv(adjusted, "compact_panel_proliferation_adjusted_model.tsv")
    write_tsv(leave_cohort_out, "compact_panel_leave_one_cohort_out.tsv")
    write_tsv(evidence, "compact_panel_gene_evidence_matrix.tsv")
    write_tsv(nested, "reference_programme_nested_size_context.tsv")

    wnt = hallmark.loc[hallmark["hallmark"].eq("HALLMARK_WNT_BETA_CATENIN_SIGNALING")].iloc[0]
    top_down = hallmark.sort_values("rank_biserial").iloc[0]
    manifest = {
        "analysis": "100-gene transparency and post hoc 10-gene candidate-panel audit",
        "primary_reference_programme": "100 genes; 50 per direction; discovery locked",
        "compact_panel_status": "post hoc, biology-guided, fixed for evaluation",
        "random_panel_seed": SEED,
        "random_panel_count": N_RANDOM_PANELS,
        "random_panel_universe": f"platform-complete top-{TOP_RANK_UNIVERSE} genes per locked arm",
        "hallmark_resource": "MSigDB Hallmark 2026.1.Hs",
        "hallmark_sha256": HALLMARK_SHA256,
        "key_results": {
            "eligible_stable_genes": int(flow.loc[flow["stage_order"].eq(4), "n_features"].iloc[0]),
            "selected_reference_genes": int(flow.loc[flow["stage_order"].eq(5), "n_features"].iloc[0]),
            "wnt_rank_auc": float(wnt["rank_auc"]),
            "wnt_q_bh": float(wnt["q_bh"]),
            "strongest_down_hallmark": str(top_down["hallmark"]),
            "strongest_down_rank_auc": float(top_down["rank_auc"]),
            "compact_panel_external_effect_sd": float(
                random_summary.loc[
                    random_summary["metric"].eq("five_cohort_fixed_effect_sd"),
                    "observed_10_gene_panel",
                ].iloc[0]
            ),
            "compact_panel_proliferation_adjusted_effect_sd": float(adjusted["adenoma_coef_sd"].iloc[0]),
        },
        "input_sha256": {
            "results/route_signature_locked/discovery_gene_stability_audit.tsv": sha256(
                ROOT / "results" / "route_signature_locked" / "discovery_gene_stability_audit.tsv"
            ),
            "results/route_signature_locked/discovery_locked_signature_genes.tsv": sha256(
                ROOT / "results" / "route_signature_locked" / "discovery_locked_signature_genes.tsv"
            ),
        },
    }
    (OUT_DIR / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    checks = {
        "selection_flow_counts": flow["n_features"].tolist() == [33698, 6127, 2496, 1504, 100],
        "balanced_reference": int(flow.iloc[-1]["n_adenoma_up"]) == 50
        and int(flow.iloc[-1]["n_adenoma_down"]) == 50,
        "panel_has_true_discovery_ranks": set(panel["discovery_rank_within_direction"].astype(int))
        == {1, 4, 5, 8, 9, 10, 12, 14, 16},
        "panel_is_platform_complete": bool(
            coverage.loc[coverage["gene"].isin(PANEL_GENES), "complete_principal_platform_coverage"].all()
        ),
        "random_panel_count": len(random_results) == N_RANDOM_PANELS,
        "random_effect_reproduces_exact_model": np.isclose(
            random_summary.loc[
                random_summary["metric"].eq("five_cohort_fixed_effect_sd"),
                "observed_10_gene_panel",
            ].iloc[0],
            pd.read_csv(
                ROOT
                / "results"
                / "translation_reduced_panel_v2_0"
                / "external_reduced_pooled_model.tsv",
                sep="\t",
            )["adenoma_coef_sd"].iloc[0],
            atol=1e-10,
        ),
        "wnt_is_top_positive_hallmark": hallmark.sort_values("rank_biserial", ascending=False).iloc[0][
            "hallmark"
        ]
        == "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    }
    qa = pd.DataFrame([{"check": key, "pass": value} for key, value in checks.items()])
    write_tsv(qa, "analysis_qa.tsv")
    if not qa["pass"].all():
        raise RuntimeError(
            "Transparency analysis QA failed: " + ", ".join(qa.loc[~qa["pass"], "check"])
        )
    print(f"Wrote {OUT_DIR} ({len(qa)} QA checks passed)")


if __name__ == "__main__":
    main()
