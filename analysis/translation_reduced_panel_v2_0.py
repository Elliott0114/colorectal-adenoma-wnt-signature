#!/usr/bin/env python3
"""Evaluate a fixed 10-gene translational reduction of the locked programme.

The 100-gene discovery-locked programme remains the primary research readout.
This script evaluates one post hoc, non-reselected 10-gene panel for score
concordance, retained effects, FFPE measurability, perturbation direction and
leave-one-gene-out stability.  It does not optimise a diagnostic threshold.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "analysis"
OUT_DIR = ROOT / "results" / "translation_reduced_panel_v2_0"
sys.path.insert(0, str(ANALYSIS_DIR))

import conventional_route_signature_transfer as transfer  # noqa: E402
import discovery_locked_route_signature as discovery  # noqa: E402
import external_sporadic_adenoma_validation as external  # noqa: E402
import gse117606_paired_route_validation as ffpe  # noqa: E402
import perturbation_validation_locked_route as perturbation  # noqa: E402


SEED = 20260710
UP_GENES = ["OLFM4", "ASCL2", "RNF43", "NKD1", "AXIN2"]
DOWN_GENES = ["FABP1", "CA2", "PCK1", "LGALS4", "AQP8"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT_DIR / name, sep="\t", index=False)


def reduced_signature() -> pd.DataFrame:
    locked = pd.read_csv(
        ROOT / "results" / "route_signature_locked" / "discovery_locked_signature_genes.tsv",
        sep="\t",
    )
    genes = UP_GENES + DOWN_GENES
    panel = locked.loc[locked["gene"].isin(genes)].copy()
    if set(panel["gene"]) != set(genes):
        missing = sorted(set(genes) - set(panel["gene"]))
        raise RuntimeError(f"Reduced-panel genes absent from locked programme: {missing}")
    panel["panel_arm"] = np.where(
        panel["gene"].isin(UP_GENES), "WNT_stem_progenitor_up", "mature_differentiation_down"
    )
    panel["route_weight"] = np.where(panel["gene"].isin(UP_GENES), 1.0, -1.0)
    panel["reduced_rank_within_arm"] = panel["gene"].map(
        {gene: rank for rank, gene in enumerate(UP_GENES, 1)}
        | {gene: rank for rank, gene in enumerate(DOWN_GENES, 1)}
    )
    panel["selection_status"] = "post_hoc_fixed_exploratory_reduction"
    panel["gene_swapping_permitted"] = False
    panel["signature_direction"] = np.where(
        panel["gene"].isin(UP_GENES), "adenoma_up", "adenoma_down"
    )
    panel["rank_within_direction"] = panel["reduced_rank_within_arm"]
    return panel.sort_values(["route_weight", "reduced_rank_within_arm"], ascending=[False, True])


def spearman_row(scope: str, cohort: str, reduced: pd.Series, full: pd.Series) -> dict[str, object]:
    keep = reduced.notna() & full.notna()
    result = stats.spearmanr(reduced.loc[keep], full.loc[keep])
    return {
        "scope": scope,
        "cohort": cohort,
        "n": int(keep.sum()),
        "spearman_rho_reduced_vs_100_gene": float(result.statistic),
        "p_value": float(result.pvalue),
    }


def auc_and_p(values: pd.DataFrame, score: str, group: str = "route_group") -> dict[str, float]:
    adenoma = values.loc[values[group].isin(["conventional_adenoma", "adenoma"]), score].dropna()
    normal = values.loc[values[group].eq("normal"), score].dropna()
    test = stats.mannwhitneyu(adenoma, normal, alternative="two-sided")
    return {
        "n_adenoma": int(len(adenoma)),
        "n_normal": int(len(normal)),
        "auc": float(test.statistic / (len(adenoma) * len(normal))),
        "p_value": float(test.pvalue),
        "median_difference": float(adenoma.median() - normal.median()),
    }


def chen_analysis(panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    meta_by_dataset: dict[str, pd.DataFrame] = {}
    expr_by_dataset: dict[str, pd.DataFrame] = {}
    for dataset in ["discovery", "validation"]:
        meta, expr = discovery.load_chen_pseudobulk(dataset)
        meta_by_dataset[dataset] = meta
        expr_by_dataset[dataset] = expr

    reduced_scores, _, reduced_paired = transfer.chen_signature_scores(
        meta_by_dataset, expr_by_dataset, panel
    )
    reduced_scores = reduced_scores.rename(
        columns={
            "score__ca_route_up": "reduced_up",
            "score__ca_route_down": "reduced_down",
            "score__ca_route_signature": "reduced_score",
            "score__ca_route_n_up_present": "reduced_n_up_present",
            "score__ca_route_n_down_present": "reduced_n_down_present",
        }
    )
    full_scores = pd.read_csv(
        ROOT / "results" / "route_signature_locked" / "chen_locked_signature_scores.tsv",
        sep="\t",
    )[["dataset", "specimen_id", "score__ca_route_signature"]]
    reduced_scores = reduced_scores.merge(
        full_scores, on=["dataset", "specimen_id"], how="left", validate="one_to_one"
    ).rename(columns={"score__ca_route_signature": "full_100_gene_score"})

    perf_rows = []
    corr_rows = []
    for dataset, frame in reduced_scores.groupby("dataset", sort=False):
        perf_rows.append({"dataset": dataset, **auc_and_p(frame, "reduced_score")})
        corr_rows.append(
            spearman_row("Chen", dataset, frame["reduced_score"], frame["full_100_gene_score"])
        )
    performance = pd.DataFrame(perf_rows)
    correlations = pd.DataFrame(corr_rows)
    return reduced_scores, performance, reduced_paired, expr_by_dataset


def external_gene_data(panel: pd.DataFrame) -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    wanted = set(panel["gene"]) | set(external.PROLIFERATION_CONTROL)
    mappings = {
        "GSE8671": external.gpl570_mapping(),
        "GSE50114": external.gpl6480_mapping(),
        "GSE41657": external.gpl6480_mapping(),
        "GSE40362": external.gpl8432_mapping(),
    }
    log2_transform = {
        "GSE8671": True,
        "GSE50114": False,
        "GSE41657": False,
        "GSE40362": False,
    }
    output: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    for accession in ["GSE8671", "GSE50114", "GSE41657", "GSE40362"]:
        meta, expression, _ = external.parse_series_matrix(accession)
        meta = external.add_design_fields(accession, meta)
        if accession == "GSE50114":
            expression = external.gse50114_raw_expression()
        if log2_transform[accession]:
            expression = np.log2(expression.clip(lower=0) + 1)
        gene_expression, _ = external.choose_features_without_labels(
            expression, mappings[accession], wanted, accession
        )
        output[accession] = (meta, gene_expression)
    meta, expression, _ = external.direct_rnaseq_expression(wanted)
    output["GSE72820"] = (meta, expression)
    return output


def score_gene_data(
    meta: pd.DataFrame,
    expression: pd.DataFrame,
    panel: pd.DataFrame,
    cohort: str,
    score_column: str = "route_score_k5",
) -> pd.DataFrame:
    z = expression.sub(expression.mean(axis=0), axis=1).div(
        expression.std(axis=0, ddof=1).replace(0, np.nan), axis=1
    )
    up = panel.loc[panel["route_weight"].eq(1), "gene"].tolist()
    down = panel.loc[panel["route_weight"].eq(-1), "gene"].tolist()
    if not up or not down:
        raise RuntimeError("Every reduced score must retain both biological arms")
    if not set(up + down).issubset(z.columns):
        missing = sorted(set(up + down) - set(z.columns))
        raise RuntimeError(f"{cohort}: missing reduced-panel genes: {missing}")
    out = meta.copy().set_index("sample_id")
    out[score_column] = z[up].mean(axis=1) - z[down].mean(axis=1)
    out = out.reset_index()
    out.insert(0, "cohort", cohort)
    out["patient_cluster_id"] = cohort + "::" + out["patient_id"].astype(str)
    return out


def external_analysis(
    panel: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, tuple[pd.DataFrame, pd.DataFrame]]]:
    gene_data = external_gene_data(panel)
    scores = pd.concat(
        [score_gene_data(meta, expr, panel, cohort) for cohort, (meta, expr) in gene_data.items()],
        ignore_index=True,
        sort=False,
    )
    current = pd.read_csv(
        ROOT / "results" / "external_sporadic_adenoma_validation" / "sample_scores.tsv",
        sep="\t",
    )[["cohort", "sample_id", "route_score_k50"]]
    scores = scores.merge(current, on=["cohort", "sample_id"], how="left", validate="one_to_one")

    paired = {"GSE8671", "GSE72820"}
    tests = pd.DataFrame(
        [
            external.cohort_comparison(
                frame,
                cohort,
                5,
                "adenoma",
                "normal",
                cohort in paired,
                "adenoma_vs_normal",
            )
            for cohort, frame in scores.groupby("cohort", sort=False)
        ]
    )
    pooled = pd.DataFrame([external.one_stage_model(scores, 5)])
    current_pooled = pd.read_csv(
        ROOT
        / "results"
        / "external_sporadic_adenoma_validation"
        / "one_stage_patient_cluster_models.tsv",
        sep="\t",
    )
    current_pooled = current_pooled.loc[
        current_pooled["signature_size_per_direction"].eq(50)
    ].iloc[0]
    pooled["full_100_gene_effect_sd"] = float(current_pooled["adenoma_coef_sd"])
    pooled["retained_effect_fraction"] = (
        pooled["adenoma_coef_sd"] / pooled["full_100_gene_effect_sd"]
    )
    correlations = pd.DataFrame(
        [
            spearman_row(
                "external_sporadic",
                cohort,
                frame["route_score_k5"],
                frame["route_score_k50"],
            )
            for cohort, frame in scores.groupby("cohort", sort=False)
        ]
    )
    return scores, tests, pooled, correlations, gene_data


def ffpe_analysis(
    panel: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    meta = ffpe.parse_metadata()
    expression = ffpe.parse_expression()
    if list(expression.index) != meta["sample_id"].tolist():
        raise RuntimeError("GSE117606 expression and metadata order mismatch")
    wanted = set(panel["gene"]) | set(external.PROLIFERATION_CONTROL)
    mapping = ffpe.map_symbols_to_features(wanted, expression)
    gene_expression = ffpe.select_gene_expression(expression, mapping)
    scores = score_gene_data(meta, gene_expression, panel, "GSE117606")
    current = pd.read_csv(
        ROOT
        / "results"
        / "gse117606_paired_route_validation"
        / "sample_scores_all_cohort_scaling.tsv",
        sep="\t",
    )[["sample_id", "route_score_k50"]]
    scores = scores.merge(current, on="sample_id", how="left", validate="one_to_one")
    test = pd.DataFrame(
        [
            external.cohort_comparison(
                scores,
                "GSE117606",
                5,
                "adenoma",
                "normal",
                True,
                "conventional_adenoma_vs_adjacent_mucosa",
            )
        ]
    )
    corr = pd.DataFrame(
        [
            spearman_row(
                "FFPE",
                "GSE117606",
                scores["route_score_k5"],
                scores["route_score_k50"],
            )
        ]
    )
    return scores, test, corr, gene_expression


def perturbation_analysis(panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(perturbation.GSE125_PATH, sep="\t", low_memory=False)
    sample_columns = [
        column
        for column in raw.columns
        if perturbation.re.match(r"^Donor\d+_(WT|APC)_(w|wo)_Wnt$", column)
    ]
    counts = raw[sample_columns].apply(pd.to_numeric, errors="raise")
    symbols = raw["hgnc_symbol"].fillna(raw["symbol"])
    counts = perturbation.aggregate_symbol_counts(counts, symbols)
    expression = perturbation.log2_cpm(counts)
    meta = perturbation.gse125_metadata(sample_columns)
    scores, coverage, _ = perturbation.score_expression(expression, panel)
    scores = meta.merge(scores, on="sample_id", validate="one_to_one")
    _, summary, _ = perturbation.gse125_contrasts(scores)
    summary = summary.loc[summary["feature"].eq("route_score")].copy()
    summary.insert(0, "dataset", "GSE125472")
    summary["direction_matches_expected"] = np.where(
        summary["expected_direction"].eq(0),
        np.nan,
        np.sign(summary["mean_difference"]) == summary["expected_direction"],
    )

    mapping = perturbation.ensembl_to_symbol()
    clone_rows = []
    sample_rows = []
    for cell_line, filename in [
        ("HCT116", "GSE135328_count_HCT116.txt.gz"),
        ("HT29", "GSE135328_count_HT29.txt.gz"),
    ]:
        counts = perturbation.read_gse135_counts(perturbation.GSE135_DIR / filename, mapping)
        expression = perturbation.log2_cpm(counts)
        meta135 = perturbation.gse135_metadata(cell_line, expression.columns.tolist())
        score135, _, _ = perturbation.score_expression(expression, panel)
        score135 = meta135.merge(score135, on="sample_id", validate="one_to_one")
        sample_rows.append(score135)
        clone = (
            score135.groupby(["cell_line", "clone_id", "clone_number", "genotype"], observed=True)[
                "route_score"
            ]
            .mean()
            .reset_index()
        )
        wt = clone.loc[clone["genotype"].eq("WT"), "route_score"].iloc[0]
        for row in clone.loc[clone["genotype"].isin(["KO", "Het"])].itertuples(index=False):
            clone_rows.append(
                {
                    "dataset": "GSE135328",
                    "cell_line": cell_line,
                    "clone_id": row.clone_id,
                    "genotype": row.genotype,
                    "reduced_route_score": row.route_score,
                    "difference_vs_WT": row.route_score - wt,
                    "expected_direction": -1,
                    "direction_matches_expected": row.route_score - wt < 0,
                }
            )
    perturbation_scores = pd.concat(
        [scores.assign(dataset="GSE125472"), *[x.assign(dataset="GSE135328") for x in sample_rows]],
        ignore_index=True,
        sort=False,
    )
    coverage["n_expected"] = coverage["feature"].map(
        {"route_up": 5, "route_down": 5}
    ).fillna(coverage["n_expected"])
    coverage["coverage_fraction"] = coverage["n_present"] / coverage["n_expected"]
    return perturbation_scores, summary, pd.DataFrame(clone_rows)


def paired_ffpe_metrics(scores: pd.DataFrame, score_column: str) -> dict[str, object]:
    wide = (
        scores.loc[scores["tissue_group"].isin(["normal", "adenoma"])]
        .pivot_table(index="patient_id", columns="tissue_group", values=score_column, aggfunc="median")
        .dropna()
    )
    delta = wide["adenoma"] - wide["normal"]
    test = stats.wilcoxon(delta, zero_method="wilcox", alternative="two-sided")
    return {
        "n_pairs": int(len(delta)),
        "median_paired_difference": float(delta.median()),
        "positive_fraction": float((delta > 0).mean()),
        "p_value": float(test.pvalue),
    }


def heldout_metrics(meta: pd.DataFrame, expression: pd.DataFrame, panel: pd.DataFrame) -> dict[str, object]:
    up = panel.loc[panel["route_weight"].eq(1), "gene"].tolist()
    down = panel.loc[panel["route_weight"].eq(-1), "gene"].tolist()
    score = transfer.route_score_from_expression(expression, up, down)["score__ca_route_signature"]
    frame = meta.copy()
    frame["score"] = score.to_numpy()
    auc = auc_and_p(frame, "score")
    donor = (
        frame.loc[frame["route_group"].isin(["normal", "conventional_adenoma"])]
        .groupby(["donor_id", "route_group"], observed=True)["score"]
        .median()
        .unstack()
        .dropna()
    )
    delta = donor["conventional_adenoma"] - donor["normal"]
    auc.update(
        {
            "n_paired_donors": int(len(delta)),
            "paired_median_difference": float(delta.median()),
            "paired_p_value": float(
                stats.wilcoxon(delta, zero_method="wilcox", alternative="two-sided").pvalue
            ),
        }
    )
    return auc


def leave_one_gene_out(
    panel: pd.DataFrame,
    external_data: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    chen_validation_meta: pd.DataFrame,
    chen_validation_expr: pd.DataFrame,
    ffpe_meta_scores: pd.DataFrame,
    ffpe_gene_expression: pd.DataFrame,
    full_reduced_external_effect: float,
) -> pd.DataFrame:
    rows = []
    for omitted in panel["gene"].tolist():
        subset = panel.loc[panel["gene"].ne(omitted)].copy()
        ext_scores = pd.concat(
            [
                score_gene_data(meta, expr, subset, cohort)
                for cohort, (meta, expr) in external_data.items()
            ],
            ignore_index=True,
            sort=False,
        )
        model = external.one_stage_model(ext_scores, 5)
        chen = heldout_metrics(chen_validation_meta, chen_validation_expr, subset)
        ffpe_score = score_gene_data(
            ffpe_meta_scores.drop(columns=[
                column
                for column in ["cohort", "route_score_k5", "route_score_k50", "patient_cluster_id"]
                if column in ffpe_meta_scores.columns
            ]),
            ffpe_gene_expression,
            subset,
            "GSE117606",
        )
        ffpe_metrics = paired_ffpe_metrics(ffpe_score, "route_score_k5")
        rows.append(
            {
                "omitted_gene": omitted,
                "omitted_arm": panel.loc[panel["gene"].eq(omitted), "panel_arm"].iloc[0],
                "n_genes_retained": len(subset),
                "external_pooled_effect_sd": model["adenoma_coef_sd"],
                "external_ci_low": model["ci_low"],
                "external_ci_high": model["ci_high"],
                "external_retained_fraction_vs_10_gene": model["adenoma_coef_sd"]
                / full_reduced_external_effect,
                "heldout_auc": chen["auc"],
                "heldout_p_value": chen["p_value"],
                "ffpe_n_pairs": ffpe_metrics["n_pairs"],
                "ffpe_median_paired_difference": ffpe_metrics["median_paired_difference"],
                "ffpe_positive_fraction": ffpe_metrics["positive_fraction"],
                "ffpe_p_value": ffpe_metrics["p_value"],
            }
        )
    return pd.DataFrame(rows)


def protein_anchor_table() -> pd.DataFrame:
    protein = pd.read_csv(
        ROOT
        / "results"
        / "public_adenoma_protein_triangulation"
        / "candidate_public_protein_evidence_matrix.tsv",
        sep="\t",
    )
    protein = protein.loc[protein["gene"].isin(["OLFM4", "CA2", "FABP1"])].copy()
    paired = pd.read_csv(
        ROOT / "results" / "pxd000445_candidate_reanalysis" / "psm_candidate_paired_tests.tsv",
        sep="\t",
    )
    paired = paired.loc[
        paired["analysis_set"].eq("author_qc_21_pairs")
        & paired["gene"].isin(["OLFM4", "CA2", "FABP1"])
    ].copy()
    keep = [
        "gene",
        "expected_direction",
        "age_sex_adjusted_log2_effect",
        "age_sex_adjusted_ci_low",
        "age_sex_adjusted_ci_high",
        "age_sex_adjusted_q_bh_candidates",
        "pxd017269_detection_fraction",
        "detected_in_nine_patient_dvp",
        "public_protein_evidence_tier",
    ]
    out = protein[keep].merge(
        paired[
            [
                "gene",
                "n_pairs",
                "median_adenoma_minus_normal_log2",
                "paired_positive_fraction",
                "p_paired_wilcoxon",
                "q_value_bh_within_analysis_set",
            ]
        ],
        on="gene",
        how="left",
    )
    out["primary_tissue_anchor"] = out["gene"].isin(["OLFM4", "CA2"])
    out["evidence_role"] = out["gene"].map(
        {
            "OLFM4": "primary_up_arm_protein_anchor",
            "CA2": "primary_down_arm_protein_anchor",
            "FABP1": "complementary_down_arm_paired_support",
        }
    )
    return out.sort_values(["primary_tissue_anchor", "gene"], ascending=[False, True])


def gate_summary(
    performance: pd.DataFrame,
    correlations: pd.DataFrame,
    pooled: pd.DataFrame,
    ffpe_test: pd.DataFrame,
    perturb_summary: pd.DataFrame,
    clone_effects: pd.DataFrame,
    loo: pd.DataFrame,
    external_scores: pd.DataFrame,
) -> pd.DataFrame:
    current_auc = pd.read_csv(
        ROOT / "results" / "route_signature_locked" / "chen_locked_signature_discrimination.tsv",
        sep="\t",
    ).set_index("dataset").loc["validation", "auc_adenoma_vs_normal"]
    reduced_auc = performance.set_index("dataset").loc["validation", "auc"]
    expected_gse125 = perturb_summary.loc[perturb_summary["expected_direction"].ne(0)]
    expected_clone = clone_effects.loc[clone_effects["genotype"].eq("KO")]
    coverage_counts = external_scores.groupby("cohort").agg(
        n_up=("route_score_k5", lambda x: 5), n_down=("route_score_k5", lambda x: 5)
    )
    rows = [
        {
            "gate": "five_cohort_complete_10_of_10_coverage",
            "criterion": "10/10 genes available in each of five external cohorts",
            "observed": f"{len(coverage_counts)}/5 cohorts; 10/10 each",
            "pass": len(coverage_counts) == 5,
        },
        {
            "gate": "score_concordance",
            "criterion": "minimum external Spearman rho >= 0.80",
            "observed": float(
                correlations.loc[correlations["scope"].eq("external_sporadic"), "spearman_rho_reduced_vs_100_gene"].min()
            ),
            "pass": bool(
                correlations.loc[correlations["scope"].eq("external_sporadic"), "spearman_rho_reduced_vs_100_gene"].min()
                >= 0.80
            ),
        },
        {
            "gate": "heldout_auc_retention",
            "criterion": "held-out AUC drop <= 0.05",
            "observed": float(current_auc - reduced_auc),
            "pass": bool(current_auc - reduced_auc <= 0.05),
        },
        {
            "gate": "external_effect_retention",
            "criterion": "pooled effect retains >= 80% of 100-gene effect",
            "observed": float(pooled["retained_effect_fraction"].iloc[0]),
            "pass": bool(pooled["retained_effect_fraction"].iloc[0] >= 0.80),
        },
        {
            "gate": "ffpe_pair_direction",
            "criterion": ">= 90% positive paired changes and two-sided P < 0.05",
            "observed": (
                f"{ffpe_test['paired_positive_fraction'].iloc[0]:.3f}; "
                f"P={ffpe_test['p_paired_wilcoxon'].iloc[0]:.3g}"
            ),
            "pass": bool(
                ffpe_test["paired_positive_fraction"].iloc[0] >= 0.90
                and ffpe_test["p_paired_wilcoxon"].iloc[0] < 0.05
            ),
        },
        {
            "gate": "prespecified_perturbation_direction",
            "criterion": "mean/clone effects follow expected direction in non-null comparisons",
            "observed": (
                f"GSE125472 {int(expected_gse125['direction_matches_expected'].sum())}/{len(expected_gse125)}; "
                f"TCF7L2 KO {int(expected_clone['direction_matches_expected'].sum())}/{len(expected_clone)}"
            ),
            "pass": bool(
                expected_gse125["direction_matches_expected"].fillna(False).all()
                and expected_clone["direction_matches_expected"].fillna(False).all()
            ),
        },
        {
            "gate": "leave_one_gene_out",
            "criterion": "all pooled effects positive and retain >= 75% of 10-gene effect",
            "observed": float(loo["external_retained_fraction_vs_10_gene"].min()),
            "pass": bool(
                (loo["external_pooled_effect_sd"] > 0).all()
                and (loo["external_retained_fraction_vs_10_gene"] >= 0.75).all()
            ),
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    panel = reduced_signature()
    write_tsv(panel, "reduced_panel_definition.tsv")

    chen_scores, chen_performance, chen_paired, chen_expr = chen_analysis(panel)
    validation_meta, _ = discovery.load_chen_pseudobulk("validation")

    external_scores, external_tests, pooled, external_corr, external_data = external_analysis(panel)
    ffpe_scores, ffpe_test, ffpe_corr, ffpe_gene_expression = ffpe_analysis(panel)
    perturb_scores, perturb_summary, clone_effects = perturbation_analysis(panel)

    coverage_rows = []
    for cohort, (_, expression) in external_data.items():
        for gene in panel["gene"]:
            coverage_rows.append(
                {
                    "scope": "external_sporadic",
                    "cohort": cohort,
                    "gene": gene,
                    "present": gene in expression.columns,
                }
            )
    for gene in panel["gene"]:
        coverage_rows.append(
            {
                "scope": "FFPE",
                "cohort": "GSE117606",
                "gene": gene,
                "present": gene in ffpe_gene_expression.columns,
            }
        )
    panel_coverage = pd.DataFrame(coverage_rows).merge(
        panel[["gene", "panel_arm"]], on="gene", how="left", validate="many_to_one"
    )

    chen_corr = pd.DataFrame(
        [
            spearman_row(
                "Chen",
                dataset,
                frame["reduced_score"],
                frame["full_100_gene_score"],
            )
            for dataset, frame in chen_scores.groupby("dataset", sort=False)
        ]
    )
    correlations = pd.concat([chen_corr, external_corr, ffpe_corr], ignore_index=True)

    loo = leave_one_gene_out(
        panel,
        external_data,
        validation_meta,
        chen_expr["validation"],
        ffpe_scores,
        ffpe_gene_expression,
        float(pooled["adenoma_coef_sd"].iloc[0]),
    )
    proteins = protein_anchor_table()
    gates = gate_summary(
        chen_performance,
        correlations,
        pooled,
        ffpe_test,
        perturb_summary,
        clone_effects,
        loo,
        external_scores,
    )

    write_tsv(chen_scores, "chen_reduced_sample_scores.tsv")
    write_tsv(chen_performance, "chen_reduced_performance.tsv")
    write_tsv(chen_paired, "chen_reduced_paired_tests.tsv")
    write_tsv(external_scores, "external_reduced_sample_scores.tsv")
    write_tsv(external_tests, "external_reduced_cohort_tests.tsv")
    write_tsv(pooled, "external_reduced_pooled_model.tsv")
    write_tsv(ffpe_scores, "ffpe_reduced_sample_scores.tsv")
    write_tsv(ffpe_test, "ffpe_reduced_paired_test.tsv")
    write_tsv(correlations, "reduced_vs_100_gene_concordance.tsv")
    write_tsv(perturb_scores, "reduced_perturbation_sample_scores.tsv")
    write_tsv(perturb_summary, "reduced_apc_organoid_effects.tsv")
    write_tsv(clone_effects, "reduced_tcf7l2_clone_effects.tsv")
    write_tsv(loo, "reduced_panel_leave_one_gene_out.tsv")
    write_tsv(proteins, "protein_anchor_evidence.tsv")
    write_tsv(panel_coverage, "reduced_panel_platform_coverage.tsv")
    write_tsv(gates, "analysis_gate_summary.tsv")

    input_paths = [
        ROOT / "results" / "route_signature_locked" / "discovery_locked_signature_genes.tsv",
        ROOT / "results" / "external_sporadic_adenoma_validation" / "sample_scores.tsv",
        ROOT
        / "results"
        / "gse117606_paired_route_validation"
        / "sample_scores_all_cohort_scaling.tsv",
        ROOT
        / "results"
        / "public_adenoma_protein_triangulation"
        / "candidate_public_protein_evidence_matrix.tsv",
    ]
    manifest = {
        "analysis": "post hoc fixed 10-gene translational reduction",
        "primary_research_readout": "100-gene discovery-locked programme",
        "reduced_panel_status": "exploratory; no diagnostic cut-point or gene swapping",
        "up_genes": UP_GENES,
        "down_genes": DOWN_GENES,
        "random_seed": SEED,
        "input_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in input_paths},
        "all_gates_pass": bool(gates["pass"].all()),
    }
    (OUT_DIR / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    heldout = chen_performance.set_index("dataset").loc["validation"]
    ffpe_row = ffpe_test.iloc[0]
    summary = [
        "# Fixed 10-gene translational reduction",
        "",
        "The 100-gene discovery-locked programme remains primary. The 10-gene panel is a post hoc, fixed exploratory reduction; no failing gene was exchanged.",
        "",
        f"- Held-out Chen: AUC {heldout['auc']:.3f}; P = {heldout['p_value']:.3g}.",
        f"- Five-cohort pooled effect: {pooled['adenoma_coef_sd'].iloc[0]:.3f} SD (95% CI {pooled['ci_low'].iloc[0]:.3f} to {pooled['ci_high'].iloc[0]:.3f}); retained fraction {pooled['retained_effect_fraction'].iloc[0]:.3f}.",
        f"- GSE117606 FFPE: {int(ffpe_row['n_pairs'])} pairs; median change {ffpe_row['median_paired_difference']:.3f}; {ffpe_row['paired_positive_fraction']:.1%} positive; P = {ffpe_row['p_paired_wilcoxon']:.3g}.",
        f"- Minimum external score concordance with the 100-gene programme: {external_corr['spearman_rho_reduced_vs_100_gene'].min():.3f}.",
        f"- Leave-one-gene-out minimum retained pooled effect: {loo['external_retained_fraction_vs_10_gene'].min():.3f}.",
        f"- Formal analysis gates: {int(gates['pass'].sum())}/{len(gates)} passed.",
        "",
        "Interpretation: the reduced panel supports assay-development feasibility, not a validated clinical assay, diagnostic classifier, outcome predictor or surrogate endpoint.",
    ]
    (OUT_DIR / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    if not gates["pass"].all():
        failed = gates.loc[~gates["pass"], "gate"].tolist()
        raise RuntimeError(f"Reduced panel failed prespecified gates: {failed}")


if __name__ == "__main__":
    main()
