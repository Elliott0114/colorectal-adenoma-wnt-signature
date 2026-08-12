#!/usr/bin/env python3
"""Patient-paired validation of the locked route in public GSE117606.

The locked 50-up/50-down signature is transferred without re-selection. The
primary comparison is conventional adenoma versus patient-matched adjacent
mucosa; sessile serrated adenomas and carcinomas are retained for context only.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

from external_sporadic_adenoma_validation import (
    PRIMARY_SIZE,
    PROLIFERATION_CONTROL,
    SEED,
    SIGNATURE_SIZES,
    bootstrap_primary,
    cohort_comparison,
    score_nested_signatures,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "data_sources" / "GSE117606"
OUT_DIR = ROOT / "results" / "gse117606_paired_route_validation"
SERIES_PATH = SOURCE_DIR / "GSE117606_series_matrix.txt.gz"
GENE_INFO_PATH = SOURCE_DIR / "Homo_sapiens.gene_info.gz"
SIGNATURE_PATH = ROOT / "results" / "route_signature_locked" / "discovery_locked_signature_genes.tsv"
PROTEIN_EVIDENCE_PATH = (
    ROOT
    / "results"
    / "public_adenoma_protein_triangulation"
    / "candidate_public_protein_evidence_matrix.tsv"
)

SOURCE_URLS = {
    "GSE117606_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117606/"
        "matrix/GSE117606_series_matrix.txt.gz"
    ),
    "GSE117606_GEO_jreumers_CRC_concurrent.rds.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117606/suppl/"
        "GSE117606_GEO_jreumers_CRC_concurrent.rds.gz"
    ),
    "Homo_sapiens.gene_info.gz": (
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/"
        "Homo_sapiens.gene_info.gz"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def parse_metadata() -> pd.DataFrame:
    metadata: dict[str, list[list[str]]] = {}
    with gzip.open(SERIES_PATH, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                continue
            parts = next(csv.reader([line], delimiter="\t"))
            metadata.setdefault(parts[0], []).append(parts[1:])

    samples = pd.DataFrame(
        {
            "sample_id": metadata["!Sample_geo_accession"][0],
            "sample_title": metadata["!Sample_title"][0],
            "source_name": metadata["!Sample_source_name_ch1"][0],
        }
    )
    for values in metadata.get("!Sample_characteristics_ch1", []):
        pairs = [value.split(":", 1) if ":" in value else ["characteristic", value] for value in values]
        key = pd.Series([clean_key(pair[0]) for pair in pairs]).mode().iloc[0]
        samples[key] = [pair[1].strip() if len(pair) > 1 else "" for pair in pairs]

    samples["patient_id"] = samples["patient_id"].astype(str)
    samples["age_years"] = pd.to_numeric(samples["age"], errors="coerce")
    samples["sex"] = samples["gender"].str.lower()
    samples["tissue_group"] = samples["tissue"].str.lower().map(
        {
            "adjacent mucosa": "normal",
            "adenoma": "adenoma",
            "ssa": "ssa",
            "tumor": "crc",
        }
    )
    if samples["patient_id"].isna().any() or samples["tissue_group"].isna().any():
        raise ValueError("GSE117606 patient or tissue labels could not be parsed")
    return samples


def parse_expression() -> pd.DataFrame:
    matrix = pd.read_csv(
        SERIES_PATH,
        sep="\t",
        comment="!",
        index_col=0,
        compression="gzip",
    )
    matrix.index = matrix.index.astype(str)
    expression = matrix.apply(pd.to_numeric, errors="coerce").T
    expression.index.name = "sample_id"
    return expression


def map_symbols_to_features(genes: set[str], expression: pd.DataFrame) -> pd.DataFrame:
    info = pd.read_csv(GENE_INFO_PATH, sep="\t", dtype=str, low_memory=False)
    direct = info.loc[info["Symbol"].isin(genes), ["GeneID", "Symbol"]].rename(
        columns={"Symbol": "gene"}
    )
    direct["mapping_method"] = "current_ncbi_symbol"
    missing = genes - set(direct["gene"])
    alias_rows = []
    if missing:
        synonym_sets = info["Synonyms"].fillna("").map(lambda value: set(value.split("|")))
        for gene in sorted(missing):
            matches = info.loc[synonym_sets.map(lambda values: gene in values), ["GeneID", "Symbol"]]
            if len(matches) == 1:
                alias_rows.append(
                    {
                        "GeneID": matches.iloc[0]["GeneID"],
                        "gene": gene,
                        "mapping_method": f"unique_ncbi_synonym_of_{matches.iloc[0]['Symbol']}",
                    }
                )
    mapping = pd.concat([direct, pd.DataFrame(alias_rows)], ignore_index=True)
    mapping["feature_id"] = mapping["GeneID"].astype(str) + "_at"
    mapping["feature_present"] = mapping["feature_id"].isin(expression.columns)
    mapping = mapping.sort_values(["gene", "mapping_method", "GeneID"]).drop_duplicates("gene")
    return mapping


def select_gene_expression(
    expression: pd.DataFrame, mapping: pd.DataFrame
) -> pd.DataFrame:
    selected = mapping.loc[mapping["feature_present"]].copy()
    output = expression[selected["feature_id"]].copy()
    output.columns = selected["gene"].tolist()
    return output


def paired_gene_tests(
    meta: pd.DataFrame,
    gene_expression: pd.DataFrame,
    candidates: pd.DataFrame,
) -> pd.DataFrame:
    data = meta[["sample_id", "patient_id", "tissue_group"]].merge(
        gene_expression.reset_index(), on="sample_id", how="left"
    )
    rows = []
    for row in candidates.itertuples(index=False):
        if row.gene not in gene_expression.columns:
            rows.append({"gene": row.gene, "n_pairs": 0})
            continue
        wide = (
            data.loc[data["tissue_group"].isin(["normal", "adenoma"])]
            .pivot_table(index="patient_id", columns="tissue_group", values=row.gene, aggfunc="median")
            .reindex(columns=["normal", "adenoma"])
            .dropna()
        )
        delta = wide["adenoma"] - wide["normal"]
        if len(delta) >= 2 and not np.allclose(delta, 0):
            wilcoxon_p = float(stats.wilcoxon(delta, alternative="two-sided").pvalue)
            paired_t_p = float(stats.ttest_rel(wide["adenoma"], wide["normal"]).pvalue)
            se = float(stats.sem(delta))
            ci = stats.t.interval(0.95, len(delta) - 1, loc=float(delta.mean()), scale=se)
        else:
            wilcoxon_p = paired_t_p = np.nan
            ci = (np.nan, np.nan)
        sd_delta = float(delta.std(ddof=1)) if len(delta) >= 2 else np.nan
        rows.append(
            {
                "gene": row.gene,
                "n_pairs": len(delta),
                "median_normal": float(wide["normal"].median()) if len(wide) else np.nan,
                "median_adenoma": float(wide["adenoma"].median()) if len(wide) else np.nan,
                "median_paired_delta": float(delta.median()) if len(delta) else np.nan,
                "mean_paired_delta": float(delta.mean()) if len(delta) else np.nan,
                "mean_delta_ci_low": float(ci[0]),
                "mean_delta_ci_high": float(ci[1]),
                "paired_cohen_dz": float(delta.mean() / sd_delta) if sd_delta > 0 else np.nan,
                "paired_positive_fraction": float((delta > 0).mean()) if len(delta) else np.nan,
                "p_paired_wilcoxon": wilcoxon_p,
                "p_paired_t": paired_t_p,
                "expected_direction": row.expected_direction,
                "gse117606_direction_match": (
                    bool(np.sign(delta.median()) == row.expected_direction) if len(delta) else False
                ),
            }
        )
    output = candidates.merge(pd.DataFrame(rows), on="gene", how="left")
    keep = output["p_paired_wilcoxon"].notna()
    output["q_paired_wilcoxon_bh_candidates"] = np.nan
    if keep.any():
        output.loc[keep, "q_paired_wilcoxon_bh_candidates"] = multipletests(
            output.loc[keep, "p_paired_wilcoxon"], method="fdr_bh"
        )[1]
    return output.sort_values(
        ["q_paired_wilcoxon_bh_candidates", "direction_stability"],
        ascending=[True, False],
        na_position="last",
    )


def proliferation_adjusted_paired_delta(scores: pd.DataFrame) -> pd.DataFrame:
    subset = scores.loc[scores["tissue_group"].isin(["normal", "adenoma"])].copy()
    route = subset.pivot_table(
        index="patient_id", columns="tissue_group", values=f"route_score_k{PRIMARY_SIZE}", aggfunc="median"
    ).dropna()
    proliferation = subset.pivot_table(
        index="patient_id", columns="tissue_group", values="score__proliferation_control", aggfunc="median"
    ).dropna()
    patients = route.index.intersection(proliferation.index)
    route_delta = route.loc[patients, "adenoma"] - route.loc[patients, "normal"]
    proliferation_delta = proliferation.loc[patients, "adenoma"] - proliferation.loc[patients, "normal"]
    design = sm.add_constant(
        pd.DataFrame({"proliferation_delta": proliferation_delta}), has_constant="add"
    )
    fit = sm.OLS(route_delta, design).fit(cov_type="HC3")
    interval = fit.conf_int().loc["const"]
    return pd.DataFrame(
        [
            {
                "model": "paired_delta_route_on_paired_delta_proliferation",
                "n_pairs": len(patients),
                "adjusted_mean_adenoma_minus_normal": float(fit.params["const"]),
                "hc3_ci_low": float(interval.iloc[0]),
                "hc3_ci_high": float(interval.iloc[1]),
                "hc3_p_value": float(fit.pvalues["const"]),
                "proliferation_delta_coefficient": float(fit.params["proliferation_delta"]),
                "proliferation_delta_p_value": float(fit.pvalues["proliferation_delta"]),
                "condition_number": float(np.linalg.cond(design.to_numpy(dtype=float))),
            }
        ]
    )


def integrated_priority(gene_tests: pd.DataFrame, protein: pd.DataFrame) -> pd.DataFrame:
    protein_columns = [
        "gene",
        "public_protein_evidence_tier",
        "age_sex_adjusted_log2_effect",
        "age_sex_adjusted_ci_low",
        "age_sex_adjusted_ci_high",
        "age_sex_adjusted_q_bh_candidates",
        "pxd017269_detection_fraction",
        "detected_in_nine_patient_dvp",
    ]
    output = gene_tests.merge(protein[protein_columns], on="gene", how="left")
    output["paired_transcript_fdr_10pct"] = output["q_paired_wilcoxon_bh_candidates"].le(0.1)
    output["public_protein_directional_support"] = output[
        "public_protein_evidence_tier"
    ].str.startswith(("A_", "B_"), na=False)
    output["ffpe_detection_gate_80pct"] = output["pxd017269_detection_fraction"].ge(0.8)
    strict = (
        output["locked_signature_member"].astype(bool)
        & output["gse117606_direction_match"].astype(bool)
        & output["paired_transcript_fdr_10pct"].astype(bool)
        & output["public_protein_directional_support"].astype(bool)
        & output["ffpe_detection_gate_80pct"].astype(bool)
    )
    partial = (
        output["locked_signature_member"].astype(bool)
        & output["gse117606_direction_match"].astype(bool)
        & output["ffpe_detection_gate_80pct"].astype(bool)
    )
    output["wet_lab_evidence_class"] = np.select(
        [strict, partial, output["locked_signature_member"].astype(bool)],
        [
            "P1_locked_paired_transcript_and_directional_protein_support",
            "P2_locked_paired_transcript_with_partial_protein_support",
            "P3_locked_but_incomplete_or_conflicting_replication",
        ],
        default="context_marker_not_locked",
    )
    order = {
        "P1_locked_paired_transcript_and_directional_protein_support": 0,
        "P2_locked_paired_transcript_with_partial_protein_support": 1,
        "P3_locked_but_incomplete_or_conflicting_replication": 2,
        "context_marker_not_locked": 3,
    }
    output["sort_order"] = output["wet_lab_evidence_class"].map(order)
    return output.sort_values(
        ["sort_order", "q_paired_wilcoxon_bh_candidates", "direction_stability"],
        ascending=[True, True, False],
        na_position="last",
    ).drop(columns="sort_order")


def source_manifest() -> pd.DataFrame:
    rows = []
    for path in sorted(SOURCE_DIR.glob("*")):
        if path.is_file():
            rows.append(
                {
                    "file": str(path.relative_to(ROOT)),
                    "analysis_role": (
                        "primary_expression_and_metadata"
                        if path.name == SERIES_PATH.name
                        else "gene_identifier_mapping"
                        if path.name == GENE_INFO_PATH.name
                        else "official_processed_object_retained_not_analyzed"
                    ),
                    "source_url": SOURCE_URLS.get(path.name, ""),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = parse_metadata()
    expression = parse_expression()
    if list(expression.index) != meta["sample_id"].tolist():
        raise ValueError("GSE117606 expression and metadata sample order differ")

    signature = pd.read_csv(SIGNATURE_PATH, sep="\t")
    protein = pd.read_csv(PROTEIN_EVIDENCE_PATH, sep="\t")
    candidates = protein[
        [
            "gene",
            "panel_role",
            "locked_signature_member",
            "selected",
            "discovery_effect_adenoma_minus_normal",
            "bootstrap_ci_low",
            "bootstrap_ci_high",
            "direction_stability",
            "expected_direction",
        ]
    ].copy()
    wanted = set(signature["gene"]) | set(PROLIFERATION_CONTROL) | set(candidates["gene"])
    mapping = map_symbols_to_features(wanted, expression)
    gene_expression = select_gene_expression(expression, mapping)

    scores_all, coverage_all = score_nested_signatures(
        meta, gene_expression, signature, "GSE117606"
    )
    scores_primary_set, coverage_primary_set = score_nested_signatures(
        meta,
        gene_expression,
        signature,
        "GSE117606",
        normalization_groups={"normal", "adenoma"},
    )
    coverage = pd.DataFrame(
        [
            {**row, "normalization_scope": "all_cohort_samples"}
            for row in coverage_all
        ]
        + [
            {**row, "normalization_scope": "normal_and_adenoma_only"}
            for row in coverage_primary_set
        ]
    )

    comparisons = []
    for scope, scores in [
        ("all_cohort_samples", scores_all),
        ("normal_and_adenoma_only", scores_primary_set),
    ]:
        for size in SIGNATURE_SIZES:
            comparison = cohort_comparison(
                scores,
                "GSE117606",
                size,
                "adenoma",
                "normal",
                paired=True,
                comparison="conventional_adenoma_vs_adjacent_mucosa",
            )
            comparison["normalization_scope"] = scope
            comparisons.append(comparison)
    comparisons = pd.DataFrame(comparisons)
    bootstrap = pd.DataFrame(
        [bootstrap_primary(scores_all, "GSE117606", paired=True, rng=np.random.default_rng(SEED))]
    )
    proliferation_model = proliferation_adjusted_paired_delta(scores_all)
    gene_tests = paired_gene_tests(meta, gene_expression, candidates)
    priority = integrated_priority(gene_tests, protein)

    meta.to_csv(OUT_DIR / "sample_metadata.tsv", sep="\t", index=False)
    mapping.to_csv(OUT_DIR / "gene_to_entrez_mapping.tsv", sep="\t", index=False)
    coverage.to_csv(OUT_DIR / "locked_signature_coverage.tsv", sep="\t", index=False)
    scores_all.to_csv(OUT_DIR / "sample_scores_all_cohort_scaling.tsv", sep="\t", index=False)
    scores_primary_set.to_csv(
        OUT_DIR / "sample_scores_normal_adenoma_scaling.tsv", sep="\t", index=False
    )
    comparisons.to_csv(OUT_DIR / "paired_route_tests.tsv", sep="\t", index=False)
    bootstrap.to_csv(OUT_DIR / "paired_route_bootstrap.tsv", sep="\t", index=False)
    proliferation_model.to_csv(
        OUT_DIR / "paired_proliferation_adjusted_model.tsv", sep="\t", index=False
    )
    gene_tests.to_csv(OUT_DIR / "candidate_gene_paired_tests.tsv", sep="\t", index=False)
    priority.to_csv(OUT_DIR / "candidate_integrated_wetlab_priority.tsv", sep="\t", index=False)
    source_manifest().to_csv(OUT_DIR / "source_file_manifest.tsv", sep="\t", index=False)
    (OUT_DIR / "analysis_manifest.json").write_text(
        json.dumps(
            {
                "analysis": "GSE117606 patient-paired locked-route validation",
                "signature_reselected": False,
                "primary_comparison": "conventional adenoma vs patient-matched adjacent mucosa",
                "python_packages": {
                    "numpy": np.__version__,
                    "pandas": pd.__version__,
                    "scipy": scipy.__version__,
                    "statsmodels": statsmodels.__version__,
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    primary = comparisons.loc[
        comparisons["signature_size_per_direction"].eq(PRIMARY_SIZE)
        & comparisons["normalization_scope"].eq("all_cohort_samples")
    ].iloc[0]
    tissue_counts = meta["tissue_group"].value_counts()
    paired_patients = (
        meta.loc[meta["tissue_group"].isin(["normal", "adenoma"])]
        .groupby("patient_id")["tissue_group"]
        .nunique()
        .eq(2)
        .sum()
    )
    panel = priority.loc[priority["panel_role"].ne("locked_signature_candidate")]
    lines = [
        "GSE117606 patient-paired locked-route validation",
        "=================================================",
        f"Samples: {len(meta)} from {meta['patient_id'].nunique()} deposited patient IDs; conventional adenoma={tissue_counts.get('adenoma', 0)}, adjacent mucosa={tissue_counts.get('normal', 0)}, SSA={tissue_counts.get('ssa', 0)}, CRC={tissue_counts.get('crc', 0)}.",
        "The GEO series summary reports 70 patients, whereas the deposited sample titles resolve to 75 unique patient IDs; this discrepancy is retained as a source-data audit finding and does not alter the 51 complete-pair primary analysis.",
        f"Complete conventional adenoma-adjacent mucosa pairs: {paired_patients}.",
        f"Locked k50 route paired median delta: {primary['median_paired_difference']:.3f}; Wilcoxon p={primary['p_paired_wilcoxon']:.3g}; positive-pair fraction={primary['paired_positive_fraction']:.3f}.",
        f"Proliferation-adjusted paired mean delta: {proliferation_model.iloc[0]['adjusted_mean_adenoma_minus_normal']:.3f} (HC3 95% CI {proliferation_model.iloc[0]['hc3_ci_low']:.3f} to {proliferation_model.iloc[0]['hc3_ci_high']:.3f}).",
        "The locked signature was transferred without re-selection; SSA and CRC were not included in the primary comparison.",
        "",
        "Current IHC panel and prespecified alternative",
    ]
    for row in panel.itertuples(index=False):
        lines.append(
            f"- {row.gene}: paired median transcript delta={row.median_paired_delta:.3f}; "
            f"q={row.q_paired_wilcoxon_bh_candidates:.3g}; direction match={row.gse117606_direction_match}; "
            f"integrated class={row.wet_lab_evidence_class}."
        )
    lines.extend(
        [
            "",
            "Interpretation boundary",
            "The integrated class is an exploratory wet-lab prioritization aid, not a change to the locked route or a clinical biomarker claim.",
        ]
    )
    (OUT_DIR / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
