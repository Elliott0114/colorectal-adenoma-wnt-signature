#!/usr/bin/env python3
"""Validate the frozen adenoma route in additional public transcriptomic cohorts.

The locked signature is transferred without gene re-selection. Array probes or
transcript rows are selected by variance across all deposited samples, without
using disease labels. Cohorts with material design confounding are retained as
descriptive stress tests and are not counted as primary inferential evidence.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from external_sporadic_adenoma_validation import (
    PRIMARY_SIZE,
    PROLIFERATION_CONTROL,
    SIGNATURE_SIZES,
    bootstrap_primary,
    choose_features_without_labels,
    cohort_comparison,
    parse_series_matrix,
    score_nested_signatures,
    split_symbols,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "data_sources" / "additional_adenoma_validation"
OUT_DIR = ROOT / "results" / "expanded_public_adenoma_validation"
SIGNATURE_PATH = (
    ROOT / "results" / "route_signature_locked" / "discovery_locked_signature_genes.tsv"
)
CANDIDATE_PATH = (
    ROOT
    / "results"
    / "public_adenoma_protein_triangulation"
    / "candidate_inventory.tsv"
)
GPL570_PATH = ROOT / "data_sources" / "GEO_bulk_recurrence" / "GPL570" / "GPL570.annot.gz"
GPL6480_PATH = (
    ROOT / "data_sources" / "GEO_sporadic_adenoma_validation" / "GPL6480.annot.gz"
)

N_BOOTSTRAPS = 5_000
SEED = 20260714
CURRENT_PANEL_ROLES = {
    "current_primary_up_arm",
    "current_primary_down_arm",
    "current_localization_readout",
    "current_context_nuclear_marker",
    "proliferation_covariate",
    "prespecified_down_arm_alternative",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_matrix(path: Path, accession: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Reuse the established parser while pointing it at an explicit directory."""
    import external_sporadic_adenoma_validation as established

    previous = established.SOURCE_DIR
    established.SOURCE_DIR = path.parent
    try:
        return parse_series_matrix(accession)
    finally:
        established.SOURCE_DIR = previous


def geo_annotation_mapping(path: Path) -> pd.DataFrame:
    rows: list[list[str]] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if parts and parts[0] == "!platform_table_begin":
                in_table = True
                continue
            if parts and parts[0] == "!platform_table_end":
                break
            if in_table:
                rows.append(parts)
    if not rows:
        raise ValueError(f"No platform table found in {path}")
    table = pd.DataFrame(rows[1:], columns=rows[0])
    output = []
    for feature_id, symbols in table[["ID", "Gene symbol"]].itertuples(index=False):
        for gene in split_symbols(symbols):
            output.append({"feature_id": str(feature_id), "gene": gene})
    return pd.DataFrame(output).drop_duplicates()


def direct_gene_selection(
    expression: pd.DataFrame,
    cohort: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = pd.DataFrame(
        {
            "cohort": cohort,
            "feature_id": expression.columns.astype(str),
            "gene": expression.columns.astype(str),
            "variance_all_samples": expression.var(axis=0, ddof=1).to_numpy(),
            "selection_uses_group_labels": False,
        }
    )
    return expression, selected


def read_gse100179(
    wanted: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    series_path = SOURCE_DIR / "GSE100179" / "GSE100179_series_matrix.txt.gz"
    supplement_path = (
        SOURCE_DIR
        / "GSE100179"
        / "GSE100179_2017-07-21-annotated-RMA-SKETCH.RMA-GENE-FULL-Group1.txt.gz"
    )
    meta, _, _ = parse_matrix(series_path, "GSE100179")
    header = pd.read_csv(supplement_path, sep="\t", nrows=0).columns.tolist()
    sample_columns = header[1:61]
    usecols = [header[0], *sample_columns, "Gene Symbol"]
    table = pd.read_csv(supplement_path, sep="\t", usecols=usecols, low_memory=False)
    table = table.rename(columns={header[0]: "feature_id"})

    expression = table[sample_columns].apply(pd.to_numeric, errors="coerce").T
    if expression.shape[0] != len(meta):
        raise ValueError("GSE100179 supplement/sample metadata size mismatch")
    expression.index = meta["sample_id"]
    expression.columns = table["feature_id"].astype(str)

    mapping_rows = []
    for feature_id, symbols in table[["feature_id", "Gene Symbol"]].itertuples(index=False):
        for gene in split_symbols(symbols):
            mapping_rows.append({"feature_id": str(feature_id), "gene": gene})
    mapping = pd.DataFrame(mapping_rows).drop_duplicates()
    gene_expression, selected = choose_features_without_labels(
        expression, mapping, wanted, "GSE100179"
    )

    tissue = meta["tissue"].str.lower()
    meta["tissue_group"] = np.select(
        [tissue.str.contains("adenoma"), tissue.str.contains("healthy"), tissue.str.contains("adenocarcinoma")],
        ["adenoma", "normal", "crc"],
        default="other",
    )
    meta["patient_id"] = meta["sample_id"]
    meta["pathologic_grade"] = meta["tissue"]
    return meta, gene_expression, selected


def read_gse37364(
    wanted: set[str], mapping: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = SOURCE_DIR / "GSE37364" / "GSE37364_series_matrix.txt.gz"
    meta, expression, _ = parse_matrix(path, "GSE37364")
    gene_expression, selected = choose_features_without_labels(
        expression, mapping, wanted, "GSE37364"
    )
    title = meta["sample_title"]
    meta["tissue_group"] = np.select(
        [title.str.startswith(("AH", "AL")), title.str.startswith("N"), title.str.startswith("CRC")],
        ["adenoma", "normal", "crc"],
        default="other",
    )
    meta["patient_id"] = meta["sample_id"]
    meta["pathologic_grade"] = np.select(
        [title.str.startswith("AH"), title.str.startswith("AL")],
        ["high_grade", "low_grade"],
        default=meta["tissue_group"],
    )
    return meta, gene_expression, selected


def read_gse164541(
    wanted: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = SOURCE_DIR / "GSE164541" / "GSE164541_ANT_count.csv.gz"
    table = pd.read_csv(path)
    sample_columns = table.columns[3:].tolist()
    counts = table[["gene_name", *sample_columns]].copy()
    counts["gene_name"] = counts["gene_name"].astype(str).str.strip()
    counts = counts.loc[counts["gene_name"].ne("") & counts["gene_name"].ne("nan")]
    counts = counts.groupby("gene_name", as_index=True)[sample_columns].sum()
    library_sizes = counts.sum(axis=0)
    log_cpm = np.log2(counts.div(library_sizes, axis=1) * 1_000_000 + 0.5).T
    expression = log_cpm.loc[:, log_cpm.columns.intersection(sorted(wanted))]
    expression.index.name = "sample_id"
    expression, selected = direct_gene_selection(expression, "GSE164541")

    meta = pd.DataFrame({"sample_id": sample_columns})
    meta["sample_title"] = meta["sample_id"]
    meta["source_name"] = "paired normal-adenoma-tumor RNA-seq"
    prefix = meta["sample_id"].str.extract(r"P([ANT])", expand=False)
    meta["tissue_group"] = prefix.map({"A": "adenoma", "N": "normal", "T": "crc"})
    meta["patient_id"] = meta["sample_id"].str.extract(r"(\d+)$", expand=False)
    meta["pathologic_grade"] = meta["tissue_group"]
    if meta["patient_id"].isna().any():
        raise ValueError("GSE164541 patient identifiers could not be parsed")
    return meta, expression, selected


def read_gse20916(
    wanted: set[str], mapping: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = SOURCE_DIR / "GSE20916" / "GSE20916_series_matrix.txt.gz"
    meta, expression, _ = parse_matrix(path, "GSE20916")
    gene_expression, selected = choose_features_without_labels(
        expression, mapping, wanted, "GSE20916"
    )
    title = meta["sample_title"].str.lower()
    meta["tissue_group"] = np.select(
        [
            title.str.contains("distant_normal"),
            title.str.contains("adenoma"),
            title.str.contains("adenocarcinoma|carcinoma", regex=True),
            title.str.contains("normal_colon"),
        ],
        ["distant_normal", "adenoma", "crc", "normal"],
        default="other",
    )
    meta["stratum"] = np.where(title.str.contains(r"\(micro\)"), "micro", "macro")
    meta["cell_compartment"] = np.where(
        title.str.startswith("colonic_crypt_epithelial_cells"), "crypt_epithelium", "mucosa_or_bulk"
    )
    meta["patient_id"] = meta["sample_id"]
    meta["pathologic_grade"] = meta["tissue_group"]
    return meta, gene_expression, selected


def read_gse226739(
    wanted: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    archive = SOURCE_DIR / "GSE226739" / "PMC11574242_supplementaryFiles.zip"
    member = "41598_2024_70455_MOESM6_ESM.xlsx"
    with zipfile.ZipFile(archive) as zf:
        content = zf.read(member)
    cra = pd.read_excel(io.BytesIO(content), sheet_name="CRA_vs_HC", header=16)
    crc = pd.read_excel(io.BytesIO(content), sheet_name="CRC_vs_HC", header=16)
    cra.columns = [str(column) for column in cra.columns]
    crc.columns = [str(column) for column in crc.columns]
    cra_sample_columns = cra.columns[11:].tolist()
    crc_sample_columns = crc.columns[11:].tolist()
    normal_columns = [column for column in cra_sample_columns if column in set(crc_sample_columns)]
    adenoma_columns = [column for column in cra_sample_columns if column not in set(normal_columns)]
    if len(normal_columns) != 5 or len(adenoma_columns) != 5:
        raise ValueError("GSE226739 supplementary comparison does not resolve to 5+5 samples")

    table = cra.loc[cra["Gene_Name"].notna()].copy()
    table["feature_id"] = table["Track_id"].astype(str)
    expression = np.log2(
        table[cra_sample_columns].apply(pd.to_numeric, errors="coerce").clip(lower=0) + 1
    ).T
    expression.columns = table["feature_id"]
    mapping = pd.DataFrame(
        {
            "feature_id": table["feature_id"].astype(str),
            "gene": table["Gene_Name"].astype(str).str.strip(),
        }
    ).drop_duplicates()
    gene_expression, selected = choose_features_without_labels(
        expression, mapping, wanted, "GSE226739"
    )
    meta = pd.DataFrame({"sample_id": cra_sample_columns})
    meta["sample_title"] = meta["sample_id"]
    meta["source_name"] = "author-deposited FPKM comparison matrix"
    meta["tissue_group"] = np.where(meta["sample_id"].isin(adenoma_columns), "adenoma", "normal")
    meta["patient_id"] = meta["sample_id"]
    meta["pathologic_grade"] = meta["tissue_group"]
    return meta, gene_expression, selected


def read_gse71187(
    wanted: set[str], mapping: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = SOURCE_DIR / "GSE71187" / "GSE71187_series_matrix.txt.gz"
    meta, expression, _ = parse_matrix(path, "GSE71187")
    gene_expression, selected = choose_features_without_labels(
        expression, mapping, wanted, "GSE71187"
    )
    source = meta["source_name"].str.lower()
    meta["tissue_group"] = np.select(
        [source.str.contains("adenoma"), source.str.contains("normal colorectal")],
        ["adenoma", "normal"],
        default="other",
    )
    meta["patient_id"] = meta["sample_id"]
    meta["pathologic_grade"] = meta["tissue_group"]
    return meta, gene_expression, selected


def candidate_gene_tests(
    cohort: str,
    meta: pd.DataFrame,
    expression: pd.DataFrame,
    candidates: pd.DataFrame,
    paired: bool,
) -> pd.DataFrame:
    joined = meta.set_index("sample_id").join(expression, how="inner")
    joined = joined.loc[joined["tissue_group"].isin(["normal", "adenoma"])]
    rows = []
    for candidate in candidates.itertuples(index=False):
        gene = candidate.gene
        if gene not in expression.columns:
            continue
        normal = joined.loc[joined["tissue_group"].eq("normal"), gene].dropna()
        adenoma = joined.loc[joined["tissue_group"].eq("adenoma"), gene].dropna()
        if len(normal) < 2 or len(adenoma) < 2:
            continue
        mw = stats.mannwhitneyu(adenoma, normal, alternative="two-sided")
        row = {
            "cohort": cohort,
            "gene": gene,
            "panel_role": candidate.panel_role,
            "expected_direction": candidate.expected_direction,
            "paired": paired,
            "n_normal": len(normal),
            "n_adenoma": len(adenoma),
            "median_normal": float(normal.median()),
            "median_adenoma": float(adenoma.median()),
            "median_effect_adenoma_minus_normal": float(adenoma.median() - normal.median()),
            "auc_adenoma_greater": float(mw.statistic / (len(adenoma) * len(normal))),
            "p_mannwhitney": float(mw.pvalue),
            "n_pairs": np.nan,
            "median_paired_effect": np.nan,
            "paired_positive_fraction": np.nan,
            "p_paired_wilcoxon": np.nan,
            "primary_p_value": float(mw.pvalue),
        }
        effect_for_direction = row["median_effect_adenoma_minus_normal"]
        if paired:
            wide = joined.pivot_table(
                index="patient_id", columns="tissue_group", values=gene, aggfunc="median"
            )
            pair = wide[["adenoma", "normal"]].dropna()
            delta = pair["adenoma"] - pair["normal"]
            if len(delta) >= 2:
                wilcoxon_p = (
                    1.0
                    if np.allclose(delta.to_numpy(dtype=float), 0.0)
                    else float(
                        stats.wilcoxon(
                            delta, alternative="two-sided", zero_method="wilcox"
                        ).pvalue
                    )
                )
                row.update(
                    {
                        "n_pairs": len(delta),
                        "median_paired_effect": float(delta.median()),
                        "paired_positive_fraction": float((delta > 0).mean()),
                        "p_paired_wilcoxon": wilcoxon_p,
                        "primary_p_value": wilcoxon_p,
                    }
                )
                effect_for_direction = float(delta.median())
        row["effect_used_for_direction"] = effect_for_direction
        row["direction_matches_prespecified"] = bool(
            np.sign(effect_for_direction) == np.sign(candidate.expected_direction)
        )
        rows.append(row)
    output = pd.DataFrame(rows)
    if not output.empty:
        output["q_value_bh_within_cohort"] = multipletests(
            output["primary_p_value"], method="fdr_bh"
        )[1]
    return output


def source_manifest(paths_to_urls: dict[Path, str]) -> pd.DataFrame:
    rows = []
    for path, url in paths_to_urls.items():
        rows.append(
            {
                "relative_path": str(path.relative_to(ROOT)),
                "source_url": url,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    signature = pd.read_csv(SIGNATURE_PATH, sep="\t")
    signature["rank_within_direction"] = pd.to_numeric(signature["rank_within_direction"])
    candidates = pd.read_csv(CANDIDATE_PATH, sep="\t")
    wanted = set(signature["gene"]) | set(candidates["gene"]) | set(PROLIFERATION_CONTROL)
    gpl570 = geo_annotation_mapping(GPL570_PATH)
    gpl6480 = geo_annotation_mapping(GPL6480_PATH)

    base_cohorts: dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]] = {
        "GSE100179": read_gse100179(wanted),
        "GSE37364": read_gse37364(wanted, gpl570),
        "GSE164541": read_gse164541(wanted),
        "GSE20916": read_gse20916(wanted, gpl570),
        "GSE226739": read_gse226739(wanted),
        "GSE71187": read_gse71187(wanted, gpl6480),
    }

    cohort_specs = [
        {
            "cohort": "GSE100179",
            "base": "GSE100179",
            "paired": False,
            "recruitment_cluster": "SOTE_Budapest",
            "primary_inference_eligible": True,
            "subset": lambda meta: pd.Series(True, index=meta.index),
            "data_type": "HTA2.0 RMA gene-level microarray",
            "design_limitation": "Shares recruitment institute and investigators with GSE37364; do not count as a fully independent center.",
        },
        {
            "cohort": "GSE37364",
            "base": "GSE37364",
            "paired": False,
            "recruitment_cluster": "SOTE_Budapest",
            "primary_inference_eligible": True,
            "subset": lambda meta: pd.Series(True, index=meta.index),
            "data_type": "Affymetrix GPL570 microarray",
            "design_limitation": "Same recruitment institute/investigator group as GSE100179; older, distinct accession and BioProject.",
        },
        {
            "cohort": "GSE164541_paired",
            "base": "GSE164541",
            "paired": True,
            "recruitment_cluster": "GSE164541_paired_set",
            "primary_inference_eligible": True,
            "subset": lambda meta: pd.Series(True, index=meta.index),
            "data_type": "RNA-seq raw counts; log2(CPM+0.5)",
            "design_limitation": "Only five patient triplets; paired estimate is precise in design but small in sample size.",
        },
        {
            "cohort": "GSE20916_micro_crypt",
            "base": "GSE20916",
            "paired": False,
            "recruitment_cluster": "GSE20916_Warsaw",
            "primary_inference_eligible": True,
            "subset": lambda meta: meta["stratum"].eq("micro") & meta["cell_compartment"].eq("crypt_epithelium") & meta["tissue_group"].isin(["normal", "adenoma"]),
            "data_type": "GPL570 microdissected crypt epithelium",
            "design_limitation": "Technical stratum of GSE20916; correlated with the mucosal and macrodissected strata. Distant normal samples excluded.",
        },
        {
            "cohort": "GSE20916_micro_mucosa",
            "base": "GSE20916",
            "paired": False,
            "recruitment_cluster": "GSE20916_Warsaw",
            "primary_inference_eligible": True,
            "subset": lambda meta: meta["stratum"].eq("micro") & meta["cell_compartment"].eq("mucosa_or_bulk") & meta["tissue_group"].isin(["normal", "adenoma"]),
            "data_type": "GPL570 microdissected mucosa",
            "design_limitation": "Technical stratum of GSE20916; correlated with the crypt and macrodissected strata. Distant normal samples excluded.",
        },
        {
            "cohort": "GSE20916_macro",
            "base": "GSE20916",
            "paired": False,
            "recruitment_cluster": "GSE20916_Warsaw",
            "primary_inference_eligible": True,
            "subset": lambda meta: meta["stratum"].eq("macro") & meta["tissue_group"].isin(["normal", "adenoma"]),
            "data_type": "GPL570 macrodissected tissue",
            "design_limitation": "Technical stratum of GSE20916 and not independent of its microdissected evidence route.",
        },
        {
            "cohort": "GSE226739",
            "base": "GSE226739",
            "paired": False,
            "recruitment_cluster": "Longhua_Shanghai",
            "primary_inference_eligible": True,
            "subset": lambda meta: pd.Series(True, index=meta.index),
            "data_type": "RNA-seq author-deposited FPKM; log2(FPKM+1)",
            "design_limitation": "Five samples per group; sample mapping was recovered from the publication comparison sheets rather than inferred from column order.",
        },
        {
            "cohort": "GSE71187_descriptive",
            "base": "GSE71187",
            "paired": False,
            "recruitment_cluster": "GSE71187_multisite",
            "primary_inference_eligible": False,
            "subset": lambda meta: meta["tissue_group"].isin(["normal", "adenoma"]),
            "data_type": "Agilent GPL6480 microarray",
            "design_limitation": "Disease group is perfectly confounded with collection source/procedure/institution; direction-only stress test, not an inferential replication.",
        },
    ]

    score_frames = []
    coverage_rows: list[dict] = []
    selected_frames = []
    test_rows = []
    bootstrap_rows = []
    candidate_frames = []
    audit_rows = []
    rng = np.random.default_rng(SEED)

    for base, (_, _, selected) in base_cohorts.items():
        selected_frames.append(selected.assign(base_accession=base))

    for spec in cohort_specs:
        meta, expression, _ = base_cohorts[spec["base"]]
        keep = spec["subset"](meta)
        cohort_meta = meta.loc[keep].copy()
        cohort_expression = expression.reindex(cohort_meta["sample_id"])
        if not {"normal", "adenoma"}.issubset(set(cohort_meta["tissue_group"])):
            raise ValueError(f"{spec['cohort']}: normal/adenoma groups are incomplete")
        scores, coverage = score_nested_signatures(
            cohort_meta,
            cohort_expression,
            signature,
            spec["cohort"],
            normalization_groups={"normal", "adenoma"},
        )
        scores["recruitment_cluster"] = spec["recruitment_cluster"]
        scores["primary_inference_eligible"] = spec["primary_inference_eligible"]
        score_frames.append(scores)
        coverage_rows.extend(coverage)

        test = cohort_comparison(
            scores,
            spec["cohort"],
            PRIMARY_SIZE,
            "adenoma",
            "normal",
            spec["paired"],
            "adenoma_vs_normal",
        )
        test.update(
            {
                "recruitment_cluster": spec["recruitment_cluster"],
                "primary_inference_eligible": spec["primary_inference_eligible"],
                "design_limitation": spec["design_limitation"],
            }
        )
        test_rows.append(test)
        bootstrap = bootstrap_primary(scores, spec["cohort"], spec["paired"], rng)
        bootstrap["primary_inference_eligible"] = spec["primary_inference_eligible"]
        bootstrap_rows.append(bootstrap)

        candidate_test = candidate_gene_tests(
            spec["cohort"], cohort_meta, cohort_expression, candidates, spec["paired"]
        )
        candidate_test["recruitment_cluster"] = spec["recruitment_cluster"]
        candidate_test["primary_inference_eligible"] = spec["primary_inference_eligible"]
        candidate_frames.append(candidate_test)

        counts = cohort_meta["tissue_group"].value_counts()
        audit_rows.append(
            {
                "cohort": spec["cohort"],
                "base_accession": spec["base"],
                "recruitment_cluster": spec["recruitment_cluster"],
                "data_type": spec["data_type"],
                "paired_primary": spec["paired"],
                "n_normal": int(counts.get("normal", 0)),
                "n_adenoma": int(counts.get("adenoma", 0)),
                "n_complete_pairs": int(
                    cohort_meta.pivot_table(
                        index="patient_id", columns="tissue_group", values="sample_id", aggfunc="first"
                    )
                    .dropna(subset=["normal", "adenoma"])
                    .shape[0]
                )
                if spec["paired"]
                else 0,
                "primary_inference_eligible": spec["primary_inference_eligible"],
                "normalization_scope": "normal_and_adenoma_only",
                "feature_selection": "highest-variance mapped feature across all deposited base-accession samples; disease labels unused",
                "design_limitation": spec["design_limitation"],
            }
        )

    scores = pd.concat(score_frames, ignore_index=True)
    coverage = pd.DataFrame(coverage_rows)
    selected = pd.concat(selected_frames, ignore_index=True)
    route_tests = pd.DataFrame(test_rows)
    bootstraps = pd.DataFrame(bootstrap_rows)
    candidate_tests = pd.concat(candidate_frames, ignore_index=True)
    audit = pd.DataFrame(audit_rows)

    cluster_route = (
        route_tests.loc[route_tests["primary_inference_eligible"]]
        .groupby("recruitment_cluster", as_index=False)
        .agg(
            n_correlated_routes=("cohort", "nunique"),
            median_route_effect=("median_difference_a_minus_b", "median"),
            all_routes_positive=("median_difference_a_minus_b", lambda x: bool((x > 0).all())),
        )
    )
    eligible_candidate = candidate_tests.loc[candidate_tests["primary_inference_eligible"]].copy()
    cluster_candidate = (
        eligible_candidate.groupby(["gene", "panel_role", "expected_direction", "recruitment_cluster"], as_index=False)
        .agg(
            n_correlated_routes=("cohort", "nunique"),
            cluster_median_effect=("effect_used_for_direction", "median"),
        )
    )
    cluster_candidate["cluster_direction_matches"] = (
        np.sign(cluster_candidate["cluster_median_effect"])
        == np.sign(cluster_candidate["expected_direction"])
    )
    candidate_cluster_summary = (
        cluster_candidate.groupby(["gene", "panel_role", "expected_direction"], as_index=False)
        .agg(
            n_recruitment_clusters_evaluable=("recruitment_cluster", "nunique"),
            n_cluster_direction_matches=("cluster_direction_matches", "sum"),
            median_cluster_effect=("cluster_median_effect", "median"),
        )
    )
    candidate_cluster_summary["cluster_direction_match_fraction"] = (
        candidate_cluster_summary["n_cluster_direction_matches"]
        / candidate_cluster_summary["n_recruitment_clusters_evaluable"]
    )

    scores.to_csv(OUT_DIR / "sample_scores.tsv", sep="\t", index=False)
    coverage.to_csv(OUT_DIR / "locked_signature_coverage.tsv", sep="\t", index=False)
    selected.to_csv(OUT_DIR / "selected_features_without_labels.tsv", sep="\t", index=False)
    route_tests.to_csv(OUT_DIR / "route_tests.tsv", sep="\t", index=False)
    bootstraps.to_csv(OUT_DIR / "route_primary_bootstrap.tsv", sep="\t", index=False)
    candidate_tests.to_csv(OUT_DIR / "candidate_gene_tests.tsv", sep="\t", index=False)
    candidate_tests.loc[candidate_tests["panel_role"].isin(CURRENT_PANEL_ROLES)].to_csv(
        OUT_DIR / "current_panel_candidate_tests.tsv", sep="\t", index=False
    )
    audit.to_csv(OUT_DIR / "cohort_design_audit.tsv", sep="\t", index=False)
    cluster_route.to_csv(OUT_DIR / "route_recruitment_cluster_summary.tsv", sep="\t", index=False)
    candidate_cluster_summary.to_csv(
        OUT_DIR / "candidate_recruitment_cluster_summary.tsv", sep="\t", index=False
    )

    manifest = source_manifest(
        {
            SOURCE_DIR / "GSE100179" / "GSE100179_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100179/matrix/GSE100179_series_matrix.txt.gz",
            SOURCE_DIR / "GSE100179" / "GSE100179_2017-07-21-annotated-RMA-SKETCH.RMA-GENE-FULL-Group1.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100179/suppl/GSE100179_2017-07-21-annotated-RMA-SKETCH.RMA-GENE-FULL-Group1.txt.gz",
            SOURCE_DIR / "GSE37364" / "GSE37364_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE37nnn/GSE37364/matrix/GSE37364_series_matrix.txt.gz",
            SOURCE_DIR / "GSE164541" / "GSE164541_ANT_count.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164541/suppl/GSE164541_ANT_count.csv.gz",
            SOURCE_DIR / "GSE164541" / "GSE164541_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164541/matrix/GSE164541_series_matrix.txt.gz",
            SOURCE_DIR / "GSE20916" / "GSE20916_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20916/matrix/GSE20916_series_matrix.txt.gz",
            SOURCE_DIR / "GSE226739" / "GSE226739_mRNA_Expression_Profiling.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226739/suppl/GSE226739_mRNA_Expression_Profiling.xlsx",
            SOURCE_DIR / "GSE226739" / "GSE226739_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226739/matrix/GSE226739_series_matrix.txt.gz",
            SOURCE_DIR / "GSE226739" / "PMC11574242_supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11574242/supplementaryFiles",
            SOURCE_DIR / "GSE71187" / "GSE71187_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71187/matrix/GSE71187_series_matrix.txt.gz",
            GPL570_PATH: "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
            GPL6480_PATH: "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz",
        }
    )
    manifest.to_csv(OUT_DIR / "source_file_manifest.tsv", sep="\t", index=False)

    primary = route_tests.loc[route_tests["primary_inference_eligible"]]
    panel_summary = candidate_cluster_summary.loc[
        candidate_cluster_summary["panel_role"].isin(CURRENT_PANEL_ROLES)
    ].sort_values(["cluster_direction_match_fraction", "n_recruitment_clusters_evaluable"], ascending=False)
    lines = [
        "Additional public adenoma transcriptomic validation",
        "==================================================",
        f"Frozen-score routes analyzed: {len(route_tests)}; inferentially eligible: {len(primary)}.",
        f"Independent recruitment clusters after correlation control: {cluster_route['recruitment_cluster'].nunique()}.",
        f"Eligible routes with positive k50 adenoma-minus-normal effect: {int((primary['median_difference_a_minus_b'] > 0).sum())}/{len(primary)}.",
        f"Recruitment clusters with positive median k50 effect: {int((cluster_route['median_route_effect'] > 0).sum())}/{len(cluster_route)}.",
        "",
        "Design guardrails",
        "- GSE100179 and GSE37364 are separate accessions/BioProjects but share the SOTE recruitment group; they count as one center cluster.",
        "- GSE20916 micro-crypt, micro-mucosa, and macro routes are technical strata of one study and count as one recruitment cluster.",
        "- GSE71187 is direction-only because disease status is confounded with collection source and procedure.",
        "- GSE164541 is paired but has only five patient triplets; exact paired inference is necessarily low-powered.",
        "",
        "Current wet-lab panel: independent-cluster transcript direction",
    ]
    for row in panel_summary.itertuples(index=False):
        lines.append(
            f"- {row.gene} ({row.panel_role}): {int(row.n_cluster_direction_matches)}/{int(row.n_recruitment_clusters_evaluable)} clusters match the prespecified direction; median effect={row.median_cluster_effect:.3f}."
        )
    (OUT_DIR / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
