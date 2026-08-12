#!/usr/bin/env python3
"""Candidate-focused reanalysis of public colorectal adenoma proteomics PXD000445.

Two evidence layers are kept separate:
1. PRIDE mzTab identification files establish experiment-level detectability.
2. The publisher-deposited PSM table supports a conservative paired iTRAQ
   reanalysis using unique, confident peptides. This is not represented as an
   exact reproduction of the authors' protein-family assembly pipeline.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "data_sources" / "public_adenoma_proteomics" / "PXD000445"
OUT_DIR = ROOT / "results" / "pxd000445_candidate_reanalysis"
CANDIDATE_DIR = ROOT / "results" / "public_adenoma_protein_triangulation"
CANDIDATE_PATH = CANDIDATE_DIR / "candidate_inventory.tsv"
EXISTING_EVIDENCE_PATH = CANDIDATE_DIR / "candidate_public_protein_evidence_matrix.tsv"
TRANSCRIPT_SUMMARY_PATH = (
    ROOT
    / "results"
    / "expanded_public_adenoma_validation"
    / "candidate_recruitment_cluster_summary.tsv"
)
HPA_PATH = ROOT / "metadata" / "hpa_candidate_assayability_2026-07-14.tsv"
SUPPLEMENT_ARCHIVE = SOURCE_DIR / "S1535947620330978_mmc1.zip"
PSM_MEMBER = "mcp.M113.035105-8.xlsx"

AUTHOR_QC_EXCLUDED_EXPERIMENTS = {3, 4, 5}
PAIR_CHANNELS = [("115", "116"), ("117", "118"), ("119", "121")]
EXPECTED_POSITIVE_CONTROLS = {"SORD": 1}
EXTRA_ACCESSIONS = {"SORD": {"Q00796"}}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_mztab_proteins(
    paths: list[Path], wanted_genes: set[str]
) -> tuple[pd.DataFrame, dict[str, str], pd.DataFrame]:
    experiment_rows = []
    accession_gene_votes: dict[str, Counter] = defaultdict(Counter)
    raw_rows = []
    for path in paths:
        experiment_match = re.search(r"(\d+)\.pride\.mztab", path.name)
        experiment_accession = experiment_match.group(1) if experiment_match else path.stem
        header: list[str] | None = None
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for parts in reader:
                if not parts:
                    continue
                if parts[0] == "PRH":
                    header = parts[1:]
                    continue
                if parts[0] != "PRT" or header is None:
                    continue
                row = dict(zip(header, parts[1:], strict=False))
                description = row.get("description", "")
                match = re.search(r"\bGN=([^\s]+)", description)
                if not match:
                    continue
                gene = match.group(1).strip()
                accession = row.get("accession", "").strip()
                accession_gene_votes[accession][gene] += 1
                if gene not in wanted_genes:
                    continue
                psms = pd.to_numeric(row.get("num_psms_ms_run[1]"), errors="coerce")
                peptides = pd.to_numeric(
                    row.get("num_peptides_distinct_ms_run[1]"), errors="coerce"
                )
                raw_rows.append(
                    {
                        "gene": gene,
                        "pride_result_accession": experiment_accession,
                        "mztab_file": path.name,
                        "protein_accession": accession,
                        "num_psms": psms,
                        "num_distinct_peptides": peptides,
                    }
                )
    raw = pd.DataFrame(raw_rows)
    if not raw.empty:
        experiment_rows = (
            raw.groupby(["gene", "pride_result_accession", "mztab_file"], as_index=False)
            .agg(
                protein_accessions=("protein_accession", lambda x: ";".join(sorted(set(x)))),
                max_num_psms=("num_psms", "max"),
                max_num_distinct_peptides=("num_distinct_peptides", "max"),
            )
        )
    experiment_table = pd.DataFrame(experiment_rows)
    accession_to_gene = {
        accession: votes.most_common(1)[0][0]
        for accession, votes in accession_gene_votes.items()
        if len(votes) > 0
    }
    return experiment_table, accession_to_gene, raw


def add_candidate_accession_mapping(
    accession_to_gene: dict[str, str], candidates: pd.DataFrame
) -> dict[str, str]:
    mapping = dict(accession_to_gene)
    ffpe = pd.read_csv(CANDIDATE_DIR / "pxd017269_ffpe_detectability.tsv", sep="\t")
    for row in ffpe[["gene", "pxd017269_protein_accessions"]].itertuples(index=False):
        if pd.isna(row.pxd017269_protein_accessions):
            continue
        for accession in str(row.pxd017269_protein_accessions).split(";"):
            accession = accession.strip()
            if accession:
                mapping.setdefault(accession, row.gene)
                mapping.setdefault(accession.split("-", 1)[0], row.gene)
    for gene, accessions in EXTRA_ACCESSIONS.items():
        for accession in accessions:
            mapping.setdefault(accession, gene)
    return mapping


def conservative_psm_reanalysis(
    accession_to_gene: dict[str, str], wanted_genes: set[str]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    with zipfile.ZipFile(SUPPLEMENT_ARCHIVE) as zf:
        workbook_bytes = io.BytesIO(zf.read(PSM_MEMBER))
    workbook = openpyxl.load_workbook(workbook_bytes, read_only=True, data_only=True)
    worksheet = workbook["Complete PSM Table"]
    rows = worksheet.iter_rows(values_only=True)
    header = [str(value) for value in next(rows)]
    index = {name: position for position, name in enumerate(header)}
    required = {
        "ExperimentInformation",
        "proteinAccession",
        "PeptideSequence",
        "ExpectationValue",
        "isBold_0_1",
        "isUniq_0_1",
    }
    if not required.issubset(index):
        raise ValueError(f"PSM table missing fields: {sorted(required - set(index))}")
    label_ratio_columns = []
    for number in range(1, 8):
        label_ratio_columns.append(
            (index[f"iTRAQ_label_{number}"], index[f"iTRAQ_ratio_{number}"])
        )

    peptide_ratios: dict[tuple[int, str, str, str], list[float]] = defaultdict(list)
    filter_counts = Counter()
    for row in rows:
        filter_counts["psm_rows_total"] += 1
        experiment_text = str(row[index["ExperimentInformation"]])
        experiment_match = re.search(r"exp0?(\d+)$", experiment_text, flags=re.IGNORECASE)
        if not experiment_match:
            filter_counts["missing_experiment_number"] += 1
            continue
        experiment = int(experiment_match.group(1))
        accession = str(row[index["proteinAccession"]]).strip()
        gene = accession_to_gene.get(accession) or accession_to_gene.get(
            accession.split("-", 1)[0]
        )
        if gene not in wanted_genes:
            continue
        filter_counts["candidate_psm_rows"] += 1
        expectation = pd.to_numeric(row[index["ExpectationValue"]], errors="coerce")
        is_bold = pd.to_numeric(row[index["isBold_0_1"]], errors="coerce")
        is_unique = pd.to_numeric(row[index["isUniq_0_1"]], errors="coerce")
        if not np.isfinite(expectation) or expectation > 0.05:
            filter_counts["candidate_rows_failed_expectation"] += 1
            continue
        if is_bold != 1:
            filter_counts["candidate_rows_not_master_assignment"] += 1
            continue
        if is_unique != 1:
            filter_counts["candidate_rows_not_unique_peptide"] += 1
            continue
        peptide = str(row[index["PeptideSequence"]]).strip()
        filter_counts["candidate_rows_retained"] += 1
        for label_index, ratio_index in label_ratio_columns:
            label = str(row[label_index])
            if "/" not in label:
                continue
            channel = label.split("/", 1)[0]
            value = pd.to_numeric(row[ratio_index], errors="coerce")
            if np.isfinite(value) and value > 0:
                peptide_ratios[(experiment, gene, peptide, channel)].append(float(np.log2(value)))

    peptide_rows = []
    for (experiment, gene, peptide, channel), values in peptide_ratios.items():
        peptide_rows.append(
            {
                "experiment": experiment,
                "gene": gene,
                "peptide": peptide,
                "channel": channel,
                "log2_ratio_vs_114": float(np.median(values)),
                "n_psms_for_peptide_channel": len(values),
            }
        )
    peptide_table = pd.DataFrame(peptide_rows)
    if peptide_table.empty:
        raise RuntimeError("No candidate peptide ratios passed the prespecified filters")
    protein_channel = (
        peptide_table.groupby(["experiment", "gene", "channel"], as_index=False)
        .agg(
            protein_log2_ratio_vs_114=("log2_ratio_vs_114", "median"),
            n_unique_peptides=("peptide", "nunique"),
            n_psms=("n_psms_for_peptide_channel", "sum"),
        )
    )
    wide = protein_channel.pivot_table(
        index=["experiment", "gene"],
        columns="channel",
        values="protein_log2_ratio_vs_114",
        aggfunc="median",
    )
    pair_rows = []
    for (experiment, gene), values in wide.iterrows():
        for pair_number, (normal_channel, adenoma_channel) in enumerate(PAIR_CHANNELS, start=1):
            if normal_channel not in values or adenoma_channel not in values:
                continue
            normal = values.get(normal_channel)
            adenoma = values.get(adenoma_channel)
            if not np.isfinite(normal) or not np.isfinite(adenoma):
                continue
            normal_peptides = protein_channel.loc[
                protein_channel["experiment"].eq(experiment)
                & protein_channel["gene"].eq(gene)
                & protein_channel["channel"].eq(normal_channel),
                "n_unique_peptides",
            ].iloc[0]
            adenoma_peptides = protein_channel.loc[
                protein_channel["experiment"].eq(experiment)
                & protein_channel["gene"].eq(gene)
                & protein_channel["channel"].eq(adenoma_channel),
                "n_unique_peptides",
            ].iloc[0]
            pair_rows.append(
                {
                    "experiment": experiment,
                    "gene": gene,
                    "pair_number_within_experiment": pair_number,
                    "pair_id": f"Exp{experiment:02d}_Pair{pair_number}",
                    "normal_channel": normal_channel,
                    "adenoma_channel": adenoma_channel,
                    "normal_log2_ratio_vs_114": float(normal),
                    "adenoma_log2_ratio_vs_114": float(adenoma),
                    "adenoma_minus_normal_log2": float(adenoma - normal),
                    "n_unique_peptides_normal": int(normal_peptides),
                    "n_unique_peptides_adenoma": int(adenoma_peptides),
                    "passes_author_batch_qc": experiment
                    not in AUTHOR_QC_EXCLUDED_EXPERIMENTS,
                }
            )
    pair_table = pd.DataFrame(pair_rows)
    filter_table = pd.DataFrame(
        [{"metric": key, "value": value} for key, value in sorted(filter_counts.items())]
    )
    return peptide_table, pair_table, filter_table


def paired_gene_summary(
    pair_table: pd.DataFrame, candidates: pd.DataFrame
) -> pd.DataFrame:
    expected = candidates.set_index("gene")["expected_direction"].to_dict()
    expected.update(EXPECTED_POSITIVE_CONTROLS)
    rows = []
    for analysis_set, subset in [
        ("all_30_pairs_deposited", pair_table),
        ("author_qc_21_pairs", pair_table.loc[pair_table["passes_author_batch_qc"]]),
    ]:
        for gene, part in subset.groupby("gene", sort=True):
            delta = part["adenoma_minus_normal_log2"].dropna()
            if delta.empty:
                continue
            if np.allclose(delta, 0):
                wilcoxon_p = 1.0
            elif len(delta) >= 2:
                wilcoxon_p = float(
                    stats.wilcoxon(delta, alternative="two-sided", zero_method="wilcox").pvalue
                )
            else:
                wilcoxon_p = np.nan
            expected_direction = expected.get(gene, np.nan)
            effect = float(delta.median())
            rows.append(
                {
                    "analysis_set": analysis_set,
                    "gene": gene,
                    "expected_direction": expected_direction,
                    "n_pairs": len(delta),
                    "n_experiments": part["experiment"].nunique(),
                    "median_adenoma_minus_normal_log2": effect,
                    "mean_adenoma_minus_normal_log2": float(delta.mean()),
                    "paired_positive_fraction": float((delta > 0).mean()),
                    "p_paired_wilcoxon": wilcoxon_p,
                    "min_unique_peptides_per_sample": int(
                        part[["n_unique_peptides_normal", "n_unique_peptides_adenoma"]]
                        .min(axis=1)
                        .min()
                    ),
                    "direction_matches_prespecified": (
                        bool(np.sign(effect) == np.sign(expected_direction))
                        if np.isfinite(expected_direction)
                        else np.nan
                    ),
                }
            )
    output = pd.DataFrame(rows)
    output["q_value_bh_within_analysis_set"] = np.nan
    for analysis_set, index_values in output.groupby("analysis_set").groups.items():
        index_values = list(index_values)
        finite = output.loc[index_values, "p_paired_wilcoxon"].notna()
        finite_indices = list(np.array(index_values)[finite.to_numpy()])
        if finite_indices:
            output.loc[finite_indices, "q_value_bh_within_analysis_set"] = multipletests(
                output.loc[finite_indices, "p_paired_wilcoxon"], method="fdr_bh"
            )[1]
    return output


def mztab_detection_summary(
    experiment_table: pd.DataFrame, genes: list[str], n_files: int
) -> pd.DataFrame:
    rows = []
    for gene in genes:
        part = experiment_table.loc[experiment_table["gene"].eq(gene)]
        rows.append(
            {
                "gene": gene,
                "pxd000445_mztab_detected": not part.empty,
                "pxd000445_n_mztab_experiments_detected": part[
                    "pride_result_accession"
                ].nunique(),
                "pxd000445_n_mztab_experiments_total": n_files,
                "pxd000445_mztab_experiment_detection_fraction": (
                    part["pride_result_accession"].nunique() / n_files
                ),
                "pxd000445_max_distinct_peptides": (
                    float(part["max_num_distinct_peptides"].max()) if not part.empty else np.nan
                ),
                "pxd000445_max_psms": (
                    float(part["max_num_psms"].max()) if not part.empty else np.nan
                ),
                "pxd000445_accessions": (
                    ";".join(
                        sorted(
                            {
                                accession
                                for value in part["protein_accessions"].dropna()
                                for accession in str(value).split(";")
                            }
                        )
                    )
                    if not part.empty
                    else ""
                ),
            }
        )
    return pd.DataFrame(rows)


def integrated_wetlab_table(
    detection: pd.DataFrame, paired_summary: pd.DataFrame
) -> pd.DataFrame:
    existing = pd.read_csv(EXISTING_EVIDENCE_PATH, sep="\t")
    transcript = pd.read_csv(TRANSCRIPT_SUMMARY_PATH, sep="\t")
    hpa = pd.read_csv(HPA_PATH, sep="\t")
    paired_qc = paired_summary.loc[
        paired_summary["analysis_set"].eq("author_qc_21_pairs")
    ].add_prefix("pxd000445_psm_")
    paired_qc = paired_qc.rename(columns={"pxd000445_psm_gene": "gene"})
    output = (
        existing.merge(detection, on="gene", how="left")
        .merge(paired_qc, on="gene", how="left")
        .merge(transcript, on=["gene", "panel_role", "expected_direction"], how="left")
        .merge(hpa, on="gene", how="left", suffixes=("", "_hpa"))
    )
    output["pxd000445_quantitatively_evaluable"] = output[
        "pxd000445_psm_n_pairs"
    ].fillna(0).gt(0)
    output["new_transcript_has_direction_conflict"] = output[
        "cluster_direction_match_fraction"
    ].fillna(0).lt(1)
    output["pxd000445_psm_has_direction_conflict"] = (
        output["pxd000445_quantitatively_evaluable"]
        & output["pxd000445_psm_direction_matches_prespecified"].eq(False)
    )
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
    candidates = pd.read_csv(CANDIDATE_PATH, sep="\t")
    wanted_genes = set(candidates["gene"]) | set(EXPECTED_POSITIVE_CONTROLS)
    mztab_paths = sorted((SOURCE_DIR / "generated").glob("*.mztab.gz"))
    if len(mztab_paths) != 11:
        raise ValueError(f"Expected 11 mzTab files, found {len(mztab_paths)}")

    experiment_table, accession_to_gene, raw_mztab = parse_mztab_proteins(
        mztab_paths, wanted_genes
    )
    accession_to_gene = add_candidate_accession_mapping(accession_to_gene, candidates)
    detection = mztab_detection_summary(
        experiment_table, sorted(wanted_genes), len(mztab_paths)
    )
    peptide_table, pair_table, filter_table = conservative_psm_reanalysis(
        accession_to_gene, wanted_genes
    )
    paired_summary = paired_gene_summary(pair_table, candidates)
    integrated = integrated_wetlab_table(detection, paired_summary)

    experiment_table.to_csv(
        OUT_DIR / "mztab_candidate_experiment_detection.tsv", sep="\t", index=False
    )
    detection.to_csv(OUT_DIR / "mztab_candidate_detection_summary.tsv", sep="\t", index=False)
    peptide_table.to_csv(
        OUT_DIR / "psm_unique_peptide_channel_ratios.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pair_table.to_csv(OUT_DIR / "psm_candidate_paired_deltas.tsv", sep="\t", index=False)
    paired_summary.to_csv(OUT_DIR / "psm_candidate_paired_tests.tsv", sep="\t", index=False)
    filter_table.to_csv(OUT_DIR / "psm_filter_audit.tsv", sep="\t", index=False)
    integrated.to_csv(OUT_DIR / "wetlab_candidate_integrated_evidence.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "experiment": number,
                "passes_author_batch_qc": number not in AUTHOR_QC_EXCLUDED_EXPERIMENTS,
                "qc_basis": (
                    "Retained in the authors' 21-pair analysis"
                    if number not in AUTHOR_QC_EXCLUDED_EXPERIMENTS
                    else "Excluded: samples formed the 18-sample experiment-specific cluster in Supplementary Figure 1"
                ),
            }
            for number in range(1, 11)
        ]
    ).to_csv(OUT_DIR / "author_batch_qc_reconstruction.tsv", sep="\t", index=False)

    manifest_paths = {
        SOURCE_DIR / "README.txt": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2014/04/PXD000445/README.txt",
        SUPPLEMENT_ARCHIVE: "https://ars.els-cdn.com/content/image/1-s2.0-S1535947620330978-mmc1.zip",
    }
    for path in mztab_paths:
        manifest_paths[path] = (
            "https://ftp.pride.ebi.ac.uk/pride/data/archive/2014/04/PXD000445/generated/"
            + path.name
        )
    source_manifest(manifest_paths).to_csv(
        OUT_DIR / "source_file_manifest.tsv", sep="\t", index=False
    )

    focus_genes = [
        "OLFM4",
        "FABP1",
        "ETHE1",
        "CA2",
        "EPHB3",
        "CTNNB1",
        "SOX9",
        "MKI67",
        "KRT20",
    ]
    focus = integrated.set_index("gene").reindex(focus_genes)
    lines = [
        "PXD000445 candidate-focused proteomics reanalysis",
        "=================================================",
        f"PRIDE identification files audited: {len(mztab_paths)}.",
        "The public mzTab files are identification-only; no differential effect is inferred from mzTab detection.",
        "Publisher PSM table: 30 paired samples in 10 experiments; the authors' QC retained 21 pairs in experiments 1, 2, and 6-10.",
        "Paired iTRAQ reanalysis uses expectation <=0.05, master protein assignment, unique peptides, peptide-level median aggregation, and within-pair adenoma-minus-normal contrasts.",
        "This candidate-focused aggregation is a sensitivity reanalysis, not an exact reconstruction of the authors' protein-family assembly model.",
        "",
        "Current and reserve panel evidence",
    ]
    for gene, row in focus.iterrows():
        if pd.isna(row.get("panel_role")):
            continue
        psm_effect = row.get("pxd000445_psm_median_adenoma_minus_normal_log2")
        psm_text = (
            f"paired effect={psm_effect:.3f}, n={int(row['pxd000445_psm_n_pairs'])}"
            if pd.notna(psm_effect)
            else "paired effect not estimable"
        )
        lines.append(
            f"- {gene}: existing tier={row['public_protein_evidence_tier']}; "
            f"mzTab detection={row['pxd000445_mztab_experiment_detection_fraction']:.3f}; "
            f"{psm_text}; transcript cluster concordance={row['cluster_direction_match_fraction']:.3f}."
        )
    sord = paired_summary.loc[
        paired_summary["analysis_set"].eq("author_qc_21_pairs")
        & paired_summary["gene"].eq("SORD")
    ]
    if not sord.empty:
        row = sord.iloc[0]
        lines.extend(
            [
                "",
                f"Positive-control SORD: paired median log2 effect={row['median_adenoma_minus_normal_log2']:.3f}, n={int(row['n_pairs'])}, Wilcoxon p={row['p_paired_wilcoxon']:.3g}.",
            ]
        )
    (OUT_DIR / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
