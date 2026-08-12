#!/usr/bin/env python3
"""Lightweight Becker multiome scATAC TSS accessibility analysis.

This is a feasibility analysis, not a full peak-calling workflow. It uses
tabix-indexed 10x fragments to count accessibility in +/-2 kb GENCODE TSS
windows for the manuscript's marker modules, then tests whether polyp samples
show higher WNT promoter accessibility than normal/unaffected mucosa.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
import shutil
import tarfile
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data_sources" / "Becker_NatGenet_2022_GEO"
RAW_TAR = DATA_DIR / "GSE201349_RAW_multiome_large.tar"
SERIES_MATRIX = DATA_DIR / "GSE201348_series_matrix.txt.gz"
GENCODE_GTF = ROOT / "references" / "genome" / "gencode.v44.annotation.gtf.gz"
EXTRACT_DIR = DATA_DIR / "multiome_tss_subset"
OUT_DIR = ROOT / "results" / "becker_multiome_tss"

WINDOW_BP = 2_000
MAX_BY_GROUP = {"normal_unaffected": 16, "polyp": 24, "crc": 4}

GENE_SETS = {
    "wnt_stem": ["LGR5", "ASCL2", "OLFM4", "AXIN2", "SOX9", "EPHB2", "SMOC2"],
    "wnt_core_ihc": ["OLFM4", "SOX9", "EPHB2"],
    "serrated_metaplasia": ["MUC5AC", "MUC6", "TFF1", "TFF2", "TFF3", "REG4", "AGR2", "SPINK4", "KRT7", "ANXA10"],
    "antigen_presentation_ifn": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "B2M", "TAP1", "STAT1", "IRF1", "CXCL10", "CXCL11"],
    "proliferation_control": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "TYMS", "UBE2C", "CENPF"],
    "housekeeping_control": ["ACTB", "GAPDH", "RPLP0", "RPL13A", "TUBB", "HPRT1", "PPIA", "RPS18", "TBP", "GUSB", "YWHAZ", "SDHA", "HMBS"],
}

PRIMARY_MODULES = ["wnt_stem", "wnt_core_ihc", "serrated_metaplasia", "antigen_presentation_ifn", "proliferation_control", "housekeeping_control"]
CONTRASTS = {
    "wnt_stem_minus_housekeeping": ("score__wnt_stem", "score__housekeeping_control"),
    "wnt_core_ihc_minus_housekeeping": ("score__wnt_core_ihc", "score__housekeeping_control"),
}


def clean(value: str) -> str:
    return value.strip().strip('"')


def normalize_column(value: str) -> str:
    value = value.lower().replace(" ", "_").replace("-", "_")
    return re.sub(r"[^a-z0-9_]+", "", value)


def patient_id_from_sample_name(value: str) -> str:
    value = str(value).strip()
    if not value:
        return "unknown"
    if value.startswith("CRC"):
        return value.split("_", 1)[0]
    if "-" in value:
        return value.split("-", 1)[0]
    return re.sub(r"(_.*)$", "", value)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    rows: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for parts in reader:
            if not parts:
                continue
            if parts[0] == "!series_matrix_table_begin":
                break
            if not parts[0].startswith("!Sample_"):
                continue
            key = parts[0].lstrip("!")
            rows.setdefault(key, []).append([clean(v) for v in parts[1:]])

    meta = pd.DataFrame(
        {
            "scrna_geo_accession": rows["Sample_geo_accession"][0],
            "sample_title": rows["Sample_title"][0],
        }
    )
    for values in rows.get("Sample_characteristics_ch1", []):
        parsed = [v.split(":", 1) if ":" in v else [v, ""] for v in values]
        keys = [normalize_column(p[0]) for p in parsed]
        if not keys:
            continue
        col = pd.Series(keys).mode().iloc[0]
        meta[col] = [p[1].strip() if len(p) > 1 else "" for p in parsed]

    meta["sample_core"] = (
        meta["sample_title"]
        .str.replace(", snRNAseq", "", regex=False)
        .str.replace(", Replicate1", "", regex=False)
        .str.replace(", Replicate2", "", regex=False)
        .str.strip()
    )
    meta["patient_id"] = meta["sample_core"].map(patient_id_from_sample_name)
    meta["disease_stage_group"] = meta["disease_stage"].map(
        {
            "Normal": "normal_unaffected",
            "Unaffected": "normal_unaffected",
            "Polyp": "polyp",
            "CRC": "crc",
        }
    )
    meta["familial_adenomatous_polyposis"] = meta["familial_adenomatous_polyposis"].replace({"": "unknown"}).fillna("unknown")
    meta["sex"] = meta["sex"].replace({"": "unknown"}).fillna("unknown")
    keep_cols = [
        "sample_core",
        "scrna_geo_accession",
        "sample_title",
        "disease_stage",
        "disease_stage_group",
        "familial_adenomatous_polyposis",
        "sex",
        "patient_id",
    ]
    return meta[keep_cols].drop_duplicates("sample_core")


def fragment_sample_core(member_name: str) -> str:
    stem = Path(member_name).name.removesuffix("_fragments.tsv.gz")
    raw = stem.split("_", 1)[1] if "_" in stem else stem
    return re.sub(r"-D_\d+$", "", raw)


def list_fragment_members(tar_path: Path) -> pd.DataFrame:
    rows = []
    index_by_prefix: dict[str, str] = {}
    with tarfile.open(tar_path, "r") as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            if name.endswith("_fragments.tsv.gz.tbi.gz"):
                prefix = name.removesuffix(".tbi.gz")
                index_by_prefix[prefix] = member.name
            elif name.endswith("_fragments.tsv.gz"):
                rows.append(
                    {
                        "multiome_geo_accession": name.split("_", 1)[0],
                        "sample_core": fragment_sample_core(name),
                        "fragment_member": member.name,
                        "fragment_file": name,
                        "fragment_size": int(member.size),
                    }
                )
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise RuntimeError(f"No fragment files found in {tar_path}")
    frame["index_member"] = frame["fragment_file"].map(index_by_prefix)
    frame["index_file_gz"] = frame["fragment_file"] + ".tbi.gz"
    frame["fragment_gb"] = frame["fragment_size"] / 1024**3
    return frame


def select_samples(fragment_meta: pd.DataFrame, sample_meta: pd.DataFrame) -> pd.DataFrame:
    merged = fragment_meta.merge(sample_meta, on="sample_core", how="left")
    merged = merged.dropna(subset=["disease_stage_group"]).copy()
    merged = merged.sort_values(["sample_core", "fragment_size", "multiome_geo_accession"])
    merged = merged.drop_duplicates("sample_core", keep="first")

    selected = []
    for group, max_n in MAX_BY_GROUP.items():
        part = merged.loc[merged["disease_stage_group"] == group].sort_values(["fragment_size", "sample_core"]).head(max_n)
        selected.append(part)
    out = pd.concat(selected, ignore_index=True)
    group_rank = {"normal_unaffected": 0, "polyp": 1, "crc": 2}
    out["_group_rank"] = out["disease_stage_group"].map(group_rank).fillna(99)
    out = out.sort_values(["_group_rank", "fragment_size", "sample_core"]).drop(columns="_group_rank").reset_index(drop=True)
    return out


def extract_selected(selected: pd.DataFrame) -> None:
    EXTRACT_DIR.mkdir(parents=True, exist_ok=True)
    needed = set(selected["fragment_member"]) | set(selected["index_member"].dropna())
    name_to_size: dict[str, int] = {}
    for _, row in selected.iterrows():
        name_to_size[row["fragment_member"]] = int(row["fragment_size"])

    with tarfile.open(RAW_TAR, "r") as tar:
        extracted = 0
        for member in tar:
            if member.name not in needed:
                continue
            target = EXTRACT_DIR / Path(member.name).name
            expected_size = int(member.size)
            if target.exists() and target.stat().st_size == expected_size:
                continue
            tmp = target.with_suffix(target.suffix + ".partial")
            print(f"Extracting {target.name} ({expected_size / 1024**3:.2f} GB)", flush=True)
            source = tar.extractfile(member)
            if source is None:
                raise FileNotFoundError(member.name)
            with source, open(tmp, "wb") as out:
                shutil.copyfileobj(source, out, length=16 * 1024 * 1024)
            tmp.replace(target)
            extracted += 1
        print(f"Extraction step complete; new files extracted: {extracted}", flush=True)

    for index_gz in sorted(EXTRACT_DIR.glob("*_fragments.tsv.gz.tbi.gz")):
        index = EXTRACT_DIR / index_gz.name.removesuffix(".gz")
        if index.exists() and index.stat().st_size > 0:
            continue
        print(f"Decompressing index {index_gz.name}", flush=True)
        tmp = index.with_suffix(index.suffix + ".partial")
        with gzip.open(index_gz, "rb") as src, open(tmp, "wb") as dst:
            shutil.copyfileobj(src, dst)
        tmp.replace(index)


def parse_gtf_attributes(value: str) -> dict[str, str]:
    attrs = {}
    for item in value.rstrip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, raw = item.split(" ", 1)
        attrs[key] = raw.strip().strip('"')
    return attrs


def load_tss_windows(gtf_path: Path, genes: Iterable[str]) -> pd.DataFrame:
    wanted = set(genes)
    chosen: dict[str, dict[str, object]] = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene = attrs.get("gene_name", "")
            if gene not in wanted:
                continue
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", ""))
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            tss_1based = start if strand == "+" else end
            row = {
                "gene": gene,
                "chrom": parts[0],
                "strand": strand,
                "gene_type": gene_type,
                "tss_1based": tss_1based,
                "window_start_0based": max(0, tss_1based - 1 - WINDOW_BP),
                "window_end_0based_exclusive": tss_1based + WINDOW_BP,
            }
            prior = chosen.get(gene)
            if prior is None or (prior["gene_type"] != "protein_coding" and gene_type == "protein_coding"):
                chosen[gene] = row

    rows = [chosen[gene] for gene in sorted(chosen)]
    frame = pd.DataFrame(rows)
    missing = sorted(wanted.difference(chosen))
    if missing:
        print("Missing genes in GENCODE:", ",".join(missing), flush=True)
    return frame


def fragment_weight(parts: list[str]) -> int:
    if len(parts) >= 5:
        try:
            return int(parts[4])
        except ValueError:
            return 1
    return 1


def count_window(tabix: pysam.TabixFile, chrom: str, start: int, end: int) -> int:
    if chrom not in tabix.contigs:
        return 0
    total = 0
    for line in tabix.fetch(chrom, start, end):
        total += fragment_weight(line.rstrip("\n").split("\t"))
    return total


def count_sample_windows(sample: pd.Series, windows: pd.DataFrame) -> list[dict[str, object]]:
    frag_path = EXTRACT_DIR / sample["fragment_file"]
    index_path = Path(str(frag_path) + ".tbi")
    if not frag_path.exists():
        raise FileNotFoundError(frag_path)
    if not index_path.exists():
        raise FileNotFoundError(index_path)

    rows = []
    with pysam.TabixFile(str(frag_path)) as tabix:
        for _, window in windows.iterrows():
            count = count_window(
                tabix,
                str(window["chrom"]),
                int(window["window_start_0based"]),
                int(window["window_end_0based_exclusive"]),
            )
            rows.append(
                {
                    "multiome_geo_accession": sample["multiome_geo_accession"],
                    "sample_core": sample["sample_core"],
                    "gene": window["gene"],
                    "fragment_count": count,
                }
            )
    return rows


def module_availability_rows(windows: pd.DataFrame) -> list[dict[str, object]]:
    present = set(windows["gene"])
    rows = []
    for module, genes in GENE_SETS.items():
        hit = [gene for gene in genes if gene in present]
        miss = [gene for gene in genes if gene not in present]
        rows.append(
            {
                "module": module,
                "n_requested": len(genes),
                "n_present": len(hit),
                "present_genes": ",".join(hit),
                "missing_genes": ",".join(miss),
            }
        )
    return rows


def score_modules(selected: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    merged = counts.merge(selected, on=["multiome_geo_accession", "sample_core"], how="left")
    totals = (
        merged.groupby(["multiome_geo_accession", "sample_core"], as_index=False)["fragment_count"]
        .sum()
        .rename(columns={"fragment_count": "tss_total_counts"})
    )
    merged = merged.merge(totals, on=["multiome_geo_accession", "sample_core"], how="left")
    merged["gene_log_tss_cpm"] = np.log1p(merged["fragment_count"] / merged["tss_total_counts"].clip(lower=1) * 1_000_000.0)

    rows = []
    for (accession, sample_core), part in merged.groupby(["multiome_geo_accession", "sample_core"], sort=False):
        row = {
            "multiome_geo_accession": accession,
            "sample_core": sample_core,
            "tss_total_counts": float(part["tss_total_counts"].iloc[0]),
        }
        for module, genes in GENE_SETS.items():
            sub = part.loc[part["gene"].isin(genes)]
            row[f"module_count__{module}"] = int(sub["fragment_count"].sum()) if not sub.empty else 0
            row[f"raw__{module}"] = float(sub["gene_log_tss_cpm"].mean()) if not sub.empty else np.nan
            row[f"n_present__{module}"] = int(sub["gene"].nunique())
        rows.append(row)

    scored = pd.DataFrame(rows).merge(selected, on=["multiome_geo_accession", "sample_core"], how="left")
    for module in GENE_SETS:
        col = f"raw__{module}"
        values = scored[col].astype(float)
        sd = values.std(skipna=True)
        scored[f"score__{module}"] = (values - values.mean(skipna=True)) / sd if sd and not pd.isna(sd) else np.nan
    for contrast, (positive, negative) in CONTRASTS.items():
        scored[f"contrast__{contrast}"] = scored[positive] - scored[negative]
    scored["log10_tss_total_counts"] = np.log10(scored["tss_total_counts"].clip(lower=1))
    return scored


def safe_mannwhitney(a: pd.Series, b: pd.Series) -> float:
    a = a.dropna().astype(float)
    b = b.dropna().astype(float)
    if len(a) < 3 or len(b) < 3:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def module_tests(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    comparisons = [
        ("all_polyp_vs_normal_unaffected", scored, "polyp", "normal_unaffected"),
        ("fap_only_polyp_vs_unaffected", scored.loc[scored["familial_adenomatous_polyposis"] == "Y"], "polyp", "normal_unaffected"),
        ("crc_vs_normal_unaffected", scored, "crc", "normal_unaffected"),
    ]
    test_items = [(module, f"score__{module}") for module in PRIMARY_MODULES]
    test_items.extend((contrast, f"contrast__{contrast}") for contrast in CONTRASTS)
    for comparison, frame, group_a, group_b in comparisons:
        for module, col in test_items:
            a = frame.loc[frame["disease_stage_group"] == group_a, col].dropna().astype(float)
            b = frame.loc[frame["disease_stage_group"] == group_b, col].dropna().astype(float)
            rows.append(
                {
                    "comparison": comparison,
                    "module": module,
                    "score_column": col,
                    "group_a": group_a,
                    "group_b": group_b,
                    "n_a": len(a),
                    "n_b": len(b),
                    "median_a": float(a.median()) if len(a) else np.nan,
                    "median_b": float(b.median()) if len(b) else np.nan,
                    "delta_median_a_minus_b": float(a.median() - b.median()) if len(a) and len(b) else np.nan,
                    "p_value": safe_mannwhitney(a, b),
                    "test": "Mann-Whitney U",
                }
            )

    out = pd.DataFrame(rows)
    out["p_adj_bh_within_comparison"] = np.nan
    for comparison, idx in out.groupby("comparison").groups.items():
        p = out.loc[idx, "p_value"]
        ok = p.notna()
        if ok.any():
            out.loc[p.index[ok], "p_adj_bh_within_comparison"] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def paired_patient_tests(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    frame = scored.loc[scored["disease_stage_group"].isin(["polyp", "normal_unaffected"])].copy()
    by_patient = frame.groupby(["patient_id", "disease_stage_group"], as_index=False)[[f"score__{m}" for m in PRIMARY_MODULES]].mean()
    for module in PRIMARY_MODULES:
        col = f"score__{module}"
        wide = by_patient.pivot(index="patient_id", columns="disease_stage_group", values=col).dropna(subset=["polyp", "normal_unaffected"])
        if len(wide) >= 3:
            p_value = float(wilcoxon(wide["polyp"], wide["normal_unaffected"]).pvalue)
        else:
            p_value = np.nan
        rows.append(
            {
                "comparison": "patient_paired_polyp_vs_normal_unaffected",
                "module": module,
                "n_pairs": int(len(wide)),
                "median_delta_polyp_minus_normal": float((wide["polyp"] - wide["normal_unaffected"]).median()) if len(wide) else np.nan,
                "p_value": p_value,
                "test": "Wilcoxon signed-rank on patient-level means",
            }
        )
    out = pd.DataFrame(rows)
    ok = out["p_value"].notna()
    out["p_adj_bh"] = np.nan
    if ok.any():
        out.loc[ok, "p_adj_bh"] = multipletests(out.loc[ok, "p_value"], method="fdr_bh")[1]
    return out


def adjusted_models(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    frame = scored.loc[scored["disease_stage_group"].isin(["polyp", "normal_unaffected"])].copy()
    frame["is_polyp"] = (frame["disease_stage_group"] == "polyp").astype(int)
    for module in ["wnt_stem", "wnt_core_ihc"]:
        col = f"score__{module}"
        model_specs = [
            ("base", col, f"{col} ~ is_polyp + score__proliferation_control + log10_tss_total_counts"),
            ("housekeeping_adjusted", col, f"{col} ~ is_polyp + score__proliferation_control + score__housekeeping_control + log10_tss_total_counts"),
            (
                "wnt_minus_housekeeping_contrast",
                f"contrast__{module}_minus_housekeeping",
                f"contrast__{module}_minus_housekeeping ~ is_polyp + score__proliferation_control + log10_tss_total_counts",
            ),
        ]
        for model_type, outcome, formula in model_specs:
            required = [outcome, "is_polyp", "score__proliferation_control", "log10_tss_total_counts"]
            if "score__housekeeping_control" in formula:
                required.append("score__housekeeping_control")
            data = frame[required].dropna()
            if len(data) < 10:
                continue
            model = smf.ols(formula, data=data).fit(cov_type="HC3")
            rows.append(
                {
                    "outcome": outcome,
                    "term": "is_polyp",
                    "model_type": model_type,
                    "n": int(model.nobs),
                    "coef": float(model.params.get("is_polyp", np.nan)),
                    "se_hc3": float(model.bse.get("is_polyp", np.nan)),
                    "p_value": float(model.pvalues.get("is_polyp", np.nan)),
                    "r_squared": float(model.rsquared),
                    "model": formula,
                    "note": "FAP status not included because polyp status is near-perfectly confounded with FAP in this subset.",
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_adj_bh"] = multipletests(out["p_value"], method="fdr_bh")[1]
    return out


def write_summary(selected: pd.DataFrame, tests: pd.DataFrame, paired: pd.DataFrame, models: pd.DataFrame, availability: pd.DataFrame) -> None:
    lines = [
        "Becker multiome TSS accessibility feasibility analysis",
        "=" * 58,
        "",
        f"Selected samples: {len(selected)}",
        selected["disease_stage_group"].value_counts().to_string(),
        "",
        f"TSS window: +/-{WINDOW_BP} bp around GENCODE v44 gene TSS.",
        "Normalization: per-gene log1p(fragment count / total assayed TSS-window fragments * 1e6), then module mean and sample z-score.",
        "Caveat: promoter/TSS-window accessibility only; no peak calling, motif deviation, or cell-type-resolved ATAC clustering.",
        "",
        "Module gene availability:",
        availability.to_string(index=False),
        "",
        "Primary all-sample polyp vs normal/unaffected tests:",
    ]
    primary = tests.loc[tests["comparison"] == "all_polyp_vs_normal_unaffected"].copy()
    if not primary.empty:
        lines.append(primary[["module", "n_a", "n_b", "delta_median_a_minus_b", "p_value", "p_adj_bh_within_comparison"]].to_string(index=False))
    lines += ["", "FAP-only polyp vs unaffected sensitivity:"]
    fap = tests.loc[tests["comparison"] == "fap_only_polyp_vs_unaffected"].copy()
    if not fap.empty:
        lines.append(fap[["module", "n_a", "n_b", "delta_median_a_minus_b", "p_value", "p_adj_bh_within_comparison"]].to_string(index=False))
    lines += ["", "Patient-paired sensitivity:"]
    if not paired.empty:
        lines.append(paired[["module", "n_pairs", "median_delta_polyp_minus_normal", "p_value", "p_adj_bh"]].to_string(index=False))
    lines += ["", "Adjusted WNT models:"]
    if not models.empty:
        lines.append(models[["outcome", "model_type", "n", "coef", "se_hc3", "p_value", "p_adj_bh", "r_squared"]].to_string(index=False))
    else:
        lines.append("No adjusted models fit.")
    (OUT_DIR / "becker_multiome_tss_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Parsing metadata and fragment tar index", flush=True)
    sample_meta = parse_series_matrix(SERIES_MATRIX)
    fragment_meta = list_fragment_members(RAW_TAR)
    selected = select_samples(fragment_meta, sample_meta)
    selected.to_csv(OUT_DIR / "becker_multiome_tss_sample_metadata.tsv", sep="\t", index=False)
    print(selected.groupby("disease_stage_group").size().to_string(), flush=True)
    print(f"Selected compressed fragment size: {selected['fragment_gb'].sum():.2f} GB", flush=True)

    extract_selected(selected)

    wanted_genes = sorted({gene for genes in GENE_SETS.values() for gene in genes})
    windows = load_tss_windows(GENCODE_GTF, wanted_genes)
    availability = pd.DataFrame(module_availability_rows(windows))
    windows.to_csv(OUT_DIR / "becker_multiome_tss_gene_windows.tsv", sep="\t", index=False)
    availability.to_csv(OUT_DIR / "becker_multiome_tss_module_gene_availability.tsv", sep="\t", index=False)

    print("Counting TSS-window fragments", flush=True)
    gene_rows: list[dict[str, object]] = []
    for i, (_, sample) in enumerate(selected.iterrows(), start=1):
        print(f"[{i}/{len(selected)}] {sample['multiome_geo_accession']} {sample['sample_core']} {sample['disease_stage_group']}", flush=True)
        gene_rows.extend(count_sample_windows(sample, windows))

    gene_counts = pd.DataFrame(gene_rows)
    scored = score_modules(selected, gene_counts)
    gene_counts = gene_counts.merge(
        scored[["multiome_geo_accession", "sample_core", "tss_total_counts"]],
        on=["multiome_geo_accession", "sample_core"],
        how="left",
    )
    gene_counts["gene_log_tss_cpm"] = np.log1p(gene_counts["fragment_count"] / gene_counts["tss_total_counts"].clip(lower=1) * 1_000_000.0)

    tests = module_tests(scored)
    paired = paired_patient_tests(scored)
    models = adjusted_models(scored)

    gene_counts.to_csv(OUT_DIR / "becker_multiome_tss_gene_counts.tsv", sep="\t", index=False)
    scored.to_csv(OUT_DIR / "becker_multiome_tss_module_scores.tsv", sep="\t", index=False)
    tests.to_csv(OUT_DIR / "becker_multiome_tss_module_tests.tsv", sep="\t", index=False)
    paired.to_csv(OUT_DIR / "becker_multiome_tss_paired_tests.tsv", sep="\t", index=False)
    models.to_csv(OUT_DIR / "becker_multiome_tss_adjusted_models.tsv", sep="\t", index=False)
    write_summary(selected, tests, paired, models, availability)
    print(f"Wrote results to {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
