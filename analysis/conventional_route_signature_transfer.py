#!/usr/bin/env python3
"""Data-driven conventional adenoma route signature and cross-cohort transfer."""

from __future__ import annotations

import csv
import gzip
import io
import re
import tarfile
import warnings
from collections import defaultdict
from io import StringIO
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from statsmodels.duration.hazard_regression import PHReg


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "route_signature"

CHEN_DIR = ROOT / "data_sources" / "Chen_Cell_2021_CELLxGENE"
CHEN_DATASETS = {
    "discovery": CHEN_DIR / "chen_discovery_epithelial.h5ad",
    "validation": CHEN_DIR / "chen_validation_epithelial.h5ad",
}

BECKER_DIR = ROOT / "data_sources" / "Becker_NatGenet_2022_GEO"
BECKER_SERIES = BECKER_DIR / "GSE201348_series_matrix.txt.gz"
BECKER_TAR = BECKER_DIR / "GSE201348_RAW_scRNA.tar"

ATLAS_PATH = ROOT / "data_sources" / "CRC_Atlas_CZI_core" / "crc_atlas_core_czi.h5ad"

GSE_DIR = ROOT / "data_sources" / "GEO_bulk_recurrence" / "GSE39582"
GPL_DIR = ROOT / "data_sources" / "GEO_bulk_recurrence" / "GPL570"
GSE39582_SERIES = GSE_DIR / "GSE39582_series_matrix.txt.gz"
GPL570_ANNOT = GPL_DIR / "GPL570.annot.gz"

ROUTE_MAP = {
    "NL": "normal",
    "TA": "conventional_adenoma",
    "TV": "conventional_adenoma",
    "TVA": "conventional_adenoma",
    "HP": "serrated",
    "SSL": "serrated",
    "UNC": "uncertain",
}

WNT_STEM_GENES = ["LGR5", "ASCL2", "OLFM4", "AXIN2", "SOX9", "EPHB2", "SMOC2"]
WNT_CORE_IHC_GENES = ["OLFM4", "SOX9", "EPHB2"]
PROLIFERATION_GENES = ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "TYMS", "UBE2C", "CENPF"]
EPITHELIAL_MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "MUC13", "TACSTD2", "CDH1"]

CELL_CYCLE_EXCLUDE = {
    "MKI67",
    "TOP2A",
    "PCNA",
    "MCM2",
    "MCM3",
    "MCM4",
    "MCM5",
    "MCM6",
    "MCM7",
    "TYMS",
    "UBE2C",
    "CENPF",
    "CDK1",
    "CCNA2",
    "CCNB1",
    "CCNB2",
    "AURKA",
    "AURKB",
    "BIRC5",
    "CDC20",
    "CDC45",
    "CDC6",
    "MKI67IP",
}

SIGNATURE_UP_N = 50
SIGNATURE_DOWN_N = 50

CARRIER_GROUPS = {
    "normal_epithelial": ("normal", "Epithelial cell"),
    "polyp_epithelial": ("polyp", "Epithelial cell"),
    "polyp_cancer": ("polyp", "Cancer cell"),
    "tumor_epithelial": ("tumor", "Epithelial cell"),
    "tumor_cancer": ("tumor", "Cancer cell"),
    "metastasis_epithelial": ("metastasis", "Epithelial cell"),
    "metastasis_cancer": ("metastasis", "Cancer cell"),
}


def decode(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


def read_categorical(group: h5py.Group) -> np.ndarray:
    codes = group["codes"][()]
    cats = np.array([decode(v) for v in group["categories"][()]], dtype=object)
    out = np.empty(len(codes), dtype=object)
    valid = codes >= 0
    out[valid] = cats[codes[valid]]
    out[~valid] = None
    return out


def read_h5_col(root: h5py.Group, key: str) -> np.ndarray:
    obj = root[key]
    if isinstance(obj, h5py.Group):
        return read_categorical(obj)
    arr = obj[()]
    if arr.dtype.kind in "SOU":
        return np.array([decode(v) for v in arr], dtype=object)
    return arr


def bh_adjust(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return out
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy()
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out.loc[order] = np.minimum(adjusted, 1.0)
    return out


def is_excluded_gene(gene: str) -> bool:
    upper = gene.upper()
    if upper.startswith("ENSG"):
        return True
    if upper in CELL_CYCLE_EXCLUDE:
        return True
    if upper.startswith(("MT-", "RPL", "RPS", "MTRNR")):
        return True
    if upper in {"MALAT1", "XIST", "NEAT1", "JUN", "FOS", "HBB", "HBA1", "HBA2"}:
        return True
    return False


def load_chen_specimen_pseudobulk(path: Path, dataset_name: str, chunk_size: int = 384) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_cache = OUT_DIR / f"chen_{dataset_name}_specimen_pseudobulk_meta.tsv"
    expr_cache = OUT_DIR / f"chen_{dataset_name}_specimen_pseudobulk_expression.tsv.gz"
    if meta_cache.exists() and expr_cache.exists():
        meta = pd.read_csv(meta_cache, sep="\t")
        expr = pd.read_csv(expr_cache, sep="\t", compression="gzip")
        return meta, expr

    with h5py.File(path, "r", rdcc_nbytes=512 * 1024 * 1024, rdcc_nslots=1_000_003) as f:
        genes = np.array(read_h5_col(f["var"], "feature_name"), dtype=object)
        specimen_id = read_h5_col(f["obs"], "HTAN Specimen ID")
        donor_id = read_h5_col(f["obs"], "donor_id")
        polyp_type = read_h5_col(f["obs"], "Polyp_Type")
        sample_classification = read_h5_col(f["obs"], "Sample_Classification")
        route_group = pd.Series(polyp_type).map(ROUTE_MAP).fillna("other").to_numpy(dtype=object)

        meta = pd.DataFrame(
            {
                "dataset": dataset_name,
                "specimen_id": specimen_id,
                "donor_id": donor_id,
                "polyp_type": polyp_type,
                "sample_classification": sample_classification,
                "route_group": route_group,
            }
        )
        specimen_meta = (
            meta.groupby(["dataset", "specimen_id"], observed=True)
            .agg(
                donor_id=("donor_id", "first"),
                polyp_type=("polyp_type", "first"),
                sample_classification=("sample_classification", "first"),
                route_group=("route_group", "first"),
                n_cells=("route_group", "size"),
            )
            .reset_index()
        )
        specimen_meta = specimen_meta[specimen_meta["n_cells"] >= 100].reset_index(drop=True)
        specimen_to_code = {sid: i for i, sid in enumerate(specimen_meta["specimen_id"])}
        keep_mask = np.array([sid in specimen_to_code for sid in specimen_id])
        kept_cells = np.flatnonzero(keep_mask)
        group_codes = np.array([specimen_to_code[sid] for sid in specimen_id[keep_mask]], dtype=np.int32)
        design = csr_matrix(
            (np.ones(len(kept_cells), dtype=np.float32), (group_codes, kept_cells)),
            shape=(len(specimen_meta), len(specimen_id)),
        )
        counts = np.asarray(design.sum(axis=1)).ravel().astype(np.float32)
        x = f["X"]
        expression = np.zeros((len(specimen_meta), len(genes)), dtype=np.float32)
        for start in range(0, len(genes), chunk_size):
            end = min(start + chunk_size, len(genes))
            block = x[:, start:end].astype(np.float32)
            expression[:, start:end] = (design @ block) / counts[:, None]
            print(f"{dataset_name}: pseudobulk genes {end:,}/{len(genes):,}", flush=True)

    expr = pd.DataFrame(expression, columns=genes)
    if expr.columns.duplicated().any():
        expr = expr.T.groupby(level=0, sort=False).mean().T
    specimen_meta.to_csv(meta_cache, sep="\t", index=False)
    expr.to_csv(expr_cache, sep="\t", index=False, compression="gzip")
    return specimen_meta, expr


def gene_level_tests(meta: pd.DataFrame, expr: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    rows = []
    conv_idx = meta.index[meta["route_group"].eq("conventional_adenoma")].to_numpy()
    normal_idx = meta.index[meta["route_group"].eq("normal")].to_numpy()
    serrated_idx = meta.index[meta["route_group"].eq("serrated")].to_numpy()
    values = expr.to_numpy(dtype=np.float32)
    for pos, gene in enumerate(expr.columns):
        conv = values[conv_idx, pos]
        normal = values[normal_idx, pos]
        serrated = values[serrated_idx, pos] if len(serrated_idx) else np.array([], dtype=np.float32)
        if np.nanmax(values[:, pos]) <= 0:
            continue
        p_conv = mannwhitneyu(conv, normal, alternative="two-sided").pvalue if len(conv) >= 3 and len(normal) >= 3 else np.nan
        p_serrated = mannwhitneyu(conv, serrated, alternative="two-sided").pvalue if len(conv) >= 3 and len(serrated) >= 3 else np.nan
        rows.append(
            {
                "dataset": dataset_name,
                "gene": gene,
                "median_conventional": float(np.nanmedian(conv)),
                "median_normal": float(np.nanmedian(normal)),
                "median_serrated": float(np.nanmedian(serrated)) if len(serrated) else np.nan,
                "delta_conventional_minus_normal": float(np.nanmedian(conv) - np.nanmedian(normal)),
                "delta_conventional_minus_serrated": float(np.nanmedian(conv) - np.nanmedian(serrated)) if len(serrated) else np.nan,
                "p_conventional_vs_normal": float(p_conv),
                "p_conventional_vs_serrated": float(p_serrated) if not pd.isna(p_serrated) else np.nan,
                "mean_expression": float(np.nanmean(values[:, pos])),
                "excluded_gene": is_excluded_gene(str(gene)),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_conventional_vs_normal"] = out.groupby("dataset")["p_conventional_vs_normal"].transform(bh_adjust)
    return out


def build_signature(stats: pd.DataFrame) -> pd.DataFrame:
    wide = stats.pivot(index="gene", columns="dataset")
    rows = []
    for gene in sorted(stats["gene"].unique()):
        try:
            disc_delta = float(wide.loc[gene, ("delta_conventional_minus_normal", "discovery")])
            val_delta = float(wide.loc[gene, ("delta_conventional_minus_normal", "validation")])
            disc_p = float(wide.loc[gene, ("p_conventional_vs_normal", "discovery")])
            val_p = float(wide.loc[gene, ("p_conventional_vs_normal", "validation")])
            disc_q = float(wide.loc[gene, ("q_conventional_vs_normal", "discovery")])
            val_q = float(wide.loc[gene, ("q_conventional_vs_normal", "validation")])
            disc_expr = float(wide.loc[gene, ("mean_expression", "discovery")])
            val_expr = float(wide.loc[gene, ("mean_expression", "validation")])
            excluded = bool(wide.loc[gene, ("excluded_gene", "discovery")]) or bool(wide.loc[gene, ("excluded_gene", "validation")])
        except Exception:
            continue
        combined_direction = np.sign(disc_delta) == np.sign(val_delta)
        score = np.sign(disc_delta + val_delta) * (abs(disc_delta) + abs(val_delta)) * (
            -np.log10(max(disc_p, 1e-300)) + -np.log10(max(val_p, 1e-300))
        )
        rows.append(
            {
                "gene": gene,
                "delta_discovery": disc_delta,
                "delta_validation": val_delta,
                "p_discovery": disc_p,
                "p_validation": val_p,
                "q_discovery": disc_q,
                "q_validation": val_q,
                "mean_expression_discovery": disc_expr,
                "mean_expression_validation": val_expr,
                "direction_consistent": combined_direction,
                "excluded_gene": excluded,
                "selection_score": float(score),
            }
        )
    combined = pd.DataFrame(rows)
    combined["selection_magnitude"] = combined["selection_score"].abs()
    combined["passes_nominal_filter"] = (combined["p_discovery"] < 0.10) & (combined["p_validation"] < 0.10)
    eligible = combined[
        combined["direction_consistent"]
        & ~combined["excluded_gene"]
        & (combined["mean_expression_discovery"] > 0.001)
        & (combined["mean_expression_validation"] > 0.001)
    ].copy()
    up = eligible[(eligible["delta_discovery"] > 0) & (eligible["delta_validation"] > 0)].nlargest(SIGNATURE_UP_N, "selection_magnitude").copy()
    down = eligible[(eligible["delta_discovery"] < 0) & (eligible["delta_validation"] < 0)].nlargest(SIGNATURE_DOWN_N, "selection_magnitude").copy()
    up["signature_direction"] = "adenoma_up"
    down["signature_direction"] = "adenoma_down"
    signature = pd.concat([up, down], ignore_index=True)
    signature["rank_within_direction"] = signature.groupby("signature_direction")["selection_magnitude"].rank(
        ascending=False,
        method="first",
    )
    return signature.sort_values(["signature_direction", "rank_within_direction", "gene"])


def zscore_frame(expr: pd.DataFrame) -> pd.DataFrame:
    values = expr.astype(float)
    sd = values.std(axis=0, skipna=True).replace(0, np.nan)
    return values.sub(values.mean(axis=0, skipna=True), axis=1).div(sd, axis=1)


def route_score_from_expression(expr: pd.DataFrame, up_genes: list[str], down_genes: list[str], prefix: str = "score__") -> pd.DataFrame:
    present_up = [gene for gene in up_genes if gene in expr.columns]
    present_down = [gene for gene in down_genes if gene in expr.columns]
    z = zscore_frame(expr[present_up + present_down])
    out = pd.DataFrame(index=expr.index)
    out[f"{prefix}ca_route_up"] = z[present_up].mean(axis=1) if present_up else np.nan
    out[f"{prefix}ca_route_down"] = z[present_down].mean(axis=1) if present_down else np.nan
    out[f"{prefix}ca_route_signature"] = out[f"{prefix}ca_route_up"] - out[f"{prefix}ca_route_down"]
    out[f"{prefix}ca_route_n_up_present"] = len(present_up)
    out[f"{prefix}ca_route_n_down_present"] = len(present_down)
    return out


def strip_score_prefix(col: str) -> str:
    if col.startswith("score__"):
        return col.removeprefix("score__")
    for prefix in ["score_all__", "score_epi__"]:
        if col.startswith(prefix):
            return col.removeprefix("score_")
    return col


def module_score(expr: pd.DataFrame, genes: list[str], name: str) -> pd.Series:
    present = [gene for gene in genes if gene in expr.columns]
    if not present:
        return pd.Series(np.nan, index=expr.index, name=name)
    z = zscore_frame(expr[present])
    return z[present].mean(axis=1).rename(name)


def compare_groups(frame: pd.DataFrame, group_col: str, score_cols: list[str], comparisons: list[tuple[str, str, str]]) -> pd.DataFrame:
    rows = []
    for label, group_a, group_b in comparisons:
        for col in score_cols:
            aval = frame.loc[frame[group_col].eq(group_a), col].dropna()
            bval = frame.loc[frame[group_col].eq(group_b), col].dropna()
            if len(aval) < 3 or len(bval) < 3:
                continue
            p_value = mannwhitneyu(aval, bval, alternative="two-sided").pvalue
            rows.append(
                {
                    "comparison": label,
                    "score": strip_score_prefix(col),
                    "group_a": group_a,
                    "group_b": group_b,
                    "n_a": len(aval),
                    "n_b": len(bval),
                    "median_a": float(aval.median()),
                    "median_b": float(bval.median()),
                    "delta_a_minus_b": float(aval.median() - bval.median()),
                    "p_mannwhitney": float(p_value),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_within_comparison"] = out.groupby("comparison")["p_mannwhitney"].transform(bh_adjust)
        out = out.sort_values(["comparison", "p_mannwhitney", "score"])
    return out


def chen_signature_scores(
    meta_by_dataset: dict[str, pd.DataFrame],
    expr_by_dataset: dict[str, pd.DataFrame],
    signature: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    up_genes = signature.loc[signature["signature_direction"].eq("adenoma_up"), "gene"].tolist()
    down_genes = signature.loc[signature["signature_direction"].eq("adenoma_down"), "gene"].tolist()
    frames = []
    for dataset, meta in meta_by_dataset.items():
        expr = expr_by_dataset[dataset].copy()
        score = route_score_from_expression(expr, up_genes, down_genes)
        score["score__wnt_stem"] = module_score(expr, WNT_STEM_GENES, "score__wnt_stem")
        score["score__wnt_core_ihc"] = module_score(expr, WNT_CORE_IHC_GENES, "score__wnt_core_ihc")
        score["score__proliferation_control"] = module_score(expr, PROLIFERATION_GENES, "score__proliferation_control")
        frame = pd.concat([meta.reset_index(drop=True), score.reset_index(drop=True)], axis=1)
        frames.append(frame)
    all_scores = pd.concat(frames, ignore_index=True)
    tests = compare_groups(
        all_scores[all_scores["route_group"].isin(["normal", "conventional_adenoma", "serrated"])],
        "route_group",
        ["score__ca_route_signature", "score__ca_route_up", "score__ca_route_down", "score__wnt_stem", "score__proliferation_control"],
        [
            ("conventional_vs_normal", "conventional_adenoma", "normal"),
            ("serrated_vs_normal", "serrated", "normal"),
            ("conventional_vs_serrated", "conventional_adenoma", "serrated"),
        ],
    )

    paired_rows = []
    donor_route = (
        all_scores.groupby(["dataset", "donor_id", "route_group"], observed=True)
        [["score__ca_route_signature", "score__wnt_stem", "score__proliferation_control"]]
        .median()
        .reset_index()
    )
    donor_route["pair_id"] = donor_route["dataset"].astype(str) + "::" + donor_route["donor_id"].astype(str)
    for dataset, frame in list(donor_route.groupby("dataset", observed=True)) + [("combined", donor_route)]:
        subset = frame[frame["route_group"].isin(["conventional_adenoma", "normal"])]
        wide = subset.pivot(index="pair_id", columns="route_group", values=["score__ca_route_signature", "score__wnt_stem", "score__proliferation_control"])
        for col in ["score__ca_route_signature", "score__wnt_stem", "score__proliferation_control"]:
            if (col, "conventional_adenoma") not in wide.columns or (col, "normal") not in wide.columns:
                continue
            pair = wide[[(col, "conventional_adenoma"), (col, "normal")]].dropna()
            if len(pair) < 3:
                continue
            diff = pair[(col, "conventional_adenoma")] - pair[(col, "normal")]
            paired_rows.append(
                {
                    "dataset": dataset,
                    "score": col.replace("score__", ""),
                    "n_pairs": len(diff),
                    "median_delta_conventional_minus_normal": float(diff.median()),
                    "p_wilcoxon": float(wilcoxon(diff, zero_method="wilcox", alternative="two-sided").pvalue),
                }
            )
    paired = pd.DataFrame(paired_rows)
    return all_scores, tests, paired


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


def parse_becker_series(path: Path) -> pd.DataFrame:
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
    accessions = rows["Sample_geo_accession"][0]
    titles = rows["Sample_title"][0]
    meta = pd.DataFrame({"geo_accession": accessions, "sample_title": titles})
    for values in rows.get("Sample_characteristics_ch1", []):
        parsed = [v.split(":", 1) if ":" in v else [v, ""] for v in values]
        keys = [normalize_column(p[0]) for p in parsed]
        if not keys:
            continue
        col = pd.Series(keys).mode().iloc[0]
        meta[col] = [p[1].strip() if len(p) > 1 else "" for p in parsed]
    meta["sample_name_from_title"] = (
        meta["sample_title"].str.replace(", snRNAseq", "", regex=False).str.replace(", Replicate1", "", regex=False).str.replace(", Replicate2", "", regex=False).str.strip()
    )
    meta["patient_id"] = meta["sample_name_from_title"].map(patient_id_from_sample_name)
    meta["disease_stage_group"] = meta["disease_stage"].map({"Normal": "normal_unaffected", "Unaffected": "normal_unaffected", "Polyp": "polyp", "CRC": "crc"})
    meta["familial_adenomatous_polyposis"] = meta["familial_adenomatous_polyposis"].replace({"": "unknown"}).fillna("unknown")
    meta["sex"] = meta["sex"].replace({"": "unknown"}).fillna("unknown")
    return meta


def becker_tar_prefixes(tar_path: Path) -> pd.DataFrame:
    rows = []
    with tarfile.open(tar_path, "r") as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            if not name.endswith("_matrix.mtx.gz"):
                continue
            prefix = name.removesuffix("_matrix.mtx.gz")
            rows.append({"geo_accession": prefix.split("_", 1)[0], "file_prefix": prefix})
    return pd.DataFrame(rows).sort_values("geo_accession")


def read_features(tar: tarfile.TarFile, member_name: str) -> list[str]:
    with tar.extractfile(member_name) as raw:
        if raw is None:
            raise FileNotFoundError(member_name)
        with gzip.GzipFile(fileobj=raw) as gz:
            genes = []
            for line in io.TextIOWrapper(gz, encoding="utf-8", errors="replace"):
                parts = line.rstrip("\n").split("\t")
                genes.append(parts[1] if len(parts) > 1 else parts[0])
    return genes


def read_matrix(tar: tarfile.TarFile, member_name: str):
    with tar.extractfile(member_name) as raw:
        if raw is None:
            raise FileNotFoundError(member_name)
        with gzip.GzipFile(fileobj=raw) as gz:
            data = io.BytesIO(gz.read())
    return mmread(data).tocsr()


def log_cpm(counts: np.ndarray, library_size: float) -> np.ndarray:
    if library_size <= 0:
        return np.full(counts.shape, np.nan)
    return np.log1p(counts / library_size * 1_000_000.0)


def score_becker(signature: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    up_genes = signature.loc[signature["signature_direction"].eq("adenoma_up"), "gene"].tolist()
    down_genes = signature.loc[signature["signature_direction"].eq("adenoma_down"), "gene"].tolist()
    wanted = sorted(set(up_genes + down_genes + WNT_STEM_GENES + WNT_CORE_IHC_GENES + PROLIFERATION_GENES + EPITHELIAL_MARKERS))
    meta = parse_becker_series(BECKER_SERIES).merge(becker_tar_prefixes(BECKER_TAR), on="geo_accession", how="left")
    rows = []
    availability = []
    with tarfile.open(BECKER_TAR, "r") as tar:
        for idx, row in meta.sort_values("geo_accession").reset_index(drop=True).iterrows():
            prefix = row["file_prefix"]
            genes = read_features(tar, f"{prefix}_features.tsv.gz")
            gene_to_indices: dict[str, list[int]] = defaultdict(list)
            for pos, gene in enumerate(genes):
                gene_to_indices[str(gene)].append(pos)
            present = [gene for gene in wanted if gene in gene_to_indices]
            missing = [gene for gene in wanted if gene not in gene_to_indices]
            availability.append({"geo_accession": row["geo_accession"], "n_present": len(present), "n_missing": len(missing), "missing_genes": ",".join(missing)})
            matrix = read_matrix(tar, f"{prefix}_matrix.mtx.gz")
            total_counts_per_cell = np.asarray(matrix.sum(axis=0)).ravel()
            total_all = float(total_counts_per_cell.sum())
            epithelial_idx = [i for gene in EPITHELIAL_MARKERS for i in gene_to_indices.get(gene, [])]
            epithelial_mask = np.asarray(matrix[epithelial_idx, :].sum(axis=0)).ravel() > 0 if epithelial_idx else np.zeros(matrix.shape[1], dtype=bool)
            total_epi = float(total_counts_per_cell[epithelial_mask].sum()) if epithelial_mask.any() else 0.0
            selected_idx = [gene_to_indices[gene][0] for gene in present]
            selected_counts_all = np.asarray(matrix[selected_idx, :].sum(axis=1)).ravel() if selected_idx else np.array([])
            selected_counts_epi = np.asarray(matrix[selected_idx, :][:, epithelial_mask].sum(axis=1)).ravel() if selected_idx and epithelial_mask.any() else np.full(len(selected_idx), np.nan)
            out = {
                "geo_accession": row["geo_accession"],
                "n_nuclei": int(matrix.shape[1]),
                "n_epithelial_marker_positive": int(epithelial_mask.sum()),
                "epithelial_marker_positive_fraction": float(epithelial_mask.mean()) if matrix.shape[1] else np.nan,
            }
            for gene, value in zip(present, log_cpm(selected_counts_all, total_all)):
                out[f"all__{gene}"] = float(value)
            for gene, value in zip(present, log_cpm(selected_counts_epi, total_epi)):
                out[f"epi__{gene}"] = float(value)
            rows.append(out)
            print(f"Becker signature scoring {idx + 1:02d}/{len(meta)} {row['geo_accession']}", flush=True)
    expr = pd.DataFrame(rows)
    frame = meta.merge(expr, on="geo_accession", how="left")
    for scope in ["all", "epi"]:
        cols = [c for c in frame.columns if c.startswith(f"{scope}__")]
        gene_expr = frame[cols].copy()
        gene_expr.columns = [c.split("__", 1)[1] for c in gene_expr.columns]
        scores = route_score_from_expression(gene_expr, up_genes, down_genes, prefix=f"score_{scope}__")
        scores[f"score_{scope}__wnt_stem"] = module_score(gene_expr, WNT_STEM_GENES, f"score_{scope}__wnt_stem")
        scores[f"score_{scope}__wnt_core_ihc"] = module_score(gene_expr, WNT_CORE_IHC_GENES, f"score_{scope}__wnt_core_ihc")
        scores[f"score_{scope}__proliferation_control"] = module_score(gene_expr, PROLIFERATION_GENES, f"score_{scope}__proliferation_control")
        frame = pd.concat([frame, scores.reset_index(drop=True)], axis=1)
    return frame, pd.DataFrame(availability)


def becker_tests_and_models(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_cols = [
        "score_all__ca_route_signature",
        "score_epi__ca_route_signature",
        "score_all__wnt_stem",
        "score_epi__wnt_stem",
        "score_all__proliferation_control",
        "score_epi__proliferation_control",
    ]
    tests = compare_groups(
        frame,
        "disease_stage_group",
        score_cols,
        [
            ("polyp_vs_normal_unaffected", "polyp", "normal_unaffected"),
            ("crc_vs_normal_unaffected", "crc", "normal_unaffected"),
            ("crc_vs_polyp", "crc", "polyp"),
        ],
    )
    rows = []
    frame = frame.copy()
    frame["log10_n_nuclei"] = np.log10(frame["n_nuclei"].clip(lower=1))
    for outcome, prolif in [
        ("score_all__ca_route_signature", "score_all__proliferation_control"),
        ("score_epi__ca_route_signature", "score_epi__proliferation_control"),
    ]:
        mf = frame[[outcome, "disease_stage_group", "familial_adenomatous_polyposis", "sex", prolif, "log10_n_nuclei", "epithelial_marker_positive_fraction"]].replace([np.inf, -np.inf], np.nan).dropna().copy()
        mf = mf[mf["disease_stage_group"].isin(["normal_unaffected", "polyp", "crc"])]
        if len(mf) < 20:
            continue
        mf["disease_stage_group"] = pd.Categorical(mf["disease_stage_group"], categories=["normal_unaffected", "polyp", "crc"])
        y = mf[outcome].astype(float)
        x = pd.get_dummies(mf[["disease_stage_group", "familial_adenomatous_polyposis", "sex"]], drop_first=True, dtype=float)
        for col in [prolif, "log10_n_nuclei", "epithelial_marker_positive_fraction"]:
            x[col] = mf[col].astype(float)
        x = sm.add_constant(x.astype(float), has_constant="add")
        fit = sm.OLS(y, x).fit(cov_type="HC1")
        for term in x.columns:
            if term == "const":
                continue
            rows.append({"outcome": outcome.replace("score_", ""), "term": term, "n": len(mf), "coef": float(fit.params[term]), "se_hc1": float(fit.bse[term]), "p_value": float(fit.pvalues[term]), "r_squared": float(fit.rsquared)})
    models = pd.DataFrame(rows)
    if not models.empty:
        models["q_value"] = bh_adjust(models["p_value"])
        models = models.sort_values(["outcome", "p_value", "term"])
    return tests, models


def atlas_carrier_group(sample_type: np.ndarray, coarse_type: np.ndarray) -> np.ndarray:
    out = np.full(len(sample_type), "not_used", dtype=object)
    for name, (sample, coarse) in CARRIER_GROUPS.items():
        out[(sample_type == sample) & (coarse_type == coarse)] = name
    return out


def mode_or_missing(values: pd.Series) -> str:
    values = values.dropna().astype(str)
    values = values[~values.isin(["", "None", "nan", "missing"])]
    if values.empty:
        return "missing"
    return str(values.value_counts().index[0])


def sample_atlas_rows(groups: np.ndarray, max_per_group: int = 20000, seed: int = 20260704) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected = []
    for group_name in CARRIER_GROUPS:
        idx = np.flatnonzero(groups == group_name)
        if len(idx) > max_per_group:
            idx = rng.choice(idx, size=max_per_group, replace=False)
        if len(idx):
            selected.append(idx)
    return np.sort(np.concatenate(selected).astype(np.int64))


def extract_atlas_selected_expression(f: h5py.File, rows: np.ndarray, gene_indices: np.ndarray) -> np.ndarray:
    x = f["X"]
    indptr = x["indptr"][()]
    indices_ds = x["indices"]
    data_ds = x["data"]
    n_vars = int(x.attrs["shape"][1])
    lookup = np.full(n_vars, -1, dtype=np.int16)
    lookup[gene_indices] = np.arange(len(gene_indices), dtype=np.int16)
    out = np.zeros((len(rows), len(gene_indices)), dtype=np.float32)
    for out_pos, row in enumerate(rows):
        start = int(indptr[row])
        end = int(indptr[row + 1])
        row_indices = indices_ds[start:end]
        selected_pos = lookup[row_indices]
        mask = selected_pos >= 0
        if mask.any():
            out[out_pos, selected_pos[mask]] = data_ds[start:end][mask]
        if (out_pos + 1) % 20000 == 0:
            print(f"Atlas signature extraction {out_pos + 1:,}/{len(rows):,}", flush=True)
    return out


def score_atlas(signature: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    up_genes = signature.loc[signature["signature_direction"].eq("adenoma_up"), "gene"].tolist()
    down_genes = signature.loc[signature["signature_direction"].eq("adenoma_down"), "gene"].tolist()
    wanted = sorted(set(up_genes + down_genes + WNT_STEM_GENES + WNT_CORE_IHC_GENES + PROLIFERATION_GENES))
    with h5py.File(ATLAS_PATH, "r", rdcc_nbytes=512 * 1024 * 1024, rdcc_nslots=1_000_003) as f:
        sample_type = read_h5_col(f["obs"], "sample_type")
        coarse_type = read_h5_col(f["obs"], "cell_type_coarse_crc_atlas")
        groups = atlas_carrier_group(sample_type, coarse_type)
        rows = sample_atlas_rows(groups)
        gene_names = read_h5_col(f["var"], "feature_name")
        gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}
        present = [gene for gene in wanted if gene in gene_to_idx]
        gene_indices = np.array([gene_to_idx[gene] for gene in present], dtype=np.int64)
        expr = extract_atlas_selected_expression(f, rows, gene_indices)
        cells = pd.DataFrame(
            {
                "row_index": rows,
                "carrier_group": groups[rows],
                "sample_type": sample_type[rows],
                "cell_type_coarse_crc_atlas": coarse_type[rows],
                "cell_type_middle_crc_atlas": read_h5_col(f["obs"], "cell_type_middle_crc_atlas")[rows],
                "cell_type_fine_crc_atlas": read_h5_col(f["obs"], "cell_type_fine_crc_atlas")[rows],
                "donor_id": read_h5_col(f["obs"], "donor_id")[rows],
                "study_id": read_h5_col(f["obs"], "study_id")[rows],
                "microsatellite_status": read_h5_col(f["obs"], "microsatellite_status")[rows],
                "CMS_type": read_h5_col(f["obs"], "CMS_type")[rows],
            }
        )
    gene_expr = pd.DataFrame(expr, columns=present)
    scores = route_score_from_expression(gene_expr, up_genes, down_genes)
    scores["score__wnt_stem"] = module_score(gene_expr, WNT_STEM_GENES, "score__wnt_stem")
    scores["score__wnt_core_ihc"] = module_score(gene_expr, WNT_CORE_IHC_GENES, "score__wnt_core_ihc")
    scores["score__proliferation_control"] = module_score(gene_expr, PROLIFERATION_GENES, "score__proliferation_control")
    cells = pd.concat([cells, scores], axis=1)
    score_cols = [c for c in cells.columns if c.startswith("score__")]
    donors = (
        cells.groupby(["donor_id", "carrier_group"], observed=True)
        .agg(
            sample_type=("sample_type", "first"),
            cell_type_coarse_crc_atlas=("cell_type_coarse_crc_atlas", "first"),
            study_id=("study_id", lambda x: ";".join(sorted(set(map(str, x))))),
            microsatellite_status=("microsatellite_status", mode_or_missing),
            CMS_type=("CMS_type", mode_or_missing),
            n_cells_sampled=("row_index", "size"),
            **{col: (col, "median") for col in score_cols},
        )
        .reset_index()
    )
    comparisons = [
        ("polyp_epithelial_vs_normal_epithelial", "polyp_epithelial", "normal_epithelial"),
        ("polyp_cancer_vs_normal_epithelial", "polyp_cancer", "normal_epithelial"),
        ("tumor_cancer_vs_normal_epithelial", "tumor_cancer", "normal_epithelial"),
        ("metastasis_cancer_vs_normal_epithelial", "metastasis_cancer", "normal_epithelial"),
        ("tumor_cancer_vs_polyp_epithelial", "tumor_cancer", "polyp_epithelial"),
        ("metastasis_cancer_vs_tumor_cancer", "metastasis_cancer", "tumor_cancer"),
    ]
    tests = compare_groups(donors[donors["n_cells_sampled"] >= 20], "carrier_group", score_cols, comparisons)
    model_rows = []
    mf = donors[donors["n_cells_sampled"] >= 20].copy()
    mf["carrier_group"] = pd.Categorical(mf["carrier_group"], categories=list(CARRIER_GROUPS))
    mf["log10_n_cells_sampled"] = np.log10(mf["n_cells_sampled"].clip(lower=1))
    for outcome in ["score__ca_route_signature", "score__wnt_stem"]:
        model_frame = mf[[outcome, "carrier_group", "study_id", "score__proliferation_control", "log10_n_cells_sampled"]].dropna().copy()
        y = model_frame[outcome].astype(float)
        x = pd.get_dummies(model_frame[["carrier_group", "study_id"]], drop_first=True, dtype=float)
        x["score__proliferation_control"] = model_frame["score__proliferation_control"].astype(float)
        x["log10_n_cells_sampled"] = model_frame["log10_n_cells_sampled"].astype(float)
        x = sm.add_constant(x.astype(float), has_constant="add")
        fit = sm.OLS(y, x).fit(cov_type="HC1")
        for term in x.columns:
            if term == "const":
                continue
            model_rows.append({"outcome": outcome.replace("score__", ""), "term": term, "n": len(model_frame), "coef": float(fit.params[term]), "se_hc1": float(fit.bse[term]), "p_value": float(fit.pvalues[term]), "r_squared": float(fit.rsquared)})
    models = pd.DataFrame(model_rows)
    if not models.empty:
        models["q_value"] = bh_adjust(models["p_value"])
        models = models.sort_values(["outcome", "p_value", "term"])
    return donors, tests, models


def clean_cell(value: object) -> str | None:
    text = str(value).strip().strip('"')
    if text in {"", "N/A", "nan", "None"}:
        return None
    return text


def clean_key(value: str) -> str:
    return value.strip().lower().replace(" ", "_").replace(".", "_").replace("(", "").replace(")", "").replace("/", "_")


def parse_gse39582_metadata(path: Path) -> pd.DataFrame:
    rows: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                rows.append(next(csv.reader([line], delimiter="\t")))
    accessions = None
    meta: dict[str, list[str | None]] = {}
    characteristics: dict[str, list[str | None]] = defaultdict(list)
    for row in rows:
        field = row[0]
        values = [clean_cell(v) for v in row[1:]]
        if field == "!Sample_geo_accession":
            accessions = values
        elif field == "!Sample_title":
            meta["sample_title"] = values
        elif field == "!Sample_source_name_ch1":
            meta["source_name"] = values
        elif field == "!Sample_characteristics_ch1":
            for value in values:
                if value and ":" in value:
                    key, val = value.split(":", 1)
                    characteristics[clean_key(key)].append(clean_cell(val))
    if accessions is None:
        raise ValueError("No GSE39582 accessions found.")
    frame = pd.DataFrame({"geo_accession": accessions})
    for key, values in meta.items():
        frame[key] = values
    for key, values in characteristics.items():
        if len(values) == len(frame):
            frame[key] = values
    frame["age_at_diagnosis_year"] = pd.to_numeric(frame.get("age_at_diagnosis_year"), errors="coerce")
    frame["rfs_event"] = pd.to_numeric(frame.get("rfs_event"), errors="coerce")
    frame["rfs_delay"] = pd.to_numeric(frame.get("rfs_delay"), errors="coerce")
    frame["stage_group"] = frame["tnm_stage"].map({"0": "0", "1": "I", "2": "II", "3": "III", "4": "IV"})
    frame["tumor_status"] = np.where(frame["dataset"].isin(["discovery", "validation"]), "tumor", "non_tumoral")
    frame["is_independent_sample"] = frame["dependancy_sample"].isna().astype(int)
    frame["analysis_base"] = (
        frame["tumor_status"].eq("tumor")
        & frame["is_independent_sample"].eq(1)
        & frame["rfs_event"].isin([0, 1])
        & frame["rfs_delay"].notna()
    ).astype(int)
    return frame


def read_gpl_mapping(path: Path, target_genes: set[str]) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    table_lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                table_lines.append(line)
    annot = pd.read_csv(StringIO("".join(table_lines)), sep="\t", dtype=str)
    probe_to_symbols: dict[str, set[str]] = {}
    rows = []
    for _, row in annot[["ID", "Gene symbol"]].dropna().iterrows():
        probe = str(row["ID"])
        symbols = {symbol.strip() for symbol in str(row["Gene symbol"]).split("///") if symbol.strip()}
        matched = sorted(symbols & target_genes)
        if matched:
            probe_to_symbols[probe] = symbols
            rows.append({"probe_id": probe, "matched_genes": ",".join(matched), "all_gene_symbols": "///".join(sorted(symbols))})
    return pd.DataFrame(rows), probe_to_symbols


def score_gse39582(signature: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    up_genes = signature.loc[signature["signature_direction"].eq("adenoma_up"), "gene"].tolist()
    down_genes = signature.loc[signature["signature_direction"].eq("adenoma_down"), "gene"].tolist()
    wanted = sorted(set(up_genes + down_genes + WNT_STEM_GENES + WNT_CORE_IHC_GENES + PROLIFERATION_GENES))
    metadata = parse_gse39582_metadata(GSE39582_SERIES)
    mapping_rows, probe_to_symbols = read_gpl_mapping(GPL570_ANNOT, set(wanted))
    expr = pd.read_csv(GSE39582_SERIES, sep="\t", compression="gzip", comment="!", index_col=0)
    expr.index = expr.index.astype(str)
    candidate_probes = [probe for probe in probe_to_symbols if probe in expr.index]
    candidate = expr.loc[candidate_probes].astype(float)
    reference_samples = [sample for sample in metadata.loc[metadata["analysis_base"].eq(1), "geo_accession"] if sample in candidate.columns]
    selected_rows = []
    gene_expr_rows = []
    for gene in wanted:
        probes = [probe for probe in candidate.index if gene in probe_to_symbols.get(probe, set())]
        if not probes:
            continue
        variances = candidate.loc[probes, reference_samples].var(axis=1)
        selected_probe = str(variances.sort_values(ascending=False).index[0])
        selected_rows.append({"gene": gene, "selected_probe_id": selected_probe, "n_candidate_probes": len(probes), "selected_probe_variance": float(variances.loc[selected_probe])})
        gene_expr_rows.append(candidate.loc[selected_probe].rename(gene))
    selected = pd.DataFrame(selected_rows)
    gene_expr = pd.DataFrame(gene_expr_rows)
    ref = [sample for sample in metadata.loc[metadata["analysis_base"].eq(1), "geo_accession"] if sample in gene_expr.columns]
    means = gene_expr[ref].mean(axis=1)
    sds = gene_expr[ref].std(axis=1).replace(0, np.nan)
    z = gene_expr.sub(means, axis=0).div(sds, axis=0).T
    scores = metadata.copy()
    route = route_score_from_expression(z, up_genes, down_genes)
    scores = pd.concat([scores, route.reindex(scores["geo_accession"]).reset_index(drop=True)], axis=1)
    for genes, name in [(WNT_STEM_GENES, "score__wnt_stem"), (WNT_CORE_IHC_GENES, "score__wnt_core_ihc"), (PROLIFERATION_GENES, "score__proliferation_control")]:
        scores[name] = module_score(z, genes, name).reindex(scores["geo_accession"]).to_numpy()
    event = gse_event_associations(scores)
    cox = gse_cox_models(scores)
    return scores, mapping_rows, selected, pd.concat({"event_associations": event, "cox_models": cox}, names=["result_type"]).reset_index(level=0)


def gse_event_associations(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    subsets = {
        "all_independent_tumors": frame["analysis_base"].eq(1),
        "stage_ii_iii_independent": frame["analysis_base"].eq(1) & frame["stage_group"].isin(["II", "III"]),
        "pmmr_independent": frame["analysis_base"].eq(1) & frame["mmr_status"].eq("pMMR"),
        "stage_ii_iii_pmmr_independent": frame["analysis_base"].eq(1) & frame["stage_group"].isin(["II", "III"]) & frame["mmr_status"].eq("pMMR"),
    }
    score_cols = ["score__ca_route_signature", "score__ca_route_up", "score__ca_route_down", "score__wnt_stem", "score__proliferation_control"]
    for subset_name, mask in subsets.items():
        subset = frame[mask & frame["rfs_event"].isin([0, 1])].copy()
        for col in score_cols:
            yes = subset.loc[subset["rfs_event"].eq(1), col].dropna()
            no = subset.loc[subset["rfs_event"].eq(0), col].dropna()
            if len(yes) < 8 or len(no) < 8:
                continue
            p_value = mannwhitneyu(yes, no, alternative="two-sided").pvalue
            rows.append({"subset": subset_name, "score": col.replace("score__", ""), "n_event": len(yes), "n_no_event": len(no), "median_event": float(yes.median()), "median_no_event": float(no.median()), "delta_event_minus_no_event": float(yes.median() - no.median()), "p_mannwhitney": float(p_value)})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_within_subset"] = out.groupby("subset")["p_mannwhitney"].transform(bh_adjust)
        out = out.sort_values(["subset", "p_mannwhitney", "score"])
    return out


def gse_design_matrix(frame: pd.DataFrame, score_col: str) -> pd.DataFrame:
    numeric = frame[[score_col, "age_at_diagnosis_year"]].astype(float).rename(columns={score_col: "score_sd"})
    sd = numeric["score_sd"].std()
    numeric["score_sd"] = (numeric["score_sd"] - numeric["score_sd"].mean()) / sd if sd and not np.isnan(sd) else np.nan
    categorical = pd.get_dummies(frame[["stage_group", "mmr_status", "sex", "tumor_location", "chemotherapy_adjuvant", "dataset"]], drop_first=True, dtype=float)
    design = pd.concat([numeric, categorical], axis=1)
    return design.loc[:, design.nunique(dropna=True) > 1].astype(float)


def gse_cox_models(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    subsets = {
        "all_independent_tumors": frame["analysis_base"].eq(1),
        "stage_ii_iii_independent": frame["analysis_base"].eq(1) & frame["stage_group"].isin(["II", "III"]),
        "pmmr_independent": frame["analysis_base"].eq(1) & frame["mmr_status"].eq("pMMR"),
        "stage_ii_iii_pmmr_independent": frame["analysis_base"].eq(1) & frame["stage_group"].isin(["II", "III"]) & frame["mmr_status"].eq("pMMR"),
    }
    score_cols = ["score__ca_route_signature", "score__ca_route_up", "score__ca_route_down", "score__wnt_stem", "score__proliferation_control"]
    for subset_name, mask in subsets.items():
        for col in score_cols:
            needed = ["rfs_event", "rfs_delay", col, "age_at_diagnosis_year", "stage_group", "mmr_status", "sex", "tumor_location", "chemotherapy_adjuvant", "dataset"]
            mf = frame.loc[mask, needed].replace({"N/A": np.nan, "": np.nan}).dropna().copy()
            mf = mf[(mf["rfs_delay"] > 0) & mf["rfs_event"].isin([0, 1])]
            if len(mf) < 80 or mf["rfs_event"].nunique() < 2:
                continue
            design = gse_design_matrix(mf, col).dropna()
            mf = mf.loc[design.index]
            if "score_sd" not in design.columns:
                continue
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fit = PHReg(mf["rfs_delay"].astype(float).to_numpy(), design.to_numpy(float), status=mf["rfs_event"].astype(int).to_numpy(), ties="breslow").fit(disp=False)
                params = pd.Series(fit.params, index=design.columns)
                pvals = pd.Series(fit.pvalues, index=design.columns)
                ci = pd.DataFrame(fit.conf_int(), index=design.columns, columns=["low", "high"])
                beta = params["score_sd"]
                rows.append({"subset": subset_name, "score": col.replace("score__", ""), "n": len(mf), "n_event": int(mf["rfs_event"].sum()), "hazard_ratio_per_1sd": float(np.exp(beta)), "ci_low": float(np.exp(ci.loc["score_sd", "low"])), "ci_high": float(np.exp(ci.loc["score_sd", "high"])), "p_value": float(pvals["score_sd"])})
            except Exception as exc:
                rows.append({"subset": subset_name, "score": col.replace("score__", ""), "n": len(mf), "n_event": int(mf["rfs_event"].sum()), "hazard_ratio_per_1sd": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan, "fit_error": str(exc)})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_within_subset"] = out.groupby("subset")["p_value"].transform(bh_adjust)
        out = out.sort_values(["subset", "p_value", "score"], na_position="last")
    return out


def write_summary(
    path: Path,
    signature: pd.DataFrame,
    chen_tests: pd.DataFrame,
    chen_paired: pd.DataFrame,
    becker_tests: pd.DataFrame,
    becker_models: pd.DataFrame,
    atlas_tests: pd.DataFrame,
    atlas_models: pd.DataFrame,
    gse_results: pd.DataFrame,
) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("Conventional adenoma data-driven route signature transfer\n")
        fh.write("=" * 70 + "\n\n")
        fh.write("Signature composition:\n")
        fh.write(signature["signature_direction"].value_counts().to_string())
        fh.write("\n\nTop adenoma-up genes:\n")
        fh.write(", ".join(signature.loc[signature["signature_direction"].eq("adenoma_up"), "gene"].head(20)))
        fh.write("\n\nTop adenoma-down genes:\n")
        fh.write(", ".join(signature.loc[signature["signature_direction"].eq("adenoma_down"), "gene"].head(20)))
        fh.write("\n\nChen focused tests:\n")
        fh.write(chen_tests[chen_tests["score"].isin(["ca_route_signature", "wnt_stem", "proliferation_control"])].to_string(index=False))
        fh.write("\n\nChen paired donor tests:\n")
        fh.write(chen_paired.to_string(index=False))
        fh.write("\n\nBecker focused tests:\n")
        fh.write(becker_tests[becker_tests["score"].isin(["all__ca_route_signature", "epi__ca_route_signature", "all__wnt_stem", "epi__wnt_stem"])].to_string(index=False))
        fh.write("\n\nBecker focused adjusted models:\n")
        focus_becker = becker_models[becker_models["term"].str.contains("disease_stage_group_polyp|proliferation", regex=True)] if not becker_models.empty else becker_models
        fh.write(focus_becker.to_string(index=False))
        fh.write("\n\nAtlas donor-level focused tests:\n")
        fh.write(atlas_tests[atlas_tests["score"].isin(["ca_route_signature", "wnt_stem", "proliferation_control"])].to_string(index=False))
        fh.write("\n\nAtlas adjusted model focused terms:\n")
        focus_atlas = atlas_models[atlas_models["term"].str.contains("carrier_group_|proliferation", regex=True)] if not atlas_models.empty else atlas_models
        fh.write(focus_atlas[focus_atlas["outcome"].isin(["ca_route_signature", "wnt_stem"])].to_string(index=False))
        fh.write("\n\nGSE39582 focused results:\n")
        if not gse_results.empty:
            focus = gse_results[gse_results["score"].isin(["ca_route_signature", "ca_route_up", "wnt_stem", "proliferation_control"])]
            fh.write(focus.to_string(index=False))
        fh.write("\n\nInterpretation boundary:\n")
        fh.write("The route signature is a data-driven premalignant epithelial-memory score. It is not a standalone recurrence-risk biomarker unless clinical cohorts support that endpoint.\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    meta_by_dataset: dict[str, pd.DataFrame] = {}
    expr_by_dataset: dict[str, pd.DataFrame] = {}
    stats_rows = []
    for dataset_name, path in CHEN_DATASETS.items():
        meta, expr = load_chen_specimen_pseudobulk(path, dataset_name)
        meta_by_dataset[dataset_name] = meta
        expr_by_dataset[dataset_name] = expr
        stats_rows.append(gene_level_tests(meta, expr, dataset_name))

    stats = pd.concat(stats_rows, ignore_index=True)
    signature = build_signature(stats)
    stats.to_csv(OUT_DIR / "chen_conventional_route_gene_level_stats.tsv", sep="\t", index=False)
    signature.to_csv(OUT_DIR / "conventional_adenoma_route_signature_genes.tsv", sep="\t", index=False)

    chen_scores, chen_tests, chen_paired = chen_signature_scores(meta_by_dataset, expr_by_dataset, signature)
    chen_scores.to_csv(OUT_DIR / "chen_conventional_route_signature_scores.tsv", sep="\t", index=False)
    chen_tests.to_csv(OUT_DIR / "chen_conventional_route_signature_tests.tsv", sep="\t", index=False)
    chen_paired.to_csv(OUT_DIR / "chen_conventional_route_signature_paired_tests.tsv", sep="\t", index=False)

    becker_scores, becker_availability = score_becker(signature)
    becker_tests, becker_models = becker_tests_and_models(becker_scores)
    becker_scores.to_csv(OUT_DIR / "becker_conventional_route_signature_scores.tsv", sep="\t", index=False)
    becker_availability.to_csv(OUT_DIR / "becker_conventional_route_signature_gene_availability.tsv", sep="\t", index=False)
    becker_tests.to_csv(OUT_DIR / "becker_conventional_route_signature_tests.tsv", sep="\t", index=False)
    becker_models.to_csv(OUT_DIR / "becker_conventional_route_signature_adjusted_models.tsv", sep="\t", index=False)

    atlas_donors, atlas_tests, atlas_models = score_atlas(signature)
    atlas_donors.to_csv(OUT_DIR / "atlas_conventional_route_signature_donor_scores.tsv", sep="\t", index=False)
    atlas_tests.to_csv(OUT_DIR / "atlas_conventional_route_signature_tests.tsv", sep="\t", index=False)
    atlas_models.to_csv(OUT_DIR / "atlas_conventional_route_signature_adjusted_models.tsv", sep="\t", index=False)

    gse_scores, gse_mapping, gse_selected, gse_results = score_gse39582(signature)
    gse_scores.to_csv(OUT_DIR / "gse39582_conventional_route_signature_scores.tsv", sep="\t", index=False)
    gse_mapping.to_csv(OUT_DIR / "gse39582_conventional_route_signature_probe_mapping_candidates.tsv", sep="\t", index=False)
    gse_selected.to_csv(OUT_DIR / "gse39582_conventional_route_signature_selected_probes.tsv", sep="\t", index=False)
    gse_results.to_csv(OUT_DIR / "gse39582_conventional_route_signature_clinical_results.tsv", sep="\t", index=False)

    write_summary(
        OUT_DIR / "conventional_route_signature_transfer_summary.txt",
        signature,
        chen_tests,
        chen_paired,
        becker_tests,
        becker_models,
        atlas_tests,
        atlas_models,
        gse_results,
    )


if __name__ == "__main__":
    main()
