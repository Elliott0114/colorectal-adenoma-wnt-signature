#!/usr/bin/env python3
"""Genome-wide Visium module analysis for Zenodo 7760264 CRC CMS ST data."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, wilcoxon


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data_sources" / "CRC_spatial_public" / "zenodo_7760264" / "extracted"
ANNOT_DIR = DATA_DIR / "Pathology_SpotAnnotations"
OUT_DIR = ROOT / "results" / "spatial_zenodo7760264"

GENE_SETS = {
    "wnt_stem": ["LGR5", "ASCL2", "OLFM4", "AXIN2", "SOX9", "EPHB2", "SMOC2"],
    "wnt_core_ihc": ["OLFM4", "SOX9", "EPHB2"],
    "serrated_metaplasia": ["MUC5AC", "MUC6", "TFF1", "TFF2", "TFF3", "REG4", "AGR2", "SPINK4", "KRT7", "ANXA10"],
    "antigen_presentation_ifn": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "B2M", "TAP1", "STAT1", "IRF1", "CXCL10", "CXCL11"],
    "proliferation_control": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "TYMS", "UBE2C", "CENPF"],
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "MUC13", "TACSTD2", "CDH1"],
    "t_cell": ["CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B", "GZMB", "NKG7"],
    "myeloid": ["LYZ", "LST1", "C1QA", "C1QB", "CD68", "CD163", "CSF1R"],
    "fibroblast_stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "PDGFRA", "FAP"],
}

SAMPLE_INFO = {
    "SN123_A551763_Rep1": {"case": "A551763", "site": "cecum", "replicate": "Rep1"},
    "SN124_A551763_Rep2": {"case": "A551763", "site": "cecum", "replicate": "Rep2"},
    "SN123_A595688_Rep1": {"case": "A595688", "site": "right_colon", "replicate": "Rep1"},
    "SN124_A595688_Rep2": {"case": "A595688", "site": "right_colon", "replicate": "Rep2"},
    "SN048_A416371_Rep1": {"case": "A416371", "site": "right_colon", "replicate": "Rep1"},
    "SN048_A416371_Rep2": {"case": "A416371", "site": "right_colon", "replicate": "Rep2"},
    "SN84_A120838_Rep1": {"case": "A120838", "site": "sigmoid_colon", "replicate": "Rep1"},
    "SN84_A120838_Rep2": {"case": "A120838", "site": "sigmoid_colon", "replicate": "Rep2"},
    "SN048_A121573_Rep1": {"case": "A121573", "site": "rectum", "replicate": "Rep1"},
    "SN048_A121573_Rep2": {"case": "A121573", "site": "rectum", "replicate": "Rep2"},
    "SN123_A938797_Rep1_X": {"case": "A938797", "site": "rectum", "replicate": "Rep1_X"},
    "SN124_A938797_Rep2": {"case": "A938797", "site": "rectum", "replicate": "Rep2"},
    "SN123_A798015_Rep1": {"case": "A798015", "site": "rectosigmoid", "replicate": "Rep1"},
    "SN124_A798015_Rep2": {"case": "A798015", "site": "rectosigmoid", "replicate": "Rep2"},
}


def decode_arr(values) -> list[str]:
    return [v.decode("utf-8", "replace") if isinstance(v, bytes) else str(v) for v in values]


def pathology_group(value: str) -> str:
    text = str(value).strip().lower()
    if not text or text == "nan":
        return "missing"
    if "exclude" in text:
        return "exclude"
    if "tumor" in text and "stroma" in text:
        return "tumor_stroma"
    if "tumor" in text:
        return "tumor"
    if "stroma" in text or "fibroblastic" in text:
        return "stroma"
    if "epithelium" in text:
        return "non_neoplastic_epithelium"
    if "submucosa" in text:
        return "submucosa"
    return "other"


def read_10x_h5(path: Path) -> tuple[list[str], list[str], csc_matrix]:
    with h5py.File(path, "r") as f:
        barcodes = decode_arr(f["matrix/barcodes"][()])
        genes = decode_arr(f["matrix/features/name"][()])
        data = f["matrix/data"][()]
        indices = f["matrix/indices"][()]
        indptr = f["matrix/indptr"][()]
        shape = tuple(f["matrix/shape"][()])
    matrix = csc_matrix((data, indices, indptr), shape=shape)
    return barcodes, genes, matrix


def selected_gene_rows(genes: list[str], sample_id: str) -> tuple[dict[str, int], pd.DataFrame]:
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    rows = []
    for module, module_genes in GENE_SETS.items():
        present = [gene for gene in module_genes if gene in gene_to_idx]
        missing = [gene for gene in module_genes if gene not in gene_to_idx]
        rows.append(
            {
                "sample_id": sample_id,
                "module": module,
                "n_requested": len(module_genes),
                "n_present": len(present),
                "present_genes": ",".join(present),
                "missing_genes": ",".join(missing),
            }
        )
    wanted = sorted({gene for module_genes in GENE_SETS.values() for gene in module_genes if gene in gene_to_idx})
    return {gene: gene_to_idx[gene] for gene in wanted}, pd.DataFrame(rows)


def log_cpm(counts: np.ndarray, library_size: np.ndarray) -> np.ndarray:
    denom = np.where(library_size > 0, library_size, np.nan)
    return np.log1p(counts / denom * 1_000_000.0)


def read_positions(path: Path) -> pd.DataFrame:
    pos = pd.read_csv(path, header=None)
    pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    return pos


def read_annotations(sample_id: str) -> pd.DataFrame:
    path = ANNOT_DIR / f"Pathologist_Annotations_{sample_id}.csv"
    ann = pd.read_csv(path)
    candidates = [c for c in ann.columns if c.lower().startswith("pathologist annotation")]
    if not candidates:
        candidates = [c for c in ann.columns if c != "Barcode"]
    annotation_col = candidates[0]
    ann = ann.rename(columns={"Barcode": "barcode", annotation_col: "pathology_annotation"})
    ann["pathology_group"] = ann["pathology_annotation"].map(pathology_group)
    return ann[["barcode", "pathology_annotation", "pathology_group"]]


def score_sample(sample_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_id = sample_dir.name
    barcodes, genes, matrix = read_10x_h5(sample_dir / "filtered_feature_bc_matrix.h5")
    selected, availability = selected_gene_rows(genes, sample_id)
    total_counts = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    selected_genes = list(selected)
    selected_idx = [selected[gene] for gene in selected_genes]
    selected_counts = matrix[selected_idx, :].toarray().astype(float) if selected_idx else np.zeros((0, len(barcodes)))
    log_values = {
        gene: log_cpm(selected_counts[pos, :], total_counts)
        for pos, gene in enumerate(selected_genes)
    }
    frame = pd.DataFrame({"sample_id": sample_id, "barcode": barcodes, "total_counts": total_counts})
    for module, module_genes in GENE_SETS.items():
        values = [log_values[gene] for gene in module_genes if gene in log_values]
        frame[f"raw__{module}"] = np.nanmean(np.vstack(values), axis=0) if values else np.nan
    positions = read_positions(sample_dir / "spatial" / "tissue_positions_list.csv")
    annotations = read_annotations(sample_id)
    frame = frame.merge(positions, on="barcode", how="left").merge(annotations, on="barcode", how="left")
    info = SAMPLE_INFO.get(sample_id, {})
    for key, value in info.items():
        frame[key] = value
    return frame, availability


def add_z_scores(spots: pd.DataFrame) -> pd.DataFrame:
    spots = spots.copy()
    raw_cols = [c for c in spots.columns if c.startswith("raw__")]
    for col in raw_cols:
        out_col = col.replace("raw__", "score__")
        spots[out_col] = np.nan
        for sample_id, idx in spots.groupby("sample_id", observed=True).groups.items():
            values = spots.loc[idx, col].astype(float)
            sd = values.std(skipna=True)
            if pd.isna(sd) or sd == 0:
                continue
            spots.loc[idx, out_col] = (values - values.mean(skipna=True)) / sd
    return spots


def module_by_pathology(spots: pd.DataFrame) -> pd.DataFrame:
    score_cols = [c for c in spots.columns if c.startswith("score__")]
    rows = []
    for keys, group in spots.groupby(["sample_id", "case", "site", "replicate", "pathology_group"], observed=True):
        sample_id, case, site, replicate, pathology = keys
        row = {
            "sample_id": sample_id,
            "case": case,
            "site": site,
            "replicate": replicate,
            "pathology_group": pathology,
            "n_spots": len(group),
        }
        for col in score_cols:
            row[col] = float(group[col].median(skipna=True))
        rows.append(row)
    return pd.DataFrame(rows)


def paired_pathology_tests(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    score_cols = [c for c in summary.columns if c.startswith("score__")]
    contrast_defs = [
        ("tumor_vs_stroma", ["tumor", "tumor_stroma"], ["stroma"]),
        ("tumor_vs_non_neoplastic_epithelium", ["tumor", "tumor_stroma"], ["non_neoplastic_epithelium"]),
        ("tumor_stroma_vs_tumor", ["tumor_stroma"], ["tumor"]),
    ]
    for label, group_a, group_b in contrast_defs:
        a = summary[summary["pathology_group"].isin(group_a)].groupby("sample_id", observed=True)[score_cols].median()
        b = summary[summary["pathology_group"].isin(group_b)].groupby("sample_id", observed=True)[score_cols].median()
        joined = a.join(b, lsuffix="__a", rsuffix="__b", how="inner")
        for col in score_cols:
            pair = joined[[f"{col}__a", f"{col}__b"]].dropna()
            if len(pair) < 4:
                continue
            diff = pair[f"{col}__a"] - pair[f"{col}__b"]
            stat = wilcoxon(diff, zero_method="wilcox", alternative="two-sided")
            rows.append(
                {
                    "unit": "sample_paired_median",
                    "comparison": label,
                    "module": col.replace("score__", ""),
                    "n_pairs": len(diff),
                    "median_delta_a_minus_b": float(diff.median()),
                    "p_wilcoxon": float(stat.pvalue),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["comparison", "p_wilcoxon", "module"])
    return out


def spot_pathology_tests(spots: pd.DataFrame) -> pd.DataFrame:
    rows = []
    score_cols = [c for c in spots.columns if c.startswith("score__")]
    contrast_defs = [
        ("tumor_vs_stroma", ["tumor", "tumor_stroma"], ["stroma"]),
        ("tumor_vs_non_neoplastic_epithelium", ["tumor", "tumor_stroma"], ["non_neoplastic_epithelium"]),
    ]
    for label, group_a, group_b in contrast_defs:
        for col in score_cols:
            aval = spots.loc[spots["pathology_group"].isin(group_a), col].dropna()
            bval = spots.loc[spots["pathology_group"].isin(group_b), col].dropna()
            if len(aval) < 20 or len(bval) < 20:
                continue
            stat = mannwhitneyu(aval, bval, alternative="two-sided")
            rows.append(
                {
                    "unit": "spot_descriptive",
                    "comparison": label,
                    "module": col.replace("score__", ""),
                    "n_a": len(aval),
                    "n_b": len(bval),
                    "median_a": float(aval.median()),
                    "median_b": float(bval.median()),
                    "delta_a_minus_b": float(aval.median() - bval.median()),
                    "p_mannwhitney": float(stat.pvalue),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["comparison", "p_mannwhitney", "module"])
    return out


def wnt_neighborhood_features(spots: pd.DataFrame, k: int = 6) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    module_cols = [
        "score__antigen_presentation_ifn",
        "score__t_cell",
        "score__myeloid",
        "score__fibroblast_stromal",
        "score__proliferation_control",
    ]
    for sample_id, group in spots.dropna(subset=["pxl_row_in_fullres", "pxl_col_in_fullres"]).groupby("sample_id", observed=True):
        tumor = group[group["pathology_group"].isin(["tumor", "tumor_stroma"])].copy()
        if len(tumor) < 20:
            continue
        coords = group[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
        tree = cKDTree(coords)
        _, idx = tree.query(tumor[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float), k=min(k + 1, len(group)))
        idx = np.atleast_2d(idx)
        if idx.shape[1] > 1:
            idx = idx[:, 1:]
        threshold = tumor["score__wnt_stem"].quantile(0.75)
        for spot_pos, (_, spot) in enumerate(tumor.iterrows()):
            row = {
                "sample_id": sample_id,
                "barcode": spot["barcode"],
                "pathology_group": spot["pathology_group"],
                "wnt_high_tumor_spot": bool(spot["score__wnt_stem"] >= threshold),
                "score__wnt_stem": spot["score__wnt_stem"],
            }
            neighbors = group.iloc[idx[spot_pos]]
            for col in module_cols:
                row[f"neighbor_mean__{col.replace('score__', '')}"] = float(neighbors[col].mean(skipna=True))
            rows.append(row)
    features = pd.DataFrame(rows)
    test_rows = []
    if not features.empty:
        value_cols = [c for c in features.columns if c.startswith("neighbor_mean__")]
        sample_medians = features.groupby(["sample_id", "wnt_high_tumor_spot"], observed=True)[value_cols].median().reset_index()
        for col in value_cols:
            wide = sample_medians.pivot(index="sample_id", columns="wnt_high_tumor_spot", values=col).dropna()
            if True not in wide.columns or False not in wide.columns or len(wide) < 4:
                continue
            diff = wide[True] - wide[False]
            stat = wilcoxon(diff, zero_method="wilcox", alternative="two-sided")
            test_rows.append(
                {
                    "unit": "sample_paired_tumor_spot_median",
                    "feature": col,
                    "n_pairs": len(diff),
                    "median_delta_wnt_high_minus_low": float(diff.median()),
                    "p_wilcoxon": float(stat.pvalue),
                }
            )
    tests = pd.DataFrame(test_rows)
    if not tests.empty:
        tests = tests.sort_values(["p_wilcoxon", "feature"])
    return features, tests


def write_summary(path: Path, spots: pd.DataFrame, paired_tests: pd.DataFrame, neigh_tests: pd.DataFrame) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"n_samples\t{spots['sample_id'].nunique()}\n")
        fh.write(f"n_cases\t{spots['case'].nunique()}\n")
        fh.write(f"n_spots\t{len(spots)}\n")
        fh.write("\nspots_by_pathology_group\n")
        fh.write(spots["pathology_group"].value_counts(dropna=False).to_string())
        focus = paired_tests[
            paired_tests["comparison"].isin(["tumor_vs_stroma", "tumor_vs_non_neoplastic_epithelium"])
            & paired_tests["module"].isin(["wnt_stem", "wnt_core_ihc", "antigen_presentation_ifn", "t_cell", "myeloid", "fibroblast_stromal"])
        ]
        fh.write("\n\nfocused_paired_pathology_tests\n")
        fh.write(focus.to_string(index=False))
        fh.write("\n\nwnt_tumor_neighborhood_tests\n")
        fh.write(neigh_tests.to_string(index=False))
        fh.write("\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sample_dirs = sorted([p for p in DATA_DIR.iterdir() if p.is_dir() and p.name.startswith("SN")])
    spot_frames = []
    availability_frames = []
    for sample_dir in sample_dirs:
        print(f"processing {sample_dir.name}", flush=True)
        spots, availability = score_sample(sample_dir)
        spot_frames.append(spots)
        availability_frames.append(availability)
    spots = pd.concat(spot_frames, ignore_index=True)
    spots = add_z_scores(spots)
    availability = pd.concat(availability_frames, ignore_index=True)
    pathology_summary = module_by_pathology(spots)
    paired_tests = paired_pathology_tests(pathology_summary)
    spot_tests = spot_pathology_tests(spots)
    neigh_features, neigh_tests = wnt_neighborhood_features(spots)

    availability.to_csv(OUT_DIR / "zenodo7760264_module_gene_availability.tsv", sep="\t", index=False)
    spots.to_csv(OUT_DIR / "zenodo7760264_spot_module_scores.tsv", sep="\t", index=False)
    pathology_summary.to_csv(OUT_DIR / "zenodo7760264_module_scores_by_pathology.tsv", sep="\t", index=False)
    paired_tests.to_csv(OUT_DIR / "zenodo7760264_sample_paired_pathology_tests.tsv", sep="\t", index=False)
    spot_tests.to_csv(OUT_DIR / "zenodo7760264_spot_pathology_tests.tsv", sep="\t", index=False)
    neigh_features.to_csv(OUT_DIR / "zenodo7760264_wnt_tumor_neighborhood_features.tsv", sep="\t", index=False)
    neigh_tests.to_csv(OUT_DIR / "zenodo7760264_wnt_tumor_neighborhood_tests.tsv", sep="\t", index=False)
    write_summary(OUT_DIR / "zenodo7760264_visium_summary.txt", spots, paired_tests, neigh_tests)


if __name__ == "__main__":
    main()
