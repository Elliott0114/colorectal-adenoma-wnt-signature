#!/usr/bin/env python3
"""Module scoring and WNT-high epithelial neighborhood analysis for Zenodo 14602110."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data_sources" / "CRC_spatial_public" / "zenodo_14602110"
OUT_DIR = ROOT / "results" / "spatial"

GENE_SETS = {
    "wnt_stem_core": ["LGR5", "ASCL2", "OLFM4"],
    "wnt_extended": ["LGR5", "ASCL2", "OLFM4", "AXIN2", "SOX9", "EPHB2", "SMOC2", "CTNNB1"],
    "antigen_presentation_ifn": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "B2M", "TAP1", "STAT1", "IRF1", "CXCL10", "CXCL11"],
    "t_cell": ["CD3D", "CD3E", "CD4", "CD8A", "CD8B", "GZMB", "NKG7", "FOXP3"],
    "myeloid_macrophage": ["CD14", "C1QA", "C1QB", "C1QC", "CD68", "CD163", "SPP1", "LYZ"],
    "fibroblast_stromal": ["COL1A1", "COL1A2", "ACTA2", "PDGFRB", "DCN", "LUM", "TAGLN"],
    "proliferation_control": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "TYMS", "UBE2C", "CENPF"],
}

LINEAGES = ["Epithelial", "Lymphoid", "Myeloid", "Fibroblast", "Endothelial"]


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, object]:
    meta = pd.read_csv(DATA_DIR / "CRC_scRNAseq-spatial_meta_data.csv")
    genes = pd.read_csv(DATA_DIR / "CRC_scRNAseq-spatial_genes.csv")["x"].astype(str).tolist()
    barcodes = pd.read_csv(DATA_DIR / "CRC_scRNAseq-spatial_barcodes.csv")["x"].astype(str).tolist()
    counts = mmread(DATA_DIR / "CRC_scRNAseq-spatial_counts.mtx").tocsr()
    if counts.shape != (len(genes), len(barcodes)):
        raise ValueError(f"Matrix shape {counts.shape} does not match genes/barcodes {len(genes), len(barcodes)}.")
    if not meta["CellID"].astype(str).reset_index(drop=True).equals(pd.Series(barcodes, name="CellID")):
        raise ValueError("Metadata CellID order does not match barcode order.")
    return meta, pd.DataFrame({"gene": genes}), pd.DataFrame({"CellID": barcodes}), counts


def selected_gene_matrix(counts, genes: list[str]) -> tuple[np.ndarray, list[str], pd.DataFrame]:
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    selected = []
    seen = set()
    rows = []
    for module, module_genes in GENE_SETS.items():
        present = [gene for gene in module_genes if gene in gene_to_idx]
        missing = [gene for gene in module_genes if gene not in gene_to_idx]
        rows.append(
            {
                "module": module,
                "n_present": len(present),
                "present_genes": ",".join(present),
                "missing_genes": ",".join(missing),
            }
        )
        for gene in present:
            if gene not in seen:
                selected.append(gene)
                seen.add(gene)
    selected_idx = [gene_to_idx[gene] for gene in selected]
    dense = counts[selected_idx, :].T.toarray().astype(np.float32)
    totals = np.asarray(counts.sum(axis=0)).ravel().astype(np.float32)
    totals[totals <= 0] = 1
    dense = np.log1p(dense / totals[:, None] * 1e4)
    return dense, selected, pd.DataFrame(rows)


def score_modules(expr: np.ndarray, selected_genes: list[str], spatial_mask: np.ndarray) -> pd.DataFrame:
    spatial_expr = expr[spatial_mask]
    mean = spatial_expr.mean(axis=0, keepdims=True)
    sd = spatial_expr.std(axis=0, keepdims=True)
    sd[sd == 0] = 1
    z = (expr - mean) / sd
    gene_to_pos = {gene: idx for idx, gene in enumerate(selected_genes)}
    scores = {}
    for module, module_genes in GENE_SETS.items():
        present = [gene for gene in module_genes if gene in gene_to_pos]
        positions = [gene_to_pos[gene] for gene in present]
        scores[f"score__{module}"] = z[:, positions].mean(axis=1) if positions else np.nan
    return pd.DataFrame(scores)


def summarize_scores(cells: pd.DataFrame, score_cols: list[str]) -> pd.DataFrame:
    grouped = (
        cells.groupby(["batch", "sample_id", "lineage", "cresc_publication_type"], observed=True)[score_cols]
        .agg(["count", "median", "mean"])
        .reset_index()
    )
    grouped.columns = ["__".join([str(x) for x in col if x]) for col in grouped.columns.to_flat_index()]
    return grouped


def neighborhood_features_for_batch(batch_cells: pd.DataFrame, k: int = 10) -> pd.DataFrame:
    coords = batch_cells[["center_x", "center_y"]].to_numpy(float)
    if len(batch_cells) <= k + 1:
        return pd.DataFrame()
    tree = cKDTree(coords)
    _, neigh_idx = tree.query(coords, k=k + 1)
    neigh_idx = neigh_idx[:, 1:]

    epithelial = batch_cells["lineage"].eq("Epithelial").to_numpy()
    wnt = batch_cells["score__wnt_stem_core"].to_numpy(float)
    if epithelial.sum() < 50:
        return pd.DataFrame()
    epithelial_wnt = wnt[epithelial]
    background = np.nanmin(epithelial_wnt)
    high_mask = epithelial & (wnt > background + 1e-8)
    low_mask = epithelial & (wnt <= background + 1e-8)
    if high_mask.sum() < 20:
        q95 = np.nanquantile(epithelial_wnt, 0.95)
        high_mask = epithelial & (wnt >= q95)
        low_mask = epithelial & (wnt < q95)
    source_mask = high_mask | low_mask
    source_rows = np.flatnonzero(source_mask)

    score_cols = [c for c in batch_cells.columns if c.startswith("score__")]
    source_neighbors = neigh_idx[source_rows]
    lineage_array = batch_cells["lineage"].astype(str).to_numpy()
    neighbor_lineages = lineage_array[source_neighbors]
    score_matrix = batch_cells[score_cols].to_numpy(float)
    neighbor_score_means = score_matrix[source_neighbors].mean(axis=1)

    out = pd.DataFrame(
        {
            "CellID": batch_cells["CellID"].to_numpy()[source_rows],
            "batch": batch_cells["batch"].to_numpy()[source_rows],
            "sample_id": batch_cells["sample_id"].to_numpy()[source_rows],
            "wnt_group": np.where(high_mask[source_rows], "wnt_high", "wnt_low"),
            "source_wnt_score": wnt[source_rows].astype(float),
            "source_ap_ifn_score": batch_cells["score__antigen_presentation_ifn"].to_numpy(float)[source_rows],
            "source_proliferation_score": batch_cells["score__proliferation_control"].to_numpy(float)[source_rows],
        }
    )
    for lineage in LINEAGES:
        out[f"neighbor_frac__{lineage}"] = (neighbor_lineages == lineage).mean(axis=1)
    for idx, col in enumerate(score_cols):
        out[f"neighbor_mean__{col.replace('score__', '')}"] = neighbor_score_means[:, idx]
    return out


def compare_wnt_neighborhoods(features: pd.DataFrame) -> pd.DataFrame:
    value_cols = [c for c in features.columns if c.startswith(("neighbor_frac__", "neighbor_mean__", "source_"))]
    rows = []
    for subset_name, frame in {"all_spatial": features, **{b: g for b, g in features.groupby("batch", observed=True)}}.items():
        for col in value_cols:
            high = frame.loc[frame["wnt_group"].eq("wnt_high"), col].dropna()
            low = frame.loc[frame["wnt_group"].eq("wnt_low"), col].dropna()
            if len(high) < 20 or len(low) < 20:
                continue
            stat = mannwhitneyu(high, low, alternative="two-sided")
            rows.append(
                {
                    "subset": subset_name,
                    "feature": col,
                    "n_wnt_high": len(high),
                    "n_wnt_low": len(low),
                    "median_wnt_high": float(high.median()),
                    "median_wnt_low": float(low.median()),
                    "delta_high_minus_low": float(high.median() - low.median()),
                    "p_mannwhitney_cell_level": float(stat.pvalue),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["subset", "p_mannwhitney_cell_level", "feature"])
    return out


def write_summary(path: Path, spatial: pd.DataFrame, features: pd.DataFrame, tests: pd.DataFrame, availability: pd.DataFrame) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"n_spatial_cells\t{len(spatial)}\n")
        fh.write(f"n_spatial_batches\t{spatial['batch'].nunique()}\n")
        fh.write(f"n_wnt_neighborhood_sources\t{len(features)}\n")
        fh.write("\nspatial_cells_by_batch\n")
        fh.write(spatial["batch"].value_counts().to_string())
        fh.write("\n\nspatial_cells_by_lineage\n")
        fh.write(spatial["lineage"].value_counts().to_string())
        fh.write("\n\nmodule_gene_availability\n")
        fh.write(availability.to_string(index=False))
        if not tests.empty:
            focus = tests[tests["feature"].isin([
                "source_ap_ifn_score",
                "source_proliferation_score",
                "neighbor_mean__antigen_presentation_ifn",
                "neighbor_mean__t_cell",
                "neighbor_mean__myeloid_macrophage",
                "neighbor_mean__fibroblast_stromal",
                "neighbor_frac__Lymphoid",
                "neighbor_frac__Myeloid",
                "neighbor_frac__Fibroblast",
            ])]
            fh.write("\n\nfocused_wnt_high_vs_low_neighborhood_tests\n")
            fh.write(focus.to_string(index=False))
        fh.write("\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta, genes_df, _, counts = load_inputs()
    spatial_mask = meta["batch"].ne("scRNAseq").to_numpy()
    expr, selected_genes, availability = selected_gene_matrix(counts, genes_df["gene"].tolist())
    scores = score_modules(expr, selected_genes, spatial_mask)
    cells = pd.concat([meta, scores], axis=1)
    spatial = cells.loc[spatial_mask].copy()
    score_cols = [c for c in spatial.columns if c.startswith("score__")]
    score_summary = summarize_scores(spatial, score_cols)

    neighborhood_frames = []
    for _, batch_cells in spatial.groupby("batch", sort=False, observed=True):
        neighborhood_frames.append(neighborhood_features_for_batch(batch_cells.reset_index(drop=True)))
    neighborhoods = pd.concat([x for x in neighborhood_frames if not x.empty], ignore_index=True)
    tests = compare_wnt_neighborhoods(neighborhoods)

    availability.to_csv(OUT_DIR / "zenodo14602110_module_gene_availability.tsv", sep="\t", index=False)
    spatial.to_csv(OUT_DIR / "zenodo14602110_spatial_cell_module_scores.tsv", sep="\t", index=False)
    score_summary.to_csv(OUT_DIR / "zenodo14602110_module_scores_by_celltype.tsv", sep="\t", index=False)
    neighborhoods.to_csv(OUT_DIR / "zenodo14602110_wnt_neighborhood_features.tsv", sep="\t", index=False)
    tests.to_csv(OUT_DIR / "zenodo14602110_wnt_neighborhood_tests.tsv", sep="\t", index=False)
    write_summary(OUT_DIR / "zenodo14602110_spatial_module_summary.txt", spatial, neighborhoods, tests, availability)


if __name__ == "__main__":
    main()
