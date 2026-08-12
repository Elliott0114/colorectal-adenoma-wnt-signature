#!/usr/bin/env python3
"""Prepare the frozen Chen-validation input for GenKI virtual knockouts.

This script is deliberately validation-only. It does not use a virtual-knockout
result, a validation phenotype, or any downstream outcome to choose cells,
genes, targets, or scoring rules.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
SOURCE_H5AD = (
    ROOT
    / "data_sources"
    / "Chen_Cell_2021_CELLxGENE"
    / "chen_validation_epithelial.h5ad"
)
PANEL_PATH = (
    ROOT
    / "results"
    / "objective_compact_panel_v2_7"
    / "objective_compact_panel_frozen.tsv"
)
CORE_PATH = (
    ROOT
    / "results"
    / "data_adaptive_panel_pilot_v2_6"
    / "stable_error_controlled_core.tsv"
)
OUT_DIR = ROOT / "results" / "virtual_knockout_validation_v2_9"
INPUT_DIR = OUT_DIR / "input"
OUTPUT_H5AD = INPUT_DIR / "chen_validation_conventional_adenoma_balanced_raw.h5ad"

ROUTE_TYPES = ("TA", "TV", "TVA")
CELLS_PER_DONOR = 128
CELL_SEED = 20260810
N_HVG = 2000
MIN_DETECTED_FRACTION = 0.01

UPSTREAM_TARGETS = ("TCF7L2", "ASCL2", "SOX4")
MODEL_SEEDS = (20260810, 20260811)


def decode(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


def read_categorical(group: h5py.Group) -> np.ndarray:
    codes = group["codes"][()]
    categories = np.array([decode(v) for v in group["categories"][()]], dtype=object)
    out = np.empty(len(codes), dtype=object)
    valid = codes >= 0
    out[valid] = categories[codes[valid]]
    out[~valid] = None
    return out


def read_h5_col(root: h5py.Group, key: str) -> np.ndarray:
    obj = root[key]
    if isinstance(obj, h5py.Group):
        return read_categorical(obj)
    values = obj[()]
    if values.dtype.kind in "SOU":
        return np.array([decode(v) for v in values], dtype=object)
    return values


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_values(values: list[str]) -> str:
    payload = "\n".join(values).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def read_selected_csr_rows(group: h5py.Group, rows: np.ndarray) -> sparse.csr_matrix:
    """Read selected CSR rows without materializing the 127-million-NNZ matrix."""
    rows = np.asarray(rows, dtype=np.int64)
    source_indptr = group["indptr"][()]
    lengths = source_indptr[rows + 1] - source_indptr[rows]
    output_indptr = np.zeros(len(rows) + 1, dtype=np.int64)
    np.cumsum(lengths, out=output_indptr[1:])
    total_nnz = int(output_indptr[-1])
    output_data = np.empty(total_nnz, dtype=np.float32)
    output_indices = np.empty(total_nnz, dtype=np.int32)

    data_ds = group["data"]
    indices_ds = group["indices"]
    for out_row, source_row in enumerate(rows):
        src_start = int(source_indptr[source_row])
        src_end = int(source_indptr[source_row + 1])
        dst_start = int(output_indptr[out_row])
        dst_end = int(output_indptr[out_row + 1])
        output_data[dst_start:dst_end] = data_ds[src_start:src_end]
        output_indices[dst_start:dst_end] = indices_ds[src_start:src_end]

    shape = tuple(int(v) for v in group.attrs["shape"])
    return sparse.csr_matrix(
        (output_data, output_indices, output_indptr),
        shape=(len(rows), shape[1]),
        dtype=np.float32,
    )


def load_frozen_sets() -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    panel = pd.read_csv(PANEL_PATH, sep="\t")
    core = pd.read_csv(CORE_PATH, sep="\t")
    panel["gene"] = panel["gene"].astype(str).str.upper()
    core["gene"] = core["gene"].astype(str).str.upper()
    panel_genes = panel.sort_values(["pair_step", "arm"])["gene"].tolist()
    core_genes = core["gene"].tolist()
    if len(panel_genes) != 12 or len(set(panel_genes)) != 12:
        raise ValueError(f"Expected 12 unique frozen panel genes, found {len(set(panel_genes))}")
    if len(core_genes) != 287 or len(set(core_genes)) != 287:
        raise ValueError(f"Expected 287 unique core genes, found {len(set(core_genes))}")
    if not set(panel_genes).issubset(core_genes):
        raise ValueError("The frozen 12-gene panel is not a subset of the frozen 287-gene core")
    return panel, core, panel_genes, core_genes


def donor_balanced_cells(obs: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    eligible = obs[obs["polyp_type"].isin(ROUTE_TYPES)].copy()
    donor_counts = eligible.groupby("donor_id", observed=True).size()
    if donor_counts.empty:
        raise ValueError("No conventional-adenoma cells found")
    if int(donor_counts.min()) < CELLS_PER_DONOR:
        raise ValueError(
            f"Smallest donor has {int(donor_counts.min())} cells; "
            f"cannot draw the frozen {CELLS_PER_DONOR} cells per donor"
        )

    selected: list[int] = []
    audit_rows: list[dict[str, object]] = []
    for donor_pos, donor in enumerate(sorted(eligible["donor_id"].unique())):
        donor_frame = eligible[eligible["donor_id"].eq(donor)]
        rng = np.random.default_rng(CELL_SEED + donor_pos)
        chosen = np.sort(
            rng.choice(
                donor_frame.index.to_numpy(dtype=np.int64),
                size=CELLS_PER_DONOR,
                replace=False,
            )
        )
        selected.extend(chosen.tolist())
        selected_frame = donor_frame.loc[chosen]
        audit_rows.append(
            {
                "donor_id": donor,
                "n_available": len(donor_frame),
                "n_selected": len(selected_frame),
                "n_specimens_available": donor_frame["specimen_id"].nunique(),
                "n_specimens_selected": selected_frame["specimen_id"].nunique(),
                "available_polyp_types": ",".join(sorted(donor_frame["polyp_type"].unique())),
                "selected_polyp_types": ",".join(sorted(selected_frame["polyp_type"].unique())),
            }
        )
    return np.array(sorted(selected), dtype=np.int64), pd.DataFrame(audit_rows)


def main() -> None:
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    panel, core, panel_genes, core_genes = load_frozen_sets()
    all_targets = list(dict.fromkeys([*UPSTREAM_TARGETS, *panel_genes]))

    with h5py.File(SOURCE_H5AD, "r", rdcc_nbytes=512 * 1024 * 1024) as handle:
        obs_all = pd.DataFrame(
            {
                "cell_id": read_h5_col(handle["obs"], "Cells"),
                "donor_id": read_h5_col(handle["obs"], "donor_id"),
                "specimen_id": read_h5_col(handle["obs"], "HTAN Specimen ID"),
                "polyp_type": read_h5_col(handle["obs"], "Polyp_Type"),
                "sample_classification": read_h5_col(
                    handle["obs"], "Sample_Classification"
                ),
                "cell_type": read_h5_col(handle["obs"], "Cell_Type"),
            }
        )
        selected_rows, donor_audit = donor_balanced_cells(obs_all)
        selected_obs = obs_all.loc[selected_rows].copy()
        selected_obs["source_row_index"] = selected_rows
        selected_obs.index = selected_obs["cell_id"].astype(str)

        genes = np.array(
            [str(v).upper() for v in read_h5_col(handle["raw"]["var"], "feature_name")],
            dtype=object,
        )
        counts = read_selected_csr_rows(handle["raw"]["X"], selected_rows)

    first_occurrence = ~pd.Index(genes).duplicated(keep="first")
    duplicate_rows = pd.DataFrame(
        {
            "gene": genes[~first_occurrence],
            "action": "dropped_duplicate_symbol_after_first_occurrence",
        }
    )
    genes = genes[first_occurrence]
    counts = counts[:, first_occurrence].tocsr()

    n_cells = counts.shape[0]
    detected_cells = np.asarray(counts.getnnz(axis=0)).ravel()
    detection_fraction = detected_cells / n_cells
    mean_count = np.asarray(counts.mean(axis=0)).ravel()
    min_detected_cells = max(3, math.ceil(MIN_DETECTED_FRACTION * n_cells))
    hvg_eligible = detected_cells >= min_detected_cells

    hvg_adata = ad.AnnData(
        X=counts[:, hvg_eligible].copy(),
        var=pd.DataFrame(index=pd.Index(genes[hvg_eligible], name="gene")),
    )
    sc.pp.normalize_total(hvg_adata, target_sum=1e4)
    sc.pp.log1p(hvg_adata)
    sc.pp.highly_variable_genes(
        hvg_adata,
        flavor="seurat",
        n_top_genes=min(N_HVG, hvg_adata.n_vars),
        inplace=True,
    )
    hvg_genes = set(hvg_adata.var_names[hvg_adata.var["highly_variable"]])

    gene_to_position = {gene: i for i, gene in enumerate(genes)}
    forced_genes = set(core_genes) | set(panel_genes) | set(all_targets)
    forced_detected = {
        gene
        for gene in forced_genes
        if gene in gene_to_position
        and detected_cells[gene_to_position[gene]] >= min_detected_cells
    }
    selected_genes = sorted(hvg_genes | forced_detected)
    selected_positions = np.array(
        [gene_to_position[gene] for gene in selected_genes], dtype=np.int64
    )

    selected_var = pd.DataFrame(index=pd.Index(selected_genes, name="gene"))
    selected_var["source_feature_index"] = selected_positions
    selected_var["detected_cells"] = detected_cells[selected_positions]
    selected_var["detection_fraction"] = detection_fraction[selected_positions]
    selected_var["mean_raw_count"] = mean_count[selected_positions]
    selected_var["hvg_background"] = selected_var.index.isin(hvg_genes)
    selected_var["fixed_287_core"] = selected_var.index.isin(core_genes)
    selected_var["fixed_12_panel"] = selected_var.index.isin(panel_genes)
    selected_var["prespecified_ko_target"] = selected_var.index.isin(all_targets)
    selected_var["panel_arm"] = selected_var.index.map(
        panel.set_index("gene")["arm"].to_dict()
    )
    selected_var["panel_route_weight"] = selected_var.index.map(
        panel.set_index("gene")["route_weight"].to_dict()
    )

    missing_targets = sorted(set(all_targets) - set(selected_genes))
    missing_panel = sorted(set(panel_genes) - set(selected_genes))
    if missing_targets:
        raise ValueError(f"Prespecified KO targets not measurable: {missing_targets}")
    if missing_panel:
        raise ValueError(f"Frozen panel genes not measurable: {missing_panel}")

    output = ad.AnnData(
        X=counts[:, selected_positions].copy(),
        obs=selected_obs,
        var=selected_var,
    )
    output.uns["validation_only"] = True
    output.uns["panel_genes"] = panel_genes
    output.uns["core_genes_measurable"] = sorted(set(core_genes) & set(selected_genes))
    output.uns["upstream_targets"] = list(UPSTREAM_TARGETS)
    output.uns["all_prespecified_targets"] = all_targets
    output.uns["cell_sampling_seed"] = CELL_SEED
    output.uns["cells_per_donor"] = CELLS_PER_DONOR
    output.uns["hvg_count_requested"] = N_HVG
    output.write_h5ad(OUTPUT_H5AD, compression="gzip")

    fixed_audit_genes = sorted(forced_genes)
    gene_audit_rows = []
    for gene in fixed_audit_genes:
        position = gene_to_position.get(gene)
        gene_audit_rows.append(
            {
                "gene": gene,
                "in_source": position is not None,
                "detected_cells": int(detected_cells[position]) if position is not None else 0,
                "detection_fraction": (
                    float(detection_fraction[position]) if position is not None else 0.0
                ),
                "passes_detection_filter": (
                    bool(detected_cells[position] >= min_detected_cells)
                    if position is not None
                    else False
                ),
                "included_in_genki_universe": gene in selected_genes,
                "fixed_287_core": gene in core_genes,
                "fixed_12_panel": gene in panel_genes,
                "prespecified_ko_target": gene in all_targets,
            }
        )
    gene_audit = pd.DataFrame(gene_audit_rows)

    donor_audit.to_csv(INPUT_DIR / "input_cell_audit.tsv", sep="\t", index=False)
    gene_audit.to_csv(INPUT_DIR / "input_fixed_gene_audit.tsv", sep="\t", index=False)
    selected_var.reset_index().to_csv(
        INPUT_DIR / "input_gene_universe.tsv", sep="\t", index=False
    )
    duplicate_rows.to_csv(
        INPUT_DIR / "input_duplicate_symbol_audit.tsv", sep="\t", index=False
    )

    contract = {
        "analysis_role": "validation_not_discovery",
        "frozen_before_virtual_knockout": True,
        "source_dataset": "Chen et al. independent validation epithelial cohort",
        "eligible_cells": "Polyp_Type in TA, TV, or TVA",
        "sampling": {
            "unit": "donor",
            "cells_per_donor": CELLS_PER_DONOR,
            "seed": CELL_SEED,
            "phenotype_or_signature_score_used": False,
        },
        "gene_universe": {
            "background": f"top {N_HVG} label-blind HVGs in balanced adenoma cells",
            "forced": "measurable frozen 287-gene core, 12-gene panel, and prespecified targets",
            "minimum_detected_fraction": MIN_DETECTED_FRACTION,
            "validation_outcome_used": False,
        },
        "frozen_sets": {
            "core_requested_n": len(core_genes),
            "core_measurable_n": int(gene_audit.query("fixed_287_core and included_in_genki_universe").shape[0]),
            "core_sha256": sha256_values(core_genes),
            "panel_n": len(panel_genes),
            "panel_genes": panel_genes,
            "panel_sha256": sha256_values(panel_genes),
        },
        "prespecified_targets": {
            "upstream_context": list(UPSTREAM_TARGETS),
            "panel_members": panel_genes,
            "all_unique": all_targets,
        },
        "model": {
            "method": "GenKI",
            "version": "0.2.1",
            "model_seeds": list(MODEL_SEEDS),
            "epochs": 100,
            "pc_components": 5,
            "edge_percentile_cutoff": 85,
            "latent_dimensions": 2,
            "primary_distance": "KL divergence",
            "distance_sensitivity": "Earth mover distance",
            "seed_handling": (
                "fixed GenKI link split; requested initialization seed restored "
                "after the package split routine"
            ),
        },
        "prespecified_interpretation": {
            "virtual_knockout_output": "unsigned network-embedding perturbation magnitude",
            "direction_inferred_from_genki": False,
            "direction_calibration": "existing empirical APC/WNT/TCF7L2/ASCL2 perturbations only",
            "new_gene_nomination_allowed": False,
            "panel_reselection_or_reweighting_allowed": False,
        },
        "prespecified_tests": [
            "fixed 287-core enrichment after each upstream-context knockout",
            "leave-target-out fixed 12-panel enrichment after each upstream-context knockout",
            "leave-target-out within-panel network coherence across the 12 panel-member knockouts",
            "cross-seed rank stability",
            "expression-and-network-degree-matched gene-set nulls",
        ],
        "matched_null_replicates": 10000,
        "multiple_testing": "Benjamini-Hochberg within the prespecified target-by-gene-set family",
    }
    (INPUT_DIR / "frozen_validation_contract.json").write_text(
        json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    manifest = {
        "source_h5ad": str(SOURCE_H5AD.relative_to(ROOT)),
        "source_h5ad_sha256": sha256_file(SOURCE_H5AD),
        "panel_path": str(PANEL_PATH.relative_to(ROOT)),
        "panel_sha256": sha256_file(PANEL_PATH),
        "core_path": str(CORE_PATH.relative_to(ROOT)),
        "core_sha256": sha256_file(CORE_PATH),
        "output_h5ad": str(OUTPUT_H5AD.relative_to(ROOT)),
        "output_h5ad_sha256": sha256_file(OUTPUT_H5AD),
        "n_source_cells": int(len(obs_all)),
        "n_eligible_adenoma_cells": int(obs_all["polyp_type"].isin(ROUTE_TYPES).sum()),
        "n_selected_cells": int(output.n_obs),
        "n_selected_donors": int(output.obs["donor_id"].nunique()),
        "n_selected_genes": int(output.n_vars),
        "n_hvg_background": int(selected_var["hvg_background"].sum()),
        "n_core_measurable": int(selected_var["fixed_287_core"].sum()),
        "n_panel_measurable": int(selected_var["fixed_12_panel"].sum()),
        "min_detected_cells": int(min_detected_cells),
        "duplicate_symbols_dropped": duplicate_rows["gene"].tolist(),
    }
    (INPUT_DIR / "input_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print(json.dumps(manifest, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
