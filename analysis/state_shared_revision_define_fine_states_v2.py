#!/usr/bin/env python3
"""Define normal-reference epithelial fine states without programme genes.

Analysis date: 2026-08-30
Random seed: 20260830
Analysis unit for downstream inference: donor; cells are measurement units only.
"""

from __future__ import annotations

import gzip
import gc
import hashlib
import json
import platform
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy import sparse
from scipy.stats import rankdata
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score


ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "results" / "state_aware_program_v1"
OUT_DIR = ROOT / "results" / "state_shared_revision_v2" / "fine_states"
CONTRACT = (
    ROOT
    / "analysis"
    / "contracts"
    / "state_shared_revision_validation_v2_2026-08-30.md"
)
H5AD = {
    "discovery": (
        ROOT
        / "data_sources"
        / "Chen_Cell_2021_CELLxGENE"
        / "chen_discovery_epithelial.h5ad"
    ),
    "validation": (
        ROOT
        / "data_sources"
        / "Chen_Cell_2021_CELLxGENE"
        / "chen_validation_epithelial.h5ad"
    ),
}
COMMON_EFFECTS = PARENT / "common_effects" / "cross_state_common_effects.tsv.gz"

STATES = ("ABS", "GOB", "TAC")
K_VALUES = (3, 4, 5)
PRIMARY_K = 4
N_HVG = 1500
N_PCS = 20
SEED = 20260830
ROW_BATCH = 512
ROUTE_MAP = {
    "NL": "normal",
    "TA": "conventional_adenoma",
    "TV": "conventional_adenoma",
    "TVA": "conventional_adenoma",
}

EXCLUDED_SYMBOLS = {
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
    "CCNB1",
    "CCNB2",
    "CDK1",
    "CDC20",
    "CDC25A",
    "CDC25B",
    "CDC25C",
    "FOS",
    "FOSB",
    "JUN",
    "JUNB",
    "JUND",
    "ATF3",
    "EGR1",
    "EGR2",
    "EGR3",
    "DUSP1",
    "DUSP2",
    "HSPA1A",
    "HSPA1B",
    "HSP90AA1",
    "HSP90AB1",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


def read_categorical(group: h5py.Group) -> np.ndarray:
    codes = group["codes"][()]
    categories = np.array([decode(value) for value in group["categories"][()]])
    output = np.full(len(codes), None, dtype=object)
    valid = codes >= 0
    output[valid] = categories[codes[valid]]
    return output


def read_column(root: h5py.Group, name: str) -> np.ndarray:
    node = root[name]
    if isinstance(node, h5py.Group):
        return read_categorical(node)
    values = node[()]
    if values.dtype.kind in "SOU":
        return np.array([decode(value) for value in values], dtype=object)
    return values


def load_metadata(path: Path, partition: str) -> tuple[pd.DataFrame, np.ndarray]:
    with h5py.File(path, "r") as handle:
        genes = read_column(handle["var"], "feature_name").astype(str)
        metadata = pd.DataFrame(
            {
                "partition": partition,
                "specimen_id": read_column(handle["obs"], "HTAN Specimen ID"),
                "donor_id": read_column(handle["obs"], "donor_id"),
                "polyp_type": read_column(handle["obs"], "Polyp_Type"),
                "broad_state": read_column(handle["obs"], "Cell_Type"),
                "tissue": read_column(handle["obs"], "tissue"),
            }
        )
    metadata["route"] = metadata["polyp_type"].map(ROUTE_MAP)
    metadata["cell_index"] = np.arange(len(metadata), dtype=np.int64)
    return metadata, genes


def excluded_gene(gene: str, programme: set[str]) -> bool:
    upper = gene.upper()
    return (
        upper in programme
        or upper in EXCLUDED_SYMBOLS
        or upper.startswith("MT-")
        or upper.startswith("RPL")
        or upper.startswith("RPS")
        or upper.startswith("IGH")
        or upper.startswith("IGK")
        or upper.startswith("IGL")
    )


def feature_statistics_all_states(
    path: Path,
    normal_indices_by_state: dict[str, np.ndarray],
    genes: np.ndarray,
    eligible_gene_names: set[str],
) -> dict[str, pd.DataFrame]:
    eligible_mask = np.array(
        [gene in eligible_gene_names for gene in genes], dtype=bool
    )
    eligible_columns = np.flatnonzero(eligible_mask).astype(np.int64)
    eligible_genes = genes[eligible_mask]
    with h5py.File(path, "r") as handle:
        raw_group = handle["raw"]["X"]
        raw_genes = read_column(handle["raw"]["var"], "feature_name").astype(str)
        if not np.array_equal(raw_genes, genes):
            raise ValueError("raw and transformed gene orders differ")
        raw = sparse.csr_matrix(
            (
                raw_group["data"][()].astype(np.float32, copy=False),
                raw_group["indices"][()].astype(np.int32, copy=False),
                raw_group["indptr"][()].astype(np.int64, copy=False),
            ),
            shape=tuple(int(value) for value in raw_group.attrs["shape"]),
        )
    row_library_sizes = np.asarray(raw.sum(axis=1)).ravel()
    if np.any(row_library_sizes <= 0):
        raise ValueError("A raw-count row has zero library size")
    outputs: dict[str, pd.DataFrame] = {}
    for state in STATES:
        rows = normal_indices_by_state[state]
        values = raw[rows, :][:, eligible_columns].astype(np.float32)
        values = sparse.diags(
            (10_000.0 / row_library_sizes[rows]).astype(np.float32)
        ).dot(values).tocsr()
        values.data = np.log1p(values.data)
        means = np.asarray(values.mean(axis=0)).ravel().astype(np.float64)
        squared = values.copy()
        squared.data = np.square(squared.data)
        sum_of_squares = np.asarray(squared.sum(axis=0)).ravel()
        variances = (
            sum_of_squares - len(rows) * np.square(means)
        ) / (len(rows) - 1)
        detection = np.asarray(values.getnnz(axis=0)).ravel() / len(rows)
        output = pd.DataFrame(
            {
                "gene": eligible_genes,
                "mean": means,
                "variance": variances,
                "detected_fraction": detection,
            }
        )
        output = output.loc[
            np.isfinite(output["variance"])
            & output["variance"].gt(0)
            & output["detected_fraction"].between(0.01, 0.99)
        ].copy()
        output = output.sort_values(
            ["variance", "detected_fraction", "gene"],
            ascending=[False, False, True],
        ).reset_index(drop=True)
        output["variance_rank"] = np.arange(1, len(output) + 1)
        output["selected_hvg"] = output["variance_rank"].le(N_HVG)
        outputs[state] = output
        del values, squared
        gc.collect()
    del raw
    gc.collect()
    return outputs


def read_matrix(
    path: Path,
    row_indices: np.ndarray,
    all_genes: np.ndarray,
    selected_genes: list[str],
) -> np.ndarray:
    gene_to_index = {gene: index for index, gene in enumerate(all_genes)}
    pairs = sorted((gene_to_index[gene], gene) for gene in selected_genes)
    columns = np.array([index for index, _ in pairs], dtype=np.int64)
    order = np.array([selected_genes.index(gene) for _, gene in pairs], dtype=np.int64)
    inverse = np.argsort(order)
    sorted_rows = np.sort(np.asarray(row_indices, dtype=np.int64))
    chunks: list[np.ndarray] = []
    with h5py.File(path, "r") as handle:
        matrix = handle["X"]
        for start in range(0, matrix.shape[0], ROW_BATCH):
            stop = min(start + ROW_BATCH, matrix.shape[0])
            local_rows = sorted_rows[(sorted_rows >= start) & (sorted_rows < stop)]
            if not len(local_rows):
                continue
            values = matrix[start:stop, columns].astype(np.float32)
            chunks.append(values[local_rows - start, :][:, inverse])
    return np.vstack(chunks)


def standardize(values: np.ndarray, mean: np.ndarray, scale: np.ndarray) -> np.ndarray:
    output = (values.astype(np.float64) - mean) / scale
    return np.clip(output, -10.0, 10.0).astype(np.float32)


def programme_rank_score(
    path: Path,
    row_indices: np.ndarray,
    all_genes: np.ndarray,
    up_genes: list[str],
    down_genes: list[str],
) -> np.ndarray:
    requested = up_genes + down_genes
    gene_to_index = {gene: index for index, gene in enumerate(all_genes)}
    present_up = [gene for gene in up_genes if gene in gene_to_index]
    present_down = [gene for gene in down_genes if gene in gene_to_index]
    requested_present = present_up + present_down
    if len(present_up) < 0.9 * len(up_genes) or len(present_down) < 0.9 * len(down_genes):
        raise ValueError("Fewer than 90% of one programme arm is measurable")
    pairs = sorted((gene_to_index[gene], gene) for gene in requested_present)
    columns = np.array([index for index, _ in pairs], dtype=np.int64)
    ordered_genes = [gene for _, gene in pairs]
    up_positions = np.array(
        [index for index, gene in enumerate(ordered_genes) if gene in set(present_up)]
    )
    down_positions = np.array(
        [index for index, gene in enumerate(ordered_genes) if gene in set(present_down)]
    )
    sorted_rows = np.sort(np.asarray(row_indices, dtype=np.int64))
    output: list[np.ndarray] = []
    with h5py.File(path, "r") as handle:
        matrix = handle["X"]
        for start in range(0, matrix.shape[0], ROW_BATCH):
            stop = min(start + ROW_BATCH, matrix.shape[0])
            local_rows = sorted_rows[(sorted_rows >= start) & (sorted_rows < stop)]
            if not len(local_rows):
                continue
            values = matrix[start:stop, columns].astype(np.float64)
            values = values[local_rows - start, :]
            ranks = rankdata(values, axis=1, method="average") / values.shape[1]
            score = ranks[:, up_positions].mean(axis=1) - ranks[:, down_positions].mean(axis=1)
            output.append(score.astype(np.float32))
    return np.concatenate(output)


def write_tsv_gzip(frame: pd.DataFrame, path: Path) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, sep="\t", index=False)


def main() -> None:
    np.random.seed(SEED)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in [*H5AD.values(), COMMON_EFFECTS, CONTRACT]:
        if not path.exists():
            raise FileNotFoundError(path)

    common = pd.read_csv(COMMON_EFFECTS, sep="\t")
    strict = common.loc[common["strict_state_shared"].astype(bool)].copy()
    up_genes = sorted(strict.loc[strict["shared_direction"].eq("up"), "gene"])
    down_genes = sorted(strict.loc[strict["shared_direction"].eq("down"), "gene"])
    programme = set(up_genes + down_genes)
    if len(programme) != 1843:
        raise ValueError("The frozen state-shared programme is not 1,843 genes")

    metadata: dict[str, pd.DataFrame] = {}
    genes: dict[str, np.ndarray] = {}
    for partition, path in H5AD.items():
        metadata[partition], genes[partition] = load_metadata(path, partition)
    common_gene_names = set(genes["discovery"]).intersection(genes["validation"])
    eligible = {
        gene
        for gene in common_gene_names
        if not excluded_gene(gene, programme)
    }

    overlap = set(metadata["discovery"]["donor_id"]).intersection(
        metadata["validation"]["donor_id"]
    )
    validation_allowed = ~metadata["validation"]["donor_id"].isin(overlap)

    normal_indices_by_state = {
        state: metadata["discovery"].loc[
            metadata["discovery"]["broad_state"].eq(state)
            & metadata["discovery"]["route"].eq("normal"),
            "cell_index",
        ].to_numpy()
        for state in STATES
    }
    feature_statistics_by_state = feature_statistics_all_states(
        H5AD["discovery"],
        normal_indices_by_state,
        genes["discovery"],
        eligible,
    )

    assignment_records: list[pd.DataFrame] = []
    feature_records: list[pd.DataFrame] = []
    model_records: list[dict[str, object]] = []
    centroid_records: list[pd.DataFrame] = []

    for state_index, state in enumerate(STATES):
        normal_indices = normal_indices_by_state[state]
        statistics = feature_statistics_by_state[state].copy()
        statistics.insert(0, "broad_state", state)
        selected = statistics.loc[statistics["selected_hvg"], "gene"].tolist()
        if len(selected) != N_HVG:
            raise ValueError(f"Only {len(selected)} eligible HVGs were selected for {state}")
        feature_records.append(statistics)

        normal_values = read_matrix(
            H5AD["discovery"], normal_indices, genes["discovery"], selected
        )
        centre = normal_values.mean(axis=0, dtype=np.float64)
        scale = normal_values.std(axis=0, ddof=1, dtype=np.float64)
        if np.any(~np.isfinite(scale)) or np.any(scale <= 0):
            raise ValueError(f"Invalid normal-reference feature scale for {state}")
        normal_z = standardize(normal_values, centre, scale)
        pca = PCA(
            n_components=N_PCS,
            svd_solver="randomized",
            random_state=SEED + state_index,
        )
        normal_pcs = pca.fit_transform(normal_z)

        fitted_clusters: dict[int, KMeans] = {}
        for k in K_VALUES:
            model = KMeans(
                n_clusters=k,
                n_init=50,
                random_state=SEED + 10 * state_index + k,
            ).fit(normal_pcs)
            fitted_clusters[k] = model
            sample_size = min(5000, len(normal_pcs))
            sample_rng = np.random.default_rng(SEED + 100 * state_index + k)
            sample_index = sample_rng.choice(
                len(normal_pcs), size=sample_size, replace=False
            )
            silhouette = silhouette_score(
                normal_pcs[sample_index], model.labels_[sample_index]
            )
            model_records.append(
                {
                    "broad_state": state,
                    "k": k,
                    "n_normal_reference_cells": len(normal_pcs),
                    "n_hvg": len(selected),
                    "n_pcs": N_PCS,
                    "explained_variance_fraction": pca.explained_variance_ratio_.sum(),
                    "inertia": model.inertia_,
                    "silhouette_sample_n": sample_size,
                    "silhouette": silhouette,
                }
            )
            centroids = pd.DataFrame(
                model.cluster_centers_,
                columns=[f"PC{index + 1}" for index in range(N_PCS)],
            )
            centroids.insert(0, "fine_state", [f"{state}_F{i + 1}" for i in range(k)])
            centroids.insert(0, "k", k)
            centroids.insert(0, "broad_state", state)
            centroid_records.append(centroids)

        for partition in ("discovery", "validation"):
            target_mask = (
                metadata[partition]["broad_state"].eq(state)
                & metadata[partition]["route"].isin(
                    ["normal", "conventional_adenoma"]
                )
            )
            if partition == "validation":
                target_mask &= validation_allowed
            target = metadata[partition].loc[target_mask].copy()
            target_indices = target["cell_index"].to_numpy()
            target_values = read_matrix(
                H5AD[partition], target_indices, genes[partition], selected
            )
            target_pcs = pca.transform(standardize(target_values, centre, scale))
            for k, model in fitted_clusters.items():
                labels = model.predict(target_pcs)
                target[f"fine_state_k{k}"] = [
                    f"{state}_F{label + 1}" for label in labels
                ]
                distances = model.transform(target_pcs)
                target[f"centroid_distance_k{k}"] = distances[
                    np.arange(len(labels)), labels
                ]
            target["programme_rank_score"] = programme_rank_score(
                H5AD[partition],
                target_indices,
                genes[partition],
                up_genes,
                down_genes,
            )
            assignment_records.append(target)

    assignments = pd.concat(assignment_records, ignore_index=True)
    features = pd.concat(feature_records, ignore_index=True)
    model_diagnostics = pd.DataFrame(model_records)
    centroids = pd.concat(centroid_records, ignore_index=True)

    for k in K_VALUES:
        counts = assignments.groupby(
            ["partition", "broad_state", f"fine_state_k{k}"], observed=True
        ).size()
        if (counts == 0).any():
            raise ValueError(f"An assigned fine state is empty at k={k}")

    write_tsv_gzip(assignments, OUT_DIR / "cell_fine_state_assignments.tsv.gz")
    write_tsv_gzip(features, OUT_DIR / "normal_reference_feature_statistics.tsv.gz")
    model_diagnostics.to_csv(
        OUT_DIR / "fine_state_model_diagnostics.tsv", sep="\t", index=False
    )
    centroids.to_csv(
        OUT_DIR / "normal_reference_cluster_centroids.tsv", sep="\t", index=False
    )

    audit = (
        assignments.groupby(
            ["partition", "route", "broad_state"], observed=True
        )
        .agg(
            n_cells=("cell_index", "size"),
            n_specimens=("specimen_id", "nunique"),
            n_donors=("donor_id", "nunique"),
            mean_programme_rank_score=("programme_rank_score", "mean"),
        )
        .reset_index()
    )
    audit.to_csv(OUT_DIR / "assignment_audit.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "state_shared_revision_define_fine_states_v2",
        "created_utc": pd.Timestamp.utcnow().isoformat(),
        "random_seed": SEED,
        "inputs_sha256": {
            "discovery_h5ad": sha256(H5AD["discovery"]),
            "validation_h5ad": sha256(H5AD["validation"]),
            "common_effects": sha256(COMMON_EFFECTS),
            "contract": sha256(CONTRACT),
        },
        "normal_reference": "discovery normal cells separately within ABS, GOB and TAC",
        "programme_genes_used_for_fine_state_definition": False,
        "hvg_selection": (
            "top 1,500 variance-ranked genes among genes detected in 1%-99% of "
            "normal-reference cells after fixed exclusions"
        ),
        "fine_state_resolutions": list(K_VALUES),
        "primary_resolution": PRIMARY_K,
        "programme_cell_score": (
            "mean within-cell percentile rank of frozen up genes minus mean rank "
            "of frozen down genes"
        ),
        "discovery_validation_overlap_removed": len(overlap),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
            "h5py": h5py.__version__,
        },
    }
    with (OUT_DIR / "analysis_manifest.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)

    print(
        f"Assigned {len(assignments):,} cells across three programme-blind "
        "normal-reference fine-state models."
    )


if __name__ == "__main__":
    main()
