#!/usr/bin/env python3
"""Outcome-blind feature-coverage audit for the frozen 1,843-gene programme."""

from __future__ import annotations

import gzip
import io
import json
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = ROOT / "results" / "state_aware_program_v1"
PROGRAM_PATH = RESULT_ROOT / "common_effects" / "cross_state_common_effects.tsv.gz"
OUT_DIR = RESULT_ROOT / "full_program_coverage_audit"


def load_programme() -> pd.DataFrame:
    frame = pd.read_csv(PROGRAM_PATH, sep="\t")
    strict = frame["strict_state_shared"]
    if strict.dtype != bool:
        strict = strict.astype(str).str.lower().eq("true")
    programme = frame.loc[strict, ["gene", "shared_direction"]].drop_duplicates("gene")
    if len(programme) != 1843:
        raise RuntimeError(f"Expected 1,843 frozen genes, found {len(programme)}")
    if programme["shared_direction"].value_counts().to_dict() != {"down": 959, "up": 884}:
        raise RuntimeError("Frozen programme arm dimensions changed")
    return programme


def coverage_row(layer: str, unit: str, genes: set[str], programme: pd.DataFrame) -> dict[str, object]:
    up = set(programme.loc[programme["shared_direction"].eq("up"), "gene"])
    down = set(programme.loc[programme["shared_direction"].eq("down"), "gene"])
    return {
        "layer": layer,
        "unit": unit,
        "expected_up": len(up),
        "present_up": len(up & genes),
        "coverage_up": len(up & genes) / len(up),
        "expected_down": len(down),
        "present_down": len(down & genes),
        "coverage_down": len(down & genes) / len(down),
        "present_total": len((up | down) & genes),
        "coverage_total": len((up | down) & genes) / len(up | down),
    }


def decode(values: np.ndarray) -> set[str]:
    return {
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    }


def audit_becker(programme: pd.DataFrame) -> list[dict[str, object]]:
    tar_path = ROOT / "data_sources" / "Becker_NatGenet_2022_GEO" / "GSE201348_RAW_scRNA.tar"
    rows: list[dict[str, object]] = []
    with tarfile.open(tar_path, "r") as archive:
        feature_members = sorted(
            (member for member in archive.getmembers() if member.name.endswith("_features.tsv.gz")),
            key=lambda member: member.name,
        )
        for member in feature_members:
            raw = archive.extractfile(member)
            if raw is None:
                raise FileNotFoundError(member.name)
            with gzip.GzipFile(fileobj=raw) as compressed:
                genes = {
                    (parts[1] if len(parts) > 1 else parts[0])
                    for line in io.TextIOWrapper(compressed, encoding="utf-8", errors="replace")
                    if (parts := line.rstrip("\n").split("\t"))
                }
            rows.append(coverage_row("Becker epithelial snRNA-seq", member.name, genes, programme))
    return rows


def audit_atlas(programme: pd.DataFrame) -> list[dict[str, object]]:
    path = ROOT / "data_sources" / "CRC_Atlas_CZI_core" / "crc_atlas_core_czi.h5ad"
    with h5py.File(path, "r") as handle:
        values = handle["var"]["feature_name"]
        if isinstance(values, h5py.Group):
            codes = values["codes"][()]
            categories = values["categories"][()]
            genes = decode(categories[codes[codes >= 0]])
        else:
            genes = decode(values[()])
    return [coverage_row("CRC Atlas scRNA-seq", path.name, genes, programme)]


def audit_perturbations(programme: pd.DataFrame) -> list[dict[str, object]]:
    import sys

    sys.path.insert(0, str(ROOT / "analysis"))
    import computational_closure_validation as closure

    homology, _ = closure.strict_human_mouse_map()
    loaders = {
        "GSE114059": closure.load_gse114059,
        "GSE67186": lambda: closure.load_gse67186(homology),
        "GSE130822": lambda: closure.load_gse130822(homology),
        "GSE171910": closure.load_gse171910,
        "GSE125472": closure.load_gse125_expression,
    }
    rows = []
    for dataset, loader in loaders.items():
        expression, _, _ = loader()
        rows.append(coverage_row("Genetic/pharmacological perturbation", dataset, set(expression.index.astype(str)), programme))
    return rows


def audit_spatial(programme: pd.DataFrame) -> list[dict[str, object]]:
    spatial_root = ROOT / "data_sources" / "Zenodo14602110_spatial" / "Data"
    if not spatial_root.exists():
        candidates = sorted((ROOT / "data_sources").glob("**/filtered_feature_bc_matrix.h5"))
    else:
        candidates = sorted(spatial_root.glob("*/filtered_feature_bc_matrix.h5"))
    rows = []
    for path in candidates:
        with h5py.File(path, "r") as handle:
            genes = decode(handle["matrix/features/name"][()])
        rows.append(coverage_row("Visium spatial transcriptomics", path.parent.name, genes, programme))
    return rows


def summarize(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame.groupby("layer", observed=True)
        .agg(
            n_units=("unit", "size"),
            minimum_up_coverage=("coverage_up", "min"),
            median_up_coverage=("coverage_up", "median"),
            minimum_down_coverage=("coverage_down", "min"),
            median_down_coverage=("coverage_down", "median"),
            minimum_total_coverage=("coverage_total", "min"),
            median_total_coverage=("coverage_total", "median"),
        )
        .reset_index()
    )


def main() -> None:
    programme = load_programme()
    rows = []
    rows.extend(audit_becker(programme))
    rows.extend(audit_atlas(programme))
    rows.extend(audit_perturbations(programme))
    rows.extend(audit_spatial(programme))
    detail = pd.DataFrame(rows)
    if detail.empty:
        raise RuntimeError("No platform feature inventories were found")
    summary = summarize(detail)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    detail.to_csv(OUT_DIR / "full_programme_feature_coverage_by_unit.tsv", sep="\t", index=False)
    summary.to_csv(OUT_DIR / "full_programme_feature_coverage_summary.tsv", sep="\t", index=False)
    manifest = {
        "analysis": "audit_state_shared_full_program_coverage_v1",
        "decision_data": "feature identities only; no phenotypes, effects, or significance results",
        "programme_genes": 1843,
        "programme_arms": {"up": 884, "down": 959},
        "n_units": int(len(detail)),
    }
    (OUT_DIR / "audit_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
