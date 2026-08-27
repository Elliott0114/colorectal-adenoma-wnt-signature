#!/usr/bin/env python3
"""Audit the seven-main-figure Communications Biology release package."""

from __future__ import annotations

import csv
import hashlib
import struct
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "figures" / "communications_biology_v1.2"
SOURCE_DIR = ROOT / "data" / "source_data"
RESULT_DIR = ROOT / "results" / "cell_state_decomposition_v1"
AUDIT_DIR = FIGURE_DIR / "source_data"

STEMS = [
    "figure1_discovery_core_and_objective_reduction",
    "figure2_independent_replication_and_ffpe",
    "figure3_cell_state_decomposition",
    "figure4_rna_atac_regulatory_support",
    "figure5_crc_atlas_cross_sectional_recurrence",
    "figure6_empirical_and_virtual_perturbation_support",
    "figure7_spatial_and_protein_context",
    "figureS1_core_composition_and_portability",
    "figureS2_external_and_ffpe_sensitivity",
    "figureS3_signature_transparency_and_random_benchmark",
    "figureS4_rna_atac_robustness",
    "figureS5_crc_atlas_source_audit",
    "figureS6_perturbation_boundaries",
    "figureS7_virtual_knockout_robustness",
    "figureS8_spatial_and_protein_assayability",
]

FIGURE3_SOURCES = [
    "figure3a_epithelial_proportions.tsv",
    "figure3b_within_state_effects.tsv",
    "figure3c_paired_donor_scores.tsv",
    "figure3d_state_contributions.tsv",
    "figure3e_exact_decomposition.tsv",
    "figure3f_core_compact_concordance.tsv",
    "figure3_source_data_manifest.tsv",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def png_metadata(path: Path) -> tuple[int, int, float | None]:
    with path.open("rb") as handle:
        if handle.read(8) != b"\x89PNG\r\n\x1a\n":
            raise ValueError(f"Invalid PNG signature: {path}")
        width = height = None
        dpi = None
        while True:
            length_raw = handle.read(4)
            if not length_raw:
                break
            length = struct.unpack(">I", length_raw)[0]
            chunk_type = handle.read(4)
            data = handle.read(length)
            handle.read(4)
            if chunk_type == b"IHDR":
                width, height = struct.unpack(">II", data[:8])
            elif chunk_type == b"pHYs" and len(data) == 9:
                pixels_x, _, unit = struct.unpack(">IIB", data)
                if unit == 1:
                    dpi = pixels_x * 0.0254
            elif chunk_type == b"IEND":
                break
    if width is None or height is None:
        raise ValueError(f"PNG lacks IHDR: {path}")
    return width, height, dpi


def tsv_all_true(path: Path, column: str) -> bool:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return bool(rows) and all(str(row[column]).lower() in {"true", "1"} for row in rows)


def main() -> None:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []
    for stem in STEMS:
        for suffix in ("pdf", "png", "svg"):
            path = FIGURE_DIR / f"{stem}.{suffix}"
            exists = path.is_file()
            width = height = ""
            dpi: float | str = ""
            valid = exists and path.stat().st_size > 10_000
            if exists and suffix == "png":
                width, height, dpi_value = png_metadata(path)
                dpi = "" if dpi_value is None else round(dpi_value, 2)
                valid = valid and width >= 1200 and height >= 800
                valid = valid and dpi_value is not None and dpi_value >= 295
            elif exists and suffix == "pdf":
                valid = valid and path.read_bytes()[:5] == b"%PDF-"
            elif exists and suffix == "svg":
                try:
                    root = ET.parse(path).getroot()
                    valid = valid and root.tag.endswith("svg")
                except ET.ParseError:
                    valid = False
            manifest.append(
                {
                    "figure": stem,
                    "format": suffix.upper(),
                    "path": str(path.relative_to(ROOT)),
                    "width_pixels": width,
                    "height_pixels": height,
                    "dpi": dpi,
                    "bytes": path.stat().st_size if exists else 0,
                    "sha256": sha256(path) if exists else "",
                    "pass": valid,
                }
            )

    manifest_path = AUDIT_DIR / "figure_export_manifest_v1_3.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(manifest[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(manifest)

    checks = [
        {
            "check": "15 figures in PDF, PNG and SVG",
            "value": len(manifest),
            "criterion": "45 files",
            "pass": len(manifest) == 45 and all(bool(row["pass"]) for row in manifest),
        },
        {
            "check": "Figure 3 panel source data",
            "value": sum((SOURCE_DIR / name).is_file() for name in FIGURE3_SOURCES),
            "criterion": "7 files",
            "pass": all((SOURCE_DIR / name).is_file() for name in FIGURE3_SOURCES),
        },
        {
            "check": "Cell-state analysis QA",
            "value": "all rows",
            "criterion": "pass",
            "pass": tsv_all_true(RESULT_DIR / "analysis_qa.tsv", "pass"),
        },
        {
            "check": "Independent R numerical audit",
            "value": "16 decompositions",
            "criterion": "pass",
            "pass": tsv_all_true(RESULT_DIR / "independent_r_audit.tsv", "pass"),
        },
    ]
    summary_path = AUDIT_DIR / "figure_package_audit_summary_v1_3.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(checks[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(checks)

    if not all(bool(row["pass"]) for row in checks):
        raise SystemExit("Figure package v1.3 audit failed")
    print(f"Figure package v1.3 audit passed: {len(checks)}/{len(checks)} checks")


if __name__ == "__main__":
    main()
