#!/usr/bin/env python3
"""Renumber the six legacy main figures after insertion of the new Figure 3."""

from __future__ import annotations

import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "figures" / "communications_biology_v1.2"
MOVES = [
    ("figure6_spatial_and_protein_context", "figure7_spatial_and_protein_context"),
    (
        "figure5_empirical_and_virtual_perturbation_support",
        "figure6_empirical_and_virtual_perturbation_support",
    ),
    (
        "figure4_crc_atlas_cross_sectional_recurrence",
        "figure5_crc_atlas_cross_sectional_recurrence",
    ),
    ("figure3_rna_atac_regulatory_support", "figure4_rna_atac_regulatory_support"),
]


def main() -> None:
    copied = 0
    for old_stem, new_stem in MOVES:
        for suffix in ("pdf", "png", "svg", "tiff"):
            source = FIGURE_DIR / f"{old_stem}.{suffix}"
            if source.is_file():
                shutil.copy2(source, FIGURE_DIR / f"{new_stem}.{suffix}")
                copied += 1
    for old_stem, _ in MOVES:
        for suffix in ("pdf", "png", "svg", "tiff"):
            obsolete = FIGURE_DIR / f"{old_stem}.{suffix}"
            if obsolete.is_file():
                obsolete.unlink()
    if copied < 12:
        raise RuntimeError(f"Expected at least 12 distributed figure files, copied {copied}")
    print(f"Renumbered downstream figures after Figure 3 insertion: {copied} files")


if __name__ == "__main__":
    main()
