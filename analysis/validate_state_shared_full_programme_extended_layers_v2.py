#!/usr/bin/env python3
"""Project the frozen 1,843-gene programme across orthogonal validation layers."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
RESULT_ROOT = ROOT / "results" / "state_aware_program_v1"
OUT = RESULT_ROOT / "extended_validation_full_programme"
COMMON_PATH = RESULT_ROOT / "common_effects" / "cross_state_common_effects.tsv.gz"
COVERAGE_PATH = (
    RESULT_ROOT
    / "full_program_coverage_audit"
    / "full_programme_feature_coverage_by_unit.tsv"
)
ADDENDUM_PATH = (
    ANALYSIS
    / "contracts"
    / "state_shared_full_programme_projection_addendum_v1_2026-08-30.md"
)
sys.path.insert(0, str(ANALYSIS))

import atlas_locked_study_influence as atlas_influence  # noqa: E402
import becker_locked_rna_atac_patient_robustness as becker_patient  # noqa: E402
import becker_multiome_regulatory_window_accessibility as regulatory  # noqa: E402
import becker_rna_atac_concordance as concordance  # noqa: E402
import computational_closure_validation as closure  # noqa: E402
import conventional_route_signature_transfer as transfer  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write(frame: pd.DataFrame, relative_path: str, compress: bool = False) -> None:
    path = OUT / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        compression="gzip" if compress else None,
    )


def load_programme() -> pd.DataFrame:
    common = pd.read_csv(COMMON_PATH, sep="\t")
    strict = common["strict_state_shared"]
    if strict.dtype != bool:
        strict = strict.astype(str).str.lower().eq("true")
    programme = common.loc[strict, ["gene", "shared_direction"]].copy()
    programme["signature_direction"] = np.where(
        programme["shared_direction"].eq("up"), "adenoma_up", "adenoma_down"
    )
    programme["route_weight"] = np.where(programme["shared_direction"].eq("up"), 1.0, -1.0)
    if len(programme) != 1843:
        raise RuntimeError(f"Expected 1,843 genes, found {len(programme)}")
    return programme


def assert_coverage_rule() -> pd.DataFrame:
    coverage = pd.read_csv(COVERAGE_PATH, sep="\t")
    failed = coverage.loc[
        (coverage["coverage_up"] < 0.75)
        | (coverage["coverage_down"] < 0.75)
        | (coverage["present_up"] < 100)
        | (coverage["present_down"] < 100)
    ]
    if not failed.empty:
        raise RuntimeError(f"Outcome-blind projection rule failed:\n{failed.to_string(index=False)}")
    return coverage


def becker_layers(programme: pd.DataFrame) -> None:
    becker_dir = OUT / "becker"
    rna_atac_dir = OUT / "becker_rna_atac"
    becker_dir.mkdir(parents=True, exist_ok=True)
    rna_atac_dir.mkdir(parents=True, exist_ok=True)

    scores, availability = transfer.score_becker(programme)
    score_path = becker_dir / "becker_full_programme_scores.tsv"
    scores.to_csv(score_path, sep="\t", index=False)
    write(availability, "becker/becker_full_programme_feature_availability.tsv")
    tests, models = transfer.becker_tests_and_models(scores)
    write(tests, "becker/becker_full_programme_tests.tsv")
    write(models, "becker/becker_full_programme_adjusted_models.tsv")

    concordance.RNA_SIGNATURE_PATH = score_path
    concordance.OUT_DIR = rna_atac_dir
    concordance.main()

    regulatory.RNA_ATAC_DIR = rna_atac_dir
    regulatory_scores = pd.read_csv(
        ROOT
        / "results"
        / "becker_multiome_regulatory_windows"
        / "becker_multiome_regulatory_window_distance_bin_scores.tsv",
        sep="\t",
    )
    regulatory_corr = regulatory.rna_atac_correlations(regulatory_scores)
    write(
        regulatory_corr,
        "becker_rna_atac/becker_full_programme_regulatory_window_rna_correlations.tsv",
    )

    becker_patient.INPUT = rna_atac_dir / "becker_rna_atac_paired_scores.tsv"
    becker_patient.OUT_DIR = rna_atac_dir
    becker_patient.main()


def atlas_layers(programme: pd.DataFrame) -> None:
    atlas_dir = OUT / "crc_atlas"
    atlas_dir.mkdir(parents=True, exist_ok=True)
    donors, tests, models = transfer.score_atlas(programme)
    donor_path = atlas_dir / "atlas_full_programme_donor_scores.tsv"
    donors.to_csv(donor_path, sep="\t", index=False)
    write(tests, "crc_atlas/atlas_full_programme_tests.tsv")
    write(models, "crc_atlas/atlas_full_programme_adjusted_models.tsv")

    atlas_influence.INPUT = donor_path
    atlas_influence.OUT_DIR = atlas_dir
    atlas_influence.main()


def perturbation_and_spatial_layers(programme: pd.DataFrame) -> None:
    signature = programme[["gene", "signature_direction", "route_weight"]].copy()
    homology_map, homology_table = closure.strict_human_mouse_map()
    datasets = {
        "GSE114059": ("human", "patient_derived_crc_organoid", *closure.load_gse114059()),
        "GSE67186": ("mouse", "in_vivo_colon_polyp", *closure.load_gse67186(homology_map)),
        "GSE130822": ("mouse", "colonic_stem_cells", *closure.load_gse130822(homology_map)),
        "GSE171910": ("human", "conditional_crc_wnt_models", *closure.load_gse171910()),
    }
    sample_parts: list[pd.DataFrame] = []
    unit_parts: list[pd.DataFrame] = []
    summary_parts: list[pd.DataFrame] = []
    coverage_parts: list[pd.DataFrame] = []
    expected_up = int((signature["route_weight"] == 1).sum())
    expected_down = int((signature["route_weight"] == -1).sum())
    for dataset, (species, model_system, expression, metadata, comparisons) in datasets.items():
        scores, _, coverage = closure.score_expression(expression, signature, dataset)
        coverage.loc[coverage["feature"].eq("route_up"), "n_expected"] = expected_up
        coverage.loc[coverage["feature"].eq("route_down"), "n_expected"] = expected_down
        merged, units, summaries = closure.comparison_effects(
            dataset, species, model_system, scores, metadata, comparisons
        )
        merged.insert(0, "dataset", dataset)
        merged.insert(1, "species", species)
        merged.insert(2, "model_system", model_system)
        sample_parts.append(merged)
        unit_parts.append(units)
        summary_parts.append(summaries)
        coverage_parts.append(coverage)

    spatial_spots, spatial_summary, spatial_units, spatial_tests = closure.spatial_locked_route(signature)
    samples = pd.concat(sample_parts, ignore_index=True, sort=False)
    units = pd.concat(unit_parts, ignore_index=True, sort=False)
    summaries = pd.concat(summary_parts, ignore_index=True, sort=False)
    coverage = pd.concat(coverage_parts, ignore_index=True, sort=False)
    coverage["coverage_fraction"] = coverage["n_present"] / coverage["n_expected"]

    write(samples, "perturbation_spatial/perturbation_full_programme_sample_scores.tsv")
    write(units, "perturbation_spatial/perturbation_full_programme_unit_effects.tsv")
    write(summaries, "perturbation_spatial/perturbation_full_programme_effect_summary.tsv")
    write(coverage, "perturbation_spatial/perturbation_full_programme_feature_coverage.tsv")
    write(spatial_spots, "perturbation_spatial/spatial_full_programme_spot_scores.tsv.gz", compress=True)
    write(spatial_summary, "perturbation_spatial/spatial_full_programme_pathology_summary.tsv")
    write(spatial_units, "perturbation_spatial/spatial_full_programme_section_effects.tsv")
    write(spatial_tests, "perturbation_spatial/spatial_full_programme_tests.tsv")
    write(homology_table, "perturbation_spatial/human_mouse_one_to_one_homology.tsv")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    coverage = assert_coverage_rule()
    programme = load_programme()
    becker_layers(programme)
    atlas_layers(programme)
    perturbation_and_spatial_layers(programme)
    manifest = {
        "analysis": "validate_state_shared_full_programme_extended_layers_v2",
        "analysis_date": "2026-08-30",
        "programme_genes": 1843,
        "gene_reselection": False,
        "weight_fitting": False,
        "cutpoint_optimization": False,
        "outcome_blind_coverage_rule_passed": True,
        "minimum_observed_up_coverage": float(coverage["coverage_up"].min()),
        "minimum_observed_down_coverage": float(coverage["coverage_down"].min()),
        "input_hashes": {
            "common_effects": sha256(COMMON_PATH),
            "coverage_audit": sha256(COVERAGE_PATH),
            "analysis_addendum": sha256(ADDENDUM_PATH),
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    (OUT / "extended_validation_full_programme_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"output_dir\t{OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
