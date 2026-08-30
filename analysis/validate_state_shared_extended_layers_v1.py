#!/usr/bin/env python3
"""Orthogonal validation of the frozen eight-gene state-shared readout."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
RESULT_ROOT = ROOT / "results" / "state_aware_program_v1"
OUT = RESULT_ROOT / "extended_validation"
PANEL_PATH = RESULT_ROOT / "panel_derivation" / "compact_state_shared_panel_frozen.tsv"
ADDENDUM_PATH = (
    ANALYSIS / "contracts" / "state_shared_extended_validation_addendum_v1_2026-08-29.md"
)
sys.path.insert(0, str(ANALYSIS))

import atlas_locked_study_influence as atlas_influence  # noqa: E402
import becker_locked_rna_atac_patient_robustness as becker_patient  # noqa: E402
import becker_multiome_regulatory_window_accessibility as regulatory  # noqa: E402
import becker_rna_atac_concordance as concordance  # noqa: E402
import computational_closure_validation as closure  # noqa: E402
import conventional_route_signature_transfer as transfer  # noqa: E402


EXPECTED_HASHES = {
    "compact_panel": "c5997e572342a72da8441df312fba4e3461cacfa5e30d0d8590a3f23ae3d96f0",
    "extended_addendum": "52d638412a71c624177af9c0a0da14e10d0d970abd66823662e25dae8747ebea",
}


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


def load_panel() -> pd.DataFrame:
    actual_hashes = {
        "compact_panel": sha256(PANEL_PATH),
        "extended_addendum": sha256(ADDENDUM_PATH),
    }
    if actual_hashes != EXPECTED_HASHES:
        raise RuntimeError(f"Frozen extended-validation input changed: {actual_hashes}")
    panel = pd.read_csv(PANEL_PATH, sep="\t").copy()
    if len(panel) != 8 or panel["arm"].value_counts().to_dict() != {"up": 4, "down": 4}:
        raise RuntimeError("Frozen compact readout has unexpected dimensions")
    if panel["validation_outcomes_used"].astype(str).str.lower().ne("false").any():
        raise RuntimeError("Frozen compact readout records use of validation outcomes")
    panel["signature_direction"] = np.where(
        panel["arm"].eq("up"), "adenoma_up", "adenoma_down"
    )
    return panel


def becker_layers(panel: pd.DataFrame) -> dict[str, object]:
    becker_dir = OUT / "becker"
    rna_atac_dir = OUT / "becker_rna_atac"
    becker_dir.mkdir(parents=True, exist_ok=True)
    rna_atac_dir.mkdir(parents=True, exist_ok=True)

    scores, availability = transfer.score_becker(panel)
    panel_genes = set(panel["gene"])
    availability = availability.copy()
    availability["panel_genes_expected"] = len(panel_genes)
    availability["panel_genes_missing"] = availability["missing_genes"].fillna("").map(
        lambda value: ";".join(sorted(panel_genes.intersection(str(value).split(","))))
    )
    availability["panel_genes_present"] = availability["panel_genes_expected"] - availability[
        "panel_genes_missing"
    ].map(lambda value: 0 if value == "" else len(value.split(";")))
    if availability["panel_genes_present"].min() != 8:
        raise RuntimeError("Becker data do not contain all eight frozen genes")
    tests, models = transfer.becker_tests_and_models(scores)
    score_path = becker_dir / "becker_compact_8_scores.tsv"
    scores.to_csv(score_path, sep="\t", index=False)
    write(availability, "becker/becker_compact_8_gene_availability.tsv")
    write(tests, "becker/becker_compact_8_tests.tsv")
    write(models, "becker/becker_compact_8_adjusted_models.tsv")

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
        "becker_rna_atac/becker_compact_8_regulatory_window_rna_correlations.tsv",
    )

    becker_patient.INPUT = rna_atac_dir / "becker_rna_atac_paired_scores.tsv"
    becker_patient.OUT_DIR = rna_atac_dir
    becker_patient.main()

    focus_test = tests.loc[
        (tests["comparison"] == "polyp_vs_normal_unaffected")
        & (tests["score"] == "epi__ca_route_signature")
    ]
    correlations = pd.read_csv(rna_atac_dir / "becker_rna_atac_correlations.tsv", sep="\t")
    focus_corr = correlations.loc[
        (correlations["subset"] == "normal_polyp")
        & (correlations["analysis_id"] == "rna_epi_ca_route__atac_wnt_tcf_ascl2_axis")
    ]
    return {
        "becker_polyp_minus_normal": float(focus_test.iloc[0]["delta_a_minus_b"]),
        "becker_polyp_p": float(focus_test.iloc[0]["p_mannwhitney"]),
        "rna_atac_axis_rho": float(focus_corr.iloc[0]["spearman_rho"]),
        "rna_atac_axis_p": float(focus_corr.iloc[0]["p_spearman"]),
    }


def atlas_layers(panel: pd.DataFrame) -> dict[str, object]:
    atlas_dir = OUT / "crc_atlas"
    atlas_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(transfer.ATLAS_PATH, "r") as handle:
        genes = set(map(str, transfer.read_h5_col(handle["var"], "feature_name")))
    missing = sorted(set(panel["gene"]).difference(genes))
    coverage = pd.DataFrame(
        [{"expected": 8, "present": 8 - len(missing), "missing_genes": ";".join(missing)}]
    )
    write(coverage, "crc_atlas/atlas_compact_8_gene_availability.tsv")
    if missing:
        raise RuntimeError(f"CRC Atlas is missing frozen genes: {missing}")

    donors, tests, models = transfer.score_atlas(panel)
    donor_path = atlas_dir / "atlas_compact_8_donor_scores.tsv"
    donors.to_csv(donor_path, sep="\t", index=False)
    write(tests, "crc_atlas/atlas_compact_8_tests.tsv")
    write(models, "crc_atlas/atlas_compact_8_adjusted_models.tsv")

    atlas_influence.INPUT = donor_path
    atlas_influence.OUT_DIR = atlas_dir
    atlas_influence.main()
    focus = tests.loc[
        (tests["comparison"] == "polyp_cancer_vs_normal_epithelial")
        & (tests["score"] == "ca_route_signature")
    ]
    influence = pd.read_csv(atlas_dir / "atlas_locked_study_influence.tsv", sep="\t")
    influence = influence.loc[
        (influence["outcome"] == "score__ca_route_signature")
        & (influence["state"] == "polyp_cancer")
    ]
    return {
        "atlas_polyp_cancer_minus_normal": float(focus.iloc[0]["delta_a_minus_b"]),
        "atlas_polyp_cancer_p": float(focus.iloc[0]["p_mannwhitney"]),
        "atlas_polyp_leave_one_study_out_positive_fraction": float(
            influence.iloc[0]["loo_positive_fraction"]
        ),
    }


def perturbation_and_spatial_layers(panel: pd.DataFrame) -> dict[str, object]:
    signature = panel[["gene", "signature_direction", "route_weight"]].copy()
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
    matched_parts: list[pd.DataFrame] = []
    null_parts: list[pd.DataFrame] = []
    for dataset, (species, model_system, expression, metadata, comparisons) in datasets.items():
        scores, z_scores, coverage = closure.score_expression(expression, signature, dataset)
        coverage.loc[coverage["feature"].eq("route_up"), "n_expected"] = 4
        coverage.loc[coverage["feature"].eq("route_down"), "n_expected"] = 4
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
        matched, null = closure.expression_matched_test(
            dataset, expression, z_scores, signature, metadata, comparisons
        )
        matched_parts.append(matched)
        null_parts.append(null)

    spatial_spots, spatial_summary, spatial_units, spatial_tests = closure.spatial_locked_route(
        signature
    )
    samples = pd.concat(sample_parts, ignore_index=True, sort=False)
    units = pd.concat(unit_parts, ignore_index=True, sort=False)
    summaries = pd.concat(summary_parts, ignore_index=True, sort=False)
    coverage = pd.concat(coverage_parts, ignore_index=True, sort=False)
    coverage["coverage_fraction"] = coverage["n_present"] / coverage["n_expected"]
    matched = pd.concat(matched_parts, ignore_index=True, sort=False)
    null = pd.concat(null_parts, ignore_index=True, sort=False)

    write(samples, "perturbation_spatial/perturbation_sample_scores.tsv")
    write(units, "perturbation_spatial/perturbation_unit_effects.tsv")
    write(summaries, "perturbation_spatial/perturbation_effect_summary.tsv")
    write(coverage, "perturbation_spatial/feature_coverage.tsv")
    write(matched, "perturbation_spatial/expression_matched_tests.tsv")
    write(null, "perturbation_spatial/expression_matched_null.tsv.gz", compress=True)
    write(spatial_spots, "perturbation_spatial/spatial_compact_8_spot_scores.tsv.gz", compress=True)
    write(spatial_summary, "perturbation_spatial/spatial_compact_8_pathology_summary.tsv")
    write(spatial_units, "perturbation_spatial/spatial_compact_8_section_effects.tsv")
    write(spatial_tests, "perturbation_spatial/spatial_compact_8_tests.tsv")
    write(homology_table, "perturbation_spatial/human_mouse_one_to_one_homology.tsv")

    route = summaries.loc[summaries["feature"].eq("route_score")].copy()
    expected = route.loc[route["expected_direction"].ne(0)].copy()
    expected["direction_supported"] = (
        expected["expected_direction"] * expected["mean_difference"] > 0
    )
    spatial_focus = spatial_tests.loc[
        (spatial_tests["comparison"] == "tumor_vs_non_neoplastic_epithelium")
        & (spatial_tests["feature"] == "route_score")
    ]
    return {
        "additional_perturbations_direction_supported": int(expected["direction_supported"].sum()),
        "additional_perturbations_total": int(len(expected)),
        "spatial_tumor_minus_epithelium": float(spatial_focus.iloc[0]["median_difference"]),
        "spatial_sections_positive": int(spatial_focus.iloc[0]["n_positive"]),
        "spatial_sections_total": int(spatial_focus.iloc[0]["n_sections"]),
        "spatial_p": float(spatial_focus.iloc[0]["p_wilcoxon"]),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    panel = load_panel()
    becker = becker_layers(panel)
    atlas = atlas_layers(panel)
    perturbation_spatial = perturbation_and_spatial_layers(panel)
    summary = pd.DataFrame(
        [
            {"layer": "becker", **becker},
            {"layer": "crc_atlas", **atlas},
            {"layer": "perturbation_spatial", **perturbation_spatial},
        ]
    )
    write(summary, "extended_validation_summary.tsv")
    manifest = {
        "analysis": "validate_state_shared_extended_layers_v1",
        "analysis_date": "2026-08-29",
        "frozen_input_hashes": {
            "compact_panel": sha256(PANEL_PATH),
            "extended_addendum": sha256(ADDENDUM_PATH),
        },
        "compact_readout_genes": len(panel),
        "gene_reselection": False,
        "weight_fitting": False,
        "cutpoint_optimization": False,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    (OUT / "extended_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"output_dir\t{OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
