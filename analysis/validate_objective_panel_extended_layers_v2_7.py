#!/usr/bin/env python3
"""Extended, label-blind validation of the frozen objective 12-gene panel.

The panel definition is read and checksummed before any validation resource is
opened.  No gene replacement, weight fitting, or cut-point optimization is
performed here.  This script adds the Becker RNA/ATAC, CRC Atlas, additional
WNT/ASCL2/pharmacologic perturbation, spatial, and leave-one-gene-out layers.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
OUT = ROOT / "results" / "objective_compact_panel_v2_7" / "extended_validation"
PANEL_PATH = (
    ROOT
    / "results"
    / "objective_compact_panel_v2_7"
    / "objective_compact_panel_frozen.tsv"
)
sys.path.insert(0, str(ANALYSIS))

import atlas_locked_study_influence as atlas_influence  # noqa: E402
import becker_locked_rna_atac_patient_robustness as becker_patient  # noqa: E402
import becker_multiome_regulatory_window_accessibility as regulatory  # noqa: E402
import becker_rna_atac_concordance as concordance  # noqa: E402
import computational_closure_validation as closure  # noqa: E402
import conventional_route_signature_transfer as transfer  # noqa: E402
import discovery_locked_route_signature as discovery  # noqa: E402
import external_sporadic_adenoma_validation as external  # noqa: E402
import translation_reduced_panel_v2_0 as reduced  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write(frame: pd.DataFrame, filename: str, compress: bool = False) -> None:
    frame.to_csv(
        OUT / filename,
        sep="\t",
        index=False,
        compression="gzip" if compress else None,
    )


def fixed_panel() -> pd.DataFrame:
    panel = pd.read_csv(PANEL_PATH, sep="\t")
    if len(panel) != 12:
        raise RuntimeError(f"Expected 12 frozen genes, found {len(panel)}")
    if panel["validation_outcomes_used"].astype(str).str.lower().ne("false").any():
        raise RuntimeError("Frozen panel records use of validation outcomes")
    if panel.groupby("arm")["gene"].size().to_dict() != {"down": 6, "up": 6}:
        raise RuntimeError("Frozen panel is not balanced at six genes per arm")
    panel = panel.copy()
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
    tests, models = transfer.becker_tests_and_models(scores)
    score_path = becker_dir / "becker_objective_panel_scores.tsv"
    scores.to_csv(score_path, sep="\t", index=False)
    availability.to_csv(
        becker_dir / "becker_objective_panel_gene_availability.tsv",
        sep="\t",
        index=False,
    )
    tests.to_csv(
        becker_dir / "becker_objective_panel_tests.tsv", sep="\t", index=False
    )
    models.to_csv(
        becker_dir / "becker_objective_panel_adjusted_models.tsv",
        sep="\t",
        index=False,
    )

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
    regulatory_corr.to_csv(
        rna_atac_dir / "becker_objective_regulatory_window_rna_correlations.tsv",
        sep="\t",
        index=False,
    )

    becker_patient.INPUT = rna_atac_dir / "becker_rna_atac_paired_scores.tsv"
    becker_patient.OUT_DIR = rna_atac_dir
    becker_patient.main()

    focus_test = tests.loc[
        (tests["comparison"] == "polyp_vs_normal_unaffected")
        & (tests["score"] == "epi__ca_route_signature")
    ]
    corr = pd.read_csv(
        rna_atac_dir / "becker_rna_atac_correlations.tsv", sep="\t"
    )
    focus_corr = corr.loc[
        (corr["subset"] == "normal_polyp")
        & (corr["analysis_id"] == "rna_epi_ca_route__atac_wnt_tcf_ascl2_axis")
    ]
    return {
        "becker_polyp_minus_normal": (
            float(focus_test.iloc[0]["delta_a_minus_b"]) if len(focus_test) else np.nan
        ),
        "becker_polyp_p": (
            float(focus_test.iloc[0]["p_mannwhitney"]) if len(focus_test) else np.nan
        ),
        "rna_atac_tcf_ascl2_rho": (
            float(focus_corr.iloc[0]["spearman_rho"]) if len(focus_corr) else np.nan
        ),
        "rna_atac_tcf_ascl2_p": (
            float(focus_corr.iloc[0]["p_spearman"]) if len(focus_corr) else np.nan
        ),
    }


def atlas_layers(panel: pd.DataFrame) -> dict[str, object]:
    atlas_dir = OUT / "crc_atlas"
    atlas_dir.mkdir(parents=True, exist_ok=True)
    donors, tests, models = transfer.score_atlas(panel)
    donor_path = atlas_dir / "atlas_objective_panel_donor_scores.tsv"
    donors.to_csv(donor_path, sep="\t", index=False)
    tests.to_csv(
        atlas_dir / "atlas_objective_panel_tests.tsv", sep="\t", index=False
    )
    models.to_csv(
        atlas_dir / "atlas_objective_panel_adjusted_models.tsv",
        sep="\t",
        index=False,
    )

    atlas_influence.INPUT = donor_path
    atlas_influence.OUT_DIR = atlas_dir
    atlas_influence.main()

    focus = tests.loc[
        (tests["comparison"] == "polyp_cancer_vs_normal_epithelial")
        & (tests["score"] == "ca_route_signature")
    ]
    influence = pd.read_csv(
        atlas_dir / "atlas_locked_study_influence.tsv", sep="\t"
    )
    influence = influence.loc[
        (influence["outcome"] == "score__ca_route_signature")
        & (influence["state"] == "polyp_cancer")
    ]
    return {
        "atlas_polyp_cancer_minus_normal": (
            float(focus.iloc[0]["delta_a_minus_b"]) if len(focus) else np.nan
        ),
        "atlas_polyp_cancer_p": (
            float(focus.iloc[0]["p_mannwhitney"]) if len(focus) else np.nan
        ),
        "atlas_polyp_loo_positive_fraction": (
            float(influence.iloc[0]["loo_positive_fraction"])
            if len(influence)
            else np.nan
        ),
    }


def perturbation_and_spatial_layers(panel: pd.DataFrame) -> dict[str, object]:
    perturb_dir = OUT / "perturbation_spatial"
    perturb_dir.mkdir(parents=True, exist_ok=True)
    signature = panel[["gene", "signature_direction", "route_weight"]].copy()
    homology_map, homology_table = closure.strict_human_mouse_map()
    datasets = {
        "GSE114059": (
            "human",
            "patient_derived_crc_organoid",
            *closure.load_gse114059(),
        ),
        "GSE67186": (
            "mouse",
            "in_vivo_colon_polyp",
            *closure.load_gse67186(homology_map),
        ),
        "GSE130822": (
            "mouse",
            "colonic_stem_cells",
            *closure.load_gse130822(homology_map),
        ),
        "GSE171910": (
            "human",
            "conditional_crc_wnt_models",
            *closure.load_gse171910(),
        ),
    }
    sample_parts: list[pd.DataFrame] = []
    unit_parts: list[pd.DataFrame] = []
    summary_parts: list[pd.DataFrame] = []
    coverage_parts: list[pd.DataFrame] = []
    matched_parts: list[pd.DataFrame] = []
    null_parts: list[pd.DataFrame] = []
    for dataset, (species, model_system, expression, metadata, comparisons) in datasets.items():
        scores, z, coverage = closure.score_expression(expression, signature, dataset)
        coverage.loc[coverage["feature"] == "route_up", "n_expected"] = 6
        coverage.loc[coverage["feature"] == "route_down", "n_expected"] = 6
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
        tests, null = closure.expression_matched_test(
            dataset, expression, z, signature, metadata, comparisons
        )
        matched_parts.append(tests)
        null_parts.append(null)

    spatial_spots, spatial_summary, spatial_units, spatial_tests = (
        closure.spatial_locked_route(signature)
    )
    sample_scores = pd.concat(sample_parts, ignore_index=True, sort=False)
    unit_effects = pd.concat(unit_parts, ignore_index=True, sort=False)
    effect_summary = pd.concat(summary_parts, ignore_index=True, sort=False)
    coverage = pd.concat(coverage_parts, ignore_index=True, sort=False)
    coverage["coverage_fraction"] = coverage["n_present"] / coverage["n_expected"]
    matched = pd.concat(matched_parts, ignore_index=True, sort=False)
    null = pd.concat(null_parts, ignore_index=True, sort=False)

    write(sample_scores, "perturbation_spatial/perturbation_sample_scores.tsv")
    write(unit_effects, "perturbation_spatial/perturbation_unit_effects.tsv")
    write(effect_summary, "perturbation_spatial/perturbation_effect_summary.tsv")
    write(coverage, "perturbation_spatial/feature_coverage.tsv")
    write(matched, "perturbation_spatial/expression_matched_tests.tsv")
    write(null, "perturbation_spatial/expression_matched_null.tsv.gz", compress=True)
    write(
        spatial_spots,
        "perturbation_spatial/spatial_objective_panel_spot_scores.tsv.gz",
        compress=True,
    )
    write(
        spatial_summary,
        "perturbation_spatial/spatial_objective_panel_pathology_summary.tsv",
    )
    write(
        spatial_units,
        "perturbation_spatial/spatial_objective_panel_section_effects.tsv",
    )
    write(spatial_tests, "perturbation_spatial/spatial_objective_panel_tests.tsv")
    write(
        homology_table,
        "perturbation_spatial/human_mouse_one_to_one_homology.tsv",
    )

    route = effect_summary.loc[effect_summary["feature"].eq("route_score")].copy()
    route["direction_supported"] = np.where(
        route["expected_direction"].eq(0),
        np.nan,
        route["expected_direction"] * route["mean_difference"] > 0,
    )
    spatial_focus = spatial_tests.loc[
        (spatial_tests["comparison"] == "tumor_vs_non_neoplastic_epithelium")
        & (spatial_tests["feature"] == "route_score")
    ]
    expected = route.loc[route["expected_direction"].ne(0)]
    return {
        "additional_perturbations_direction_supported": int(
            expected["direction_supported"].fillna(False).sum()
        ),
        "additional_perturbations_total": int(len(expected)),
        "spatial_tumor_minus_epithelium": (
            float(spatial_focus.iloc[0]["median_difference"])
            if len(spatial_focus)
            else np.nan
        ),
        "spatial_sections_positive": (
            int(spatial_focus.iloc[0]["n_positive"]) if len(spatial_focus) else 0
        ),
        "spatial_sections_total": (
            int(spatial_focus.iloc[0]["n_sections"]) if len(spatial_focus) else 0
        ),
        "spatial_p": (
            float(spatial_focus.iloc[0]["p_wilcoxon"])
            if len(spatial_focus)
            else np.nan
        ),
    }


def leave_one_gene_out(panel: pd.DataFrame) -> pd.DataFrame:
    working = panel.copy()
    working["panel_arm"] = working["arm"]
    external_data = reduced.external_gene_data(working)
    chen_meta, chen_expr = discovery.load_chen_pseudobulk("validation")
    _, _, full_external, _, _ = reduced.external_analysis(working)
    ffpe_scores, _, _, ffpe_expression = reduced.ffpe_analysis(working)
    result = reduced.leave_one_gene_out(
        working,
        external_data,
        chen_meta,
        chen_expr,
        ffpe_scores,
        ffpe_expression,
        float(full_external["adenoma_coef_sd"].iloc[0]),
    )
    return result.rename(
        columns={
            "external_retained_fraction_vs_10_gene":
                "external_retained_fraction_vs_primary_12_gene"
        }
    )


def qa_summary(
    becker: dict[str, object],
    atlas: dict[str, object],
    closure_result: dict[str, object],
    loo: pd.DataFrame,
) -> pd.DataFrame:
    rows = [
        {
            "check": "becker_polyp_effect_positive",
            "passed": bool(becker["becker_polyp_minus_normal"] > 0),
            "observed": becker["becker_polyp_minus_normal"],
        },
        {
            "check": "becker_rna_atac_axis_positive",
            "passed": bool(becker["rna_atac_tcf_ascl2_rho"] > 0),
            "observed": becker["rna_atac_tcf_ascl2_rho"],
        },
        {
            "check": "atlas_polyp_cancer_effect_positive",
            "passed": bool(atlas["atlas_polyp_cancer_minus_normal"] > 0),
            "observed": atlas["atlas_polyp_cancer_minus_normal"],
        },
        {
            "check": "atlas_polyp_leave_one_study_out_all_positive",
            "passed": bool(atlas["atlas_polyp_loo_positive_fraction"] == 1),
            "observed": atlas["atlas_polyp_loo_positive_fraction"],
        },
        {
            "check": "all_additional_non_null_perturbations_expected_direction",
            "passed": bool(
                closure_result["additional_perturbations_direction_supported"]
                == closure_result["additional_perturbations_total"]
            ),
            "observed": (
                f"{closure_result['additional_perturbations_direction_supported']}/"
                f"{closure_result['additional_perturbations_total']}"
            ),
        },
        {
            "check": "spatial_all_sections_positive",
            "passed": bool(
                closure_result["spatial_sections_positive"]
                == closure_result["spatial_sections_total"]
                and closure_result["spatial_sections_total"] > 0
            ),
            "observed": (
                f"{closure_result['spatial_sections_positive']}/"
                f"{closure_result['spatial_sections_total']}"
            ),
        },
        {
            "check": "leave_one_gene_out_external_effect_positive",
            "passed": bool((loo["external_pooled_effect_sd"] > 0).all()),
            "observed": float(loo["external_pooled_effect_sd"].min()),
        },
        {
            "check": "leave_one_gene_out_heldout_auc_at_least_0_85",
            "passed": bool((loo["heldout_auc"] >= 0.85).all()),
            "observed": float(loo["heldout_auc"].min()),
        },
        {
            "check": "leave_one_gene_out_ffpe_positive_fraction_at_least_0_85",
            "passed": bool((loo["ffpe_positive_fraction"] >= 0.85).all()),
            "observed": float(loo["ffpe_positive_fraction"].min()),
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checksum_before_validation = sha256(PANEL_PATH)
    panel = fixed_panel()

    becker = becker_layers(panel)
    atlas = atlas_layers(panel)
    closure_result = perturbation_and_spatial_layers(panel)
    loo = leave_one_gene_out(panel)
    write(loo, "objective_panel_leave_one_gene_out.tsv")
    qa = qa_summary(becker, atlas, closure_result, loo)
    write(qa, "extended_validation_qa.tsv")

    summary = pd.DataFrame(
        [
            {"layer": "becker", **becker},
            {"layer": "crc_atlas", **atlas},
            {"layer": "perturbation_spatial", **closure_result},
        ]
    )
    write(summary, "extended_validation_summary.tsv")
    manifest = {
        "analysis": "extended validation of frozen objective 12-gene panel",
        "panel_sha256_before_validation_read": checksum_before_validation,
        "panel_frozen_before_validation": True,
        "gene_reselection": False,
        "weight_fitting": False,
        "cutpoint_optimization": False,
        "validation_layers": [
            "Becker snRNA",
            "Becker matched RNA-ATAC",
            "CRC Atlas donor-level recurrence",
            "CRC Atlas leave-one-study-out",
            "additional WNT/ASCL2/pharmacologic perturbations",
            "six-section spatial pathology",
            "leave-one-gene-out robustness",
        ],
        "all_strict_qa_checks_passed": bool(qa["passed"].all()),
    }
    (OUT / "extended_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
