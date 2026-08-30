#!/usr/bin/env python3
"""Fast integrity checks for the public identity-remodelling release."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SIGNATURE_DIR = ROOT / "data" / "signatures"
FIGURE_DIR = ROOT / "figures" / "communications_biology_v3.0"
MAX_PUBLIC_FILE_BYTES = 50 * 1024 * 1024

MAIN_FIGURES = [
    "figure1_study_design_and_programme_derivation",
    "figure2_donor_disjoint_identity_remodelling",
    "figure3_external_recurrence_and_archival_transfer",
    "figure4_functional_architecture",
    "figure5_full_programme_orthogonal_context",
    "figure6_full_programme_perturbation_support",
]
SUPP_FIGURES = [
    "figureS1_sampling_and_donor_stability",
    "figureS2_fine_state_and_composition_sensitivities",
    "figureS3_functional_and_regulatory_structure",
    "figureS4_external_and_ffpe_sensitivities",
    "figureS5_reduced_readout_audit",
    "figureS6_multiomic_atlas_spatial_protein",
    "figureS7_empirical_perturbation_context",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def close(observed: float, expected: float, tolerance: float = 1e-10) -> bool:
    return abs(observed - expected) <= tolerance


def check_definitions() -> None:
    ranking = pd.read_csv(SIGNATURE_DIR / "common_effect_ranking_8221.tsv.gz", sep="\t")
    programme = pd.read_csv(SIGNATURE_DIR / "state_aware_high_confidence_1843.tsv", sep="\t")
    candidates = pd.read_csv(SIGNATURE_DIR / "portable_candidates_53.tsv", sep="\t")
    reduced = pd.read_csv(SIGNATURE_DIR / "compact_candidate_8_genes.tsv", sep="\t")
    require(len(ranking) == 8221, "Expected an 8,221-gene common-effect ranking")
    require(len(programme) == 1843, "Expected 1,843 high-confidence genes")
    require(
        programme["shared_direction"].value_counts().to_dict()
        == {"down": 959, "up": 884},
        "Unexpected programme arms",
    )
    require(len(candidates) == 53, "Expected 53 measurable candidates")
    require(len(reduced) == 8, "Expected eight reduced-readout genes")
    require(
        reduced["arm"].value_counts().to_dict() == {"up": 4, "down": 4},
        "Reduced candidate is not balanced",
    )
    require(
        set(reduced["gene"])
        == {"EPHB2", "REG1A", "LTBP1", "RNF43", "CALM2", "COX6C", "B2M", "ACAA2"},
        "Reduced membership changed",
    )
    require(set(reduced["gene"]).issubset(set(candidates["gene"])), "Reduced genes fall outside the measurable universe")
    require(not reduced["validation_outcomes_used"].astype(bool).any(), "Validation leakage flag changed")
    require(reduced["panel_frozen_before_validation"].astype(bool).all(), "Freeze flag changed")


def check_principal_results() -> None:
    donor = pd.read_csv(
        ROOT / "results/state_shared_revision_v2/donor_site/donor_disjoint_replication_summary.tsv",
        sep="\t",
    ).iloc[0]
    require(int(donor["n_testable_strict_genes"]) == 1646, "Donor-disjoint testable count changed")
    require(close(donor["common_direction_match_fraction"], 0.986026731470231), "Direction agreement changed")
    require(close(donor["discovery_validation_effect_spearman"], 0.912185029960394), "Discovery-validation fidelity changed")

    decomposition = pd.read_csv(
        ROOT / "results/state_shared_revision_v2/fine_state_models/programme_composition_decomposition.tsv",
        sep="\t",
    )
    primary = decomposition.loc[
        decomposition["partition"].eq("validation")
        & decomposition["k"].eq(4)
        & decomposition["scope"].eq("all")
        & decomposition["decomposition_type"].eq("fine_state")
    ].set_index("component")
    total = primary.loc["total", "estimate"]
    composition = primary.loc["composition", "estimate"]
    within = primary.loc["within", "estimate"]
    require(close(total, composition + within), "Primary decomposition does not close")
    require(close(within / total, 0.7911514513877733, 1e-9), "Primary within-state fraction changed")

    pathway = pd.read_csv(
        ROOT / "results/state_aware_program_v1/functional_architecture_v1/pathway_replication/pathway_replication_summary.tsv",
        sep="\t",
    )
    replicated = pathway.loc[pathway["replicated"].astype(bool)].copy()
    require(len(replicated) == 58, "Expected 58 replicated pathways")
    require(replicated["discovery_common_direction"].value_counts().to_dict() == {"Down": 54, "Up": 4}, "Unexpected pathway directions")
    require(replicated["community_id"].nunique() == 28, "Expected 28 leading-edge communities")

    all_pathways = pd.read_csv(
        ROOT / "results/state_aware_program_v1/functional_architecture_v1/pathway_replication.tsv",
        sep="\t",
    )
    wnt = all_pathways.loc[all_pathways["gene_set"].eq("HALLMARK_WNT_BETA_CATENIN_SIGNALING")]
    require(len(wnt) == 8 and wnt["direction"].eq("Up").all(), "WNT direction context changed")
    require(not wnt["replicated"].astype(bool).any(), "WNT should not be labelled pathway-replicated")

    meta = pd.read_csv(
        ROOT / "results/state_shared_revision_v2/external_meta/random_effects_meta_summary.tsv",
        sep="\t",
    )
    full = meta.loc[
        meta["signature_id"].eq("state_shared_1843")
        & meta["excluded_cohort"].eq("__NONE__")
    ].iloc[0]
    require(close(full["pooled_standardized_effect"], 1.90288083576639), "External pooled effect changed")


def check_figures() -> None:
    for stem in MAIN_FIGURES + SUPP_FIGURES:
        for suffix in ("pdf", "png", "svg", "tiff"):
            path = FIGURE_DIR / f"{stem}.{suffix}"
            require(path.is_file() and path.stat().st_size > 0, f"Missing figure: {path.name}")
    source_files = list((FIGURE_DIR / "source_data").glob("figure*.tsv"))
    require(len(source_files) >= 68, "Too few panel source-data files")
    require(
        not (FIGURE_DIR / "source_data/figureS2g_dslab_common_composition.tsv").exists(),
        "Governed DSLab patient-level source is public",
    )
    for name in (
        "figure4a_pathway_replication_heatmap.tsv",
        "figure4b_running_enrichment.tsv",
        "figure4c_pathway_network_nodes.tsv",
        "figure4d_leading_edge_gene_heatmap.tsv",
    ):
        require((FIGURE_DIR / "source_data" / name).is_file(), f"Missing Figure 4 source: {name}")


def check_release_hygiene() -> None:
    required = [
        "README.md",
        "CITATION.cff",
        "environment.yml",
        "data/datasets.tsv",
        "workflow/pipeline.tsv",
        "analysis/state_aware_pathway_replication_v1.R",
        "analysis/build_functional_architecture_figure4_v1.R",
        "analysis/validate_state_shared_full_programme_extended_layers_v2.py",
    ]
    missing = [path for path in required if not (ROOT / path).is_file()]
    require(not missing, "Missing required files: " + ", ".join(missing))

    forbidden_path = re.compile(
        r"287|signature_12|genki|virtual_knockout|wgcna|legacy_audit|objective_compact|stability_consensus",
        flags=re.IGNORECASE,
    )
    forbidden = [
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if ".git" not in path.parts and forbidden_path.search(path.as_posix())
    ]
    require(not forbidden, "Retired or audit-only paths remain: " + ", ".join(forbidden))

    oversized = [
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.is_file() and ".git" not in path.parts and path.stat().st_size > MAX_PUBLIC_FILE_BYTES
    ]
    require(not oversized, "Files exceed 50 MiB: " + ", ".join(oversized))

    local_home = "/home/" + "elliottlv"
    forbidden_text = re.compile(re.escape(local_home) + r"|ghp_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}")
    hits: list[str] = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or ".git" in path.parts or path.suffix.lower() not in {
            ".md", ".py", ".r", ".json", ".txt", ".tsv", ".yml", ".yaml", ".cff"
        }:
            continue
        if forbidden_text.search(path.read_text(encoding="utf-8", errors="ignore")):
            hits.append(str(path.relative_to(ROOT)))
    require(not hits, "Machine path or credential-like text found: " + ", ".join(hits))


def main() -> None:
    checks = [
        ("frozen definitions", check_definitions),
        ("principal numerical results", check_principal_results),
        ("figures and source data", check_figures),
        ("release hygiene", check_release_hygiene),
    ]
    for label, function in checks:
        function()
        print(f"[PASS] {label}")
    print(f"Release verification passed: {len(checks)}/{len(checks)} checks")


if __name__ == "__main__":
    main()
