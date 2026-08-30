#!/usr/bin/env python3
"""Fast integrity checks for the public identity-remodelling release."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
SIGNATURE_DIR = ROOT / "data" / "signatures"
FIGURE_DIR = ROOT / "figures" / "communications_biology_v2.0"
MAX_PUBLIC_FILE_BYTES = 50 * 1024 * 1024

MAIN_FIGURES = [
    "figure1_study_design_and_programme_derivation",
    "figure2_donor_disjoint_identity_remodelling",
    "figure3_external_recurrence_and_archival_transfer",
    "figure4_regulatory_and_epithelial_context",
    "figure5_genetic_perturbation_support",
    "figure6_reduced_measurement_candidate",
]
SUPP_FIGURES = [
    "figureS1_sampling_and_donor_stability",
    "figureS2_historical_gene_set_audit",
    "figureS3_fine_state_and_composition_sensitivities",
    "figureS4_functional_and_regulatory_structure",
    "figureS5_external_and_ffpe_sensitivities",
    "figureS6_compact_derivation_and_benchmarks",
    "figureS7_multiomic_atlas_spatial_protein",
    "figureS8_perturbation_and_virtual_context",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def close(observed: float, expected: float, tolerance: float = 1e-10) -> bool:
    return abs(observed - expected) <= tolerance


def check_definitions() -> None:
    ranking = pd.read_csv(SIGNATURE_DIR / "common_effect_ranking_8221.tsv.gz", sep="\t")
    strict = pd.read_csv(SIGNATURE_DIR / "state_aware_high_confidence_1843.tsv", sep="\t")
    candidates = pd.read_csv(SIGNATURE_DIR / "portable_candidates_53.tsv", sep="\t")
    compact = pd.read_csv(SIGNATURE_DIR / "compact_candidate_8_genes.tsv", sep="\t")
    require(len(ranking) == 8221, "Expected an 8,221-gene common-effect ranking")
    require(len(strict) == 1843, "Expected 1,843 high-confidence genes")
    require(strict["shared_direction"].value_counts().to_dict() == {"down": 959, "up": 884}, "Unexpected programme arms")
    require(len(candidates) == 53, "Expected 53 portable candidates")
    require(len(compact) == 8, "Expected eight compact genes")
    require(compact["arm"].value_counts().to_dict() == {"up": 4, "down": 4}, "Compact candidate is not balanced")
    require(set(compact["gene"]) == {"EPHB2", "REG1A", "LTBP1", "RNF43", "CALM2", "COX6C", "B2M", "ACAA2"}, "Compact membership changed")
    require(set(compact["gene"]).issubset(set(candidates["gene"])), "Compact genes are not portable candidates")
    require(not compact["validation_outcomes_used"].astype(bool).any(), "Validation leakage flag changed")
    require(compact["panel_frozen_before_validation"].astype(bool).all(), "Freeze flag changed")


def check_principal_results() -> None:
    donor = pd.read_csv(ROOT / "results/state_shared_revision_v2/donor_site/donor_disjoint_replication_summary.tsv", sep="\t").iloc[0]
    require(int(donor["n_testable_strict_genes"]) == 1646, "Donor-disjoint testable count changed")
    require(close(donor["common_direction_match_fraction"], 0.986026731470231), "Direction agreement changed")
    require(close(donor["discovery_validation_effect_spearman"], 0.912185029960394), "Discovery-validation fidelity changed")

    decomposition = pd.read_csv(ROOT / "results/state_shared_revision_v2/fine_state_models/programme_composition_decomposition.tsv", sep="\t")
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

    meta = pd.read_csv(ROOT / "results/state_shared_revision_v2/external_meta/random_effects_meta_summary.tsv", sep="\t")
    full = meta.loc[meta["signature_id"].eq("state_shared_1843") & meta["excluded_cohort"].eq("__NONE__")].iloc[0]
    require(close(full["pooled_standardized_effect"], 1.90288083576639), "External pooled effect changed")

    internal = pd.read_csv(ROOT / "results/state_shared_revision_v2/compact_rank/random_eight_gene_benchmark_summary.tsv", sep="\t").iloc[0]
    external = pd.read_csv(ROOT / "results/state_shared_revision_v2/external_rank/random_eight_gene_benchmark_summary.tsv", sep="\t").iloc[0]
    require(internal["observed_compact_spearman"] > internal["random_q95"], "Internal compact benchmark boundary changed")
    require(external["observed_median_cohort_spearman"] < external["random_q95"], "External non-uniqueness boundary changed")


def check_figures_and_workbooks() -> None:
    for stem in MAIN_FIGURES + SUPP_FIGURES:
        for suffix in ("pdf", "png", "svg"):
            require((FIGURE_DIR / f"{stem}.{suffix}").is_file(), f"Missing figure: {stem}.{suffix}")
    source = load_workbook(ROOT / "data/source_data/Source_Data.xlsx", read_only=True)
    supplement = load_workbook(ROOT / "data/supplementary/Supplementary_Tables_1-12.xlsx", read_only=True)
    require(len(source.sheetnames) == 86, "Unexpected Source Data worksheet count")
    require(len(supplement.sheetnames) == 72, "Unexpected supplementary worksheet count")
    require(not any("dslab_common" in name.lower() for name in source.sheetnames), "Governed DSLab patient-level source was included")
    require(all(any(name.startswith(f"T{number}") for name in supplement.sheetnames) for number in range(1, 13)), "A supplementary table group is missing")


def check_required_files_and_hygiene() -> None:
    required = [
        "README.md",
        "CITATION.cff",
        "environment.yml",
        "data/datasets.tsv",
        "workflow/pipeline.tsv",
        "analysis/state_aware_build_discovery_pseudobulk_v1.R",
        "analysis/state_aware_integrate_common_effects_v1.R",
        "analysis/state_shared_revision_define_fine_states_v2.py",
        "analysis/state_shared_revision_fine_state_models_v2.R",
        "analysis/derive_state_shared_compact_panel_v1.R",
        "analysis/build_state_shared_revision_figure6_v3.R",
    ]
    missing = [path for path in required if not (ROOT / path).is_file()]
    require(not missing, "Missing required files: " + ", ".join(missing))
    oversized = [str(path.relative_to(ROOT)) for path in ROOT.rglob("*") if path.is_file() and ".git" not in path.parts and path.stat().st_size > MAX_PUBLIC_FILE_BYTES]
    require(not oversized, "Files exceed 50 MiB: " + ", ".join(oversized))
    local_home = "/home/" + "elliottlv"
    forbidden_text = re.compile(re.escape(local_home) + r"|ghp_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}")
    hits: list[str] = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or ".git" in path.parts or path.suffix.lower() not in {".md", ".py", ".r", ".json", ".txt", ".tsv", ".yml", ".yaml", ".cff"}:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if forbidden_text.search(text):
            hits.append(str(path.relative_to(ROOT)))
    require(not hits, "Machine path or credential-like text found: " + ", ".join(hits))
    require(not (ROOT / "data/source_data/figureS3g_dslab_common_composition.tsv").exists(), "Governed DSLab source file is public")


def main() -> None:
    checks = [
        ("frozen definitions", check_definitions),
        ("principal numerical results", check_principal_results),
        ("figures and workbooks", check_figures_and_workbooks),
        ("required files and hygiene", check_required_files_and_hygiene),
    ]
    for label, function in checks:
        function()
        print(f"[PASS] {label}")
    print(f"Release verification passed: {len(checks)}/{len(checks)} checks")


if __name__ == "__main__":
    main()
