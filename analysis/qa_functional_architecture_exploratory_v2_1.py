#!/usr/bin/env python3
"""Audit corrected exploratory WGCNA integration outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RUN = (
    ROOT
    / "results"
    / "state_aware_program_v1"
    / "functional_architecture_exploratory_v2_1"
)
REPLAY = RUN.parent / "functional_architecture_exploratory_v2_1_replay"
FIGURE = (
    ROOT
    / "figures"
    / "communications_biology_v3.0"
    / "functional_architecture_wgcna_candidate_v2_1"
)
EXPECTED_MODULES = {"M02", "M03", "M04", "M05", "M06", "M09", "M10"}
EXPECTED_COHORTS = {"GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return values.astype(str).str.lower().eq("true")


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: object) -> None:
        checks.append({"check": name, "pass": bool(passed), "detail": detail})

    validation = pd.read_csv(RUN / "module_validation.tsv", sep="\t")
    selected = set(validation.loc[as_bool(validation["analysis_route_pass"]), "module"])
    check("selected_modules", selected == EXPECTED_MODULES, sorted(selected))
    check(
        "validation_columns_unique",
        not validation.columns.duplicated().any(),
        validation.columns[validation.columns.duplicated()].tolist(),
    )
    check(
        "no_external_gate_pass",
        int(as_bool(validation["external_gate_pass"]).sum()) == 0,
        int(as_bool(validation["external_gate_pass"]).sum()),
    )
    selected_validation = validation.loc[validation["module"].isin(EXPECTED_MODULES)]
    check(
        "all_selected_technical_gate_pass",
        as_bool(selected_validation["technical_gate_pass"]).all(),
        int(as_bool(selected_validation["technical_gate_pass"]).sum()),
    )

    effects = pd.read_csv(
        RUN / "module_external_validation/module_external_cohort_effects.tsv", sep="\t"
    )
    cohorts = set(effects["cohort"])
    check("five_external_cohorts", cohorts == EXPECTED_COHORTS, sorted(cohorts))
    coverage = pd.read_csv(
        RUN / "module_external_validation/module_platform_coverage.tsv", sep="\t"
    )
    gse8671 = coverage.loc[
        coverage["cohort"].eq("GSE8671") & coverage["module"].isin(EXPECTED_MODULES)
    ]
    check(
        "gse8671_full_annotation_coverage",
        len(gse8671) == 7 and gse8671["coverage_fraction"].min() >= 0.80,
        float(gse8671["coverage_fraction"].min()),
    )

    perturbation_gate = pd.read_csv(
        RUN / "module_perturbation_protein/module_perturbation_gate.tsv", sep="\t"
    )
    check(
        "all_selected_modules_projected_to_perturbations",
        set(perturbation_gate["module"]) == EXPECTED_MODULES,
        sorted(perturbation_gate["module"].tolist()),
    )
    priorities = pd.read_csv(RUN / "protein_priorities.tsv", sep="\t")
    check("sentinel_count", len(priorities) == 16, len(priorities))
    check(
        "sentinels_exploratory_only",
        set(priorities["final_status"].dropna()) == {"exploratory_only"},
        sorted(set(priorities["final_status"].dropna())),
    )
    check(
        "no_regulatory_node",
        not priorities["priority_role"].eq("regulatory_node").any(),
        int(priorities["priority_role"].eq("regulatory_node").sum()),
    )

    replay = pd.read_csv(REPLAY / "replay_hash_audit.tsv", sep="\t")
    check(
        "external_replay_identical",
        as_bool(replay["identical"]).all(),
        int(as_bool(replay["identical"]).sum()),
    )

    figure_stems = [
        "figure4_functional_architecture_wgcna_candidate_v2_1",
        "figureS_wgcna_structure_and_context_candidate_v2_1",
    ]
    figure_files = [FIGURE / f"{stem}.{extension}" for stem in figure_stems for extension in ["pdf", "svg", "png", "tiff"]]
    check(
        "figure_assets_complete",
        all(path.exists() and path.stat().st_size > 10_000 for path in figure_files),
        {path.name: path.stat().st_size if path.exists() else 0 for path in figure_files},
    )

    manifest_path = RUN / "integration_summary/exploratory_wgcna_integration_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    output_hashes = manifest["output_sha256"]
    integration_outputs = {
        "integrated_summary": RUN / "integration_summary/module_integrated_evidence.tsv",
        "perturbation_correlations": RUN / "integration_summary/module_perturbation_profile_correlations.tsv",
        "figure_pdf": (
            ROOT
            / "figures/communications_biology_v3.0/functional_architecture_exploratory_v2_1/exploratory_wgcna_module_integration_v2_1.pdf"
        ),
        "figure_svg": (
            ROOT
            / "figures/communications_biology_v3.0/functional_architecture_exploratory_v2_1/exploratory_wgcna_module_integration_v2_1.svg"
        ),
    }
    observed_hashes = {key: sha256(path) for key, path in integration_outputs.items()}
    check(
        "integration_manifest_hashes",
        observed_hashes == output_hashes,
        {"observed": observed_hashes, "manifest": output_hashes},
    )

    required_text = [
        RUN / "integration_summary/exploratory_wgcna_integration_report_v2_1.md",
        ROOT
        / "manuscripts/communications_biology_2026-08-12/wgcna_manuscript_insert_candidate_v2_1_2026-08-31.md",
    ]
    check(
        "reports_present",
        all(path.exists() and path.stat().st_size > 2_000 for path in required_text),
        {path.name: path.stat().st_size if path.exists() else 0 for path in required_text},
    )

    passed = all(item["pass"] for item in checks)
    output = {
        "analysis": "qa_functional_architecture_exploratory_v2_1",
        "all_checks_pass": passed,
        "n_checks": len(checks),
        "checks": checks,
    }
    (RUN / "functional_architecture_exploratory_v2_1_QA.json").write_text(
        json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if not passed:
        failed = [item["check"] for item in checks if not item["pass"]]
        raise RuntimeError(f"QA failed: {failed}")
    print(f"Exploratory WGCNA v2.1 QA passed: {len(checks)} checks")


if __name__ == "__main__":
    main()
