#!/usr/bin/env python3
"""Run fast integrity checks for the public manuscript repository."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SIGNATURE_DIR = ROOT / "data" / "signatures"
EXPECTED_SIGNATURE_SHA256 = (
    "91e564f42bdc5fb0188605b76116bcfb126e5102d62a33d6890fa67018928380"
)
MAX_PUBLIC_FILE_BYTES = 50 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_signatures() -> None:
    core = pd.read_csv(SIGNATURE_DIR / "core_287_genes.tsv", sep="\t")
    portable = pd.read_csv(
        SIGNATURE_DIR / "portable_candidates_62_genes.tsv", sep="\t"
    )
    signature_path = SIGNATURE_DIR / "signature_12_genes.tsv"
    signature = pd.read_csv(signature_path, sep="\t")

    require(len(core) == 287, f"Expected 287 core genes, observed {len(core)}")
    require(len(portable) == 62, f"Expected 62 portable genes, observed {len(portable)}")
    require(len(signature) == 12, f"Expected 12 signature genes, observed {len(signature)}")
    require(core["gene"].is_unique, "Core contains duplicate symbols")
    require(portable["gene"].is_unique, "Portable set contains duplicate symbols")
    require(signature["gene"].is_unique, "Signature contains duplicate symbols")
    require(set(signature["arm"]) == {"up", "down"}, "Signature is not two-armed")
    require(
        signature["arm"].value_counts().to_dict() == {"up": 6, "down": 6},
        "Signature is not balanced 6 up / 6 down",
    )
    require(
        set(signature["gene"]).issubset(set(portable["gene"])),
        "Signature is not a subset of the portable candidates",
    )
    require(
        set(portable["gene"]).issubset(set(core["gene"])),
        "Portable candidates are not a subset of the 287-gene core",
    )
    require(
        not signature["fixed_gene_count_used"].astype(bool).any(),
        "Signature metadata indicates a fixed count",
    )
    require(
        not signature["validation_outcomes_used"].astype(bool).any(),
        "Signature metadata indicates validation leakage",
    )
    require(
        sha256(signature_path) == EXPECTED_SIGNATURE_SHA256,
        "Frozen signature checksum changed",
    )


def check_required_files() -> None:
    required = [
        "README.md",
        "LICENSE",
        "CITATION.cff",
        "environment.yml",
        "environment-virtual-knockout.yml",
        "data/datasets.tsv",
        "data/supplementary/Supplementary_Tables_1-12.xlsx",
        "workflow/pipeline.tsv",
        "analysis/plot_jtm_submission_figures_v2_8.R",
        "results/virtual_knockout_validation_v2_9/genki_analysis_manifest.json",
        "results/virtual_knockout_validation_v2_9/genki_reproducibility_summary.json",
    ]
    required.extend(
        f"figures/manuscript/figure{number}_{stem}.svg"
        for number, stem in [
            (1, "discovery_core_and_objective_reduction"),
            (2, "independent_replication_and_ffpe"),
            (3, "rna_atac_regulatory_support"),
            (4, "crc_atlas_cross_sectional_recurrence"),
            (5, "empirical_and_virtual_perturbation_support"),
            (6, "spatial_and_protein_context"),
        ]
    )
    missing = [path for path in required if not (ROOT / path).is_file()]
    require(not missing, "Missing required files: " + ", ".join(missing))


def check_virtual_deletion_manifest() -> None:
    result_dir = ROOT / "results" / "virtual_knockout_validation_v2_9"
    manifest = json.loads((result_dir / "genki_analysis_manifest.json").read_text(encoding="utf-8"))
    audit = json.loads(
        (result_dir / "genki_reproducibility_summary.json").read_text(encoding="utf-8")
    )
    require(manifest["analysis_role"] == "validation_not_discovery", "GenKI role changed")
    require(manifest["n_cells"] == 1664, "Unexpected GenKI cell count")
    require(manifest["n_donors"] == 13, "Unexpected GenKI donor count")
    require(manifest["model_seeds"] == [20260810, 20260811], "Unexpected GenKI seeds")
    require(manifest["matched_null_replicates"] == 10000, "Unexpected null count")
    require(manifest["qa_all_passed"] is True, "GenKI QA did not pass")
    require(
        audit["all_reported_outputs_exactly_reproduced"] is True,
        "GenKI exact-output audit did not pass",
    )


def check_repository_hygiene() -> None:
    oversized = []
    forbidden_names = []
    sensitive_hits = []
    text_suffixes = {".md", ".txt", ".tsv", ".csv", ".py", ".r", ".yml", ".yaml", ".json", ".cff"}
    secret_patterns = [
        re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
        re.compile(r"ghp_[A-Za-z0-9]{30,}"),
        re.compile(r"github_pat_[A-Za-z0-9_]{30,}"),
        re.compile(r"AKIA[0-9A-Z]{16}"),
    ]

    for path in ROOT.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        relative = path.relative_to(ROOT)
        if path.stat().st_size > MAX_PUBLIC_FILE_BYTES:
            oversized.append(f"{relative} ({path.stat().st_size} bytes)")
        if path.name in {".env", "id_rsa", "id_ed25519"}:
            forbidden_names.append(str(relative))
        if path.suffix.lower() in text_suffixes and path.stat().st_size <= 5 * 1024 * 1024:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if any(pattern.search(text) for pattern in secret_patterns):
                sensitive_hits.append(str(relative))

    require(not oversized, "Files exceed 50 MiB: " + ", ".join(oversized))
    require(not forbidden_names, "Credential-like files found: " + ", ".join(forbidden_names))
    require(not sensitive_hits, "Sensitive-token patterns found: " + ", ".join(sensitive_hits))
    require(not (ROOT / "data_sources").exists(), "Raw data_sources directory must not be public")


def main() -> None:
    checks = [
        ("required files", check_required_files),
        ("frozen signatures", check_signatures),
        ("virtual-deletion manifest", check_virtual_deletion_manifest),
        ("repository hygiene", check_repository_hygiene),
    ]
    for label, function in checks:
        function()
        print(f"[PASS] {label}")
    print(f"Release verification passed: {len(checks)}/{len(checks)} checks")


if __name__ == "__main__":
    main()

