#!/usr/bin/env python3
"""Deterministically re-run GenKI and compare every reported result table.

Compressed TSV files are compared after decompression so gzip container
timestamps cannot create a false mismatch. The analysis manifest is compared
after removing wall-clock duration, which is not a scientific result.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "virtual_knockout_validation_v2_9"
RUN_SCRIPT = ROOT / "analysis" / "run_genki_virtual_knockout_validation_v2_9.py"
LOG_PATH = RESULT_DIR / "genki_reproducibility_rerun.log"
AUDIT_PATH = RESULT_DIR / "genki_reproducibility_audit.tsv"
SUMMARY_PATH = RESULT_DIR / "genki_reproducibility_summary.json"

RESULT_FILES = (
    "genki_training_metrics.tsv",
    "genki_gene_impact_by_seed.tsv.gz",
    "genki_gene_impact_consensus.tsv.gz",
    "genki_seed_stability.tsv",
    "genki_distance_metric_sensitivity.tsv",
    "genki_fixed_gene_set_tests.tsv",
    "genki_aggregate_validation_endpoints.tsv",
    "genki_panel_knockout_impact_matrix.tsv",
    "empirical_direction_calibration.tsv",
    "prespecified_knockout_targets.tsv",
    "genki_network_gene_metrics.tsv",
    "genki_validation_qa.tsv",
    "genki_analysis_manifest.json",
)


def normalized_bytes(path: Path) -> bytes:
    if path.name.endswith(".tsv.gz"):
        with gzip.open(path, "rb") as handle:
            return handle.read()
    if path.name == "genki_analysis_manifest.json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload.pop("elapsed_minutes", None)
        return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return path.read_bytes()


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    missing = [name for name in RESULT_FILES if not (RESULT_DIR / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing original result files: {missing}")

    before = {name: normalized_bytes(RESULT_DIR / name) for name in RESULT_FILES}
    start = time.time()
    clean_env = os.environ.copy()
    for key in (
        "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY", "http_proxy", "https_proxy", "all_proxy"
    ):
        clean_env.pop(key, None)

    with LOG_PATH.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            [sys.executable, str(RUN_SCRIPT)],
            cwd=ROOT,
            env=clean_env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    log_text = LOG_PATH.read_text(encoding="utf-8")
    log_text = log_text.replace(str(ROOT), ".")
    env_root = str(Path(sys.executable).resolve().parents[1])
    log_text = log_text.replace(env_root, "<CONDA_ENV>")
    LOG_PATH.write_text(log_text, encoding="utf-8")
    elapsed_seconds = time.time() - start

    rows: list[dict[str, object]] = []
    for name in RESULT_FILES:
        path = RESULT_DIR / name
        after = normalized_bytes(path) if path.is_file() else b""
        rows.append(
            {
                "file": name,
                "normalization": (
                    "decompressed_content"
                    if name.endswith(".tsv.gz")
                    else "json_without_elapsed_minutes"
                    if name.endswith(".json")
                    else "raw_bytes"
                ),
                "before_sha256": digest(before[name]),
                "after_sha256": digest(after),
                "before_bytes": len(before[name]),
                "after_bytes": len(after),
                "exact_match": before[name] == after,
            }
        )

    audit = pd.DataFrame(rows)
    audit.to_csv(AUDIT_PATH, sep="\t", index=False)
    all_match = bool(completed.returncode == 0 and audit["exact_match"].all())
    summary = {
        "command": "python analysis/run_genki_virtual_knockout_validation_v2_9.py",
        "exit_code": completed.returncode,
        "elapsed_seconds": elapsed_seconds,
        "files_compared": len(audit),
        "files_exactly_matched": int(audit["exact_match"].sum()),
        "all_reported_outputs_exactly_reproduced": all_match,
        "verdict": "REPRODUCIBLE" if all_match else "NOT_REPRODUCIBLE",
        "notes": (
            "All reported tabular outputs and the timing-normalized manifest were compared. "
            "Serialized model binaries were not used as scientific endpoints and were not compared."
        ),
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if not all_match:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
