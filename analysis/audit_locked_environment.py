#!/usr/bin/env python3
"""Verify the software versions used by the locked manuscript release."""

from __future__ import annotations

import importlib.metadata
import platform
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_TSV = ROOT / "reports" / "locked_environment_audit_2026-07-10.tsv"
OUT_TXT = ROOT / "reports" / "locked_environment_audit_2026-07-10.txt"

EXPECTED_PYTHON = {
    "Python": "3.11.15",
    "numpy": "1.26.4",
    "pandas": "2.3.3",
    "scipy": "1.17.1",
    "statsmodels": "0.14.6",
    "h5py": "3.13.0",
    "pysam": "0.24.0",
    "openpyxl": "3.1.5",
    "python-docx": "1.2.0",
    "matplotlib": "3.11.0",
    "scikit-learn": "1.9.0",
}
EXPECTED_R = {
    "R": "4.5.3",
    "ggplot2": "4.0.3",
    "dplyr": "1.2.1",
    "tidyr": "1.3.2",
    "patchwork": "1.3.2",
    "scales": "1.4.0",
    "ggrepel": "0.9.8",
    "svglite": "2.2.2",
    "ragg": "1.5.2",
    "survival": "3.8.6",
    "limma": "3.66.0",
}


def command_output(arguments: list[str]) -> str:
    return subprocess.run(
        arguments,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout.strip()


def r_versions() -> dict[str, str]:
    packages = [name for name in EXPECTED_R if name != "R"]
    expression = (
        'cat("R\\t", as.character(getRversion()), "\\n", sep=""); '
        f'pkgs <- c({", ".join(repr(package) for package in packages)}); '
        'for (p in pkgs) cat(p, "\\t", '
        'if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) '
        'else "MISSING", "\\n", sep="")'
    )
    output = command_output(["Rscript", "-e", expression])
    return dict(line.split("\t", 1) for line in output.splitlines())


def main() -> None:
    rows: list[dict[str, object]] = []
    observed_python = {"Python": platform.python_version()}
    for package in EXPECTED_PYTHON:
        if package == "Python":
            continue
        try:
            observed_python[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            observed_python[package] = "MISSING"
    for component, expected in EXPECTED_PYTHON.items():
        observed = observed_python[component]
        rows.append(
            {
                "layer": "Python",
                "component": component,
                "expected_version": expected,
                "observed_version": observed,
                "pass": observed == expected,
            }
        )

    observed_r = r_versions()
    for component, expected in EXPECTED_R.items():
        observed = observed_r.get(component, "MISSING")
        rows.append(
            {
                "layer": "R",
                "component": component,
                "expected_version": expected,
                "observed_version": observed,
                "pass": observed == expected,
            }
        )

    pandoc_line = command_output(["pandoc", "--version"]).splitlines()[0]
    pandoc_version = pandoc_line.split(maxsplit=1)[1]
    rows.append(
        {
            "layer": "document",
            "component": "pandoc",
            "expected_version": "3.10",
            "observed_version": pandoc_version,
            "pass": pandoc_version == "3.10",
        }
    )

    audit = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(OUT_TSV, sep="\t", index=False)
    lines = [
        "Locked computational environment audit",
        "======================================",
        f"Platform: {platform.platform()}",
        f"Checks: {len(audit)}",
        f"Passed: {int(audit['pass'].sum())}",
        f"Failed: {int((~audit['pass']).sum())}",
        "",
    ]
    lines.extend(
        f"[{'PASS' if row.pass_ else 'FAIL'}] {row.layer}/{row.component}: "
        f"expected {row.expected_version}; observed {row.observed_version}"
        for row in audit.rename(columns={"pass": "pass_"}).itertuples(index=False)
    )
    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if not audit["pass"].all():
        failed = audit.loc[~audit["pass"], "component"].tolist()
        raise SystemExit("Environment audit failed: " + ", ".join(failed))
    print(f"Locked environment audit passed: {len(audit)}/{len(audit)}")


if __name__ == "__main__":
    main()
