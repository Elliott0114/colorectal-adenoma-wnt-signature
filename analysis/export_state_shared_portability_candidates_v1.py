#!/usr/bin/env python3
"""Audit label-blind platform availability for the state-shared programme.

This script performs no statistical modelling and reads no validation outcomes.
It reuses the project's checked platform parsers solely to determine whether a
gene is represented as a feature on each intended transfer platform.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "analysis"
OUT_DIR = ROOT / "results" / "state_aware_program_v1" / "panel_derivation"
sys.path.insert(0, str(ANALYSIS_DIR))

import external_sporadic_adenoma_validation as external  # noqa: E402
import gse117606_paired_route_validation as ffpe  # noqa: E402


COMMON_PATH = (
    ROOT
    / "results"
    / "state_aware_program_v1"
    / "common_effects"
    / "cross_state_common_effects.tsv.gz"
)
EXTERNAL_COHORTS = ["GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820"]
ALL_PLATFORMS = [*EXTERNAL_COHORTS, "GSE117606"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ncbi_annotation(genes: list[str]) -> pd.DataFrame:
    info = pd.read_csv(
        ffpe.GENE_INFO_PATH,
        sep="\t",
        dtype=str,
        low_memory=False,
    )
    direct = info.loc[
        info["Symbol"].isin(genes),
        ["GeneID", "Symbol", "type_of_gene", "Nomenclature_status"],
    ].copy()
    direct = direct.rename(columns={"Symbol": "gene"})
    direct["ncbi_current_symbol"] = direct["gene"]
    direct["annotation_method"] = "current_ncbi_symbol"

    missing = sorted(set(genes) - set(direct["gene"]))
    alias_rows: list[dict[str, object]] = []
    if missing:
        synonym_sets = info["Synonyms"].fillna("").map(lambda value: set(value.split("|")))
        for gene in missing:
            matches = info.loc[
                synonym_sets.map(lambda values: gene in values),
                ["GeneID", "Symbol", "type_of_gene", "Nomenclature_status"],
            ]
            if len(matches) == 1:
                row = matches.iloc[0]
                alias_rows.append(
                    {
                        "GeneID": row["GeneID"],
                        "gene": gene,
                        "type_of_gene": row["type_of_gene"],
                        "Nomenclature_status": row["Nomenclature_status"],
                        "ncbi_current_symbol": row["Symbol"],
                        "annotation_method": "unique_ncbi_synonym",
                    }
                )
    annotation = pd.concat([direct, pd.DataFrame(alias_rows)], ignore_index=True)
    annotation = annotation.sort_values(
        ["gene", "annotation_method", "GeneID"], kind="mergesort"
    ).drop_duplicates("gene")
    return pd.DataFrame({"gene": genes}).merge(annotation, on="gene", how="left")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    common = pd.read_csv(COMMON_PATH, sep="\t")
    core = common.loc[
        common["strict_state_shared"],
        [
            "gene",
            "shared_direction",
            "common_effect",
            "common_z",
            "common_q_value",
            "max_lfsr",
        ],
    ].copy()
    core = core.rename(columns={"shared_direction": "arm"})
    if set(core["arm"]) != {"up", "down"} or core["gene"].duplicated().any():
        raise RuntimeError("Strict state-shared programme is malformed")
    genes = core["gene"].astype(str).tolist()
    wanted = set(genes)

    mappings = {
        "GSE8671": external.gpl570_mapping(),
        "GSE50114": external.gpl6480_mapping(),
        "GSE41657": external.gpl6480_mapping(),
        "GSE40362": external.gpl8432_mapping(),
    }
    presence_rows: list[dict[str, object]] = []
    for accession in ["GSE8671", "GSE50114", "GSE41657", "GSE40362"]:
        _, expression, _ = external.parse_series_matrix(accession)
        if accession == "GSE50114":
            expression = external.gse50114_raw_expression()
        present_features = set(expression.columns.astype(str))
        mapping = mappings[accession]
        measured_genes = set(
            mapping.loc[mapping["feature_id"].isin(present_features), "gene"]
        )
        presence_rows.extend(
            {
                "gene": gene,
                "platform": accession,
                "feature_present": gene in measured_genes,
                "eligibility_information": "feature availability only",
            }
            for gene in genes
        )
        del expression

    _, direct_expression, _ = external.direct_rnaseq_expression(wanted)
    measured_genes = set(direct_expression.columns.astype(str))
    presence_rows.extend(
        {
            "gene": gene,
            "platform": "GSE72820",
            "feature_present": gene in measured_genes,
            "eligibility_information": "feature availability only",
        }
        for gene in genes
    )
    del direct_expression

    ffpe_expression = ffpe.parse_expression()
    ffpe_mapping = ffpe.map_symbols_to_features(wanted, ffpe_expression)
    ffpe_present = set(ffpe_mapping.loc[ffpe_mapping["feature_present"], "gene"])
    presence_rows.extend(
        {
            "gene": gene,
            "platform": "GSE117606",
            "feature_present": gene in ffpe_present,
            "eligibility_information": "feature availability only",
        }
        for gene in genes
    )

    presence = pd.DataFrame(presence_rows)
    wide = presence.pivot(
        index="gene", columns="platform", values="feature_present"
    ).reset_index()
    annotation = ncbi_annotation(genes)
    audit = (
        core.merge(annotation, on="gene", how="left", validate="one_to_one")
        .merge(wide, on="gene", how="left", validate="one_to_one")
    )
    audit["n_platforms_present"] = audit[ALL_PLATFORMS].sum(axis=1).astype(int)
    audit["all_six_platforms_present"] = audit[ALL_PLATFORMS].all(axis=1)
    audit["protein_coding"] = audit["type_of_gene"].eq("protein-coding")
    audit["objective_selection_eligible"] = (
        audit["all_six_platforms_present"] & audit["protein_coding"]
    )
    audit["eligibility_rule"] = (
        "strict state-shared core; NCBI protein-coding; feature present on all "
        "five external platforms and GSE117606 FFPE"
    )
    candidates = audit.loc[audit["objective_selection_eligible"]].copy()
    if candidates.groupby("arm")["gene"].nunique().min() < 10:
        raise RuntimeError("Fewer than ten portable candidates remain in one arm")

    presence.to_csv(OUT_DIR / "gene_platform_presence_long.tsv", sep="\t", index=False)
    audit.to_csv(OUT_DIR / "state_shared_gene_portability_audit.tsv", sep="\t", index=False)
    candidates.to_csv(
        OUT_DIR / "portable_state_shared_candidate_universe.tsv",
        sep="\t",
        index=False,
    )
    summary = (
        audit.groupby("arm", observed=True)
        .agg(
            strict_state_shared_genes=("gene", "nunique"),
            all_six_platforms_present=("all_six_platforms_present", "sum"),
            portable_protein_coding=("objective_selection_eligible", "sum"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_DIR / "candidate_universe_summary.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "label-blind portability audit for the strict state-shared programme",
        "created_utc": pd.Timestamp.utcnow().isoformat(),
        "input_path": str(COMMON_PATH.resolve()),
        "input_sha256": sha256(COMMON_PATH),
        "platforms": ALL_PLATFORMS,
        "gene_annotation": "NCBI Homo_sapiens.gene_info",
        "outcome_labels_used": False,
        "validation_effects_used": False,
        "eligibility_rule": (
            "strict state-shared core; NCBI protein-coding; feature present on all "
            "five external platforms and GSE117606 FFPE"
        ),
        "summary": summary.to_dict(orient="records"),
    }
    (OUT_DIR / "portability_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        "Portable state-shared candidates: "
        + ", ".join(
            f"{row.arm}={int(row.portable_protein_coding)}"
            for row in summary.itertuples(index=False)
        )
    )


if __name__ == "__main__":
    main()
