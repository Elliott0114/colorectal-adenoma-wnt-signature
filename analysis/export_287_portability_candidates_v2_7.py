#!/usr/bin/env python3
"""Create a label-blind portability audit for the 287-gene discovery core.

Only platform feature availability and NCBI gene type are inspected.  No tissue
labels, score effects, validation outcomes or perturbation results enter this
eligibility step.  The resulting candidate universe is consumed by the R-only
donor-grouped sparse-panel selection workflow.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "analysis"
OUT_DIR = ROOT / "results" / "objective_compact_panel_v2_7"
sys.path.insert(0, str(ANALYSIS_DIR))

import gse117606_paired_route_validation as ffpe  # noqa: E402
import translation_reduced_panel_v2_0 as reduced  # noqa: E402


EXTERNAL_COHORTS = ["GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820"]
ALL_PLATFORMS = [*EXTERNAL_COHORTS, "GSE117606"]


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
    core = pd.read_csv(
        ROOT
        / "results"
        / "data_adaptive_panel_pilot_v2_6"
        / "stable_error_controlled_core.tsv",
        sep="\t",
    )
    if len(core) != 287 or set(core["arm"]) != {"up", "down"}:
        raise RuntimeError("Expected the fixed 287-gene two-arm discovery core")
    genes = core["gene"].astype(str).tolist()

    external_data = reduced.external_gene_data(pd.DataFrame({"gene": genes}))
    presence_rows: list[dict[str, object]] = []
    for cohort in EXTERNAL_COHORTS:
        expression = external_data[cohort][1]
        measured = set(expression.columns)
        presence_rows.extend(
            {
                "gene": gene,
                "platform": cohort,
                "feature_present": gene in measured,
                "eligibility_information": "feature availability only; labels not accessed",
            }
            for gene in genes
        )

    ffpe_expression = ffpe.parse_expression()
    ffpe_mapping = ffpe.map_symbols_to_features(set(genes), ffpe_expression)
    ffpe_present = set(ffpe_mapping.loc[ffpe_mapping["feature_present"], "gene"])
    presence_rows.extend(
        {
            "gene": gene,
            "platform": "GSE117606",
            "feature_present": gene in ffpe_present,
            "eligibility_information": "feature availability only; labels not accessed",
        }
        for gene in genes
    )
    presence = pd.DataFrame(presence_rows)

    wide = presence.pivot(index="gene", columns="platform", values="feature_present").reset_index()
    annotation = ncbi_annotation(genes)
    audit = (
        core[
            [
                "gene",
                "arm",
                "logFC",
                "adj.P.Val",
                "direction_stability",
                "bootstrap_ci_low",
                "bootstrap_ci_high",
            ]
        ]
        .merge(annotation, on="gene", how="left", validate="one_to_one")
        .merge(wide, on="gene", how="left", validate="one_to_one")
    )
    audit["n_platforms_present"] = audit[ALL_PLATFORMS].sum(axis=1).astype(int)
    audit["all_six_platforms_present"] = audit[ALL_PLATFORMS].all(axis=1)
    audit["protein_coding"] = audit["type_of_gene"].eq("protein-coding")
    audit["objective_selection_eligible"] = (
        audit["all_six_platforms_present"] & audit["protein_coding"]
    )
    audit["eligibility_rule"] = (
        "287-core member; protein-coding; feature present on all five external platforms and FFPE"
    )

    candidates = audit.loc[audit["objective_selection_eligible"]].copy()
    if candidates.groupby("arm")["gene"].nunique().min() < 10:
        raise RuntimeError("Fewer than 10 portable protein-coding genes remain in one arm")

    presence.to_csv(OUT_DIR / "gene_platform_presence_long.tsv", sep="\t", index=False)
    audit.to_csv(OUT_DIR / "gene_portability_audit.tsv", sep="\t", index=False)
    candidates.to_csv(
        OUT_DIR / "portable_protein_coding_candidate_universe.tsv",
        sep="\t",
        index=False,
    )
    summary = (
        audit.groupby("arm", observed=True)
        .agg(
            core_genes=("gene", "nunique"),
            portable_all_six=("all_six_platforms_present", "sum"),
            portable_protein_coding=("objective_selection_eligible", "sum"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_DIR / "candidate_universe_summary.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "label-blind portability audit for the 287-gene core",
        "platforms": ALL_PLATFORMS,
        "gene_annotation": "NCBI Homo_sapiens.gene_info",
        "outcome_labels_used": False,
        "expression_effects_used": False,
        "eligibility_rule": (
            "member of the 287-gene stable core; NCBI protein-coding; feature present "
            "on all five external platforms and GSE117606 FFPE"
        ),
    }
    (OUT_DIR / "portability_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
