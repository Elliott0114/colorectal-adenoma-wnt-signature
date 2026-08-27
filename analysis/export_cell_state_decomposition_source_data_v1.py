#!/usr/bin/env python3
"""Export exact panel-level source tables for the cell-state decomposition figure."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "cell_state_decomposition_v1"
OUT_DIR = RESULT_DIR / "source_data"


def read(name: str, compressed: bool = False) -> pd.DataFrame:
    suffix = ".tsv.gz" if compressed else ".tsv"
    return pd.read_csv(RESULT_DIR / f"{name}{suffix}", sep="\t")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    abundance = read("differential_abundance_propeller_style")
    panel_a = abundance.loc[
        abundance["dataset"].eq("validation") & abundance["transformation"].eq("arcsine_sqrt")
    ].copy()
    panel_a.to_csv(OUT_DIR / "figure3a_epithelial_proportions.tsv", sep="\t", index=False)

    within = read("within_cell_state_effects")
    panel_b = within.loc[
        within["dataset"].eq("validation")
        & within["score"].eq("core_287")
        & within["min_cells_per_specimen_state"].eq(20)
        & within["comparison"].eq("conventional_vs_normal")
    ].copy()
    panel_b.to_csv(OUT_DIR / "figure3b_within_state_effects.tsv", sep="\t", index=False)

    state = read("specimen_cell_state_scores", compressed=True)
    panel_c = state.loc[
        state["dataset"].eq("validation")
        & state["cell_type"].isin(["ABS", "GOB", "TAC"])
        & state["n_cells"].ge(20)
        & state["route_group"].isin(["normal", "conventional_adenoma"])
    ].copy()
    panel_c = (
        panel_c.groupby(
            ["donor_id", "route_group", "cell_type", "cell_type_label"], observed=True
        )["score__core_287__mean"]
        .mean()
        .rename("core_287_score")
        .reset_index()
    )
    paired = (
        panel_c.groupby(["donor_id", "cell_type"], observed=True)["route_group"]
        .nunique()
        .rename("n_routes")
        .reset_index()
    )
    panel_c = panel_c.merge(
        paired.loc[paired["n_routes"].eq(2), ["donor_id", "cell_type"]],
        on=["donor_id", "cell_type"],
        how="inner",
    )
    panel_c.to_csv(OUT_DIR / "figure3c_paired_donor_scores.tsv", sep="\t", index=False)

    contributions = read("decomposition_state_contributions")
    panel_d = contributions.loc[
        contributions["dataset"].eq("validation")
        & contributions["score"].eq("core_287")
        & contributions["state_set"].eq("all_states")
        & contributions["comparison"].eq("conventional_vs_normal")
    ].copy()
    panel_d.to_csv(OUT_DIR / "figure3d_state_contributions.tsv", sep="\t", index=False)

    decomposition = read("decomposition_summary")
    panel_e = decomposition.loc[
        decomposition["score"].eq("core_287")
        & decomposition["comparison"].eq("conventional_vs_normal")
        & (
            (decomposition["dataset"].eq("validation") & decomposition["state_set"].isin(["all_states", "canonical_states"]))
            | (decomposition["dataset"].eq("discovery") & decomposition["state_set"].eq("all_states"))
        )
    ].copy()
    panel_e.to_csv(OUT_DIR / "figure3e_exact_decomposition.tsv", sep="\t", index=False)

    donor = read("donor_route_decomposition_inputs", compressed=True)
    panel_f = donor.loc[
        donor["dataset"].eq("validation")
        & donor["state_set"].eq("all_states")
        & donor["route_group"].isin(["normal", "conventional_adenoma"])
    ].copy()
    panel_f = (
        panel_f.groupby(["donor_id", "route_group", "score"], observed=True)["contribution"]
        .sum()
        .rename("total_score")
        .reset_index()
        .pivot(index=["donor_id", "route_group"], columns="score", values="total_score")
        .reset_index()
    )
    panel_f.columns.name = None
    panel_f.to_csv(OUT_DIR / "figure3f_core_compact_concordance.tsv", sep="\t", index=False)

    manifest = []
    for path in sorted(OUT_DIR.glob("figure3*.tsv")):
        frame = pd.read_csv(path, sep="\t")
        manifest.append({"file": path.name, "rows": len(frame), "columns": len(frame.columns)})
    pd.DataFrame(manifest).to_csv(OUT_DIR / "figure3_source_data_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
