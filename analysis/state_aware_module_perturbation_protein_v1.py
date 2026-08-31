#!/usr/bin/env python3
"""Project retained modules into perturbations and nominate experimental proteins."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
RESULT_ROOT = ROOT / "results" / "state_aware_program_v1"
SOURCE_ROOT = RESULT_ROOT / "functional_architecture_v1"
OUT_ROOT = Path(
    os.environ.get("STATE_AWARE_MODULE_RUN_ROOT", str(SOURCE_ROOT))
).resolve()
OUT = OUT_ROOT / "module_perturbation_protein"
MODULE_PATH = OUT_ROOT / "consensus_modules.tsv"
VALIDATION_PATH = OUT_ROOT / "module_validation.tsv"
CANDIDATE_PATH = (
    OUT_ROOT / "module_technical_validation" / "sentinel_candidate_base.tsv"
)
REGULATOR_PATH = (
    RESULT_ROOT
    / "identity_reversal_target_prioritization_v1"
    / "regulator_replication_and_programme_projection.tsv"
)
TWO_COMPONENT_PATH = (
    RESULT_ROOT
    / "identity_reversal_target_prioritization_v1"
    / "perturbation_two_component_summary.tsv"
)
MODULE_COMMUNITY_OVERLAP_PATH = (
    SOURCE_ROOT / "consensus_wgcna" / "module_pathway_community_overlap.tsv"
)
PATHWAY_SUMMARY_PATH = (
    SOURCE_ROOT / "pathway_replication" / "pathway_replication_summary.tsv"
)
ADDENDUM_PATH = (
    ANALYSIS
    / "contracts"
    / "state_aware_functional_architecture_downstream_addendum_v1_2026-08-30.md"
)
EXPECTED_ADDENDUM_HASH = (
    "32343afe117d09007066fbe01f8fbe7cf4a11ee628dbe6a121f0f814968da3bb"
)
ROUTING_ADDENDUM_PATH = (
    Path(
        os.environ.get(
            "STATE_AWARE_MODULE_ROUTING_ADDENDUM_PATH",
            str(
                ANALYSIS
                / "contracts"
                / "state_aware_functional_architecture_routing_addendum_v1_2026-08-30.md"
            ),
        )
    ).resolve()
)
EXPECTED_ROUTING_ADDENDUM_HASH = os.environ.get(
    "STATE_AWARE_MODULE_ROUTING_ADDENDUM_SHA256",
    "45120db0d56cc31e0610e63e1b28f1e09638b5fff78b404f34a0347a7cb62ea2",
)
ROUTE_COLUMN = os.environ.get(
    "STATE_AWARE_MODULE_ROUTE_COLUMN", "internal_gate_pass"
)
MAIN_LABEL = os.environ.get("STATE_AWARE_MODULE_MAIN_LABEL", "Main")
SUPPLEMENT_LABEL = os.environ.get(
    "STATE_AWARE_MODULE_SUPPLEMENT_LABEL", "Supplement"
)
DESCRIPTIVE_DOWNSTREAM = os.environ.get(
    "STATE_AWARE_MODULE_DESCRIPTIVE_DOWNSTREAM", "false"
).lower() == "true"

sys.path.insert(0, str(ANALYSIS))
import computational_closure_validation as closure  # noqa: E402
import public_adenoma_protein_triangulation as protein  # noqa: E402


SEED = 20260830
N_BOOTSTRAPS = 10_000


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


def route_pass(frame: pd.DataFrame) -> pd.Series:
    if ROUTE_COLUMN not in frame.columns:
        raise RuntimeError(f"Routing column is missing: {ROUTE_COLUMN}")
    return as_bool(frame[ROUTE_COLUMN])


def functional_community_increment(validation: pd.DataFrame) -> pd.Series:
    """Apply the frozen direction-matched community-overlap routing rule."""
    community_overlap = pd.read_csv(MODULE_COMMUNITY_OVERLAP_PATH, sep="\t")
    pathway_summary = pd.read_csv(PATHWAY_SUMMARY_PATH, sep="\t")
    replicated_pathways = pathway_summary.loc[
        as_bool(pathway_summary["replicated"]),
        ["community_id", "discovery_common_direction"],
    ].drop_duplicates()
    direction_count = replicated_pathways.groupby(
        "community_id", observed=True
    )["discovery_common_direction"].nunique()
    if (direction_count > 1).any():
        raise RuntimeError("A frozen pathway community is not direction coherent")
    direction = (
        replicated_pathways.drop_duplicates("community_id")
        .set_index("community_id")["discovery_common_direction"]
    )
    output: dict[str, bool] = {}
    for row in validation.itertuples(index=False):
        local = community_overlap.loc[
            community_overlap["module"].eq(row.module)
            & (community_overlap["fdr"] <= 0.05)
        ].copy()
        local["community_direction"] = local["community_id"].map(direction)
        output[row.module] = bool(
            local["community_direction"].eq(row.heldout_direction).any()
        )
    return validation["module"].map(output).fillna(False).astype(bool)


def module_scores(
    expression: pd.DataFrame,
    module_gene_sets: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = closure.aggregate_gene_rows(expression)
    z = closure.gene_zscores(expression)
    scores = pd.DataFrame(index=z.columns)
    coverage_rows: list[dict[str, object]] = []
    for module, genes in module_gene_sets.items():
        present = [
            gene
            for gene in genes
            if gene in z.index and np.isfinite(z.loc[gene].std(ddof=1))
        ]
        coverage = len(present) / len(genes)
        passed = len(present) >= 15 and coverage >= 0.50
        scores[module] = z.loc[present].mean(axis=0) if passed else np.nan
        coverage_rows.append(
            {
                "module": module,
                "module_size": len(genes),
                "n_present": len(present),
                "coverage_fraction": coverage,
                "coverage_gate_pass": passed,
                "present_genes": ";".join(present),
            }
        )
    scores.index.name = "sample_id"
    return scores.reset_index(), pd.DataFrame(coverage_rows)


def comparison_units(
    metadata: pd.DataFrame,
    scores: pd.DataFrame,
    comparisons: list[dict[str, object]],
    dataset: str,
    species: str,
    model_system: str,
    direction_sign: dict[str, int],
) -> pd.DataFrame:
    merged = metadata.merge(scores, on="sample_id", validate="one_to_one")
    records: list[dict[str, object]] = []
    for comparison in comparisons:
        subset = merged
        if comparison.get("units") is not None:
            subset = subset.loc[subset["unit_id"].isin(comparison["units"])]
        for module in direction_sign:
            wide = subset.groupby(["unit_id", "condition"], observed=True)[module].mean().unstack()
            target = comparison["target_condition"]
            reference = comparison["reference_condition"]
            if target not in wide or reference not in wide:
                continue
            difference = (wide[target] - wide[reference]).dropna()
            for unit_id, value in difference.items():
                records.append(
                    {
                        "module": module,
                        "dataset": dataset,
                        "species": species,
                        "model_system": model_system,
                        "comparison": comparison["comparison"],
                        "unit_id": unit_id,
                        "target_minus_reference": float(value),
                        "module_adenoma_direction_sign": direction_sign[module],
                        "module_reversal": float(-direction_sign[module] * value),
                    }
                )
    return pd.DataFrame(records)


def gse135_units(
    module_gene_sets: dict[str, list[str]],
    direction_sign: dict[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    unit_parts: list[pd.DataFrame] = []
    coverage_parts: list[pd.DataFrame] = []
    for cell_line, (expression, metadata) in closure.load_gse135_expression().items():
        scores, coverage = module_scores(expression, module_gene_sets)
        coverage.insert(0, "dataset", f"GSE135328_{cell_line}")
        coverage_parts.append(coverage)
        merged = metadata.merge(scores, on="sample_id", validate="one_to_one")
        for module in direction_sign:
            clone = merged.groupby(["clone_id", "condition"], observed=True)[module].mean().reset_index()
            wild_type = clone.loc[clone["condition"].eq("WT"), module]
            if wild_type.empty:
                continue
            reference = float(wild_type.iloc[0])
            for genotype, comparison in [
                ("Het", "tcf7l2_heterozygous_vs_WT"),
                ("KO", "tcf7l2_KO_vs_WT"),
            ]:
                local = clone.loc[clone["condition"].eq(genotype)]
                for row in local.itertuples(index=False):
                    value = float(getattr(row, module) - reference)
                    unit_parts.append(
                        pd.DataFrame(
                            [
                                {
                                    "module": module,
                                    "dataset": "GSE135328",
                                    "species": "human",
                                    "model_system": "crc_cell_line_clone",
                                    "comparison": comparison,
                                    "unit_id": row.clone_id,
                                    "target_minus_reference": value,
                                    "module_adenoma_direction_sign": direction_sign[module],
                                    "module_reversal": -direction_sign[module] * value,
                                }
                            ]
                        )
                    )
    return (
        pd.concat(unit_parts, ignore_index=True, sort=False),
        pd.concat(coverage_parts, ignore_index=True, sort=False),
    )


def bootstrap_interval(values: np.ndarray, seed: int) -> tuple[float, float]:
    if len(values) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    bootstrap = rng.choice(values, size=(N_BOOTSTRAPS, len(values)), replace=True).mean(axis=1)
    return float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))


def summarize_perturbations(units: pd.DataFrame) -> pd.DataFrame:
    role = pd.read_csv(TWO_COMPONENT_PATH, sep="\t")[
        ["dataset", "comparison", "interpretation_role"]
    ].drop_duplicates()
    rows: list[dict[str, object]] = []
    for keys, frame in units.groupby(
        ["module", "dataset", "species", "model_system", "comparison"],
        observed=True,
    ):
        module, dataset, species, model_system, comparison = keys
        values = frame["module_reversal"].to_numpy(float)
        ci_low, ci_high = bootstrap_interval(values, SEED + len(rows))
        nonzero = values[values != 0]
        sign_p = (
            stats.binomtest(
                int((nonzero > 0).sum()), len(nonzero), 0.5, alternative="two-sided"
            ).pvalue
            if len(nonzero)
            else np.nan
        )
        rows.append(
            {
                "module": module,
                "dataset": dataset,
                "species": species,
                "model_system": model_system,
                "comparison": comparison,
                "n_units": len(values),
                "mean_module_reversal": float(np.mean(values)),
                "median_module_reversal": float(np.median(values)),
                "bootstrap_mean_ci_low": ci_low,
                "bootstrap_mean_ci_high": ci_high,
                "n_positive_reversal": int((values > 0).sum()),
                "all_units_positive_reversal": bool((values > 0).all()),
                "p_exact_sign_two_sided": float(sign_p),
            }
        )
    output = pd.DataFrame(rows).merge(
        role, on=["dataset", "comparison"], how="left", validate="many_to_one"
    )
    output["interpretation_role"] = output["interpretation_role"].fillna(
        "causal or pathway perturbation"
    )
    return output


def perturbation_gate(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module, frame in summary.groupby("module", observed=True):
        causal = frame.loc[
            frame["interpretation_role"].eq("causal or pathway perturbation")
        ]
        positive_datasets = causal.loc[
            causal["mean_module_reversal"] > 0, "dataset"
        ].nunique()
        multi_unit_all_positive = (
            (causal["n_units"] >= 3) & causal["all_units_positive_reversal"]
        ).any()
        rows.append(
            {
                "module": module,
                "n_independent_datasets_positive_reversal": int(positive_datasets),
                "one_model_n3_all_positive": bool(multi_unit_all_positive),
                "perturbation_gate_pass": positive_datasets >= 2
                or multi_unit_all_positive,
            }
        )
    return pd.DataFrame(rows)


def spatial_support(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = candidates["gene"].tolist()
    section_rows: list[dict[str, object]] = []
    for sample_dir in sorted(
        path
        for path in closure.SPATIAL_DIR.iterdir()
        if path.is_dir() and (path / "filtered_feature_bc_matrix.h5").exists()
    ):
        sample_id = sample_dir.name
        barcodes, all_genes, matrix = closure.read_10x_h5(
            sample_dir / "filtered_feature_bc_matrix.h5"
        )
        gene_index = {gene: index for index, gene in enumerate(all_genes)}
        present = [gene for gene in genes if gene in gene_index]
        total = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
        selected = matrix[[gene_index[gene] for gene in present], :].toarray().astype(float)
        expression = pd.DataFrame(
            np.log2(selected / np.where(total > 0, total, np.nan) * 1_000_000 + 1),
            index=present,
            columns=barcodes,
        )
        annotation_path = (
            closure.SPATIAL_ANNOT_DIR / f"Pathologist_Annotations_{sample_id}.csv"
        )
        annotation = pd.read_csv(annotation_path)
        annotation_column = next(
            (
                column
                for column in annotation.columns
                if column.lower().startswith("pathologist annotation")
            ),
            None,
        )
        if annotation_column is None:
            annotation_column = next(column for column in annotation.columns if column != "Barcode")
        annotation = annotation.rename(
            columns={"Barcode": "barcode", annotation_column: "pathology_annotation"}
        )
        annotation["pathology_group"] = annotation["pathology_annotation"].map(
            closure.pathology_group
        )
        for gene in genes:
            if gene not in expression.index:
                continue
            values = pd.DataFrame(
                {"barcode": expression.columns, "expression": expression.loc[gene].to_numpy()}
            ).merge(annotation[["barcode", "pathology_group"]], on="barcode", how="inner")
            tumour = values.loc[
                values["pathology_group"].isin(["tumor", "tumor_stroma"]),
                "expression",
            ]
            normal = values.loc[
                values["pathology_group"].eq("non_neoplastic_epithelium"),
                "expression",
            ]
            if tumour.empty or normal.empty:
                continue
            section_rows.append(
                {
                    "gene": gene,
                    "sample_id": sample_id,
                    "n_tumour_spots": len(tumour),
                    "n_nonneoplastic_epithelial_spots": len(normal),
                    "tumour_minus_nonneoplastic_median": float(
                        tumour.median() - normal.median()
                    ),
                }
            )
    sections = pd.DataFrame(section_rows)
    summary_rows = []
    expected = candidates.set_index("gene")["expected_direction"]
    for gene in genes:
        local = sections.loc[sections["gene"].eq(gene)]
        values = local["tumour_minus_nonneoplastic_median"].to_numpy(float)
        direction = int(expected[gene])
        oriented = direction * values
        test = (
            stats.wilcoxon(oriented, zero_method="wilcox", alternative="two-sided")
            if len(oriented) >= 4 and (oriented != 0).any()
            else None
        )
        summary_rows.append(
            {
                "gene": gene,
                "expected_direction": direction,
                "n_informative_sections": len(oriented),
                "median_direction_oriented_difference": float(np.median(oriented))
                if len(oriented)
                else np.nan,
                "direction_concordant_fraction": float((oriented > 0).mean())
                if len(oriented)
                else np.nan,
                "p_paired_wilcoxon": float(test.pvalue) if test is not None else np.nan,
                "spatial_directional_support": len(oriented) >= 4
                and np.median(oriented) > 0
                and (oriented > 0).mean() >= 2 / 3
                and test is not None
                and test.pvalue <= 0.10,
            }
        )
    return sections, pd.DataFrame(summary_rows)


def protein_and_spatial_candidates(
    candidate_base: pd.DataFrame,
    validation: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    candidates = candidate_base.loc[as_bool(candidate_base["pre_orthogonal_gate_pass"])].copy()
    if candidates.empty:
        return candidates, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    module_direction = validation.set_index("module")["heldout_direction"]
    candidates["expected_direction"] = candidates["module"].map(
        lambda module: 1 if module_direction[module] == "Up" else -1
    )
    inventory = candidates[["gene", "expected_direction"]].drop_duplicates()
    pxd2 = protein.pxd002137_tests(inventory)
    pxd17 = protein.pxd017269_detectability(inventory)
    pxd46 = protein.pxd046999_presence(inventory)
    spatial_sections, spatial = spatial_support(inventory)
    evidence = (
        candidates.merge(pxd2, on="gene", how="left")
        .merge(pxd17, on="gene", how="left")
        .merge(pxd46, on="gene", how="left")
        .merge(spatial, on=["gene", "expected_direction"], how="left")
    )
    evidence["protein_directional_support"] = (
        evidence["protein_gene_present"].fillna(False)
        & (np.sign(evidence["age_sex_adjusted_log2_effect"]) == evidence["expected_direction"])
        & (evidence["age_sex_adjusted_p"] <= 0.10)
    )
    evidence["spatial_directional_support"] = evidence[
        "spatial_directional_support"
    ].fillna(False)
    evidence["orthogonal_evidence_count"] = (
        evidence["protein_directional_support"].astype(int)
        + evidence["spatial_directional_support"].astype(int)
    )
    evidence["orthogonal_gate_pass"] = evidence["orthogonal_evidence_count"] >= 1
    return evidence, pxd2, pxd17, pxd46, spatial_sections


def non_dominated_sentinels(evidence: pd.DataFrame) -> pd.DataFrame:
    eligible = evidence.loc[evidence["orthogonal_gate_pass"]].copy()
    records = []
    for module, frame in eligible.groupby("module", observed=True):
        dominated = pd.Series(False, index=frame.index)
        for index, row in frame.iterrows():
            other = frame.drop(index=index)
            weakly_better = (
                (other["bootstrap_kme_states_ci_low_gt_0.40"] >= row["bootstrap_kme_states_ci_low_gt_0.40"])
                & (other["heldout_common_fdr"] <= row["heldout_common_fdr"])
                & (other["external_cohorts_measurable"] >= row["external_cohorts_measurable"])
                & (other["orthogonal_evidence_count"] >= row["orthogonal_evidence_count"])
            )
            strictly_better = (
                (other["bootstrap_kme_states_ci_low_gt_0.40"] > row["bootstrap_kme_states_ci_low_gt_0.40"])
                | (other["heldout_common_fdr"] < row["heldout_common_fdr"])
                | (other["external_cohorts_measurable"] > row["external_cohorts_measurable"])
                | (other["orthogonal_evidence_count"] > row["orthogonal_evidence_count"])
            )
            dominated.loc[index] = bool((weakly_better & strictly_better).any())
        frame["pareto_non_dominated"] = ~dominated
        selected = (
            frame.loc[frame["pareto_non_dominated"]]
            .sort_values(
                [
                    "bootstrap_kme_states_ci_low_gt_0.40",
                    "heldout_common_fdr",
                    "external_cohorts_measurable",
                    "orthogonal_evidence_count",
                    "gene",
                ],
                ascending=[False, True, False, False, True],
            )
            .head(3)
            .copy()
        )
        selected["priority_role"] = "measurement_sentinel"
        selected["final_status"] = "nominated"
        records.append(selected)
    return pd.concat(records, ignore_index=True, sort=False) if records else pd.DataFrame()


def regulatory_nodes(
    modules: pd.DataFrame,
    validation: pd.DataFrame,
    perturbation_summary: pd.DataFrame,
    perturbation_gate_table: pd.DataFrame,
) -> pd.DataFrame:
    retained = validation.loc[
        route_pass(validation)
        & validation["external_gate_pass"]
        & validation["technical_gate_pass"],
        "module",
    ]
    gene_module = modules.loc[modules["module"].isin(retained)].set_index("gene")["module"]
    gate_by_module = perturbation_gate_table.set_index("module")["perturbation_gate_pass"]
    regulator = pd.read_csv(REGULATOR_PATH, sep="\t")
    regulator["replicated"] = as_bool(regulator["replicated"])
    rows: list[dict[str, object]] = []
    for row in regulator.loc[regulator["replicated"]].itertuples(index=False):
        gene = row.source
        if gene not in gene_module.index:
            continue
        module = gene_module[gene]
        gene_row = modules.loc[modules["gene"].eq(gene)].iloc[0]
        expression_direction = np.sign(gene_row["heldout_common_z"])
        activity_direction = -1 if row.activity_direction == "suppressed_in_adenoma" else 1
        no_conflict = expression_direction == activity_direction
        if (
            gate_by_module.get(module, False)
            and gene_row["heldout_common_fdr"] <= 0.10
            and no_conflict
        ):
            rows.append(
                {
                    "gene": gene,
                    "module": module,
                    "priority_role": "regulatory_node",
                    "regulatory_evidence": "replicated_signed_regulon",
                    "proposed_intervention": row.proposed_intervention,
                    "expression_activity_direction_consistent": True,
                    "final_status": "nominated",
                }
            )

    direct_map = {
        "ASCL2": ("GSE130822", "ascl2_ko_vs_resting_wt", "suppress"),
        "TCF7L2": ("GSE135328", "tcf7l2_KO_vs_WT", "suppress"),
    }
    for gene, (dataset, comparison, intervention) in direct_map.items():
        if gene not in gene_module.index:
            continue
        module = gene_module[gene]
        gene_row = modules.loc[modules["gene"].eq(gene)].iloc[0]
        specific = perturbation_summary.loc[
            perturbation_summary["module"].eq(module)
            & perturbation_summary["dataset"].eq(dataset)
            & perturbation_summary["comparison"].eq(comparison)
        ]
        no_conflict = gene_row["heldout_common_z"] > 0 and intervention == "suppress"
        if (
            not specific.empty
            and specific["mean_module_reversal"].iloc[0] > 0
            and gate_by_module.get(module, False)
            and gene_row["heldout_common_fdr"] <= 0.10
            and no_conflict
        ):
            rows.append(
                {
                    "gene": gene,
                    "module": module,
                    "priority_role": "regulatory_node",
                    "regulatory_evidence": f"independent_perturbation:{dataset}:{comparison}",
                    "proposed_intervention": intervention,
                    "expression_activity_direction_consistent": True,
                    "final_status": "nominated",
                }
            )
    if not rows:
        return pd.DataFrame(
            columns=[
                "gene",
                "module",
                "priority_role",
                "regulatory_evidence",
                "proposed_intervention",
                "expression_activity_direction_consistent",
                "final_status",
            ]
        )
    return pd.DataFrame(rows).drop_duplicates(["gene", "module", "priority_role"])


def main() -> None:
    required = [
        MODULE_PATH,
        VALIDATION_PATH,
        CANDIDATE_PATH,
        REGULATOR_PATH,
        TWO_COMPONENT_PATH,
        MODULE_COMMUNITY_OVERLAP_PATH,
        PATHWAY_SUMMARY_PATH,
        ADDENDUM_PATH,
        ROUTING_ADDENDUM_PATH,
    ]
    if not all(path.exists() for path in required):
        raise RuntimeError("At least one perturbation/protein input is missing")
    if sha256(ADDENDUM_PATH) != EXPECTED_ADDENDUM_HASH:
        raise RuntimeError("The frozen downstream validation addendum changed")
    if sha256(ROUTING_ADDENDUM_PATH) != EXPECTED_ROUTING_ADDENDUM_HASH:
        raise RuntimeError("The frozen final-routing addendum changed")
    OUT.mkdir(parents=True, exist_ok=True)
    modules = pd.read_csv(MODULE_PATH, sep="\t")
    validation_input_hash = sha256(VALIDATION_PATH)
    validation = pd.read_csv(VALIDATION_PATH, sep="\t")
    validation["functional_community_increment_pass"] = (
        functional_community_increment(validation)
    )
    route_mask = route_pass(validation)
    if DESCRIPTIVE_DOWNSTREAM:
        retained_mask = route_mask
    else:
        retained_mask = (
            route_mask
            & as_bool(validation["external_gate_pass"])
            & as_bool(validation["technical_gate_pass"])
        )
    retained = validation.loc[retained_mask, "module"].tolist()
    if not retained:
        validation["perturbation_gate_pass"] = False
        validation["n_sentinel_proteins"] = 0
        validation["n_regulatory_nodes"] = 0
        validation["interpretive_increment_pass"] = validation[
            "functional_community_increment_pass"
        ]
        validation["routing_status"] = "Audit"
        validation.to_csv(VALIDATION_PATH, sep="\t", index=False)
        priority_columns = [
            "gene", "module", "priority_role", "regulatory_evidence",
            "proposed_intervention", "final_status",
        ]
        pd.DataFrame(columns=priority_columns).to_csv(
            OUT_ROOT / "protein_priorities.tsv", sep="\t", index=False
        )
        empty_outputs = {
            "module_perturbation_unit_effects.tsv": [
                "module", "dataset", "species", "model_system", "comparison",
                "unit_id", "target_minus_reference", "module_reversal",
            ],
            "module_perturbation_summary.tsv": [
                "module", "dataset", "comparison", "n_units",
                "mean_module_reversal", "interpretation_role",
            ],
            "module_perturbation_gate.tsv": [
                "module", "n_independent_datasets_positive_reversal",
                "one_model_n3_all_positive", "perturbation_gate_pass",
            ],
            "module_perturbation_coverage.tsv": [
                "dataset", "module", "module_size", "n_present",
                "coverage_fraction", "coverage_gate_pass",
            ],
            "sentinel_orthogonal_evidence.tsv": [
                "gene", "module", "orthogonal_gate_pass",
            ],
            "sentinel_pxd002137_tests.tsv": ["gene"],
            "sentinel_pxd017269_detectability.tsv": ["gene"],
            "sentinel_pxd046999_presence.tsv": ["gene"],
            "sentinel_spatial_section_effects.tsv": ["gene", "sample_id"],
        }
        for filename, columns in empty_outputs.items():
            pd.DataFrame(columns=columns).to_csv(OUT / filename, sep="\t", index=False)
        manifest = {
            "analysis": "state_aware_module_perturbation_protein_v1",
            "created_utc": pd.Timestamp.utcnow().isoformat(),
            "random_seed": SEED,
            "modules_entering_analysis": [],
            "modules_routed_main": [],
            "sentinel_proteins": [],
            "regulatory_nodes": [],
            "input_sha256": {
                "consensus_modules": sha256(MODULE_PATH),
                "module_validation_before_final_routing": validation_input_hash,
                "module_community_overlap": sha256(MODULE_COMMUNITY_OVERLAP_PATH),
                "pathway_summary": sha256(PATHWAY_SUMMARY_PATH),
                "downstream_addendum": sha256(ADDENDUM_PATH),
                "routing_addendum": sha256(ROUTING_ADDENDUM_PATH),
            },
            "output_sha256": {
                "module_validation": sha256(VALIDATION_PATH),
                "protein_priorities": sha256(OUT_ROOT / "protein_priorities.tsv"),
            },
        }
        (OUT / "module_perturbation_protein_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        print("No retained module entered perturbation/protein nomination")
        return
    candidate_base = pd.read_csv(CANDIDATE_PATH, sep="\t")

    module_gene_sets = {
        module: modules.loc[modules["module"].eq(module), "gene"].tolist()
        for module in retained
    }
    direction_sign = {
        row.module: 1 if row.heldout_direction == "Up" else -1
        for row in validation.loc[validation["module"].isin(retained)].itertuples()
    }
    homology, _ = closure.strict_human_mouse_map()
    datasets = {
        "GSE114059": (
            "human",
            "patient_derived_crc_organoid",
            *closure.load_gse114059(),
        ),
        "GSE67186": (
            "mouse",
            "in_vivo_colon_polyp",
            *closure.load_gse67186(homology),
        ),
        "GSE130822": (
            "mouse",
            "colonic_stem_cells",
            *closure.load_gse130822(homology),
        ),
        "GSE171910": (
            "human",
            "conditional_crc_wnt_models",
            *closure.load_gse171910(),
        ),
    }
    unit_parts: list[pd.DataFrame] = []
    coverage_parts: list[pd.DataFrame] = []
    for dataset, (species, model_system, expression, metadata, comparisons) in datasets.items():
        scores, coverage = module_scores(expression, module_gene_sets)
        coverage.insert(0, "dataset", dataset)
        coverage_parts.append(coverage)
        unit_parts.append(
            comparison_units(
                metadata,
                scores,
                comparisons,
                dataset,
                species,
                model_system,
                direction_sign,
            )
        )

    gse125_expression, gse125_metadata, _ = closure.load_gse125_expression()
    gse125_comparisons = [
        {
            "comparison": "wnt_rspo_withdrawal_in_WT",
            "target_condition": "WT_without",
            "reference_condition": "WT_with",
            "units": None,
        },
        {
            "comparison": "wnt_rspo_withdrawal_in_APC_KO",
            "target_condition": "APC_without",
            "reference_condition": "APC_with",
            "units": None,
        },
    ]
    scores, coverage = module_scores(gse125_expression, module_gene_sets)
    coverage.insert(0, "dataset", "GSE125472")
    coverage_parts.append(coverage)
    unit_parts.append(
        comparison_units(
            gse125_metadata,
            scores,
            gse125_comparisons,
            "GSE125472",
            "human",
            "isogenic_human_colon_organoid",
            direction_sign,
        )
    )
    gse135_unit_table, gse135_coverage = gse135_units(
        module_gene_sets, direction_sign
    )
    unit_parts.append(gse135_unit_table)
    coverage_parts.append(gse135_coverage)

    units = pd.concat(unit_parts, ignore_index=True, sort=False)
    coverage = pd.concat(coverage_parts, ignore_index=True, sort=False)
    perturbation_summary = summarize_perturbations(units)
    perturbation_gate_table = perturbation_gate(perturbation_summary)
    evidence, pxd2, pxd17, pxd46, spatial_sections = protein_and_spatial_candidates(
        candidate_base, validation
    )
    sentinels = non_dominated_sentinels(evidence) if not evidence.empty else pd.DataFrame()
    if DESCRIPTIVE_DOWNSTREAM and not sentinels.empty:
        fully_validated_modules = set(
            validation.loc[
                route_pass(validation)
                & as_bool(validation["external_gate_pass"])
                & as_bool(validation["technical_gate_pass"]),
                "module",
            ]
        )
        sentinels.loc[
            ~sentinels["module"].isin(fully_validated_modules), "final_status"
        ] = "exploratory_only"
    regulators = regulatory_nodes(
        modules, validation, perturbation_summary, perturbation_gate_table
    )
    priorities = pd.concat([sentinels, regulators], ignore_index=True, sort=False)

    validation = validation.merge(
        perturbation_gate_table, on="module", how="left"
    )
    validation["perturbation_gate_pass"] = validation[
        "perturbation_gate_pass"
    ].fillna(False)
    sentinel_counts = (
        sentinels.groupby("module").size().rename("n_sentinel_proteins")
        if not sentinels.empty
        else pd.Series(dtype=int, name="n_sentinel_proteins")
    )
    regulator_counts = (
        regulators.groupby("module").size().rename("n_regulatory_nodes")
        if not regulators.empty
        else pd.Series(dtype=int, name="n_regulatory_nodes")
    )
    validation = validation.merge(
        sentinel_counts, left_on="module", right_index=True, how="left"
    ).merge(regulator_counts, left_on="module", right_index=True, how="left")
    validation[["n_sentinel_proteins", "n_regulatory_nodes"]] = validation[
        ["n_sentinel_proteins", "n_regulatory_nodes"]
    ].fillna(0).astype(int)
    validation["interpretive_increment_pass"] = (
        validation["functional_community_increment_pass"]
        | validation["perturbation_gate_pass"]
        | (validation["n_sentinel_proteins"] > 0)
        | (validation["n_regulatory_nodes"] > 0)
    )
    core = (
        route_pass(validation)
        & as_bool(validation["external_gate_pass"])
        & as_bool(validation["technical_gate_pass"])
    )
    validation["routing_status"] = np.select(
        [core & validation["interpretive_increment_pass"], core],
        [MAIN_LABEL, SUPPLEMENT_LABEL],
        default="Audit",
    )

    units.to_csv(OUT / "module_perturbation_unit_effects.tsv", sep="\t", index=False)
    perturbation_summary.to_csv(
        OUT / "module_perturbation_summary.tsv", sep="\t", index=False
    )
    perturbation_gate_table.to_csv(
        OUT / "module_perturbation_gate.tsv", sep="\t", index=False
    )
    coverage.to_csv(OUT / "module_perturbation_coverage.tsv", sep="\t", index=False)
    evidence.to_csv(OUT / "sentinel_orthogonal_evidence.tsv", sep="\t", index=False)
    pxd2.to_csv(OUT / "sentinel_pxd002137_tests.tsv", sep="\t", index=False)
    pxd17.to_csv(OUT / "sentinel_pxd017269_detectability.tsv", sep="\t", index=False)
    pxd46.to_csv(OUT / "sentinel_pxd046999_presence.tsv", sep="\t", index=False)
    spatial_sections.to_csv(
        OUT / "sentinel_spatial_section_effects.tsv", sep="\t", index=False
    )
    priorities.to_csv(OUT_ROOT / "protein_priorities.tsv", sep="\t", index=False)
    validation.to_csv(VALIDATION_PATH, sep="\t", index=False)

    manifest = {
        "analysis": "state_aware_module_perturbation_protein_v1",
        "created_utc": pd.Timestamp.utcnow().isoformat(),
        "random_seed": SEED,
        "modules_entering_analysis": retained,
        "route_column": ROUTE_COLUMN,
        "descriptive_downstream": DESCRIPTIVE_DOWNSTREAM,
        "modules_routed_main": validation.loc[
            validation["routing_status"].eq(MAIN_LABEL), "module"
        ].tolist(),
        "sentinel_proteins": sentinels.get("gene", pd.Series(dtype=str)).tolist(),
        "regulatory_nodes": regulators.get("gene", pd.Series(dtype=str)).tolist(),
        "input_sha256": {
            "consensus_modules": sha256(MODULE_PATH),
            "module_validation_before_final_routing": validation_input_hash,
            "sentinel_candidate_base": sha256(CANDIDATE_PATH),
            "regulator_evidence": sha256(REGULATOR_PATH),
            "two_component_summary": sha256(TWO_COMPONENT_PATH),
            "module_community_overlap": sha256(MODULE_COMMUNITY_OVERLAP_PATH),
            "pathway_summary": sha256(PATHWAY_SUMMARY_PATH),
            "downstream_addendum": sha256(ADDENDUM_PATH),
            "routing_addendum": sha256(ROUTING_ADDENDUM_PATH),
        },
        "output_sha256": {
            "module_validation": sha256(VALIDATION_PATH),
            "protein_priorities": sha256(OUT_ROOT / "protein_priorities.tsv"),
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": __import__("scipy").__version__,
        },
    }
    (OUT / "module_perturbation_protein_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Module perturbation/protein analysis complete: "
        f"{sum(validation['routing_status'].eq(MAIN_LABEL))} modules routed {MAIN_LABEL}; "
        f"{len(sentinels)} sentinels; {len(regulators)} regulatory nodes"
    )


if __name__ == "__main__":
    main()
