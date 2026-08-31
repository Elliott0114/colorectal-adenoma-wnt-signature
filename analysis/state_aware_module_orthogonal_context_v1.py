#!/usr/bin/env python3
"""Project retained consensus modules into orthogonal tissue contexts."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
RESULT_ROOT = ROOT / "results" / "state_aware_program_v1"
SOURCE_ROOT = RESULT_ROOT / "functional_architecture_v1"
ARCHITECTURE_ROOT = Path(
    os.environ.get("STATE_AWARE_MODULE_RUN_ROOT", str(SOURCE_ROOT))
).resolve()
OUT = ARCHITECTURE_ROOT / "module_orthogonal_context"
MODULE_PATH = ARCHITECTURE_ROOT / "consensus_modules.tsv"
VALIDATION_PATH = ARCHITECTURE_ROOT / "module_validation.tsv"
ROUTE_COLUMN = os.environ.get(
    "STATE_AWARE_MODULE_ROUTE_COLUMN", "internal_gate_pass"
)
DESCRIPTIVE_DOWNSTREAM = os.environ.get(
    "STATE_AWARE_MODULE_DESCRIPTIVE_DOWNSTREAM", "false"
).lower() == "true"
PAIRED_ATAC_PATH = (
    RESULT_ROOT
    / "extended_validation_full_programme"
    / "becker_rna_atac"
    / "becker_rna_atac_paired_scores.tsv"
)
CONTRACT_PATH = (
    ANALYSIS
    / "contracts"
    / "state_aware_functional_architecture_context_addendum_v1_2026-08-30.md"
)
EXPECTED_CONTRACT_HASH = (
    "e24e2fcbcc1ce99b145287eb6b1d0f09f27a1853c75f6152ca8df337fcdf3073"
)

sys.path.insert(0, str(ANALYSIS))
import computational_closure_validation as closure  # noqa: E402
import conventional_route_signature_transfer as route  # noqa: E402
import public_adenoma_protein_triangulation as protein  # noqa: E402


SEED = 20260830
MIN_GENES = 15
MIN_COVERAGE = 0.50


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


def score_modules(
    expression: pd.DataFrame,
    module_genes: dict[str, list[str]],
    direction: dict[str, int],
    context: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = expression.astype(float)
    standard_deviation = expression.std(axis=0, ddof=1)
    variable = set(standard_deviation.index[np.isfinite(standard_deviation) & (standard_deviation > 0)])
    z = expression.sub(expression.mean(axis=0), axis=1).div(
        standard_deviation.replace(0, np.nan), axis=1
    )
    scores = pd.DataFrame(index=expression.index)
    coverage_rows: list[dict[str, object]] = []
    for module, genes in module_genes.items():
        present = [gene for gene in genes if gene in expression.columns and gene in variable]
        coverage = len(present) / len(genes)
        passed = len(present) >= MIN_GENES and coverage >= MIN_COVERAGE
        scores[module] = direction[module] * z[present].mean(axis=1) if passed else np.nan
        coverage_rows.append(
            {
                "context": context,
                "module": module,
                "module_size": len(genes),
                "n_measurable": len(present),
                "coverage_fraction": coverage,
                "coverage_gate_pass": passed,
                "measurable_genes": ";".join(present),
            }
        )
    scores.index.name = expression.index.name or "sample_id"
    return scores, pd.DataFrame(coverage_rows)


def clustered_binary_effect(
    frame: pd.DataFrame,
    score_column: str,
    group_column: str,
    cluster_column: str,
    positive_group: str,
    reference_group: str,
    covariates: list[str] | None = None,
) -> dict[str, object]:
    covariates = covariates or []
    columns = [score_column, group_column, cluster_column, *covariates]
    data = frame.loc[
        frame[group_column].isin([positive_group, reference_group]), columns
    ].replace([np.inf, -np.inf], np.nan).dropna().copy()
    data["is_positive"] = data[group_column].eq(positive_group).astype(float)
    score_sd = data[score_column].std(ddof=1)
    if len(data) < 8 or not np.isfinite(score_sd) or score_sd <= 0:
        return {
            "n": len(data),
            "n_clusters": data[cluster_column].nunique(),
            "effect_sd": np.nan,
            "standard_error": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p_value": np.nan,
        }
    y = (data[score_column] - data[score_column].mean()) / score_sd
    x = pd.DataFrame({"is_positive": data["is_positive"]}, index=data.index)
    if covariates:
        numeric = [column for column in covariates if pd.api.types.is_numeric_dtype(data[column])]
        categorical = [column for column in covariates if column not in numeric]
        for column in numeric:
            x[column] = data[column].astype(float)
        if categorical:
            x = pd.concat(
                [x, pd.get_dummies(data[categorical], drop_first=True, dtype=float)],
                axis=1,
            )
    x = sm.add_constant(x.astype(float), has_constant="add")
    fit = sm.OLS(y, x).fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data[cluster_column].astype(str),
            "use_correction": True,
            "df_correction": True,
        },
        use_t=True,
    )
    interval = fit.conf_int().loc["is_positive"]
    return {
        "n": len(data),
        "n_clusters": data[cluster_column].nunique(),
        "effect_sd": float(fit.params["is_positive"]),
        "standard_error": float(fit.bse["is_positive"]),
        "ci_low": float(interval.iloc[0]),
        "ci_high": float(interval.iloc[1]),
        "p_value": float(fit.pvalues["is_positive"]),
    }


def load_becker_module_expression(wanted: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = route.parse_becker_series(route.BECKER_SERIES).merge(
        route.becker_tar_prefixes(route.BECKER_TAR),
        on="geo_accession",
        how="left",
        validate="one_to_one",
    )
    rows: list[dict[str, object]] = []
    with tarfile.open(route.BECKER_TAR, "r") as archive:
        for row in metadata.sort_values("geo_accession").itertuples(index=False):
            genes = route.read_features(archive, f"{row.file_prefix}_features.tsv.gz")
            indices: dict[str, list[int]] = defaultdict(list)
            for index, gene in enumerate(genes):
                indices[str(gene)].append(index)
            matrix = route.read_matrix(archive, f"{row.file_prefix}_matrix.mtx.gz")
            total_per_nucleus = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
            epithelial_indices = [
                index
                for gene in route.EPITHELIAL_MARKERS
                for index in indices.get(gene, [])
            ]
            epithelial_mask = (
                np.asarray(matrix[epithelial_indices, :].sum(axis=0)).ravel() > 0
                if epithelial_indices
                else np.zeros(matrix.shape[1], dtype=bool)
            )
            library_size = float(total_per_nucleus[epithelial_mask].sum())
            output: dict[str, object] = {
                "geo_accession": row.geo_accession,
                "n_nuclei": int(matrix.shape[1]),
                "n_epithelial_marker_positive": int(epithelial_mask.sum()),
                "epithelial_marker_positive_fraction": float(epithelial_mask.mean()),
            }
            present = sorted(gene for gene in wanted if gene in indices)
            for gene in present:
                count = float(matrix[indices[gene][0], :][:, epithelial_mask].sum())
                output[gene] = float(route.log_cpm(np.array([count]), library_size)[0])
            rows.append(output)
    expression = pd.DataFrame(rows).set_index("geo_accession")
    technical_columns = [
        "n_nuclei",
        "n_epithelial_marker_positive",
        "epithelial_marker_positive_fraction",
    ]
    technical = expression[technical_columns].reset_index()
    metadata = (
        metadata.set_index("geo_accession")
        .loc[expression.index]
        .reset_index()
        .merge(technical, on="geo_accession", validate="one_to_one")
    )
    return metadata, expression.drop(columns=technical_columns)


def becker_context(
    module_genes: dict[str, list[str]],
    direction: dict[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    wanted = set().union(*map(set, module_genes.values()))
    metadata, expression = load_becker_module_expression(wanted)
    scores, coverage = score_modules(expression, module_genes, direction, "Becker_epithelial_RNA")
    frame = metadata.merge(scores.reset_index(), on="geo_accession", validate="one_to_one")
    frame["log10_n_nuclei"] = np.log10(frame["n_nuclei"].clip(lower=1))
    effect_rows = []
    for module in module_genes:
        value = clustered_binary_effect(
            frame,
            module,
            "disease_stage_group",
            "patient_id",
            "polyp",
            "normal_unaffected",
            [
                "log10_n_nuclei",
                "epithelial_marker_positive_fraction",
                "familial_adenomatous_polyposis",
                "sex",
            ],
        )
        effect_rows.append({"module": module, **value})

    paired = pd.read_csv(PAIRED_ATAC_PATH, sep="\t")
    paired = paired.loc[paired["analysis_set"].eq("normal_polyp")].merge(
        scores.reset_index().rename(columns={"geo_accession": "scrna_geo_accession"}),
        on="scrna_geo_accession",
        how="inner",
        validate="one_to_one",
    )
    correlation_rows = []
    for module in module_genes:
        sample = paired[[module, "atac_tss__wnt_stem"]].dropna()
        sample_correlation = stats.spearmanr(
            sample[module], sample["atac_tss__wnt_stem"]
        )
        patient = (
            paired.groupby("patient_id", observed=True)[[module, "atac_tss__wnt_stem"]]
            .median()
            .dropna()
        )
        patient_correlation = stats.spearmanr(
            patient[module], patient["atac_tss__wnt_stem"]
        )
        correlation_rows.append(
            {
                "module": module,
                "n_samples": len(sample),
                "n_patients": len(patient),
                "sample_spearman_rho": float(sample_correlation.statistic),
                "sample_spearman_p": float(sample_correlation.pvalue),
                "patient_median_spearman_rho": float(patient_correlation.statistic),
                "patient_median_spearman_p": float(patient_correlation.pvalue),
            }
        )
    return frame, coverage, pd.DataFrame(effect_rows), pd.DataFrame(correlation_rows)


def atlas_context(
    module_genes: dict[str, list[str]],
    direction: dict[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    wanted = sorted(set().union(*map(set, module_genes.values())))
    with h5py.File(
        route.ATLAS_PATH,
        "r",
        rdcc_nbytes=512 * 1024 * 1024,
        rdcc_nslots=1_000_003,
    ) as archive:
        sample_type = route.read_h5_col(archive["obs"], "sample_type")
        coarse_type = route.read_h5_col(archive["obs"], "cell_type_coarse_crc_atlas")
        carrier = route.atlas_carrier_group(sample_type, coarse_type)
        rows = route.sample_atlas_rows(carrier)
        donor = route.read_h5_col(archive["obs"], "donor_id")[rows]
        study = route.read_h5_col(archive["obs"], "study_id")[rows]
        carrier_selected = carrier[rows]
        group_frame = pd.DataFrame(
            {"donor_id": donor, "carrier_group": carrier_selected, "study_id": study}
        )
        group_frame["group_id"] = (
            group_frame["donor_id"].astype(str)
            + "__"
            + group_frame["carrier_group"].astype(str)
        )
        group_ids = pd.Index(group_frame["group_id"].drop_duplicates())
        group_lookup = {value: index for index, value in enumerate(group_ids)}
        row_group = group_frame["group_id"].map(group_lookup).to_numpy(int)

        gene_names = route.read_h5_col(archive["var"], "feature_name")
        gene_to_index = {str(gene): index for index, gene in enumerate(gene_names)}
        present = [gene for gene in wanted if gene in gene_to_index]
        selected_indices = np.array([gene_to_index[gene] for gene in present], dtype=np.int64)
        feature_lookup = np.full(len(gene_names), -1, dtype=np.int32)
        feature_lookup[selected_indices] = np.arange(len(selected_indices), dtype=np.int32)
        sums = np.zeros((len(group_ids), len(present)), dtype=np.float64)
        counts = np.zeros(len(group_ids), dtype=np.int64)
        x = archive["X"]
        indptr = x["indptr"][()]
        indices = x["indices"]
        data = x["data"]
        for output_index, row_index in enumerate(rows):
            start = int(indptr[row_index])
            end = int(indptr[row_index + 1])
            local_index = indices[start:end]
            selected_position = feature_lookup[local_index]
            keep = selected_position >= 0
            if keep.any():
                sums[row_group[output_index], selected_position[keep]] += data[start:end][keep]
            counts[row_group[output_index]] += 1
    expression = pd.DataFrame(
        sums / np.maximum(counts[:, None], 1), index=group_ids, columns=present
    )
    group_metadata = (
        group_frame.groupby("group_id", observed=True)
        .agg(
            donor_id=("donor_id", "first"),
            carrier_group=("carrier_group", "first"),
            study_id=("study_id", lambda values: ";".join(sorted(set(map(str, values))))),
        )
        .loc[group_ids]
        .reset_index()
    )
    group_metadata["n_cells_sampled"] = counts
    scores, coverage = score_modules(expression, module_genes, direction, "CRC_Atlas")
    frame = group_metadata.merge(
        scores.reset_index().rename(columns={scores.index.name: "group_id"}),
        on="group_id",
        validate="one_to_one",
    )
    effect_rows = []
    for module in module_genes:
        value = clustered_binary_effect(
            frame.loc[frame["n_cells_sampled"] >= 20],
            module,
            "carrier_group",
            "donor_id",
            "polyp_epithelial",
            "normal_epithelial",
            ["study_id", "n_cells_sampled"],
        )
        effect_rows.append({"module": module, **value})
    return frame, coverage, pd.DataFrame(effect_rows)


def spatial_context(
    module_genes: dict[str, list[str]],
    direction: dict[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    wanted = sorted(set().union(*map(set, module_genes.values())))
    section_rows: list[dict[str, object]] = []
    coverage_parts: list[pd.DataFrame] = []
    for sample_dir in sorted(
        path
        for path in closure.SPATIAL_DIR.iterdir()
        if path.is_dir() and (path / "filtered_feature_bc_matrix.h5").exists()
    ):
        sample_id = sample_dir.name
        barcodes, genes, matrix = closure.read_10x_h5(
            sample_dir / "filtered_feature_bc_matrix.h5"
        )
        gene_index = {gene: index for index, gene in enumerate(genes)}
        present = [gene for gene in wanted if gene in gene_index]
        total = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
        selected = matrix[[gene_index[gene] for gene in present], :].toarray().astype(float)
        expression = pd.DataFrame(
            np.log2(selected / np.where(total > 0, total, np.nan) * 1_000_000 + 1).T,
            index=barcodes,
            columns=present,
        )
        expression.index.name = "barcode"
        scores, coverage = score_modules(
            expression, module_genes, direction, f"spatial_{sample_id}"
        )
        coverage_parts.append(coverage)
        annotation = pd.read_csv(
            closure.SPATIAL_ANNOT_DIR / f"Pathologist_Annotations_{sample_id}.csv"
        )
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
        scored = scores.reset_index().merge(
            annotation[["barcode", "pathology_group"]],
            on="barcode",
            how="inner",
        )
        for module in module_genes:
            tumour = scored.loc[
                scored["pathology_group"].isin(["tumor", "tumor_stroma"]), module
            ].dropna()
            normal = scored.loc[
                scored["pathology_group"].eq("non_neoplastic_epithelium"), module
            ].dropna()
            if tumour.empty or normal.empty:
                continue
            section_rows.append(
                {
                    "module": module,
                    "sample_id": sample_id,
                    "n_tumour_spots": len(tumour),
                    "n_nonneoplastic_epithelial_spots": len(normal),
                    "direction_oriented_tumour_minus_normal_median": float(
                        tumour.median() - normal.median()
                    ),
                }
            )
    sections = pd.DataFrame(section_rows)
    summary_rows = []
    for module in module_genes:
        values = sections.loc[
            sections["module"].eq(module),
            "direction_oriented_tumour_minus_normal_median",
        ].to_numpy(float)
        test = (
            stats.wilcoxon(values, zero_method="wilcox", alternative="two-sided")
            if len(values) >= 4 and (values != 0).any()
            else None
        )
        summary_rows.append(
            {
                "module": module,
                "n_informative_sections": len(values),
                "median_direction_oriented_effect": float(np.median(values))
                if len(values)
                else np.nan,
                "positive_section_fraction": float((values > 0).mean())
                if len(values)
                else np.nan,
                "p_paired_wilcoxon": float(test.pvalue) if test is not None else np.nan,
            }
        )
    return (
        sections,
        pd.concat(coverage_parts, ignore_index=True),
        pd.DataFrame(summary_rows),
    )


def protein_context(
    modules: pd.DataFrame,
    retained: list[str],
    direction: dict[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    inventory = modules.loc[modules["module"].isin(retained), ["gene", "module"]].copy()
    inventory["expected_direction"] = inventory["module"].map(direction)
    pxd2 = protein.pxd002137_tests(inventory[["gene", "expected_direction"]])
    pxd17 = protein.pxd017269_detectability(
        inventory[["gene", "expected_direction"]]
    ).rename(
        columns={
            "protein_group_present": "pxd017269_protein_group_present",
            "detection_fraction": "pxd017269_detection_fraction",
        }
    )
    pxd46 = protein.pxd046999_presence(
        inventory[["gene", "expected_direction"]]
    ).rename(
        columns={"detected_in_nine_patient_dvp": "pxd046999_detected"}
    )
    evidence = (
        inventory.merge(pxd2, on="gene", how="left")
        .merge(pxd17, on="gene", how="left")
        .merge(pxd46, on="gene", how="left")
    )
    evidence["directionally_concordant_pxd002137"] = (
        evidence["protein_gene_present"].fillna(False)
        & (np.sign(evidence["age_sex_adjusted_log2_effect"]) == evidence["expected_direction"])
        & (evidence["age_sex_adjusted_p"] <= 0.10)
    )
    summary = (
        evidence.groupby("module", observed=True)
        .agg(
            n_module_genes=("gene", "size"),
            n_pxd002137_present=("protein_gene_present", "sum"),
            n_pxd002137_directionally_concordant=(
                "directionally_concordant_pxd002137",
                "sum",
            ),
            n_pxd017269_detected=("pxd017269_protein_group_present", "sum"),
            n_pxd046999_detected=("pxd046999_detected", "sum"),
        )
        .reset_index()
    )
    return evidence, summary


def main() -> None:
    required = [MODULE_PATH, VALIDATION_PATH, PAIRED_ATAC_PATH, CONTRACT_PATH]
    if not all(path.exists() for path in required):
        raise RuntimeError("At least one orthogonal-context input is missing")
    if sha256(CONTRACT_PATH) != EXPECTED_CONTRACT_HASH:
        raise RuntimeError("The frozen orthogonal-context addendum changed")
    OUT.mkdir(parents=True, exist_ok=True)
    modules = pd.read_csv(MODULE_PATH, sep="\t")
    validation = pd.read_csv(VALIDATION_PATH, sep="\t")
    if ROUTE_COLUMN not in validation.columns:
        raise RuntimeError(f"Routing column is missing: {ROUTE_COLUMN}")
    route_mask = as_bool(validation[ROUTE_COLUMN])
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
        (OUT / "module_orthogonal_context_manifest.json").write_text(
            json.dumps(
                {
                    "analysis": "state_aware_module_orthogonal_context_v1",
                    "created_utc": pd.Timestamp.utcnow().isoformat(),
                    "retained_modules": [],
                    "projection_performed": False,
                    "contract_sha256": sha256(CONTRACT_PATH),
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        print("No retained modules entered orthogonal-context projection")
        return

    module_genes = {
        module: modules.loc[modules["module"].eq(module), "gene"].astype(str).tolist()
        for module in retained
    }
    direction = {
        row.module: 1 if row.heldout_direction == "Up" else -1
        for row in validation.loc[validation["module"].isin(retained)].itertuples()
    }

    becker_scores, becker_coverage, becker_effects, rna_atac = becker_context(
        module_genes, direction
    )
    atlas_scores, atlas_coverage, atlas_effects = atlas_context(module_genes, direction)
    spatial_sections, spatial_coverage, spatial_summary = spatial_context(
        module_genes, direction
    )
    protein_evidence, protein_summary = protein_context(modules, retained, direction)

    becker_scores.to_csv(OUT / "becker_module_scores.tsv", sep="\t", index=False)
    becker_coverage.to_csv(OUT / "becker_module_coverage.tsv", sep="\t", index=False)
    becker_effects.to_csv(OUT / "becker_module_effects.tsv", sep="\t", index=False)
    rna_atac.to_csv(OUT / "module_rna_atac_correlations.tsv", sep="\t", index=False)
    atlas_scores.to_csv(OUT / "atlas_module_scores.tsv", sep="\t", index=False)
    atlas_coverage.to_csv(OUT / "atlas_module_coverage.tsv", sep="\t", index=False)
    atlas_effects.to_csv(OUT / "atlas_module_effects.tsv", sep="\t", index=False)
    spatial_sections.to_csv(OUT / "spatial_module_section_effects.tsv", sep="\t", index=False)
    spatial_coverage.to_csv(OUT / "spatial_module_coverage.tsv", sep="\t", index=False)
    spatial_summary.to_csv(OUT / "spatial_module_summary.tsv", sep="\t", index=False)
    protein_evidence.to_csv(OUT / "module_protein_evidence.tsv", sep="\t", index=False)
    protein_summary.to_csv(OUT / "module_protein_summary.tsv", sep="\t", index=False)

    context_summary = (
        pd.DataFrame({"module": retained})
        .merge(becker_effects.add_prefix("becker_").rename(columns={"becker_module": "module"}), on="module", how="left")
        .merge(rna_atac.add_prefix("rna_atac_").rename(columns={"rna_atac_module": "module"}), on="module", how="left")
        .merge(atlas_effects.add_prefix("atlas_").rename(columns={"atlas_module": "module"}), on="module", how="left")
        .merge(spatial_summary.add_prefix("spatial_").rename(columns={"spatial_module": "module"}), on="module", how="left")
        .merge(protein_summary.add_prefix("protein_").rename(columns={"protein_module": "module"}), on="module", how="left")
    )
    context_summary.to_csv(OUT / "module_orthogonal_context_summary.tsv", sep="\t", index=False)
    validation = validation.merge(context_summary, on="module", how="left")
    validation.to_csv(VALIDATION_PATH, sep="\t", index=False)

    manifest = {
        "analysis": "state_aware_module_orthogonal_context_v1",
        "created_utc": pd.Timestamp.utcnow().isoformat(),
        "random_seed": SEED,
        "route_column": ROUTE_COLUMN,
        "descriptive_downstream": DESCRIPTIVE_DOWNSTREAM,
        "retained_modules": retained,
        "interpretation": "descriptive context; cannot rescue a failed module gate",
        "input_sha256": {
            "consensus_modules": sha256(MODULE_PATH),
            "context_contract": sha256(CONTRACT_PATH),
            "paired_rna_atac": sha256(PAIRED_ATAC_PATH),
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "statsmodels": __import__("statsmodels").__version__,
            "h5py": h5py.__version__,
        },
    }
    (OUT / "module_orthogonal_context_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Orthogonal module context complete for {len(retained)} modules")


if __name__ == "__main__":
    main()
