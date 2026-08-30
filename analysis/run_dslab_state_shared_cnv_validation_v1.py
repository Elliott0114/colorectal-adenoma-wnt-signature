#!/usr/bin/env python3
"""Validate the frozen state-shared programme against CNV composition in DSLab.

Private identifiers and cell barcodes are used only in memory for linkage.
Every written result is patient-tokenised and aggregate.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

import decompose_dslab_cnv_composition_v1 as decomposition
import run_dslab_identity_validation_v1 as dslab


ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = ROOT / "results" / "state_aware_program_v1"
OUT = RESULT_ROOT / "dslab_cnv_validation"
COMMON_FILE = RESULT_ROOT / "common_effects" / "cross_state_common_effects.tsv.gz"
PANEL_FILE = RESULT_ROOT / "panel_derivation" / "compact_state_shared_panel_frozen.tsv"
ADDENDUM_FILE = (
    ROOT
    / "analysis"
    / "contracts"
    / "state_aware_dslab_cnv_validation_addendum_v1_2026-08-29.md"
)
CLINICAL_FILE = (
    ROOT
    / "results"
    / "dslab_identity_validation_v1"
    / "pathology_validation"
    / "dslab_deidentified_clinical_groups.tsv"
)
BOOTSTRAPS = 5000
SEED = 20260829
MIN_CNV_COMPARTMENT_CELLS = 50

EXPECTED_HASHES = {
    "common_effects": "a1ac4b7b67ac279782e04e386971d7463e169cbdbc24d7da4c314e34a4f3e946",
    "compact_panel": "c5997e572342a72da8441df312fba4e3461cacfa5e30d0d8590a3f23ae3d96f0",
    "implementation_addendum": "433bc2867f738beef2fb2ace5a00e8de936490ca9e9291bd361b75c46960c2e2",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_frozen_gene_sets() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, list[str]]]]:
    common = pd.read_csv(COMMON_FILE, sep="\t")
    strict_flag = common["strict_state_shared"]
    if strict_flag.dtype != bool:
        strict_flag = strict_flag.astype(str).str.lower().eq("true")
    strict = common.loc[strict_flag, ["gene", "shared_direction"]]
    strict = strict.rename(columns={"shared_direction": "arm"})
    panel = pd.read_csv(PANEL_FILE, sep="\t", usecols=["gene", "arm", "pair_step"])
    if len(strict) != 1843 or strict["arm"].value_counts().to_dict() != {"down": 959, "up": 884}:
        raise RuntimeError("Frozen strict programme does not match 1,843 genes (884 up, 959 down)")
    if len(panel) != 8 or panel["arm"].value_counts().to_dict() != {"up": 4, "down": 4}:
        raise RuntimeError("Frozen compact readout does not contain four genes per arm")
    if not set(panel["gene"]).issubset(set(strict["gene"])):
        raise RuntimeError("Compact readout is not a subset of the strict programme")
    definitions = {
        "state_shared_1843": {
            arm: sorted(strict.loc[strict["arm"].eq(arm), "gene"].astype(str))
            for arm in ("up", "down")
        },
        "compact_8": {
            arm: panel.loc[panel["arm"].eq(arm)].sort_values("pair_step")["gene"].astype(str).tolist()
            for arm in ("up", "down")
        },
    }
    return strict, panel, definitions


def bootstrap_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    observed = spearmanr(x, y)
    rng = np.random.default_rng(SEED)
    estimates: list[float] = []
    for _ in range(BOOTSTRAPS):
        index = rng.integers(0, len(x), size=len(x))
        estimate = spearmanr(x[index], y[index]).statistic
        if np.isfinite(estimate):
            estimates.append(float(estimate))
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(observed.statistic), float(low), float(high), float(observed.pvalue)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    actual_hashes = {
        "common_effects": sha256(COMMON_FILE),
        "compact_panel": sha256(PANEL_FILE),
        "implementation_addendum": sha256(ADDENDUM_FILE),
    }
    if actual_hashes != EXPECTED_HASHES:
        raise RuntimeError("A frozen state-aware input changed before DSLab validation")

    strict, panel, definitions = load_frozen_gene_sets()
    all_genes = definitions["state_shared_1843"]["up"] + definitions["state_shared_1843"]["down"]
    gene_to_arm = strict.set_index("gene")["arm"].to_dict()
    panel_genes = set(panel["gene"])

    token_map = json.loads(dslab.MATRIX_MAP.read_text(encoding="utf-8"))
    matrix_tokens = sorted(token for token in token_map if token.startswith("M"))
    matrix_dirs = {token: Path(token_map[token]) for token in matrix_tokens}
    patient_tokens = {token: f"P{index:02d}" for index, token in enumerate(matrix_tokens, start=1)}
    cnv_path = dslab.find_cnv_table()
    matrix_patient, raw_barcodes = dslab.resolve_private_linkage(
        matrix_dirs,
        cnv_path,
    )
    cnv_metadata = dslab.load_private_cnv_metadata(
        cnv_path,
        matrix_patient,
        raw_barcodes,
    )

    quantile_rows: list[dict[str, object]] = []
    detection_rows: list[dict[str, object]] = []
    sample_audit_rows: list[dict[str, object]] = []

    for token in matrix_tokens:
        patient_token = patient_tokens[token]
        features = dslab.read_features(matrix_dirs[token])
        matrix = dslab.read_matrix(matrix_dirs[token])
        barcodes = raw_barcodes[token]
        if matrix.shape != (len(features), len(barcodes)):
            raise RuntimeError(f"Dimension mismatch for opaque matrix {token}")

        metadata = cnv_metadata[token]
        cell_columns: list[int] = []
        cnv_values: list[float] = []
        cnv_thresholds: list[float] = []
        for column, barcode in enumerate(barcodes):
            record = metadata.get(barcode)
            if record is None or record["annotation"] != "Polyp_Epi":
                continue
            cell_columns.append(column)
            cnv_values.append(float(record["cnv_fraction"]))
            cnv_thresholds.append(float(record["cnv_threshold"]))
        if not cell_columns:
            raise RuntimeError(f"No source-confirmed epithelial cells for opaque matrix {token}")

        selected_matrix = matrix[:, np.asarray(cell_columns, dtype=int)]
        library_sizes = np.asarray(selected_matrix.sum(axis=0)).ravel().astype(np.float64)
        aggregation, present_genes, duplicate_counts = dslab.build_gene_aggregation(
            features,
            all_genes,
        )
        missing_panel = sorted(panel_genes.difference(present_genes))
        if missing_panel:
            raise RuntimeError(
                f"Opaque matrix {token} lacks frozen compact genes: {','.join(missing_panel)}"
            )
        gene_counts_sparse = aggregation @ selected_matrix
        present_index = np.asarray(
            [index for index, gene in enumerate(all_genes) if gene in set(present_genes)],
            dtype=int,
        )
        measurable_genes = [all_genes[index] for index in present_index]
        gene_counts = gene_counts_sparse[present_index, :].toarray().astype(np.float32)
        with np.errstate(divide="ignore", invalid="ignore"):
            log_expression = np.log1p(
                gene_counts
                / np.maximum(library_sizes, 1.0)[None, :]
                * np.float32(10_000.0)
            )

        full_cell_scores = dslab.cell_rank_scores(
            log_expression,
            measurable_genes,
            definitions["state_shared_1843"],
        )
        compact_cell_scores = dslab.cell_rank_scores(
            log_expression,
            measurable_genes,
            definitions["compact_8"],
        )
        cnv_values_array = np.asarray(cnv_values, dtype=float)
        cnv_threshold_array = np.asarray(cnv_thresholds, dtype=float)
        numeric_cnv = np.isfinite(cnv_values_array) & np.isfinite(cnv_threshold_array)
        high = numeric_cnv & (cnv_values_array >= cnv_threshold_array)
        low = numeric_cnv & ~high
        compartments = {
            "all_epithelium": np.ones(len(cell_columns), dtype=bool),
            "cnv_low": low,
            "cnv_high": high,
        }
        for compartment, mask in compartments.items():
            if mask.sum() == 0:
                continue
            quantile_rows.append(
                dslab.summarise_quantiles(
                    patient_token,
                    compartment,
                    full_cell_scores[mask],
                    "state_shared_1843",
                )
            )
            quantile_rows.append(
                dslab.summarise_quantiles(
                    patient_token,
                    compartment,
                    compact_cell_scores[mask],
                    "compact_8",
                )
            )

        measurable_lookup = {gene: index for index, gene in enumerate(measurable_genes)}
        for gene in sorted(panel_genes):
            index = measurable_lookup[gene]
            detection_rows.append(
                {
                    "patient_token": patient_token,
                    "gene": gene,
                    "arm": gene_to_arm[gene],
                    "detection_fraction": float((gene_counts[index, :] > 0).mean()),
                    "feature_rows_aggregated": duplicate_counts[gene],
                }
            )
        sample_audit_rows.append(
            {
                "patient_token": patient_token,
                "raw_cells": len(barcodes),
                "source_metadata_linked_cells": len(metadata),
                "source_confirmed_epithelial_cells": len(cell_columns),
                "cnv_low_cells": int(low.sum()),
                "cnv_high_cells": int(high.sum()),
                "cnv_missing_cells": int((~numeric_cnv).sum()),
                "strict_genes_measurable": len(measurable_genes),
                "strict_up_genes_measurable": sum(
                    gene_to_arm[gene] == "up" for gene in measurable_genes
                ),
                "strict_down_genes_measurable": sum(
                    gene_to_arm[gene] == "down" for gene in measurable_genes
                ),
                "compact_genes_measurable": sum(gene in measurable_genes for gene in panel_genes),
            }
        )
        del matrix, selected_matrix, gene_counts_sparse, gene_counts, log_expression
        print(
            f"processed\t{patient_token}\t{len(cell_columns)} epithelial cells\t"
            f"{len(measurable_genes)}/1843 genes",
            flush=True,
        )

    quantiles = pd.DataFrame(quantile_rows)
    sample_audit = pd.DataFrame(sample_audit_rows)
    detection = pd.DataFrame(detection_rows)
    if not (
        (sample_audit["cnv_low_cells"] >= MIN_CNV_COMPARTMENT_CELLS)
        & (sample_audit["cnv_high_cells"] >= MIN_CNV_COMPARTMENT_CELLS)
    ).all():
        raise RuntimeError("At least one patient fails the frozen 50-cell CNV-compartment rule")
    quantiles.to_csv(OUT / "dslab_state_shared_cell_score_quantiles.tsv", sep="\t", index=False)
    sample_audit.to_csv(OUT / "dslab_state_shared_sample_audit.tsv", sep="\t", index=False)
    detection.to_csv(OUT / "dslab_compact_8_detection.tsv", sep="\t", index=False)

    clinical = pd.read_csv(CLINICAL_FILE, sep="\t")
    populations = {
        "all_polyps": set(clinical["patient_token"]),
        "adenomatous_polyps": set(
            clinical.loc[
                clinical["adenomatous_binary"].eq("adenomatous"),
                "patient_token",
            ]
        ),
        "conventional_adenomas": set(
            clinical.loc[
                clinical["pathology_group"].eq("conventional_adenoma"),
                "patient_token",
            ]
        ),
    }

    decomposition.SEED = SEED
    decomposition.BOOTSTRAPS = BOOTSTRAPS
    decompositions = [
        decomposition.decompose_score(
            quantiles,
            sample_audit,
            score,
            population,
            patients,
        )
        for population, patients in populations.items()
        for score in ("state_shared_1843", "compact_8")
    ]
    combined = pd.concat(decompositions, ignore_index=True)
    combined.to_csv(
        OUT / "dslab_state_shared_cnv_composition_decomposition.tsv",
        sep="\t",
        index=False,
    )

    statistic_rows: list[dict[str, object]] = []
    for (population, score), frame in combined.groupby(
        ["population", "score"],
        sort=False,
    ):
        observed = decomposition.metrics(frame)
        intervals = decomposition.bootstrap_intervals(frame)
        for metric, estimate in observed.items():
            low, high = intervals.get(metric, (math.nan, math.nan))
            statistic_rows.append(
                {
                    "analysis": "cnv_composition_decomposition",
                    "population": population,
                    "score": score,
                    "metric": metric,
                    "estimate": estimate,
                    "ci_low": low,
                    "ci_high": high,
                    "p_value": math.nan,
                    "n_patients": len(frame),
                }
            )
        paired = wilcoxon(
            frame["high_mean"],
            frame["low_mean"],
            alternative="two-sided",
            method="auto",
        )
        statistic_rows.append(
            {
                "analysis": "paired_cnv_compartment",
                "population": population,
                "score": score,
                "metric": "median_high_minus_low_cell_rank_score",
                "estimate": float(np.median(frame["high_mean"] - frame["low_mean"])),
                "ci_low": math.nan,
                "ci_high": math.nan,
                "p_value": float(paired.pvalue),
                "n_patients": len(frame),
            }
        )

    for population, frame in combined.groupby("population", sort=False):
        wide = frame.pivot(
            index="patient_token",
            columns="score",
            values="common_composition_mean",
        ).dropna()
        estimate, low, high, p_value = bootstrap_spearman(
            wide["state_shared_1843"].to_numpy(dtype=float),
            wide["compact_8"].to_numpy(dtype=float),
        )
        statistic_rows.append(
            {
                "analysis": "common_composition_compact_fidelity",
                "population": population,
                "score": "compact_8_vs_state_shared_1843",
                "metric": "spearman_rho",
                "estimate": estimate,
                "ci_low": low,
                "ci_high": high,
                "p_value": p_value,
                "n_patients": len(wide),
            }
        )

    all_epithelium = quantiles.loc[quantiles["compartment"].eq("all_epithelium")]
    score_wide = all_epithelium.pivot(
        index="patient_token",
        columns="score",
        values="mean",
    ).reset_index()
    score_wide = score_wide.merge(
        clinical[["patient_token", "pathology_group", "adenomatous_binary"]],
        on="patient_token",
        validate="one_to_one",
    )
    for score in ("state_shared_1843", "compact_8"):
        conventional = score_wide.loc[
            score_wide["pathology_group"].eq("conventional_adenoma"),
            score,
        ].to_numpy(dtype=float)
        non_adenomatous = score_wide.loc[
            score_wide["adenomatous_binary"].eq("non_adenomatous"),
            score,
        ].to_numpy(dtype=float)
        test = mannwhitneyu(conventional, non_adenomatous, alternative="two-sided")
        statistic_rows.append(
            {
                "analysis": "descriptive_pathology_negative_control",
                "population": "conventional_vs_non_adenomatous_polyps",
                "score": score,
                "metric": "rank_auc",
                "estimate": float(test.statistic / (len(conventional) * len(non_adenomatous))),
                "ci_low": math.nan,
                "ci_high": math.nan,
                "p_value": float(test.pvalue),
                "n_patients": len(conventional) + len(non_adenomatous),
            }
        )

    statistics = pd.DataFrame(statistic_rows)
    statistics.to_csv(
        OUT / "dslab_state_shared_cnv_validation_statistics.tsv",
        sep="\t",
        index=False,
    )

    qc = {
        "patient_count": len(matrix_tokens),
        "epithelial_cell_count": int(sample_audit["source_confirmed_epithelial_cells"].sum()),
        "cnv_low_cell_count": int(sample_audit["cnv_low_cells"].sum()),
        "cnv_high_cell_count": int(sample_audit["cnv_high_cells"].sum()),
        "all_compartments_at_least_50_cells": bool(
            (
                (sample_audit["cnv_low_cells"] >= MIN_CNV_COMPARTMENT_CELLS)
                & (sample_audit["cnv_high_cells"] >= MIN_CNV_COMPARTMENT_CELLS)
            ).all()
        ),
        "minimum_strict_genes_measurable": int(sample_audit["strict_genes_measurable"].min()),
        "compact_8_present_in_every_sample": bool(
            sample_audit["compact_genes_measurable"].eq(8).all()
        ),
        "source_patient_ids_exported": False,
        "cell_ids_or_barcodes_exported": False,
    }
    manifest = {
        "analysis": "run_dslab_state_shared_cnv_validation_v1",
        "analysis_date": "2026-08-29",
        "frozen_input_hashes": actual_hashes,
        "restricted_linkage_map_sha256": sha256(dslab.MATRIX_MAP),
        "clinical_group_file_sha256": sha256(CLINICAL_FILE),
        "minimum_cnv_compartment_cells": MIN_CNV_COMPARTMENT_CELLS,
        "bootstrap_replicates": BOOTSTRAPS,
        "seed": SEED,
        "score": "within-cell percentile-rank, mean(up) minus mean(down), equal arm weights",
        "validation_outcomes_used_for_gene_selection": False,
        "fine_state_projection_used": False,
        "interpretation": "Composition accounting and compact-readout fidelity; not clinical validation.",
        "qc": qc,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }
    (OUT / "dslab_state_shared_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(qc, indent=2), flush=True)
    print(f"output_dir\t{OUT.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
