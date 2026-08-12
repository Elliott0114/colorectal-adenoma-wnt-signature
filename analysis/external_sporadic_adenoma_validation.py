#!/usr/bin/env python3
"""Validate the discovery-locked route in independent sporadic adenoma cohorts."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "data_sources" / "GEO_sporadic_adenoma_validation"
ROUTE_DIR = ROOT / "results" / "route_signature_locked"
OUT_DIR = ROOT / "results" / "external_sporadic_adenoma_validation"
GPL570_MAP = (
    ROOT
    / "results"
    / "geo_bulk"
    / "multicohort_recurrence"
    / "gpl570_target_probe_mapping.tsv"
)

SIGNATURE_SIZES = [10, 20, 30, 50]
PRIMARY_SIZE = 50
MIN_COMPONENT_COVERAGE = 0.80
N_BOOTSTRAPS = 5_000
SEED = 20260710
PROLIFERATION_CONTROL = [
    "MKI67",
    "TOP2A",
    "PCNA",
    "MCM2",
    "MCM5",
    "TYMS",
    "UBE2C",
    "CENPF",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def split_symbols(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text.lower() in {"nan", "na", "---"}:
        return []
    values = re.split(r"\s*(?:///|//|;|\||,)\s*", text)
    return [item.strip() for item in values if item.strip() and item.strip() != "---"]


def parse_series_matrix(
    accession: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[list[str]]]]:
    path = SOURCE_DIR / f"{accession}_series_matrix.txt.gz"
    metadata: dict[str, list[list[str]]] = {}
    table_rows: list[list[str]] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_rows.append(next(csv.reader([line], delimiter="\t")))
            elif line.startswith("!"):
                parts = next(csv.reader([line], delimiter="\t"))
                metadata.setdefault(parts[0], []).append(parts[1:])

    accessions = metadata["!Sample_geo_accession"][0]
    sample_meta = pd.DataFrame(
        {
            "sample_id": accessions,
            "sample_title": metadata["!Sample_title"][0],
            "source_name": metadata["!Sample_source_name_ch1"][0],
        }
    )
    for values in metadata.get("!Sample_characteristics_ch1", []):
        parsed = [value.split(":", 1) if ":" in value else ["characteristic", value] for value in values]
        keys = [clean_key(pair[0]) for pair in parsed]
        key = pd.Series(keys).mode().iloc[0]
        sample_meta[key] = [pair[1].strip() if len(pair) > 1 else "" for pair in parsed]

    if len(table_rows) <= 1:
        return sample_meta, pd.DataFrame(index=sample_meta["sample_id"]), metadata
    header = table_rows[0]
    if header[1:] != accessions:
        raise ValueError(f"{accession}: series matrix sample order does not match metadata")
    expression = pd.DataFrame(table_rows[1:], columns=header)
    expression = expression.rename(columns={header[0]: "feature_id"}).set_index("feature_id")
    expression = expression.apply(pd.to_numeric, errors="coerce").T
    expression.index.name = "sample_id"
    return sample_meta, expression, metadata


def gpl570_mapping() -> pd.DataFrame:
    mapping = pd.read_csv(GPL570_MAP, sep="\t")
    rows = []
    for row in mapping.itertuples(index=False):
        for gene in split_symbols(row.matched_module_genes):
            rows.append({"feature_id": str(row.probe_id), "gene": gene})
    return pd.DataFrame(rows).drop_duplicates()


def gpl6480_mapping() -> pd.DataFrame:
    path = SOURCE_DIR / "GPL6480.annot.gz"
    rows: list[list[str]] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if parts and parts[0] == "!platform_table_begin":
                in_table = True
                continue
            if parts and parts[0] == "!platform_table_end":
                break
            if in_table:
                rows.append(parts)
    table = pd.DataFrame(rows[1:], columns=rows[0])
    output = []
    for row in table[["ID", "Gene symbol"]].itertuples(index=False):
        for gene in split_symbols(row[1]):
            output.append({"feature_id": str(row[0]), "gene": gene})
    return pd.DataFrame(output).drop_duplicates()


def gpl8432_mapping() -> pd.DataFrame:
    path = SOURCE_DIR / "GPL8432_HUMANREF-8_V3_0_R1_11282963_A_WGDASL.txt.gz"
    lines: list[str] = []
    in_probes = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.rstrip("\n") == "[Probes]":
                in_probes = True
                continue
            if in_probes and line.startswith("["):
                break
            if in_probes:
                lines.append(line)
    table = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
    table["gene_for_mapping"] = table["Symbol"].fillna(table["ILMN_Gene"])
    output = []
    for row in table[["Probe_Id", "gene_for_mapping"]].itertuples(index=False):
        for gene in split_symbols(row[1]):
            output.append({"feature_id": str(row[0]), "gene": gene})
    return pd.DataFrame(output).drop_duplicates()


def choose_features_without_labels(
    expression: pd.DataFrame,
    mapping: pd.DataFrame,
    wanted_genes: set[str],
    cohort: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = mapping.loc[mapping["gene"].isin(wanted_genes)].copy()
    mapping = mapping.loc[mapping["feature_id"].isin(expression.columns)].copy()
    variances = expression.var(axis=0, skipna=True, ddof=1)
    mapping["variance_all_samples"] = mapping["feature_id"].map(variances)
    selected = (
        mapping.sort_values(
            ["gene", "variance_all_samples", "feature_id"],
            ascending=[True, False, True],
        )
        .drop_duplicates("gene")
        .copy()
    )
    selected.insert(0, "cohort", cohort)
    selected["selection_uses_group_labels"] = False
    gene_expression = expression[selected["feature_id"]].copy()
    gene_expression.columns = selected["gene"].tolist()
    return gene_expression, selected


def direct_rnaseq_expression(wanted_genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    sample_meta = []
    for path in sorted(SOURCE_DIR.glob("GSE72820_*RPKM.xlsx.gz")):
        match = re.search(r"MDA(\d+)N-", path.name)
        if not match:
            raise ValueError(f"Cannot parse patient from {path.name}")
        patient = match.group(1)
        with gzip.open(path, "rb") as handle:
            content = handle.read()
        table = pd.read_excel(io.BytesIO(content))
        table["gene"] = table["gene"].astype(str)
        table = table.loc[table["gene"].isin(wanted_genes), ["gene", "value_1", "value_2"]]
        table = table.groupby("gene", as_index=False)[["value_1", "value_2"]].median()
        for tissue, column in [("normal", "value_1"), ("adenoma", "value_2")]:
            sample_id = f"GSE72820_{patient}_{tissue}"
            row = {"sample_id": sample_id}
            row.update(dict(zip(table["gene"], np.log2(table[column].clip(lower=0) + 1), strict=True)))
            rows.append(row)
            sample_meta.append(
                {
                    "sample_id": sample_id,
                    "sample_title": sample_id,
                    "source_name": tissue,
                    "tissue_group": tissue,
                    "patient_id": patient,
                    "pathologic_grade": tissue,
                }
            )
    expression = pd.DataFrame(rows).set_index("sample_id")
    selected = pd.DataFrame(
        {
            "cohort": "GSE72820",
            "feature_id": expression.columns,
            "gene": expression.columns,
            "variance_all_samples": expression.var(axis=0, ddof=1).to_numpy(),
            "selection_uses_group_labels": False,
        }
    )
    return pd.DataFrame(sample_meta), expression, selected


def gse50114_raw_expression() -> pd.DataFrame:
    path = OUT_DIR / "GSE50114_raw_limma_log2_quantile.tsv.gz"
    if not path.exists():
        raise FileNotFoundError(
            f"Run analysis/preprocess_gse50114_raw.R before scoring: {path}"
        )
    table = pd.read_csv(path, sep="\t")
    expression = table.set_index("feature_id").T
    expression.index.name = "sample_id"
    return expression


def add_design_fields(accession: str, meta: pd.DataFrame) -> pd.DataFrame:
    out = meta.copy()
    if accession == "GSE8671":
        out["tissue_group"] = np.where(out["source_name"].str.contains("adenoma", case=False), "adenoma", "normal")
        out["patient_id"] = out["sample_title"].str.extract(r"patient #(\d+)", expand=False)
        out["pathologic_grade"] = out["tissue_group"]
    elif accession == "GSE50114":
        tissue = out["tissue"].str.lower()
        out["tissue_group"] = np.select(
            [tissue.str.contains("adenoma"), tissue.str.contains("normal")],
            ["adenoma", "normal"],
            default="other",
        )
        out["patient_id"] = out["sample_title"].str.extract(r"patient\s+(\d+)", expand=False)
        out["pathologic_grade"] = out["tissue_group"]
    elif accession == "GSE41657":
        grade = out["pathologic_grade"].str.lower()
        out["tissue_group"] = np.select(
            [grade.str.contains("normal"), grade.str.contains("dysplasia"), grade.str.contains("adenocarcinoma")],
            ["normal", "adenoma", "crc"],
            default="other",
        )
        out["patient_id"] = out["sample_title"]
        out["grade_group"] = np.select(
            [
                grade.str.contains("normal"),
                grade.str.contains("low-grade"),
                grade.str.contains("high-grade"),
                grade.str.contains("adenocarcinoma"),
            ],
            ["normal", "low_grade", "high_grade", "crc"],
            default="other",
        )
    elif accession == "GSE40362":
        source = out["source_name"].str.lower()
        out["tissue_group"] = np.select(
            [source.str.contains("adenomatous"), source.str.contains("hyperplastic"), source.str.contains("normal")],
            ["adenoma", "hyperplastic", "normal"],
            default="other",
        )
        out["patient_id"] = out["sample_title"]
        out["pathologic_grade"] = out["tissue_group"]
    else:
        raise ValueError(accession)
    if out["patient_id"].isna().any():
        raise ValueError(f"{accession}: patient ID could not be recovered for every sample")
    return out


def score_nested_signatures(
    meta: pd.DataFrame,
    expression: pd.DataFrame,
    signature: pd.DataFrame,
    cohort: str,
    normalization_groups: set[str] | None = None,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    output = meta.copy().set_index("sample_id")
    expression = expression.reindex(output.index)
    reference_expression = expression
    normalization_scope = "all_cohort_samples"
    if normalization_groups is not None:
        reference_mask = output["tissue_group"].isin(normalization_groups)
        reference_expression = expression.loc[reference_mask]
        normalization_scope = "normal_and_adenoma_only"
    mean = reference_expression.mean(axis=0, skipna=True)
    sd = reference_expression.std(axis=0, skipna=True, ddof=1).replace(0, np.nan)
    z = expression.sub(mean, axis=1).div(sd, axis=1)
    present_proliferation = [gene for gene in PROLIFERATION_CONTROL if gene in z.columns]
    proliferation_coverage = len(present_proliferation) / len(PROLIFERATION_CONTROL)
    if proliferation_coverage < MIN_COMPONENT_COVERAGE:
        raise RuntimeError(
            f"{cohort}: proliferation-control coverage {proliferation_coverage:.1%}"
        )
    output["score__proliferation_control"] = z[present_proliferation].mean(
        axis=1, skipna=True
    )
    output["n_proliferation_control_present"] = z[present_proliferation].notna().sum(
        axis=1
    )
    coverage_rows = []
    for size in SIGNATURE_SIZES:
        up = signature.loc[
            signature["signature_direction"].eq("adenoma_up")
            & signature["rank_within_direction"].le(size),
            "gene",
        ].tolist()
        down = signature.loc[
            signature["signature_direction"].eq("adenoma_down")
            & signature["rank_within_direction"].le(size),
            "gene",
        ].tolist()
        present_up = [gene for gene in up if gene in z.columns]
        present_down = [gene for gene in down if gene in z.columns]
        coverage_up = len(present_up) / size
        coverage_down = len(present_down) / size
        if min(coverage_up, coverage_down) < MIN_COMPONENT_COVERAGE:
            raise RuntimeError(
                f"{cohort} size {size}: component coverage {coverage_up:.1%}/{coverage_down:.1%}"
            )
        score = z[present_up].mean(axis=1, skipna=True) - z[present_down].mean(axis=1, skipna=True)
        available_up = z[present_up].notna().sum(axis=1)
        available_down = z[present_down].notna().sum(axis=1)
        if ((available_up / len(present_up) < MIN_COMPONENT_COVERAGE) | (available_down / len(present_down) < MIN_COMPONENT_COVERAGE)).any():
            raise RuntimeError(f"{cohort} size {size}: sample-level gene coverage below threshold")
        output[f"route_score_k{size}"] = score
        output[f"n_up_present_k{size}"] = available_up
        output[f"n_down_present_k{size}"] = available_down
        coverage_rows.append(
            {
                "cohort": cohort,
                "signature_size_per_direction": size,
                "requested_up": size,
                "present_up": len(present_up),
                "coverage_up": coverage_up,
                "requested_down": size,
                "present_down": len(present_down),
                "coverage_down": coverage_down,
                "proliferation_requested": len(PROLIFERATION_CONTROL),
                "proliferation_present": len(present_proliferation),
                "proliferation_coverage": proliferation_coverage,
                "passes_80pct_component_rule": True,
            }
        )
    output = output.reset_index()
    output.insert(0, "cohort", cohort)
    output["normalization_scope"] = normalization_scope
    output["patient_cluster_id"] = cohort + "::" + output["patient_id"].astype(str)
    return output, coverage_rows


def cohort_comparison(
    scores: pd.DataFrame,
    cohort: str,
    size: int,
    group_a: str,
    group_b: str,
    paired: bool,
    comparison: str,
) -> dict[str, object]:
    column = f"route_score_k{size}"
    subset = scores.loc[scores["tissue_group"].isin([group_a, group_b])].copy()
    a = subset.loc[subset["tissue_group"].eq(group_a), column].dropna()
    b = subset.loc[subset["tissue_group"].eq(group_b), column].dropna()
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    row: dict[str, object] = {
        "cohort": cohort,
        "comparison": comparison,
        "group_a": group_a,
        "group_b": group_b,
        "signature_size_per_direction": size,
        "paired_primary": paired,
        "n_a": len(a),
        "n_b": len(b),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "median_difference_a_minus_b": float(a.median() - b.median()),
        "auc_a_vs_b": float(mw.statistic / (len(a) * len(b))),
        "rank_biserial_a_vs_b": float(2 * mw.statistic / (len(a) * len(b)) - 1),
        "p_mannwhitney": float(mw.pvalue),
        "n_pairs": np.nan,
        "median_paired_difference": np.nan,
        "p_paired_wilcoxon": np.nan,
        "paired_positive_fraction": np.nan,
        "primary_p_value": float(mw.pvalue),
    }
    cluster_data = subset[["patient_cluster_id", "tissue_group", column]].dropna().copy()
    cluster_data["is_group_a"] = cluster_data["tissue_group"].eq(group_a).astype(float)
    design = sm.add_constant(cluster_data[["is_group_a"]], has_constant="add")
    cluster_fit = sm.OLS(cluster_data[column].astype(float), design).fit(
        cov_type="cluster",
        cov_kwds={
            "groups": cluster_data["patient_cluster_id"],
            "use_correction": True,
            "df_correction": True,
        },
        use_t=True,
    )
    row.update(
        {
            "clustered_mean_difference": float(cluster_fit.params["is_group_a"]),
            "clustered_ci_low": float(cluster_fit.conf_int().loc["is_group_a", 0]),
            "clustered_ci_high": float(cluster_fit.conf_int().loc["is_group_a", 1]),
            "p_patient_clustered_ols": float(cluster_fit.pvalues["is_group_a"]),
            "n_patient_clusters": cluster_data["patient_cluster_id"].nunique(),
        }
    )
    cluster_data["score_z_analysis_set"] = (
        cluster_data[column] - cluster_data[column].mean()
    ) / cluster_data[column].std(ddof=1)
    standardized_fit = sm.OLS(cluster_data["score_z_analysis_set"], design).fit(
        cov_type="cluster",
        cov_kwds={
            "groups": cluster_data["patient_cluster_id"],
            "use_correction": True,
            "df_correction": True,
        },
        use_t=True,
    )
    row.update(
        {
            "clustered_standardized_mean_difference": float(
                standardized_fit.params["is_group_a"]
            ),
            "clustered_standardized_ci_low": float(
                standardized_fit.conf_int().loc["is_group_a", 0]
            ),
            "clustered_standardized_ci_high": float(
                standardized_fit.conf_int().loc["is_group_a", 1]
            ),
            "p_patient_clustered_standardized_ols": float(
                standardized_fit.pvalues["is_group_a"]
            ),
        }
    )
    if not paired:
        row["primary_p_value"] = float(cluster_fit.pvalues["is_group_a"])
    if paired:
        wide = subset.pivot_table(index="patient_id", columns="tissue_group", values=column, aggfunc="median")
        pair = wide[[group_a, group_b]].dropna()
        difference = pair[group_a] - pair[group_b]
        wilcoxon = stats.wilcoxon(difference, zero_method="wilcox", alternative="two-sided")
        row.update(
            {
                "n_pairs": len(pair),
                "median_paired_difference": float(difference.median()),
                "p_paired_wilcoxon": float(wilcoxon.pvalue),
                "paired_positive_fraction": float((difference > 0).mean()),
                "primary_p_value": float(wilcoxon.pvalue),
            }
        )
    return row


def bootstrap_primary(
    scores: pd.DataFrame,
    cohort: str,
    paired: bool,
    rng: np.random.Generator,
) -> dict[str, object]:
    column = f"route_score_k{PRIMARY_SIZE}"
    subset = scores.loc[scores["tissue_group"].isin(["adenoma", "normal"])].copy()
    auc_values = []
    delta_values = []
    if paired:
        wide = subset.pivot_table(index="patient_id", columns="tissue_group", values=column, aggfunc="median")
        wide = wide[["adenoma", "normal"]].dropna()
        for _ in range(N_BOOTSTRAPS):
            indices = rng.integers(0, len(wide), size=len(wide))
            sampled = wide.iloc[indices]
            adenoma = sampled["adenoma"].to_numpy()
            normal = sampled["normal"].to_numpy()
            auc_values.append(stats.mannwhitneyu(adenoma, normal).statistic / (len(adenoma) * len(normal)))
            delta_values.append(float(np.median(adenoma - normal)))
    else:
        patients = subset["patient_id"].unique()
        blocks = {
            patient: part for patient, part in subset.groupby("patient_id", sort=False)
        }
        for _ in range(N_BOOTSTRAPS):
            sampled_patients = rng.choice(patients, size=len(patients), replace=True)
            sampled = pd.concat([blocks[patient] for patient in sampled_patients], ignore_index=True)
            a = sampled.loc[sampled["tissue_group"].eq("adenoma"), column].to_numpy()
            n = sampled.loc[sampled["tissue_group"].eq("normal"), column].to_numpy()
            if len(a) == 0 or len(n) == 0:
                continue
            auc_values.append(stats.mannwhitneyu(a, n).statistic / (len(a) * len(n)))
            delta_values.append(float(np.median(a) - np.median(n)))
    return {
        "cohort": cohort,
        "signature_size_per_direction": PRIMARY_SIZE,
        "bootstrap_unit": "patient_pair" if paired else "patient_cluster",
        "n_bootstraps": N_BOOTSTRAPS,
        "seed": SEED,
        "auc_bootstrap_median": float(np.median(auc_values)),
        "auc_ci_low": float(np.quantile(auc_values, 0.025)),
        "auc_ci_high": float(np.quantile(auc_values, 0.975)),
        "median_delta_bootstrap_median": float(np.median(delta_values)),
        "median_delta_ci_low": float(np.quantile(delta_values, 0.025)),
        "median_delta_ci_high": float(np.quantile(delta_values, 0.975)),
    }


def one_stage_model(scores: pd.DataFrame, size: int, excluded_cohort: str = "__NONE__") -> dict[str, object]:
    column = f"route_score_k{size}"
    data = scores.loc[scores["tissue_group"].isin(["normal", "adenoma"])].copy()
    if excluded_cohort != "__NONE__":
        data = data.loc[data["cohort"].ne(excluded_cohort)].copy()
    data["score_sd_within_cohort"] = data.groupby("cohort")[column].transform(
        lambda values: (values - values.mean()) / values.std(ddof=1)
    )
    data["is_adenoma"] = data["tissue_group"].eq("adenoma").astype(float)
    cohort_design = pd.get_dummies(data["cohort"], prefix="cohort", drop_first=True, dtype=float)
    design = pd.concat([data[["is_adenoma"]].reset_index(drop=True), cohort_design.reset_index(drop=True)], axis=1)
    design = sm.add_constant(design, has_constant="add")
    fit = sm.OLS(data["score_sd_within_cohort"].reset_index(drop=True), design).fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data["patient_cluster_id"].reset_index(drop=True),
            "use_correction": True,
            "df_correction": True,
        },
        use_t=True,
    )
    return {
        "signature_size_per_direction": size,
        "excluded_cohort": excluded_cohort,
        "n_samples": len(data),
        "n_patient_clusters": data["patient_cluster_id"].nunique(),
        "n_cohorts": data["cohort"].nunique(),
        "adenoma_coef_sd": float(fit.params["is_adenoma"]),
        "se_patient_cluster": float(fit.bse["is_adenoma"]),
        "ci_low": float(fit.conf_int().loc["is_adenoma", 0]),
        "ci_high": float(fit.conf_int().loc["is_adenoma", 1]),
        "p_value": float(fit.pvalues["is_adenoma"]),
        "r_squared": float(fit.rsquared),
    }


def adjusted_cohort_model(
    scores: pd.DataFrame,
    cohort: str,
    age_column: str,
    sex_column: str,
) -> dict[str, object]:
    column = f"route_score_k{PRIMARY_SIZE}"
    data = scores.loc[
        scores["cohort"].eq(cohort)
        & scores["tissue_group"].isin(["normal", "adenoma"]),
        [column, "tissue_group", "patient_cluster_id", age_column, sex_column],
    ].dropna().copy()
    data["score_sd"] = (data[column] - data[column].mean()) / data[column].std(ddof=1)
    data["is_adenoma"] = data["tissue_group"].eq("adenoma").astype(float)
    data["age_centered"] = pd.to_numeric(data[age_column]) - pd.to_numeric(data[age_column]).mean()
    data["is_male"] = data[sex_column].str.lower().eq("male").astype(float)
    design = sm.add_constant(data[["is_adenoma", "age_centered", "is_male"]], has_constant="add")
    fit = sm.OLS(data["score_sd"], design).fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data["patient_cluster_id"],
            "use_correction": True,
            "df_correction": True,
        },
        use_t=True,
    )
    return {
        "cohort": cohort,
        "signature_size_per_direction": PRIMARY_SIZE,
        "covariates": "age + sex",
        "n_samples": len(data),
        "n_patient_clusters": data["patient_cluster_id"].nunique(),
        "adjusted_adenoma_coef_sd": float(fit.params["is_adenoma"]),
        "adjusted_ci_low": float(fit.conf_int().loc["is_adenoma", 0]),
        "adjusted_ci_high": float(fit.conf_int().loc["is_adenoma", 1]),
        "adjusted_p_value": float(fit.pvalues["is_adenoma"]),
    }


def proliferation_adjusted_model(
    scores: pd.DataFrame,
    excluded_cohort: str = "__NONE__",
) -> dict[str, object]:
    column = f"route_score_k{PRIMARY_SIZE}"
    data = scores.loc[scores["tissue_group"].isin(["normal", "adenoma"])].copy()
    if excluded_cohort != "__NONE__":
        data = data.loc[data["cohort"].ne(excluded_cohort)].copy()
    data["score_sd_within_cohort"] = data.groupby("cohort")[column].transform(
        lambda values: (values - values.mean()) / values.std(ddof=1)
    )
    data["proliferation_sd_within_cohort"] = data.groupby("cohort")[
        "score__proliferation_control"
    ].transform(lambda values: (values - values.mean()) / values.std(ddof=1))
    data["is_adenoma"] = data["tissue_group"].eq("adenoma").astype(float)
    cohort_design = pd.get_dummies(
        data["cohort"], prefix="cohort", drop_first=True, dtype=float
    )
    design = pd.concat(
        [
            data[["is_adenoma", "proliferation_sd_within_cohort"]].reset_index(
                drop=True
            ),
            cohort_design.reset_index(drop=True),
        ],
        axis=1,
    )
    design = sm.add_constant(design, has_constant="add")
    fit = sm.OLS(
        data["score_sd_within_cohort"].reset_index(drop=True), design
    ).fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data["patient_cluster_id"].reset_index(drop=True),
            "use_correction": True,
            "df_correction": True,
        },
        use_t=True,
    )
    return {
        "signature_size_per_direction": PRIMARY_SIZE,
        "excluded_cohort": excluded_cohort,
        "n_samples": len(data),
        "n_patient_clusters": data["patient_cluster_id"].nunique(),
        "n_cohorts": data["cohort"].nunique(),
        "covariates": "cohort fixed effects + proliferation control",
        "adenoma_coef_sd": float(fit.params["is_adenoma"]),
        "se_patient_cluster": float(fit.bse["is_adenoma"]),
        "ci_low": float(fit.conf_int().loc["is_adenoma", 0]),
        "ci_high": float(fit.conf_int().loc["is_adenoma", 1]),
        "p_value": float(fit.pvalues["is_adenoma"]),
        "proliferation_coef_sd": float(
            fit.params["proliferation_sd_within_cohort"]
        ),
        "proliferation_p_value": float(
            fit.pvalues["proliferation_sd_within_cohort"]
        ),
        "r_squared": float(fit.rsquared),
        "condition_number": float(np.linalg.cond(design.to_numpy(dtype=float))),
    }


def cohort_proliferation_adjusted_model(
    scores: pd.DataFrame,
    cohort: str,
) -> dict[str, object]:
    column = f"route_score_k{PRIMARY_SIZE}"
    data = scores.loc[
        scores["cohort"].eq(cohort)
        & scores["tissue_group"].isin(["normal", "adenoma"])
    ].copy()
    data["score_sd"] = (data[column] - data[column].mean()) / data[column].std(ddof=1)
    data["proliferation_sd"] = (
        data["score__proliferation_control"]
        - data["score__proliferation_control"].mean()
    ) / data["score__proliferation_control"].std(ddof=1)
    data["is_adenoma"] = data["tissue_group"].eq("adenoma").astype(float)
    design = sm.add_constant(
        data[["is_adenoma", "proliferation_sd"]], has_constant="add"
    )
    fit = sm.OLS(data["score_sd"], design).fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data["patient_cluster_id"],
            "use_correction": True,
            "df_correction": True,
        },
        use_t=True,
    )
    return {
        "cohort": cohort,
        "signature_size_per_direction": PRIMARY_SIZE,
        "n_samples": len(data),
        "n_patient_clusters": data["patient_cluster_id"].nunique(),
        "proliferation_genes_present": int(
            data["n_proliferation_control_present"].min()
        ),
        "route_proliferation_spearman_rho": float(
            stats.spearmanr(data[column], data["score__proliferation_control"]).statistic
        ),
        "adjusted_adenoma_coef_sd": float(fit.params["is_adenoma"]),
        "adjusted_se_patient_cluster": float(fit.bse["is_adenoma"]),
        "adjusted_ci_low": float(fit.conf_int().loc["is_adenoma", 0]),
        "adjusted_ci_high": float(fit.conf_int().loc["is_adenoma", 1]),
        "adjusted_p_value": float(fit.pvalues["is_adenoma"]),
        "proliferation_coef_sd": float(fit.params["proliferation_sd"]),
        "proliferation_p_value": float(fit.pvalues["proliferation_sd"]),
        "condition_number": float(np.linalg.cond(design.to_numpy(dtype=float))),
    }


def read_all_cohorts(
    signature: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    overlap = set(signature["gene"]) & set(PROLIFERATION_CONTROL)
    if overlap:
        raise RuntimeError(
            "Proliferation-control genes overlap the locked route: "
            + ", ".join(sorted(overlap))
        )
    wanted = set(signature["gene"]) | set(PROLIFERATION_CONTROL)
    mappings = {
        "GSE8671": gpl570_mapping(),
        "GSE50114": gpl6480_mapping(),
        "GSE41657": gpl6480_mapping(),
        "GSE40362": gpl8432_mapping(),
    }
    log2_transform = {"GSE8671": True, "GSE50114": False, "GSE41657": False, "GSE40362": False}
    all_scores = []
    all_primary_set_scores = []
    all_selected = []
    coverage_rows = []
    audit_rows = []

    for accession in ["GSE8671", "GSE50114", "GSE41657", "GSE40362"]:
        meta, expression, metadata = parse_series_matrix(accession)
        meta = add_design_fields(accession, meta)
        if accession == "GSE50114":
            expression = gse50114_raw_expression()
        if log2_transform[accession]:
            expression = np.log2(expression.clip(lower=0) + 1)
        gene_expression, selected = choose_features_without_labels(
            expression,
            mappings[accession],
            wanted,
            accession,
        )
        scores, coverage = score_nested_signatures(meta, gene_expression, signature, accession)
        primary_set_scores, _ = score_nested_signatures(
            meta,
            gene_expression,
            signature,
            accession,
            normalization_groups={"normal", "adenoma"},
        )
        all_scores.append(scores)
        all_primary_set_scores.append(primary_set_scores)
        all_selected.append(selected)
        coverage_rows.extend(coverage)
        counts = scores["tissue_group"].value_counts()
        paired_patients = (
            scores.loc[scores["tissue_group"].isin(["normal", "adenoma"])]
            .groupby("patient_id")["tissue_group"]
            .nunique()
            .eq(2)
            .sum()
        )
        audit_rows.append(
            {
                "cohort": accession,
                "platform": metadata["!Series_platform_id"][0][0],
                "n_samples": len(scores),
                "n_normal": int(counts.get("normal", 0)),
                "n_adenoma": int(counts.get("adenoma", 0)),
                "n_hyperplastic": int(counts.get("hyperplastic", 0)),
                "n_crc": int(counts.get("crc", 0)),
                "n_patient_ids": scores["patient_id"].nunique(),
                "n_complete_normal_adenoma_pairs": int(paired_patients),
                "processed_expression_transform": (
                    "limma normexp(offset=50) + quantile from complete raw arrays"
                    if accession == "GSE50114"
                    else ("log2(x+1)" if log2_transform[accession] else "as deposited")
                ),
                "feature_selection": "highest variance probe per locked gene; group labels unused",
                "eligible_for_locked_score_validation": True,
                "exclusion_reason": "",
            }
        )

    meta, expression, selected = direct_rnaseq_expression(wanted)
    scores, coverage = score_nested_signatures(meta, expression, signature, "GSE72820")
    primary_set_scores, _ = score_nested_signatures(
        meta,
        expression,
        signature,
        "GSE72820",
        normalization_groups={"normal", "adenoma"},
    )
    all_scores.append(scores)
    all_primary_set_scores.append(primary_set_scores)
    all_selected.append(selected)
    coverage_rows.extend(coverage)
    audit_rows.append(
        {
            "cohort": "GSE72820",
            "platform": "GPL11154 RNA-seq",
            "n_samples": len(scores),
            "n_normal": int(scores["tissue_group"].eq("normal").sum()),
            "n_adenoma": int(scores["tissue_group"].eq("adenoma").sum()),
            "n_hyperplastic": 0,
            "n_crc": 0,
            "n_patient_ids": scores["patient_id"].nunique(),
            "n_complete_normal_adenoma_pairs": scores["patient_id"].nunique(),
            "processed_expression_transform": "log2(RPKM+1)",
            "feature_selection": "deposited gene symbols; group labels unused",
            "eligible_for_locked_score_validation": True,
            "exclusion_reason": "",
        }
    )
    return (
        pd.concat(all_scores, ignore_index=True),
        pd.concat(all_selected, ignore_index=True),
        pd.DataFrame(coverage_rows),
        pd.DataFrame(audit_rows),
        pd.concat(all_primary_set_scores, ignore_index=True),
    )


def source_manifest() -> pd.DataFrame:
    urls = {
        "GSE8671_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8671/matrix/GSE8671_series_matrix.txt.gz",
        "GSE50114_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50114/matrix/GSE50114_series_matrix.txt.gz",
        "GSE50114_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50114/suppl/GSE50114_RAW.tar",
        "GSE41657_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41657/matrix/GSE41657_series_matrix.txt.gz",
        "GSE40362_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40362/matrix/GSE40362_series_matrix.txt.gz",
        "GSE40362_non-normalized.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40362/suppl/GSE40362_non-normalized.txt.gz",
        "GPL6480.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz",
        "GPL8432_HUMANREF-8_V3_0_R1_11282963_A_WGDASL.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8432/suppl/GPL8432_HUMANREF-8_V3_0_R1_11282963_A_WGDASL.txt.gz",
    }
    rows = []
    for path in sorted(SOURCE_DIR.iterdir()):
        if not path.is_file() or path.name == "SHA256SUMS.txt":
            continue
        url = urls.get(path.name)
        if path.name.startswith("GSE72820_") and path.name.endswith("RPKM.xlsx.gz"):
            url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72820/suppl/{path.name}"
        elif path.name == "GSE72820_series_matrix.txt.gz":
            url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72820/matrix/GSE72820_series_matrix.txt.gz"
        rows.append(
            {
                "file": str(path.relative_to(ROOT)),
                "source_url": url,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    signature = pd.read_csv(ROUTE_DIR / "discovery_locked_signature_genes.tsv", sep="\t")
    signature["rank_within_direction"] = pd.to_numeric(signature["rank_within_direction"])
    scores, selected, coverage, audit, primary_set_scores = read_all_cohorts(signature)

    paired_cohorts = {"GSE8671", "GSE72820"}
    test_rows = []
    for size in SIGNATURE_SIZES:
        for cohort in scores["cohort"].unique():
            cohort_scores = scores.loc[scores["cohort"].eq(cohort)]
            test_rows.append(
                cohort_comparison(
                    cohort_scores,
                    cohort,
                    size,
                    "adenoma",
                    "normal",
                    cohort in paired_cohorts,
                    "adenoma_vs_normal",
                )
            )
        for group_a, group_b, name in [
            ("hyperplastic", "normal", "hyperplastic_vs_normal"),
            ("adenoma", "hyperplastic", "adenoma_vs_hyperplastic"),
        ]:
            test_rows.append(
                cohort_comparison(
                    scores.loc[scores["cohort"].eq("GSE40362")],
                    "GSE40362",
                    size,
                    group_a,
                    group_b,
                    False,
                    name,
                )
            )
    tests = pd.DataFrame(test_rows)
    tests["q_value_within_comparison_and_size"] = tests.groupby(
        ["comparison", "signature_size_per_direction"]
    )["primary_p_value"].transform(lambda values: multipletests(values, method="fdr_bh")[1])

    rng = np.random.default_rng(SEED)
    bootstrap = pd.DataFrame(
        [
            bootstrap_primary(
                scores.loc[scores["cohort"].eq(cohort)],
                cohort,
                cohort in paired_cohorts,
                rng,
            )
            for cohort in scores["cohort"].unique()
        ]
    )

    models = pd.DataFrame([one_stage_model(scores, size) for size in SIGNATURE_SIZES])
    loo = pd.DataFrame(
        [
            one_stage_model(scores, PRIMARY_SIZE, excluded_cohort=cohort)
            for cohort in scores["cohort"].unique()
        ]
    )
    proliferation_adjusted = pd.DataFrame([proliferation_adjusted_model(scores)])
    proliferation_adjusted_loo = pd.DataFrame(
        [
            proliferation_adjusted_model(scores, excluded_cohort=cohort)
            for cohort in scores["cohort"].unique()
        ]
    )
    proliferation_adjusted_by_cohort = pd.DataFrame(
        [
            cohort_proliferation_adjusted_model(scores, cohort)
            for cohort in scores["cohort"].unique()
        ]
    )

    concordance_rows = []
    for cohort, part in scores.groupby("cohort", sort=False):
        for size in SIGNATURE_SIZES[:-1]:
            rho, p_value = stats.spearmanr(part[f"route_score_k{size}"], part[f"route_score_k{PRIMARY_SIZE}"])
            concordance_rows.append(
                {
                    "cohort": cohort,
                    "nested_size_per_direction": size,
                    "reference_size_per_direction": PRIMARY_SIZE,
                    "n_samples": len(part),
                    "spearman_rho": float(rho),
                    "p_value": float(p_value),
                }
            )
    concordance = pd.DataFrame(concordance_rows)

    grade_scores = scores.loc[scores["cohort"].eq("GSE41657")].copy()
    grade_scores["tissue_group"] = grade_scores["grade_group"]
    grade_tests = pd.DataFrame(
        [
            cohort_comparison(
                grade_scores,
                "GSE41657",
                PRIMARY_SIZE,
                group_a,
                group_b,
                False,
                comparison,
            )
            for group_a, group_b, comparison in [
                ("low_grade", "normal", "low_grade_vs_normal"),
                ("high_grade", "normal", "high_grade_vs_normal"),
                ("high_grade", "low_grade", "high_grade_vs_low_grade"),
            ]
        ]
    )
    grade_tests["q_value_bh"] = multipletests(grade_tests["primary_p_value"], method="fdr_bh")[1]

    adjusted_models = pd.DataFrame(
        [
            adjusted_cohort_model(scores, "GSE41657", "age_years", "gender"),
            adjusted_cohort_model(scores, "GSE40362", "age", "gender"),
        ]
    )

    normalization_sensitivity_rows = []
    for cohort in scores["cohort"].unique():
        original = scores.loc[scores["cohort"].eq(cohort)]
        restricted = primary_set_scores.loc[primary_set_scores["cohort"].eq(cohort)]
        original_test = cohort_comparison(
            original,
            cohort,
            PRIMARY_SIZE,
            "adenoma",
            "normal",
            cohort in paired_cohorts,
            "adenoma_vs_normal",
        )
        restricted_test = cohort_comparison(
            restricted,
            cohort,
            PRIMARY_SIZE,
            "adenoma",
            "normal",
            cohort in paired_cohorts,
            "adenoma_vs_normal",
        )
        rho, rho_p = stats.spearmanr(
            original[f"route_score_k{PRIMARY_SIZE}"],
            restricted[f"route_score_k{PRIMARY_SIZE}"],
        )
        normalization_sensitivity_rows.append(
            {
                "cohort": cohort,
                "reference_normalization": "all_cohort_samples",
                "sensitivity_normalization": "normal_and_adenoma_only",
                "n_samples_compared": len(original),
                "spearman_rho_scores": float(rho),
                "spearman_p_value": float(rho_p),
                "reference_auc": original_test["auc_a_vs_b"],
                "sensitivity_auc": restricted_test["auc_a_vs_b"],
                "reference_standardized_effect": original_test[
                    "clustered_standardized_mean_difference"
                ],
                "sensitivity_standardized_effect": restricted_test[
                    "clustered_standardized_mean_difference"
                ],
                "reference_primary_p": original_test["primary_p_value"],
                "sensitivity_primary_p": restricted_test["primary_p_value"],
            }
        )
    normalization_sensitivity = pd.DataFrame(normalization_sensitivity_rows)

    source_manifest().to_csv(OUT_DIR / "source_file_manifest.tsv", sep="\t", index=False)
    audit.to_csv(OUT_DIR / "cohort_audit.tsv", sep="\t", index=False)
    selected.to_csv(OUT_DIR / "selected_features_without_labels.tsv", sep="\t", index=False)
    coverage.to_csv(OUT_DIR / "locked_signature_coverage.tsv", sep="\t", index=False)
    scores.to_csv(OUT_DIR / "sample_scores.tsv", sep="\t", index=False)
    tests.to_csv(OUT_DIR / "cohort_tests.tsv", sep="\t", index=False)
    bootstrap.to_csv(OUT_DIR / "primary_cluster_bootstrap.tsv", sep="\t", index=False)
    models.to_csv(OUT_DIR / "one_stage_patient_cluster_models.tsv", sep="\t", index=False)
    loo.to_csv(OUT_DIR / "one_stage_leave_one_cohort_out.tsv", sep="\t", index=False)
    proliferation_adjusted.to_csv(
        OUT_DIR / "one_stage_proliferation_adjusted_model.tsv", sep="\t", index=False
    )
    proliferation_adjusted_loo.to_csv(
        OUT_DIR / "one_stage_proliferation_adjusted_leave_one_cohort_out.tsv",
        sep="\t",
        index=False,
    )
    proliferation_adjusted_by_cohort.to_csv(
        OUT_DIR / "cohort_proliferation_adjusted_models.tsv", sep="\t", index=False
    )
    concordance.to_csv(OUT_DIR / "nested_score_concordance.tsv", sep="\t", index=False)
    grade_tests.to_csv(OUT_DIR / "gse41657_grade_sensitivity.tsv", sep="\t", index=False)
    adjusted_models.to_csv(OUT_DIR / "age_sex_adjusted_models.tsv", sep="\t", index=False)
    normalization_sensitivity.to_csv(
        OUT_DIR / "normalization_scope_sensitivity.tsv", sep="\t", index=False
    )

    primary = tests.loc[
        tests["comparison"].eq("adenoma_vs_normal")
        & tests["signature_size_per_direction"].eq(PRIMARY_SIZE)
    ]
    lines = [
        "Independent sporadic colorectal adenoma validation",
        "==================================================",
        f"Candidate cohorts audited: {len(audit)}; eligible frozen-score cohorts: {scores['cohort'].nunique()}",
        f"Eligible normal/adenoma samples: {int(scores['tissue_group'].eq('normal').sum())}/{int(scores['tissue_group'].eq('adenoma').sum())}",
        f"Complete paired cohorts: GSE8671 n=32; GSE72820 n=7",
        "",
        "Locked 50-up/50-down results",
    ]
    for row in primary.itertuples(index=False):
        paired_text = (
            f"paired median delta={row.median_paired_difference:.3f}; paired P={row.p_paired_wilcoxon:.3g}"
            if row.paired_primary
            else f"patient-clustered mean difference={row.clustered_mean_difference:.3f} "
            f"(95% CI {row.clustered_ci_low:.3f} to {row.clustered_ci_high:.3f}; "
            f"P={row.p_patient_clustered_ols:.3g})"
        )
        lines.append(
            f"- {row.cohort}: n={row.n_a}/{row.n_b}; AUC={row.auc_a_vs_b:.3f}; {paired_text}."
        )
    full = models.loc[models["signature_size_per_direction"].eq(PRIMARY_SIZE)].iloc[0]
    full_proliferation_adjusted = proliferation_adjusted.iloc[0]
    lines.extend(
        [
            "",
            f"One-stage cohort-fixed, patient-clustered effect: beta={full['adenoma_coef_sd']:.3f} SD "
            f"(95% CI {full['ci_low']:.3f} to {full['ci_high']:.3f}; P={full['p_value']:.3g}; "
            f"n={int(full['n_samples'])} samples/{int(full['n_patient_clusters'])} patient clusters).",
            f"Leave-one-cohort-out beta range: {loo['adenoma_coef_sd'].min():.3f} to {loo['adenoma_coef_sd'].max():.3f} SD.",
            f"All nested-size one-stage effects positive: {bool((models['adenoma_coef_sd'] > 0).all())}.",
            f"After prespecified proliferation-control adjustment: beta={full_proliferation_adjusted['adenoma_coef_sd']:.3f} SD "
            f"(95% CI {full_proliferation_adjusted['ci_low']:.3f} to {full_proliferation_adjusted['ci_high']:.3f}; "
            f"P={full_proliferation_adjusted['p_value']:.3g}).",
            f"Proliferation-adjusted leave-one-cohort-out beta range: "
            f"{proliferation_adjusted_loo['adenoma_coef_sd'].min():.3f} to "
            f"{proliferation_adjusted_loo['adenoma_coef_sd'].max():.3f} SD.",
            "",
            "Interpretation boundary",
            "Processed public cohorts validate frozen score transport; they do not define an absolute diagnostic threshold.",
        ]
    )
    (OUT_DIR / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
