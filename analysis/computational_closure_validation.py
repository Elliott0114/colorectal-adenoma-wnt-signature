#!/usr/bin/env python3
"""Build a multi-layer computational validation closure for the locked route."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import binomtest, mannwhitneyu, wilcoxon


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "computational_closure_validation"
OUT.mkdir(parents=True, exist_ok=True)

SIGNATURE_PATH = ROOT / "results" / "route_signature_locked" / "discovery_locked_signature_genes.tsv"
MGI_PATH = ROOT / "data_sources" / "perturbation_closure" / "HOM_MouseHumanSequence.rpt"
PRIOR_PATH = ROOT / "data_sources" / "regulatory_priors" / "collectri_tcf_ascl2_wnt_controls_2026-08-08.tsv"
GSE114_DIR = ROOT / "data_sources" / "perturbation_closure" / "GSE114059"
GSE671_DIR = ROOT / "data_sources" / "perturbation_closure" / "GSE67186"
GSE130_DIR = ROOT / "data_sources" / "perturbation_closure" / "GSE130822"
GSE171_DIR = ROOT / "data_sources" / "perturbation_closure" / "GSE171910"
GSE125_PATH = ROOT / "data_sources" / "GSE125472_apc_ko_organoids" / "GSE125472_20181220_Results_Wnt_Signature_ALL.txt.gz"
GSE135_DIR = ROOT / "data_sources" / "GSE135328_tcf7l2_ko_crc"
GENE_INFO_PATH = ROOT / "data_sources" / "GSE117606" / "Homo_sapiens.gene_info.gz"
SPATIAL_DIR = ROOT / "data_sources" / "CRC_spatial_public" / "zenodo_7760264" / "extracted"
SPATIAL_ANNOT_DIR = SPATIAL_DIR / "Pathology_SpotAnnotations"

SEED = 20260808
N_PERMUTATIONS = 10_000
N_BOOTSTRAPS = 10_000

WNT_STEM = ["LGR5", "ASCL2", "OLFM4", "AXIN2", "SOX9", "EPHB2", "SMOC2"]
PROLIFERATION = ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "TYMS", "UBE2C", "CENPF"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "MUC13", "TACSTD2", "CDH1"]
FEATURES = ["route_score", "route_up", "route_down", "wnt_stem", "proliferation_control"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_signature() -> pd.DataFrame:
    signature = pd.read_csv(SIGNATURE_PATH, sep="\t")
    if len(signature) != 100:
        raise ValueError(f"Expected 100 locked genes, found {len(signature)}")
    if signature["validation_used_for_selection"].astype(str).str.lower().ne("false").any():
        raise ValueError("Validation-selected genes found in locked signature")
    signature["route_weight"] = np.where(signature["signature_direction"].eq("adenoma_up"), 1, -1)
    return signature[["gene", "signature_direction", "route_weight"]].copy()


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    library = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(library, axis=1) * 1_000_000 + 1)


def aggregate_gene_rows(expression: pd.DataFrame) -> pd.DataFrame:
    expression = expression.copy()
    expression.index = expression.index.astype(str).str.strip()
    valid_index = expression.index.notna() & (expression.index != "") & (expression.index != "nan")
    expression = expression[valid_index]
    return expression.groupby(level=0, sort=True).mean()


def gene_zscores(expression: pd.DataFrame) -> pd.DataFrame:
    means = expression.mean(axis=1)
    sds = expression.std(axis=1, ddof=1).replace(0, np.nan)
    return expression.sub(means, axis=0).div(sds, axis=0)


def mean_present(z: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, int]:
    present = [gene for gene in genes if gene in z.index]
    if not present:
        return pd.Series(np.nan, index=z.columns), 0
    return z.loc[present].mean(axis=0), len(present)


def score_expression(expression: pd.DataFrame, signature: pd.DataFrame, dataset: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expression = aggregate_gene_rows(expression)
    z = gene_zscores(expression)
    up = signature.loc[signature["route_weight"].eq(1), "gene"].tolist()
    down = signature.loc[signature["route_weight"].eq(-1), "gene"].tolist()
    route_up, n_up = mean_present(z, up)
    route_down, n_down = mean_present(z, down)
    wnt, n_wnt = mean_present(z, WNT_STEM)
    proliferation, n_prolif = mean_present(z, PROLIFERATION)
    epithelial, n_epi = mean_present(z, EPITHELIAL)
    scores = pd.DataFrame(
        {
            "sample_id": z.columns,
            "route_up": route_up.reindex(z.columns).to_numpy(),
            "route_down": route_down.reindex(z.columns).to_numpy(),
            "route_score": (route_up - route_down).reindex(z.columns).to_numpy(),
            "wnt_stem": wnt.reindex(z.columns).to_numpy(),
            "proliferation_control": proliferation.reindex(z.columns).to_numpy(),
            "epithelial_control": epithelial.reindex(z.columns).to_numpy(),
        }
    )
    coverage = pd.DataFrame(
        [
            {"dataset": dataset, "feature": "route_up", "n_expected": 50, "n_present": n_up},
            {"dataset": dataset, "feature": "route_down", "n_expected": 50, "n_present": n_down},
            {"dataset": dataset, "feature": "wnt_stem", "n_expected": len(WNT_STEM), "n_present": n_wnt},
            {"dataset": dataset, "feature": "proliferation_control", "n_expected": len(PROLIFERATION), "n_present": n_prolif},
            {"dataset": dataset, "feature": "epithelial_control", "n_expected": len(EPITHELIAL), "n_present": n_epi},
        ]
    )
    return scores, z, coverage


def bootstrap_mean(values: np.ndarray, seed: int) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = rng.choice(values, size=(N_BOOTSTRAPS, len(values)), replace=True).mean(axis=1)
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def exact_sign_p(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    nonzero = values[values != 0]
    if len(nonzero) == 0:
        return np.nan
    return float(binomtest(int((nonzero > 0).sum()), len(nonzero), 0.5, alternative="two-sided").pvalue)


def feature_expected_direction(feature: str, route_direction: int) -> int:
    if feature in {"route_score", "route_up", "wnt_stem"}:
        return route_direction
    if feature == "route_down":
        return -route_direction
    return 0


def contrast_weights(metadata: pd.DataFrame, comparison: dict[str, object]) -> tuple[np.ndarray, list[str]]:
    metadata = metadata.reset_index(drop=True)
    weights = np.zeros(len(metadata), dtype=float)
    eligible = metadata.copy()
    units = comparison.get("units")
    if units is not None:
        eligible = eligible[eligible["unit_id"].isin(units)]
    complete = []
    for unit_id, group in eligible.groupby("unit_id", observed=True):
        target_idx = group.index[group["condition"].eq(comparison["target_condition"])].to_numpy()
        reference_idx = group.index[group["condition"].eq(comparison["reference_condition"])].to_numpy()
        if len(target_idx) and len(reference_idx):
            complete.append(unit_id)
            weights[target_idx] = 1 / len(target_idx)
            weights[reference_idx] = -1 / len(reference_idx)
    if complete:
        weights /= len(complete)
    return weights, complete


def comparison_effects(
    dataset: str,
    species: str,
    model_system: str,
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    comparisons: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    merged = metadata.merge(scores, on="sample_id", validate="one_to_one")
    unit_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for comparison in comparisons:
        subset = merged
        if comparison.get("units") is not None:
            subset = subset[subset["unit_id"].isin(comparison["units"])]
        for feature in FEATURES:
            wide = subset.groupby(["unit_id", "condition"], observed=True)[feature].mean().unstack()
            if comparison["target_condition"] not in wide or comparison["reference_condition"] not in wide:
                continue
            diff = (wide[comparison["target_condition"]] - wide[comparison["reference_condition"]]).dropna()
            expected = feature_expected_direction(feature, int(comparison["expected_route_direction"]))
            for unit_id, value in diff.items():
                unit_rows.append(
                    {
                        "dataset": dataset,
                        "species": species,
                        "model_system": model_system,
                        "comparison": comparison["comparison"],
                        "unit_id": unit_id,
                        "feature": feature,
                        "expected_direction": expected,
                        "difference": float(value),
                    }
                )
            ci_low, ci_high = bootstrap_mean(diff.to_numpy(), SEED + len(summary_rows))
            summary_rows.append(
                {
                    "dataset": dataset,
                    "species": species,
                    "model_system": model_system,
                    "comparison": comparison["comparison"],
                    "feature": feature,
                    "expected_direction": expected,
                    "n_units": len(diff),
                    "mean_difference": float(diff.mean()) if len(diff) else np.nan,
                    "median_difference": float(diff.median()) if len(diff) else np.nan,
                    "bootstrap_mean_ci_low": ci_low,
                    "bootstrap_mean_ci_high": ci_high,
                    "n_expected_direction": int((expected * diff.to_numpy() > 0).sum()) if expected else np.nan,
                    "p_exact_sign_two_sided": exact_sign_p(diff.to_numpy()),
                }
            )
    return merged, pd.DataFrame(unit_rows), pd.DataFrame(summary_rows)


def expression_matched_test(
    dataset: str,
    expression: pd.DataFrame,
    z: pd.DataFrame,
    signature: pd.DataFrame,
    metadata: pd.DataFrame,
    comparisons: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    valid = z.index[z.notna().sum(axis=1).ge(2) & z.std(axis=1, ddof=1).gt(0)]
    z = z.loc[valid]
    expression = expression.reindex(valid)
    signature_present = signature[signature["gene"].isin(valid)].copy()
    signature_genes = set(signature_present["gene"])
    means = expression.mean(axis=1)
    ranks = means.rank(method="first", pct=True)
    bins = np.minimum((ranks * 10).astype(int), 9)
    candidate_by_bin = {
        int(group): np.array([gene for gene in genes.index if gene not in signature_genes], dtype=object)
        for group, genes in bins.groupby(bins, observed=True)
    }
    up = signature_present.loc[signature_present["route_weight"].eq(1), "gene"].tolist()
    down = signature_present.loc[signature_present["route_weight"].eq(-1), "gene"].tolist()
    up_bins = bins.reindex(up).value_counts().to_dict()
    down_bins = bins.reindex(down).value_counts().to_dict()
    gene_to_idx = {gene: idx for idx, gene in enumerate(z.index)}
    z_values = z.to_numpy(float)
    metadata = metadata.set_index("sample_id").loc[z.columns].reset_index()
    dataset_seed = int.from_bytes(hashlib.sha256(dataset.encode("utf-8")).digest()[:4], "little")
    rng = np.random.default_rng(SEED + dataset_seed % 10_000)
    null_rows = []
    test_rows = []
    for comparison in comparisons:
        weights, complete_units = contrast_weights(metadata, comparison)
        if not complete_units:
            continue
        observed_score = z.loc[up].mean(axis=0).to_numpy() - z.loc[down].mean(axis=0).to_numpy()
        observed = float(observed_score @ weights)
        null = np.empty(N_PERMUTATIONS, dtype=float)
        for permutation in range(N_PERMUTATIONS):
            random_up: list[str] = []
            random_down: list[str] = []
            for expression_bin in sorted(set(up_bins) | set(down_bins)):
                n_up = int(up_bins.get(expression_bin, 0))
                n_down = int(down_bins.get(expression_bin, 0))
                pool = candidate_by_bin[int(expression_bin)]
                chosen = rng.choice(pool, size=n_up + n_down, replace=False)
                random_up.extend(chosen[:n_up].tolist())
                random_down.extend(chosen[n_up:].tolist())
            up_idx = [gene_to_idx[gene] for gene in random_up]
            down_idx = [gene_to_idx[gene] for gene in random_down]
            random_score = np.nanmean(z_values[up_idx], axis=0) - np.nanmean(z_values[down_idx], axis=0)
            null[permutation] = float(random_score @ weights)
        expected = int(comparison["expected_route_direction"])
        if expected > 0:
            p_value = (1 + int((null >= observed).sum())) / (N_PERMUTATIONS + 1)
        elif expected < 0:
            p_value = (1 + int((null <= observed).sum())) / (N_PERMUTATIONS + 1)
        else:
            p_value = (1 + int((np.abs(null) >= abs(observed)).sum())) / (N_PERMUTATIONS + 1)
        test_rows.append(
            {
                "dataset": dataset,
                "comparison": comparison["comparison"],
                "expected_direction": expected,
                "n_units": len(complete_units),
                "observed_mean_difference": observed,
                "null_mean": float(null.mean()),
                "null_sd": float(null.std(ddof=1)),
                "null_ci_low": float(np.quantile(null, 0.025)),
                "null_ci_high": float(np.quantile(null, 0.975)),
                "p_expression_matched_one_sided": p_value,
                "n_permutations": N_PERMUTATIONS,
                "n_up": len(up),
                "n_down": len(down),
            }
        )
        null_rows.extend(
            {"dataset": dataset, "comparison": comparison["comparison"], "permutation": i + 1, "null_effect": value}
            for i, value in enumerate(null)
        )
    return pd.DataFrame(test_rows), pd.DataFrame(null_rows)


def strict_human_mouse_map() -> tuple[dict[str, str], pd.DataFrame]:
    homology = pd.read_csv(MGI_PATH, sep="\t", dtype=str)
    rows = []
    mapping = {}
    for class_key, group in homology.groupby("DB Class Key", observed=True):
        human = group.loc[group["NCBI Taxon ID"].eq("9606"), "Symbol"].dropna().unique()
        mouse = group.loc[group["NCBI Taxon ID"].eq("10090"), "Symbol"].dropna().unique()
        if len(human) == 1 and len(mouse) == 1:
            mapping[str(human[0])] = str(mouse[0])
            rows.append({"homology_class": class_key, "human_gene": human[0], "mouse_gene": mouse[0]})
    return mapping, pd.DataFrame(rows)


def mouse_to_human_expression(expression: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    inverse = {mouse: human for human, mouse in mapping.items()}
    keep = [gene for gene in expression.index.astype(str) if gene in inverse]
    out = expression.loc[keep].copy()
    out.index = [inverse[str(gene)] for gene in out.index]
    return aggregate_gene_rows(out)


def load_gse114059() -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    matrix_path = GSE114_DIR / "GSE114059_series_matrix.txt.gz"
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    title_line = next(line for line in lines if line.startswith("!Sample_title"))
    accession_line = next(line for line in lines if line.startswith("!Sample_geo_accession"))
    titles = [x.strip().strip('"') for x in title_line.rstrip().split("\t")[1:]]
    accessions = [x.strip().strip('"') for x in accession_line.rstrip().split("\t")[1:]]
    table_start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin")) + 1
    table_end = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_end"))
    matrix = pd.read_csv(
        io.StringIO("".join(lines[table_start:table_end])), sep="\t", index_col=0
    )
    matrix.index = matrix.index.astype(str).str.strip('"')
    matrix.columns = [str(c).strip('"') for c in matrix.columns]

    annotation = pd.read_csv(GSE114_DIR / "GPL10558.annot.gz", sep="\t", skiprows=28, dtype=str)
    annotation = annotation[["ID", "Gene symbol"]].dropna()
    annotation["gene"] = annotation["Gene symbol"].str.split(" /// ").str[0].str.strip()
    annotation = annotation[annotation["gene"].ne("")]
    probe = annotation.merge(matrix.var(axis=1).rename("variance"), left_on="ID", right_index=True, how="inner")
    probe = probe.sort_values(["gene", "variance", "ID"], ascending=[True, False, True]).drop_duplicates("gene")
    expression = matrix.loc[probe["ID"]].copy()
    expression.index = probe["gene"].to_numpy()
    metadata_rows = []
    for accession, title in zip(accessions, titles, strict=True):
        organoid_match = re.search(r"PTO\s+(\d+T)", title)
        replicate_match = re.search(r"R(\d+)$", title)
        if "PRI-724" in title:
            condition = "trametinib_plus_pri724"
        elif "Trametinib" in title:
            condition = "trametinib"
        else:
            condition = "dmso"
        metadata_rows.append(
            {
                "sample_id": accession,
                "unit_id": f"PTO{organoid_match.group(1)}",
                "condition": condition,
                "replicate": int(replicate_match.group(1)),
                "sample_title": title,
            }
        )
    comparisons = [
        {
            "comparison": "trametinib_vs_dmso",
            "target_condition": "trametinib",
            "reference_condition": "dmso",
            "expected_route_direction": 1,
            "units": None,
        },
        {
            "comparison": "pri724_reversal_of_trametinib",
            "target_condition": "trametinib_plus_pri724",
            "reference_condition": "trametinib",
            "expected_route_direction": -1,
            "units": ["PTO19T"],
        },
    ]
    return expression.astype(float), pd.DataFrame(metadata_rows), comparisons


def load_gse67186(mapping: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    raw = pd.read_csv(GSE671_DIR / "GSE67186_counts_scaled_DESeq.xls.gz", sep="\t")
    expression = raw.set_index("GeneSymbol").drop(columns=["GeneID"]).astype(float)
    expression = mouse_to_human_expression(np.log2(expression + 1), mapping)
    group_map = {
        "s_AD": ("shApc", "on_dox"),
        "s_AOD": ("shApc", "off_dox"),
        "s_AKD": ("shApc_Kras", "on_dox"),
        "s_AKOD": ("shApc_Kras", "off_dox"),
        "s_RD": ("shRenilla", "on_dox"),
        "s_ROD": ("shRenilla", "off_dox"),
    }
    metadata_rows = []
    for sample in expression.columns:
        prefix = next(prefix for prefix in group_map if sample.startswith(prefix + "_"))
        unit_id, condition = group_map[prefix]
        metadata_rows.append({"sample_id": sample, "unit_id": unit_id, "condition": condition})
    comparisons = [
        {
            "comparison": "apc_restoration_shApc",
            "target_condition": "off_dox",
            "reference_condition": "on_dox",
            "expected_route_direction": -1,
            "units": ["shApc"],
        },
        {
            "comparison": "apc_restoration_shApc_Kras",
            "target_condition": "off_dox",
            "reference_condition": "on_dox",
            "expected_route_direction": -1,
            "units": ["shApc_Kras"],
        },
        {
            "comparison": "doxycycline_control_shRenilla",
            "target_condition": "off_dox",
            "reference_condition": "on_dox",
            "expected_route_direction": 0,
            "units": ["shRenilla"],
        },
    ]
    return expression, pd.DataFrame(metadata_rows), comparisons


def load_gse130822(mapping: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    raw = pd.read_csv(GSE130_DIR / "GSE130822_Read_Counts_All_Cells.txt.gz", sep="\t")
    relevant = [c for c in raw.columns if c.startswith(("Ascl2_KO_CSC", "Resting_CSC"))]
    counts = raw.set_index("Gene")[relevant].astype(float)
    expression = mouse_to_human_expression(log2_cpm(counts), mapping)
    metadata_rows = []
    for sample in relevant:
        condition = "ascl2_ko" if sample.startswith("Ascl2_KO") else "resting_wt"
        metadata_rows.append({"sample_id": sample, "unit_id": "colonic_stem_cells", "condition": condition})
    comparisons = [
        {
            "comparison": "ascl2_ko_vs_resting_wt",
            "target_condition": "ascl2_ko",
            "reference_condition": "resting_wt",
            "expected_route_direction": -1,
            "units": ["colonic_stem_cells"],
        }
    ]
    return expression, pd.DataFrame(metadata_rows), comparisons


def load_gse171910() -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    raw = pd.read_excel(GSE171_DIR / "GSE171910_RNA_seq_results.xlsx", sheet_name="Reads FPKM")
    count_cols = [c for c in raw.columns if c.endswith("_readcount")]
    counts = raw.set_index("GeneName")[count_cols].astype(float)
    counts.columns = [c.removesuffix("_readcount") for c in count_cols]
    expression = log2_cpm(aggregate_gene_rows(counts))
    metadata_rows = []
    for sample in expression.columns:
        if sample.startswith("BC_"):
            unit = "SIBC"
            condition = "wnt_off" if "DOX" in sample else "control"
        elif sample.startswith("DL_"):
            unit = "DLD1"
            condition = "wnt_off" if "DOX" in sample else "control"
        elif sample.startswith("HT_"):
            unit = "HT29"
            condition = "wnt_off" if "ZN" in sample else "control"
        elif sample.startswith("L8_"):
            unit = "L8"
            condition = "wnt_off" if "DOX" in sample else "control"
        else:
            raise ValueError(f"Unrecognized GSE171910 sample {sample}")
        metadata_rows.append({"sample_id": sample, "unit_id": unit, "condition": condition})
    comparisons = [
        {
            "comparison": "conditional_wnt_silencing",
            "target_condition": "wnt_off",
            "reference_condition": "control",
            "expected_route_direction": -1,
            "units": None,
        }
    ]
    return expression, pd.DataFrame(metadata_rows), comparisons


def signed_collectri(path: Path) -> pd.DataFrame:
    prior = pd.read_csv(path, sep="\t")
    for col in ["consensus_direction", "consensus_stimulation", "consensus_inhibition"]:
        prior[col] = prior[col].astype(str).str.lower().eq("true")
    prior = prior[prior["consensus_direction"]].copy()
    prior["edge_sign"] = np.select(
        [prior["consensus_stimulation"] & ~prior["consensus_inhibition"], prior["consensus_inhibition"] & ~prior["consensus_stimulation"]],
        [1.0, -1.0],
        default=np.nan,
    )
    prior = prior.dropna(subset=["edge_sign"])
    grouped = prior.groupby(["source_genesymbol", "target_genesymbol"], observed=True)["edge_sign"]
    consensus = grouped.agg(["min", "max"]).reset_index()
    consensus = consensus[consensus["min"].eq(consensus["max"])].copy()
    consensus["edge_sign"] = consensus["min"]
    return consensus[["source_genesymbol", "target_genesymbol", "edge_sign"]]


def ulm_activity(z: pd.DataFrame, prior: pd.DataFrame, dataset: str, min_targets: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    activity_rows = []
    coverage_rows = []
    for tf, regulon in prior.groupby("source_genesymbol", observed=True):
        regulon = regulon[regulon["target_genesymbol"].isin(z.index)].drop_duplicates("target_genesymbol")
        n_targets = len(regulon)
        coverage_rows.append({"dataset": dataset, "tf": tf, "n_targets": n_targets})
        if n_targets < min_targets:
            continue
        x = regulon["edge_sign"].to_numpy(float)
        x = x - x.mean()
        ssx = float(np.sum(x**2))
        if ssx == 0:
            continue
        y_matrix = z.loc[regulon["target_genesymbol"]].to_numpy(float)
        for sample_idx, sample_id in enumerate(z.columns):
            y = y_matrix[:, sample_idx]
            valid = np.isfinite(y)
            if valid.sum() < min_targets:
                continue
            xv = x[valid]
            yv = y[valid]
            yv = yv - yv.mean()
            ssx_v = float(np.sum(xv**2))
            beta = float(np.sum(xv * yv) / ssx_v)
            residual = yv - beta * xv
            dof = len(yv) - 2
            se = float(np.sqrt(np.sum(residual**2) / dof / ssx_v)) if dof > 0 else np.nan
            t_value = beta / se if se and np.isfinite(se) else np.nan
            activity_rows.append(
                {"dataset": dataset, "sample_id": sample_id, "tf": tf, "ulm_activity": t_value, "n_targets": int(valid.sum())}
            )
    return pd.DataFrame(activity_rows), pd.DataFrame(coverage_rows)


def ulm_comparison_effects(
    dataset: str,
    activity: pd.DataFrame,
    metadata: pd.DataFrame,
    comparisons: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if activity.empty:
        return pd.DataFrame(), pd.DataFrame()
    wide_activity = activity.pivot(index="sample_id", columns="tf", values="ulm_activity")
    unit_rows = []
    summary_rows = []
    for tf in wide_activity.columns:
        score = wide_activity[[tf]].rename(columns={tf: "tf_activity"}).reset_index()
        merged = metadata.merge(score, on="sample_id", how="inner")
        for comparison in comparisons:
            subset = merged
            if comparison.get("units") is not None:
                subset = subset[subset["unit_id"].isin(comparison["units"])]
            pivot = subset.groupby(["unit_id", "condition"], observed=True)["tf_activity"].mean().unstack()
            if comparison["target_condition"] not in pivot or comparison["reference_condition"] not in pivot:
                continue
            diff = (pivot[comparison["target_condition"]] - pivot[comparison["reference_condition"]]).dropna()
            expected = int(comparison["expected_route_direction"])
            for unit, value in diff.items():
                unit_rows.append(
                    {"dataset": dataset, "comparison": comparison["comparison"], "tf": tf, "unit_id": unit, "difference": float(value)}
                )
            summary_rows.append(
                {
                    "dataset": dataset,
                    "comparison": comparison["comparison"],
                    "tf": tf,
                    "expected_direction": expected,
                    "n_units": len(diff),
                    "mean_difference": float(diff.mean()) if len(diff) else np.nan,
                    "n_expected_direction": int((expected * diff.to_numpy() > 0).sum()) if expected else np.nan,
                }
            )
    return pd.DataFrame(unit_rows), pd.DataFrame(summary_rows)


def decode(values) -> list[str]:
    return [v.decode("utf-8", "replace") if isinstance(v, bytes) else str(v) for v in values]


def read_10x_h5(path: Path) -> tuple[list[str], list[str], csc_matrix]:
    with h5py.File(path, "r") as handle:
        barcodes = decode(handle["matrix/barcodes"][()])
        genes = decode(handle["matrix/features/name"][()])
        matrix = csc_matrix(
            (handle["matrix/data"][()], handle["matrix/indices"][()], handle["matrix/indptr"][()]),
            shape=tuple(handle["matrix/shape"][()]),
        )
    return barcodes, genes, matrix


def pathology_group(value: str) -> str:
    text = str(value).strip().lower()
    if not text or text == "nan":
        return "missing"
    if "exclude" in text:
        return "exclude"
    if "tumor" in text and "stroma" in text:
        return "tumor_stroma"
    if "tumor" in text:
        return "tumor"
    if "stroma" in text or "fibroblastic" in text:
        return "stroma"
    if "epithelium" in text:
        return "non_neoplastic_epithelium"
    if "submucosa" in text:
        return "submucosa"
    return "other"


def spatial_locked_route(signature: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected_genes = sorted(set(signature["gene"]) | set(WNT_STEM) | set(PROLIFERATION) | set(EPITHELIAL))
    spot_frames = []
    coverage_rows = []
    for sample_dir in sorted(path for path in SPATIAL_DIR.iterdir() if path.is_dir() and (path / "filtered_feature_bc_matrix.h5").exists()):
        sample_id = sample_dir.name
        barcodes, genes, matrix = read_10x_h5(sample_dir / "filtered_feature_bc_matrix.h5")
        gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
        present = [gene for gene in selected_genes if gene in gene_to_idx]
        selected = matrix[[gene_to_idx[g] for g in present], :].toarray().astype(float)
        total = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
        expression = pd.DataFrame(np.log2(selected / np.where(total > 0, total, np.nan) * 1_000_000 + 1), index=present, columns=barcodes)
        scores, _, coverage = score_expression(expression, signature, f"spatial_{sample_id}")
        positions = pd.read_csv(sample_dir / "spatial" / "tissue_positions_list.csv", header=None)
        positions.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
        ann_path = SPATIAL_ANNOT_DIR / f"Pathologist_Annotations_{sample_id}.csv"
        annotation = pd.read_csv(ann_path)
        ann_col = next((c for c in annotation.columns if c.lower().startswith("pathologist annotation")), None)
        if ann_col is None:
            ann_col = next(c for c in annotation.columns if c != "Barcode")
        annotation = annotation.rename(columns={"Barcode": "barcode", ann_col: "pathology_annotation"})
        annotation["pathology_group"] = annotation["pathology_annotation"].map(pathology_group)
        frame = scores.rename(columns={"sample_id": "barcode"}).merge(positions, on="barcode", how="left")
        frame = frame.merge(annotation[["barcode", "pathology_annotation", "pathology_group"]], on="barcode", how="left")
        frame.insert(0, "sample_id", sample_id)
        relevant = frame["pathology_group"].isin(["tumor", "tumor_stroma", "non_neoplastic_epithelium", "stroma"])
        design = np.column_stack(
            [
                np.ones(relevant.sum()),
                frame.loc[relevant, "proliferation_control"].to_numpy(float),
                frame.loc[relevant, "epithelial_control"].to_numpy(float),
            ]
        )
        outcome = frame.loc[relevant, "route_score"].to_numpy(float)
        valid = np.isfinite(design).all(axis=1) & np.isfinite(outcome)
        residual = np.full(relevant.sum(), np.nan)
        if valid.sum() > 10:
            beta = np.linalg.lstsq(design[valid], outcome[valid], rcond=None)[0]
            residual[valid] = outcome[valid] - design[valid] @ beta
        frame["route_residual_prolif_epithelial"] = np.nan
        frame.loc[relevant, "route_residual_prolif_epithelial"] = residual
        spot_frames.append(frame)
        coverage_rows.append(coverage)
    spots = pd.concat(spot_frames, ignore_index=True)
    coverage = pd.concat(coverage_rows, ignore_index=True)
    summaries = (
        spots.groupby(["sample_id", "pathology_group"], observed=True)[
            ["route_score", "route_residual_prolif_epithelial", "wnt_stem", "proliferation_control", "epithelial_control"]
        ]
        .median()
        .reset_index()
    )
    unit_rows = []
    test_rows = []
    for comparison, target_groups, reference_groups in [
        ("tumor_vs_non_neoplastic_epithelium", ["tumor", "tumor_stroma"], ["non_neoplastic_epithelium"]),
        ("tumor_vs_stroma", ["tumor", "tumor_stroma"], ["stroma"]),
    ]:
        target = summaries[summaries["pathology_group"].isin(target_groups)].groupby("sample_id", observed=True).median(numeric_only=True)
        reference = summaries[summaries["pathology_group"].isin(reference_groups)].groupby("sample_id", observed=True).median(numeric_only=True)
        joined = target.join(reference, lsuffix="__target", rsuffix="__reference", how="inner")
        for feature in ["route_score", "route_residual_prolif_epithelial", "wnt_stem", "proliferation_control"]:
            diff = joined[f"{feature}__target"] - joined[f"{feature}__reference"]
            for sample_id, value in diff.items():
                unit_rows.append({"comparison": comparison, "sample_id": sample_id, "feature": feature, "difference": float(value)})
            p_value = float(wilcoxon(diff, zero_method="wilcox", alternative="two-sided").pvalue) if len(diff) >= 4 and (diff != 0).any() else np.nan
            test_rows.append(
                {
                    "comparison": comparison,
                    "feature": feature,
                    "n_sections": len(diff),
                    "median_difference": float(diff.median()) if len(diff) else np.nan,
                    "mean_difference": float(diff.mean()) if len(diff) else np.nan,
                    "n_positive": int((diff > 0).sum()),
                    "p_wilcoxon": p_value,
                }
            )
    return spots, summaries, pd.DataFrame(unit_rows), pd.DataFrame(test_rows)


def load_gse125_expression() -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    raw = pd.read_csv(GSE125_PATH, sep="\t", decimal=",")
    sample_cols = [c for c in raw.columns if re.match(r"^Donor\d+_(WT|APC)_(w|wo)_Wnt$", c)]
    counts = raw[sample_cols].apply(pd.to_numeric, errors="raise")
    counts.index = raw["hgnc_symbol"].fillna(raw["symbol"]).astype(str)
    expression = log2_cpm(aggregate_gene_rows(counts))
    rows = []
    for sample in sample_cols:
        match = re.match(r"^(Donor\d+)_(WT|APC)_(w|wo)_Wnt$", sample)
        genotype = match.group(2)
        medium = "with" if match.group(3) == "w" else "without"
        rows.append({"sample_id": sample, "unit_id": match.group(1), "condition": f"{genotype}_{medium}"})
    comparisons = [
        {
            "comparison": "apc_ko_vs_wt_with_wnt",
            "target_condition": "APC_with",
            "reference_condition": "WT_with",
            "expected_route_direction": 1,
            "units": None,
        }
    ]
    return expression, pd.DataFrame(rows), comparisons


def ensembl_to_symbol() -> dict[str, str]:
    mapping = {}
    with gzip.open(GENE_INFO_PATH, "rt", encoding="utf-8", errors="replace") as handle:
        info = pd.read_csv(handle, sep="\t", dtype=str)
    for _, row in info.iterrows():
        symbol = str(row["Symbol"])
        for token in str(row.get("dbXrefs", "")).split("|"):
            if token.startswith("Ensembl:"):
                mapping[token.split(":", 1)[1]] = symbol
    return mapping


def load_gse135_expression() -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    mapping = ensembl_to_symbol()
    out = {}
    for cell_line in ["HCT116", "HT29"]:
        raw = pd.read_csv(GSE135_DIR / f"GSE135328_count_{cell_line}.txt.gz", sep="\t", index_col=0)
        raw.index = [mapping.get(str(gene).split(".")[0], str(gene)) for gene in raw.index]
        expression = log2_cpm(aggregate_gene_rows(raw.astype(float)))
        rows = []
        for sample in expression.columns:
            clone = sample.rsplit("_", 1)[0]
            clone_num = clone.split("_")[-1]
            if cell_line == "HCT116":
                genotype = "WT" if clone_num == "3" else "KO"
            else:
                genotype = "WT" if clone_num == "56" else ("Het" if clone_num == "62" else "KO")
            rows.append({"sample_id": sample, "clone_id": clone, "condition": genotype})
        out[cell_line] = (expression, pd.DataFrame(rows))
    return out


def gse135_ulm_effects(activity: pd.DataFrame, metadata: pd.DataFrame, cell_line: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = activity.merge(metadata, on="sample_id", how="inner")
    unit_rows = []
    summary_rows = []
    for tf, tf_frame in merged.groupby("tf", observed=True):
        clone_activity = tf_frame.groupby(["clone_id", "condition"], observed=True)["ulm_activity"].mean().reset_index()
        wt = clone_activity.loc[clone_activity["condition"].eq("WT"), "ulm_activity"]
        if wt.empty:
            continue
        reference = float(wt.iloc[0])
        ko = clone_activity[clone_activity["condition"].eq("KO")]
        diff = ko["ulm_activity"] - reference
        for (_, row), value in zip(ko.iterrows(), diff, strict=True):
            unit_rows.append(
                {"dataset": f"GSE135328_{cell_line}", "comparison": "tcf7l2_ko_vs_wt", "tf": tf, "unit_id": row["clone_id"], "difference": float(value)}
            )
        summary_rows.append(
            {
                "dataset": f"GSE135328_{cell_line}",
                "comparison": "tcf7l2_ko_vs_wt",
                "tf": tf,
                "expected_direction": -1,
                "n_units": len(diff),
                "mean_difference": float(diff.mean()) if len(diff) else np.nan,
                "n_expected_direction": int((diff < 0).sum()),
            }
        )
    return pd.DataFrame(unit_rows), pd.DataFrame(summary_rows)


def append_existing_perturbations() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = ROOT / "results" / "perturbation_validation_locked_route"
    gse125_units = pd.read_csv(base / "gse125472_donor_contrasts.tsv", sep="\t")
    gse125_units = gse125_units[gse125_units["comparison"].isin(["APC_vs_WT_with_Wnt", "APC_vs_WT_without_Wnt"])]
    gse125_units = gse125_units.rename(columns={"donor_id": "unit_id"})
    gse125_units.insert(0, "dataset", "GSE125472")
    gse125_units.insert(1, "species", "human")
    gse125_units.insert(2, "model_system", "human_colon_organoid")

    gse125_summary = pd.read_csv(base / "gse125472_contrast_summary.tsv", sep="\t")
    gse125_summary = gse125_summary[gse125_summary["comparison"].isin(["APC_vs_WT_with_Wnt", "APC_vs_WT_without_Wnt"])]
    gse125_summary.insert(0, "dataset", "GSE125472")
    gse125_summary.insert(1, "species", "human")
    gse125_summary.insert(2, "model_system", "human_colon_organoid")

    gse135 = pd.read_csv(base / "gse135328_clone_contrasts.tsv", sep="\t")
    gse135 = gse135[gse135["genotype"].eq("KO") & gse135["feature"].isin(FEATURES)].copy()
    gse135_units = gse135.rename(columns={"clone_id": "unit_id", "difference_vs_WT": "difference"})
    gse135_units["dataset"] = "GSE135328_" + gse135_units["cell_line"]
    gse135_units["species"] = "human"
    gse135_units["model_system"] = "crc_cell_line"
    gse135_units["comparison"] = "tcf7l2_ko_vs_wt"
    gse135_units["expected_direction"] = gse135_units["feature"].map(
        {"route_score": -1, "route_up": -1, "route_down": 1, "wnt_stem": -1, "proliferation_control": 0}
    )
    keep_cols = ["dataset", "species", "model_system", "comparison", "unit_id", "feature", "expected_direction", "difference"]
    gse135_units = gse135_units[keep_cols]
    summary_rows = []
    for keys, group in gse135_units.groupby(["dataset", "species", "model_system", "comparison", "feature", "expected_direction"], observed=True):
        dataset, species, model_system, comparison, feature, expected = keys
        values = group["difference"].to_numpy(float)
        ci_low, ci_high = bootstrap_mean(values, SEED + len(summary_rows))
        summary_rows.append(
            {
                "dataset": dataset,
                "species": species,
                "model_system": model_system,
                "comparison": comparison,
                "feature": feature,
                "expected_direction": expected,
                "n_units": len(values),
                "mean_difference": float(np.mean(values)),
                "median_difference": float(np.median(values)),
                "bootstrap_mean_ci_low": ci_low,
                "bootstrap_mean_ci_high": ci_high,
                "n_expected_direction": int((expected * values > 0).sum()) if expected else np.nan,
                "p_exact_sign_two_sided": exact_sign_p(values),
            }
        )
    gse125_matched = pd.read_csv(base / "gse125472_matched_signature_tests.tsv", sep="\t")
    gse125_matched = gse125_matched[
        gse125_matched["comparison"].isin(["APC_vs_WT_with_Wnt", "APC_vs_WT_without_Wnt"])
    ].copy()
    gse135_matched = pd.read_csv(base / "gse135328_matched_signature_tests.tsv", sep="\t").rename(
        columns={
            "observed_ko_minus_wt": "observed_mean_difference",
            "p_expression_matched_lower_tail": "p_expression_matched_one_sided",
        }
    )
    gse135_matched["dataset"] = "GSE135328_" + gse135_matched["cell_line"]
    gse135_matched["comparison"] = "tcf7l2_ko_vs_wt"
    gse135_matched["expected_direction"] = -1
    matched = pd.concat([gse125_matched, gse135_matched], ignore_index=True, sort=False)
    return pd.concat([gse125_units[keep_cols], gse135_units], ignore_index=True), pd.concat([gse125_summary, pd.DataFrame(summary_rows)], ignore_index=True, sort=False), matched


def main() -> None:
    signature = read_signature()
    homology_map, homology_table = strict_human_mouse_map()
    prior = signed_collectri(PRIOR_PATH)

    datasets = {
        "GSE114059": ("human", "patient_derived_crc_organoid", *load_gse114059()),
        "GSE67186": ("mouse", "in_vivo_colon_polyp", *load_gse67186(homology_map)),
        "GSE130822": ("mouse", "colonic_stem_cells", *load_gse130822(homology_map)),
        "GSE171910": ("human", "conditional_crc_wnt_models", *load_gse171910()),
    }

    all_samples = []
    all_units = []
    all_summaries = []
    all_coverage = []
    all_matched_tests = []
    all_matched_null = []
    ulm_samples = []
    ulm_coverage = []
    ulm_units = []
    ulm_summaries = []

    for dataset, (species, model_system, expression, metadata, comparisons) in datasets.items():
        scores, z, coverage = score_expression(expression, signature, dataset)
        merged, units, summaries = comparison_effects(dataset, species, model_system, scores, metadata, comparisons)
        merged.insert(0, "dataset", dataset)
        merged.insert(1, "species", species)
        merged.insert(2, "model_system", model_system)
        all_samples.append(merged)
        all_units.append(units)
        all_summaries.append(summaries)
        all_coverage.append(coverage)
        tests, null = expression_matched_test(dataset, expression, z, signature, metadata, comparisons)
        all_matched_tests.append(tests)
        all_matched_null.append(null)
        if species == "human":
            activity, activity_coverage = ulm_activity(z, prior, dataset)
            units_activity, summaries_activity = ulm_comparison_effects(dataset, activity, metadata, comparisons)
            ulm_samples.append(activity)
            ulm_coverage.append(activity_coverage)
            ulm_units.append(units_activity)
            ulm_summaries.append(summaries_activity)

    # Add ULM calibration for the primary APC organoids.
    gse125_expression, gse125_metadata, gse125_comparisons = load_gse125_expression()
    _, gse125_z, gse125_coverage = score_expression(gse125_expression, signature, "GSE125472")
    gse125_activity, gse125_activity_coverage = ulm_activity(gse125_z, prior, "GSE125472")
    gse125_ulm_units, gse125_ulm_summary = ulm_comparison_effects(
        "GSE125472", gse125_activity, gse125_metadata, gse125_comparisons
    )
    ulm_samples.append(gse125_activity)
    ulm_coverage.append(gse125_activity_coverage)
    ulm_units.append(gse125_ulm_units)
    ulm_summaries.append(gse125_ulm_summary)

    # Add ULM calibration for TCF7L2 knockout, standardized within cell line.
    for cell_line, (expression, metadata) in load_gse135_expression().items():
        _, z, coverage = score_expression(expression, signature, f"GSE135328_{cell_line}")
        activity, activity_coverage = ulm_activity(z, prior, f"GSE135328_{cell_line}")
        units_activity, summaries_activity = gse135_ulm_effects(activity, metadata, cell_line)
        ulm_samples.append(activity)
        ulm_coverage.append(activity_coverage)
        ulm_units.append(units_activity)
        ulm_summaries.append(summaries_activity)
        all_coverage.append(coverage)
    all_coverage.append(gse125_coverage)

    existing_units, existing_summaries, existing_matched = append_existing_perturbations()
    all_units.append(existing_units)
    all_summaries.append(existing_summaries)
    all_matched_tests.append(existing_matched)

    spatial_spots, spatial_summary, spatial_units, spatial_tests = spatial_locked_route(signature)

    sample_scores = pd.concat(all_samples, ignore_index=True, sort=False)
    unit_effects = pd.concat(all_units, ignore_index=True, sort=False)
    effect_summary = pd.concat(all_summaries, ignore_index=True, sort=False)
    coverage = pd.concat(all_coverage, ignore_index=True, sort=False)
    coverage["coverage_fraction"] = coverage["n_present"] / coverage["n_expected"]
    route_coverage = coverage[coverage["feature"].isin(["route_up", "route_down"])].pivot(
        index="dataset", columns="feature", values="coverage_fraction"
    )
    route_coverage["route_reportable"] = route_coverage[["route_up", "route_down"]].ge(0.80).all(axis=1)
    route_reportable = route_coverage["route_reportable"].rename_axis("dataset").reset_index()
    unit_effects = unit_effects.merge(route_reportable, on="dataset", how="left")
    effect_summary = effect_summary.merge(route_reportable, on="dataset", how="left")
    matched_tests = pd.concat(all_matched_tests, ignore_index=True, sort=False)
    matched_tests = matched_tests.merge(route_reportable, on="dataset", how="left")
    matched_null = pd.concat(all_matched_null, ignore_index=True, sort=False)
    ulm_sample_activity = pd.concat(ulm_samples, ignore_index=True, sort=False)
    ulm_target_coverage = pd.concat(ulm_coverage, ignore_index=True, sort=False)
    ulm_unit_effects = pd.concat(ulm_units, ignore_index=True, sort=False)
    ulm_effect_summary = pd.concat(ulm_summaries, ignore_index=True, sort=False)

    virtual_panel = pd.read_csv(
        ROOT / "results" / "perturbation_validation_locked_route" / "virtual_tf_knockout_panel.tsv", sep="\t"
    )
    signed_mean_summary = pd.read_csv(
        ROOT / "results" / "perturbation_validation_locked_route" / "gse135328_clone_contrasts.tsv", sep="\t"
    )
    signed_mean_summary = signed_mean_summary[
        signed_mean_summary["feature"].eq("tf_activity__TCF7L2") & signed_mean_summary["genotype"].eq("KO")
    ].copy()

    # Post hoc sensitivity only: remove the observed OFF-DOX minus ON-DOX shift
    # in shRenilla from each Apc-restoration contrast. The frozen raw contrasts
    # remain primary and the non-null control is retained.
    dox_rows = []
    gse671_summary = effect_summary[effect_summary["dataset"].eq("GSE67186")]
    for feature, feature_rows in gse671_summary.groupby("feature", observed=True):
        control = feature_rows[feature_rows["comparison"].eq("doxycycline_control_shRenilla")]
        if control.empty:
            continue
        control_effect = float(control.iloc[0]["mean_difference"])
        for comparison in ["apc_restoration_shApc", "apc_restoration_shApc_Kras"]:
            treatment = feature_rows[feature_rows["comparison"].eq(comparison)]
            if treatment.empty:
                continue
            treatment_effect = float(treatment.iloc[0]["mean_difference"])
            dox_rows.append(
                {
                    "dataset": "GSE67186",
                    "comparison": comparison,
                    "feature": feature,
                    "raw_off_minus_on": treatment_effect,
                    "shRenilla_off_minus_on": control_effect,
                    "control_adjusted_difference_in_differences": treatment_effect - control_effect,
                    "analysis_role": "post_hoc_control_adjusted_sensitivity",
                }
            )
    dox_control_sensitivity = pd.DataFrame(dox_rows)

    evidence_rows = []
    route_summary = effect_summary[effect_summary["feature"].eq("route_score")]
    for _, row in route_summary.iterrows():
        expected = row.get("expected_direction", np.nan)
        specificity = matched_tests[
            matched_tests["dataset"].eq(row["dataset"]) & matched_tests["comparison"].eq(row["comparison"])
        ]
        specificity_p = (
            float(specificity.iloc[0]["p_expression_matched_one_sided"])
            if not specificity.empty and pd.notna(specificity.iloc[0]["p_expression_matched_one_sided"])
            else np.nan
        )
        reportable = bool(row.get("route_reportable", False))
        if pd.isna(expected) or expected == 0:
            status = "control_nonzero" if pd.notna(specificity_p) and specificity_p <= 0.05 else "control_null"
        elif not reportable:
            status = "exploratory_low_coverage"
        elif expected * row["mean_difference"] > 0:
            status = "supportive_specific" if pd.notna(specificity_p) and specificity_p <= 0.05 else "supportive_direction_only"
        else:
            status = "discordant"
        evidence_rows.append(
            {
                "layer": "empirical_perturbation",
                "dataset": row["dataset"],
                "comparison": row["comparison"],
                "effect": row["mean_difference"],
                "expected_direction": expected,
                "status": status,
                "n_units": row["n_units"],
                "coverage_reportable": reportable,
                "specificity_p": specificity_p,
            }
        )
    for _, row in spatial_tests[spatial_tests["feature"].isin(["route_score", "route_residual_prolif_epithelial"])].iterrows():
        evidence_rows.append(
            {
                "layer": "spatial_recapitulation",
                "dataset": "Zenodo7760264",
                "comparison": f"{row['comparison']}__{row['feature']}",
                "effect": row["median_difference"],
                "expected_direction": 1,
                "status": (
                    "supportive_paired"
                    if row["median_difference"] > 0 and row["p_wilcoxon"] <= 0.05
                    else "supportive_direction_only"
                    if row["median_difference"] > 0
                    else "discordant"
                ),
                "n_units": row["n_sections"],
                "coverage_reportable": True,
                "specificity_p": row["p_wilcoxon"],
            }
        )
    evidence_matrix = pd.DataFrame(evidence_rows)

    tables = {
        "perturbation_sample_scores.tsv": sample_scores,
        "perturbation_unit_effects.tsv": unit_effects,
        "perturbation_effect_summary.tsv": effect_summary,
        "feature_coverage.tsv": coverage,
        "expression_matched_signature_tests.tsv": matched_tests,
        "collectri_ulm_sample_activity.tsv": ulm_sample_activity,
        "collectri_ulm_target_coverage.tsv": ulm_target_coverage,
        "collectri_ulm_unit_effects.tsv": ulm_unit_effects,
        "collectri_ulm_effect_summary.tsv": ulm_effect_summary,
        "collectri_signed_mean_tcf7l2_clone_effects.tsv": signed_mean_summary,
        "virtual_tf_knockout_panel.tsv": virtual_panel,
        "spatial_locked_route_spot_scores.tsv.gz": spatial_spots,
        "spatial_locked_route_pathology_summary.tsv": spatial_summary,
        "spatial_locked_route_section_effects.tsv": spatial_units,
        "spatial_locked_route_tests.tsv": spatial_tests,
        "human_mouse_one_to_one_homology.tsv": homology_table,
        "gse67186_doxycycline_control_sensitivity.tsv": dox_control_sensitivity,
        "evidence_closure_matrix.tsv": evidence_matrix,
    }
    for filename, frame in tables.items():
        compression = "gzip" if filename.endswith(".gz") else None
        frame.to_csv(OUT / filename, sep="\t", index=False, compression=compression)
    matched_null.to_csv(OUT / "expression_matched_signature_null.tsv.gz", sep="\t", index=False, compression="gzip")

    existing_base = ROOT / "results" / "perturbation_validation_locked_route"
    spatial_source_paths = []
    for sample_dir in sorted(
        path for path in SPATIAL_DIR.iterdir()
        if path.is_dir() and (path / "filtered_feature_bc_matrix.h5").exists()
    ):
        sample_id = sample_dir.name
        spatial_source_paths.extend(
            [
                sample_dir / "filtered_feature_bc_matrix.h5",
                sample_dir / "spatial" / "tissue_positions_list.csv",
                SPATIAL_ANNOT_DIR / f"Pathologist_Annotations_{sample_id}.csv",
            ]
        )
    source_paths = [
        SIGNATURE_PATH,
        MGI_PATH,
        PRIOR_PATH,
        GSE114_DIR / "GSE114059_series_matrix.txt.gz",
        GSE114_DIR / "GPL10558.annot.gz",
        GSE671_DIR / "GSE67186_counts_scaled_DESeq.xls.gz",
        GSE130_DIR / "GSE130822_Read_Counts_All_Cells.txt.gz",
        GSE171_DIR / "GSE171910_RNA_seq_results.xlsx",
        GSE125_PATH,
        GSE135_DIR / "GSE135328_count_HCT116.txt.gz",
        GSE135_DIR / "GSE135328_count_HT29.txt.gz",
        GENE_INFO_PATH,
        existing_base / "gse125472_donor_contrasts.tsv",
        existing_base / "gse125472_contrast_summary.tsv",
        existing_base / "gse125472_matched_signature_tests.tsv",
        existing_base / "gse135328_clone_contrasts.tsv",
        existing_base / "gse135328_matched_signature_tests.tsv",
        existing_base / "virtual_tf_knockout_panel.tsv",
    ] + spatial_source_paths
    source_manifest = pd.DataFrame(
        [{"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in source_paths]
    )
    source_manifest.to_csv(OUT / "source_file_manifest.tsv", sep="\t", index=False)

    route_view = evidence_matrix[evidence_matrix["layer"].eq("empirical_perturbation")][
        ["dataset", "comparison", "n_units", "effect", "expected_direction", "status", "specificity_p"]
    ]
    spatial_view = spatial_tests[spatial_tests["feature"].isin(["route_score", "route_residual_prolif_epithelial"])]
    summary = "# Computational closure validation\n\n## Empirical perturbation route effects\n\n```text\n"
    summary += route_view.to_string(index=False)
    summary += "\n```\n\n## Spatial locked-route effects\n\n```text\n"
    summary += spatial_view.to_string(index=False)
    summary += "\n```\n\nAll contexts are retained, including discordant and low-information results.\n"
    (OUT / "summary.md").write_text(summary, encoding="utf-8")

    manifest = {
        "analysis": "computational_closure_validation",
        "date": "2026-08-08",
        "random_seed": SEED,
        "n_permutations": N_PERMUTATIONS,
        "n_bootstraps": N_BOOTSTRAPS,
        "signature_genes": 100,
        "drawing_backend": "R",
        "claim_boundary": "multi-layer computational validation; not local experimental validation",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote computational closure validation to {OUT}")


if __name__ == "__main__":
    main()
