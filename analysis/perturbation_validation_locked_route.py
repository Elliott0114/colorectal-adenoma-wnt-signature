#!/usr/bin/env python3
"""Validate the frozen adenoma route in independent genetic perturbation data.

The primary analysis uses donor-matched APC-knockout human colon organoids
(GSE125472). TCF7L2-knockout CRC cell lines (GSE135328) are retained as a
context-dependent stress test. A signed CollecTRI edge-deletion screen provides
topology-based virtual perturbation support without reselecting route genes.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu, wilcoxon


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "perturbation_validation_locked_route"
OUT.mkdir(parents=True, exist_ok=True)

SIGNATURE_PATH = ROOT / "results" / "route_signature_locked" / "discovery_locked_signature_genes.tsv"
GSE125_PATH = (
    ROOT
    / "data_sources"
    / "GSE125472_apc_ko_organoids"
    / "GSE125472_20181220_Results_Wnt_Signature_ALL.txt.gz"
)
GSE135_DIR = ROOT / "data_sources" / "GSE135328_tcf7l2_ko_crc"
GENE_INFO_PATH = ROOT / "data_sources" / "GSE117606" / "Homo_sapiens.gene_info.gz"
COLLECTRI_PANEL_PATH = (
    ROOT
    / "data_sources"
    / "regulatory_priors"
    / "collectri_tcf_ascl2_wnt_controls_2026-08-08.tsv"
)
COLLECTRI_ROUTE_PATH = (
    ROOT
    / "data_sources"
    / "regulatory_priors"
    / "collectri_locked_route_targets_2026-08-08.tsv"
)
CHEN_DIR = ROOT / "results" / "route_signature"

RANDOM_SEED = 20260808
N_PERMUTATIONS = 10_000
N_BOOTSTRAPS = 10_000

WNT_STEM_GENES = ["LGR5", "ASCL2", "OLFM4", "AXIN2", "SOX9", "EPHB2", "SMOC2"]
PROLIFERATION_GENES = ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "TYMS", "UBE2C", "CENPF"]
CANDIDATE_GENES = ["OLFM4", "FABP1", "ASCL2", "AXIN2", "LGR5", "SOX9"]
TF_PANEL = ["TCF7L2", "ASCL2", "MYC", "SOX9", "HNF4A", "KLF4"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bool_col(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def read_signature() -> pd.DataFrame:
    signature = pd.read_csv(SIGNATURE_PATH, sep="\t")
    if len(signature) != 100:
        raise ValueError(f"Expected 100 locked genes, found {len(signature)}")
    if signature["validation_used_for_selection"].astype(str).str.lower().ne("false").any():
        raise ValueError("Locked signature contains validation-selected genes")
    signature["route_weight"] = np.where(
        signature["signature_direction"].eq("adenoma_up"), 1.0, -1.0
    )
    return signature[["gene", "signature_direction", "rank_within_direction", "route_weight"]].copy()


def aggregate_symbol_counts(counts: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    frame = counts.copy()
    frame.insert(0, "gene", symbols.astype(str).str.strip().to_numpy())
    frame = frame[frame["gene"].notna() & frame["gene"].ne("") & frame["gene"].ne("nan")]
    return frame.groupby("gene", sort=True, observed=True)[counts.columns].sum()


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    library = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(library, axis=1) * 1_000_000 + 1)


def gene_zscores(expression: pd.DataFrame) -> pd.DataFrame:
    means = expression.mean(axis=1)
    sds = expression.std(axis=1, ddof=0).replace(0, np.nan)
    return expression.sub(means, axis=0).div(sds, axis=0)


def mean_present(z: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, int]:
    present = [gene for gene in genes if gene in z.index]
    if not present:
        return pd.Series(np.nan, index=z.columns), 0
    return z.loc[present].mean(axis=0), len(present)


def score_expression(
    expression: pd.DataFrame,
    signature: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    z = gene_zscores(expression)
    up = signature.loc[signature["route_weight"].eq(1), "gene"].tolist()
    down = signature.loc[signature["route_weight"].eq(-1), "gene"].tolist()
    up_score, n_up = mean_present(z, up)
    down_score, n_down = mean_present(z, down)
    wnt_score, n_wnt = mean_present(z, WNT_STEM_GENES)
    proliferation_score, n_proliferation = mean_present(z, PROLIFERATION_GENES)
    scores = pd.DataFrame(
        {
            "sample_id": expression.columns,
            "route_up": up_score.reindex(expression.columns).to_numpy(),
            "route_down": down_score.reindex(expression.columns).to_numpy(),
            "route_score": (up_score - down_score).reindex(expression.columns).to_numpy(),
            "wnt_stem": wnt_score.reindex(expression.columns).to_numpy(),
            "proliferation_control": proliferation_score.reindex(expression.columns).to_numpy(),
        }
    )
    for gene in CANDIDATE_GENES + ["TCF7L2", "MYC", "HNF4A", "KLF4"]:
        scores[f"gene_z__{gene}"] = z.loc[gene].reindex(expression.columns).to_numpy() if gene in z.index else np.nan
    coverage = pd.DataFrame(
        [
            {"feature": "route_up", "n_expected": 50, "n_present": n_up},
            {"feature": "route_down", "n_expected": 50, "n_present": n_down},
            {"feature": "wnt_stem", "n_expected": len(WNT_STEM_GENES), "n_present": n_wnt},
            {
                "feature": "proliferation_control",
                "n_expected": len(PROLIFERATION_GENES),
                "n_present": n_proliferation,
            },
        ]
    )
    return scores, coverage, z


def signed_collectri(path: Path) -> pd.DataFrame:
    prior = pd.read_csv(path, sep="\t")
    stimulation = bool_col(prior["consensus_stimulation"])
    inhibition = bool_col(prior["consensus_inhibition"])
    prior = prior.loc[stimulation ^ inhibition].copy()
    prior["edge_sign"] = np.where(stimulation.loc[prior.index], 1.0, -1.0)
    prior = prior[["source_genesymbol", "target_genesymbol", "edge_sign"]].dropna()
    grouped = prior.groupby(["source_genesymbol", "target_genesymbol"], observed=True)["edge_sign"].agg(
        lambda x: x.iloc[0] if x.nunique() == 1 else np.nan
    )
    return grouped.dropna().reset_index()


def regulon_activity(z: pd.DataFrame, prior: pd.DataFrame, min_targets: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    coverage = []
    for tf in TF_PANEL:
        edges = prior[prior["source_genesymbol"].eq(tf)].drop_duplicates("target_genesymbol")
        edges = edges[edges["target_genesymbol"].isin(z.index)]
        coverage.append({"tf": tf, "n_targets": len(edges), "targets": ",".join(sorted(edges["target_genesymbol"]))})
        if len(edges) < min_targets:
            continue
        target_z = z.loc[edges["target_genesymbol"]].copy()
        target_z.index = edges["target_genesymbol"].to_numpy()
        signs = edges.set_index("target_genesymbol")["edge_sign"].reindex(target_z.index)
        activity = target_z.mul(signs, axis=0).sum(axis=0) / math.sqrt(len(edges))
        rows.extend({"sample_id": sample, "tf": tf, "activity": value, "n_targets": len(edges)} for sample, value in activity.items())
    return pd.DataFrame(rows), pd.DataFrame(coverage)


def bootstrap_mean(values: np.ndarray, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(N_BOOTSTRAPS, len(values)), replace=True).mean(axis=1)
    return tuple(np.quantile(draws, [0.025, 0.975]))


def exact_sign_p(values: np.ndarray) -> float:
    nonzero = values[values != 0]
    if not len(nonzero):
        return np.nan
    return float(binomtest(int((nonzero > 0).sum()), len(nonzero), 0.5).pvalue)


def summarize_differences(values: np.ndarray, seed: int) -> dict[str, object]:
    values = np.asarray(values, dtype=float)
    ci_low, ci_high = bootstrap_mean(values, seed)
    sd = float(values.std(ddof=1)) if len(values) > 1 else np.nan
    return {
        "n_units": len(values),
        "mean_difference": float(values.mean()),
        "median_difference": float(np.median(values)),
        "sd_difference": sd,
        "paired_dz": float(values.mean() / sd) if np.isfinite(sd) and sd > 0 else np.nan,
        "bootstrap_mean_ci_low": ci_low,
        "bootstrap_mean_ci_high": ci_high,
        "min_difference": float(values.min()),
        "max_difference": float(values.max()),
        "n_positive": int((values > 0).sum()),
        "n_negative": int((values < 0).sum()),
        "p_exact_sign_two_sided": exact_sign_p(values),
    }


def gse125_metadata(sample_columns: list[str]) -> pd.DataFrame:
    rows = []
    pattern = re.compile(r"^Donor(?P<donor>\d+)_(?P<genotype>WT|APC)_(?P<medium>w|wo)_Wnt$")
    for sample in sample_columns:
        match = pattern.match(sample)
        if not match:
            raise ValueError(f"Unrecognized GSE125472 sample: {sample}")
        rows.append(
            {
                "sample_id": sample,
                "donor_id": f"Donor{match.group('donor')}",
                "genotype": match.group("genotype"),
                "wnt_rspo": "with" if match.group("medium") == "w" else "without",
            }
        )
    return pd.DataFrame(rows)


GSE125_COMPARISONS = [
    ("APC_vs_WT_with_Wnt", ("APC", "with"), ("WT", "with"), 1),
    ("APC_vs_WT_without_Wnt", ("APC", "without"), ("WT", "without"), 1),
    ("WT_withdrawal", ("WT", "without"), ("WT", "with"), -1),
    ("APC_withdrawal", ("APC", "without"), ("APC", "with"), 0),
]


def gse125_expected_direction(feature: str, comparison: str) -> int:
    """Return the preregistered direction for each program component.

    The down arm is stored as its raw mean expression, so its expected direction
    is opposite to the composite route score. Proliferation is a control rather
    than a directional endpoint.
    """
    if comparison in {"APC_vs_WT_with_Wnt", "APC_vs_WT_without_Wnt", "genotype_by_Wnt_interaction"}:
        return {"route_score": 1, "route_up": 1, "route_down": -1, "wnt_stem": 1}.get(feature, 0)
    if comparison == "WT_withdrawal":
        return {"route_score": -1, "route_up": -1, "route_down": 1, "wnt_stem": -1}.get(feature, 0)
    return 0


def gse125_contrasts(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    features = ["route_score", "route_up", "route_down", "wnt_stem", "proliferation_control"]
    donor_rows = []
    summaries = []
    loo_rows = []
    for feature in features:
        wide = scores.pivot(index="donor_id", columns=["genotype", "wnt_rspo"], values=feature)
        for comparison_id, target, reference, _ in GSE125_COMPARISONS:
            expected_direction = gse125_expected_direction(feature, comparison_id)
            diff = (wide[target] - wide[reference]).dropna()
            for donor, value in diff.items():
                donor_rows.append(
                    {
                        "feature": feature,
                        "comparison": comparison_id,
                        "expected_direction": expected_direction,
                        "donor_id": donor,
                        "difference": value,
                    }
                )
            summary = {
                "feature": feature,
                "comparison": comparison_id,
                "expected_direction": expected_direction,
                **summarize_differences(diff.to_numpy(), RANDOM_SEED + len(summaries)),
            }
            if expected_direction:
                summary["n_expected_direction"] = int((expected_direction * diff.to_numpy() > 0).sum())
            else:
                summary["n_expected_direction"] = np.nan
            summaries.append(summary)
            for omitted in diff.index:
                retained = diff.drop(omitted)
                loo_rows.append(
                    {
                        "feature": feature,
                        "comparison": comparison_id,
                        "omitted_donor": omitted,
                        "n_retained": len(retained),
                        "mean_difference": retained.mean(),
                        "median_difference": retained.median(),
                    }
                )
    # Difference-in-differences: withdrawal response in APC minus WT.
    for feature in features:
        wide = scores.pivot(index="donor_id", columns=["genotype", "wnt_rspo"], values=feature)
        diff = ((wide[("APC", "without")] - wide[("APC", "with")]) - (wide[("WT", "without")] - wide[("WT", "with")])).dropna()
        expected_direction = gse125_expected_direction(feature, "genotype_by_Wnt_interaction")
        for donor, value in diff.items():
            donor_rows.append(
                {
                    "feature": feature,
                    "comparison": "genotype_by_Wnt_interaction",
                    "expected_direction": expected_direction,
                    "donor_id": donor,
                    "difference": value,
                }
            )
        summaries.append(
            {
                "feature": feature,
                "comparison": "genotype_by_Wnt_interaction",
                "expected_direction": expected_direction,
                **summarize_differences(diff.to_numpy(), RANDOM_SEED + len(summaries)),
                "n_expected_direction": (
                    int((expected_direction * diff.to_numpy() > 0).sum()) if expected_direction else np.nan
                ),
            }
        )
        for omitted in diff.index:
            retained = diff.drop(omitted)
            loo_rows.append(
                {
                    "feature": feature,
                    "comparison": "genotype_by_Wnt_interaction",
                    "omitted_donor": omitted,
                    "n_retained": len(retained),
                    "mean_difference": retained.mean(),
                    "median_difference": retained.median(),
                }
            )
    return pd.DataFrame(donor_rows), pd.DataFrame(summaries), pd.DataFrame(loo_rows)


def donor_contrast_from_gene_score(
    score: pd.Series,
    metadata: pd.DataFrame,
    target: tuple[str, str],
    reference: tuple[str, str],
) -> float:
    frame = metadata.copy()
    frame["score"] = frame["sample_id"].map(score)
    wide = frame.pivot(index="donor_id", columns=["genotype", "wnt_rspo"], values="score")
    return float((wide[target] - wide[reference]).mean())


def matched_signature_permutations(
    expression: pd.DataFrame,
    z: pd.DataFrame,
    metadata: pd.DataFrame,
    signature: pd.DataFrame,
    comparison_map: list[tuple[str, tuple[str, str], tuple[str, str], int]],
    dataset: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    signature_present = signature[signature["gene"].isin(z.index)].copy()
    mean_expression = expression.mean(axis=1)
    bins = pd.qcut(mean_expression.rank(method="first"), q=10, labels=False)
    signature_present["expression_bin"] = signature_present["gene"].map(bins)
    candidate_pool = pd.DataFrame({"gene": z.index, "expression_bin": bins.reindex(z.index).to_numpy()})
    candidate_pool = candidate_pool[~candidate_pool["gene"].isin(signature_present["gene"])]
    grouped_candidates = {
        int(group): frame["gene"].to_numpy()
        for group, frame in candidate_pool.groupby("expression_bin", observed=True)
    }
    rng = np.random.default_rng(RANDOM_SEED + sum(map(ord, dataset)))
    observed_up = signature_present.loc[signature_present["route_weight"].eq(1), "gene"].tolist()
    observed_down = signature_present.loc[signature_present["route_weight"].eq(-1), "gene"].tolist()
    observed_score = z.loc[observed_up].mean(axis=0) - z.loc[observed_down].mean(axis=0)
    observed = {
        comparison: donor_contrast_from_gene_score(observed_score, metadata, target, reference)
        for comparison, target, reference, _ in comparison_map
    }
    null_rows = []
    members_by_bin = {
        int(group): frame.sort_values(["route_weight", "gene"])
        for group, frame in signature_present.groupby("expression_bin", observed=True)
    }
    for permutation_id in range(1, N_PERMUTATIONS + 1):
        up_genes = []
        down_genes = []
        for expression_bin, members in members_by_bin.items():
            pool = grouped_candidates[expression_bin]
            sampled = rng.choice(pool, size=len(members), replace=False)
            directions = members["route_weight"].to_numpy()
            up_genes.extend(sampled[directions == 1])
            down_genes.extend(sampled[directions == -1])
        random_score = z.loc[up_genes].mean(axis=0) - z.loc[down_genes].mean(axis=0)
        for comparison, target, reference, _ in comparison_map:
            null_rows.append(
                {
                    "dataset": dataset,
                    "permutation_id": permutation_id,
                    "comparison": comparison,
                    "null_mean_difference": donor_contrast_from_gene_score(random_score, metadata, target, reference),
                }
            )
    null = pd.DataFrame(null_rows)
    tests = []
    for comparison, _, _, expected_direction in comparison_map:
        values = null.loc[null["comparison"].eq(comparison), "null_mean_difference"].to_numpy()
        obs = observed[comparison]
        if expected_direction > 0:
            p = (1 + np.sum(values >= obs)) / (len(values) + 1)
        elif expected_direction < 0:
            p = (1 + np.sum(values <= obs)) / (len(values) + 1)
        else:
            p = (1 + np.sum(np.abs(values) >= abs(obs))) / (len(values) + 1)
        tests.append(
            {
                "dataset": dataset,
                "comparison": comparison,
                "expected_direction": expected_direction,
                "observed_mean_difference": obs,
                "null_mean": values.mean(),
                "null_sd": values.std(ddof=1),
                "null_ci_low": np.quantile(values, 0.025),
                "null_ci_high": np.quantile(values, 0.975),
                "p_expression_matched_one_sided": p,
                "n_permutations": len(values),
                "n_up": len(observed_up),
                "n_down": len(observed_down),
            }
        )
    return null, pd.DataFrame(tests)


def analyze_gse125(signature: pd.DataFrame, prior_panel: pd.DataFrame) -> dict[str, pd.DataFrame]:
    raw = pd.read_csv(GSE125_PATH, sep="\t", low_memory=False)
    sample_columns = [column for column in raw.columns if re.match(r"^Donor\d+_(WT|APC)_(w|wo)_Wnt$", column)]
    counts = raw[sample_columns].apply(pd.to_numeric, errors="raise")
    symbols = raw["hgnc_symbol"].fillna(raw["symbol"])
    counts = aggregate_symbol_counts(counts, symbols)
    expression = log2_cpm(counts)
    metadata = gse125_metadata(sample_columns)
    scores, coverage, z = score_expression(expression, signature)
    scores = metadata.merge(scores, on="sample_id", validate="one_to_one")
    activity, activity_coverage = regulon_activity(z, prior_panel)
    activity.insert(0, "dataset", "GSE125472")
    scores = scores.merge(
        activity.pivot(index="sample_id", columns="tf", values="activity").add_prefix("tf_activity__").reset_index(),
        on="sample_id",
        how="left",
    )
    donor_contrasts, contrast_summary, loo = gse125_contrasts(scores)
    null, matched_tests = matched_signature_permutations(
        expression,
        z,
        metadata,
        signature,
        GSE125_COMPARISONS,
        "GSE125472",
    )
    gene_rows = []
    gene_effect_rows = []
    gene_route_direction = dict(
        zip(
            signature["gene"],
            signature["signature_direction"].map({"adenoma_up": 1, "adenoma_down": -1}),
            strict=True,
        )
    )
    # SOX9 is a prespecified context marker but is not one of the 100 locked genes.
    gene_route_direction["SOX9"] = 1
    for gene in sorted(set(signature["gene"]) | set(CANDIDATE_GENES) | {"TCF7L2", "MYC", "HNF4A", "KLF4"}):
        if gene not in expression.index:
            continue
        frame = metadata.copy()
        frame["expression"] = frame["sample_id"].map(expression.loc[gene])
        wide = frame.pivot(index="donor_id", columns=["genotype", "wnt_rspo"], values="expression")
        for comparison, target, reference, route_direction in GSE125_COMPARISONS:
            expected_direction = int(gene_route_direction.get(gene, 0) * route_direction)
            diff = (wide[target] - wide[reference]).dropna()
            for donor, value in diff.items():
                gene_rows.append(
                    {
                        "gene": gene,
                        "comparison": comparison,
                        "donor_id": donor,
                        "log2cpm_difference": value,
                    }
                )
            gene_effect_rows.append(
                {
                    "gene": gene,
                    "comparison": comparison,
                    "expected_direction": expected_direction,
                    **summarize_differences(diff.to_numpy(), RANDOM_SEED + len(gene_effect_rows)),
                }
            )
    coverage.insert(0, "dataset", "GSE125472")
    activity_coverage.insert(0, "dataset", "GSE125472")
    return {
        "sample_scores": scores,
        "coverage": coverage,
        "activity": activity,
        "activity_coverage": activity_coverage,
        "donor_contrasts": donor_contrasts,
        "contrast_summary": contrast_summary,
        "leave_one_donor_out": loo,
        "matched_null": null,
        "matched_tests": matched_tests,
        "gene_donor_contrasts": pd.DataFrame(gene_rows),
        "gene_effects": pd.DataFrame(gene_effect_rows),
        "expression": expression,
        "z": z,
    }


def ensembl_to_symbol() -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(GENE_INFO_PATH, "rt", encoding="utf-8") as handle:
        header = next(handle).rstrip("\n").lstrip("#").split("\t")
        idx_symbol = header.index("Symbol")
        idx_xref = header.index("dbXrefs")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            symbol = fields[idx_symbol]
            for xref in fields[idx_xref].split("|"):
                if xref.startswith("Ensembl:ENSG"):
                    mapping[xref.split(":", 1)[1].split(".", 1)[0]] = symbol
    return mapping


def read_gse135_counts(path: Path, mapping: dict[str, str]) -> pd.DataFrame:
    counts = pd.read_csv(path, sep="\t", index_col=0)
    symbols = pd.Series(counts.index.astype(str).str.split(".").str[0].map(mapping), index=counts.index)
    return aggregate_symbol_counts(counts, symbols)


def gse135_metadata(cell_line: str, columns: list[str]) -> pd.DataFrame:
    genotype_map = {
        "HCT116": {"3": "WT", "18": "KO"},
        "HT29": {"56": "WT", "62": "Het", "57": "KO", "83": "KO", "86": "KO"},
    }
    rows = []
    pattern = re.compile(r"^(?:HCT|HT29)_(?P<clone>\d+)_(?P<replicate>[AB])$")
    for sample in columns:
        match = pattern.match(sample)
        if not match:
            raise ValueError(f"Unrecognized GSE135328 sample: {sample}")
        clone = match.group("clone")
        rows.append(
            {
                "sample_id": sample,
                "cell_line": cell_line,
                "clone_id": f"{cell_line}_{clone}",
                "clone_number": clone,
                "genotype": genotype_map[cell_line][clone],
                "replicate": match.group("replicate"),
            }
        )
    return pd.DataFrame(rows)


def simple_matched_test(
    expression: pd.DataFrame,
    z: pd.DataFrame,
    signature: pd.DataFrame,
    metadata: pd.DataFrame,
    cell_line: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    clone_meta = metadata[["sample_id", "clone_id", "genotype"]].drop_duplicates()
    signature_present = signature[signature["gene"].isin(z.index)].copy()
    mean_expression = expression.mean(axis=1)
    bins = pd.qcut(mean_expression.rank(method="first"), 10, labels=False)
    signature_present["expression_bin"] = signature_present["gene"].map(bins)
    pool = pd.DataFrame({"gene": z.index, "expression_bin": bins.reindex(z.index).to_numpy()})
    pool = pool[~pool["gene"].isin(signature_present["gene"])]
    pool_by_bin = {int(k): v["gene"].to_numpy() for k, v in pool.groupby("expression_bin", observed=True)}
    members_by_bin = {int(k): v.sort_values(["route_weight", "gene"]) for k, v in signature_present.groupby("expression_bin", observed=True)}
    rng = np.random.default_rng(RANDOM_SEED + sum(map(ord, cell_line)))

    def clone_effect(score: pd.Series) -> float:
        frame = clone_meta.copy()
        frame["score"] = frame["sample_id"].map(score)
        clone_scores = frame.groupby(["clone_id", "genotype"], observed=True)["score"].mean().reset_index()
        return float(clone_scores.loc[clone_scores["genotype"].eq("KO"), "score"].mean() - clone_scores.loc[clone_scores["genotype"].eq("WT"), "score"].mean())

    up = signature_present.loc[signature_present["route_weight"].eq(1), "gene"].tolist()
    down = signature_present.loc[signature_present["route_weight"].eq(-1), "gene"].tolist()
    observed = clone_effect(z.loc[up].mean(axis=0) - z.loc[down].mean(axis=0))
    values = np.empty(N_PERMUTATIONS, dtype=float)
    for i in range(N_PERMUTATIONS):
        random_up = []
        random_down = []
        for expression_bin, members in members_by_bin.items():
            sampled = rng.choice(pool_by_bin[expression_bin], size=len(members), replace=False)
            directions = members["route_weight"].to_numpy()
            random_up.extend(sampled[directions == 1])
            random_down.extend(sampled[directions == -1])
        values[i] = clone_effect(z.loc[random_up].mean(axis=0) - z.loc[random_down].mean(axis=0))
    null = pd.DataFrame(
        {"dataset": "GSE135328", "cell_line": cell_line, "permutation_id": np.arange(1, N_PERMUTATIONS + 1), "null_ko_minus_wt": values}
    )
    test = pd.DataFrame(
        [
            {
                "dataset": "GSE135328",
                "cell_line": cell_line,
                "observed_ko_minus_wt": observed,
                "null_mean": values.mean(),
                "null_sd": values.std(ddof=1),
                "null_ci_low": np.quantile(values, 0.025),
                "null_ci_high": np.quantile(values, 0.975),
                "p_expression_matched_lower_tail": (1 + np.sum(values <= observed)) / (N_PERMUTATIONS + 1),
                "n_permutations": N_PERMUTATIONS,
                "n_up": len(up),
                "n_down": len(down),
            }
        ]
    )
    return null, test


def analyze_gse135(signature: pd.DataFrame, prior_panel: pd.DataFrame) -> dict[str, pd.DataFrame]:
    mapping = ensembl_to_symbol()
    data = []
    clone_scores_all = []
    clone_contrasts = []
    gene_effects_all = []
    null_all = []
    tests_all = []
    activities_all = []
    coverage_all = []
    expression_by_cell_line: dict[str, pd.DataFrame] = {}
    for cell_line, filename in [("HCT116", "GSE135328_count_HCT116.txt.gz"), ("HT29", "GSE135328_count_HT29.txt.gz")]:
        counts = read_gse135_counts(GSE135_DIR / filename, mapping)
        expression = log2_cpm(counts)
        expression_by_cell_line[cell_line] = expression
        metadata = gse135_metadata(cell_line, expression.columns.tolist())
        scores, coverage, z = score_expression(expression, signature)
        scores = metadata.merge(scores, on="sample_id", validate="one_to_one")
        activity, activity_coverage = regulon_activity(z, prior_panel)
        activity.insert(0, "dataset", "GSE135328")
        activity.insert(1, "cell_line", cell_line)
        activities_all.append(activity)
        activity_wide = activity.pivot(index="sample_id", columns="tf", values="activity").add_prefix("tf_activity__").reset_index()
        scores = scores.merge(activity_wide, on="sample_id", how="left")
        data.append(scores)
        coverage.insert(0, "dataset", f"GSE135328_{cell_line}")
        coverage_all.append(coverage)
        activity_coverage.insert(0, "dataset", f"GSE135328_{cell_line}")
        coverage_all.append(activity_coverage.rename(columns={"tf": "feature", "n_targets": "n_present"}).assign(n_expected=np.nan))
        score_columns = [column for column in scores.columns if column not in metadata.columns]
        clone_scores = scores.groupby(["cell_line", "clone_id", "clone_number", "genotype"], observed=True)[score_columns].mean().reset_index()
        clone_scores_all.append(clone_scores)
        wt = clone_scores[clone_scores["genotype"].eq("WT")].iloc[0]
        for row in clone_scores[clone_scores["genotype"].isin(["KO", "Het"])].itertuples(index=False):
            for feature in score_columns:
                clone_contrasts.append(
                    {
                        "cell_line": cell_line,
                        "clone_id": row.clone_id,
                        "genotype": row.genotype,
                        "reference_clone": wt["clone_id"],
                        "feature": feature,
                        "difference_vs_WT": getattr(row, feature) - wt[feature],
                    }
                )
        clone_expression = expression.T.join(metadata.set_index("sample_id")[["clone_id", "genotype"]]).groupby(["clone_id", "genotype"], observed=True).mean()
        wt_expression = clone_expression.xs("WT", level="genotype").iloc[0]
        ko_expression = clone_expression.xs("KO", level="genotype")
        for clone_id, row in ko_expression.iterrows():
            delta = row - wt_expression
            tmp = pd.DataFrame(
                {
                    "cell_line": cell_line,
                    "ko_clone": clone_id,
                    "reference_clone": clone_expression.xs("WT", level="genotype").index[0],
                    "gene": delta.index,
                    "ko_minus_wt_log2cpm": delta.to_numpy(),
                }
            )
            gene_effects_all.append(tmp)
        null, test = simple_matched_test(expression, z, signature, metadata, cell_line)
        null_all.append(null)
        tests_all.append(test)
    return {
        "sample_scores": pd.concat(data, ignore_index=True),
        "clone_scores": pd.concat(clone_scores_all, ignore_index=True),
        "clone_contrasts": pd.DataFrame(clone_contrasts),
        "gene_effects": pd.concat(gene_effects_all, ignore_index=True),
        "matched_null": pd.concat(null_all, ignore_index=True),
        "matched_tests": pd.concat(tests_all, ignore_index=True),
        "activity": pd.concat(activities_all, ignore_index=True),
        "coverage": pd.concat(coverage_all, ignore_index=True, sort=False),
        "expression_by_cell_line": expression_by_cell_line,
    }


def chen_regulon_analysis(prior_panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    activity_rows = []
    test_rows = []
    coverage_rows = []
    for dataset in ["discovery", "validation"]:
        metadata = pd.read_csv(CHEN_DIR / f"chen_{dataset}_specimen_pseudobulk_meta.tsv", sep="\t")
        expression = pd.read_csv(CHEN_DIR / f"chen_{dataset}_specimen_pseudobulk_expression.tsv.gz", sep="\t")
        expression.index = metadata["specimen_id"].astype(str)
        z = gene_zscores(expression.T)
        activity, coverage = regulon_activity(z, prior_panel)
        activity = activity.merge(
            metadata[["specimen_id", "donor_id", "route_group"]],
            left_on="sample_id",
            right_on="specimen_id",
            how="left",
            validate="many_to_one",
        )
        activity.insert(0, "dataset", dataset)
        activity_rows.append(activity.drop(columns="specimen_id"))
        coverage.insert(0, "dataset", dataset)
        coverage_rows.append(coverage)
        for tf, frame in activity.groupby("tf", observed=True):
            adenoma = frame.loc[frame["route_group"].eq("conventional_adenoma"), "activity"].dropna()
            normal = frame.loc[frame["route_group"].eq("normal"), "activity"].dropna()
            if len(adenoma) >= 3 and len(normal) >= 3:
                mw = mannwhitneyu(adenoma, normal, alternative="two-sided")
                test_rows.append(
                    {
                        "dataset": dataset,
                        "tf": tf,
                        "analysis": "specimen",
                        "n_adenoma": len(adenoma),
                        "n_normal": len(normal),
                        "median_adenoma_minus_normal": adenoma.median() - normal.median(),
                        "auc": mw.statistic / (len(adenoma) * len(normal)),
                        "p_value": mw.pvalue,
                    }
                )
            donor = frame.groupby(["donor_id", "route_group"], observed=True)["activity"].median().unstack()
            if {"conventional_adenoma", "normal"}.issubset(donor.columns):
                diff = (donor["conventional_adenoma"] - donor["normal"]).dropna()
                if len(diff) >= 3:
                    test_rows.append(
                        {
                            "dataset": dataset,
                            "tf": tf,
                            "analysis": "paired_donor",
                            "n_adenoma": len(diff),
                            "n_normal": len(diff),
                            "median_adenoma_minus_normal": diff.median(),
                            "auc": np.nan,
                            "p_value": wilcoxon(diff, zero_method="wilcox", alternative="two-sided").pvalue,
                        }
                    )
    return pd.concat(activity_rows, ignore_index=True), pd.DataFrame(test_rows), pd.concat(coverage_rows, ignore_index=True)


def virtual_knockout(
    signature: pd.DataFrame,
    route_prior: pd.DataFrame,
    panel_prior: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    route_weights = signature.set_index("gene")["route_weight"]
    edges = route_prior[route_prior["target_genesymbol"].isin(route_weights.index)].copy()
    edges = edges[edges["source_genesymbol"].astype(str).str.match(r"^[A-Z0-9-]+$")]
    edges["route_weight"] = edges["target_genesymbol"].map(route_weights)
    edges["virtual_ko_component"] = -edges["edge_sign"] * edges["route_weight"]
    rng = np.random.default_rng(RANDOM_SEED)
    target_list = route_weights.index.to_numpy()
    weight_values = route_weights.to_numpy()
    rows = []
    for tf, frame in edges.groupby("source_genesymbol", observed=True):
        frame = frame.drop_duplicates("target_genesymbol")
        observed = frame["virtual_ko_component"].mean()
        null = np.empty(N_PERMUTATIONS, dtype=float)
        target_positions = pd.Index(target_list).get_indexer(frame["target_genesymbol"])
        edge_sign = frame["edge_sign"].to_numpy()
        for i in range(N_PERMUTATIONS):
            shuffled = rng.permutation(weight_values)
            null[i] = np.mean(-edge_sign * shuffled[target_positions])
        rows.append(
            {
                "tf": tf,
                "n_direct_locked_targets": len(frame),
                "n_route_up_targets": int((frame["route_weight"] > 0).sum()),
                "n_route_down_targets": int((frame["route_weight"] < 0).sum()),
                "direct_virtual_ko_route_impact": observed,
                "fraction_components_predicting_route_collapse": float((frame["virtual_ko_component"] < 0).mean()),
                "null_mean": null.mean(),
                "null_sd": null.std(ddof=1),
                "p_route_collapse_lower_tail": (1 + np.sum(null <= observed)) / (N_PERMUTATIONS + 1),
            }
        )
    ranking = pd.DataFrame(rows)
    eligible = ranking[ranking["n_direct_locked_targets"] >= 3].copy()
    eligible["collapse_rank"] = eligible["direct_virtual_ko_route_impact"].rank(method="min", ascending=True)
    ranking = ranking.merge(eligible[["tf", "collapse_rank"]], on="tf", how="left")
    ranking = ranking.sort_values(["direct_virtual_ko_route_impact", "n_direct_locked_targets"], ascending=[True, False])

    # Candidate-specific two-hop sensitivity: TF -> intermediate TF -> locked target.
    route_by_tf = {tf: frame for tf, frame in edges.groupby("source_genesymbol", observed=True)}
    panel_rows = []
    path_rows = []
    for tf in TF_PANEL:
        direct = edges[edges["source_genesymbol"].eq(tf)].drop_duplicates("target_genesymbol")
        direct_sum = direct["virtual_ko_component"].sum()
        outgoing = panel_prior[panel_prior["source_genesymbol"].eq(tf)]
        two_hop_components = []
        for first in outgoing.itertuples(index=False):
            intermediate = first.target_genesymbol
            if intermediate not in route_by_tf:
                continue
            for second in route_by_tf[intermediate].drop_duplicates("target_genesymbol").itertuples(index=False):
                route_weight = route_weights[second.target_genesymbol]
                component = -first.edge_sign * second.edge_sign * route_weight * 0.5
                two_hop_components.append(component)
                path_rows.append(
                    {
                        "tf": tf,
                        "intermediate_tf": intermediate,
                        "locked_target": second.target_genesymbol,
                        "path_sign": first.edge_sign * second.edge_sign,
                        "route_weight": route_weight,
                        "attenuated_virtual_ko_component": component,
                    }
                )
        total_n = len(direct) + len(two_hop_components)
        panel_rows.append(
            {
                "tf": tf,
                "n_direct_paths": len(direct),
                "n_two_hop_paths": len(two_hop_components),
                "direct_component_sum": direct_sum,
                "two_hop_component_sum": sum(two_hop_components),
                "combined_virtual_ko_route_impact": (direct_sum + sum(two_hop_components)) / total_n if total_n else np.nan,
            }
        )
    return ranking, pd.DataFrame(panel_rows), pd.DataFrame(path_rows)


def tcf_target_empirical_check(
    signature: pd.DataFrame,
    panel_prior: pd.DataFrame,
    gse135_gene_effects: pd.DataFrame,
) -> pd.DataFrame:
    tcf = panel_prior[panel_prior["source_genesymbol"].eq("TCF7L2")].copy()
    tcf = tcf[tcf["target_genesymbol"].isin(signature["gene"])]
    tcf["predicted_ko_direction"] = -tcf["edge_sign"]
    effects = (
        gse135_gene_effects.groupby(["cell_line", "gene"], observed=True)["ko_minus_wt_log2cpm"]
        .mean()
        .reset_index()
    )
    out = tcf.merge(effects, left_on="target_genesymbol", right_on="gene", how="left")
    out["empirical_direction"] = np.sign(out["ko_minus_wt_log2cpm"])
    out["direction_concordant"] = out["predicted_ko_direction"].eq(out["empirical_direction"])
    return out


def write_outputs(
    gse125: dict[str, pd.DataFrame],
    gse135: dict[str, pd.DataFrame],
    chen_activity: pd.DataFrame,
    chen_tests: pd.DataFrame,
    chen_coverage: pd.DataFrame,
    ranking: pd.DataFrame,
    panel: pd.DataFrame,
    paths: pd.DataFrame,
    tcf_check: pd.DataFrame,
) -> None:
    tables = {
        "gse125472_sample_scores.tsv": gse125["sample_scores"],
        "gse125472_feature_coverage.tsv": gse125["coverage"],
        "gse125472_regulon_activity.tsv": gse125["activity"],
        "gse125472_regulon_coverage.tsv": gse125["activity_coverage"],
        "gse125472_donor_contrasts.tsv": gse125["donor_contrasts"],
        "gse125472_contrast_summary.tsv": gse125["contrast_summary"],
        "gse125472_leave_one_donor_out.tsv": gse125["leave_one_donor_out"],
        "gse125472_matched_signature_tests.tsv": gse125["matched_tests"],
        "gse125472_gene_donor_contrasts.tsv": gse125["gene_donor_contrasts"],
        "gse125472_gene_effects.tsv": gse125["gene_effects"],
        "gse135328_sample_scores.tsv": gse135["sample_scores"],
        "gse135328_clone_scores.tsv": gse135["clone_scores"],
        "gse135328_clone_contrasts.tsv": gse135["clone_contrasts"],
        "gse135328_gene_effects.tsv": gse135["gene_effects"],
        "gse135328_matched_signature_tests.tsv": gse135["matched_tests"],
        "gse135328_regulon_activity.tsv": gse135["activity"],
        "gse135328_feature_coverage.tsv": gse135["coverage"],
        "chen_regulon_activity.tsv": chen_activity,
        "chen_regulon_activity_tests.tsv": chen_tests,
        "chen_regulon_coverage.tsv": chen_coverage,
        "virtual_tf_knockout_ranking.tsv": ranking,
        "virtual_tf_knockout_panel.tsv": panel,
        "virtual_tf_knockout_two_hop_paths.tsv": paths,
        "tcf7l2_predicted_empirical_target_check.tsv": tcf_check,
    }
    for filename, frame in tables.items():
        frame.to_csv(OUT / filename, sep="\t", index=False)
    gse125["matched_null"].to_csv(OUT / "gse125472_matched_signature_null.tsv.gz", sep="\t", index=False, compression="gzip")
    gse135["matched_null"].to_csv(OUT / "gse135328_matched_signature_null.tsv.gz", sep="\t", index=False, compression="gzip")

    source_paths = [
        SIGNATURE_PATH,
        GSE125_PATH,
        GSE135_DIR / "GSE135328_count_HCT116.txt.gz",
        GSE135_DIR / "GSE135328_count_HT29.txt.gz",
        GENE_INFO_PATH,
        COLLECTRI_PANEL_PATH,
        COLLECTRI_ROUTE_PATH,
    ]
    source_manifest = pd.DataFrame(
        [
            {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in source_paths
        ]
    )
    source_manifest.to_csv(OUT / "source_file_manifest.tsv", sep="\t", index=False)

    primary = gse125["contrast_summary"].query(
        "feature == 'route_score' and comparison == 'APC_vs_WT_with_Wnt'"
    ).iloc[0]
    without = gse125["contrast_summary"].query(
        "feature == 'route_score' and comparison == 'APC_vs_WT_without_Wnt'"
    ).iloc[0]
    withdrawal = gse125["contrast_summary"].query(
        "feature == 'route_score' and comparison == 'WT_withdrawal'"
    ).iloc[0]
    matched_primary = gse125["matched_tests"].query("comparison == 'APC_vs_WT_with_Wnt'").iloc[0]
    tcf_route = gse135["clone_contrasts"].query("feature == 'route_score' and genotype == 'KO'")
    tcf_activity = gse135["clone_contrasts"].query("feature == 'tf_activity__TCF7L2' and genotype == 'KO'")
    tcf_route_text = tcf_route[["cell_line", "clone_id", "difference_vs_WT"]].to_string(index=False)
    tcf_activity_text = tcf_activity[["cell_line", "clone_id", "difference_vs_WT"]].to_string(index=False)
    summary = f"""# Perturbation validation summary

## Primary empirical perturbation

- GSE125472 contains three independent donor-matched isogenic human colon organoid systems.
- APC-KO minus WT route-score change with Wnt/R-spondin: mean {primary.mean_difference:.4f}, median {primary.median_difference:.4f}, range {primary.min_difference:.4f} to {primary.max_difference:.4f}; expected direction in {int(primary.n_expected_direction)}/{int(primary.n_units)} donors.
- Expression-matched random-signature test: P={matched_primary.p_expression_matched_one_sided:.6g} ({int(matched_primary.n_permutations):,} permutations).
- APC-KO minus WT without Wnt/R-spondin: mean {without.mean_difference:.4f}; expected direction in {int(without.n_expected_direction)}/{int(without.n_units)} donors.
- WT Wnt/R-spondin withdrawal: mean {withdrawal.mean_difference:.4f}; expected negative direction in {int(withdrawal.n_expected_direction)}/{int(withdrawal.n_units)} donors.

## TCF7L2 stress test

```text
{tcf_route_text}
```

TCF7L2 regulon-activity calibration:

```text
{tcf_activity_text}
```

## Evidential interpretation

The APC organoid analysis is the main perturbational validation because it uses the relevant human colonic epithelium, isogenic editing and donor matching. The TCF7L2 analysis is a context-dependent falsification layer and must not be pooled across cell lines as if clones were patient replicates. The signed virtual knockout is topology based and is supportive only.
"""
    (OUT / "summary.md").write_text(summary, encoding="utf-8")

    manifest = {
        "analysis": "locked_route_perturbation_validation",
        "date": "2026-08-08",
        "random_seed": RANDOM_SEED,
        "n_permutations": N_PERMUTATIONS,
        "n_bootstraps": N_BOOTSTRAPS,
        "signature_genes": 100,
        "primary_dataset": "GSE125472",
        "stress_test_dataset": "GSE135328",
        "virtual_prior": "CollecTRI",
        "claim_boundary": "computational and public perturbation support; not a replacement for local TMA-IHC/mIF",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    signature = read_signature()
    prior_panel = signed_collectri(COLLECTRI_PANEL_PATH)
    prior_route = signed_collectri(COLLECTRI_ROUTE_PATH)
    gse125 = analyze_gse125(signature, prior_panel)
    gse135 = analyze_gse135(signature, prior_panel)
    chen_activity, chen_tests, chen_coverage = chen_regulon_analysis(prior_panel)
    ranking, panel, paths = virtual_knockout(signature, prior_route, prior_panel)
    tcf_check = tcf_target_empirical_check(signature, prior_panel, gse135["gene_effects"])
    write_outputs(
        gse125,
        gse135,
        chen_activity,
        chen_tests,
        chen_coverage,
        ranking,
        panel,
        paths,
        tcf_check,
    )
    print(f"Wrote perturbation validation outputs to {OUT}")


if __name__ == "__main__":
    main()
