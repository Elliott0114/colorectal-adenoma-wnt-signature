#!/usr/bin/env python3
"""Candidate regulatory-window accessibility in Becker multiome scATAC data.

The Becker GEO archive provides tabix-indexed fragments, but not a peak matrix
or motif-deviation matrix. This script therefore adds a locus-resolved layer
without claiming chromVAR or TF binding: it counts fragments in fixed 1 kb
windows across +/-20 kb around candidate route loci, then tests whether WNT
route and WNT/TCF/ASCL2-axis windows are more accessible in polyps.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

from becker_multiome_tss_accessibility import (
    DATA_DIR,
    EXTRACT_DIR,
    GENCODE_GTF,
    RAW_TAR,
    SERIES_MATRIX,
    count_window,
    extract_selected,
    list_fragment_members,
    load_tss_windows,
    parse_series_matrix,
    select_samples,
)


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "becker_multiome_regulatory_windows"
TSS_OUT_DIR = ROOT / "results" / "becker_multiome_tss"
RNA_ATAC_DIR = ROOT / "results" / "becker_rna_atac_concordance"

WINDOW_SIZE_BP = 1_000
FLANK_BP = 20_000
WINDOW_SET_ID = f"tss_pm{FLANK_BP}_bin{WINDOW_SIZE_BP}_v1"

LOCUS_SETS = {
    "wnt_route_loci": ["LGR5", "ASCL2", "OLFM4", "AXIN2", "SOX9", "EPHB2", "SMOC2"],
    "wnt_tcf_ascl2_axis_loci": ["ASCL2", "TCF7L2", "TCF7", "LEF1", "MYC"],
    "immune_ifn_loci": ["HLA-DRA", "HLA-DRB1", "B2M", "TAP1", "STAT1", "IRF1", "CXCL10", "CXCL11"],
    "proliferation_control_loci": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "TYMS"],
    "serrated_metaplasia_loci": ["MUC5AC", "MUC6", "TFF1", "TFF2", "TFF3", "REG4", "AGR2", "SPINK4", "KRT7", "ANXA10"],
    "housekeeping_control_loci": ["ACTB", "GAPDH", "RPLP0", "RPL13A", "TUBB", "HPRT1", "PPIA", "RPS18", "TBP", "GUSB", "YWHAZ", "SDHA", "HMBS"],
}

PRIMARY_SETS = [
    "wnt_route_loci",
    "wnt_tcf_ascl2_axis_loci",
    "immune_ifn_loci",
    "proliferation_control_loci",
    "serrated_metaplasia_loci",
    "housekeeping_control_loci",
]


def distance_bin(distance: float) -> str:
    value = abs(float(distance))
    if value <= 500:
        return "tss_core_1kb"
    if value <= 2_500:
        return "promoter_proximal_2_5kb"
    if value <= 10_000:
        return "proximal_regulatory_10kb"
    return "distal_flank_20kb"


def locus_membership(gene: str) -> str:
    hits = [name for name, genes in LOCUS_SETS.items() if gene in genes]
    return ";".join(hits)


def build_windows() -> tuple[pd.DataFrame, pd.DataFrame]:
    wanted_genes = sorted({gene for genes in LOCUS_SETS.values() for gene in genes})
    tss = load_tss_windows(GENCODE_GTF, wanted_genes)
    rows = []
    for _, gene in tss.iterrows():
        tss0 = int(gene["tss_1based"]) - 1
        strand = str(gene["strand"])
        for offset in range(-FLANK_BP, FLANK_BP, WINDOW_SIZE_BP):
            start = max(0, tss0 + offset)
            end = max(start + 1, tss0 + offset + WINDOW_SIZE_BP)
            midpoint = (start + end) / 2
            genomic_distance = midpoint - tss0
            oriented_distance = genomic_distance if strand == "+" else -genomic_distance
            rows.append(
                {
                    "window_set_id": WINDOW_SET_ID,
                    "window_id": f"{gene['gene']}|{offset:+06d}|{WINDOW_SIZE_BP}",
                    "gene": gene["gene"],
                    "locus_set_membership": locus_membership(str(gene["gene"])),
                    "chrom": gene["chrom"],
                    "strand": strand,
                    "tss_1based": int(gene["tss_1based"]),
                    "window_start_0based": int(start),
                    "window_end_0based_exclusive": int(end),
                    "window_size_bp": int(end - start),
                    "genomic_midpoint_distance_to_tss": float(genomic_distance),
                    "oriented_midpoint_distance_to_tss": float(oriented_distance),
                    "distance_bin": distance_bin(genomic_distance),
                }
            )
    windows = pd.DataFrame(rows)
    present = set(tss["gene"])
    availability = []
    for locus_set, genes in LOCUS_SETS.items():
        hit = [gene for gene in genes if gene in present]
        miss = [gene for gene in genes if gene not in present]
        availability.append(
            {
                "locus_set": locus_set,
                "n_requested": len(genes),
                "n_present": len(hit),
                "present_genes": ",".join(hit),
                "missing_genes": ",".join(miss),
                "n_windows": int(windows["locus_set_membership"].str.contains(locus_set, regex=False).sum()),
            }
        )
    return windows, pd.DataFrame(availability)


def existing_counts_are_current(path: Path, selected: pd.DataFrame, windows: pd.DataFrame) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        preview = pd.read_csv(path, sep="\t", usecols=["window_set_id", "multiome_geo_accession", "window_id"])
    except Exception:
        return False
    if set(preview["window_set_id"].dropna().unique()) != {WINDOW_SET_ID}:
        return False
    return (
        set(preview["multiome_geo_accession"].unique()) == set(selected["multiome_geo_accession"])
        and set(preview["window_id"].unique()) == set(windows["window_id"])
    )


def count_sample(sample: pd.Series, windows: pd.DataFrame) -> list[dict[str, object]]:
    import pysam

    frag_path = EXTRACT_DIR / sample["fragment_file"]
    rows = []
    with pysam.TabixFile(str(frag_path)) as tabix:
        for _, window in windows.iterrows():
            count = count_window(
                tabix,
                str(window["chrom"]),
                int(window["window_start_0based"]),
                int(window["window_end_0based_exclusive"]),
            )
            rows.append(
                {
                    "window_set_id": WINDOW_SET_ID,
                    "multiome_geo_accession": sample["multiome_geo_accession"],
                    "sample_core": sample["sample_core"],
                    "window_id": window["window_id"],
                    "fragment_count": int(count),
                }
            )
    return rows


def load_or_count(selected: pd.DataFrame, windows: pd.DataFrame) -> pd.DataFrame:
    count_path = OUT_DIR / "becker_multiome_regulatory_window_counts.tsv"
    if existing_counts_are_current(count_path, selected, windows):
        return pd.read_csv(count_path, sep="\t", low_memory=False)

    extract_selected(selected)
    rows: list[dict[str, object]] = []
    for i, (_, sample) in enumerate(selected.iterrows(), start=1):
        print(f"[{i}/{len(selected)}] {sample['multiome_geo_accession']} {sample['sample_core']} {sample['disease_stage_group']}", flush=True)
        rows.extend(count_sample(sample, windows))
    counts = pd.DataFrame(rows)
    counts.to_csv(count_path, sep="\t", index=False)
    return counts


def add_normalized_values(counts: pd.DataFrame, selected: pd.DataFrame, windows: pd.DataFrame) -> pd.DataFrame:
    frame = counts.merge(windows, on=["window_set_id", "window_id"], how="left")
    frame = frame.merge(selected, on=["multiome_geo_accession", "sample_core"], how="left")
    totals = (
        frame.groupby(["multiome_geo_accession", "sample_core"], as_index=False)["fragment_count"]
        .sum()
        .rename(columns={"fragment_count": "regulatory_window_total_counts"})
    )
    frame = frame.merge(totals, on=["multiome_geo_accession", "sample_core"], how="left")
    frame["window_log_cpm"] = np.log1p(frame["fragment_count"] / frame["regulatory_window_total_counts"].clip(lower=1) * 1_000_000.0)
    frame["window_z"] = frame.groupby("window_id")["window_log_cpm"].transform(
        lambda x: (x - x.mean()) / x.std(ddof=1) if x.std(ddof=1) and not pd.isna(x.std(ddof=1)) else np.nan
    )
    frame["log10_regulatory_window_total_counts"] = np.log10(frame["regulatory_window_total_counts"].clip(lower=1))
    return frame


def safe_mannwhitney(a: pd.Series, b: pd.Series) -> float:
    a = a.dropna().astype(float)
    b = b.dropna().astype(float)
    if len(a) < 3 or len(b) < 3:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def safe_wilcoxon(a: pd.Series, b: pd.Series) -> float:
    data = pd.DataFrame({"a": a, "b": b}).dropna().astype(float)
    if len(data) < 3:
        return np.nan
    delta = data["a"] - data["b"]
    if np.isclose(delta.abs().sum(), 0):
        return np.nan
    return float(wilcoxon(data["a"], data["b"]).pvalue)


def window_tests(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    main = frame.loc[frame["disease_stage_group"].isin(["polyp", "normal_unaffected"])]
    for window_id, part in main.groupby("window_id", sort=False):
        a = part.loc[part["disease_stage_group"].eq("polyp"), "window_z"]
        b = part.loc[part["disease_stage_group"].eq("normal_unaffected"), "window_z"]
        first = part.iloc[0]
        rows.append(
            {
                "comparison": "polyp_vs_normal_unaffected",
                "window_id": window_id,
                "gene": first["gene"],
                "locus_set_membership": first["locus_set_membership"],
                "chrom": first["chrom"],
                "window_start_0based": int(first["window_start_0based"]),
                "window_end_0based_exclusive": int(first["window_end_0based_exclusive"]),
                "distance_bin": first["distance_bin"],
                "oriented_midpoint_distance_to_tss": float(first["oriented_midpoint_distance_to_tss"]),
                "n_polyp": int(a.notna().sum()),
                "n_normal_unaffected": int(b.notna().sum()),
                "median_polyp": float(a.median()) if a.notna().any() else np.nan,
                "median_normal_unaffected": float(b.median()) if b.notna().any() else np.nan,
                "delta_median_polyp_minus_normal": float(a.median() - b.median()) if a.notna().any() and b.notna().any() else np.nan,
                "p_value": safe_mannwhitney(a, b),
            }
        )
    out = pd.DataFrame(rows)
    ok = out["p_value"].notna()
    out["q_value_bh_all_windows"] = np.nan
    if ok.any():
        out.loc[ok, "q_value_bh_all_windows"] = multipletests(out.loc[ok, "p_value"], method="fdr_bh")[1]
    return out.sort_values(["q_value_bh_all_windows", "p_value", "gene"], na_position="last")


def distance_bin_scores(frame: pd.DataFrame) -> pd.DataFrame:
    memberships = []
    for locus_set in PRIMARY_SETS:
        part = frame.loc[frame["locus_set_membership"].str.contains(locus_set, regex=False, na=False)]
        grouped = (
            part.groupby(
                [
                    "multiome_geo_accession",
                    "sample_core",
                    "disease_stage_group",
                    "patient_id",
                    "familial_adenomatous_polyposis",
                    "sex",
                    "distance_bin",
                ],
                as_index=False,
            )
            .agg(
                regulatory_window_score=("window_z", "mean"),
                mean_window_log_cpm=("window_log_cpm", "mean"),
                n_windows=("window_id", "nunique"),
                regulatory_window_total_counts=("regulatory_window_total_counts", "first"),
                log10_regulatory_window_total_counts=("log10_regulatory_window_total_counts", "first"),
            )
            .assign(locus_set=locus_set)
        )
        memberships.append(grouped)
    return pd.concat(memberships, ignore_index=True)


def distance_bin_tests(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    main = scores.loc[scores["disease_stage_group"].isin(["polyp", "normal_unaffected"])]
    for (locus_set, bin_name), part in main.groupby(["locus_set", "distance_bin"], sort=False):
        a = part.loc[part["disease_stage_group"].eq("polyp"), "regulatory_window_score"]
        b = part.loc[part["disease_stage_group"].eq("normal_unaffected"), "regulatory_window_score"]
        rows.append(
            {
                "comparison": "polyp_vs_normal_unaffected",
                "locus_set": locus_set,
                "distance_bin": bin_name,
                "n_polyp": int(a.notna().sum()),
                "n_normal_unaffected": int(b.notna().sum()),
                "median_polyp": float(a.median()) if a.notna().any() else np.nan,
                "median_normal_unaffected": float(b.median()) if b.notna().any() else np.nan,
                "delta_median_polyp_minus_normal": float(a.median() - b.median()) if a.notna().any() and b.notna().any() else np.nan,
                "p_value": safe_mannwhitney(a, b),
            }
        )
    out = pd.DataFrame(rows)
    out["q_value_bh_within_distance_bin"] = np.nan
    for bin_name, idx in out.groupby("distance_bin").groups.items():
        p = out.loc[idx, "p_value"]
        ok = p.notna()
        if ok.any():
            out.loc[p.index[ok], "q_value_bh_within_distance_bin"] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def fap_only_distance_bin_tests(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    main = scores.loc[
        scores["disease_stage_group"].isin(["polyp", "normal_unaffected"])
        & scores["familial_adenomatous_polyposis"].eq("Y")
    ]
    for (locus_set, bin_name), part in main.groupby(["locus_set", "distance_bin"], sort=False):
        a = part.loc[part["disease_stage_group"].eq("polyp"), "regulatory_window_score"]
        b = part.loc[part["disease_stage_group"].eq("normal_unaffected"), "regulatory_window_score"]
        rows.append(
            {
                "comparison": "fap_only_polyp_vs_normal_unaffected",
                "locus_set": locus_set,
                "distance_bin": bin_name,
                "n_polyp": int(a.notna().sum()),
                "n_normal_unaffected": int(b.notna().sum()),
                "median_polyp": float(a.median()) if a.notna().any() else np.nan,
                "median_normal_unaffected": float(b.median()) if b.notna().any() else np.nan,
                "delta_median_polyp_minus_normal": float(a.median() - b.median()) if a.notna().any() and b.notna().any() else np.nan,
                "p_value": safe_mannwhitney(a, b),
            }
        )
    out = pd.DataFrame(rows)
    out["q_value_bh_within_distance_bin"] = np.nan
    for bin_name, idx in out.groupby("distance_bin").groups.items():
        p = out.loc[idx, "p_value"]
        ok = p.notna()
        if ok.any():
            out.loc[p.index[ok], "q_value_bh_within_distance_bin"] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def paired_patient_distance_bin_tests(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    main = scores.loc[scores["disease_stage_group"].isin(["polyp", "normal_unaffected"])].copy()
    patient_scores = (
        main.groupby(["patient_id", "disease_stage_group", "locus_set", "distance_bin"], as_index=False)["regulatory_window_score"]
        .mean()
    )
    for (locus_set, bin_name), part in patient_scores.groupby(["locus_set", "distance_bin"], sort=False):
        wide = part.pivot_table(index="patient_id", columns="disease_stage_group", values="regulatory_window_score")
        wide = wide.dropna(subset=["polyp", "normal_unaffected"])
        delta = wide["polyp"] - wide["normal_unaffected"] if len(wide) else pd.Series(dtype=float)
        rows.append(
            {
                "comparison": "patient_paired_polyp_vs_normal_unaffected",
                "locus_set": locus_set,
                "distance_bin": bin_name,
                "n_pairs": int(len(wide)),
                "median_delta_polyp_minus_normal": float(delta.median()) if len(delta) else np.nan,
                "mean_delta_polyp_minus_normal": float(delta.mean()) if len(delta) else np.nan,
                "p_value": safe_wilcoxon(wide["polyp"], wide["normal_unaffected"]) if len(wide) else np.nan,
                "paired_patient_ids": ",".join(map(str, wide.index.tolist())),
            }
        )
    out = pd.DataFrame(rows)
    out["q_value_bh_within_distance_bin"] = np.nan
    for bin_name, idx in out.groupby("distance_bin").groups.items():
        p = out.loc[idx, "p_value"]
        ok = p.notna()
        if ok.any():
            out.loc[p.index[ok], "q_value_bh_within_distance_bin"] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def adjusted_models(scores: pd.DataFrame) -> pd.DataFrame:
    covar_path = TSS_OUT_DIR / "becker_multiome_tss_module_scores.tsv"
    if not covar_path.exists():
        return pd.DataFrame()
    covars = pd.read_csv(
        covar_path,
        sep="\t",
        usecols=[
            "multiome_geo_accession",
            "sample_core",
            "score__proliferation_control",
            "score__housekeeping_control",
            "log10_tss_total_counts",
        ],
    )
    frame = scores.merge(covars, on=["multiome_geo_accession", "sample_core"], how="left")
    frame = frame.loc[frame["disease_stage_group"].isin(["polyp", "normal_unaffected"])].copy()
    frame["is_polyp"] = frame["disease_stage_group"].eq("polyp").astype(int)
    rows = []
    for locus_set in ["wnt_route_loci", "wnt_tcf_ascl2_axis_loci", "immune_ifn_loci", "serrated_metaplasia_loci"]:
        for bin_name in ["tss_core_1kb", "promoter_proximal_2_5kb", "proximal_regulatory_10kb", "distal_flank_20kb"]:
            data = frame.loc[(frame["locus_set"].eq(locus_set)) & (frame["distance_bin"].eq(bin_name))].copy()
            data = data[
                [
                    "regulatory_window_score",
                    "is_polyp",
                    "score__proliferation_control",
                    "score__housekeeping_control",
                    "log10_tss_total_counts",
                ]
            ].dropna()
            if len(data) < 10 or data["is_polyp"].nunique() < 2:
                continue
            formula = "regulatory_window_score ~ is_polyp + score__proliferation_control + score__housekeeping_control + log10_tss_total_counts"
            fit = smf.ols(formula, data=data).fit(cov_type="HC3")
            rows.append(
                {
                    "locus_set": locus_set,
                    "distance_bin": bin_name,
                    "term": "is_polyp",
                    "n": int(fit.nobs),
                    "coef": float(fit.params.get("is_polyp", np.nan)),
                    "se_hc3": float(fit.bse.get("is_polyp", np.nan)),
                    "p_value": float(fit.pvalues.get("is_polyp", np.nan)),
                    "r_squared": float(fit.rsquared),
                    "model": formula,
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value_bh"] = multipletests(out["p_value"], method="fdr_bh")[1]
    return out


def rna_atac_correlations(scores: pd.DataFrame) -> pd.DataFrame:
    paired_path = RNA_ATAC_DIR / "becker_rna_atac_paired_scores.tsv"
    if not paired_path.exists():
        return pd.DataFrame()
    paired = pd.read_csv(
        paired_path,
        sep="\t",
        usecols=[
            "multiome_geo_accession",
            "sample_core",
            "disease_stage_group",
            "rna_epi__wnt_stem",
            "rna_epi__ca_route_signature",
            "rna_epi__proliferation_control",
        ],
    )
    frame = scores.merge(paired, on=["multiome_geo_accession", "sample_core", "disease_stage_group"], how="inner")
    frame = frame.loc[frame["disease_stage_group"].isin(["polyp", "normal_unaffected"])].copy()
    rows = []
    for locus_set in ["wnt_route_loci", "wnt_tcf_ascl2_axis_loci"]:
        for bin_name in ["tss_core_1kb", "promoter_proximal_2_5kb", "proximal_regulatory_10kb", "distal_flank_20kb"]:
            part = frame.loc[(frame["locus_set"].eq(locus_set)) & (frame["distance_bin"].eq(bin_name))]
            for rna_feature in ["rna_epi__wnt_stem", "rna_epi__ca_route_signature", "rna_epi__proliferation_control"]:
                valid = part[["regulatory_window_score", rna_feature]].dropna()
                if len(valid) < 10:
                    continue
                stat = spearmanr(valid["regulatory_window_score"], valid[rna_feature])
                rows.append(
                    {
                        "subset": "normal_polyp",
                        "locus_set": locus_set,
                        "distance_bin": bin_name,
                        "rna_feature": rna_feature,
                        "n": len(valid),
                        "spearman_rho": float(stat.statistic),
                        "p_spearman": float(stat.pvalue),
                    }
                )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_spearman"] = multipletests(out["p_spearman"], method="fdr_bh")[1]
    return out


def write_summary(
    selected: pd.DataFrame,
    windows: pd.DataFrame,
    availability: pd.DataFrame,
    tests: pd.DataFrame,
    bin_tests: pd.DataFrame,
    fap_tests: pd.DataFrame,
    paired_tests: pd.DataFrame,
    models: pd.DataFrame,
    correlations: pd.DataFrame,
) -> None:
    lines = [
        "Becker multiome candidate regulatory-window accessibility",
        "=" * 66,
        "",
        f"Selected samples: {len(selected)}",
        selected["disease_stage_group"].value_counts().to_string(),
        "",
        f"Window set: {WINDOW_SET_ID}",
        f"Windows counted: {len(windows)} across {windows['gene'].nunique()} genes.",
        "Interpretation boundary: fixed regulatory windows around candidate loci; not peak calling, chromVAR, motif deviation, enhancer-promoter looping, or TF binding proof.",
        "",
        "Locus availability:",
        availability.to_string(index=False),
        "",
        "Primary distance-bin tests:",
    ]
    primary = bin_tests.loc[
        bin_tests["locus_set"].isin(["wnt_route_loci", "wnt_tcf_ascl2_axis_loci", "immune_ifn_loci"])
    ].sort_values(["distance_bin", "q_value_bh_within_distance_bin", "locus_set"])
    lines.append(
        primary[
            [
                "locus_set",
                "distance_bin",
                "n_polyp",
                "n_normal_unaffected",
                "delta_median_polyp_minus_normal",
                "p_value",
                "q_value_bh_within_distance_bin",
            ]
        ].to_string(index=False)
    )
    lines += ["", "FAP-only distance-bin sensitivity:"]
    fap_primary = fap_tests.loc[
        fap_tests["locus_set"].isin(["wnt_route_loci", "wnt_tcf_ascl2_axis_loci", "immune_ifn_loci"])
    ].sort_values(["distance_bin", "q_value_bh_within_distance_bin", "locus_set"])
    lines.append(
        fap_primary[
            [
                "locus_set",
                "distance_bin",
                "n_polyp",
                "n_normal_unaffected",
                "delta_median_polyp_minus_normal",
                "p_value",
                "q_value_bh_within_distance_bin",
            ]
        ].to_string(index=False)
    )
    lines += ["", "Patient-paired distance-bin sensitivity:"]
    paired_primary = paired_tests.loc[
        paired_tests["locus_set"].isin(["wnt_route_loci", "wnt_tcf_ascl2_axis_loci", "immune_ifn_loci"])
    ].sort_values(["distance_bin", "locus_set"])
    lines.append(
        paired_primary[
            [
                "locus_set",
                "distance_bin",
                "n_pairs",
                "median_delta_polyp_minus_normal",
                "p_value",
                "q_value_bh_within_distance_bin",
            ]
        ].to_string(index=False)
    )
    lines += ["", "Top positive WNT/axis windows:"]
    top = tests.loc[
        tests["locus_set_membership"].str.contains("wnt_route_loci|wnt_tcf_ascl2_axis_loci", regex=True, na=False)
        & tests["delta_median_polyp_minus_normal"].gt(0)
    ].head(20)
    if not top.empty:
        lines.append(top[["gene", "distance_bin", "oriented_midpoint_distance_to_tss", "delta_median_polyp_minus_normal", "p_value", "q_value_bh_all_windows"]].to_string(index=False))
    if not models.empty:
        lines += ["", "Adjusted models:"]
        lines.append(models[["locus_set", "distance_bin", "n", "coef", "se_hc3", "p_value", "q_value_bh", "r_squared"]].to_string(index=False))
    if not correlations.empty:
        lines += ["", "RNA-regulatory-window Spearman correlations:"]
        lines.append(correlations.to_string(index=False))
    (OUT_DIR / "becker_multiome_regulatory_window_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sample_meta = parse_series_matrix(SERIES_MATRIX)
    fragment_meta = list_fragment_members(RAW_TAR)
    selected = select_samples(fragment_meta, sample_meta)
    selected.to_csv(OUT_DIR / "becker_multiome_regulatory_window_sample_metadata.tsv", sep="\t", index=False)

    windows, availability = build_windows()
    windows.to_csv(OUT_DIR / "becker_multiome_regulatory_windows.tsv", sep="\t", index=False)
    availability.to_csv(OUT_DIR / "becker_multiome_regulatory_window_locus_availability.tsv", sep="\t", index=False)

    counts = load_or_count(selected, windows)
    normalized = add_normalized_values(counts, selected, windows)
    tests = window_tests(normalized)
    scores = distance_bin_scores(normalized)
    bin_tests = distance_bin_tests(scores)
    fap_tests = fap_only_distance_bin_tests(scores)
    paired_tests = paired_patient_distance_bin_tests(scores)
    models = adjusted_models(scores)
    correlations = rna_atac_correlations(scores)

    normalized.to_csv(OUT_DIR / "becker_multiome_regulatory_window_normalized_counts.tsv", sep="\t", index=False)
    tests.to_csv(OUT_DIR / "becker_multiome_regulatory_window_tests.tsv", sep="\t", index=False)
    tests.head(100).to_csv(OUT_DIR / "becker_multiome_regulatory_window_top100.tsv", sep="\t", index=False)
    scores.to_csv(OUT_DIR / "becker_multiome_regulatory_window_distance_bin_scores.tsv", sep="\t", index=False)
    bin_tests.to_csv(OUT_DIR / "becker_multiome_regulatory_window_distance_bin_tests.tsv", sep="\t", index=False)
    fap_tests.to_csv(OUT_DIR / "becker_multiome_regulatory_window_fap_only_distance_bin_tests.tsv", sep="\t", index=False)
    paired_tests.to_csv(OUT_DIR / "becker_multiome_regulatory_window_paired_patient_tests.tsv", sep="\t", index=False)
    models.to_csv(OUT_DIR / "becker_multiome_regulatory_window_adjusted_models.tsv", sep="\t", index=False)
    correlations.to_csv(OUT_DIR / "becker_multiome_regulatory_window_rna_correlations.tsv", sep="\t", index=False)
    write_summary(selected, windows, availability, tests, bin_tests, fap_tests, paired_tests, models, correlations)
    print(f"Wrote results to {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
