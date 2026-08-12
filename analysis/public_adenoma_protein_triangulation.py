#!/usr/bin/env python3
"""Triangulate locked transcript candidates in public colorectal adenoma proteomes.

This analysis does not re-select or modify the locked 50-up/50-down route. It
tests protein-level direction in PXD002137, FFPE detectability in PXD017269,
and tissue-region detectability in PXD046999. These evidence roles are kept
separate because only PXD002137 contains a normal-versus-adenoma comparison.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "data_sources" / "public_adenoma_proteomics"
LOCKED_DIR = ROOT / "results" / "route_signature_locked"
OUT_DIR = ROOT / "results" / "public_adenoma_protein_triangulation"

PANEL_ROLES = {
    "OLFM4": "current_primary_up_arm",
    "FABP1": "current_primary_down_arm",
    "SOX9": "current_context_nuclear_marker",
    "CTNNB1": "current_localization_readout",
    "KRT20": "prespecified_down_arm_alternative",
    "MKI67": "proliferation_covariate",
}

SOURCE_URLS = {
    "PXD002137/E-PROT-104-MappedToGeneID.txt": (
        "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/"
        "E-PROT-104/E-PROT-104-MappedToGeneID.txt"
    ),
    "PXD002137/E-PROT-104.condensed-sdrf.tsv": (
        "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/"
        "E-PROT-104/E-PROT-104.condensed-sdrf.tsv"
    ),
    "PXD002137/E-PROT-104-analysis-methods.tsv": (
        "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/"
        "E-PROT-104/E-PROT-104-analysis-methods.tsv"
    ),
    "PXD017269/Supplement_AdenomaProteomes.xlsx": (
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/04/PXD017269/"
        "Supplement_AdenomaProteomes.xlsx"
    ),
    "PXD017269/Rawfile_Information_AdenomaPatients.xlsx": (
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/04/PXD017269/"
        "Rawfile_Information_AdenomaPatients.xlsx"
    ),
    "PXD046999/PMC11381895_supplementaryFiles.zip": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11381895/"
        "supplementaryFiles?includeInlineImage=false"
    ),
    "PXD046999/mmc2.xlsx": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11381895/"
        "supplementaryFiles?includeInlineImage=false"
    ),
    "PXD046999/CAA__MBR_LibraryFree_DIANN181.pg_matrix.tsv": (
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/08/PXD046999/"
        "CAA__MBR_LibraryFree_DIANN181.pg_matrix.tsv"
    ),
    "PXD046999/CRA_revision_timstofultra_specificcohort_DIANN181_libraryfree.pg_matrix.tsv": (
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/08/PXD046999/"
        "CRA_revision_timstofultra_specificcohort_DIANN181_libraryfree.pg_matrix.tsv"
    ),
    "PXD046999/Metadata_DVPonCRA_Kabatnik_et_al.csv": (
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/08/PXD046999/"
        "Metadata_DVPonCRA_Kabatnik_et_al.csv"
    ),
    "PXD046999/Metadata_DVPonCRA_Kabatnik_et_al_revision_updated.csv": (
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/08/PXD046999/"
        "Metadata_DVPonCRA_Kabatnik_et_al_revision_updated.csv"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_genes(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [item.strip() for item in re.split(r"[;,]", str(value)) if item.strip()]


def bh(values: pd.Series) -> pd.Series:
    output = pd.Series(np.nan, index=values.index, dtype=float)
    keep = values.notna()
    if keep.any():
        output.loc[keep] = multipletests(values.loc[keep], method="fdr_bh")[1]
    return output


def candidate_inventory() -> pd.DataFrame:
    signature = pd.read_csv(
        LOCKED_DIR / "discovery_locked_signature_genes.tsv", sep="\t"
    )
    stability = pd.read_csv(
        LOCKED_DIR / "discovery_gene_stability_audit.tsv", sep="\t"
    )
    wanted = set(signature["gene"].astype(str)) | set(PANEL_ROLES)
    candidates = stability.loc[stability["gene"].isin(wanted)].copy()
    candidates["locked_signature_member"] = candidates["gene"].isin(
        set(signature["gene"].astype(str))
    )
    candidates["panel_role"] = candidates["gene"].map(PANEL_ROLES).fillna(
        "locked_signature_candidate"
    )
    candidates["expected_direction"] = np.sign(
        candidates["discovery_effect_adenoma_minus_normal"]
    ).astype(int)
    return candidates[
        [
            "gene",
            "panel_role",
            "locked_signature_member",
            "selected",
            "discovery_effect_adenoma_minus_normal",
            "bootstrap_ci_low",
            "bootstrap_ci_high",
            "direction_stability",
            "expected_direction",
        ]
    ].sort_values(
        ["locked_signature_member", "discovery_effect_adenoma_minus_normal"],
        ascending=[False, False],
    )


def pxd002137_metadata() -> pd.DataFrame:
    path = SOURCE_DIR / "PXD002137" / "E-PROT-104.condensed-sdrf.tsv"
    long = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["experiment", "blank", "sample_id", "kind", "field", "value", "ontology"],
    )
    meta = (
        long.pivot_table(
            index="sample_id", columns="field", values="value", aggfunc="first"
        )
        .reset_index()
        .rename_axis(columns=None)
    )
    meta["tissue_group"] = meta["sampling site"].map(
        {
            "adenoma": "adenoma",
            "normal tissue adjacent to tumor": "normal",
            "tumor": "crc",
        }
    )
    meta["age_years"] = pd.to_numeric(
        meta["age"].str.extract(r"(\d+)", expand=False), errors="coerce"
    )
    meta["is_male"] = meta["sex"].eq("male").astype(int)
    return meta


def adjusted_pxd002137_effect(values: pd.Series, meta: pd.DataFrame) -> dict[str, float]:
    data = meta.loc[meta["tissue_group"].isin(["normal", "adenoma"])].copy()
    data["abundance"] = values.reindex(data["sample_id"]).to_numpy(dtype=float)
    data = data.dropna(subset=["abundance", "age_years", "is_male"])
    group_counts = data["tissue_group"].value_counts()
    if len(data) < 8 or any(group_counts.get(group, 0) < 3 for group in ["normal", "adenoma"]):
        return {
            "age_sex_adjusted_log2_effect": np.nan,
            "age_sex_adjusted_ci_low": np.nan,
            "age_sex_adjusted_ci_high": np.nan,
            "age_sex_adjusted_p": np.nan,
        }
    data["log2_abundance"] = np.log2(data["abundance"].clip(lower=0) + 1)
    data["is_adenoma"] = data["tissue_group"].eq("adenoma").astype(int)
    design = sm.add_constant(data[["is_adenoma", "age_years", "is_male"]])
    fit = sm.OLS(data["log2_abundance"], design).fit(cov_type="HC3")
    interval = fit.conf_int().loc["is_adenoma"]
    return {
        "age_sex_adjusted_log2_effect": float(fit.params["is_adenoma"]),
        "age_sex_adjusted_ci_low": float(interval.iloc[0]),
        "age_sex_adjusted_ci_high": float(interval.iloc[1]),
        "age_sex_adjusted_p": float(fit.pvalues["is_adenoma"]),
    }


def pxd002137_tests(candidates: pd.DataFrame) -> pd.DataFrame:
    path = SOURCE_DIR / "PXD002137" / "E-PROT-104-MappedToGeneID.txt"
    matrix = pd.read_csv(path, sep="\t", na_values=["NA"])
    sample_cols = [column for column in matrix.columns if column.startswith("ppb.iBAQ.")]
    matrix["Gene.Symbol"] = matrix["Gene.Symbol"].astype(str)
    expression = matrix.groupby("Gene.Symbol", sort=False)[sample_cols].sum(min_count=1)
    meta = pxd002137_metadata()
    group_samples = {
        group: meta.loc[meta["tissue_group"].eq(group), "sample_id"].tolist()
        for group in ["normal", "adenoma", "crc"]
    }

    rows = []
    for gene in candidates["gene"]:
        present = gene in expression.index
        values = (
            expression.loc[gene].astype(float)
            if present
            else pd.Series(np.nan, index=sample_cols, dtype=float)
        )
        normal_raw = values.reindex(group_samples["normal"])
        adenoma_raw = values.reindex(group_samples["adenoma"])
        normal = np.log2(normal_raw.clip(lower=0) + 1)
        adenoma = np.log2(adenoma_raw.clip(lower=0) + 1)
        valid_normal = normal.dropna()
        valid_adenoma = adenoma.dropna()
        if valid_normal.empty or valid_adenoma.empty:
            u_value = p_value = auc = np.nan
        else:
            test = stats.mannwhitneyu(valid_adenoma, valid_normal, alternative="two-sided")
            u_value = float(test.statistic)
            p_value = float(test.pvalue)
            auc = u_value / (len(valid_adenoma) * len(valid_normal))
        pooled_sd = np.sqrt(
            (
                (len(valid_adenoma) - 1) * valid_adenoma.var(ddof=1)
                + (len(valid_normal) - 1) * valid_normal.var(ddof=1)
            )
            / max(len(valid_adenoma) + len(valid_normal) - 2, 1)
        )
        standardized = (
            (valid_adenoma.mean() - valid_normal.mean()) / pooled_sd
            if np.isfinite(pooled_sd) and pooled_sd > 0
            else np.nan
        )
        row = {
            "gene": gene,
            "protein_gene_present": present,
            "n_normal": len(valid_normal),
            "n_adenoma": len(valid_adenoma),
            "normal_detection_fraction": float((normal_raw.fillna(0) > 0).mean()),
            "adenoma_detection_fraction": float((adenoma_raw.fillna(0) > 0).mean()),
            "median_log2_normal": float(valid_normal.median()) if not valid_normal.empty else np.nan,
            "median_log2_adenoma": float(valid_adenoma.median()) if not valid_adenoma.empty else np.nan,
            "adenoma_minus_normal_log2_median": (
                float(valid_adenoma.median() - valid_normal.median())
                if not valid_normal.empty and not valid_adenoma.empty
                else np.nan
            ),
            "standardized_mean_difference": float(standardized),
            "auc_adenoma_greater_than_normal": float(auc),
            "mann_whitney_u": u_value,
            "p_value": p_value,
        }
        if present:
            row.update(adjusted_pxd002137_effect(values, meta))
        else:
            row.update(
                {
                    "age_sex_adjusted_log2_effect": np.nan,
                    "age_sex_adjusted_ci_low": np.nan,
                    "age_sex_adjusted_ci_high": np.nan,
                    "age_sex_adjusted_p": np.nan,
                }
            )
        rows.append(row)
    output = pd.DataFrame(rows)
    output["q_value_bh_candidates"] = bh(output["p_value"])
    output["age_sex_adjusted_q_bh_candidates"] = bh(output["age_sex_adjusted_p"])
    return output


def pxd017269_detectability(candidates: pd.DataFrame) -> pd.DataFrame:
    path = SOURCE_DIR / "PXD017269" / "Supplement_AdenomaProteomes.xlsx"
    matrix = pd.read_excel(path, sheet_name="Fig4A 118 adenomas")
    annotation = pd.read_excel(path, sheet_name="Sample annotation")
    sample_cols = list(matrix.columns[2:])
    annotation["sample_key"] = annotation["Raw file identifier"].astype(str)
    archival = annotation.set_index("sample_key")["Archival group"].to_dict()
    rows = []
    for gene in candidates["gene"]:
        matches = matrix["Genes"].map(lambda value: gene in split_genes(value))
        if not matches.any():
            rows.append(
                {
                    "gene": gene,
                    "protein_group_present": False,
                    "protein_accessions": "",
                    "n_samples": len(sample_cols),
                    "n_detected": 0,
                    "detection_fraction": 0.0,
                    "median_log2_intensity_detected": np.nan,
                    "archival_group_kruskal_p": np.nan,
                }
            )
            continue
        options = matrix.loc[matches].copy()
        options["detection_fraction"] = options[sample_cols].notna().mean(axis=1)
        options["median_intensity"] = options[sample_cols].median(axis=1, skipna=True)
        selected = options.sort_values(
            ["detection_fraction", "median_intensity"], ascending=False
        ).iloc[0]
        values = pd.to_numeric(selected[sample_cols], errors="coerce")
        groups = {}
        for column, value in values.items():
            key = str(int(column)) if isinstance(column, (int, float, np.integer, np.floating)) else str(column)
            group = archival.get(key)
            if group is not None and pd.notna(value):
                groups.setdefault(group, []).append(float(value))
        usable = [group for group in groups.values() if len(group) >= 2]
        archival_p = (
            float(stats.kruskal(*usable).pvalue) if len(usable) >= 2 else np.nan
        )
        rows.append(
            {
                "gene": gene,
                "protein_group_present": True,
                "protein_accessions": selected["ProteinAccessions"],
                "n_samples": len(sample_cols),
                "n_detected": int(values.notna().sum()),
                "detection_fraction": float(values.notna().mean()),
                "median_log2_intensity_detected": float(values.median(skipna=True)),
                "archival_group_kruskal_p": archival_p,
            }
        )
    output = pd.DataFrame(rows)
    output["archival_group_kruskal_q_bh_candidates"] = bh(
        output["archival_group_kruskal_p"]
    )
    return output


def gene_set(values: pd.Series) -> set[str]:
    genes: set[str] = set()
    for value in values.dropna():
        genes.update(split_genes(value))
    return genes


def pxd046999_presence(candidates: pd.DataFrame) -> pd.DataFrame:
    path = SOURCE_DIR / "PXD046999" / "mmc2.xlsx"
    comparisons = {
        "9CRAstudy_Volcano_CvsNMN": ("C", "NMN"),
        "9CRAstudy_Volcano_CvsHDA": ("C", "HDA"),
        "9CRAstudy_Volcano_HDAvsNMN": ("HDA", "NMN"),
    }
    detected: set[str] = set()
    enriched = {"C": set(), "HDA": set(), "NMN": set()}
    for sheet, (left_group, right_group) in comparisons.items():
        table = pd.read_excel(path, sheet_name=sheet, header=None)
        detected.update(gene_set(table.iloc[2:, 0]))
        enriched[left_group].update(gene_set(table.iloc[2:, 3]))
        enriched[right_group].update(gene_set(table.iloc[2:, 5]))
    return pd.DataFrame(
        {
            "gene": candidates["gene"],
            "detected_in_nine_patient_dvp": candidates["gene"].isin(detected),
            "enriched_in_metachronous_crc_group": candidates["gene"].isin(enriched["C"]),
            "enriched_in_metachronous_hda_group": candidates["gene"].isin(enriched["HDA"]),
            "enriched_in_non_metachronous_group": candidates["gene"].isin(enriched["NMN"]),
        }
    )


def evidence_tier(row: pd.Series) -> str:
    if not bool(row["locked_signature_member"]):
        return "context_only_not_locked"
    if not bool(row["protein_gene_present"]):
        return "C_transcript_locked_not_detected_in_pxd002137"
    protein_direction = np.sign(row["age_sex_adjusted_log2_effect"])
    direction_match = protein_direction == row["expected_direction"]
    well_detected = min(
        row["normal_detection_fraction"], row["adenoma_detection_fraction"]
    ) >= 0.5
    ffpe_detected = row["pxd017269_detection_fraction"] >= 0.8
    adjusted_q = row["age_sex_adjusted_q_bh_candidates"]
    if not direction_match and well_detected:
        return "D_direction_conflict"
    if direction_match and well_detected and ffpe_detected and adjusted_q <= 0.1:
        return "A_directional_protein_replication_and_ffpe_detectability"
    if direction_match and well_detected and ffpe_detected:
        return "B_directional_support_and_ffpe_detectability"
    return "C_incomplete_public_protein_support"


def source_manifest() -> pd.DataFrame:
    rows = []
    for path in sorted(SOURCE_DIR.rglob("*")):
        if not path.is_file():
            continue
        relative = str(path.relative_to(SOURCE_DIR))
        rows.append(
            {
                "file": str(path.relative_to(ROOT)),
                "source_url": SOURCE_URLS.get(relative, ""),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    candidates = candidate_inventory()
    pxd2 = pxd002137_tests(candidates)
    pxd17 = pxd017269_detectability(candidates).rename(
        columns={
            "protein_group_present": "pxd017269_protein_group_present",
            "protein_accessions": "pxd017269_protein_accessions",
            "n_samples": "pxd017269_n_samples",
            "n_detected": "pxd017269_n_detected",
            "detection_fraction": "pxd017269_detection_fraction",
            "median_log2_intensity_detected": "pxd017269_median_log2_intensity_detected",
            "archival_group_kruskal_p": "pxd017269_archival_group_kruskal_p",
            "archival_group_kruskal_q_bh_candidates": "pxd017269_archival_group_kruskal_q_bh_candidates",
        }
    )
    pxd46 = pxd046999_presence(candidates)
    evidence = candidates.merge(pxd2, on="gene", how="left").merge(
        pxd17, on="gene", how="left"
    ).merge(pxd46, on="gene", how="left")
    evidence["pxd002137_adjusted_direction_match"] = (
        np.sign(evidence["age_sex_adjusted_log2_effect"])
        == evidence["expected_direction"]
    )
    evidence["public_protein_evidence_tier"] = evidence.apply(evidence_tier, axis=1)
    tier_order = {
        "A_directional_protein_replication_and_ffpe_detectability": 0,
        "B_directional_support_and_ffpe_detectability": 1,
        "C_incomplete_public_protein_support": 2,
        "C_transcript_locked_not_detected_in_pxd002137": 3,
        "D_direction_conflict": 4,
        "context_only_not_locked": 5,
    }
    evidence["tier_sort"] = evidence["public_protein_evidence_tier"].map(tier_order)
    evidence = evidence.sort_values(
        ["tier_sort", "age_sex_adjusted_q_bh_candidates", "direction_stability"],
        ascending=[True, True, False],
        na_position="last",
    ).drop(columns="tier_sort")

    candidates.to_csv(OUT_DIR / "candidate_inventory.tsv", sep="\t", index=False)
    pxd2.to_csv(OUT_DIR / "pxd002137_candidate_tests.tsv", sep="\t", index=False)
    pxd17.to_csv(OUT_DIR / "pxd017269_ffpe_detectability.tsv", sep="\t", index=False)
    pxd46.to_csv(OUT_DIR / "pxd046999_dvp_presence.tsv", sep="\t", index=False)
    evidence.to_csv(OUT_DIR / "candidate_public_protein_evidence_matrix.tsv", sep="\t", index=False)
    source_manifest().to_csv(OUT_DIR / "source_file_manifest.tsv", sep="\t", index=False)

    tier_counts = evidence["public_protein_evidence_tier"].value_counts()
    panel = evidence.loc[evidence["panel_role"].ne("locked_signature_candidate")]
    lines = [
        "Public colorectal adenoma protein triangulation",
        "================================================",
        f"Candidates audited: {len(evidence)} (locked signature: {int(evidence['locked_signature_member'].sum())}).",
        "PXD002137: 16 adenomas, 8 tumor-adjacent normal mucosae, 8 CRCs; adenoma-versus-normal is unpaired.",
        "PXD017269: FFPE detectability only (118 adenomas); no normal comparator.",
        "PXD046999: nine-patient DVP detectability/risk-group enrichment only; no normal comparator.",
        "",
        "Evidence tiers",
    ]
    for tier, count in tier_counts.items():
        lines.append(f"- {tier}: {count}")
    lines.extend(["", "Current panel and prespecified alternative"])
    for row in panel.itertuples(index=False):
        lines.append(
            f"- {row.gene} ({row.panel_role}): tier={row.public_protein_evidence_tier}; "
            f"PXD002137 adjusted effect={row.age_sex_adjusted_log2_effect:.3f}; "
            f"PXD017269 detection={row.pxd017269_detection_fraction:.3f}; "
            f"PXD046999 DVP detected={row.detected_in_nine_patient_dvp}."
        )
    lines.extend(
        [
            "",
            "Interpretation boundary",
            "These results prioritize assay candidates; they do not establish causal function, diagnostic thresholds, or clinical utility.",
            "PXD017269 and PXD046999 contribute detectability evidence, not independent normal-versus-adenoma replication.",
        ]
    )
    (OUT_DIR / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
