#!/usr/bin/env python3
"""Decompose the frozen epithelial programme into composition and within-state effects.

The held-out Chen validation partition is primary. Cells are measurement units;
all estimates and resampling operate on specimens or whole donors.
"""

from __future__ import annotations

import hashlib
import json
import platform
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data_sources" / "Chen_Cell_2021_CELLxGENE"
STANDALONE_SIGNATURE_DIR = ROOT / "data" / "signatures"
WORKSPACE_SIGNATURE_DIR = (
    ROOT / "release" / "colorectal-adenoma-wnt-signature" / "data" / "signatures"
)
SIGNATURE_DIR = (
    STANDALONE_SIGNATURE_DIR
    if STANDALONE_SIGNATURE_DIR.is_dir()
    else WORKSPACE_SIGNATURE_DIR
)
OUT_DIR = ROOT / "results" / "cell_state_decomposition_v1"
CONTRACT_PATH = ROOT / "analysis" / "contracts" / "cell_state_decomposition_v1_2026-08-27.md"

DATASETS = {
    "validation": DATA_DIR / "chen_validation_epithelial.h5ad",
    "discovery": DATA_DIR / "chen_discovery_epithelial.h5ad",
}
SIGNATURE_FILES = {
    "core_287": SIGNATURE_DIR / "core_287_genes.tsv",
    "signature_12": SIGNATURE_DIR / "signature_12_genes.tsv",
}

ROUTE_MAP = {
    "NL": "normal",
    "TA": "conventional_adenoma",
    "TV": "conventional_adenoma",
    "TVA": "conventional_adenoma",
    "HP": "serrated",
    "SSL": "serrated",
    "UNC": "uncertain",
}

CELL_LABELS = {
    "ABS": "Absorptive",
    "GOB": "Goblet",
    "ASC": "Neoplastic",
    "TAC": "Transit-amplifying",
    "SSC": "Abnormal",
    "STM": "Crypt stem",
    "CT": "Other colonic epithelium",
    "TUF": "Tuft",
    "EE": "Enteroendocrine",
}
ALL_STATES = ["ABS", "GOB", "ASC", "TAC", "SSC", "STM", "CT", "TUF", "EE"]
CANONICAL_STATES = ["ABS", "GOB", "TAC", "STM", "CT", "TUF", "EE"]
PROLIFERATION_GENES = ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "TYMS", "UBE2C", "CENPF"]
SCORES = ["core_287", "signature_12", "proliferation_control"]
DECOMPOSITION_SCORES = ["core_287", "signature_12"]
WITHIN_THRESHOLDS = [10, 20, 50]
PRIMARY_THRESHOLD = 20
BOOTSTRAPS = 5000
SEED = 20260827
COMPARISONS = [
    ("conventional_vs_normal", "conventional_adenoma", "normal"),
    ("conventional_vs_serrated", "conventional_adenoma", "serrated"),
]


def decode(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


def read_categorical(group: h5py.Group) -> np.ndarray:
    codes = group["codes"][()]
    categories = np.array([decode(v) for v in group["categories"][()]], dtype=object)
    values = np.empty(len(codes), dtype=object)
    valid = codes >= 0
    values[valid] = categories[codes[valid]]
    values[~valid] = None
    return values


def read_col(root: h5py.Group, key: str) -> np.ndarray:
    obj = root[key]
    if isinstance(obj, h5py.Group):
        return read_categorical(obj)
    values = obj[()]
    if values.dtype.kind in "SOU":
        return np.array([decode(v) for v in values], dtype=object)
    return values


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> pd.Series:
    output = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return output
    order = np.argsort(valid.to_numpy())
    ranked = valid.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    output.loc[valid.index[order]] = adjusted
    return output


def load_gene_sets() -> dict[str, dict[str, list[str]]]:
    core = pd.read_csv(SIGNATURE_FILES["core_287"], sep="\t")
    compact = pd.read_csv(SIGNATURE_FILES["signature_12"], sep="\t")
    if core["gene"].duplicated().any() or compact["gene"].duplicated().any():
        raise ValueError("Signature files must contain unique gene symbols")
    if set(core["arm"].dropna()) != {"up", "down"}:
        raise ValueError("The 287-gene core must contain up and down arms")
    if set(compact["arm"].dropna()) != {"up", "down"}:
        raise ValueError("The 12-gene signature must contain up and down arms")
    return {
        "core_287": {
            "up": core.loc[core["arm"].eq("up"), "gene"].astype(str).tolist(),
            "down": core.loc[core["arm"].eq("down"), "gene"].astype(str).tolist(),
        },
        "signature_12": {
            "up": compact.loc[compact["arm"].eq("up"), "gene"].astype(str).tolist(),
            "down": compact.loc[compact["arm"].eq("down"), "gene"].astype(str).tolist(),
        },
        "proliferation_control": {"up": PROLIFERATION_GENES, "down": []},
    }


def load_partition(
    path: Path, dataset: str, gene_sets: dict[str, dict[str, list[str]]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    requested = sorted(
        {
            gene
            for definition in gene_sets.values()
            for arm in ("up", "down")
            for gene in definition[arm]
        }
    )
    with h5py.File(path, "r") as handle:
        genes = read_col(handle["var"], "feature_name")
        gene_to_index = {str(gene): index for index, gene in enumerate(genes)}
        present = [gene for gene in requested if gene in gene_to_index]
        ordered_pairs = sorted((gene_to_index[gene], gene) for gene in present)
        indices = np.array([index for index, _ in ordered_pairs], dtype=np.int64)
        ordered_genes = [gene for _, gene in ordered_pairs]
        expression = handle["X"][:, indices].astype(np.float64)

        means = expression.mean(axis=0)
        standard_deviations = expression.std(axis=0, ddof=0)
        standard_deviations[standard_deviations == 0] = 1.0
        z = (expression - means) / standard_deviations

        obs = pd.DataFrame(
            {
                "dataset": dataset,
                "specimen_id": read_col(handle["obs"], "HTAN Specimen ID"),
                "donor_id": read_col(handle["obs"], "donor_id"),
                "polyp_type": read_col(handle["obs"], "Polyp_Type"),
                "sample_classification": read_col(handle["obs"], "Sample_Classification"),
                "cell_type": read_col(handle["obs"], "Cell_Type"),
            }
        )

    obs["route_group"] = obs["polyp_type"].map(ROUTE_MAP).fillna("other")
    gene_to_position = {gene: position for position, gene in enumerate(ordered_genes)}
    availability = []
    for score_name, definition in gene_sets.items():
        up = [gene for gene in definition["up"] if gene in gene_to_position]
        down = [gene for gene in definition["down"] if gene in gene_to_position]
        missing_up = [gene for gene in definition["up"] if gene not in gene_to_position]
        missing_down = [gene for gene in definition["down"] if gene not in gene_to_position]
        up_positions = [gene_to_position[gene] for gene in up]
        down_positions = [gene_to_position[gene] for gene in down]
        if not up_positions:
            raise ValueError(f"No up-arm genes available for {score_name} in {dataset}")
        score = z[:, up_positions].mean(axis=1)
        if down_positions:
            score = score - z[:, down_positions].mean(axis=1)
        obs[f"score__{score_name}"] = score
        availability.append(
            {
                "dataset": dataset,
                "score": score_name,
                "n_up_requested": len(definition["up"]),
                "n_up_present": len(up),
                "n_down_requested": len(definition["down"]),
                "n_down_present": len(down),
                "missing_up": ",".join(missing_up),
                "missing_down": ",".join(missing_down),
            }
        )
    return obs, pd.DataFrame(availability)


def summarize_specimen_states(cells: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_columns = [f"score__{score}" for score in SCORES]
    state = (
        cells.groupby(
            ["dataset", "specimen_id", "donor_id", "polyp_type", "route_group", "cell_type"],
            observed=True,
        )[score_columns]
        .agg(["size", "mean", "median"])
        .reset_index()
    )
    state.columns = ["__".join([str(x) for x in col if x]) for col in state.columns.to_flat_index()]
    state = state.rename(columns={"score__core_287__size": "n_cells"})
    for score in SCORES[1:]:
        redundant = f"score__{score}__size"
        if redundant in state:
            state = state.drop(columns=redundant)
    state["cell_type_label"] = state["cell_type"].map(CELL_LABELS).fillna(state["cell_type"])

    specimen = (
        cells.groupby(["dataset", "specimen_id"], observed=True)
        .agg(
            donor_id=("donor_id", "first"),
            polyp_type=("polyp_type", "first"),
            route_group=("route_group", "first"),
            n_cells_total=("cell_type", "size"),
        )
        .reset_index()
    )
    return state, specimen


def donor_route_state_values(state: pd.DataFrame, score: str, threshold: int) -> pd.DataFrame:
    value = f"score__{score}__mean"
    eligible = state.loc[state["n_cells"].ge(threshold)].copy()
    return (
        eligible.groupby(["dataset", "donor_id", "route_group", "cell_type"], observed=True)
        .agg(
            mean_score=(value, "mean"),
            n_specimens=("specimen_id", "nunique"),
            total_cells=("n_cells", "sum"),
        )
        .reset_index()
    )


def bootstrap_group_difference(
    frame: pd.DataFrame, group_a: str, group_b: str, rng: np.random.Generator
) -> np.ndarray:
    donors = sorted(frame["donor_id"].astype(str).unique())
    matrix = (
        frame.assign(donor_id=frame["donor_id"].astype(str))
        .pivot_table(index="donor_id", columns="route_group", values="mean_score", aggfunc="mean")
        .reindex(index=donors, columns=[group_a, group_b])
    )
    values = matrix.to_numpy(dtype=float)
    present = np.isfinite(values).astype(float)
    filled = np.nan_to_num(values, nan=0.0)
    weights = rng.multinomial(len(donors), np.repeat(1 / len(donors), len(donors)), size=BOOTSTRAPS)
    denominator = weights @ present
    numerator = weights @ filled
    with np.errstate(invalid="ignore", divide="ignore"):
        means = numerator / denominator
    differences = means[:, 0] - means[:, 1]
    differences[(denominator == 0).any(axis=1)] = np.nan
    return differences


def paired_wilcoxon(frame: pd.DataFrame, group_a: str, group_b: str) -> tuple[int, float, float]:
    paired = frame.pivot_table(index="donor_id", columns="route_group", values="mean_score", aggfunc="mean")
    if group_a not in paired or group_b not in paired:
        return 0, np.nan, np.nan
    paired = paired[[group_a, group_b]].dropna()
    if len(paired) < 3:
        return len(paired), np.nan, np.nan
    difference = paired[group_a] - paired[group_b]
    if np.allclose(difference, 0):
        return len(paired), float(np.median(difference)), 1.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = wilcoxon(difference, alternative="two-sided", zero_method="wilcox", method="auto")
    return len(paired), float(np.median(difference)), float(result.pvalue)


def within_state_effects(state: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for score in SCORES:
        for threshold in WITHIN_THRESHOLDS:
            donor_values = donor_route_state_values(state, score, threshold)
            for dataset, dataset_frame in donor_values.groupby("dataset", observed=True):
                for comparison, group_a, group_b in COMPARISONS:
                    for cell_type in ALL_STATES:
                        frame = dataset_frame.loc[
                            dataset_frame["cell_type"].eq(cell_type)
                            & dataset_frame["route_group"].isin([group_a, group_b])
                        ].copy()
                        counts = frame.groupby("route_group", observed=True)["donor_id"].nunique()
                        if counts.get(group_a, 0) == 0 or counts.get(group_b, 0) == 0:
                            continue
                        observed = (
                            frame.loc[frame["route_group"].eq(group_a), "mean_score"].mean()
                            - frame.loc[frame["route_group"].eq(group_b), "mean_score"].mean()
                        )
                        rng = np.random.default_rng(SEED)
                        boot = bootstrap_group_difference(frame, group_a, group_b, rng)
                        boot_valid = boot[np.isfinite(boot)]
                        n_pairs, paired_median, paired_p = paired_wilcoxon(frame, group_a, group_b)
                        rows.append(
                            {
                                "dataset": dataset,
                                "score": score,
                                "min_cells_per_specimen_state": threshold,
                                "primary_threshold": threshold == PRIMARY_THRESHOLD,
                                "comparison": comparison,
                                "group_a": group_a,
                                "group_b": group_b,
                                "cell_type": cell_type,
                                "cell_type_label": CELL_LABELS.get(cell_type, cell_type),
                                "n_donors_a": int(counts.get(group_a, 0)),
                                "n_donors_b": int(counts.get(group_b, 0)),
                                "mean_a": float(frame.loc[frame["route_group"].eq(group_a), "mean_score"].mean()),
                                "mean_b": float(frame.loc[frame["route_group"].eq(group_b), "mean_score"].mean()),
                                "effect_a_minus_b": float(observed),
                                "bootstrap_ci_low": float(np.quantile(boot_valid, 0.025)),
                                "bootstrap_ci_high": float(np.quantile(boot_valid, 0.975)),
                                "bootstrap_p_directional": float(
                                    2 * min(np.mean(boot_valid <= 0), np.mean(boot_valid >= 0))
                                ),
                                "n_paired_donors": n_pairs,
                                "paired_median_difference": paired_median,
                                "paired_wilcoxon_p": paired_p,
                            }
                        )
    output = pd.DataFrame(rows)
    if output.empty:
        return output
    group_columns = ["dataset", "score", "min_cells_per_specimen_state", "comparison"]
    output["paired_wilcoxon_q"] = output.groupby(group_columns, observed=True, group_keys=False)[
        "paired_wilcoxon_p"
    ].apply(bh_adjust)
    return output.sort_values(group_columns + ["cell_type"])


def specimen_state_grid(state: pd.DataFrame, specimen: pd.DataFrame) -> pd.DataFrame:
    state_counts = state[["dataset", "specimen_id", "cell_type", "n_cells"]].copy()
    grids = []
    for dataset, specimen_frame in specimen.groupby("dataset", observed=True):
        grid = pd.MultiIndex.from_product(
            [specimen_frame["specimen_id"].astype(str), ALL_STATES],
            names=["specimen_id", "cell_type"],
        ).to_frame(index=False)
        grid.insert(0, "dataset", dataset)
        grids.append(grid)
    full = pd.concat(grids, ignore_index=True)
    full = full.merge(state_counts, on=["dataset", "specimen_id", "cell_type"], how="left")
    full["n_cells"] = full["n_cells"].fillna(0).astype(int)
    full = full.merge(specimen, on=["dataset", "specimen_id"], how="left")
    full["proportion"] = full["n_cells"] / full["n_cells_total"]
    full["cell_type_label"] = full["cell_type"].map(CELL_LABELS)
    return full


def fit_clustered_ols(frame: pd.DataFrame, outcome: str) -> dict[str, float]:
    model_frame = frame[[outcome, "route_group", "donor_id"]].dropna().copy()
    model_frame["adenoma"] = model_frame["route_group"].eq("conventional_adenoma").astype(float)
    design = sm.add_constant(model_frame[["adenoma"]], has_constant="add")
    try:
        fit = sm.OLS(model_frame[outcome].astype(float), design.astype(float)).fit(
            cov_type="cluster", cov_kwds={"groups": model_frame["donor_id"].astype(str)}
        )
        interval = fit.conf_int().loc["adenoma"]
        return {
            "coefficient": float(fit.params["adenoma"]),
            "standard_error": float(fit.bse["adenoma"]),
            "ci_low": float(interval.iloc[0]),
            "ci_high": float(interval.iloc[1]),
            "p_value": float(fit.pvalues["adenoma"]),
        }
    except (ValueError, np.linalg.LinAlgError):
        return {key: np.nan for key in ["coefficient", "standard_error", "ci_low", "ci_high", "p_value"]}


def differential_abundance(grid: pd.DataFrame) -> pd.DataFrame:
    rows = []
    primary = grid.loc[grid["route_group"].isin(["normal", "conventional_adenoma"])].copy()
    for dataset, dataset_frame in primary.groupby("dataset", observed=True):
        for cell_type in ALL_STATES:
            frame = dataset_frame.loc[dataset_frame["cell_type"].eq(cell_type)].copy()
            frame["arcsine_sqrt"] = np.arcsin(np.sqrt(frame["proportion"].clip(0, 1)))
            frame["empirical_logit"] = np.log(
                (frame["n_cells"] + 0.5) / (frame["n_cells_total"] - frame["n_cells"] + 0.5)
            )
            donor_route = (
                frame.groupby(["donor_id", "route_group"], observed=True)["proportion"].mean().reset_index()
            )
            donor_means = donor_route.groupby("route_group", observed=True)["proportion"].mean()
            for transformation in ["arcsine_sqrt", "empirical_logit"]:
                estimate = fit_clustered_ols(frame, transformation)
                rows.append(
                    {
                        "dataset": dataset,
                        "cell_type": cell_type,
                        "cell_type_label": CELL_LABELS.get(cell_type, cell_type),
                        "transformation": transformation,
                        "n_specimens": frame["specimen_id"].nunique(),
                        "n_donors": frame["donor_id"].nunique(),
                        "n_donors_adenoma": donor_route.loc[
                            donor_route["route_group"].eq("conventional_adenoma"), "donor_id"
                        ].nunique(),
                        "n_donors_normal": donor_route.loc[
                            donor_route["route_group"].eq("normal"), "donor_id"
                        ].nunique(),
                        "mean_proportion_adenoma": float(donor_means.get("conventional_adenoma", np.nan)),
                        "mean_proportion_normal": float(donor_means.get("normal", np.nan)),
                        "absolute_proportion_difference": float(
                            donor_means.get("conventional_adenoma", np.nan) - donor_means.get("normal", np.nan)
                        ),
                        **estimate,
                    }
                )
    output = pd.DataFrame(rows)
    output["q_value"] = output.groupby(["dataset", "transformation"], observed=True, group_keys=False)[
        "p_value"
    ].apply(bh_adjust)
    return output.sort_values(["dataset", "transformation", "cell_type"])


def decomposition_donor_base(
    state: pd.DataFrame, specimen: pd.DataFrame, score: str, states: list[str]
) -> pd.DataFrame:
    value_column = f"score__{score}__mean"
    selected = state.loc[state["cell_type"].isin(states)].copy()
    totals = (
        selected.groupby(["dataset", "specimen_id"], observed=True)["n_cells"]
        .sum()
        .rename("included_cells")
        .reset_index()
    )
    selected = selected.merge(totals, on=["dataset", "specimen_id"], how="left")
    selected["proportion"] = selected["n_cells"] / selected["included_cells"]
    selected["contribution"] = selected["proportion"] * selected[value_column]

    metadata = specimen[["dataset", "specimen_id", "donor_id", "route_group", "n_cells_total"]].merge(
        totals, on=["dataset", "specimen_id"], how="inner"
    )
    grid_parts = []
    for dataset, meta in metadata.groupby("dataset", observed=True):
        grid = pd.MultiIndex.from_product(
            [meta["specimen_id"].astype(str), states], names=["specimen_id", "cell_type"]
        ).to_frame(index=False)
        grid.insert(0, "dataset", dataset)
        grid_parts.append(grid)
    grid = pd.concat(grid_parts, ignore_index=True)
    values = selected[
        ["dataset", "specimen_id", "cell_type", "proportion", "contribution"]
    ]
    grid = grid.merge(values, on=["dataset", "specimen_id", "cell_type"], how="left")
    grid[["proportion", "contribution"]] = grid[["proportion", "contribution"]].fillna(0.0)
    grid = grid.merge(metadata, on=["dataset", "specimen_id"], how="left")
    grid["retained_fraction"] = grid["included_cells"] / grid["n_cells_total"]
    donor = (
        grid.groupby(["dataset", "donor_id", "route_group", "cell_type"], observed=True)
        .agg(
            proportion=("proportion", "mean"),
            contribution=("contribution", "mean"),
            n_specimens=("specimen_id", "nunique"),
            retained_fraction=("retained_fraction", "mean"),
        )
        .reset_index()
    )
    return donor


def arrays_from_donor_base(
    donor: pd.DataFrame, dataset: str, states: list[str], group_a: str, group_b: str
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    frame = donor.loc[
        donor["dataset"].eq(dataset) & donor["route_group"].isin([group_a, group_b])
    ].copy()
    donors = sorted(frame["donor_id"].astype(str).unique())
    groups = [group_a, group_b]
    p = np.full((len(donors), 2, len(states)), np.nan, dtype=float)
    c = np.full_like(p, np.nan)
    present = np.zeros((len(donors), 2), dtype=float)
    donor_index = {donor_id: i for i, donor_id in enumerate(donors)}
    group_index = {group: i for i, group in enumerate(groups)}
    state_index = {state: i for i, state in enumerate(states)}
    for row in frame.itertuples(index=False):
        i = donor_index[str(row.donor_id)]
        j = group_index[row.route_group]
        k = state_index[row.cell_type]
        p[i, j, k] = float(row.proportion)
        c[i, j, k] = float(row.contribution)
        present[i, j] = 1.0
    for i in range(len(donors)):
        for j in range(2):
            if present[i, j]:
                if not np.isclose(np.nansum(p[i, j]), 1.0, atol=1e-10):
                    raise AssertionError("Donor-route state proportions do not sum to one")
                p[i, j] = np.nan_to_num(p[i, j], nan=0.0)
                c[i, j] = np.nan_to_num(c[i, j], nan=0.0)
    return donors, p, c, present


def group_arrays(
    weights: np.ndarray, p: np.ndarray, c: np.ndarray, present: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    denominator = weights @ present
    p_group = np.empty((weights.shape[0], 2, p.shape[2]), dtype=float)
    c_group = np.empty_like(p_group)
    for group in range(2):
        p_numerator = weights @ np.nan_to_num(p[:, group, :], nan=0.0)
        c_numerator = weights @ np.nan_to_num(c[:, group, :], nan=0.0)
        p_group[:, group, :] = np.divide(
            p_numerator,
            denominator[:, group, None],
            out=np.full_like(p_numerator, np.nan),
            where=denominator[:, group, None] > 0,
        )
        c_group[:, group, :] = np.divide(
            c_numerator,
            denominator[:, group, None],
            out=np.full_like(c_numerator, np.nan),
            where=denominator[:, group, None] > 0,
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        mu_group = c_group / p_group
    return p_group, c_group, mu_group


def decompose_arrays(
    p_group: np.ndarray, c_group: np.ndarray, mu_group: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if np.isnan(mu_group).any():
        raise ValueError("A cell state is absent from one comparison group; exact named-state decomposition is undefined")
    composition_by_state = (p_group[:, 0, :] - p_group[:, 1, :]) * (
        mu_group[:, 0, :] + mu_group[:, 1, :]
    ) / 2
    within_by_state = (mu_group[:, 0, :] - mu_group[:, 1, :]) * (
        p_group[:, 0, :] + p_group[:, 1, :]
    ) / 2
    composition = composition_by_state.sum(axis=1)
    within = within_by_state.sum(axis=1)
    total = c_group[:, 0, :].sum(axis=1) - c_group[:, 1, :].sum(axis=1)
    closure = total - composition - within
    return total, composition, within, composition_by_state, within_by_state, closure


def decomposition_analysis(
    state: pd.DataFrame, specimen: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summaries = []
    contributions = []
    bootstraps = []
    donor_bases = []
    for score in DECOMPOSITION_SCORES:
        for state_set, states in [("all_states", ALL_STATES), ("canonical_states", CANONICAL_STATES)]:
            donor = decomposition_donor_base(state, specimen, score, states)
            donor["score"] = score
            donor["state_set"] = state_set
            donor_bases.append(donor)
            for dataset in DATASETS:
                for comparison, group_a, group_b in COMPARISONS:
                    donors, p, c, present = arrays_from_donor_base(donor, dataset, states, group_a, group_b)
                    if present[:, 0].sum() == 0 or present[:, 1].sum() == 0:
                        continue
                    observed_weights = np.ones((1, len(donors)), dtype=int)
                    observed_p, observed_c, observed_mu = group_arrays(observed_weights, p, c, present)
                    total, composition, within, comp_state, within_state, closure = decompose_arrays(
                        observed_p, observed_c, observed_mu
                    )

                    rng = np.random.default_rng(SEED)
                    weights = rng.multinomial(
                        len(donors), np.repeat(1 / len(donors), len(donors)), size=BOOTSTRAPS
                    )
                    boot_p, boot_c, boot_mu = group_arrays(weights, p, c, present)
                    valid = np.isfinite(boot_mu).all(axis=(1, 2))
                    boot_total, boot_comp, boot_within, boot_comp_state, boot_within_state, boot_closure = (
                        decompose_arrays(boot_p[valid], boot_c[valid], boot_mu[valid])
                    )
                    if np.max(np.abs(boot_closure)) > 1e-10 or abs(float(closure[0])) > 1e-10:
                        raise AssertionError("Kitagawa decomposition failed numerical closure")

                    pooled_p = (observed_p[0, 0] + observed_p[0, 1]) / 2
                    standardised_a = float(np.sum(pooled_p * observed_mu[0, 0]))
                    standardised_b = float(np.sum(pooled_p * observed_mu[0, 1]))
                    boot_within_share = np.divide(
                        boot_within,
                        boot_total,
                        out=np.full_like(boot_within, np.nan),
                        where=np.abs(boot_total) > 1e-12,
                    )
                    summaries.append(
                        {
                            "dataset": dataset,
                            "score": score,
                            "state_set": state_set,
                            "comparison": comparison,
                            "group_a": group_a,
                            "group_b": group_b,
                            "n_unique_donors": len(donors),
                            "n_donors_a": int(present[:, 0].sum()),
                            "n_donors_b": int(present[:, 1].sum()),
                            "total_difference": float(total[0]),
                            "total_ci_low": float(np.quantile(boot_total, 0.025)),
                            "total_ci_high": float(np.quantile(boot_total, 0.975)),
                            "composition_component": float(composition[0]),
                            "composition_ci_low": float(np.quantile(boot_comp, 0.025)),
                            "composition_ci_high": float(np.quantile(boot_comp, 0.975)),
                            "within_state_component": float(within[0]),
                            "within_state_ci_low": float(np.quantile(boot_within, 0.025)),
                            "within_state_ci_high": float(np.quantile(boot_within, 0.975)),
                            "within_share_of_total": float(within[0] / total[0]) if abs(total[0]) > 1e-12 else np.nan,
                            "within_share_ci_low": float(np.nanquantile(boot_within_share, 0.025)),
                            "within_share_ci_high": float(np.nanquantile(boot_within_share, 0.975)),
                            "group_a_total_mean": float(observed_c[0, 0].sum()),
                            "group_b_total_mean": float(observed_c[0, 1].sum()),
                            "group_a_pooled_composition_mean": standardised_a,
                            "group_b_pooled_composition_mean": standardised_b,
                            "pooled_composition_difference": standardised_a - standardised_b,
                            "normal_means_composition_shift": float(
                                np.sum((observed_p[0, 0] - observed_p[0, 1]) * observed_mu[0, 1])
                            ),
                            "numerical_closure_error": float(closure[0]),
                            "n_valid_bootstraps": int(valid.sum()),
                        }
                    )
                    for index, cell_type in enumerate(states):
                        contributions.append(
                            {
                                "dataset": dataset,
                                "score": score,
                                "state_set": state_set,
                                "comparison": comparison,
                                "cell_type": cell_type,
                                "cell_type_label": CELL_LABELS.get(cell_type, cell_type),
                                "proportion_a": float(observed_p[0, 0, index]),
                                "proportion_b": float(observed_p[0, 1, index]),
                                "mean_score_a": float(observed_mu[0, 0, index]),
                                "mean_score_b": float(observed_mu[0, 1, index]),
                                "composition_contribution": float(comp_state[0, index]),
                                "composition_ci_low": float(np.quantile(boot_comp_state[:, index], 0.025)),
                                "composition_ci_high": float(np.quantile(boot_comp_state[:, index], 0.975)),
                                "within_state_contribution": float(within_state[0, index]),
                                "within_state_ci_low": float(np.quantile(boot_within_state[:, index], 0.025)),
                                "within_state_ci_high": float(np.quantile(boot_within_state[:, index], 0.975)),
                            }
                        )
                    for replicate, values in enumerate(zip(boot_total, boot_comp, boot_within, boot_closure), start=1):
                        bootstraps.append(
                            {
                                "dataset": dataset,
                                "score": score,
                                "state_set": state_set,
                                "comparison": comparison,
                                "bootstrap": replicate,
                                "total_difference": float(values[0]),
                                "composition_component": float(values[1]),
                                "within_state_component": float(values[2]),
                                "closure_error": float(values[3]),
                            }
                        )
    return (
        pd.DataFrame(summaries),
        pd.DataFrame(contributions),
        pd.DataFrame(bootstraps),
        pd.concat(donor_bases, ignore_index=True),
    )


def score_concordance(donor_base: pd.DataFrame) -> pd.DataFrame:
    totals = (
        donor_base.groupby(
            ["dataset", "score", "state_set", "donor_id", "route_group"], observed=True
        )["contribution"]
        .sum()
        .rename("total_score")
        .reset_index()
    )
    rows = []
    scopes = {
        "normal_and_conventional": ["normal", "conventional_adenoma"],
        "conventional_and_serrated": ["conventional_adenoma", "serrated"],
        "all_interpretable_routes": ["normal", "conventional_adenoma", "serrated"],
    }
    for (dataset, state_set), frame in totals.groupby(["dataset", "state_set"], observed=True):
        for scope, routes in scopes.items():
            selected = frame.loc[frame["route_group"].isin(routes)]
            wide = selected.pivot_table(
                index=["donor_id", "route_group"], columns="score", values="total_score", aggfunc="mean"
            ).dropna(subset=DECOMPOSITION_SCORES)
            if len(wide) < 3:
                continue
            rho, p_value = scipy.stats.spearmanr(wide["core_287"], wide["signature_12"])
            rows.append(
                {
                    "dataset": dataset,
                    "state_set": state_set,
                    "scope": scope,
                    "n_donor_route_units": len(wide),
                    "spearman_rho": float(rho),
                    "p_value": float(p_value),
                }
            )
    return pd.DataFrame(rows).sort_values(["dataset", "state_set", "scope"])


def make_summary(
    availability: pd.DataFrame,
    specimen: pd.DataFrame,
    within: pd.DataFrame,
    abundance: pd.DataFrame,
    decomposition: pd.DataFrame,
    concordance: pd.DataFrame,
) -> str:
    primary = decomposition.loc[
        decomposition["dataset"].eq("validation")
        & decomposition["score"].eq("core_287")
        & decomposition["state_set"].eq("all_states")
        & decomposition["comparison"].eq("conventional_vs_normal")
    ].iloc[0]
    canonical = decomposition.loc[
        decomposition["dataset"].eq("validation")
        & decomposition["score"].eq("core_287")
        & decomposition["state_set"].eq("canonical_states")
        & decomposition["comparison"].eq("conventional_vs_normal")
    ].iloc[0]
    compact = decomposition.loc[
        decomposition["dataset"].eq("validation")
        & decomposition["score"].eq("signature_12")
        & decomposition["state_set"].eq("all_states")
        & decomposition["comparison"].eq("conventional_vs_normal")
    ].iloc[0]
    discovery = decomposition.loc[
        decomposition["dataset"].eq("discovery")
        & decomposition["score"].eq("core_287")
        & decomposition["state_set"].eq("all_states")
        & decomposition["comparison"].eq("conventional_vs_normal")
    ].iloc[0]

    primary_within = within.loc[
        within["dataset"].eq("validation")
        & within["score"].eq("core_287")
        & within["min_cells_per_specimen_state"].eq(PRIMARY_THRESHOLD)
        & within["comparison"].eq("conventional_vs_normal")
    ].copy()
    well_supported_within = primary_within.loc[
        primary_within["n_donors_a"].ge(5) & primary_within["n_donors_b"].ge(5)
    ].sort_values("effect_a_minus_b", ascending=False)
    da = abundance.loc[
        abundance["dataset"].eq("validation")
        & abundance["transformation"].eq("arcsine_sqrt")
    ].sort_values("absolute_proportion_difference", ascending=False)
    concordance_row = concordance.loc[
        concordance["dataset"].eq("validation")
        & concordance["state_set"].eq("all_states")
        & concordance["scope"].eq("normal_and_conventional")
    ].iloc[0]

    def interval(row: pd.Series, value: str, low: str, high: str) -> str:
        return f"{row[value]:.3f} ({row[low]:.3f} to {row[high]:.3f})"

    lines = [
        "# Epithelial-state decomposition v1 — result summary",
        "",
        "## Primary held-out result",
        "",
        f"- Donor-balanced conventional-adenoma minus normal 287-gene score: {interval(primary, 'total_difference', 'total_ci_low', 'total_ci_high')}.",
        f"- Composition component: {interval(primary, 'composition_component', 'composition_ci_low', 'composition_ci_high')}.",
        f"- Within-state component: {interval(primary, 'within_state_component', 'within_state_ci_low', 'within_state_ci_high')}.",
        f"- Observed within-state share of the signed total: {primary['within_share_of_total']:.1%}. This ratio is descriptive because bootstrap ratios become unstable when a resampled total approaches zero.",
        f"- Exact numerical closure error: {primary['numerical_closure_error']:.3e}.",
        "",
        "## Prespecified robustness checks",
        "",
        f"- Canonical seven-state within-state component: {interval(canonical, 'within_state_component', 'within_state_ci_low', 'within_state_ci_high')}.",
        f"- Frozen 12-gene within-state component: {interval(compact, 'within_state_component', 'within_state_ci_low', 'within_state_ci_high')}.",
        f"- Donor-route 287-gene and 12-gene scores were concordant in held-out normal/conventional samples (Spearman rho={concordance_row['spearman_rho']:.3f}, n={int(concordance_row['n_donor_route_units'])}).",
        f"- Discovery-partition 287-gene within-state component: {interval(discovery, 'within_state_component', 'within_state_ci_low', 'within_state_ci_high')} (selection-aware replication).",
        "",
        "## Most informative state-resolved effects",
        "",
    ]
    for row in well_supported_within.itertuples(index=False):
        lines.append(
            f"- {row.cell_type_label}: effect {row.effect_a_minus_b:.3f}; 95% donor-bootstrap interval {row.bootstrap_ci_low:.3f} to {row.bootstrap_ci_high:.3f}; donors {row.n_donors_a}/{row.n_donors_b} (adenoma/normal)."
        )
    lines.extend(["", "## Largest abundance shifts in held-out validation", ""])
    for row in da.head(5).itertuples(index=False):
        lines.append(
            f"- {row.cell_type_label}: {row.mean_proportion_adenoma:.1%} versus {row.mean_proportion_normal:.1%}; absolute difference {row.absolute_proportion_difference:+.1%}; q={row.q_value:.3g}."
        )
    lines.extend(
        [
            "",
            "## Interpretation rule",
            "",
            "The analysis distinguishes an accounting contribution from altered epithelial-state abundance and an accounting contribution from different scores within deposited states. It does not establish causality, lineage conversion or unbiased in situ cell abundance.",
            "",
            "## Data audit",
            "",
            f"- Specimens analysed: {len(specimen)}.",
            f"- Frozen score availability rows: {len(availability)}.",
            "- Cells were never used as inferential replicates.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    gene_sets = load_gene_sets()
    cells = []
    availability = []
    for dataset, path in DATASETS.items():
        frame, available = load_partition(path, dataset, gene_sets)
        cells.append(frame)
        availability.append(available)
    all_cells = pd.concat(cells, ignore_index=True)
    availability_table = pd.concat(availability, ignore_index=True)
    state, specimen = summarize_specimen_states(all_cells)
    del all_cells

    grid = specimen_state_grid(state, specimen)
    within = within_state_effects(state)
    abundance = differential_abundance(grid)
    decomposition, contributions, bootstrap, donor_base = decomposition_analysis(state, specimen)
    concordance = score_concordance(donor_base)

    availability_table.to_csv(OUT_DIR / "gene_availability.tsv", sep="\t", index=False)
    specimen.to_csv(OUT_DIR / "specimen_inventory.tsv", sep="\t", index=False)
    state.to_csv(OUT_DIR / "specimen_cell_state_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    grid.to_csv(OUT_DIR / "specimen_cell_state_proportions.tsv.gz", sep="\t", index=False, compression="gzip")
    donor_base.to_csv(OUT_DIR / "donor_route_decomposition_inputs.tsv.gz", sep="\t", index=False, compression="gzip")
    within.to_csv(OUT_DIR / "within_cell_state_effects.tsv", sep="\t", index=False)
    abundance.to_csv(OUT_DIR / "differential_abundance_propeller_style.tsv", sep="\t", index=False)
    decomposition.to_csv(OUT_DIR / "decomposition_summary.tsv", sep="\t", index=False)
    contributions.to_csv(OUT_DIR / "decomposition_state_contributions.tsv", sep="\t", index=False)
    bootstrap.to_csv(OUT_DIR / "decomposition_bootstrap.tsv.gz", sep="\t", index=False, compression="gzip")
    concordance.to_csv(OUT_DIR / "core_compact_score_concordance.tsv", sep="\t", index=False)

    qa = pd.DataFrame(
        [
            {
                "check": "all_decompositions_close",
                "value": float(decomposition["numerical_closure_error"].abs().max()),
                "criterion": "<=1e-10",
                "pass": bool(decomposition["numerical_closure_error"].abs().max() <= 1e-10),
            },
            {
                "check": "all_bootstraps_close",
                "value": float(bootstrap["closure_error"].abs().max()),
                "criterion": "<=1e-10",
                "pass": bool(bootstrap["closure_error"].abs().max() <= 1e-10),
            },
            {
                "check": "primary_dataset_present",
                "value": int((specimen["dataset"] == "validation").sum()),
                "criterion": ">0",
                "pass": bool((specimen["dataset"] == "validation").any()),
            },
            {
                "check": "frozen_core_gene_count",
                "value": int(
                    availability_table.loc[
                        availability_table["score"].eq("core_287"), ["n_up_present", "n_down_present"]
                    ].iloc[0].sum()
                ),
                "criterion": "287",
                "pass": bool(
                    availability_table.loc[
                        availability_table["score"].eq("core_287"), ["n_up_present", "n_down_present"]
                    ].iloc[0].sum()
                    == 287
                ),
            },
        ]
    )
    qa.to_csv(OUT_DIR / "analysis_qa.tsv", sep="\t", index=False)
    if not qa["pass"].all():
        raise AssertionError("One or more prespecified QA gates failed")

    manifest = {
        "analysis": "cell_state_composition_decomposition_v1",
        "date": "2026-08-27",
        "contract": str(CONTRACT_PATH.relative_to(ROOT)),
        "seed": SEED,
        "bootstrap_replicates": BOOTSTRAPS,
        "python": platform.python_version(),
        "packages": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
            "h5py": h5py.__version__,
        },
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [*DATASETS.values(), *SIGNATURE_FILES.values(), Path(__file__)]
        },
        "outputs": sorted(path.name for path in OUT_DIR.iterdir() if path.is_file()),
    }
    (OUT_DIR / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (OUT_DIR / "summary.md").write_text(
        make_summary(availability_table, specimen, within, abundance, decomposition, concordance), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
