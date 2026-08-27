# Frozen analysis contract: epithelial-state decomposition v1

**Frozen on:** 27 August 2026, before outcome computation  
**Status:** post hoc mechanistic analysis with prospectively frozen implementation choices  
**Purpose:** distinguish epithelial-state composition from expression changes within matched deposited epithelial states.

## 1. Immutable inputs

| Input | SHA-256 |
|---|---|
| `data_sources/Chen_Cell_2021_CELLxGENE/chen_validation_epithelial.h5ad` | `f51bf02885d637073743093844c1db00a173196f0e1184b37f39a5118468668d` |
| `data_sources/Chen_Cell_2021_CELLxGENE/chen_discovery_epithelial.h5ad` | `2e2b5aa4909fd314ea0c6dba17f6ced16de95c1bb3317e69a3edfe2caff17440` |
| `release/colorectal-adenoma-wnt-signature/data/signatures/core_287_genes.tsv` | `c649ea65e6f32c667d347872bd16ec4712ef12202e359d680052b44def9cddf2` |
| `release/colorectal-adenoma-wnt-signature/data/signatures/signature_12_genes.tsv` | `91e564f42bdc5fb0188605b76116bcfb126e5102d62a33d6890fa67018928380` |

No gene membership, arm, weight, cell-state label or route mapping may be altered after outcome inspection.

## 2. Analysis hierarchy

- **Primary dataset:** Chen held-out validation epithelial partition.
- **Replication dataset:** Chen discovery epithelial partition; explicitly selection-aware because the core was defined there.
- **Primary contrast:** conventional adenoma minus normal.
- **Specificity contrast:** conventional adenoma minus serrated lesion.
- **Primary score:** frozen 287-gene core.
- **Compact replication score:** frozen 12-gene signature.
- **Inferential unit:** donor; specimens are nested measurements and cells are measurement units only.

## 3. Fixed mappings

### Lesion route

- Normal: `NL`.
- Conventional adenoma: `TA`, `TV`, `TVA`.
- Serrated: `HP`, `SSL`.
- `UNC` and unrecognised labels are excluded from group contrasts but retained in the data audit.

### Cell states

All deposited `Cell_Type` values are retained in the primary all-state analysis. The canonical-state sensitivity excludes only disease-labelled `ASC` (neoplastic) and `SSC` (abnormal), retaining `ABS`, `GOB`, `TAC`, `STM`, `CT`, `TUF` and `EE`.

## 4. Score calculation

For each Chen partition separately:

1. Match available genes by exact symbol.
2. Standardise each matched gene across all epithelial cells in that partition without using lesion labels: \(z_{ij}=(x_{ij}-\bar{x}_j)/s_j\).
3. Calculate each cell score as the unweighted mean of up-arm z scores minus the unweighted mean of down-arm z scores.
4. Report matched gene counts for both scores; do not replace missing genes.
5. Aggregate scores to specimen-by-cell-state means. A stratum is eligible for within-state analysis at the prespecified minimum cell threshold.

Primary minimum: 20 cells per specimen-state. Sensitivities: 10 and 50 cells.

The minimum-cell rule applies to separately reported **within-state effect estimates**. The exact decomposition uses every scored epithelial cell: a state absent from a specimen has proportion and contribution zero, whereas an observed state contributes its observed mean regardless of stratum size. This preserves the identity between the donor-balanced overall score and the sum of state contributions. As a descriptive stability audit, decomposition is additionally repeated after retaining only specimen-state strata meeting 10, 20 or 50 cells, renormalising proportions among retained cells and reporting the retained-cell fraction; these restricted-population results cannot replace the all-cell primary decomposition.

## 5. Donor balancing

Within each donor, lesion route and cell state, eligible specimens are averaged with equal specimen weight. For decomposition, each specimen contributes a cell-state proportion \(p_{dsk}\), mean score \(\mu_{dsk}\), and contribution \(c_{dsk}=p_{dsk}\mu_{dsk}\). Specimens are averaged equally within donor-route before group aggregation.

No cell-level P value is permitted.

## 6. Within-state effects

For every state and score:

- effect = donor-balanced group mean in conventional adenoma minus normal;
- uncertainty = percentile 95% interval from 5,000 whole-donor bootstrap samples;
- seed = `20260827`;
- donors are sampled with replacement and all observations belonging to a sampled donor are retained together;
- paired donor differences and an exact/approximate Wilcoxon signed-rank P value are reported when at least three donors contain both routes and the eligible cell state;
- Benjamini–Hochberg adjustment is applied across states within one dataset, score and threshold.

Effect sizes and intervals, not dichotomised significance, govern interpretation.

## 7. Differential abundance

For each specimen, state proportion = state cell count / total epithelial-cell count.

- Primary transformation: \(\arcsin(\sqrt p)\).
- Sensitivity transformation: empirical logit \(\log[(n_k+0.5)/(n-n_k+0.5)]\).
- Model: transformed proportion regressed on conventional-adenoma status versus normal, with donor-clustered sandwich standard errors.
- Absolute mean proportion difference and donor/specimen counts are always reported.
- BH adjustment is performed across states within each dataset and transformation.

This is described as a propeller-style transformed-proportion analysis, not as use of the propeller software package.

## 8. Exact decomposition

For group \(g\) and state \(k\):

- \(p_{gk}\) is the mean donor-route proportion;
- \(c_{gk}\) is the mean donor-route contribution;
- \(\mu_{gk}=c_{gk}/p_{gk}\).

The conventional-adenoma-minus-normal contrast is:

\[
\Delta = \sum_k(p_{Ak}-p_{Nk})(\mu_{Ak}+\mu_{Nk})/2
       + \sum_k(\mu_{Ak}-\mu_{Nk})(p_{Ak}+p_{Nk})/2.
\]

The first term is labelled **composition** and the second **within-state expression**. Signed state-specific contributions are retained. Whole-donor bootstrap resampling (5,000 samples, seed `20260827`) is used for intervals.

Numerical closure gate: `abs(total - composition - within_state) <= 1e-10` for every non-bootstrap and bootstrap calculation, allowing only a documented machine-precision exception.

Counterfactual summaries are also reported:

- composition-only shift holding normal state means fixed;
- expression-only shift using the pooled mean composition;
- group means under the pooled composition.

These are descriptive standardisations, not causal counterfactuals.

## 9. Prespecified sensitivities

1. all nine versus seven canonical states;
2. minimum 10, 20 and 50 cells per specimen-state for within-state estimates;
3. 287-gene versus 12-gene scores;
4. held-out validation versus selection-aware discovery;
5. conventional adenoma versus serrated specificity;
6. presence or absence of a fixed proliferation-control score using `MKI67`, `TOP2A`, `PCNA`, `MCM2`, `MCM5`, `TYMS`, `UBE2C` and `CENPF` where available.

## 10. Success and failure gates

The analysis may support a main-text statement that the programme is not purely compositional only if:

1. the held-out validation within-state component is positive;
2. its whole-donor bootstrap interval excludes zero or the same qualitative conclusion is independently supported by multiple major cell states with compatible intervals;
3. the canonical-state sensitivity preserves direction;
4. the 12-gene score is directionally concordant;
5. discovery and validation decompositions are directionally concordant;
6. exact numerical closure passes.

If these conditions are not met, the result is reported as mixed/equivocal and is not elevated to a headline biological claim.

## 11. Prohibited claims

- No causal mediation.
- No inferred lineage transition.
- No assertion that dissociation-derived proportions reproduce in situ abundance.
- No pharmacodynamic, diagnostic or clinical-performance claim.
- No result-dependent gene, cell-state or threshold selection.
