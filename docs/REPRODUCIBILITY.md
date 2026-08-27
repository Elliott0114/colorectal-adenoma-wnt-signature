# Reproducibility and claim boundaries

## Inferential units

Human analyses use the patient, donor or patient cluster as the independent
unit. Perturbation analyses use donors, organoid models, edited clones, cell
lines or oncogenic backgrounds. Spatial analyses use tissue sections. Cells,
nuclei, technical libraries and spots do not add inferential degrees of
freedom.

## Information flow

Validation outcomes were unavailable during all three definition steps:

```text
33,698 assayed features
  -> 6,127 expressed and non-technical genes
  -> 287 error-controlled genes
  -> 62 label-blind portable candidates
  -> 12 frozen genes selected from discovery-only fidelity
  -> held-out and external evaluation
```

The compact signature is benchmarked against 10,000 balanced matched
alternatives. It is reproducible but not uniquely superior; that negative
boundary is retained in the manuscript and supplementary source data.

## Evidence vocabulary

Use the following terms consistently when interpreting repository outputs:

| Evidence | Supported wording | Wording not supported |
|---|---|---|
| 287-gene discovery result | threshold-defined epithelial state or core | complete biological programme boundary |
| 12-gene result | frozen research signature or compact representation | diagnostic panel or validated biomarker |
| Held-out Chen data | held-out subset from the source study | independent cohort |
| Five GEO cohorts | independent transcriptomic replication | prospective clinical validation |
| Cell-state decomposition | accounting contributions from recovered-cell composition and within-state expression | causal mediation, lineage conversion or unbiased in situ abundance |
| RNA-ATAC | regulatory association or concordance | direct transcription-factor binding |
| CRC Atlas | cross-sectional recurrence across states | longitudinal progression or lineage continuity |
| Empirical perturbations | signed response support | universal causal mechanism |
| GenKI virtual deletion | unsigned network-coupling support | experimental knockout or causal proof |
| Spatial and protein layers | tissue localisation and protein anchors | protein implementation of the 12-gene signature |

## Virtual deletion

The GenKI analysis was prespecified as a validation layer. It used 1,664
donor-balanced held-out adenoma epithelial cells from 13 donors, 2,198 genes,
two model seeds and 10,000 matched null sets. The reported aggregate endpoints
and target-level results are stored in
`results/virtual_knockout_validation_v2_9/`.

GenKI distances are unsigned embedding perturbation magnitudes. Consequently,
virtual deletion is never used to infer the direction of biological change.
The dual-seed rerun reproduced all 13 reported tabular outputs exactly after
normalising non-scientific timing fields.

## Cell-state decomposition

The primary analysis uses the isolated Chen held-out partition. Cell scores are
summarised first within specimen–state strata; donors are the independent units
for effects and whole-donor bootstrap intervals. Differential abundance uses a
propeller-style arcsine-square-root model with donor-clustered inference. The
exact symmetric decomposition is evaluated for all deposited epithelial states
and after exclusion of the disease-labelled neoplastic and abnormal states.

The decomposition is an accounting identity. It does not identify causal
mediation or reconstruct lineage transitions. Dissociation-derived proportions
refer to recovered cells and are not treated as unbiased in situ abundance.
An independent base-R implementation reproduced all 16 reported identities
with maximum absolute error 2.22 × 10^-16.

## Reproduction tiers

1. **Release integrity**: `make verify` checks gene counts, arms, checksums,
   required outputs, sensitive-content patterns and file-size limits.
2. **Figure reproduction**: `make figures` regenerates all figures from included
   derived result tables without downloading raw data and audits the final
   seven-main-figure package.
3. **Full reanalysis**: acquire public raw data, create the two environments and
   run the stages listed in `workflow/pipeline.tsv`.
4. **Virtual-deletion rerun**: `make virtual-deletion-audit` reruns GenKI and
   compares the scientific outputs with the frozen files.
