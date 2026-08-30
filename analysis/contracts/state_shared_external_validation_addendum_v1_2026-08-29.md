# Implementation addendum: external validation of the frozen state-shared programme and compact readout (v1)

**Date:** 29 August 2026
**Status:** FROZEN before calculating any external-cohort, FFPE or perturbation outcome for the new 1,843- or eight-gene definitions
**Parent contract:** `state_aware_program_rederivation_v1_2026-08-29.md`

## Frozen inputs

- The biological programme is the 1,843-gene strict state-shared core (884 up, 959 down).
- The compact readout is the discovery-only eight-gene balanced panel (four up, four down).
- Membership, direction and equal-arm weights cannot be changed using any result below.
- Existing external matrices, metadata mappings, patient clusters, contrast definitions and proliferation-control genes are reused unchanged from the audited v2.7 validation pipeline.

## Score construction

1. Within each cohort or experimental dataset, transform each measured gene to a sample-level z score without using disease labels.
2. Calculate each score as mean z(up) minus mean z(down), so unequal arm sizes cannot dominate.
3. The full-programme score uses every measurable frozen core gene. Report up- and down-arm coverage and do not impute missing genes.
4. The compact score requires all eight frozen genes. Missingness causes a dataset-level measurement failure; no substitution is permitted.

## Validation layers and outcomes

### Five independent sporadic-adenoma cohorts

- GSE8671, GSE50114, GSE41657, GSE40362 and GSE72820.
- Report per-cohort adenoma-versus-normal effects under the existing paired or unpaired design.
- Fit the existing one-stage patient-cluster model across cohorts, before and after proliferation adjustment, and repeat leave-one-cohort-out.
- Report compact-versus-full Spearman correlation within every cohort.

### Paired FFPE series

- GSE117606, using the fixed 51 adenoma–adjacent mucosa pairs.
- Report paired score difference, positive-pair fraction and two-sided Wilcoxon P value for full and compact scores.
- Report compact-versus-full correlation and gene-level paired directions for the eight frozen genes.

### Empirical perturbations

- Reuse the fixed APC-organoid and TCF7L2-knockout comparisons from the audited validation pipeline.
- Report compact and measurable full-programme score changes under the existing contrast directions.
- These analyses test responsiveness and concordance, not causal sufficiency of an eight-gene circuit.

## Interpretation boundary

- The compact panel is an RNA reconstruction readout if it preserves programme ordering and group effects across platforms.
- No AUC threshold, diagnostic cut point, prognostic claim, pharmacodynamic claim or clinical context of use is evaluated.
- External resources cannot be called independent if they were used previously to evaluate an older panel; for the new panel they are locked outcome-validation datasets, but not prospective clinical validation.
