# Implementation addendum: orthogonal validation of the frozen state-shared compact readout (v1)

**Date:** 29 August 2026
**Status:** FROZEN before opening Becker RNA–ATAC, CRC Atlas, additional perturbation, or spatial outcomes for the new eight-gene readout
**Parent contract:** `state_aware_program_rederivation_v1_2026-08-29.md`

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: execution
- Origin Date: 2026-08-29
- Verification Status: UNVERIFIED
- Version Label: state_shared_extended_validation_v1

## Frozen input and boundary

- The compact readout is fixed at eight genes before this analysis: four adenoma-up genes and four adenoma-down genes.
- Membership, direction, equal-arm scoring, and panel size cannot be changed using any result below.
- The full 1,843-gene programme has already undergone held-out, bulk-cohort, and paired-FFPE validation. Orthogonal layers below evaluate portability and biological concordance of the compact reconstruction; they do not redefine the full programme.
- No diagnostic cut point, prognostic endpoint, clinical utility, or pharmacodynamic claim is evaluated.

## Locked validation layers

1. **Becker snRNA and matched RNA–ATAC:** score the fixed readout in all nuclei and epithelial-marker-positive nuclei; test polyp versus unaffected normal tissue; correlate the RNA score with the pre-existing WNT/TCF/ASCL2 accessibility axis and repeat the existing patient-level robustness analyses.
2. **CRC Atlas:** calculate donor-level scores without cell or gene reselection; test polyp/cancer epithelial groups against normal epithelium and repeat the existing leave-one-study-out influence analysis.
3. **Additional perturbations:** reuse the fixed contrasts in GSE114059, GSE67186, GSE130822, and GSE171910. Human–mouse conversion is restricted to the existing one-to-one homology map. Report every non-null contrast, including discordant directions.
4. **Spatial sections:** score the six fixed sections and compare pathologist-annotated tumour regions with non-neoplastic epithelium. Report raw and proliferation/epithelial-control-adjusted score differences.

## Scoring and missingness

- Within each dataset, each measurable gene is standardised without using outcome labels.
- Score = mean z(four up genes) minus mean z(four down genes).
- Human datasets require all eight genes unless the source platform makes a gene structurally unavailable; any failure is reported without substitution.
- Mouse datasets report one-to-one homologue coverage. No paralogue substitution is permitted.
- Existing comparison units, metadata mappings, controls, null construction, and statistical tests are reused unchanged.

## Interpretation

- Concordance across these layers supports the eight-gene set as a compact RNA reconstruction of the state-shared epithelial programme.
- Perturbation movement supports responsiveness of the readout, not causal sufficiency of the eight genes.
- Spatial enrichment supports tissue localisation, not diagnostic performance.
