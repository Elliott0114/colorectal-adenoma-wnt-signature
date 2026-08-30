# Implementation note: efficient whole-donor leave-out stability audit (v1)

**Date:** 29 August 2026
**Status:** FROZEN before any donor leave-out result was generated
**Parent contract:** `state_aware_program_rederivation_v1_2026-08-29.md`

The initial implementation attempted to repeat the full per-gene `dream` mixed model for every donor-by-state fold. It was interrupted after more than six minutes without completing or writing a single fold result; at the observed rate, the audit would require several hours. No output was inspected and no gene, threshold or expected direction was changed.

For computational tractability, the descriptive stability audit uses the following fixed donor-aware approximation:

1. retain the frozen 1,843 strict state-shared genes and the complete 8,221-gene discovery normalisation universe;
2. leave out the complete donor across all specimens and states;
3. recalculate TMM factors and voom precision weights within each fold;
4. fit `limma::lmFit` with donor as the repeated-measures block and the consensus within-donor correlation estimated by `limma::duplicateCorrelation`;
5. integrate ABS, GOB and TAC effects by GLS using the frozen full-discovery null-correlation matrix;
6. report fold-to-full rank correlation and sign stability only.

The primary effect estimates, programme definition, significance tests and biological conclusions continue to use the prespecified `dream` mixed models. This approximation is a robustness audit and cannot add, remove or reweight a gene.
