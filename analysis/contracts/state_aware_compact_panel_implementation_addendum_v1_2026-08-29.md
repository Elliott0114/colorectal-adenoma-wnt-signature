# Implementation addendum: compact RNA reconstruction of the state-shared programme (v1)

**Date:** 29 August 2026
**Status:** FROZEN before any new compact-panel membership, path or fidelity result was calculated
**Parent contract:** `state_aware_program_rederivation_v1_2026-08-29.md`
**Scope:** implementation details for Section 10 of the parent contract; this addendum does not modify the state-shared programme.

## Outcomes already known at freeze

- The strict state-shared core contains 1,843 genes (884 up and 959 down).
- The existing 12-gene panel has nine genes in the shared modelling universe; all nine have the expected direction and meet the strict shared rule.
- SLC7A5, SDC3 and DPEP1 do not pass the three-state shared expression universe and therefore trigger objective rederivation under Section 10.2 of the parent contract.
- No compact-panel path, selected size, gene membership, held-out fidelity or external validation result has been calculated.

## Fixed input and eligibility

1. The target is the complete strict state-shared programme, not the original 287-gene score.
2. Candidate genes must be members of the strict core, NCBI protein-coding genes and present as features on all five prespecified external transcriptomic platforms plus GSE117606 FFPE.
3. Platform eligibility uses feature identifiers only. Tissue labels, group effects and validation scores are prohibited.
4. Candidate direction is the frozen shared direction; directions cannot be swapped during reconstruction.

## Discovery profiles and target score

1. Use the 134 raw-count specimen-by-state pseudobulks already constructed for ABS, GOB and TAC.
2. Apply TMM normalisation and log-CPM transformation.
3. Within every donor-held-out fold, centre and scale each gene separately within each epithelial state using training profiles only; apply those training parameters to the held-out donor.
4. The full-programme target is the equal-arm score: mean z expression of strict up genes minus mean z expression of strict down genes. Equal-arm averaging prevents the 959-gene down arm from dominating by gene count.
5. Leave out the complete donor, including every specimen and state, in each outer fold.

## Balanced reconstruction path

1. At each step add exactly one frozen-up and one frozen-down candidate gene.
2. Select the pair that maximises Pearson correlation with the full-programme target in the training profiles; ties are broken alphabetically.
3. The implementation may use an algebraically equivalent vectorised covariance calculation, but must reproduce exhaustive pair evaluation.
4. The path is evaluated up to the smaller of the two candidate-arm sizes and `number of training donors minus one`. This prevents a reconstruction path longer than the independent biological training units.
5. For each step, calculate a single donor-grouped out-of-fold score and its Spearman correlation with the out-of-fold full-programme target.

## Size and stability rules

1. The primary panel size is the Kneedle maximum distance from the diagonal of the monotone envelope of the donor-grouped out-of-fold fidelity curve.
2. Report the one-standard-error size as a sensitivity, using 2,000 whole-donor bootstrap replicates of the out-of-fold fidelity curve.
3. Refit the balanced path on all discovery profiles only after the primary size is fixed.
4. Report gene selection frequencies from the donor-held-out paths and from 2,000 whole-donor bootstrap refits at the fixed number of pairs.
5. Bootstrap resampling preserves all profiles belonging to each sampled donor. Duplicate sampled donors are treated as distinct bootstrap clusters.
6. Biological-module coverage may break an exact numerical tie between otherwise identical training correlations, but may not override a higher reconstruction correlation. No exact tie is anticipated.

## Freeze and validation boundary

- Freeze genes, directions, equal-arm scoring and state-specific standardisation parameters before opening Chen held-out data.
- The held-out and external datasets may evaluate or falsify the compact readout but may not replace genes, change signs, alter weights or change its size.
- The output is termed a **compact RNA readout of the state-shared programme**, not a diagnostic, prognostic or pharmacodynamic biomarker.
