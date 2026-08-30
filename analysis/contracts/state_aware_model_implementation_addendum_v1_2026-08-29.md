# Implementation addendum: state-aware discovery models v1

**Date:** 29 August 2026
**Status:** FROZEN before the full-genome differential-state models were run
**Parent contract SHA-256:** `0f39a03154408d14bb0bbe0c1e4f55d498e94d403bd6b45d4be5dd8918555f69`

## Reason for this addendum

Two alphabetically selected, non-inferential smoke tests were used only to verify software compatibility and coefficient extraction. No held-out data were opened and no gene, state, threshold or scientific hypothesis was changed.

The first smoke test identified incompatibility between `lme4 2.0-6` and the installed `lmerTest`/`variancePartition` interface. The environment was therefore pinned to `lme4 1.1-38` and `Matrix 1.7-4`. The second test established that `MArrayLM2` objects must be moderated with the `variancePartition::eBayes` S4 method rather than `limma::eBayes`.

## Frozen implementation clarifications

1. Per-state models use `variancePartition::voomWithDreamWeights` followed by `variancePartition::dream` with Satterthwaite degrees of freedom.
2. The model remains `expression ~ route + (1 | donor_id)`.
3. Empirical-Bayes moderation uses `variancePartition::eBayes(..., robust = TRUE)`.
4. Model-based unmoderated standard errors are retained for multivariate adaptive shrinkage; moderated statistics are used for gene-level reporting.
5. Four SOCK workers are used. Parallelisation changes runtime only.
6. The shared gene universe is the intersection of genes passing default `edgeR::filterByExpr` in ABS, GOB and TAC.
7. Paired-donor sensitivity sums multiple specimens within donor-route-state and fits `~ donor_id + route` using voom/limma without redefining the gene universe.
