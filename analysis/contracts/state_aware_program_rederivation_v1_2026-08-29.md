# Analysis contract: state-aware redefinition of the adenoma epithelial programme (v1)

**Date:** 29 August 2026
**Status:** FROZEN before any new gene-level state-aware outcome was calculated; execution authorised by the user on 29 August 2026
**Nature of analysis:** post hoc biological redefinition with a prospectively fixed gene-level implementation
**Primary purpose:** determine whether the adenoma-associated epithelial signal represents coordinated expression remodeling within comparable epithelial states, rather than only a change in the abundance of lesion-labelled epithelial cells.

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: plan
- Origin Date: 2026-08-29
- Verification Status: UNVERIFIED
- Version Label: state_aware_program_plan_v1

## 1. Known information and separation from new outcomes

The following information was already known when this contract was written:

- the original 287-gene set was derived from donor-aware conventional-adenoma versus normal differential expression in the Chen discovery partition;
- score-level decomposition in the Chen held-out partition showed that the frozen 287-gene score differed within absorptive, goblet and transit-amplifying states;
- the 12-gene score was broadly concordant with the 287-gene score;
- current self-generated DSLab data support composition independence but do not provide reliable transferred cell-state labels.

The following outcomes have **not** yet been calculated and must remain unopened until this contract is frozen:

- raw-count, gene-level differential-state effects in individual canonical epithelial states;
- the cross-state common-effect ranking and any state-shared gene set;
- pathway and regulator enrichment of the new common-effect ranking;
- overlap and enrichment of the current 287 genes in that ranking;
- fidelity of the existing 12-gene panel to the new state-shared programme;
- corresponding held-out gene-level validation results.

No held-out result may be used to change discovery-stage eligibility, filtering, model, gene membership, gene direction, weights or compact-panel size.

## 2. Research questions and hypotheses

### Primary question

Does conventional adenoma induce a coordinated transcriptional displacement within multiple comparable epithelial states?

### Primary hypothesis

A common adenoma-associated effect is present across adequately sampled absorptive, goblet and transit-amplifying epithelial states, with concordant direction at gene and programme levels.

### Secondary hypotheses

1. The original 287-gene set is enriched for, but need not be identical to, the cross-state common-effect programme.
2. Biological interpretation of the common-effect ranking resolves into coordinated regulatory and functional modules rather than an unstructured differential-expression list.
3. A compact RNA panel can reconstruct the common-effect programme without fixing its gene count in advance.
4. Frozen multi-omic, spatial, perturbational and self-generated datasets support different parts of this programme without participating in gene selection.

## 3. Analysis hierarchy and immutable inputs

### 3.1 Discovery only

`data_sources/Chen_Cell_2021_CELLxGENE/chen_discovery_epithelial.h5ad`
SHA-256: `2e2b5aa4909fd314ea0c6dba17f6ced16de95c1bb3317e69a3edfe2caff17440`

The discovery partition alone may be used for gene-level modelling, programme definition, biological module discovery and compact-panel construction.

### 3.2 First held-out validation

`data_sources/Chen_Cell_2021_CELLxGENE/chen_validation_epithelial.h5ad`
SHA-256: `f51bf02885d637073743093844c1db00a173196f0e1184b37f39a5118468668d`

This partition may be opened only after discovery outputs and compact-panel rules are frozen.

### 3.3 Existing reference signatures

- `release/colorectal-adenoma-wnt-signature/data/signatures/core_287_genes.tsv`
  SHA-256: `c649ea65e6f32c667d347872bd16ec4712ef12202e359d680052b44def9cddf2`
- `release/colorectal-adenoma-wnt-signature/data/signatures/signature_12_genes.tsv`
  SHA-256: `91e564f42bdc5fb0188605b76116bcfb126e5102d62a33d6890fa67018928380`

These sets are comparators, not priors for the new gene-level model.

## 4. Population, contrast and state eligibility

### 4.1 Primary contrast

- Normal: `NL`
- Conventional adenoma: `TA`, `TV`, `TVA`
- Serrated lesions (`HP`, `SSL`) are excluded from derivation and reserved for specificity analysis.
- `UNC` and unrecognised labels are excluded.

### 4.2 State eligibility is based on sampling only

Disease-labelled states `ASC` (neoplastic) and `SSC` (abnormal) are excluded from the primary derivation. A deposited canonical state is eligible if, among specimen-state strata containing at least 20 cells, it has:

1. at least 10 donors in each route;
2. at least 6 donors represented in both routes; and
3. non-zero expression information sufficient for pseudobulk normalisation.

Applied without reference to gene-expression outcomes, these rules retain:

- `ABS` — absorptive;
- `GOB` — goblet;
- `TAC` — transit-amplifying.

Crypt stem, other colonic epithelium, tuft and enteroendocrine states are descriptive sensitivities only because of sparse donor support. No state may be added or removed after gene-level outcomes are inspected.

## 5. Raw-count pseudobulk construction

1. Use integer counts from `raw/X`, never log-normalised `X`, for differential-state modelling.
2. Preserve deposited gene symbols and state labels; record duplicate-symbol resolution before modelling.
3. Sum counts within `specimen_id × state` strata.
4. Retain strata containing at least 20 epithelial cells.
5. Retain multiple specimens from one donor; donor is the biological unit and specimen is a repeated measurement.
6. Record cells, total counts, detected genes, specimens and donors for every route-state stratum.
7. Exclude genes using expression-only filtering (`edgeR::filterByExpr`) within each eligible state. The primary shared universe is the intersection of genes passing filtering in all three states; the union is a sensitivity analysis.
8. Use TMM library normalisation and precision weights. Cell count and library size affect weights, not the biological sample size.

Numerical and metadata audits must pass before any differential model is fitted.

## 6. Donor-aware differential-state models

### 6.1 State-specific primary models

For each eligible state and gene, fit a precision-weighted mixed model:

`expression ~ route + (1 | donor_id)`

The route coefficient is conventional adenoma minus normal. Multiple specimens per donor remain in the model. Extract log-fold change, standard error, moderated test statistic and nominal P value.

### 6.2 Paired-donor sensitivity

Repeat each state-specific analysis using only donors represented in both routes, with donor blocking. This analysis tests whether the direction depends on between-donor imbalance.

### 6.3 Other fixed sensitivities

- minimum specimen-state sizes of 10 and 50 cells;
- exclusion of curated cell-cycle genes from downstream programme scoring;
- conventional adenoma versus serrated lesion, used only after the primary contrast is frozen;
- all adequately observed canonical states as a descriptive extension, with no effect on membership.

No cell-level P values and no route-dependent re-clustering are permitted.

## 7. Cross-state integration without a fixed gene count

### 7.1 Continuous common-effect ranking is the primary biological object

Combine the three state-specific effect estimates using their estimated null correlation. A correlation-adjusted generalized Stouffer statistic provides one signed, genome-wide common-effect ranking. Inconsistent state effects cancel rather than being hidden by pooling.

### 7.2 Multivariate shrinkage

Use multivariate adaptive shrinkage (`mashr`) with canonical and data-driven covariance components to estimate:

- posterior effects in each state;
- local false sign rates (LFSR);
- the degree and pattern of effect sharing across states.

The strict state-shared core is defined without a target size:

1. BH-adjusted common-effect P value ≤ 0.05;
2. identical posterior effect direction in ABS, GOB and TAC; and
3. LFSR ≤ 0.05 in all three states.

A relaxed sensitivity uses LFSR ≤ 0.10. The continuous ranking, not the number of genes passing either threshold, drives pathway analysis and the main biological conclusions.

### 7.3 Stability audit

Repeat effect estimation after leaving out each donor in turn. Report gene-level sign stability and rank stability. Stability is descriptive and cannot be used to tune the main thresholds after inspection.

## 8. Relationship to the original 287-gene set

The current 287 genes are not automatically retained or discarded. They undergo four fixed audits:

1. separate competitive enrichment tests for the 89-gene up arm and 198-gene down arm against the complete common-effect ranking;
2. state-specific and common-effect sign concordance;
3. overlap with the strict and relaxed state-shared cores;
4. comparison with 10,000 expression- and detectability-matched random gene sets of identical arm sizes.

Interpretation follows the evidence:

- **Outcome A — strong cross-state support:** retain the 287-gene set as a donor-stable lesion core that is substantially supported across canonical states; do not call 287 a natural biological boundary.
- **Outcome B — partial cross-state support:** retain 287 only as a broad lesion-associated core and promote the objectively defined state-shared subset/ranking as the biological programme.
- **Outcome C — weak cross-state support:** replace the 287-gene main claim and report the failure transparently; the prior score-level finding cannot rescue gene-level inconsistency.

No percentage overlap is designated as a success threshold. Enrichment effect sizes, confidence intervals, matched-null probabilities and held-out replication jointly determine the interpretation.

## 9. Biological knowledge integration

Biological knowledge is applied **after**, not before, unbiased common-effect estimation.

1. Run competitive rank-based enrichment on the complete common-effect statistic using:
   - Hallmark v2026.1 (already frozen locally);
   - version-frozen Reactome pathways;
   - a trimmed Gene Ontology Biological Process collection.
2. Correct within each collection using BH FDR.
3. Extract leading-edge genes for supported pathways.
4. Build a pathway-overlap graph from leading-edge Jaccard similarity and group redundant terms into data-driven modules.
5. Assign concise biological labels only after reviewing the member pathways and genes.
6. Estimate pathway and transcription-factor activities with fixed PROGENy and CollecTRI priors as orthogonal interpretation. These activities cannot alter gene membership.

Expected hypotheses include WNT/stem-regulatory activation and loss of mature absorptive/metabolic functions, but all curated pathways remain eligible and discordant findings must be retained.

## 10. Compact-panel audit and possible rederivation

### 10.1 Test the existing 12-gene panel first

The existing panel is evaluated against the new state-shared programme using:

- arm-direction agreement for every gene;
- donor-grouped leave-one-donor-out correlation with the full common-programme score;
- fidelity within ABS, GOB and TAC separately;
- held-out correlation after the panel is frozen;
- preservation of both activation and loss-of-function biological modules.

### 10.2 Objective rederivation if the existing panel fails

If any gene has the wrong shared-effect direction or the score does not preserve the programme across states, reconstruct the panel from the strict state-shared core intersected with a label-blind, protein-coding, cross-platform availability set.

1. Add one up/down gene pair at each step using discovery-only donor-grouped out-of-fold reconstruction of the full common-programme score.
2. Choose panel size at the reproducible knee of the out-of-fold fidelity curve; do not specify 10, 12 or any other number in advance.
3. Use 2,000 whole-donor bootstraps to report selection stability.
4. Use biological-module coverage only as a prespecified tie-breaker between panels with indistinguishable out-of-fold fidelity.
5. Freeze genes, arms and scoring rule before opening the held-out partition.

The output is called a compact RNA readout, not a diagnostic or pharmacodynamic biomarker.

## 11. Held-out and orthogonal validation sequence

After discovery freeze, validation proceeds in this order:

1. **Chen held-out epithelial partition:** reproduce state-specific effects, common-ranking enrichment and compact/full-score concordance without tuning.
2. **Independent adenoma transcriptomic cohorts and paired FFPE:** test only the frozen programme/readout at patient or specimen level.
3. **Becker RNA–ATAC:** assess whether programme direction is supported by concordant accessibility, without selecting genes from ATAC results.
4. **APC/WNT/ASCL2/TCF7L2 perturbations and virtual perturbation:** test predicted bidirectional movement of the frozen score; no new exploratory target search.
5. **CRC atlas and spatial datasets:** assess cross-study state support and in-tissue localisation.
6. **Public protein evidence:** use reciprocal tissue anchors as protein-level triangulation, not as proof that the RNA panel is an assay.
7. **DSLab self-generated data:** retain only composition-standardised and CNV-composition-independence analyses unless a reliable state-mapping method is independently established.

External datasets may falsify or support the frozen programme but cannot redefine it.

## 12. Primary success, mixed and failure interpretations

### Supports the revised main line

- the discovery common-effect ranking is non-null and biologically coherent;
- the frozen programme is directionally enriched in the held-out state-aware analysis;
- held-out ABS, GOB and TAC programme effects are concordant;
- the 287 audit or its objective subset shows enrichment beyond matched random sets;
- a frozen compact readout preserves the full programme without outcome-guided tuning.

### Mixed result

The average programme replicates but gene sharing is incomplete or restricted to two states. The manuscript must then describe state-predominant remodeling rather than a universal epithelial identity shift.

### Failure

Held-out state-aware effects are inconsistent, biological modules do not replicate, or apparent support is driven by disease-labelled states. The current 287-gene programme must not be presented as a distributed epithelial remodeling programme.

## 13. Planned figures and manuscript consequences

### Revised Figure 1

- **a:** study resources and discovery/validation separation;
- **b:** specimen-by-state raw-count pseudobulk design and objective state eligibility;
- **c:** ABS/GOB/TAC state-specific effects and cross-state sharing map;
- **d:** common-effect pathway/module network;
- **e:** relationship between the common programme and the original 287 genes;
- **f:** full programme versus objectively sized compact RNA readout.

The old multi-threshold gene-count waterfall moves to Supplementary Methods/Figure S1.

### Revised Figure 2

Lead with held-out within-state replication and composition-versus-expression decomposition. This figure establishes that the signal is not reducible to an increased fraction of lesion-labelled epithelial cells.

### Later main figures

Retain external cohort/FFPE replication, RNA–ATAC/regulatory support, perturbation/virtual perturbation, and spatial/protein/self-data triangulation. The virtual perturbation result should remain in the main text if it directly tests the frozen programme rather than opening a new target-discovery branch.

## 14. Reproducibility and implementation

- Statistical analysis and all manuscript figures will be implemented in R.
- Create a dedicated, version-locked state-aware environment; do not alter `crc-premalignant-locked`.
- Required packages include `edgeR`, `dreamlet`/`variancePartition`, `mashr`, `limma`, `decoupleR`, `msigdbr`, `igraph`, `ggplot2` and `patchwork`.
- Record package versions, random seeds, input hashes, output hashes and session information.
- Every figure must have machine-readable source data and a manifest linking panel to script and input.
- All new outputs will be written to `results/state_aware_program_v1/`; no existing result is overwritten.

## 15. Execution stages and gates

1. **Environment and input audit** → gate: raw counts, metadata, hashes and state eligibility reproduce exactly.
2. **Discovery pseudobulk and state-specific models** → gate: model convergence, valid residual diagnostics and no donor leakage.
3. **Cross-state integration and biological modules** → gate: contract-compliant outputs generated without held-out access.
4. **287 audit and compact-panel decision** → gate: discovery membership and scoring rule frozen.
5. **Held-out validation** → gate: all planned outcomes reported, including negative findings.
6. **Orthogonal evidence re-scoring** → gate: no dataset has contributed to selection.
7. **Figure/manuscript rebuild** → gate: every quantitative statement maps to source data.
8. **Independent reproducibility audit** → gate: clean rerun reproduces tables, figures and manifests.
