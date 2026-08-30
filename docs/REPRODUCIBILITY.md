# Reproducibility and claim boundaries

## Inferential units

Human analyses use the patient, donor or patient cluster as the independent
unit. Perturbation analyses use donors, organoid models, edited clones, cell
lines or oncogenic backgrounds. Spatial analyses use tissue sections. Cells,
nuclei, technical libraries and spots do not add inferential degrees of freedom.

## Information flow

```text
33,703 assayed discovery features
  -> three sampling-eligible epithelial states
  -> 8,221 genes testable in all three states
  -> continuous cross-state common-effect ranking
  -> 1,843 genes passing fixed confidence rules
  -> freeze programme membership, direction and scoring
  -> donor-disjoint and independent validation
  -> 53 label-independent portable protein-coding candidates
  -> eight-gene donor-held-out reconstruction
  -> compact-readout fidelity and random-panel benchmarks
```

No fixed target gene count is used for the 1,843-gene subset. The eight-gene
size is selected at the data-defined knee of the discovery donor-held-out
fidelity curve; validation outcomes are not used.

## Evidence vocabulary

| Evidence | Supported wording | Wording not supported |
|---|---|---|
| 8,221-gene ranking | continuous common-effect programme used for functional interpretation | exhaustive biological boundary |
| 1,843-gene subset | high-confidence state-aware subset | top-1,843 signature |
| Donor-disjoint Chen analysis | isolated source-study validation with donor overlap removed | external independent cohort |
| Fine-state decomposition | exact accounting of measured composition and expression within projected fine states | causal mediation, lineage conversion or unbiased in situ abundance |
| Five transcriptomic cohorts | independent cross-platform recurrence | prospective clinical validation |
| Eight-gene result | candidate reduced measurement of the full programme | uniquely optimal signature, diagnostic or validated biomarker |
| RNA-ATAC | regulatory association after patient-aware adjustment | direct transcription-factor binding |
| CRC Atlas | recurrence across cross-sectional epithelial contexts | longitudinal progression |
| Empirical genetic perturbations | signed responsiveness to APC-WNT regulation | universal perturbation effect size |
| GenKI virtual deletion | unsigned historical regulator context | experimental knockout or causal proof |
| Spatial and protein layers | localisation and reciprocal tissue anchors | protein implementation of the eight-gene RNA score |

## Fine-state analysis

Fine states are learned separately within absorptive, goblet and
transit-amplifying epithelium from normal discovery cells. All 1,843 programme
genes, lesion labels and prespecified technical or stress-related feature sets
are excluded before feature selection. The reference is frozen before every
validation cell is projected.

At the primary four-cluster resolution, exact donor-balanced decomposition
attributes 79.1% of the validation difference to expression within fine states
and 20.9% to changed fine-state proportions. Three- and five-cluster solutions
retain the same ordering. The identity reconstructs the observed total to
numerical tolerance, and uncertainty is estimated by resampling whole donors.

## Compact-readout benchmark

The compact candidate is evaluated against 2,000 direction-balanced portable
eight-gene panels. It modestly exceeds the internal 95th percentile but its
median external fidelity does not exceed the external 95th percentile. This
boundary is retained: the eight genes are tractable, not uniquely optimal.

## Virtual deletion

GenKI uses a donor-balanced held-out adenoma network, fixed targets, two model
seeds and expression-matched null genes. Distances are unsigned embedding
perturbation magnitudes. Direction is taken only from empirical perturbation
datasets. Virtual-deletion outputs are supplementary and do not participate in
programme or compact-readout definition.

## Data-governance boundary

Raw third-party data are not redistributed. Patient-level DSLab matrices,
metadata and plotted source values are also excluded from the public release;
they are available from the corresponding author subject to institutional
governance. Public figure source data and Supplementary Tables contain only
public-data outputs and aggregate DSLab statistics.

## Reproduction tiers

1. **Release integrity:** `make verify` checks frozen definitions, required
   outputs, figure inventories, sensitive-content patterns and file-size limits.
2. **Public figure reproduction:** `make figures` regenerates six main figures
   and Supplementary Figs. 5–8 from distributed derived tables. Supplementary
   Figs. 1–4 are distributed in final form; Fig. 3g requires governed DSLab
   patient-level source values.
3. **Full public-data reanalysis:** acquire the public raw inputs, create the
   locked environment and run the stages in `workflow/pipeline.tsv`.
4. **Virtual-deletion audit:** create the separate GenKI environment and run
   `make virtual-deletion-audit`.
