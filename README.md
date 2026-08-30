# Colorectal adenoma epithelial identity remodelling

Reproducible code, frozen gene definitions, derived results and publication
figures for the Communications Biology manuscript:

> **Within-state epithelial identity remodelling predominates over compositional change in colorectal adenomas**

Release v2.1.0 reorganises the public workflow around the final biological
object: a state-shared epithelial programme derived in absorptive, goblet and
transit-amplifying cells.

## Main result

The analysis separates two sources of an adenoma–normal expression difference:

1. altered proportions of recoverable epithelial states; and
2. transcriptional change within those states.

The continuous discovery analysis contains 8,221 jointly testable genes.
Fixed false-discovery, posterior-direction and local false-sign criteria retain
1,843 high-confidence genes (884 adenoma-up and 959 adenoma-down). In the
donor-disjoint partition, 98.6% of testable programme genes retain the common
direction and discovery–validation effects correlate at Spearman ρ = 0.912.
At the primary fine-state resolution, exact donor-balanced decomposition
attributes 79.1% of the score difference to expression within fine states and
20.9% to changed fine-state proportions.

The programme combines focused WNT/stem regulation with broad loss of mature
metabolic, transport and absorptive functions. Independent cohorts, paired
FFPE tissue, epithelial atlases, matched RNA–ATAC, spatial sections and
genetic perturbations test this frozen full programme.

## Candidate reduced readout

An eight-gene candidate was derived only after the 1,843-gene programme had
been frozen:

- adenoma-up: **EPHB2, REG1A, LTBP1, RNF43**
- adenoma-down: **CALM2, COX6C, B2M, ACAA2**

It is one tractable reduced representation of a transcriptionally redundant
programme. Internal donor-held-out fidelity is high, whereas external
direction-balanced random-panel benchmarking does not identify this
combination as uniquely optimal. It is not presented as a clinically
validated biomarker.

## Repository map

- `analysis/` — ordered R and Python analysis scripts.
- `analysis/contracts/` — frozen analysis contracts, including the
  outcome-blind full-programme projection rule.
- `workflow/pipeline.tsv` — manuscript-facing execution order.
- `data/signatures/` — frozen ranking, high-confidence programme,
  platform-measurable candidate universe and eight-gene candidate.
- `data/datasets.tsv` — accession, primary-publication and analysis-role
  registry.
- `results/state_aware_program_v1/` — derived state-aware results.
- `results/state_shared_revision_v2/` — donor-disjoint, fine-state,
  meta-analysis and reduced-readout benchmarks.
- `figures/communications_biology_v2.1/` — six main and seven
  supplementary figures in PDF, SVG, TIFF and PNG, with panel source data.
- `tests/verify_release.py` — release-integrity and numerical checks.

## Reproduce the public analyses

The project is pinned to the `crc-premalignant-locked` environment.

```bash
conda env create -f environment.yml
make verify
make figures
```

The complete ordered workflow is in `workflow/pipeline.tsv`. Several early
steps require public source matrices that are not redistributed; the frozen
derived inputs used for manuscript figures and audits are included where
licensing permits.

## Full-programme downstream projection

The full programme is projected into the Becker, CRC Atlas, perturbation and
spatial layers only when both directional arms have at least 75% feature
coverage and neither arm has fewer than 100 measurable genes. Eligibility is
determined from platform feature inventories before phenotype labels or
effect estimates are inspected.

Relevant entry points:

```bash
conda run -n crc-premalignant-locked python \
  analysis/audit_state_shared_full_program_coverage_v1.py

conda run -n crc-premalignant-locked python \
  analysis/validate_state_shared_full_programme_extended_layers_v2.py
```

## Data governance

Public accession-level inputs are listed in `data/datasets.tsv`.
Patient-level values from the self-generated DSLab series are not distributed.
Only aggregate, non-identifying statistics are included; governed values are
available from the corresponding author on reasonable request and subject to
institutional requirements.

## Citation

Citation metadata are provided in `CITATION.cff`. The immutable release tag
and commit should be cited with the associated manuscript.
