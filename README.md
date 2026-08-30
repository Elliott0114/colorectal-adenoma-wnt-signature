# Colorectal adenoma epithelial identity remodelling

Reproducible code, frozen gene definitions, derived results and publication
figures for the Communications Biology manuscript:

> **Within-state epithelial identity remodelling predominates over compositional change in colorectal adenomas**

Release v3.0.0 follows the final manuscript structure and contains no retired
gene-set definitions or exploratory analyses that did not meet the
prespecified reporting gates.

## Main result

The study separates an adenoma–normal expression difference into altered
proportions of recoverable epithelial states and transcriptional change within
those states. Donor-aware models provide a continuous ranking of 8,221 jointly
testable genes. Fixed false-discovery, posterior-direction and local
false-sign criteria retain 1,843 high-confidence genes: 884 adenoma-up and 959
adenoma-down. In the donor-disjoint partition, 98.6% of testable programme
genes retain direction and discovery–validation effects correlate at
Spearman rho = 0.912. Exact donor-balanced decomposition assigns 79.1% of the
primary score difference to expression within fine states and 20.9% to changed
fine-state proportions.

## Functional architecture

Competitive gene-set tests use the complete signed 8,221-gene ranking in the
discovery and donor-disjoint partitions. Fifty-eight pathways satisfy the
fixed cross-partition and cross-state replication rule. Fifty-four are
enriched toward adenoma-down expression, showing broad loss of oxidative,
metabolic, transport and mature epithelial functions; four are enriched
toward adenoma-up expression. Leading-edge overlap reduces the 58 pathways to
28 nonredundant communities. Hallmark WNT signalling is directionally positive
in all eight contexts but does not pass the common-effect FDR thresholds, so it
is interpreted with gene anchors, regulatory, chromatin and perturbation
evidence rather than as a replicated pathway-level result.

## Candidate reduced readout

An eight-gene candidate was derived only after the full programme was frozen:

- adenoma-up: **EPHB2, REG1A, LTBP1, RNF43**
- adenoma-down: **CALM2, COX6C, B2M, ACAA2**

It is retained in the Supplementary Information as one tractable reduced
representation. External direction-balanced random-panel benchmarking does
not establish unique optimality, and it is not presented as a clinical
biomarker.

## Repository map

- `analysis/` — manuscript-facing R and Python scripts.
- `analysis/contracts/` — frozen implementation contracts.
- `workflow/pipeline.tsv` — ordered analysis map.
- `data/signatures/` — complete ranking, high-confidence programme,
  platform-measurable universe and eight-gene candidate.
- `data/datasets.tsv` — accession, primary-publication and analysis-role
  registry.
- `results/state_aware_program_v1/functional_architecture_v1/` — complete
  pathway replication results and leading-edge communities.
- `figures/communications_biology_v3.0/` — six main and seven supplementary
  figures in PDF, SVG, TIFF and PNG, with public panel source data.
- `tests/verify_release.py` — release-integrity and numerical checks.

## Reproduce or verify

The project is pinned to the `crc-premalignant-locked` environment.

```bash
conda env create -f environment.yml
make verify
make figures-public
```

Early analytical steps require source matrices downloaded from the public
repositories in `data/datasets.tsv`; these raw matrices are not redistributed.
The current figures and non-identifying derived tables are committed so that
the reported estimates can be audited without repeating every download.

## Data governance

Patient-level values from the self-generated DSLab series are not distributed.
Only aggregate, non-identifying statistics are included. Governed values and
additional code/data details are available from the corresponding author on
reasonable request and subject to institutional requirements.

## Citation

Citation metadata are in `CITATION.cff`. Cite the manuscript together with the
immutable v3.0.0 release used in its Code Availability statement.
