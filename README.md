# Colorectal adenoma epithelial remodelling

Reproducible code, frozen definitions, derived results and publication figures
for the Communications Biology manuscript:

> **Within-state epithelial remodelling in colorectal adenomas comprises separable regulatory and mature-function responses**

Release v4.0.0 follows the state-shared-remodelling manuscript structure. It
includes the donor-aware consensus-network and two-coordinate perturbation
analyses added after v3.0.0. Historical 287-gene and 12-gene development sets
are not included. Non-replicating screens, where retained for reproducibility,
remain audit-only and are not used as manuscript claims.

## Main result

The study separates the adenoma–normal epithelial transcriptomic difference
into altered proportions of recoverable normal-reference substates and
expression change within those substates. Donor-aware models provide a
continuous ranking of 8,221 jointly testable genes. Fixed false-discovery,
posterior-direction and local false-sign criteria retain 1,843 high-confidence
genes (884 adenoma-up and 959 adenoma-down); the count is an output of the
criteria rather than a prespecified programme size.

In the donor-disjoint partition, 98.6% of testable high-confidence genes retain
their common direction and discovery–validation effects correlate at Spearman
rho = 0.912. Exact donor-balanced decomposition assigns 79.1% of the primary
score difference to expression within normal-reference epithelial substates
and 20.9% to their changed proportions.

## Functional and perturbational architecture

Competitive tests use the complete signed 8,221-gene rank. Fifty-eight pathways
satisfy the fixed cross-partition replication rule; 54 map to adenoma-down
expression and form 28 leading-edge communities dominated by oxidative,
metabolic and mature epithelial functions. Hallmark WNT signalling is
directionally positive in all state-level contexts but does not pass the
common-effect FDR thresholds.

An exploratory consensus WGCNA uses all 8,221 genes after removal of the
adenoma–normal main effect from independent donor profiles. Seven modules are
direction-consistent and preserved in all three held-out state networks. They
organise the continuous disease-effect rank into a focused adenoma-up group and
multiple adenoma-down mature-function modules; they do not replace the rank or
define portable clinical scores.

Prespecified genetic and pharmacological perturbations are represented by two
response coordinates: WNT/stem suppression and mature-function restoration.
These coordinates can separate. For example, ASCL2 knockout favours the
WNT/stem coordinate while moving the mature-function coordinate adversely,
whereas trametinib shows the reciprocal partial response.

## Candidate reduced readout

An eight-gene candidate is retained only in the Supplementary Information as
one tractable approximation of a transcriptionally redundant response:

- adenoma-up: **EPHB2, REG1A, LTBP1, RNF43**
- adenoma-down: **CALM2, COX6C, B2M, ACAA2**

External direction-balanced random-panel benchmarking does not establish its
unique optimality, and it is not presented as a diagnostic, prognostic or
clinically deployable biomarker.

## Repository map

- `analysis/` — manuscript-facing R and Python scripts.
- `analysis/contracts/` — frozen implementation contracts.
- `workflow/pipeline.tsv` — ordered analysis map.
- `data/signatures/` — complete rank, high-confidence genes,
  platform-measurable universe and candidate reduced readout.
- `data/datasets.tsv` — accession, primary-publication and analysis-role
  registry.
- `results/state_aware_program_v1/functional_architecture_v1/` — pathway
  replication and leading-edge communities.
- `results/state_aware_program_v1/functional_architecture_exploratory_v2_1/`
  — consensus modules, preservation, cross-context projection, perturbation
  responses and sentinel routing.
- `results/state_aware_program_v1/identity_reversal_target_prioritization_v1/`
  — PROGENy, CollecTRI and two-coordinate perturbation outputs.
- `figures/communications_biology_v5.0/` — six main and eight supplementary
  figures with panel-level public source data.
- `tests/verify_release.py` — release-integrity and numerical checks.

## Reproduce or verify

The project is pinned to the `crc-premalignant-locked` environment.

```bash
conda env create -f environment.yml
make verify
make figures-public
```

Early analytical steps require source matrices downloaded from the repositories
listed in `data/datasets.tsv`; raw third-party matrices are not redistributed.
Committed non-identifying results allow the manuscript estimates and figures to
be audited without repeating every download.

## Data governance

Patient-level values from the self-generated DSLab series are not distributed.
Only aggregate, non-identifying statistics are included. Governed values are
available from the corresponding author on reasonable request and subject to
institutional requirements.

## Citation

Citation metadata are in `CITATION.cff`. Cite the manuscript together with the
immutable v4.0.0 release used in its Code Availability statement.
