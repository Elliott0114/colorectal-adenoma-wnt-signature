# Colorectal adenoma epithelial identity remodelling

This repository contains analysis code, frozen gene definitions, derived
results and figure source data supporting the manuscript:

> Within-state epithelial identity remodelling predominates over compositional
> change in colorectal adenomas

The main analysis integrates public, de-identified datasets. A self-generated
de-identified colorectal-polyp single-cell series is used only for a
copy-number-composition sensitivity analysis; its patient-level data are not
distributed in this repository.

## Study question

A lesion-level expression difference can reflect changed proportions of
epithelial states, changed expression within comparable states, or both. The
analysis therefore asks whether conventional colorectal adenomas undergo a
compartment-wide change in epithelial identity after measured cell-state
composition is accounted for.

## Analysis architecture

1. Specimen-by-state pseudobulk profiles were built for absorptive, goblet and
   transit-amplifying epithelium. Donors, rather than cells, were the biological
   replicates.
2. State-specific mixed-model effects were integrated into a continuous
   8,221-gene common-effect ranking. Fixed false-discovery, posterior-direction
   and local false-sign criteria yielded 1,843 high-confidence genes: 884
   adenoma-up and 959 adenoma-down. No target gene count was supplied.
3. The deposited validation partition was made donor-disjoint. Of the 1,646
   high-confidence genes testable in all three validation states, 98.6%
   retained the common direction; discovery and validation effects correlated
   at rho = 0.912.
4. Fine states were learned from normal discovery cells after exclusion of all
   programme genes and lesion labels, then frozen before validation projection.
   Exact donor-balanced decomposition assigned 79.1% of the primary validation
   difference to expression within fine states and 20.9% to altered fine-state
   proportions.
5. Independent transcriptomic cohorts, paired FFPE tissue, matched RNA-ATAC,
   epithelial atlases, spatial tissue sections and empirical genetic
   perturbations tested transfer and biological context.
6. Only after the full programme was frozen, platform availability and
   donor-held-out reconstruction reduced it to a balanced eight-gene candidate
   measurement. Random-panel benchmarking is retained to show that the compact
   implementation is tractable but not uniquely optimal.

The 8,221-gene ranking and confidence-defined 1,843-gene subset describe the
biological finding. The eight genes provide one reduced measurement; they are
not presented as a diagnostic, a clinically validated biomarker or a closed
regulatory circuit. The historical 287- and 12-gene sets are retained only for
transparent supplementary comparison.

## Frozen definitions

The full programme score uses equal-arm averaging:

```text
mean(z[884 adenoma-up genes]) - mean(z[959 adenoma-down genes])
```

The compact candidate score is:

```text
mean(z[EPHB2, REG1A, LTBP1, RNF43])
  - mean(z[CALM2, COX6C, B2M, ACAA2])
```

Authoritative files are in [`data/signatures`](data/signatures). The compact
set was fixed after discovery-only donor-held-out size selection; validation
outcomes were not used for membership or panel size.

## Repository layout

```text
analysis/               analysis and R figure-generation scripts
analysis/contracts/     frozen analysis contracts
data/datasets.tsv       accession and primary-publication registry
data/signatures/        8,221-, 1,843-, 53- and 8-gene definitions
data/source_data/       panel-level source data for public analyses
data/supplementary/     Supplementary Tables 1-12 workbook
figures/communications_biology_v2.0/
                         six main and eight supplementary figures
results/state_aware_program_v1/
                         programme derivation and validation outputs
results/state_shared_revision_v2/
                         donor-disjoint, fine-state and benchmark outputs
tests/                   release-integrity checks
workflow/                ordered analysis map
```

Versioned script names are retained to preserve the audit trail. Names that
refer to an earlier manuscript target do not alter the scientific content.

## Quick verification

From the repository root:

```bash
conda env create -f environment.yml
conda run -n crc-premalignant-locked python tests/verify_release.py
```

To regenerate the six main figures and Supplementary Figs. 5–8 from the
distributed derived tables:

```bash
make figures
```

The R renderers export PDF and SVG vector files, 600-dpi TIFF files and 300-dpi
PNG files. Supplementary Figs. 1–4 are distributed in final form; their renderer
is retained for provenance, but Fig. 3g requires the governed DSLab patient-level
source. Full reanalysis requires the public raw data listed in
[`data/datasets.tsv`](data/datasets.tsv); these third-party inputs are not
redistributed.

## Self-generated data boundary

Patient-level DSLab matrices, metadata and source values are available from the
corresponding author on reasonable request and subject to institutional
data-governance requirements. Only aggregate statistics are distributed here.
The source values for the patient-level panel in Supplementary Fig. 3g are
therefore omitted from the public Source Data workbook.

## Virtual deletion

GenKI virtual deletion is supplementary, unsigned regulator-context evidence.
It tests whether fixed target deletions preferentially affect historical gene
sets relative to matched genes. Biological direction is supplied by empirical
perturbation data. The analysis uses a separate environment because its PyTorch
requirements differ from the main workflow:

```bash
conda env create -f environment-virtual-knockout.yml
conda run -n crc-premalignant-virtual-ko python \
  analysis/audit_genki_reproducibility_v2_9.py
```

## Reproducibility

- Patients, donors, patient clusters, models, clones or tissue sections are the
  inferential units; cells, nuclei and spots are nested measurements.
- Fine-state uncertainty uses 5,000 whole-donor bootstrap samples.
- The compact candidate is compared with 2,000 direction-balanced portable
  panels in donor-disjoint and external data.
- Random seeds, input hashes, environment versions and quality gates are
  recorded in result manifests.
- Figure source data are distributed independently of raw third-party data.

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for the supported
interpretations and analysis boundaries.

## Licence and citation

Code is released under the MIT License. Derived tables are provided for
research reproducibility; source-dataset licences and terms remain in force.
Cite the associated manuscript, this software release and the primary
publication for every reused dataset. Citation metadata are provided in
[`CITATION.cff`](CITATION.cff).
