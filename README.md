# Colorectal adenoma WNT-associated epithelial state

This repository contains the analysis code, frozen gene definitions, derived
results and figure source data supporting the manuscript:

> Multicohort analysis defines a reproducible epithelial state associated with
> WNT in colorectal adenomas

The study uses previously published, de-identified datasets. No new patient or
experimental data were generated.

## Study architecture

The analysis deliberately separates biological definition from compact
measurement and from downstream validation.

1. A donor-aware discovery analysis identified a **287-gene core**. Membership
   required a Benjamini-Hochberg false-discovery rate of at most 0.05,
   directional stability of at least 0.90, a donor-bootstrap 95% interval that
   excluded zero, and concordant model and bootstrap directions. No fixed
   top-gene count was used.
2. A label-blind portability audit retained **62 protein-coding genes** present
   on five external expression platforms and one FFPE platform.
3. A balanced, leave-one-donor-out reconstruction path represented the
   287-gene score. The knee of its monotone fidelity curve selected six gene
   pairs, producing the frozen **12-gene signature**.
4. Only after the signature was frozen were held-out, external, FFPE,
   RNA-ATAC, atlas, perturbation, spatial and protein results evaluated.

The 287-gene core defines the analysed epithelial state. The 12 genes are an
unweighted research signature for measuring that state; they are not presented
as a uniquely optimal classifier, a self-contained causal circuit or a
clinically validated biomarker.

## Frozen signature

The score is calculated within a dataset as:

```text
mean(z[ASCL2, SOX4, SLC7A5, AXIN2, SDC3, DPEP1])
  - mean(z[LGALS3, COX6C, NDUFA1, CYCS, LRRC19, ETHE1])
```

All 12 genes are required in the reported analyses. The authoritative tables
are in [`data/signatures`](data/signatures).

## Repository layout

```text
analysis/              analysis and R figure-generation scripts
data/datasets.tsv      public dataset and primary-publication index
data/signatures/       frozen 287-, 62- and 12-gene definitions
data/source_data/      source data underlying manuscript figures
data/supplementary/    supplementary tables workbook
figures/communications_biology_v1.2/
                        final submission figures in SVG, PDF and PNG
figures/manuscript/    frozen pre-layout-refinement figures retained for provenance
results/               locked derived results and reproducibility manifests
tests/                  release-integrity and portability checks
workflow/               execution order and reproducibility documentation
```

The versioned script names are retained to preserve the exact analysis trail.
Names containing `jtm` refer to an earlier manuscript target and do not change
the scientific content. The Communications Biology submission uses the same
frozen v2.8 analysis outputs, followed only by journal-specific workflow,
typography and panel-alignment refinements.

## Quick verification

From the repository root:

```bash
conda env create -f environment.yml
conda run -n crc-premalignant-locked python tests/verify_release.py
conda run -n crc-premalignant-locked python analysis/audit_locked_environment.py
```

To regenerate and audit the complete Communications Biology figure package
from the included derived results:

```bash
conda run -n crc-premalignant-locked Rscript \
  analysis/refine_communications_biology_workflows_v1_1.R
conda run -n crc-premalignant-locked Rscript \
  analysis/revise_communications_biology_figure1_v1_2.R
conda run -n crc-premalignant-locked Rscript \
  analysis/refine_communications_biology_figure_audit_fixes_v1_2.R
conda run -n crc-premalignant-locked Rscript \
  analysis/refine_communications_biology_alignment_v1_3.R
conda run -n crc-premalignant-locked Rscript \
  analysis/audit_communications_biology_figures_v1_2.R
```

The renderers export SVG and PDF vector files, 600-dpi TIFF files, 300-dpi PNG
files and panel-level source data. TIFF files are generated locally but omitted
from Git because the corresponding vector and 300-dpi files are distributed.
The expected execution order for full reanalysis is listed in
[`workflow/pipeline.tsv`](workflow/pipeline.tsv).

## Full reanalysis and raw data

Raw third-party datasets are not redistributed. Download them from the public
repositories listed in [`data/datasets.tsv`](data/datasets.tsv), preserve the
directory names documented in [`docs/DATA.md`](docs/DATA.md), and then run the
stages in [`workflow/pipeline.tsv`](workflow/pipeline.tsv). Large single-cell,
ATAC, spatial and proteomic inputs make the full workflow substantially more
computationally intensive than the quick verification.

The GenKI virtual-deletion analysis uses a separate environment because its
PyTorch and NumPy requirements differ from the main analysis:

```bash
conda env create -f environment-virtual-knockout.yml
conda run -n crc-premalignant-virtual-ko python \
  analysis/audit_genki_reproducibility_v2_9.py
```

This analysis is confirmatory and unsigned: it tests whether prespecified
target deletions preferentially affect the frozen gene sets relative to matched
genes. Direction is supplied only by the empirical perturbation datasets.

## Reproducibility

- The signature file has SHA256
  `91e564f42bdc5fb0188605b76116bcfb126e5102d62a33d6890fa67018928380`.
- Patients, donors, models, clones or tissue sections are the inferential units;
  cells, nuclei and spots do not inflate sample size.
- Random seeds, environment versions, input hashes and output checks are
  recorded in the result manifests.
- Figure source data are distributed independently of the raw third-party data.

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for the evidence and
claim boundaries.

## Licence and citation

Code in this repository is released under the MIT License. Derived tables are
provided for research reproducibility; source-dataset licences and terms remain
in force. Cite both the eventual manuscript and the primary publications for
all reused datasets. Repository citation metadata are provided in
[`CITATION.cff`](CITATION.cff).
