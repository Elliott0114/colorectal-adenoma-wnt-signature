# Data acquisition and layout

## Redistribution policy

This repository does not redistribute raw third-party datasets. Public inputs
must be obtained from their original repositories and remain subject to the
original licences and terms. The accession, primary publication and public URL
for every reported evidence layer are indexed in [`../data/datasets.tsv`](../data/datasets.tsv).

The repository does include:

- frozen gene definitions;
- analysis-ready, non-identifying derived result tables;
- figure source data;
- a supplementary workbook;
- checksums and reproducibility manifests; and
- the balanced, de-identified input used for the reported GenKI audit.

Patient-level DSLab matrices, metadata and plotted source values are not
redistributed. Aggregate DSLab statistics are included; governed patient-level
data are available from the corresponding author on reasonable request.

## Expected top-level raw-data directory

Full reanalysis expects a `data_sources/` directory at the repository root. It
is ignored by Git. The principal paths used by the frozen scripts are:

```text
data_sources/
├── Chen_Cell_2021_CELLxGENE/
├── Becker_NatGenet_2022_GEO/
├── CRC_Atlas_CZI_core/
├── external_sporadic_adenoma/
├── GSE117606/
├── GSE125472_apc_ko_organoids/
├── GSE135328_tcf7l2_ko_crc/
├── perturbation_closure/
├── CRC_spatial_public/zenodo_7760264/extracted/
├── public_adenoma_proteomics/
├── reference_gene_sets_v2_5/
└── regulatory_priors/
```

Some download and preprocessing scripts create more specific subdirectories.
The exact paths are constants near the top of each script and are intentionally
kept explicit for auditability.

The state-aware discovery and donor-disjoint validation workflows expect these
two Chen files:

```text
data_sources/Chen_Cell_2021_CELLxGENE/chen_discovery_epithelial.h5ad
data_sources/Chen_Cell_2021_CELLxGENE/chen_validation_epithelial.h5ad
```

Their hashes and versions are recorded in the frozen analysis manifests.

## Licensed reference resources

MSigDB Hallmark release 2026.1.Hs is not redistributed. Obtain the gene-set
file under the applicable MSigDB terms and place it in
`data_sources/reference_gene_sets_v2_5/`. The analysis records the release and
SHA256 checksum. NCBI gene information and mouse-human orthology tables should
likewise be downloaded from their authoritative sources.

## Integrity

Where possible, the frozen scripts write source URLs, byte counts and SHA256
hashes. A changed upstream file should be treated as a new input version rather
than silently substituted for the reported file.
