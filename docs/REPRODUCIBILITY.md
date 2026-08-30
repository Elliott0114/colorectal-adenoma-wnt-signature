# Reproducibility notes

## Frozen biological object

The primary object is a continuous 8,221-gene common-effect ranking derived
from donor-aware models in absorptive, goblet and transit-amplifying
epithelium. Fixed confidence criteria retain 1,843 genes:

```text
common-effect FDR <= 0.05
+ same posterior direction in all three states
+ local false-sign rate <= 0.05 in every state
= 884 adenoma-up + 959 adenoma-down genes
```

The observed gene count is an output of these rules, not a prespecified top-N
cutoff.

## Independent units

- Discovery and held-out single-cell analyses: donor.
- External transcriptomics: patient cluster or matched patient pair.
- Becker multiome: sample with patient-aware robustness analyses.
- Perturbations: donor, model, edited clone, cell line or oncogenic background.
- Spatial transcriptomics: paired tissue section.

Cells, nuclei and spots are nested measurements rather than independent
replicates.

## Identity-versus-composition analysis

Fine states are learned only from normal discovery epithelium after excluding
all programme genes, lesion labels and predefined nuisance families. All
remaining cells are projected without refitting. Exact donor-balanced
decomposition reconstructs the observed score difference as:

```text
total difference = fine-state composition + expression within fine states
```

The primary four-state solution assigns 79.1% to within-state expression and
20.9% to fine-state composition. Three- and five-state solutions are fixed
sensitivities.

## Outcome-blind full-programme projection

For Becker, CRC Atlas, perturbation and spatial analyses, platform feature
inventories are audited before outcomes are read. A dataset is eligible when
both programme arms have at least 75% coverage and neither arm has fewer than
100 measurable genes. No gene selection, coefficient fitting or
outcome-derived weighting is performed in these layers.

## Candidate reduced readout

The eight-gene candidate is reconstructed from discovery profiles after the
full programme is frozen. Platform availability reduces the candidate universe
to 53 protein-coding genes. Donor-held-out balanced-pair reconstruction and
the data-defined knee select four up/down pairs.

Direction-balanced random panels from the same measurable universe provide a
direct benchmark. The selected candidate modestly exceeds the internal 95th
percentile but not the external 95th percentile. It is therefore interpreted
as one tractable representation of a redundant programme, not as a uniquely
optimal signature.

## Execution

```bash
conda env create -f environment.yml
make verify
make figures
```

The full step order, environment and principal output of each analysis are in
`workflow/pipeline.tsv`. Public accession metadata and primary publications
are in `data/datasets.tsv`.

Supplementary Figure 2g uses governed, pseudonymised patient-level values from
the self-generated DSLab series. The rendered figure and aggregate statistics
are included, but the public `make figures` command retains the committed
Supplementary Figure 2 asset when those governed values are unavailable.

## Governed data

Patient-level DSLab values are excluded from this public release. Aggregate
statistics used for the copy-number-composition sensitivity are included.
