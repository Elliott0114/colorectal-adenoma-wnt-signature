# Reproducibility notes

## Frozen biological object

The primary object is a continuous 8,221-gene common-effect ranking from
donor-aware models in absorptive, goblet and transit-amplifying epithelium.
Fixed confidence criteria retain 1,843 genes:

```text
common-effect FDR <= 0.05
+ same posterior direction in all three states
+ local false-sign rate <= 0.05 in every state
= 884 adenoma-up + 959 adenoma-down genes
```

The gene count is the result of fixed confidence rules, not a top-N cutoff.

## Independent units

- Discovery and held-out single-cell analyses: donor.
- External transcriptomics: patient cluster or matched patient pair.
- Becker multiome: sample, with patient-aware robustness analyses.
- Perturbations: donor, model, edited clone, cell line or oncogenic background.
- Spatial transcriptomics: paired tissue section.

Cells, nuclei and spots are nested observations rather than independent
replicates.

## Identity-versus-composition analysis

Fine states are learned only from normal discovery epithelium after excluding
programme genes, lesion labels and predefined nuisance families. Remaining
cells are projected without refitting. Exact donor-balanced decomposition
reconstructs:

```text
total difference = fine-state composition + expression within fine states
```

The primary four-state solution assigns 79.1% to within-state expression and
20.9% to fine-state composition. Three- and five-state solutions are fixed
sensitivities.

## Replicated pathway architecture

`cameraPR` tests the complete signed 8,221-gene ranking in discovery common,
discovery ABS/GOB/TAC, held-out common and held-out ABS/GOB/TAC contexts.
Replication requires discovery-common FDR <= 0.05, held-out-common FDR <= 0.10,
agreement in at least five of six state contexts, and no state-level opposite
direction at FDR <= 0.10. `fgsea` supplies running-enrichment and leading-edge
coordinates only. A fixed leading-edge Jaccard threshold of 0.25 and seeded
Leiden clustering reduce redundant pathways to communities.

## Outcome-blind full-programme projection

For Becker, CRC Atlas, perturbation and spatial analyses, platform feature
inventories are audited before outcomes are read. A dataset is eligible when
both programme arms have at least 75% coverage and neither arm has fewer than
100 measurable genes. No gene selection, coefficient fitting or
outcome-derived weighting is performed in these layers.

## Candidate reduced readout

The supplementary eight-gene candidate is reconstructed from discovery
profiles after the full programme is frozen. Platform availability reduces the
candidate universe to 53 protein-coding genes. Donor-held-out balanced-pair
reconstruction and the data-defined knee select four up/down pairs. External
random-panel benchmarking does not establish unique optimality.

## Execution and governed data

```bash
conda env create -f environment.yml
make verify
make pathway
make figures-public
```

The full order is in `workflow/pipeline.tsv`. Supplementary Figure 2g uses
governed, pseudonymised patient-level DSLab values. Its public release contains
the rendered panel context and aggregate statistics, but not patient-level
source values. Consequently `figures-public` rebuilds only figures whose
complete non-governed derived inputs are public.
