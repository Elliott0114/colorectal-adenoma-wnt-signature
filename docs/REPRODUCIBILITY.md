# Reproducibility notes

## Primary transcriptomic response

The primary analytical substrate is a continuous 8,221-gene common-effect
ranking from donor-aware models in absorptive, goblet and transit-amplifying
epithelium. Fixed confidence criteria identify 1,843 high-confidence,
state-shared differentially expressed genes:

```text
common-effect FDR <= 0.05
+ same posterior direction in all three states
+ local false-sign rate <= 0.05 in every state
= 884 adenoma-up + 959 adenoma-down genes
```

The count is the consequence of the fixed rules, not a top-N cutoff. Pathway
and network analyses use the complete 8,221-gene universe where specified.

## Independent units

- Discovery and held-out single-cell analyses: donor.
- External transcriptomics: patient cluster or matched patient pair.
- Becker multiome: sample, with patient-aware robustness analyses.
- Perturbations: donor, model, edited clone, cell line or oncogenic background.
- Spatial transcriptomics: paired tissue section.

Cells, nuclei and spots are nested observations rather than independent
replicates.

## Within-substate versus compositional analysis

Normal-reference epithelial substates are learned only from normal discovery
epithelium after excluding the 1,843 high-confidence genes, lesion labels and
predefined nuisance families. Remaining cells are projected without refitting.
Exact donor-balanced decomposition reconstructs:

```text
total difference = substate composition + expression within substates
```

The primary four-substate solution assigns 79.1% to within-substate expression
and 20.9% to substate composition. Three- and five-substate solutions are fixed
sensitivities.

## Replicated pathway architecture

`cameraPR` tests the complete signed 8,221-gene ranking in discovery common,
discovery ABS/GOB/TAC, held-out common and held-out ABS/GOB/TAC contexts.
Replication requires discovery-common FDR <= 0.05, held-out-common FDR <= 0.10,
agreement in at least five of six state contexts, and no state-level opposite
direction at FDR <= 0.10. `fgsea` supplies running-enrichment and leading-edge
coordinates only. A fixed leading-edge Jaccard threshold of 0.25 and seeded
Leiden clustering reduce redundant pathways to communities.

## Exploratory consensus modules

Consensus WGCNA is fitted to all 8,221 testable genes after residualising the
normal/adenoma main effect within each state. The network therefore asks how the
same transcriptomic response is functionally organised across independent
donors; it is not a second differential-expression screen. Modules are
described only when their direction recurs in held-out rankings, external
cohorts and orthogonal contexts. The corresponding scripts retain preservation,
technical-correlation and parameter-sensitivity audits whether or not a module
is routed to the manuscript.

## Two-coordinate perturbation projection

Genetic and pharmacological contrasts are projected without refitting onto two
prespecified biological coordinates: suppression of the adenoma-up WNT/stem
component and restoration of adenoma-down mature epithelial functions. The
module-by-perturbation matrix records all included contrasts, including
discordant and null responses. The analysis distinguishes selective regulatory
suppression from coordinated reversal; it does not infer clinical efficacy.

## Candidate reduced readout

The supplementary eight-gene candidate is reconstructed from discovery
profiles after the high-confidence response is frozen. Platform availability
reduces the candidate universe to 53 protein-coding genes. Donor-held-out
balanced-pair reconstruction and the data-defined knee select four up/down
pairs. External random-panel benchmarking does not establish unique
optimality.

## Execution and governed data

```bash
conda env create -f environment.yml
make verify
make pathway
make network
make perturbation
make figures-public
```

The full order is in `workflow/pipeline.tsv`. Supplementary Figure 2g uses
governed, pseudonymised patient-level DSLab values. Its public release contains
the rendered panel context and aggregate statistics, but not patient-level
source values. Consequently `figures-public` rebuilds only figures whose
complete non-governed derived inputs are public.
