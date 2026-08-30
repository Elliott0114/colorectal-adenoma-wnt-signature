# State-shared epithelial identity remodelling: revision validation contract

Date locked: 2026-08-30

Target manuscript: Communications Biology major-revision candidate

Parent programme: `state_aware_program_v1`

## Scope

This revision does not rederive the 1,843-gene state-shared programme or the frozen eight-gene compact readout. It tests the central interpretation that conventional adenoma is associated with a coordinated change in epithelial identity across comparable epithelial states, rather than a signal explained only by the proportions of recovered cell populations.

No validation outcome may alter programme membership, gene direction, compact-panel membership, panel size or the full-programme scoring rule.

## Primary estimands and gates

### 1. Donor-disjoint held-out replication

- Remove every validation donor represented in the discovery partition before any validation estimate is calculated.
- Recalculate frozen full- and compact-score effects, gene-direction replication and discovery-versus-validation common-effect correlation.
- Gate: the full score must remain positive overall and within ABS, GOB and TAC; at least 90% of testable strict genes must retain the common direction; discovery-versus-validation Spearman correlation must remain at least 0.80.

### 2. Anatomical-site sensitivity

- Attach specimen-level anatomical site from the deposited CELLxGENE metadata without using expression values.
- Primary adjustment: distal/proximal site plus deposited broad state, with donor and specimen repeated measures retained.
- Paired sensitivity: donors contributing conventional adenoma and normal tissue from the same anatomical site.
- Gate: the adjusted full-score coefficient must remain positive; all three state-specific same-site paired estimates must remain positive where at least five paired donors are available.

### 3. Fine-state analysis independent of programme genes

- Fine states are defined separately within ABS, GOB and TAC from discovery normal cells only.
- The feature space excludes all 1,843 programme genes, mitochondrial and ribosomal genes, immunoglobulin genes, and a fixed cell-cycle/stress exclusion list.
- Highly variable features are selected from discovery normal cells without route labels. A normal-reference PCA and k-means model are then frozen and used to assign discovery adenoma and all validation cells to the nearest normal-reference centroid.
- Primary resolution: four fine states per broad state. Sensitivities: three and five fine states.
- A specimen-by-broad-state-by-fine-state stratum requires at least 20 cells.
- Gate: the donor-disjoint validation full-score route coefficient remains positive after fine-state adjustment at all three resolutions; the within-fine-state component is positive and constitutes more than half of the total standardised contrast at the primary resolution.

### 4. Current-programme composition decomposition

- Decompose the frozen 1,843-gene score defined by the final state-aware confidence rule.
- Use donor-route means and normal-reference fine states. Report the observed total difference, composition component, within-fine-state component, and their donor-bootstrap 95% intervals.
- Broad-state-only decomposition is a secondary comparator.
- Gate: the within-state component is positive and larger than the composition component in the donor-disjoint validation partition.

### 5. Independent transcriptomic synthesis and specificity

- Retain cohort-specific effects as the evidence unit.
- Pool standardised adenoma-versus-normal effects with a random-effects model using REML and Hartung-Knapp intervals; report tau-squared, I-squared and a 95% prediction interval.
- Refit leave-one-cohort-out estimates.
- Evaluate available serrated/hyperplastic, carcinoma and grade-defined groups as prespecified boundary analyses. These analyses test context, not diagnostic classification.
- Gate: the pooled adenoma-versus-normal effect and every leave-one-cohort-out estimate remain positive. No claim of adenoma specificity is permitted unless adenoma differs consistently from all available non-adenoma inflammatory/serrated comparators.

### 6. Compact-readout portability

- The eight genes remain frozen.
- Add a fixed single-sample rank score computed from each sample's transcriptome: mean centred rank of the four up genes minus mean centred rank of the four down genes. No disease-label-dependent centring or cohort-wise z-standardisation is permitted.
- Compare this score with the frozen full-programme score in held-out and external datasets. Benchmark against 2,000 arm-balanced random eight-gene panels sampled from the portable 1,843-gene universe.
- Gate: held-out correlation with the full score is at least 0.75 and exceeds the 95th percentile of random panels. This supports reconstruction, not clinical diagnostic use.

### 7. Small-cluster inference

- For Becker and organoid analyses, retain patient/donor as the inference unit.
- Use patient-clustered small-sample corrections where estimable. For three-donor organoid contrasts, report all donor-specific contrasts and the range; do not present a narrow bootstrap interval as population precision.

## Claim discipline

Allowed central claim if all primary gates pass:

> Conventional colorectal adenomas show a donor-reproducible epithelial identity displacement that is shared across canonical epithelial states and remains after accounting for anatomical site and finer within-state composition.

Not allowed from these data:

- independence from every possible epithelial substate;
- adenoma-specific diagnostic performance;
- prospective clinical utility;
- a validated pharmacodynamic endpoint;
- causal mediation of progression;
- population-level precision inferred from three organoid donors.

## Reproducibility

- Random seed: 20260830.
- All new analyses write to `results/state_shared_revision_v2/`.
- Every output directory must contain a manifest with input hashes, package versions, analysis-unit definitions and any failed quality gate.
- The previous manuscript and results remain unchanged; the revised manuscript is a new version.
