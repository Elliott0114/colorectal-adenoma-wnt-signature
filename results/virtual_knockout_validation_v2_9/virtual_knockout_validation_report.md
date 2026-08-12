# Validation-only virtual-knockout report

## Material Passport

- **Origin Skill:** experiment-agent
- **Origin Mode:** validate
- **Origin Date:** 2026-08-10
- **Verification Status:** VERIFIED
- **Version Label:** virtual_knockout_validation_v2_9

## Validation question

Do prespecified WNT/stem regulatory targets and the frozen 12-gene panel show virtual-knockout coupling to the already-defined 287-gene adenoma state, without using virtual-knockout results to nominate, replace or reweight genes?

## Answer

The virtual-knockout analysis supports the **broader 287-gene programme** as an ASCL2/SOX4-coupled network-responsive state. It gives **modest, sensitivity-dependent support** to the 12-gene panel as a compact readout of that state, but it does **not** support interpreting the 12 genes as a self-contained causal or self-regulating circuit.

This distinction is the appropriate closing claim for the manuscript:

> The 12 genes form a reproducible two-arm expression signature because their fixed, equal-weight score reconstructs and transfers the broader state. Virtual knockout adds orthogonal evidence that this broader state is coupled to prespecified WNT/stem regulators, while the absence of within-panel knockout coherence argues against treating the 12 genes themselves as a causal module.

## Frozen design

- **Input:** conventional-adenoma epithelial cells from the independent Chen validation cohort; 13 donors balanced to 128 cells each (1,664 cells total).
- **Gene universe:** 2,000 label-blind highly variable genes plus measurable frozen genes and prespecified targets (2,198 genes total).
- **Coverage:** 257/287 core genes and 12/12 panel genes passed the prespecified ≥1% cell-detection rule.
- **Targets:** TCF7L2, ASCL2 and SOX4 as upstream-context targets, plus all 12 frozen panel genes (13 unique knockouts).
- **Model:** GenKI 0.2.1, 100 epochs, two initialization seeds, KL divergence as the primary unsigned impact metric and earth mover distance (EMD) as sensitivity analysis. GenKI learns gene embeddings from wild-type single-cell expression and an inferred regulatory network, then removes all edges to and from the target gene for the virtual knockout [as described in the primary method paper](https://academic.oup.com/nar/article/51/13/6578/7184155).
- **Null model:** 10,000 gene sets matched for mean expression, detection fraction and network degree for every target/set/metric test.
- **Multiplicity:** Benjamini–Hochberg correction across the three prespecified primary aggregate KL endpoints and within the prespecified target-level families.
- **Direction:** not inferred from GenKI. The method authors explicitly identify the lack of directional prediction as a limitation; direction was calibrated only against eight existing real APC/WNT/TCF7L2/ASCL2 perturbation contrasts.

## Statistical findings

| Finding | Test and estimate | Interpretation | Confidence |
|---|---|---|---|
| Upstream targets → measurable 287-gene core | KL observed mean impact percentile 0.6484 versus matched-null 0.6302 (95% null interval 0.6222–0.6385); matched-null z = 4.36; empirical P = 9.999 × 10⁻⁵; BH q = 3.000 × 10⁻⁴; 2/3 targets above their matched means | Strong aggregate enrichment relative to abundance- and degree-matched genes | SOLID |
| Upstream targets → leave-target-out 12-gene panel | KL 0.8110 versus 0.7853 (95% null interval 0.7569–0.8122); z = 1.80; P = 0.0308; q = 0.0462; 3/3 targets above matched means | Modest aggregate enrichment, close to the correction boundary | CAUTION |
| Distance-metric sensitivity for core | EMD 0.6230 versus 0.6096; z = 3.04; P = 0.00110 | Confirms the core-level result | SOLID |
| Distance-metric sensitivity for panel | EMD 0.7663 versus 0.7436; z = 1.57; P = 0.0531 | Directionally concordant but does not cross 0.05 | CAUTION |
| Panel-member knockout → leave-target-out panel | KL 0.8126 versus 0.8098; z = 0.338; P = 0.3749; q = 0.3749; EMD P = 0.2890 | No evidence that the 12 genes form an internally coherent knockout-responsive circuit | SOLID negative result; does not prove equivalence |
| ASCL2 knockout → core | KL z = 4.16; empirical P = 0.000200; upstream-family q = 0.000600 | Strong matched enrichment | SOLID |
| SOX4 knockout → core | KL z = 5.10; empirical P = 9.999 × 10⁻⁵; upstream-family q = 0.000600 | Strong matched enrichment | SOLID |
| TCF7L2 knockout → core | KL z = −1.48; empirical P = 0.932; upstream-family q = 0.932 | Not enriched after matching, despite positive unadjusted rank summaries | CAUTION; report the matched result |
| Individual upstream targets → panel | TCF7L2 q = 0.282; ASCL2 q = 0.258; SOX4 q = 0.162 | None is individually significant after family correction | CAUTION |
| Model reconstruction | Edge AUC 0.9483/0.9477 and average precision 0.9410/0.9407 across the two seeds | The learned graph autoencoder reconstructed held-out edges consistently; these metrics are technical QA, not biological validation | SOLID technical QA |
| Rank robustness | KL across-seed Spearman ρ median 0.964 (range 0.539–0.999); KL-versus-EMD ρ median 0.982 (range 0.964–0.996) | Metric choice is highly stable; seed stability is high overall but only moderate for SOX4 | CAUTION for target-specific top-rank interpretation |
| Empirical direction calibration | All 8/8 existing contrasts shifted in the prespecified direction | Supplies directional context that the unsigned GenKI analysis cannot provide | CAUTION because three contrasts contain one model only |

## Warnings that must remain in the manuscript

| Type | Detail | Affected claim |
|---|---|---|
| Unsigned output | GenKI ranks perturbation magnitude in gene-embedding space and cannot predict whether expression increases or decreases | Do not describe GenKI as reproducing the signed 12-gene score response |
| Partial core coverage | Thirty of 287 fixed core genes failed the prespecified detection threshold in this input | Refer to the “257 measurable members of the fixed core,” not the complete core |
| Borderline compact-panel evidence | The primary KL aggregate reaches q = 0.0462, whereas the EMD sensitivity is P = 0.0531 | Use “modest” or “suggestive” support, not definitive validation |
| Target heterogeneity | ASCL2 and SOX4 support the core result; TCF7L2 does not after matched-null adjustment | Do not state that all three upstream nodes independently validate the programme |
| Negative internal coherence | Aggregate within-panel coherence is null and no individual panel-member result survives family correction | Do not present the 12 genes as a causal circuit |
| Unit of inference | Cells were balanced by donor to fit the network; no cell-level P value or donor-level clinical effect is inferred | Do not count 1,664 cells as biological replicates |
| Study/platform scope | The input is a donor-independent validation cohort from the same source study, not a new prospective or cross-platform perturbation cohort | Call this orthogonal computational support, not external experimental validation |
| Low-unit direction contrasts | ASCL2 knockout and two APC-restoration contrasts each contain one model | Treat these directions as descriptive |
| Secondary engine feasibility | scTenifoldKnk was installed but did not complete a scientifically valid fixed-universe run within the feasibility window | Do not cite an incomplete second-engine result |

## Statistical fallacy scan

- **Coverage:** 11/11 fallacy types checked

| Fallacy | Severity | Assessment | Reporting control |
|---|---|---|---|
| Simpson’s paradox | CAUTION | The three-target aggregate masks target heterogeneity: ASCL2 and SOX4 are positive, whereas matched TCF7L2 is not | Show the target-level decomposition beside the aggregate result |
| Ecological fallacy | NOTE | The model operates on a pooled, donor-balanced gene network; it does not estimate individual-patient responses | Restrict inference to fixed gene-set coupling in this epithelial dataset |
| Berkson’s paradox | CAUTION | Selection is restricted to conventional-adenoma epithelial cells and therefore cannot represent all colorectal tissues or disease stages | State the lesion/cell-type scope explicitly |
| Collider bias | NOTE | No outcome-conditioned covariate adjustment was used in the virtual-knockout test; cell-type restriction is part of the estimand rather than a causal adjustment set | Avoid causal language and do not generalize outside the restricted epithelial state |
| Base-rate neglect | NOTE | No sensitivity, specificity, PPV or NPV is estimated in this analysis | Do not use the virtual-knockout result as a diagnostic-performance claim |
| Regression to the mean | NOTE | No group was selected by an extreme pre-perturbation score and no pre/post clinical comparison was made | No specific correction required |
| Survivorship bias | NOTE | All 13 eligible donors met the fixed cell requirement and each contributed the same number of cells | Report the donor-balancing rule and available/selected cell audit |
| Look-elsewhere effect | SOLID control | Targets and gene sets were frozen, all target results are retained and multiple-testing families are explicit | Keep the negative TCF7L2 and within-panel results |
| Garden of forking paths | CAUTION | The analysis was not prospectively registered, but the frozen contract, two metrics, two seeds and no-reselection rule bound the main researcher degrees of freedom | Retain the contract, manifests and feasibility audit |
| Correlation ≠ causation | CAUTION | A virtual edge deletion is a model perturbation, not a biological knockout | Use “virtual-knockout coupling/support,” not “causal validation” |
| Reverse causality | CAUTION | GenKI is unsigned and cannot establish regulatory direction | Obtain direction only from the existing real perturbation datasets |

## Reproducibility

- **Method:** deterministic full re-run with the same input, environment and two seeds.
- **Verdict:** REPRODUCIBLE.
- **Comparison:** 13/13 reported output files matched exactly. Gzip tables were compared after decompression; the manifest was compared after excluding wall-clock duration. Serialized model binaries were not scientific endpoints and were not part of the equality criterion.
- **Runtime:** 107.2 seconds for the verification re-run on the current CPU environment.

## Manuscript-level conclusion

The defensible manuscript conclusion is not that virtual knockout “proves” the signature. It is that a frozen, empirically transferable 12-gene signature represents a broader 287-gene state for which a validation-only GenKI analysis provides independent network-perturbation support, strongest for ASCL2 and SOX4. The null internal-coherence result clarifies the panel’s role: it is a compact measurement construct, not a proposed 12-node mechanism.

