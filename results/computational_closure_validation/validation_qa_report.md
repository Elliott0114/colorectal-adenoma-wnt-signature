# Validation QA report

Date: 2026-08-08  
Mode: validate  
Overall confidence: **CAUTION** — coherent multi-layer computational support,
with small contextual-unit counts, one discordant pharmacological endpoint,
cross-species limitations and no validation experiment in the authors' material.

## Statistical findings audit

| Finding | Effect and evidence | Statistical interpretation | QA status |
|---|---|---|---|
| APC-KO, +WNT | mean Δ 0.742; 2/3 donors; matched-null P 0.000100 | Route is more responsive than expression-matched random signatures, but the donor bootstrap interval crosses zero | CAUTION |
| APC-KO, −WNT | mean Δ 2.094; 3/3 donors; matched-null P 0.000100 | Large and directionally consistent within this three-donor organoid dataset | SOLID within model |
| Conditional WNT-off | mean Δ −0.806; 4/4 models; matched-null P 0.000100 | Reciprocal direction across four distinct CRC models; exact sign P is 0.125 because n=4 | SOLID direction; CAUTION inference |
| Ascl2-KO | Δ −0.554; matched-null P 0.00510 | Specific versus random signatures, but only one contextual contrast (2 KO versus 3 WT libraries) | CAUTION |
| TCF7L2-KO, HT29 | mean Δ −1.341; 3/3 clones; matched-null P 0.000100 | Strong cell-line-specific support | SOLID within model |
| TCF7L2-KO, HCT116 | Δ −0.089; matched-null P 0.276 | Non-specific and weak; not evidence for a universal dependency | RED_FLAG if presented as positive validation |
| Trametinib | route mean Δ 0.073; 2/4 PDOs; matched-null P 0.255 | WNT/stemness rises, but the full route has only direction-level support | CAUTION |
| PRI-724 reversal | route Δ +0.176; WNT/stemness Δ −0.455 | Pathway control reverses; locked route does not | Boundary / discordant |
| Apc restoration | raw Δ −0.794 and −0.955 | Expected direction, but route-arm coverage is 74%/76% and the doxycycline control is non-null | Exploratory only |
| Spatial tumor versus non-neoplastic epithelium | median Δ 0.474; 6/6 sections; exact Wilcoxon P 0.03125 | Paired section-level spatial recapitulation | SOLID within dataset |
| Adjusted spatial route | median Δ 0.449; 6/6 sections; exact Wilcoxon P 0.03125 | Persists after two prespecified nuisance controls; residualization is a sensitivity analysis | SOLID sensitivity |

### Assumption and multiplicity notes

- No t tests or asymptotic normality claims were used for the small contextual
  samples. Exact sign/Wilcoxon tests and unit-level effects are primary.
- A bootstrap interval with one contextual unit is not an inferential interval;
  these entries are treated as descriptive.
- Expression-matched empirical P values quantify signature specificity, not
  population-level biological replication.
- Primary directions and endpoints were frozen before outcomes. Secondary TF
  and control panels involve multiple related readouts and are descriptive; no
  family-wise claim is made from them.
- The two spatial tests shown are the raw primary endpoint and its prespecified
  sensitivity, not independent discoveries.

## Fallacy scan

Coverage: **11/11 checked**.

| Fallacy | Severity | Audit finding | Required handling |
|---|---|---|---|
| 1. Simpson's paradox | NOTE | Cell lines, organoids and donors were not pooled to hide subgroup reversals. APC +WNT includes one discordant donor and trametinib includes two discordant PDOs; both remain shown. | Keep contextual units visible. |
| 2. Ecological fallacy | CAUTION | Section-level spatial differences cannot be interpreted as single-cell effects or patient-level diagnostic performance. | Claim spatial recapitulation only. |
| 3. Berkson's paradox | CAUTION | Public perturbation models are selected experimental systems and do not represent an unselected adenoma population. | Limit generalization to tested models. |
| 4. Collider bias | CAUTION | Proliferation and epithelial-content scores could lie downstream of tumor state; adjusting for them could induce bias. | Keep the unadjusted spatial comparison primary and residualization secondary. |
| 5. Base-rate neglect | NOT APPLICABLE | No sensitivity, specificity, PPV or NPV is reported. | None. |
| 6. Regression to the mean | NOTE | Models were not selected on extreme route scores for pre/post analysis. | None beyond retaining controls. |
| 7. Survivorship bias | NOTE | GEO/Zenodo processed deposits do not expose a study-level attrition analysis; no additional outcome-based sample exclusion was made here. | State reliance on deposited eligible samples. |
| 8. Look-elsewhere effect | CAUTION | The full route and directions were frozen and all tested contexts are reported, but many secondary TF/control outputs exist. | Do not treat secondary panels as independent hypothesis tests. |
| 9. Garden of forking paths | NOTE | A dated contract froze scoring, directions, coverage and units. The doxycycline difference-in-differences was added only after observing a non-null control and is labelled post hoc. | Preserve the contract and label the sensitivity. |
| 10. Correlation ≠ causation | CAUTION | Public genetic interventions are causal manipulations in their source models, but projection of this route and the spatial comparison do not prove causal dependence in human adenoma. | Use “supports” and “recapitulates”, not “proves”. |
| 11. Reverse causality | CAUTION | Temporal direction exists for intervention datasets but not for cross-sectional spatial tissue. | Do not infer route activation caused the spatial tumor state. |

## Reproducibility verification

- Classification: deterministic conditional on fixed inputs, software versions
  and seed 20260808.
- Command: `conda run -n crc-premalignant-locked python analysis/computational_closure_validation.py`
- Final rerun verdict: **REPRODUCIBLE**.
- Non-gzip result tables: exact SHA256 match before versus after rerun.
- Expression-matched null table: exact decompressed-content SHA256 match.
- Spatial spot-score table: exact decompressed-content SHA256 match.
- Figure source tables, TIFF and PNG: exact SHA256 match after an independent R
  rerun.
- All actual analysis inputs, including prior-result tables and the 14 Visium
  HDF5/position/annotation triplets, are enumerated with SHA256 in
  `source_file_manifest.tsv`.

## Environment

- Python 3.11.15; pandas 2.3.3; NumPy 1.26.4; SciPy 1.17.1; h5py 3.13.0;
  openpyxl 3.1.5.
- R 4.5.3; ggplot2 4.0.3; dplyr 1.2.1; tidyr 1.3.2; patchwork 1.3.2;
  svglite 2.2.2; ragg 1.5.2.
- No package was installed, upgraded or modified in the protected environment.

## Final wording verdict

Acceptable: “independent perturbation and spatial analyses support a
WNT–TCF/ASCL2-responsive epithelial route.”

Not acceptable: “the route was experimentally validated”, “PRI-724 rescued the
route”, or “the route causally drives human adenoma progression.”
