# Outcome-blind addendum: full-programme projection across orthogonal layers

## Purpose

This addendum prespecifies projection of the frozen 1,843-gene state-shared programme into downstream multi-omic, atlas, spatial and perturbation analyses.

## Frozen object

- Programme: 1,843 genes defined before the analyses in this addendum.
- Directional arms: 884 adenoma-up genes and 959 adenoma-down genes.
- Score: mean sample-wise z score of all measurable up-arm genes minus the corresponding mean for all measurable down-arm genes.
- No gene selection, coefficient fitting, outcome-derived weighting or cut-point selection is permitted.

## Outcome-blind feasibility rule

The decision to project the full programme is based only on platform feature inventories. A dataset is eligible when both directional arms have at least 75% feature coverage and neither arm contains fewer than 100 measurable genes. Phenotype labels, effect estimates and P values are not inspected for this decision.

## Prespecified analyses

1. Recalculate the Becker epithelial RNA score and repeat the locked RNA–ATAC associations without changing the accessibility features or covariates.
2. Recalculate CRC Atlas donor–study carrier scores and repeat the locked study-adjusted and leave-one-study-out models.
3. Recalculate the four independent genetic/pharmacological perturbation datasets with the full-programme score and the original unit-level contrasts.
4. Recalculate the Visium spot and section contrasts with the full-programme score and the original proliferation/epithelial adjustment.

The eight-gene candidate remains a downstream measurement-reduction analysis and is not used to define or validate the biological programme in these layers.
