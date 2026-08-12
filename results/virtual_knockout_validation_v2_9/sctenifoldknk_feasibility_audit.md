# scTenifoldKnk feasibility audit

## Scope

scTenifoldKnk was evaluated only as an optional second in-silico perturbation engine. It was not required for the prespecified GenKI validation and no scTenifoldKnk output was used in statistical inference or manuscript conclusions.

## Environment

- R package: `scTenifoldKnk` 1.0.3
- Network engine: `scTenifoldNet` 1.4
- Isolated environment: `crc-premalignant-virtual-ko`

## Attempts

1. **Full fixed input (2,198 genes; 1,664 balanced cells):** the official workflow did not complete its first network stage within the monitored feasibility window and produced no analyzable output.
2. **Compact label-blind input (1,025 genes; 1,664 balanced cells):** the universe consisted of the top 800 label-blind highly variable genes plus every measurable prespecified target, all 12 panel genes and all 257 measurable members of the fixed 287-gene core. With the package defaults (`nNet = 10`, `nCells = 500`, `nComp = 3`, `q = 0.9`, `K = 3`, manifold dimension 2), the first seed reached approximately 10% after about 45 minutes and reported an estimated duration of roughly 6 hours per seed.

## Decision

The compact run was stopped before producing a scientific result. Further downsizing was not pursued because removing fixed-core genes or sharply reducing the label-blind background would change the prespecified target set or the reference universe and would therefore weaken, rather than strengthen, a validation-only analysis. The completed dual-seed GenKI analysis, its matched nulls and existing real perturbation calibration remain the sole virtual-perturbation evidence.

## Status

- Installation: completed
- Reproducible input exports: completed
- Complete scTenifoldKnk result: not available
- Included in inferential results: no

