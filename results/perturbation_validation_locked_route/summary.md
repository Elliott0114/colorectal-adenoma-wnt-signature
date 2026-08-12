# Perturbation validation summary

## Primary empirical perturbation

- GSE125472 contains three independent donor-matched isogenic human colon organoid systems.
- APC-KO minus WT route-score change with Wnt/R-spondin: mean 0.7420, median 1.0976, range -0.4632 to 1.5916; expected direction in 2/3 donors.
- Expression-matched random-signature test: P=9.999e-05 (10,000 permutations).
- APC-KO minus WT without Wnt/R-spondin: mean 2.0942; expected direction in 3/3 donors.
- WT Wnt/R-spondin withdrawal: mean -1.5394; expected negative direction in 3/3 donors.

## TCF7L2 stress test

```text
cell_line  clone_id  difference_vs_WT
   HCT116 HCT116_18         -0.089076
     HT29   HT29_57         -1.209350
     HT29   HT29_83         -0.914394
     HT29   HT29_86         -1.900408
```

TCF7L2 regulon-activity calibration:

```text
cell_line  clone_id  difference_vs_WT
   HCT116 HCT116_18          2.545304
     HT29   HT29_57         -4.151301
     HT29   HT29_83         -2.558806
     HT29   HT29_86         -2.666616
```

## Evidential interpretation

The APC organoid analysis is the main perturbational validation because it uses the relevant human colonic epithelium, isogenic editing and donor matching. The TCF7L2 analysis is a context-dependent falsification layer and must not be pooled across cell lines as if clones were patient replicates. The signed virtual knockout is topology based and is supportive only.
