# Computational closure validation

## Empirical perturbation route effects

```text
         dataset                    comparison  n_units    effect  expected_direction                    status  specificity_p
       GSE114059            trametinib_vs_dmso        4  0.072666                   1 supportive_direction_only       0.255474
       GSE114059 pri724_reversal_of_trametinib        1  0.176235                  -1                discordant       0.942206
        GSE67186         apc_restoration_shApc        1 -0.793567                  -1  exploratory_low_coverage       0.002600
        GSE67186    apc_restoration_shApc_Kras        1 -0.954564                  -1  exploratory_low_coverage       0.018398
        GSE67186 doxycycline_control_shRenilla        1  0.712480                   0           control_nonzero       0.003400
       GSE130822        ascl2_ko_vs_resting_wt        1 -0.553918                  -1       supportive_specific       0.005099
       GSE171910     conditional_wnt_silencing        4 -0.806095                  -1       supportive_specific       0.000100
       GSE125472            APC_vs_WT_with_Wnt        3  0.741998                   1       supportive_specific       0.000100
       GSE125472         APC_vs_WT_without_Wnt        3  2.094181                   1       supportive_specific       0.000100
GSE135328_HCT116               tcf7l2_ko_vs_wt        1 -0.089076                  -1 supportive_direction_only       0.275572
  GSE135328_HT29               tcf7l2_ko_vs_wt        3 -1.341384                  -1       supportive_specific       0.000100
```

## Spatial locked-route effects

```text
                        comparison                          feature  n_sections  median_difference  mean_difference  n_positive  p_wilcoxon
tumor_vs_non_neoplastic_epithelium                      route_score           6           0.473659         0.592297           6    0.031250
tumor_vs_non_neoplastic_epithelium route_residual_prolif_epithelial           6           0.449142         0.464468           6    0.031250
                   tumor_vs_stroma                      route_score          12          -0.031977        -0.056611           1    0.020996
                   tumor_vs_stroma route_residual_prolif_epithelial          12           0.014717         0.003962           8    0.791016
```

All contexts are retained, including discordant and low-information results.
