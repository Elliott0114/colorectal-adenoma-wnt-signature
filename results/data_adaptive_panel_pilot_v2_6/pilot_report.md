# Data-adaptive gene-panel pilot (v2.6)

## Scope

This pilot was run in parallel with the v2.5 submission package. No manuscript, figure or table in v2.5 was overwritten. Gene eligibility, forward reconstruction and panel size were determined from Chen discovery data only. Chen validation data were read after the adaptive panel was fixed.

## Main results

- Stable error-controlled core: 287 genes (89 up; 198 down).
- One-SE adaptive panel: 28 genes (21 up; 7 down).
- Discovery grouped out-of-fold fidelity to the full stable core: Spearman rho = 0.956.
- Chen held-out AUC: adaptive 0.917; 100-gene reference 0.929; biology-guided 10-gene 0.936.
- Chen held-out adaptive fidelity to the full stable core: Spearman rho = 0.954.
- Chen held-out paired positive fraction for the adaptive panel: 0.857 (7 donors).
- Cluster-bootstrap AUC difference, adaptive minus 100-gene: median -0.010, 95% CI -0.049 to 0.000.

## Preliminary decision

The discovery-only adaptive method passes the prespecified pilot screen for expansion to the five external cohorts. This is not yet sufficient to replace the manuscript panel.

## Boundary

This is a small-sample exploratory feature-selection analysis. The stable core was fixed from the complete discovery dataset before grouped reconstruction, so grouped out-of-fold fidelity is less conservative than a fully nested re-selection pipeline. Replacement would require re-running the complete selector inside each discovery fold and confirming transportability in all external transcriptomic and FFPE resources.
