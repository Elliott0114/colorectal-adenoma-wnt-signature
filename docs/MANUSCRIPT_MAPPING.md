# Manuscript-to-code map

| Manuscript component | Principal code | Locked outputs |
|---|---|---|
| Figure 1 and core definition | `data_adaptive_gene_panel_pilot_v2_6.R`, `export_287_portability_candidates_v2_7.py`, `select_objective_compact_panel_v2_7.R` | `results/data_adaptive_panel_pilot_v2_6/`, `results/objective_compact_panel_v2_7/` |
| Held-out and external validation | `validate_objective_compact_panel_v2_7.py`, `external_sporadic_adenoma_validation.py`, `gse117606_paired_route_validation.py` | `results/objective_compact_panel_v2_7/`, `results/external_sporadic_adenoma_validation/`, `results/gse117606_paired_route_validation/` |
| RNA-ATAC | `becker_rna_atac_concordance.py`, `becker_locked_rna_atac_patient_robustness.py`, `becker_multiome_regulatory_window_accessibility.py` | `results/becker_multiome_regulatory_windows/` and compact-panel extended validation |
| CRC Atlas | `atlas_locked_study_influence.py`, `validate_objective_panel_extended_layers_v2_7.py` | compact-panel extended validation and `results/figure_data_locked/` |
| Empirical perturbations | `perturbation_validation_locked_route.py`, `computational_closure_validation.py` | `results/perturbation_validation_locked_route/`, `results/computational_closure_validation/` |
| Virtual deletion | `prepare_genki_virtual_knockout_validation_v2_9.py`, `run_genki_virtual_knockout_validation_v2_9.py`, `audit_genki_reproducibility_v2_9.py` | `results/virtual_knockout_validation_v2_9/` |
| Spatial and protein context | `spatial_zenodo7760264_visium_analysis.py`, `public_adenoma_protein_triangulation.py`, `pxd000445_candidate_reanalysis.py` | `results/computational_closure_validation/`, `results/public_adenoma_protein_triangulation/`, `results/pxd000445_candidate_reanalysis/` |
| All manuscript figures | `plot_jtm_submission_figures_v2_8.R` | `figures/jtm_submission_v2.8/` and `data/source_data/` |

Earlier versioned scripts are retained only where the final renderer imports
their audited visual components or where they document a necessary provenance
step. The v2.8 gene definitions and outputs are authoritative.

