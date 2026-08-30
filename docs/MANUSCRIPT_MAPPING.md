# Manuscript-to-code map

| Manuscript component | Principal code | Frozen outputs |
|---|---|---|
| Fig. 1: state-aware derivation and interpretation | `state_aware_build_discovery_pseudobulk_v1.R`, `state_aware_fit_discovery_models_v1.R`, `state_aware_integrate_common_effects_v1.R`, `state_aware_leave_one_donor_out_stability_v1.R`, `state_aware_interpret_common_program_v1.R` | `results/state_aware_program_v1/discovery_*`, `common_effects/`, `donor_leaveout_stability/`, `interpretation/` |
| Fig. 2: donor-disjoint replication and fine-state decomposition | `state_aware_validate_frozen_program_v1.R`, `state_shared_revision_donor_site_v2.R`, `state_shared_revision_define_fine_states_v2.py`, `state_shared_revision_fine_state_models_v2.R` | `results/state_aware_program_v1/heldout_validation/`, `results/state_shared_revision_v2/donor_site/`, `fine_states/`, `fine_state_models/` |
| Fig. 3: independent cohorts and FFPE | `validate_state_shared_external_layers_v1.py`, `state_shared_revision_external_meta_v2.R` | `results/state_aware_program_v1/external_validation/`, `results/state_shared_revision_v2/external_meta/` |
| Fig. 4: RNA-ATAC and epithelial atlas | `validate_state_shared_extended_layers_v1.py` | `results/state_aware_program_v1/extended_validation/becker*`, `crc_atlas/` |
| Fig. 5: APC-WNT genetic perturbations | `validate_state_shared_external_layers_v1.py`, `validate_state_shared_extended_layers_v1.py` | `results/state_aware_program_v1/external_validation/`, `extended_validation/perturbation_spatial/` |
| Fig. 6: compact eight-gene measurement | `export_state_shared_portability_candidates_v1.py`, `derive_state_shared_compact_panel_v1.R`, `state_shared_revision_compact_rank_v2.R`, `state_shared_revision_external_rank_v2.py` | `results/state_aware_program_v1/panel_derivation/`, `results/state_shared_revision_v2/compact_rank/`, `external_rank/` |
| DSLab composition sensitivity | `run_dslab_state_shared_cnv_validation_v1.py` | aggregate statistics in `results/state_aware_program_v1/dslab_cnv_validation/`; patient-level outputs are not public |
| Historical 287/12-gene audit | `state_aware_audit_legacy_signatures_v1.R` | `results/state_aware_program_v1/legacy_audit/` |
| Virtual-deletion context | `prepare_genki_virtual_knockout_validation_v2_9.py`, `run_genki_virtual_knockout_validation_v2_9.py`, `audit_genki_reproducibility_v2_9.py` | `results/virtual_knockout_validation_v2_9/` |
| Main figure rendering | `build_state_shared_revision_figure1_v2.R`, `build_state_shared_revision_figure2_v3.R` through `build_state_shared_revision_figure6_v3.R` | `figures/communications_biology_v2.0/figure1*` through `figure6*` |
| Supplementary figure rendering | `build_state_shared_revision_supplementary_figures_s1_s4_v3.R`, `build_state_shared_revision_supplementary_figures_s5_s8_v3.R` | `figures/communications_biology_v2.0/figureS1*` through `figureS8*` |

The continuous 8,221-gene ranking is used for pathway interpretation. The
1,843-gene subset is defined by fixed confidence rules, and the eight-gene
candidate is derived only after the biological programme has been frozen.
Historical 287- and 12-gene files are maintained solely for auditability.
