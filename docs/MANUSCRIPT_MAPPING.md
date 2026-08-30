# Manuscript-to-code mapping

| Manuscript element | Primary scripts | Principal outputs |
|---|---|---|
| Fig. 1: study design and state-aware programme derivation | `state_aware_build_discovery_pseudobulk_v1.R`, `state_aware_fit_discovery_models_v1.R`, `state_aware_integrate_common_effects_v1.R`, `state_aware_leave_one_donor_out_stability_v1.R` | `results/state_aware_program_v1/common_effects/`, `donor_leaveout_stability/` |
| Fig. 2: donor-disjoint identity versus composition | `state_shared_revision_donor_site_v2.R`, `state_shared_revision_define_fine_states_v2.py`, `state_shared_revision_fine_state_models_v2.R` | `results/state_shared_revision_v2/donor_site/`, `fine_states/`, `fine_state_models/` |
| Fig. 3: external recurrence and paired FFPE transfer | `validate_state_shared_external_layers_v1.py`, `state_shared_revision_external_meta_v2.R` | `results/state_aware_program_v1/external_validation/`, `results/state_shared_revision_v2/external_meta/` |
| Fig. 4: full-programme regulatory and epithelial contexts | `audit_state_shared_full_program_coverage_v1.py`, `validate_state_shared_full_programme_extended_layers_v2.py` | `results/state_aware_program_v1/full_program_coverage_audit/`, `extended_validation_full_programme/becker*/`, `crc_atlas/` |
| Fig. 5: full-programme genetic perturbations | `validate_state_shared_full_programme_extended_layers_v2.py` | `results/state_aware_program_v1/extended_validation_full_programme/perturbation_spatial/` |
| Fig. 6: candidate reduced readout | `export_state_shared_portability_candidates_v1.py`, `derive_state_shared_compact_panel_v1.R`, `state_shared_revision_compact_rank_v2.R`, `state_shared_revision_external_rank_v2.py` | `results/state_aware_program_v1/panel_derivation/`, `results/state_shared_revision_v2/compact_rank/`, `external_rank/` |
| Supplementary Figs. 1–7 | `build_state_shared_revision_supplementary_figures_s1_s4_v3.R`, `build_state_shared_revision_supplementary_figures_s5_s8_v3.R` | `figures/communications_biology_v2.1/` |

The full 1,843-gene programme defines the biological observation. The
eight-gene candidate appears only after full-programme discovery, validation
and perturbation analyses and is evaluated as one reduced representation.
