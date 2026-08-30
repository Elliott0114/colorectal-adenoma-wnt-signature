# Manuscript-to-code mapping

| Manuscript element | Primary scripts | Principal outputs |
|---|---|---|
| Fig. 1: resources and programme derivation | `state_aware_build_discovery_pseudobulk_v1.R`, `state_aware_fit_discovery_models_v1.R`, `state_aware_integrate_common_effects_v1.R`, `state_aware_leave_one_donor_out_stability_v1.R` | `results/state_aware_program_v1/common_effects/`, `donor_leaveout_stability/` |
| Fig. 2: identity remodelling versus composition | `state_shared_revision_donor_site_v2.R`, `state_shared_revision_define_fine_states_v2.py`, `state_shared_revision_fine_state_models_v2.R` | `results/state_shared_revision_v2/donor_site/`, `fine_states/`, `fine_state_models/` |
| Fig. 3: external recurrence and FFPE transfer | `validate_state_shared_external_layers_v1.py`, `state_shared_revision_external_meta_v2.R` | `results/state_aware_program_v1/external_validation/`, `results/state_shared_revision_v2/external_meta/` |
| Fig. 4: replicated functional architecture | `state_aware_pathway_replication_v1.R` | `results/state_aware_program_v1/functional_architecture_v1/pathway_replication*` |
| Fig. 5: RNA–ATAC and atlas context | `audit_state_shared_full_program_coverage_v1.py`, `validate_state_shared_full_programme_extended_layers_v2.py` | `results/state_aware_program_v1/extended_validation_full_programme/becker*/`, `crc_atlas/` |
| Fig. 6: empirical perturbations | `validate_state_shared_full_programme_extended_layers_v2.py` | `results/state_aware_program_v1/extended_validation_full_programme/perturbation_spatial/` |
| Supplementary Fig. 5: candidate reduced readout | `export_state_shared_portability_candidates_v1.py`, `derive_state_shared_compact_panel_v1.R`, `state_shared_revision_compact_rank_v2.R`, `state_shared_revision_external_rank_v2.py` | `results/state_aware_program_v1/panel_derivation/`, `results/state_shared_revision_v2/compact_rank/`, `external_rank/` |

The full 1,843-gene programme is the biological object. The eight-gene
candidate is a supplementary measurement implementation and does not define
the programme or the manuscript's primary conclusion.
