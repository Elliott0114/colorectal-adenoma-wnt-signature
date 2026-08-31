# Manuscript-to-code mapping

| Manuscript element | Primary scripts | Principal outputs |
|---|---|---|
| Fig. 1: resources and derivation of the state-shared response | `state_aware_build_discovery_pseudobulk_v1.R`, `state_aware_fit_discovery_models_v1.R`, `state_aware_integrate_common_effects_v1.R`, `state_aware_leave_one_donor_out_stability_v1.R` | `results/state_aware_program_v1/common_effects/`, `donor_leaveout_stability/` |
| Fig. 2: within-substate remodelling versus composition | `state_shared_revision_donor_site_v2.R`, `state_shared_revision_define_fine_states_v2.py`, `state_shared_revision_fine_state_models_v2.R` | `results/state_shared_revision_v2/donor_site/`, `fine_states/`, `fine_state_models/` |
| Fig. 3: external recurrence and FFPE transfer | `validate_state_shared_external_layers_v1.py`, `state_shared_revision_external_meta_v2.R` | `results/state_aware_program_v1/external_validation/`, `results/state_shared_revision_v2/external_meta/` |
| Fig. 4: replicated pathways and exploratory module architecture | `run_functional_architecture_v1.R`, `state_aware_consensus_wgcna_v1.R`, `state_aware_module_exploratory_routing_v2.R`, `build_exploratory_wgcna_integration_v2_1.R` | `results/state_aware_program_v1/functional_architecture_v1/`, `functional_architecture_exploratory_v2_1/` |
| Fig. 5: multiomic and tissue context | `state_aware_module_orthogonal_context_v1.py`, `state_aware_module_external_validation_v1.py`, `state_aware_module_meta_analysis_v1.R` | `functional_architecture_exploratory_v2_1/module_external_validation/`, `module_orthogonal_context/` |
| Fig. 6: separable perturbation responses | `state_aware_module_perturbation_protein_v1.py`, `adenoma_identity_reversal_target_prioritization_v1.R` | `functional_architecture_exploratory_v2_1/module_perturbation_protein/`, `identity_reversal_target_prioritization_v1/` |
| Supplementary Fig. 6: candidate reduced readout | `export_state_shared_portability_candidates_v1.py`, `derive_state_shared_compact_panel_v1.R`, `state_shared_revision_compact_rank_v2.R`, `state_shared_revision_external_rank_v2.py` | `results/state_aware_program_v1/panel_derivation/`, `results/state_shared_revision_v2/compact_rank/`, `external_rank/` |

The continuous 8,221-gene signed ranking is the primary substrate for pathway
and network analyses. The 1,843 genes retained by fixed directional-confidence
criteria are a high-confidence, state-shared differential-expression subset,
not a natural biological boundary and not a clinical signature. Consensus
modules organise the same response into interpretable functional components;
they do not constitute a second independent discovery track. The supplementary
eight-gene candidate is one compact measurement implementation and does not
define the manuscript's primary conclusion.
