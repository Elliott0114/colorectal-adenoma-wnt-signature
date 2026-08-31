.PHONY: verify environment-audit pathway network perturbation figures-public

verify:
	conda run -n crc-premalignant-locked python tests/verify_release.py

environment-audit:
	conda run -n crc-premalignant-locked python analysis/audit_locked_environment.py

pathway:
	conda run -n crc-premalignant-locked Rscript analysis/state_aware_pathway_replication_v1.R

network:
	conda run -n crc-premalignant-locked Rscript analysis/run_functional_architecture_exploratory_v2_1.R

perturbation:
	conda run -n crc-premalignant-locked Rscript analysis/adenoma_identity_reversal_target_prioritization_v1.R

figures-public:
	CB_FIGURE_OUT_DIR=figures/communications_biology_v5.0 CB_FIGURE_STEM=figure1_state_shared_derivation conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure1_v3.R
	CB_FIGURE_OUT_DIR=figures/communications_biology_v5.0 CB_FIGURE_STEM=figure2_within_substate_remodelling conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure2_v3.R
	CB_FIGURE_OUT_DIR=figures/communications_biology_v5.0 CB_FIGURE_STEM=figure3_external_tissue_recurrence conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure3_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_communications_biology_v5_figure4.R
	conda run -n crc-premalignant-locked Rscript analysis/build_communications_biology_v5_figure5.R
	conda run -n crc-premalignant-locked Rscript analysis/build_communications_biology_v5_figure6.R
	CB_OUT_DIR=figures/communications_biology_v5.0 conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_supplementary_figures_s1_s4_v3.R
	CB_OUT_DIR=figures/communications_biology_v5.0 CB_SUPPLEMENT_NUMBER_OFFSET=1 conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_supplementary_figures_s5_s8_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_communications_biology_v5_supplementary_figure8.R
