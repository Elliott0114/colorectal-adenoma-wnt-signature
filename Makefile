.PHONY: verify environment-audit figures cell-state-audit virtual-deletion-audit

verify:
	conda run -n crc-premalignant-locked python analysis/audit_communications_biology_figures_v1_3.py
	conda run -n crc-premalignant-locked python tests/verify_release.py

environment-audit:
	conda run -n crc-premalignant-locked python analysis/audit_locked_environment.py

figures:
	conda run -n crc-premalignant-locked Rscript analysis/plot_jtm_submission_figures_v2_8.R
	conda run -n crc-premalignant-locked Rscript analysis/refine_communications_biology_workflows_v1_1.R
	conda run -n crc-premalignant-locked Rscript analysis/revise_communications_biology_figure1_v1_2.R
	conda run -n crc-premalignant-locked Rscript analysis/refine_communications_biology_figure_audit_fixes_v1_2.R
	conda run -n crc-premalignant-locked Rscript analysis/refine_communications_biology_alignment_v1_3.R
	conda run -n crc-premalignant-locked python analysis/assemble_communications_biology_v1_3.py
	CELL_STATE_FIGURE_DIR=figures/communications_biology_v1.2 CELL_STATE_FIGURE_STEM=figure3_cell_state_decomposition conda run -n crc-premalignant-locked Rscript analysis/plot_cell_state_decomposition_v1.R
	conda run -n crc-premalignant-locked python analysis/audit_communications_biology_figures_v1_3.py

cell-state-audit:
	conda run -n crc-premalignant-locked Rscript analysis/audit_cell_state_decomposition_v1.R

virtual-deletion-audit:
	conda run -n crc-premalignant-virtual-ko python analysis/audit_genki_reproducibility_v2_9.py
