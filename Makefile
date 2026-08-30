.PHONY: verify environment-audit figures

verify:
	conda run -n crc-premalignant-locked python tests/verify_release.py

environment-audit:
	conda run -n crc-premalignant-locked python analysis/audit_locked_environment.py

figures:
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure1_v2.R
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure2_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure3_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure4_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure5_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure6_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_supplementary_figures_s1_s4_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_supplementary_figures_s5_s8_v3.R
