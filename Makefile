.PHONY: verify environment-audit figures virtual-deletion-audit

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
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_supplementary_figures_s5_s8_v3.R

virtual-deletion-audit:
	conda run -n crc-premalignant-virtual-ko python analysis/audit_genki_reproducibility_v2_9.py
