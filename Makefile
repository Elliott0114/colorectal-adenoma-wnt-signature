.PHONY: verify environment-audit figures virtual-deletion-audit

verify:
	conda run -n crc-premalignant-locked python tests/verify_release.py

environment-audit:
	conda run -n crc-premalignant-locked python analysis/audit_locked_environment.py

figures:
	conda run -n crc-premalignant-locked Rscript analysis/plot_jtm_submission_figures_v2_8.R

virtual-deletion-audit:
	conda run -n crc-premalignant-virtual-ko python analysis/audit_genki_reproducibility_v2_9.py

