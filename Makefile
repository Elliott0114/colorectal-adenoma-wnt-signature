.PHONY: verify environment-audit pathway figures-public

verify:
	conda run -n crc-premalignant-locked python tests/verify_release.py

environment-audit:
	conda run -n crc-premalignant-locked python analysis/audit_locked_environment.py

pathway:
	conda run -n crc-premalignant-locked Rscript analysis/state_aware_pathway_replication_v1.R

figures-public:
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure1_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure2_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_state_shared_revision_figure3_v3.R
	conda run -n crc-premalignant-locked Rscript analysis/build_functional_architecture_figure4_v1.R
	conda run -n crc-premalignant-locked Rscript analysis/build_full_programme_figure5_v1.R
	conda run -n crc-premalignant-locked Rscript analysis/build_full_programme_figure6_v1.R
	conda run -n crc-premalignant-locked Rscript analysis/build_functional_architecture_supplement_s3_v1.R
	conda run -n crc-premalignant-locked Rscript analysis/build_full_programme_supplement_s4_v1.R
	conda run -n crc-premalignant-locked Rscript analysis/build_reduced_readout_supplement_v4.R
