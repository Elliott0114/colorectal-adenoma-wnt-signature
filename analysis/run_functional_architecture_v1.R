#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork = TRUE)
analysis_dir <- file.path(root, "analysis")
rscript <- file.path(Sys.getenv("CONDA_PREFIX"), "bin", "Rscript")
if (!file.exists(rscript) || basename(Sys.getenv("CONDA_PREFIX")) !=
    "crc-premalignant-locked") {
  stop("Run this workflow inside the crc-premalignant-locked Conda environment")
}

stages <- c(
  "state_aware_pathway_replication_v1.R",
  "state_aware_consensus_wgcna_v1.R",
  "audit_consensus_wgcna_replay_v1.R",
  "run_functional_architecture_downstream_v1.R",
  "audit_functional_architecture_v1.R"
)

for (stage in stages) {
  message("Running functional-architecture stage: ", stage)
  stage_environment <- if (stage %in% c(
    "state_aware_consensus_wgcna_v1.R",
    "audit_consensus_wgcna_replay_v1.R"
  )) {
    c(
      "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1",
      "MKL_NUM_THREADS=1", "NUMEXPR_NUM_THREADS=1"
    )
  } else {
    character()
  }
  status <- system2(
    rscript,
    file.path(analysis_dir, stage),
    env = stage_environment,
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L)) {
    stop(stage, " failed with exit status ", status)
  }
}

message("Functional-architecture workflow completed and passed audit")
