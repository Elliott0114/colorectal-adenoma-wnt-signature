#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork = TRUE)
Sys.setenv(
  STATE_AWARE_MODULE_RUN_ROOT = file.path(
    root, "results", "state_aware_program_v1",
    "functional_architecture_exploratory_v2_1"
  ),
  STATE_AWARE_MODULE_ROUTING_ADDENDUM_PATH = file.path(
    root, "analysis", "contracts",
    "state_aware_functional_architecture_exploratory_v2_1_full_gpl570_2026-08-31.md"
  )
)
source("analysis/build_exploratory_wgcna_integration_v2.R")
