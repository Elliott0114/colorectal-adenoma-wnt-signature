#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(reticulate)
})

root <- normalizePath(".", mustWork = TRUE)
analysis_dir <- file.path(root, "analysis")
result_root <- file.path(root, "results", "state_aware_program_v1")
source_run_root <- file.path(
  result_root, "functional_architecture_exploratory_v2"
)
run_root <- file.path(
  result_root, "functional_architecture_exploratory_v2_1"
)
selection_path <- file.path(source_run_root, "module_exploratory_routing.tsv")
routing_addendum <- file.path(
  analysis_dir, "contracts",
  "state_aware_functional_architecture_exploratory_v2_1_full_gpl570_2026-08-31.md"
)
if (!all(file.exists(c(selection_path, routing_addendum)))) {
  stop("The v2 route table or v2.1 correction addendum is missing")
}
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

python_binary <- file.path(Sys.getenv("CONDA_PREFIX"), "bin", "python")
if (!file.exists(python_binary)) {
  stop("The active crc-premalignant-locked Python interpreter was not found")
}
use_python(python_binary, required = TRUE)

run_python_stage <- function(filename) {
  message("Running exploratory v2.1 Python stage: ", filename)
  py_run_file(file.path(analysis_dir, filename), local = FALSE, convert = TRUE)
}

run_r_stage <- function(filename) {
  message("Running exploratory v2.1 R stage: ", filename)
  status <- system2(
    file.path(Sys.getenv("CONDA_PREFIX"), "bin", "Rscript"),
    file.path(analysis_dir, filename),
    stdout = "", stderr = ""
  )
  if (!identical(status, 0L)) {
    stop(filename, " failed with exit status ", status)
  }
}

Sys.setenv(
  STATE_AWARE_MODULE_RUN_ROOT = normalizePath(run_root, mustWork = TRUE),
  STATE_AWARE_MODULE_SELECTION_PATH = normalizePath(
    selection_path, mustWork = TRUE
  ),
  STATE_AWARE_MODULE_ROUTE_COLUMN = "analysis_route_pass",
  STATE_AWARE_MODULE_DESCRIPTIVE_DOWNSTREAM = "true",
  STATE_AWARE_MODULE_ROUTING_ADDENDUM_PATH = normalizePath(
    routing_addendum, mustWork = TRUE
  ),
  STATE_AWARE_MODULE_ROUTING_ADDENDUM_SHA256 = digest(
    routing_addendum, algo = "sha256", file = TRUE
  ),
  STATE_AWARE_MODULE_MAIN_LABEL = "exploratory_main_candidate",
  STATE_AWARE_MODULE_SUPPLEMENT_LABEL = "exploratory_supplement_candidate"
)

run_python_stage("state_aware_module_external_validation_v1.py")
run_r_stage("state_aware_module_meta_analysis_v1.R")
run_r_stage("state_aware_module_technical_sentinel_v1.R")
run_python_stage("state_aware_module_orthogonal_context_v1.py")
run_python_stage("state_aware_module_perturbation_protein_v1.py")

message("All exploratory functional-architecture v2.1 stages completed")
