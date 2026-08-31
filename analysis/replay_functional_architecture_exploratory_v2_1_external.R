#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
  library(reticulate)
})

root <- normalizePath(".", mustWork = TRUE)
analysis_dir <- file.path(root, "analysis")
result_root <- file.path(root, "results", "state_aware_program_v1")
canonical_root <- file.path(
  result_root, "functional_architecture_exploratory_v2_1"
)
replay_root <- file.path(
  result_root, "functional_architecture_exploratory_v2_1_replay"
)
selection_path <- file.path(
  result_root, "functional_architecture_exploratory_v2",
  "module_exploratory_routing.tsv"
)
dir.create(replay_root, recursive = TRUE, showWarnings = FALSE)

python_binary <- file.path(Sys.getenv("CONDA_PREFIX"), "bin", "python")
if (!file.exists(python_binary)) {
  stop("The active crc-premalignant-locked Python interpreter was not found")
}
use_python(python_binary, required = TRUE)

Sys.setenv(
  STATE_AWARE_MODULE_RUN_ROOT = normalizePath(replay_root, mustWork = TRUE),
  STATE_AWARE_MODULE_SELECTION_PATH = normalizePath(
    selection_path, mustWork = TRUE
  ),
  STATE_AWARE_MODULE_ROUTE_COLUMN = "analysis_route_pass"
)

py_run_file(
  file.path(analysis_dir, "state_aware_module_external_validation_v1.py"),
  local = FALSE, convert = TRUE
)
status <- system2(
  file.path(Sys.getenv("CONDA_PREFIX"), "bin", "Rscript"),
  file.path(analysis_dir, "state_aware_module_meta_analysis_v1.R"),
  stdout = "", stderr = ""
)
if (!identical(status, 0L)) {
  stop("The replay meta-analysis failed with exit status ", status)
}

relative_outputs <- file.path(
  "module_external_validation",
  c(
    "module_external_cohort_effects.tsv",
    "module_external_random_effects.tsv",
    "module_external_leave_one_cohort_out.tsv",
    "module_ffpe_paired_effects.tsv",
    "module_platform_coverage.tsv",
    "module_gene_platform_measurability.tsv"
  )
)
sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)
comparison <- do.call(rbind, lapply(relative_outputs, function(relative_path) {
  canonical <- file.path(canonical_root, relative_path)
  replay <- file.path(replay_root, relative_path)
  if (!file.exists(canonical) || !file.exists(replay)) {
    stop("A canonical or replay output is missing: ", relative_path)
  }
  canonical_hash <- sha256(canonical)
  replay_hash <- sha256(replay)
  data.frame(
    output = relative_path,
    canonical_sha256 = canonical_hash,
    replay_sha256 = replay_hash,
    identical = identical(canonical_hash, replay_hash),
    stringsAsFactors = FALSE
  )
}))
write.table(
  comparison, file.path(replay_root, "replay_hash_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write_json(
  list(
    analysis = "replay_functional_architecture_exploratory_v2_1_external",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    all_outputs_identical = all(comparison$identical),
    outputs = split(comparison, seq_len(nrow(comparison)))
  ),
  file.path(replay_root, "replay_hash_audit.json"),
  pretty = TRUE, auto_unbox = TRUE
)
if (!all(comparison$identical)) {
  stop("At least one replay output differs from the canonical v2.1 run")
}
message("External-validation replay passed: all six result tables are identical")
