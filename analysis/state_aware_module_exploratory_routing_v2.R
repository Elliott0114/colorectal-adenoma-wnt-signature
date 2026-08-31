#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
source_root <- file.path(result_root, "functional_architecture_v1")
run_root <- file.path(result_root, "functional_architecture_exploratory_v2")
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  module_summary = file.path(
    source_root, "consensus_wgcna", "module_internal_gate_summary.tsv"
  ),
  parameter_stability = file.path(
    source_root, "consensus_wgcna", "audit", "module_parameter_stability.tsv"
  ),
  contract = file.path(
    root, "analysis", "contracts",
    "state_aware_functional_architecture_exploratory_routing_v2_2026-08-31.md"
  )
)
if (!all(file.exists(unlist(paths)))) {
  stop("At least one frozen v1 routing input is missing")
}

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)
summary <- read.delim(paths$module_summary, check.names = FALSE)
stability <- read.delim(paths$parameter_stability, check.names = FALSE)

stability_wide <- reshape(
  stability[, c("module", "sensitivity", "maximum_jaccard")],
  idvar = "module", timevar = "sensitivity", direction = "wide"
)
colnames(stability_wide) <- sub(
  "maximum_jaccard.consensus_quantile_", "jaccard_q",
  colnames(stability_wide), fixed = TRUE
)

routing <- merge(summary, stability_wide, by = "module", all.x = TRUE, sort = FALSE)
routing$original_internal_gate_pass <- routing$internal_gate_pass
routing$analysis_route_pass <- routing$preservation_pass &
  routing$heldout_enrichment_pass
routing$programme_linked <- routing$biological_overlap_pass
routing$parameter_boundary_stable_v1 <- routing$parameter_stable
routing$analysis_route <- ifelse(
  routing$analysis_route_pass,
  ifelse(
    routing$programme_linked,
    "exploratory_programme_linked",
    "exploratory_broader_adenoma_associated"
  ),
  "audit_no_direction_consistent_heldout_association"
)
routing$analysis_version <- "post_result_exploratory_v2"

expected_modules <- c("M02", "M03", "M04", "M05", "M06", "M09", "M10")
selected_modules <- routing$module[routing$analysis_route_pass]
if (!identical(sort(selected_modules), sort(expected_modules))) {
  stop(
    "The frozen v2 entry rule selected an unexpected module set: ",
    paste(sort(selected_modules), collapse = ",")
  )
}

output_path <- file.path(run_root, "module_exploratory_routing.tsv")
write.table(
  routing, output_path, sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_module_exploratory_routing_v2",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  status = "post-result exploratory routing frozen before v2 downstream outcomes",
  route_column = "analysis_route_pass",
  selected_modules = selected_modules,
  input_sha256 = lapply(paths, sha256),
  output_sha256 = list(module_exploratory_routing = sha256(output_path)),
  package_versions = list(
    R = as.character(getRversion()),
    digest = as.character(packageVersion("digest")),
    jsonlite = as.character(packageVersion("jsonlite"))
  )
)
write_json(
  manifest,
  file.path(run_root, "module_exploratory_routing_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)

message(
  "Exploratory v2 routing frozen: ", length(selected_modules),
  " modules selected (", paste(selected_modules, collapse = ", "), ")"
)
