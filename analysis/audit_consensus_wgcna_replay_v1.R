#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
set.seed(20260830)
allowWGCNAThreads(nThreads = 8)

root <- normalizePath(".", mustWork = TRUE)
architecture_root <- file.path(
  root, "results", "state_aware_program_v1", "functional_architecture_v1"
)
wgcna_root <- file.path(architecture_root, "consensus_wgcna")
audit_root <- file.path(wgcna_root, "audit")
qc_root <- file.path(root, "qc", "functional_architecture_v1")
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

projection_path <- file.path(wgcna_root, "network_projection_object.rds")
audit_object_path <- file.path(audit_root, "consensus_network_audit_object.rds")
if (!file.exists(projection_path) || !file.exists(audit_object_path)) {
  stop("The completed consensus-network objects are required for replay")
}

projection <- readRDS(projection_path)
reference <- readRDS(audit_object_path)
multi_expr <- lapply(
  projection$discovery$states[c("ABS", "GOB", "TAC")],
  function(value) list(data = as.data.frame(
    value$residual_donor_expression[, projection$network_genes, drop = FALSE]
  ))
)

replay <- blockwiseConsensusModules(
  multi_expr,
  maxBlockSize = 9000,
  randomSeed = 20260830,
  individualTOMInfo = NULL,
  corType = "bicor",
  maxPOutliers = 0.05,
  power = projection$selected_power,
  networkType = "signed",
  TOMType = "signed",
  saveIndividualTOMs = FALSE,
  individualTOMFileNames = file.path(
    qc_root, "replay-set%s-block%b.RData"
  ),
  saveConsensusTOMs = FALSE,
  networkCalibration = "single quantile",
  calibrationQuantile = 0.95,
  consensusQuantile = 0.25,
  useMean = FALSE,
  deepSplit = 2,
  minModuleSize = 30,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  nThreads = 8,
  verbose = 2
)

reference_colors <- as.integer(reference$primary$colors)
replay_colors <- as.integer(replay$colors)
exact_assignment <- identical(reference_colors, replay_colors)
reference_non_grey <- sum(reference_colors != 0L)
replay_non_grey <- sum(replay_colors != 0L)
result <- data.frame(
  check_id = c(
    "same_gene_order", "exact_numeric_module_assignment",
    "same_module_count", "same_non_grey_gene_count"
  ),
  passed = c(
    identical(names(reference$primary$colors), names(replay$colors)),
    exact_assignment,
    length(setdiff(unique(reference_colors), 0L)) ==
      length(setdiff(unique(replay_colors), 0L)),
    reference_non_grey == replay_non_grey
  ),
  detail = c(
    paste("genes", length(reference_colors)),
    paste("exact", exact_assignment),
    paste(
      "reference", length(setdiff(unique(reference_colors), 0L)),
      "replay", length(setdiff(unique(replay_colors), 0L))
    ),
    paste("reference", reference_non_grey, "replay", replay_non_grey)
  ),
  stringsAsFactors = FALSE
)
write.table(
  result, file.path(qc_root, "wgcna_primary_replay.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)
manifest <- list(
  analysis = "audit_consensus_wgcna_replay_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = 20260830,
  passed = all(result$passed),
  input_sha256 = list(
    projection = sha256(projection_path),
    reference_network = sha256(audit_object_path)
  ),
  parameters = list(
    power = projection$selected_power,
    consensus_quantile = 0.25,
    network_type = "signed",
    tom_type = "signed",
    min_module_size = 30,
    deep_split = 2,
    merge_cut_height = 0.25
  ),
  package_versions = list(
    R = as.character(getRversion()),
    WGCNA = as.character(packageVersion("WGCNA"))
  ),
  session_info = capture.output(sessionInfo())
)
write_json(
  manifest, file.path(qc_root, "wgcna_primary_replay_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)
if (!all(result$passed)) {
  stop("The deterministic primary-network replay did not reproduce exactly")
}
message("Primary consensus network replay reproduced the exact assignment")
