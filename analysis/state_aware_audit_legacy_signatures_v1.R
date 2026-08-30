#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
set.seed(20260829)

root <- normalizePath(".", mustWork = TRUE)
common_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "common_effects",
  "cross_state_common_effects.tsv.gz"
)
state_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "discovery_models",
  "state_specific_primary_effects.tsv.gz"
)
pseudobulk_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "discovery_pseudobulk",
  "discovery_state_pseudobulk.rds"
)
core_path <- file.path(
  root,
  "release",
  "colorectal-adenoma-wnt-signature",
  "data",
  "signatures",
  "core_287_genes.tsv"
)
panel_path <- file.path(
  root,
  "release",
  "colorectal-adenoma-wnt-signature",
  "data",
  "signatures",
  "signature_12_genes.tsv"
)
contract_path <- file.path(
  root,
  "analysis",
  "contracts",
  "state_aware_program_rederivation_v1_2026-08-29.md"
)
out_dir <- file.path(root, "results", "state_aware_program_v1", "legacy_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_contract_sha256 <-
  "0f39a03154408d14bb0bbe0c1e4f55d498e94d403bd6b45d4be5dd8918555f69"
if (!identical(
  digest::digest(contract_path, algo = "sha256", file = TRUE),
  expected_contract_sha256
)) {
  stop("The frozen parent contract changed before the legacy-signature audit")
}

required_paths <- c(common_path, state_path, pseudobulk_path, core_path, panel_path)
if (any(!file.exists(required_paths))) {
  stop("One or more required inputs are missing")
}

common <- read.delim(common_path, check.names = FALSE)
state_effects <- read.delim(state_path, check.names = FALSE)
pseudobulk <- readRDS(pseudobulk_path)
core <- read.delim(core_path, check.names = FALSE)
panel <- read.delim(panel_path, check.names = FALSE)

states <- c("ABS", "GOB", "TAC")
universe <- common$gene
if (anyDuplicated(universe) || !setequal(universe, unique(state_effects$gene))) {
  stop("The common-effect and state-model gene universes are inconsistent")
}
if (nrow(core) != 287L || nrow(panel) != 12L) {
  stop("The frozen comparator signatures have unexpected sizes")
}
if (!all(core$arm %in% c("up", "down")) || !all(panel$arm %in% c("up", "down"))) {
  stop("Unexpected signature arm label")
}

state_ave_expr <- tapply(
  state_effects$AveExpr,
  state_effects$gene,
  mean,
  na.rm = TRUE
)
count_index <- match(universe, rownames(pseudobulk$counts))
if (anyNA(count_index)) {
  stop("Common-effect genes are absent from the discovery pseudobulk object")
}
detection_fraction <- Matrix::rowMeans(
  pseudobulk$counts[count_index, , drop = FALSE] > 0
)
gene_covariates <- data.frame(
  gene = universe,
  mean_ave_expr = as.numeric(state_ave_expr[universe]),
  detection_fraction = as.numeric(detection_fraction),
  stringsAsFactors = FALSE
)

quantile_bin <- function(value, groups) {
  ranked <- rank(value, ties.method = "average", na.last = "keep")
  pmin(groups, pmax(1L, ceiling(groups * ranked / sum(is.finite(value)))))
}
gene_covariates$expression_bin <- quantile_bin(gene_covariates$mean_ave_expr, 10L)
gene_covariates$detection_bin <- quantile_bin(
  gene_covariates$detection_fraction,
  5L
)
gene_covariates$match_stratum <- paste(
  gene_covariates$expression_bin,
  gene_covariates$detection_bin,
  sep = "_"
)

common_by_gene <- common[match(universe, common$gene), , drop = FALSE]
rownames(common_by_gene) <- common_by_gene$gene

audit_signature <- function(signature, label) {
  expected_sign <- ifelse(signature$arm == "up", 1, -1)
  audit <- data.frame(
    signature = label,
    gene = signature$gene,
    expected_arm = signature$arm,
    expected_sign = expected_sign,
    in_common_universe = signature$gene %in% universe,
    stringsAsFactors = FALSE
  )
  index <- match(audit$gene, common_by_gene$gene)
  audit$common_effect <- common_by_gene$common_effect[index]
  audit$common_z <- common_by_gene$common_z[index]
  audit$common_q_value <- common_by_gene$common_q_value[index]
  audit$max_lfsr <- common_by_gene$max_lfsr[index]
  audit$strict_state_shared <- common_by_gene$strict_state_shared[index]
  audit$relaxed_state_shared <- common_by_gene$relaxed_state_shared[index]
  audit$common_direction_matches <-
    sign(audit$common_effect) == audit$expected_sign
  for (state in states) {
    audit[[paste0("logFC_", state)]] <- common_by_gene[
      index,
      paste0("logFC_", state)
    ]
    audit[[paste0("posterior_mean_", state)]] <- common_by_gene[
      index,
      paste0("posterior_mean_", state)
    ]
    audit[[paste0("lfsr_", state)]] <- common_by_gene[
      index,
      paste0("lfsr_", state)
    ]
    audit[[paste0("direction_matches_", state)]] <-
      sign(audit[[paste0("posterior_mean_", state)]]) == audit$expected_sign
  }
  direction_columns <- paste0("direction_matches_", states)
  audit$all_state_directions_match <- apply(
    audit[direction_columns],
    1,
    function(value) all(value %in% TRUE)
  )
  audit
}

core_audit <- audit_signature(core, "original_287")
panel_audit <- audit_signature(panel, "existing_12")

statistics <- common$common_z
names(statistics) <- common$gene
camera_index <- list(
  original_287_up = match(
    core$gene[core$arm == "up" & core$gene %in% universe],
    names(statistics)
  ),
  original_287_down = match(
    core$gene[core$arm == "down" & core$gene %in% universe],
    names(statistics)
  )
)
camera_result <- cameraPR(
  statistic = statistics,
  index = camera_index,
  use.ranks = TRUE,
  inter.gene.cor = 0.01,
  sort = FALSE,
  directional = TRUE
)
camera_result$gene_set <- rownames(camera_result)
rownames(camera_result) <- NULL
camera_result <- camera_result[, c(
  "gene_set",
  setdiff(colnames(camera_result), "gene_set")
)]

testable_core <- core[core$gene %in% universe, c("gene", "arm"), drop = FALSE]
testable_core$expected_sign <- ifelse(testable_core$arm == "up", 1, -1)
target_covariates <- gene_covariates[
  match(testable_core$gene, gene_covariates$gene),
  ,
  drop = FALSE
]
candidate_covariates <- gene_covariates[
  !gene_covariates$gene %in% core$gene,
  ,
  drop = FALSE
]

target_groups <- split(seq_len(nrow(target_covariates)), target_covariates$match_stratum)
candidate_groups <- split(
  candidate_covariates$gene,
  candidate_covariates$match_stratum
)
target_groups <- target_groups[order(vapply(
  target_groups,
  function(index) {
    length(candidate_groups[[target_covariates$match_stratum[index[1L]]]]) /
      length(index)
  },
  numeric(1)
))]

scaled_covariates <- scale(gene_covariates[c("mean_ave_expr", "detection_fraction")])
rownames(scaled_covariates) <- gene_covariates$gene

sample_matched_set <- function() {
  selected <- rep(NA_character_, nrow(target_covariates))
  used <- character()
  for (indices in target_groups) {
    stratum <- target_covariates$match_stratum[indices[1L]]
    pool <- setdiff(candidate_groups[[stratum]], used)
    if (length(pool) < length(indices)) {
      pool <- setdiff(
        candidate_covariates$gene[
          candidate_covariates$expression_bin ==
            target_covariates$expression_bin[indices[1L]]
        ],
        used
      )
    }
    if (length(pool) < length(indices)) {
      target_centre <- colMeans(
        scaled_covariates[target_covariates$gene[indices], , drop = FALSE]
      )
      available <- setdiff(candidate_covariates$gene, used)
      distance <- rowSums(
        (scaled_covariates[available, , drop = FALSE] -
          matrix(
            target_centre,
            nrow = length(available),
            ncol = length(target_centre),
            byrow = TRUE
          ))^2
      )
      pool <- available[order(distance)][seq_len(max(length(indices), 50L))]
    }
    sampled <- sample(pool, length(indices), replace = FALSE)
    selected[indices] <- sampled
    used <- c(used, sampled)
  }
  if (anyNA(selected) || anyDuplicated(selected)) {
    stop("Matched-set sampler produced missing or duplicated genes")
  }
  selected
}

score_set <- function(genes, arms) {
  z <- statistics[genes]
  up <- z[arms == "up"]
  down <- -z[arms == "down"]
  c(
    up_signed_mean_z = mean(up),
    down_signed_mean_z = mean(down),
    joint_signed_mean_z = mean(c(up, down))
  )
}

observed_metrics <- score_set(testable_core$gene, testable_core$arm)
n_random <- 10000L
message("Generating ", n_random, " expression/detectability-matched random sets")
null_metrics <- matrix(
  NA_real_,
  nrow = n_random,
  ncol = length(observed_metrics),
  dimnames = list(NULL, names(observed_metrics))
)
matching_quality <- matrix(
  NA_real_,
  nrow = n_random,
  ncol = 2L,
  dimnames = list(NULL, c("smd_mean_ave_expr", "smd_detection_fraction"))
)
target_mean <- colMeans(target_covariates[c(
  "mean_ave_expr",
  "detection_fraction"
)])
target_sd <- apply(target_covariates[c(
  "mean_ave_expr",
  "detection_fraction"
)], 2, sd)

for (iteration in seq_len(n_random)) {
  sampled <- sample_matched_set()
  null_metrics[iteration, ] <- score_set(sampled, testable_core$arm)
  sampled_covariates <- gene_covariates[
    match(sampled, gene_covariates$gene),
    c("mean_ave_expr", "detection_fraction"),
    drop = FALSE
  ]
  matching_quality[iteration, ] <-
    (colMeans(sampled_covariates) - target_mean) / target_sd
}

matched_summary <- do.call(rbind, lapply(
  names(observed_metrics),
  function(metric) {
    null <- null_metrics[, metric]
    data.frame(
      metric = metric,
      observed = observed_metrics[[metric]],
      null_mean = mean(null),
      null_sd = sd(null),
      null_ci_low = unname(quantile(null, 0.025)),
      null_ci_high = unname(quantile(null, 0.975)),
      matched_z = (observed_metrics[[metric]] - mean(null)) / sd(null),
      empirical_p_one_sided =
        (1 + sum(null >= observed_metrics[[metric]])) / (n_random + 1),
      stringsAsFactors = FALSE
    )
  }
))

null_output <- data.frame(
  replicate = seq_len(n_random),
  null_metrics,
  matching_quality,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

coverage_summary <- data.frame(
  signature = c("original_287", "original_287_up", "original_287_down", "existing_12"),
  frozen_genes = c(
    nrow(core),
    sum(core$arm == "up"),
    sum(core$arm == "down"),
    nrow(panel)
  ),
  testable_genes = c(
    sum(core$gene %in% universe),
    sum(core$arm == "up" & core$gene %in% universe),
    sum(core$arm == "down" & core$gene %in% universe),
    sum(panel$gene %in% universe)
  ),
  expected_direction_matches = c(
    sum(core_audit$common_direction_matches %in% TRUE),
    sum(core_audit$expected_arm == "up" & core_audit$common_direction_matches %in% TRUE),
    sum(core_audit$expected_arm == "down" & core_audit$common_direction_matches %in% TRUE),
    sum(panel_audit$common_direction_matches %in% TRUE)
  ),
  strict_state_shared = c(
    sum(core_audit$strict_state_shared %in% TRUE),
    sum(core_audit$expected_arm == "up" & core_audit$strict_state_shared %in% TRUE),
    sum(core_audit$expected_arm == "down" & core_audit$strict_state_shared %in% TRUE),
    sum(panel_audit$strict_state_shared %in% TRUE)
  ),
  stringsAsFactors = FALSE
)
coverage_summary$coverage_fraction <-
  coverage_summary$testable_genes / coverage_summary$frozen_genes

write.table(
  core_audit,
  gzfile(file.path(out_dir, "original_287_gene_level_audit.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  panel_audit,
  file.path(out_dir, "existing_12_gene_level_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  camera_result,
  file.path(out_dir, "original_287_competitive_enrichment.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  matched_summary,
  file.path(out_dir, "original_287_matched_null_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  null_output,
  gzfile(file.path(out_dir, "original_287_matched_null_10000.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  coverage_summary,
  file.path(out_dir, "legacy_signature_coverage_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  gene_covariates,
  gzfile(file.path(out_dir, "matched_null_gene_covariates.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_audit_legacy_signatures_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  contract_sha256 = expected_contract_sha256,
  inputs = as.list(setNames(
    vapply(
      required_paths,
      digest::digest,
      character(1),
      algo = "sha256",
      file = TRUE
    ),
    basename(required_paths)
  )),
  universe_genes = length(universe),
  matched_random_replicates = n_random,
  matching_variables = c("mean state-specific AveExpr", "pseudobulk detection fraction"),
  matching_strata = "expression decile by detection quintile; documented nearest-stratum fallback",
  competitive_test = "limma::cameraPR, ranks, inter-gene correlation 0.01",
  random_seed = 20260829,
  coverage = coverage_summary,
  matched_summary = matched_summary,
  package_versions = as.list(vapply(
    c("limma", "Matrix", "digest", "jsonlite"),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "legacy_audit_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Legacy audit complete: ",
  sum(core_audit$common_direction_matches %in% TRUE),
  "/",
  sum(core_audit$in_common_universe),
  " testable 287 genes match the expected common-effect direction"
)
