#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(Matrix)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
pseudobulk_path <- file.path(
  result_root,
  "discovery_pseudobulk",
  "discovery_state_pseudobulk.rds"
)
state_effect_path <- file.path(
  result_root,
  "discovery_models",
  "state_specific_primary_effects.tsv.gz"
)
common_effect_path <- file.path(
  result_root,
  "common_effects",
  "cross_state_common_effects.tsv.gz"
)
null_correlation_path <- file.path(
  result_root,
  "common_effects",
  "estimated_null_correlation.tsv"
)
contract_path <- file.path(
  root,
  "analysis",
  "contracts",
  "state_aware_program_rederivation_v1_2026-08-29.md"
)
implementation_note_path <- file.path(
  root,
  "analysis",
  "contracts",
  "state_aware_leaveout_stability_implementation_note_v1_2026-08-29.md"
)
out_dir <- file.path(result_root, "donor_leaveout_stability")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_hashes <- c(
  pseudobulk =
    "2c63422dd5cf8124098974fec62648bb24b546adb407ec41d794c8709cdb4c96",
  state_effects =
    "9f2ed71d3f0b01a18a07ae168abaac6249cf49c0fbe9229aaad797c3dbfc5d88",
  common_effects =
    "a1ac4b7b67ac279782e04e386971d7463e169cbdbc24d7da4c314e34a4f3e946",
  null_correlation =
    "b6e0ce7b6cc8350c18b3547883ad1ab5abe8e4f661455ca53423a0ca06d84c78",
  contract =
    "0f39a03154408d14bb0bbe0c1e4f55d498e94d403bd6b45d4be5dd8918555f69",
  implementation_note =
    "32fa283540f7971a1d5431e9e000c77f73b205afbc8d8f73b10c33ecb5344ff8"
)
actual_hashes <- c(
  pseudobulk = digest::digest(pseudobulk_path, algo = "sha256", file = TRUE),
  state_effects = digest::digest(state_effect_path, algo = "sha256", file = TRUE),
  common_effects = digest::digest(
    common_effect_path,
    algo = "sha256",
    file = TRUE
  ),
  null_correlation = digest::digest(
    null_correlation_path,
    algo = "sha256",
    file = TRUE
  ),
  contract = digest::digest(contract_path, algo = "sha256", file = TRUE),
  implementation_note = digest::digest(
    implementation_note_path,
    algo = "sha256",
    file = TRUE
  )
)
if (!identical(expected_hashes, actual_hashes)) {
  stop("A frozen discovery input changed before donor leave-out analysis")
}

pseudobulk <- readRDS(pseudobulk_path)
metadata <- as.data.frame(pseudobulk$metadata)
full_state <- read.delim(state_effect_path, check.names = FALSE)
full_common <- read.delim(common_effect_path, check.names = FALSE)
null_table <- read.delim(null_correlation_path, check.names = FALSE)
states <- c("ABS", "GOB", "TAC")
strict <- full_common[full_common$strict_state_shared, , drop = FALSE]
strict_genes <- sort(strict$gene)
model_genes <- Reduce(
  intersect,
  lapply(states, function(state) {
    full_state$gene[full_state$cell_type == state]
  })
)
model_genes <- rownames(pseudobulk$counts)[
  rownames(pseudobulk$counts) %in% model_genes
]
donors <- sort(unique(as.character(metadata$donor_id)))
coefficient_name <- "routeconventional_adenoma"
if (length(strict_genes) != 1843L ||
    length(model_genes) != 8221L ||
    length(donors) != 27L) {
  stop("Frozen strict-gene or discovery-donor count changed")
}
if (!all(strict_genes %in% rownames(pseudobulk$counts))) {
  stop("A strict gene is absent from the discovery pseudobulk matrix")
}

null_correlation <- as.matrix(null_table[, states, drop = FALSE])
rownames(null_correlation) <- null_table$state
null_correlation <- null_correlation[states, states, drop = FALSE]
if (min(eigen(null_correlation, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) {
  stop("Frozen null-correlation matrix is not positive definite")
}

fit_state_without_donor <- function(state, held_out_donor) {
  index <- which(
    metadata$cell_type == state &
      as.character(metadata$donor_id) != held_out_donor
  )
  local_metadata <- droplevels(metadata[index, , drop = FALSE])
  local_metadata$route <- relevel(factor(local_metadata$route), ref = "normal")
  local_metadata$donor_id <- factor(local_metadata$donor_id)
  counts_all <- as.matrix(
    pseudobulk$counts[model_genes, index, drop = FALSE]
  )
  rownames(local_metadata) <- colnames(counts_all)
  dge_all <- DGEList(counts = round(counts_all))
  dge_all <- calcNormFactors(dge_all, method = "TMM")
  dge <- dge_all[strict_genes, , keep.lib.sizes = TRUE]
  design <- model.matrix(~ route, data = local_metadata)
  voom_object <- voom(dge, design, plot = FALSE)
  correlation <- duplicateCorrelation(
    voom_object,
    design,
    block = local_metadata$donor_id
  )
  fit <- lmFit(
    voom_object,
    design,
    block = local_metadata$donor_id,
    correlation = correlation$consensus.correlation
  )
  coefficient_index <- match(coefficient_name, colnames(fit$coefficients))
  if (is.na(coefficient_index)) {
    stop("Route coefficient is absent in donor leave-out fit")
  }
  standard_error <-
    fit$stdev.unscaled[, coefficient_index] * fit$sigma
  output <- data.frame(
    held_out_donor = held_out_donor,
    cell_type = state,
    gene = rownames(fit$coefficients),
    logFC = fit$coefficients[, coefficient_index],
    raw_se = standard_error,
    n_profiles = nrow(local_metadata),
    n_donors = nlevels(local_metadata$donor_id),
    consensus_donor_correlation = correlation$consensus.correlation,
    held_out_donor_present_in_state = held_out_donor %in%
      as.character(metadata$donor_id[metadata$cell_type == state]),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(output$logFC)) ||
      any(!is.finite(output$raw_se)) ||
      any(output$raw_se <= 0)) {
    stop("Invalid donor leave-out estimate for ", held_out_donor, " / ", state)
  }
  output
}

integrate_fold <- function(state_effects, held_out_donor) {
  bhat <- vapply(
    states,
    function(state) {
      local <- state_effects[state_effects$cell_type == state, , drop = FALSE]
      local$logFC[match(strict_genes, local$gene)]
    },
    numeric(length(strict_genes))
  )
  shat <- vapply(
    states,
    function(state) {
      local <- state_effects[state_effects$cell_type == state, , drop = FALSE]
      local$raw_se[match(strict_genes, local$gene)]
    },
    numeric(length(strict_genes))
  )
  rownames(bhat) <- strict_genes
  rownames(shat) <- strict_genes
  colnames(bhat) <- states
  colnames(shat) <- states
  ones <- rep(1, length(states))
  result <- t(vapply(
    seq_along(strict_genes),
    function(index) {
      covariance <- diag(shat[index, ]) %*%
        null_correlation %*%
        diag(shat[index, ])
      precision <- solve(covariance)
      denominator <- as.numeric(t(ones) %*% precision %*% ones)
      effect <- as.numeric(
        t(ones) %*% precision %*% bhat[index, ] / denominator
      )
      standard_error <- sqrt(1 / denominator)
      c(common_effect = effect, common_se = standard_error, common_z = effect / standard_error)
    },
    numeric(3L)
  ))
  data.frame(
    held_out_donor = held_out_donor,
    gene = strict_genes,
    result,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

message(
  "Running 27 whole-donor leave-out refits for 1,843 frozen genes using three processes"
)
fold_results <- parallel::mclapply(
  donors,
  function(donor) {
    state_effects <- do.call(
      rbind,
      lapply(states, fit_state_without_donor, held_out_donor = donor)
    )
    common_effects <- integrate_fold(state_effects, donor)
    message("completed donor leave-out ", donor)
    list(state = state_effects, common = common_effects)
  },
  mc.cores = 3L,
  mc.preschedule = FALSE,
  mc.set.seed = FALSE
)
failed <- vapply(fold_results, inherits, logical(1), what = "try-error")
if (any(failed)) {
  stop("At least one donor leave-out refit failed")
}
leaveout_state <- do.call(rbind, lapply(fold_results, `[[`, "state"))
leaveout_common <- do.call(rbind, lapply(fold_results, `[[`, "common"))

full_state_strict <- full_state[
  full_state$gene %in% strict_genes,
  c("gene", "cell_type", "logFC"),
  drop = FALSE
]
colnames(full_state_strict)[3L] <- "full_logFC"
leaveout_state <- merge(
  leaveout_state,
  full_state_strict,
  by = c("gene", "cell_type"),
  all.x = TRUE,
  sort = FALSE
)
leaveout_state$direction_match <-
  sign(leaveout_state$logFC) == sign(leaveout_state$full_logFC)

full_common_strict <- strict[, c("gene", "common_effect", "common_z", "shared_direction")]
colnames(full_common_strict)[2:3] <- c("full_common_effect", "full_common_z")
leaveout_common <- merge(
  leaveout_common,
  full_common_strict,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
leaveout_common$direction_match <-
  sign(leaveout_common$common_effect) == sign(leaveout_common$full_common_effect)

state_gene_stability <- aggregate(
  direction_match ~ gene + cell_type,
  data = leaveout_state[leaveout_state$held_out_donor_present_in_state, ],
  FUN = mean
)
state_gene_stability <- tidyr::pivot_wider(
  state_gene_stability,
  names_from = cell_type,
  values_from = direction_match,
  names_prefix = "state_sign_stability_"
)
common_gene_stability <- aggregate(
  direction_match ~ gene,
  data = leaveout_common,
  FUN = mean
)
colnames(common_gene_stability)[2L] <- "common_sign_stability"
common_ranges <- aggregate(
  common_effect ~ gene,
  data = leaveout_common,
  FUN = function(value) c(
    minimum = min(value),
    median = median(value),
    maximum = max(value)
  )
)
common_ranges <- data.frame(
  gene = common_ranges$gene,
  leaveout_common_effect_min = common_ranges$common_effect[, "minimum"],
  leaveout_common_effect_median = common_ranges$common_effect[, "median"],
  leaveout_common_effect_max = common_ranges$common_effect[, "maximum"],
  stringsAsFactors = FALSE
)
gene_stability <- Reduce(
  function(x, y) merge(x, y, by = "gene", all = TRUE, sort = FALSE),
  list(
    full_common_strict,
    state_gene_stability,
    common_gene_stability,
    common_ranges
  )
)
gene_stability$minimum_state_sign_stability <- apply(
  gene_stability[, paste0("state_sign_stability_", states), drop = FALSE],
  1,
  min,
  na.rm = TRUE
)

rank_records <- list()
record_index <- 0L
for (donor in donors) {
  for (state in states) {
    local <- leaveout_state[
      leaveout_state$held_out_donor == donor & leaveout_state$cell_type == state,
      ,
      drop = FALSE
    ]
    record_index <- record_index + 1L
    rank_records[[record_index]] <- data.frame(
      held_out_donor = donor,
      scope = state,
      n_genes = nrow(local),
      effect_spearman = cor(local$logFC, local$full_logFC, method = "spearman"),
      direction_match_fraction = mean(local$direction_match),
      donor_present_in_scope = unique(local$held_out_donor_present_in_state),
      stringsAsFactors = FALSE
    )
  }
  local <- leaveout_common[leaveout_common$held_out_donor == donor, , drop = FALSE]
  record_index <- record_index + 1L
  rank_records[[record_index]] <- data.frame(
    held_out_donor = donor,
    scope = "common_GLS",
    n_genes = nrow(local),
    effect_spearman = cor(
      local$common_effect,
      local$full_common_effect,
      method = "spearman"
    ),
    direction_match_fraction = mean(local$direction_match),
    donor_present_in_scope = TRUE,
    stringsAsFactors = FALSE
  )
}
rank_stability <- do.call(rbind, rank_records)

summary_table <- data.frame(
  metric = c(
    "strict_genes",
    "donor_leaveout_folds",
    "genes_common_sign_stability_1",
    "genes_common_sign_stability_ge_0.90",
    "median_common_sign_stability",
    "minimum_common_sign_stability",
    "median_minimum_state_sign_stability",
    "minimum_minimum_state_sign_stability",
    "median_common_effect_spearman",
    "minimum_common_effect_spearman"
  ),
  value = c(
    nrow(gene_stability),
    length(donors),
    sum(gene_stability$common_sign_stability == 1),
    sum(gene_stability$common_sign_stability >= 0.90),
    median(gene_stability$common_sign_stability),
    min(gene_stability$common_sign_stability),
    median(gene_stability$minimum_state_sign_stability),
    min(gene_stability$minimum_state_sign_stability),
    median(rank_stability$effect_spearman[rank_stability$scope == "common_GLS"]),
    min(rank_stability$effect_spearman[rank_stability$scope == "common_GLS"])
  ),
  stringsAsFactors = FALSE
)

write.table(
  leaveout_state,
  gzfile(file.path(out_dir, "donor_leaveout_state_effects.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  leaveout_common,
  gzfile(file.path(out_dir, "donor_leaveout_common_effects.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  gene_stability,
  file.path(out_dir, "strict_gene_leaveout_stability.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  rank_stability,
  file.path(out_dir, "donor_leaveout_rank_stability.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  summary_table,
  file.path(out_dir, "donor_leaveout_stability_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_leave_one_donor_out_stability_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_hashes = as.list(actual_hashes),
  genes = length(strict_genes),
  donors = length(donors),
  states = states,
  refit = list(
    formula = "voom expression ~ route; donor repeated-measures block",
    normalisation = "TMM",
    weights = "limma::voom",
    donor_correlation = "limma::duplicateCorrelation consensus correlation",
    common_integration = "GLS using frozen full-discovery null correlation",
    parallel_folds = 3L
  ),
  role = "descriptive robustness audit; primary effects remain dream mixed models",
  membership_changed_from_stability = FALSE,
  package_versions = as.list(vapply(
    c("edgeR", "limma", "Matrix"),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "donor_leaveout_stability_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Donor leave-out stability completed; median common sign stability=",
  sprintf("%.3f", median(gene_stability$common_sign_stability))
)
