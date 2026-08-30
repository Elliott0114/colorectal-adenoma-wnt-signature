#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BiocParallel)
  library(edgeR)
  library(limma)
  library(variancePartition)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
pseudobulk_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "discovery_pseudobulk",
  "discovery_state_pseudobulk.rds"
)
common_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "common_effects",
  "cross_state_common_effects.tsv.gz"
)
primary_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "discovery_models",
  "state_specific_primary_effects.tsv.gz"
)
purity_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "purity_audit",
  "specimen_state_purity_metrics.tsv"
)
out_dir <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "purity_adjusted_sensitivity"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(file.exists(c(
  pseudobulk_path,
  common_path,
  primary_path,
  purity_path
)))) {
  stop("One or more purity-sensitivity inputs are missing")
}

pseudobulk <- readRDS(pseudobulk_path)
common <- read.delim(common_path, check.names = FALSE)
primary <- read.delim(primary_path, check.names = FALSE)
purity <- read.delim(purity_path, check.names = FALSE)
strict_genes <- common$gene[common$strict_state_shared]
states <- c("ABS", "GOB", "TAC")
coefficient <- "routeconventional_adenoma"

metadata <- as.data.frame(pseudobulk$metadata)
metadata <- merge(
  metadata,
  purity[c(
    "specimen_id",
    "cell_type",
    "any_nonepithelial_positive_fraction"
  )],
  by = c("specimen_id", "cell_type"),
  all.x = TRUE,
  sort = FALSE
)
metadata <- metadata[
  match(pseudobulk$metadata$pseudobulk_id, metadata$pseudobulk_id),
  ,
  drop = FALSE
]
if (anyNA(metadata$any_nonepithelial_positive_fraction)) {
  stop("Purity covariate could not be mapped to every pseudobulk")
}

parallel_param <- SnowParam(
  workers = 4L,
  type = "SOCK",
  progressbar = TRUE,
  stop.on.error = FALSE
)
results <- list()
diagnostics <- list()

for (state in states) {
  message("Fitting purity-adjusted sensitivity for ", state)
  started <- Sys.time()
  column_index <- which(metadata$cell_type == state)
  state_metadata <- droplevels(metadata[column_index, , drop = FALSE])
  state_metadata$route <- relevel(factor(state_metadata$route), ref = "normal")
  state_metadata$donor_id <- factor(state_metadata$donor_id)
  state_metadata$nonepithelial_marker_z <- as.numeric(scale(
    state_metadata$any_nonepithelial_positive_fraction
  ))
  state_counts <- as.matrix(
    pseudobulk$counts[strict_genes, column_index, drop = FALSE]
  )
  rownames(state_metadata) <- colnames(state_counts)
  dge <- DGEList(counts = round(state_counts))
  dge <- calcNormFactors(dge, method = "TMM")
  formula <- ~ route + nonepithelial_marker_z + (1 | donor_id)
  voom_object <- voomWithDreamWeights(
    dge,
    formula,
    state_metadata,
    plot = FALSE,
    BPPARAM = parallel_param
  )
  raw_fit <- dream(
    voom_object,
    formula,
    state_metadata,
    ddf = "Satterthwaite",
    BPPARAM = parallel_param
  )
  coefficient_index <- match(coefficient, colnames(raw_fit$coefficients))
  if (is.na(coefficient_index)) {
    stop("Route coefficient is absent for ", state)
  }
  raw_se <- raw_fit$stdev.unscaled[, coefficient_index] * raw_fit$sigma
  moderated_fit <- variancePartition::eBayes(raw_fit, robust = TRUE)
  table <- variancePartition::topTable(
    moderated_fit,
    coef = coefficient,
    number = Inf,
    sort.by = "none"
  )
  table$gene <- rownames(table)
  table$cell_type <- state
  table$raw_se <- raw_se[table$gene]
  table <- table[, c(
    "gene",
    "cell_type",
    "logFC",
    "raw_se",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B",
    "z.std"
  )]
  if (anyNA(table$logFC) || anyNA(table$raw_se)) {
    stop("Purity-adjusted model produced missing estimates for ", state)
  }
  results[[state]] <- table
  diagnostics[[state]] <- data.frame(
    cell_type = state,
    n_genes = nrow(table),
    n_specimens = nrow(state_metadata),
    n_donors = nlevels(state_metadata$donor_id),
    contamination_fraction_min = min(
      state_metadata$any_nonepithelial_positive_fraction
    ),
    contamination_fraction_median = median(
      state_metadata$any_nonepithelial_positive_fraction
    ),
    contamination_fraction_max = max(
      state_metadata$any_nonepithelial_positive_fraction
    ),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    stringsAsFactors = FALSE
  )
}

adjusted <- do.call(rbind, results)
comparison <- merge(
  primary[primary$gene %in% strict_genes, c("gene", "cell_type", "logFC")],
  adjusted[c("gene", "cell_type", "logFC", "adj.P.Val")],
  by = c("gene", "cell_type"),
  suffixes = c("_primary", "_purity_adjusted"),
  sort = FALSE
)
comparison$direction_retained <-
  sign(comparison$logFC_primary) == sign(comparison$logFC_purity_adjusted)
comparison$attenuation_ratio <-
  comparison$logFC_purity_adjusted / comparison$logFC_primary

state_summary <- do.call(rbind, lapply(states, function(state) {
  frame <- comparison[comparison$cell_type == state, , drop = FALSE]
  data.frame(
    cell_type = state,
    strict_genes = nrow(frame),
    direction_retained_fraction = mean(frame$direction_retained),
    spearman_primary_vs_adjusted = cor(
      frame$logFC_primary,
      frame$logFC_purity_adjusted,
      method = "spearman"
    ),
    median_attenuation_ratio = median(frame$attenuation_ratio),
    adjusted_fdr_le_0.05 = sum(frame$adj.P.Val <= 0.05),
    stringsAsFactors = FALSE
  )
}))

gene_summary <- reshape(
  comparison[c(
    "gene",
    "cell_type",
    "logFC_purity_adjusted",
    "adj.P.Val",
    "direction_retained"
  )],
  idvar = "gene",
  timevar = "cell_type",
  direction = "wide"
)
direction_columns <- paste0("direction_retained.", states)
gene_summary$direction_retained_all_states <- apply(
  gene_summary[direction_columns],
  1,
  function(value) all(value %in% TRUE)
)

write.table(
  adjusted,
  gzfile(file.path(out_dir, "purity_adjusted_state_effects.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  comparison,
  gzfile(file.path(out_dir, "primary_vs_purity_adjusted.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  gene_summary,
  gzfile(file.path(out_dir, "purity_adjusted_gene_summary.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  state_summary,
  file.path(out_dir, "purity_adjusted_state_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  do.call(rbind, diagnostics),
  file.path(out_dir, "purity_adjusted_model_diagnostics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_purity_adjusted_sensitivity_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  genes = length(strict_genes),
  formula = "expression ~ route + z(any non-epithelial marker-positive fraction) + (1|donor)",
  scope = "post hoc contamination sensitivity; does not redefine gene membership",
  input_hashes = as.list(setNames(
    vapply(
      c(pseudobulk_path, common_path, primary_path, purity_path),
      digest::digest,
      character(1),
      algo = "sha256",
      file = TRUE
    ),
    basename(c(pseudobulk_path, common_path, primary_path, purity_path))
  )),
  summary = state_summary,
  package_versions = as.list(vapply(
    c("dreamlet", "edgeR", "limma", "variancePartition", "lme4"),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "purity_adjusted_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Purity-adjusted sensitivity complete: ",
  sum(gene_summary$direction_retained_all_states),
  "/",
  nrow(gene_summary),
  " strict genes retained direction in all states"
)
