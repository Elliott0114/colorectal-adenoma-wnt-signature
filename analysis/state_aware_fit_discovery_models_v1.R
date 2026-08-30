#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BiocParallel)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(variancePartition)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
arguments <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(arguments)) arguments[1L] else "discovery"
if (!dataset %in% c("discovery", "validation")) {
  stop("Dataset must be either discovery or validation")
}
input_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  paste0(dataset, "_pseudobulk"),
  paste0(dataset, "_state_pseudobulk.rds")
)
contract_path <- file.path(
  root,
  "analysis",
  "contracts",
  "state_aware_program_rederivation_v1_2026-08-29.md"
)
addendum_path <- file.path(
  root,
  "analysis",
  "contracts",
  "state_aware_model_implementation_addendum_v1_2026-08-29.md"
)
out_dir <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  paste0(dataset, "_models")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_contract_sha256 <-
  "0f39a03154408d14bb0bbe0c1e4f55d498e94d403bd6b45d4be5dd8918555f69"
actual_contract_sha256 <- digest::digest(contract_path, algo = "sha256", file = TRUE)
if (!identical(actual_contract_sha256, expected_contract_sha256)) {
  stop("The parent analysis contract changed after pseudobulk construction")
}

input <- readRDS(input_path)
if (!identical(input$contract_sha256, expected_contract_sha256)) {
  stop("Pseudobulk object was not generated under the frozen parent contract")
}

counts <- input$counts
metadata <- as.data.frame(input$metadata)
if (!inherits(counts, "sparseMatrix")) {
  stop("Pseudobulk counts are not stored as a sparse Matrix object")
}
if (!identical(colnames(counts), metadata$pseudobulk_id)) {
  stop("Count columns and metadata rows are misaligned")
}
if (any(abs(counts@x - round(counts@x)) > 1e-8)) {
  stop("Pseudobulk matrix contains non-integer values")
}

states <- c("ABS", "GOB", "TAC")
coefficient_name <- "routeconventional_adenoma"
workers <- 4L
parallel_param <- SnowParam(
  workers = workers,
  type = "SOCK",
  progressbar = TRUE,
  stop.on.error = FALSE
)

filter_records <- list()
genes_by_state <- list()
for (state in states) {
  column_index <- which(metadata$cell_type == state)
  state_metadata <- droplevels(metadata[column_index, , drop = FALSE])
  state_metadata$route <- relevel(factor(state_metadata$route), ref = "normal")
  state_metadata$donor_id <- factor(state_metadata$donor_id)
  state_counts <- as.matrix(counts[, column_index, drop = FALSE])
  rownames(state_metadata) <- colnames(state_counts)

  dge <- DGEList(counts = round(state_counts))
  filter_design <- model.matrix(~ route, data = state_metadata)
  keep <- filterByExpr(dge, design = filter_design)
  genes_by_state[[state]] <- rownames(dge)[keep]
  filter_records[[state]] <- data.frame(
    gene = rownames(dge),
    cell_type = state,
    passes_filter = keep,
    stringsAsFactors = FALSE
  )
  message(
    state,
    ": ",
    sum(keep),
    " of ",
    length(keep),
    " genes pass expression-only filtering"
  )
}

shared_genes <- Reduce(intersect, genes_by_state)
union_genes <- Reduce(union, genes_by_state)
if (length(shared_genes) < 1000L) {
  stop("Unexpectedly few genes pass filtering in all three states")
}
shared_genes <- rownames(counts)[rownames(counts) %in% shared_genes]
message("Shared primary universe: ", length(shared_genes), " genes")

filter_table <- do.call(rbind, filter_records)
filter_wide <- tidyr::pivot_wider(
  filter_table,
  names_from = cell_type,
  values_from = passes_filter,
  names_prefix = "passes_"
)
filter_wide$in_shared_universe <- filter_wide$gene %in% shared_genes
filter_wide$in_union_universe <- filter_wide$gene %in% union_genes

state_results <- list()
state_fits <- list()
model_diagnostics <- list()
normalization_records <- list()

for (state in states) {
  message("Fitting primary mixed model for ", state)
  started <- Sys.time()
  column_index <- which(metadata$cell_type == state)
  state_metadata <- droplevels(metadata[column_index, , drop = FALSE])
  state_metadata$route <- relevel(factor(state_metadata$route), ref = "normal")
  state_metadata$donor_id <- factor(state_metadata$donor_id)
  state_counts <- as.matrix(counts[shared_genes, column_index, drop = FALSE])
  rownames(state_metadata) <- colnames(state_counts)

  dge <- DGEList(counts = round(state_counts))
  dge <- calcNormFactors(dge, method = "TMM")
  formula <- ~ route + (1 | donor_id)
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

  if (!coefficient_name %in% colnames(raw_fit$coefficients)) {
    stop("Expected route coefficient is absent for ", state)
  }
  coefficient_index <- match(coefficient_name, colnames(raw_fit$coefficients))
  raw_standard_error <-
    raw_fit$stdev.unscaled[, coefficient_index] * raw_fit$sigma

  moderated_fit <- variancePartition::eBayes(raw_fit, robust = TRUE)
  table <- variancePartition::topTable(
    moderated_fit,
    coef = coefficient_name,
    number = Inf,
    sort.by = "none"
  )
  table$gene <- rownames(table)
  table$cell_type <- state
  table$raw_se <- as.numeric(raw_standard_error[table$gene])
  table$raw_t <- as.numeric(raw_fit$t[table$gene, coefficient_index])
  table$raw_p_value <- as.numeric(raw_fit$p.value[table$gene, coefficient_index])
  table$moderated_se <- abs(table$logFC / table$t)
  table$moderated_se[!is.finite(table$moderated_se)] <- NA_real_
  table <- table[, c(
    "gene",
    "cell_type",
    "logFC",
    "raw_se",
    "raw_t",
    "raw_p_value",
    "moderated_se",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B",
    "z.std"
  )]

  model_errors <- attr(raw_fit, "errors")
  initial_errors <- attr(raw_fit, "error.initial")
  n_errors <- if (is.null(model_errors)) 0L else length(model_errors)
  n_initial_errors <- if (is.null(initial_errors)) 0L else length(initial_errors)
  if (anyNA(table$logFC) || anyNA(table$raw_se) || anyNA(table$P.Value)) {
    stop("Primary model produced missing route estimates for ", state)
  }
  if (any(table$raw_se <= 0)) {
    stop("Primary model produced non-positive standard errors for ", state)
  }

  state_results[[state]] <- table
  state_fits[[state]] <- list(
    raw_fit = raw_fit,
    moderated_fit = moderated_fit,
    voom_object = voom_object
  )
  normalization_records[[state]] <- data.frame(
    pseudobulk_id = colnames(dge),
    cell_type = state,
    library_size = dge$samples$lib.size,
    normalization_factor = dge$samples$norm.factors,
    effective_library_size = dge$samples$lib.size * dge$samples$norm.factors,
    stringsAsFactors = FALSE
  )
  model_diagnostics[[state]] <- data.frame(
    cell_type = state,
    n_genes = nrow(table),
    n_specimens = nrow(state_metadata),
    n_donors = nlevels(state_metadata$donor_id),
    n_normal_specimens = sum(state_metadata$route == "normal"),
    n_adenoma_specimens = sum(state_metadata$route == "conventional_adenoma"),
    n_model_errors = n_errors,
    n_initial_model_errors = n_initial_errors,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    stringsAsFactors = FALSE
  )
  message(
    "Completed ",
    state,
    " in ",
    round(model_diagnostics[[state]]$elapsed_seconds, 1),
    " seconds"
  )
}

state_effects <- do.call(rbind, state_results)
normalization_table <- do.call(rbind, normalization_records)
diagnostics_table <- do.call(rbind, model_diagnostics)

paired_results <- list()
for (state in states) {
  message("Fitting paired-donor sensitivity for ", state)
  column_index <- which(metadata$cell_type == state)
  state_metadata <- droplevels(metadata[column_index, , drop = FALSE])
  donor_route_presence <- table(state_metadata$donor_id, state_metadata$route) > 0
  paired_donors <- rownames(donor_route_presence)[rowSums(donor_route_presence) == 2L]
  keep_columns <- column_index[metadata$donor_id[column_index] %in% paired_donors]
  paired_metadata <- droplevels(metadata[keep_columns, , drop = FALSE])
  paired_metadata$route <- relevel(factor(paired_metadata$route), ref = "normal")
  paired_metadata$donor_id <- factor(paired_metadata$donor_id)
  paired_counts <- counts[shared_genes, keep_columns, drop = FALSE]

  donor_route <- interaction(
    paired_metadata$donor_id,
    paired_metadata$route,
    drop = TRUE,
    lex.order = TRUE,
    sep = "__"
  )
  aggregation <- sparseMatrix(
    i = seq_along(donor_route),
    j = as.integer(donor_route),
    x = 1,
    dims = c(length(donor_route), nlevels(donor_route))
  )
  donor_route_counts <- as.matrix(paired_counts %*% aggregation)
  colnames(donor_route_counts) <- levels(donor_route)
  donor_route_metadata <- unique(data.frame(
    donor_route = as.character(donor_route),
    donor_id = as.character(paired_metadata$donor_id),
    route = as.character(paired_metadata$route),
    stringsAsFactors = FALSE
  ))
  donor_route_metadata <- donor_route_metadata[
    match(colnames(donor_route_counts), donor_route_metadata$donor_route),
    ,
    drop = FALSE
  ]
  donor_route_metadata$donor_id <- factor(donor_route_metadata$donor_id)
  donor_route_metadata$route <- relevel(factor(donor_route_metadata$route), ref = "normal")
  rownames(donor_route_metadata) <- donor_route_metadata$donor_route

  dge <- DGEList(counts = round(donor_route_counts))
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~ donor_id + route, data = donor_route_metadata)
  voom_object <- voom(dge, design, plot = FALSE)
  fit <- lmFit(voom_object, design)
  fit <- limma::eBayes(fit, robust = TRUE)
  coefficient_index <- match(coefficient_name, colnames(design))
  if (is.na(coefficient_index)) {
    stop("Paired model route coefficient is absent for ", state)
  }
  table <- limma::topTable(
    fit,
    coef = coefficient_index,
    number = Inf,
    sort.by = "none"
  )
  table$gene <- rownames(table)
  table$cell_type <- state
  table$n_paired_donors <- length(paired_donors)
  table <- table[, c(
    "gene",
    "cell_type",
    "n_paired_donors",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )]
  paired_results[[state]] <- table
}
paired_effects <- do.call(rbind, paired_results)

write.table(
  state_effects,
  gzfile(file.path(out_dir, "state_specific_primary_effects.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  paired_effects,
  gzfile(file.path(out_dir, "state_specific_paired_sensitivity.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter_wide,
  gzfile(file.path(out_dir, "expression_filter_audit.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  normalization_table,
  file.path(out_dir, "tmm_normalization_factors.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  diagnostics_table,
  file.path(out_dir, "model_diagnostics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  list(
    state_fits = state_fits,
    shared_genes = shared_genes,
    union_genes = union_genes,
    coefficient = coefficient_name,
    parent_contract_sha256 = actual_contract_sha256,
    addendum_sha256 = digest::digest(addendum_path, algo = "sha256", file = TRUE)
  ),
  file.path(out_dir, paste0(dataset, "_state_model_fits.rds")),
  compress = "xz"
)

manifest <- list(
  analysis = paste0("state_aware_fit_", dataset, "_models_v1"),
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input = list(
    path = normalizePath(input_path),
    sha256 = digest::digest(input_path, algo = "sha256", file = TRUE),
    pseudobulks = ncol(counts),
    genes_before_filtering = nrow(counts)
  ),
  contracts = list(
    parent_path = normalizePath(contract_path),
    parent_sha256 = actual_contract_sha256,
    addendum_path = normalizePath(addendum_path),
    addendum_sha256 = digest::digest(addendum_path, algo = "sha256", file = TRUE)
  ),
  model = list(
    formula = "expression ~ route + (1 | donor_id)",
    coefficient = coefficient_name,
    normalisation = "TMM",
    expression_filter = "edgeR::filterByExpr default",
    shared_universe_genes = length(shared_genes),
    union_universe_genes = length(union_genes),
    degrees_of_freedom = "Satterthwaite",
    empirical_bayes = "variancePartition::eBayes robust=TRUE",
    workers = workers
  ),
  diagnostics = diagnostics_table,
  package_versions = as.list(vapply(
    c(
      "edgeR",
      "limma",
      "variancePartition",
      "dreamlet",
      "lme4",
      "Matrix",
      "mashr"
    ),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  random_seed = NULL,
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "model_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "All ",
  dataset,
  " state models completed for ",
  length(shared_genes),
  " shared genes"
)
