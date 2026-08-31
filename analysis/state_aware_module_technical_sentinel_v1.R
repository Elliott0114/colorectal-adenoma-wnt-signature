#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
  library(nlme)
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
set.seed(20260830)
allowWGCNAThreads(nThreads = 8)

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
source_root <- file.path(result_root, "functional_architecture_v1")
run_root_env <- Sys.getenv("STATE_AWARE_MODULE_RUN_ROOT", unset = "")
out_root <- if (nzchar(run_root_env)) {
  normalizePath(run_root_env, mustWork = FALSE)
} else {
  source_root
}
wgcna_root <- file.path(source_root, "consensus_wgcna")
external_root <- file.path(out_root, "module_external_validation")
technical_root <- file.path(out_root, "module_technical_validation")
audit_root <- file.path(technical_root, "audit")
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  modules = file.path(source_root, "consensus_modules.tsv"),
  modules_output = file.path(out_root, "consensus_modules.tsv"),
  validation = file.path(out_root, "module_validation.tsv"),
  projection = file.path(wgcna_root, "network_projection_object.rds"),
  discovery_pseudobulk = file.path(
    result_root, "discovery_pseudobulk", "discovery_state_pseudobulk.rds"
  ),
  heldout_pseudobulk = file.path(
    result_root, "validation_pseudobulk", "validation_state_pseudobulk.rds"
  ),
  measurability = file.path(
    external_root, "module_gene_platform_measurability.tsv"
  ),
  gene_info = file.path(
    root, "data_sources", "GSE117606", "Homo_sapiens.gene_info.gz"
  ),
  addendum = file.path(
    root, "analysis", "contracts",
    "state_aware_functional_architecture_downstream_addendum_v1_2026-08-30.md"
  )
)
base_inputs <- paths[setdiff(
  names(paths), c("modules_output", "measurability", "gene_info")
)]
if (!all(file.exists(unlist(base_inputs)))) {
  stop("At least one core technical/sentinel input is missing")
}
sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)
if (sha256(paths$addendum) !=
    "32343afe117d09007066fbe01f8fbe7cf4a11ee628dbe6a121f0f814968da3bb") {
  stop("The frozen downstream validation addendum changed")
}

modules <- read.delim(paths$modules, check.names = FALSE)
validation <- read.delim(paths$validation, check.names = FALSE)
route_column <- Sys.getenv(
  "STATE_AWARE_MODULE_ROUTE_COLUMN", unset = "internal_gate_pass"
)
if (!route_column %in% colnames(validation)) {
  stop("Routing column is missing: ", route_column)
}
as_bool <- function(value) {
  if (is.logical(value)) value else tolower(as.character(value)) == "true"
}
route_pass <- as_bool(validation[[route_column]])
descriptive_downstream <- tolower(Sys.getenv(
  "STATE_AWARE_MODULE_DESCRIPTIVE_DOWNSTREAM", unset = "false"
)) == "true"
retained_modules <- validation$module[
  route_pass & (descriptive_downstream | validation$external_gate_pass)
]

empty_output <- function() {
  write.table(
    data.frame(), file.path(technical_root, "module_technical_covariate_audit.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(), file.path(technical_root, "sentinel_candidate_base.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  validation$technical_gate_pass <- FALSE
  validation$routing_status <- ifelse(
    route_pass & validation$external_gate_pass,
    "audit_technical_gate_not_evaluable",
    validation$routing_status
  )
  write.table(
    validation, paths$validation,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
if (!length(retained_modules)) {
  empty_output()
  message("No module entered technical/sentinel analysis under the active routing rule")
  quit(save = "no", status = 0)
}
if (!file.exists(paths$measurability) || !file.exists(paths$gene_info)) {
  stop("Platform measurability or gene annotation is missing for retained modules")
}

projection <- readRDS(paths$projection)
module_labels <- projection$module_labels
retained_genes <- modules$gene[modules$module %in% retained_modules]
retained_labels <- module_labels[retained_genes]

orient_eigengene <- function(expression, genes, module) {
  local_genes <- intersect(genes, colnames(expression))
  if (length(local_genes) < 2L) {
    return(rep(NA_real_, nrow(expression)))
  }
  local_expression <- expression[, local_genes, drop = FALSE]
  eigengene <- moduleEigengenes(
    local_expression,
    colors = rep(module, length(local_genes)),
    excludeGrey = TRUE, align = "along average"
  )$eigengenes[[1L]]
  average <- rowMeans(scale(local_expression), na.rm = TRUE)
  if (is.finite(cor(eigengene, average, use = "pairwise.complete.obs")) &&
      cor(eigengene, average, use = "pairwise.complete.obs") < 0) {
    eigengene <- -eigengene
  }
  eigengene
}

technical_records <- list()
profile_eigengenes <- list()
technical_counter <- 0L
for (state in c("ABS", "GOB", "TAC")) {
  expression <- projection$discovery$states[[state]]$profile_expression
  metadata <- projection$discovery$states[[state]]$profile_metadata
  stopifnot(identical(rownames(expression), metadata$profile_id))
  metadata$log10_library_size <- log10(metadata$library_size_input + 1)
  metadata$log10_cell_count <- log10(metadata$n_cells + 1)
  metadata$is_adenoma <- as.integer(metadata$route == "conventional_adenoma")
  for (module in retained_modules) {
    genes <- modules$gene[modules$module == module]
    eigengene <- orient_eigengene(expression, genes, module)
    metadata$module_eigengene <- eigengene
    rho_library <- suppressWarnings(cor(
      eigengene, metadata$log10_library_size,
      method = "spearman", use = "pairwise.complete.obs"
    ))
    rho_cells <- suppressWarnings(cor(
      eigengene, metadata$log10_cell_count,
      method = "spearman", use = "pairwise.complete.obs"
    ))
    fit <- tryCatch(
      lme(
        module_eigengene ~ is_adenoma + log10_library_size + log10_cell_count,
        random = ~1 | donor_id, data = metadata, method = "REML",
        control = lmeControl(returnObject = TRUE)
      ),
      error = function(error) NULL
    )
    coefficient <- p_value <- NA_real_
    if (!is.null(fit)) {
      coefficient_table <- summary(fit)$tTable
      coefficient <- coefficient_table["is_adenoma", "Value"]
      p_value <- coefficient_table["is_adenoma", "p-value"]
    }
    expected_sign <- ifelse(
      validation$heldout_direction[match(module, validation$module)] == "Up",
      1, -1
    )
    technical_counter <- technical_counter + 1L
    technical_records[[technical_counter]] <- data.frame(
      module = module, state = state,
      n_profiles = nrow(metadata), n_donors = length(unique(metadata$donor_id)),
      spearman_library_size = rho_library,
      spearman_cell_count = rho_cells,
      adjusted_route_coefficient = coefficient,
      expected_route_sign = expected_sign,
      adjusted_route_direction_concordant = is.finite(coefficient) &&
        sign(coefficient) == expected_sign,
      adjusted_route_p_value = p_value,
      stringsAsFactors = FALSE
    )
    profile_eigengenes[[paste(state, module, sep = "__")]] <- data.frame(
      module = module, state = state,
      donor_id = metadata$donor_id, route = metadata$route,
      n_cells = metadata$n_cells,
      library_size_input = metadata$library_size_input,
      module_eigengene = eigengene,
      stringsAsFactors = FALSE
    )
  }
}
technical_table <- do.call(rbind, technical_records)
technical_gate <- do.call(rbind, lapply(retained_modules, function(module) {
  local <- technical_table[technical_table$module == module, ]
  high_technical_correlation <- any(
    abs(local$spearman_library_size) >= 0.50 |
      abs(local$spearman_cell_count) >= 0.50,
    na.rm = TRUE
  )
  no_adjusted_support <- all(
    !local$adjusted_route_direction_concordant |
      !is.finite(local$adjusted_route_p_value) |
      local$adjusted_route_p_value > 0.10
  )
  data.frame(
    module = module,
    high_technical_correlation = high_technical_correlation,
    no_adjusted_state_support = no_adjusted_support,
    technical_gate_pass = !(high_technical_correlation && no_adjusted_support),
    stringsAsFactors = FALSE
  )
}))

message("Bootstrapping kME uncertainty for observed hub candidates")
observed_candidate <- do.call(rbind, lapply(retained_modules, function(module) {
  local <- modules[modules$module == module, ]
  observed_passes <- rowSums(cbind(
    local$abs_kme_discovery_ABS >= 0.70,
    local$abs_kme_discovery_GOB >= 0.70,
    local$abs_kme_discovery_TAC >= 0.70
  ), na.rm = TRUE)
  local$observed_kme_states_0.70 <- observed_passes
  local[observed_passes >= 2L, , drop = FALSE]
}))

bootstrap_records <- list()
bootstrap_counter <- 0L
for (state in c("ABS", "GOB", "TAC")) {
  expression <- projection$discovery$states[[state]]$residual_donor_expression
  for (module in retained_modules) {
    genes <- observed_candidate$gene[observed_candidate$module == module]
    genes <- intersect(genes, colnames(expression))
    if (!length(genes)) next
    eigengene <- projection$discovery_kme[[state]]$eigengenes[[
      paste0("ME", modules$module_color[
        match(module, modules$module)
      ])
    ]]
    if (is.null(eigengene)) {
      eigengene <- orient_eigengene(
        expression,
        modules$gene[modules$module == module],
        module
      )
    }
    if (cor(eigengene, rowMeans(scale(expression[, genes, drop = FALSE])),
            use = "pairwise.complete.obs") < 0) {
      eigengene <- -eigengene
    }
    n_donors <- nrow(expression)
    bootstrap_matrix <- matrix(
      NA_real_, nrow = length(genes), ncol = 1000L,
      dimnames = list(genes, NULL)
    )
    for (bootstrap_index in seq_len(1000L)) {
      indices <- sample.int(n_donors, n_donors, replace = TRUE)
      correlation <- bicor(
        expression[indices, genes, drop = FALSE],
        eigengene[indices],
        use = "p", maxPOutliers = 0.05
      )
      bootstrap_matrix[, bootstrap_index] <- abs(as.numeric(correlation))
    }
    lower <- apply(bootstrap_matrix, 1L, quantile, probs = 0.025, na.rm = TRUE)
    median <- apply(bootstrap_matrix, 1L, median, na.rm = TRUE)
    upper <- apply(bootstrap_matrix, 1L, quantile, probs = 0.975, na.rm = TRUE)
    for (gene in genes) {
      bootstrap_counter <- bootstrap_counter + 1L
      bootstrap_records[[bootstrap_counter]] <- data.frame(
        gene = gene, module = module, state = state,
        bootstrap_absolute_kme_median = median[gene],
        bootstrap_absolute_kme_ci_low = lower[gene],
        bootstrap_absolute_kme_ci_high = upper[gene],
        bootstrap_ci_low_gt_0.40 = lower[gene] > 0.40,
        n_bootstraps = 1000L, seed = 20260830,
        stringsAsFactors = FALSE
      )
    }
  }
}
bootstrap_kme <- if (length(bootstrap_records)) {
  do.call(rbind, bootstrap_records)
} else {
  data.frame()
}
for (state in c("ABS", "GOB", "TAC")) {
  if (nrow(bootstrap_kme)) {
    local <- bootstrap_kme[bootstrap_kme$state == state, ]
    key <- paste(local$gene, local$module, sep = "__")
    module_key <- paste(modules$gene, modules$module, sep = "__")
    modules[[paste0("bootstrap_abs_kme_ci_low_", state)]] <-
      local$bootstrap_absolute_kme_ci_low[match(module_key, key)]
    modules[[paste0("bootstrap_abs_kme_median_", state)]] <-
      local$bootstrap_absolute_kme_median[match(module_key, key)]
    modules[[paste0("bootstrap_abs_kme_ci_high_", state)]] <-
      local$bootstrap_absolute_kme_ci_high[match(module_key, key)]
  } else {
    modules[[paste0("bootstrap_abs_kme_ci_low_", state)]] <- NA_real_
    modules[[paste0("bootstrap_abs_kme_median_", state)]] <- NA_real_
    modules[[paste0("bootstrap_abs_kme_ci_high_", state)]] <- NA_real_
  }
}

discovery_object <- readRDS(paths$discovery_pseudobulk)
heldout_object <- readRDS(paths$heldout_pseudobulk)
discovery_metadata <- as.data.frame(discovery_object$metadata)
heldout_metadata <- as.data.frame(heldout_object$metadata)
overlap_donors <- intersect(
  unique(as.character(discovery_metadata$donor_id)),
  unique(as.character(heldout_metadata$donor_id))
)
heldout_keep <- !heldout_metadata$donor_id %in% overlap_donors
state_keep_discovery <- discovery_metadata$cell_type %in% c("ABS", "GOB", "TAC")
state_keep_heldout <- heldout_keep &
  heldout_metadata$cell_type %in% c("ABS", "GOB", "TAC")
detection <- data.frame(gene = observed_candidate$gene, stringsAsFactors = FALSE)
detection$discovery_detection_fraction <- vapply(detection$gene, function(gene) {
  mean(discovery_object$counts[gene, state_keep_discovery] > 0)
}, numeric(1))
detection$heldout_detection_fraction <- vapply(detection$gene, function(gene) {
  if (!gene %in% rownames(heldout_object$counts)) return(0)
  mean(heldout_object$counts[gene, state_keep_heldout] > 0)
}, numeric(1))

measurability <- read.delim(paths$measurability, check.names = FALSE)
measurability$measurable <- as.logical(measurability$measurable)
platform_summary <- do.call(rbind, lapply(observed_candidate$gene, function(gene) {
  local <- measurability[measurability$gene == gene, ]
  data.frame(
    gene = gene,
    external_cohorts_measurable = sum(
      local$measurable & local$platform != "GSE117606", na.rm = TRUE
    ),
    ffpe_measurable = any(
      local$measurable & local$platform == "GSE117606", na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}))

gene_info <- read.delim(
  gzfile(paths$gene_info), check.names = FALSE, quote = "", comment.char = ""
)
gene_type_column <- if ("type_of_gene" %in% colnames(gene_info)) {
  "type_of_gene"
} else {
  "type.of.gene"
}
protein_coding <- unique(as.character(gene_info$Symbol[
  gene_info[[gene_type_column]] == "protein-coding"
]))

bootstrap_summary <- if (nrow(bootstrap_kme)) {
  do.call(rbind, lapply(unique(bootstrap_kme$gene), function(gene) {
    local <- bootstrap_kme[bootstrap_kme$gene == gene, ]
    data.frame(
      gene = gene,
      bootstrap_kme_states_ci_low_gt_0.40 = sum(
        local$bootstrap_ci_low_gt_0.40, na.rm = TRUE
      ),
      minimum_bootstrap_kme_ci_low = min(
        local$bootstrap_absolute_kme_ci_low, na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame(gene = character(), stringsAsFactors = FALSE)
}

candidate_base <- merge(
  observed_candidate, bootstrap_summary, by = "gene", all.x = TRUE,
  suffixes = c("", ".bootstrap")
)
candidate_base <- merge(candidate_base, detection, by = "gene", all.x = TRUE)
candidate_base <- merge(candidate_base, platform_summary, by = "gene", all.x = TRUE)
candidate_base$protein_coding <- candidate_base$gene %in% protein_coding
candidate_base$heldout_gene_direction_concordant <- ifelse(
  validation$heldout_direction[match(candidate_base$module, validation$module)] == "Up",
  candidate_base$heldout_common_z > 0,
  candidate_base$heldout_common_z < 0
)
candidate_base$kme_gate_pass <-
  candidate_base$bootstrap_kme_states_ci_low_gt_0.40 >= 2L
candidate_base$detection_gate_pass <-
  candidate_base$discovery_detection_fraction >= 0.80 &
  candidate_base$heldout_detection_fraction >= 0.80
candidate_base$heldout_gene_gate_pass <-
  candidate_base$heldout_gene_direction_concordant &
  is.finite(candidate_base$heldout_common_fdr) &
  candidate_base$heldout_common_fdr <= 0.10
candidate_base$platform_gate_pass <-
  candidate_base$external_cohorts_measurable >= 4L &
  candidate_base$ffpe_measurable
candidate_base$pre_orthogonal_gate_pass <-
  candidate_base$protein_coding & candidate_base$kme_gate_pass &
  candidate_base$detection_gate_pass & candidate_base$heldout_gene_gate_pass &
  candidate_base$platform_gate_pass

technical_gate_columns <- c(
  "high_technical_correlation",
  "no_adjusted_state_support",
  "technical_gate_pass"
)
validation <- validation[, setdiff(
  colnames(validation), technical_gate_columns
), drop = FALSE]
validation <- merge(validation, technical_gate, by = "module", all.x = TRUE)
validation$technical_gate_pass[is.na(validation$technical_gate_pass)] <- FALSE
route_pass <- as_bool(validation[[route_column]])
validation$routing_status <- ifelse(
  route_pass & validation$external_gate_pass &
    !validation$technical_gate_pass,
  "audit_technical_gate_failed",
  validation$routing_status
)

write.table(
  technical_table, file.path(technical_root, "module_technical_covariate_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  technical_gate, file.path(technical_root, "module_technical_gate.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, profile_eigengenes),
  file.path(audit_root, "module_profile_eigengenes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  bootstrap_kme, file.path(technical_root, "sentinel_kme_bootstrap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  candidate_base, file.path(technical_root, "sentinel_candidate_base.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  modules, paths$modules_output,
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  validation, paths$validation,
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_module_technical_sentinel_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = 20260830,
  route_column = route_column,
  descriptive_downstream = descriptive_downstream,
  retained_modules_entering_technical_gate = retained_modules,
  modules_passing_technical_gate = technical_gate$module[
    technical_gate$technical_gate_pass
  ],
  pre_orthogonal_sentinel_candidates = candidate_base$gene[
    candidate_base$pre_orthogonal_gate_pass
  ],
  input_sha256 = as.list(vapply(
    paths[setdiff(names(paths), "modules_output")], sha256, character(1)
  )),
  output_sha256 = list(consensus_modules = sha256(paths$modules_output)),
  package_versions = list(
    R = as.character(getRversion()),
    WGCNA = as.character(packageVersion("WGCNA")),
    nlme = as.character(packageVersion("nlme")),
    digest = as.character(packageVersion("digest")),
    jsonlite = as.character(packageVersion("jsonlite"))
  ),
  session_info = capture.output(sessionInfo())
)
write_json(
  manifest, file.path(technical_root, "module_technical_sentinel_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)
message(
  "Technical/sentinel base analysis complete: ",
  sum(technical_gate$technical_gate_pass), " modules passed technical gate; ",
  sum(candidate_base$pre_orthogonal_gate_pass),
  " candidates await orthogonal evidence"
)
