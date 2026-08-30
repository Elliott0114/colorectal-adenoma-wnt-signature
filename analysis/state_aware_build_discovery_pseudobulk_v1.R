#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(rhdf5)
  library(SummarizedExperiment)
  library(zellkonverter)
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
  "data_sources",
  "Chen_Cell_2021_CELLxGENE",
  paste0("chen_", dataset, "_epithelial.h5ad")
)
contract_path <- file.path(
  root,
  "analysis",
  "contracts",
  "state_aware_program_rederivation_v1_2026-08-29.md"
)
out_dir <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  paste0(dataset, "_pseudobulk")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_input_sha256 <- c(
  discovery = "2e2b5aa4909fd314ea0c6dba17f6ced16de95c1bb3317e69a3edfe2caff17440",
  validation = "f51bf02885d637073743093844c1db00a173196f0e1184b37f39a5118468668d"
)[[dataset]]
actual_input_sha256 <- digest::digest(input_path, algo = "sha256", file = TRUE)
if (!identical(actual_input_sha256, expected_input_sha256)) {
  stop("Discovery H5AD SHA-256 does not match the frozen contract")
}

route_map <- c(
  NL = "normal",
  TA = "conventional_adenoma",
  TV = "conventional_adenoma",
  TVA = "conventional_adenoma"
)
canonical_states <- c("ABS", "GOB", "TAC", "STM", "CT", "TUF", "EE")
minimum_cells <- 20L
minimum_donors_per_route <- 10L
minimum_paired_donors <- 6L
block_size <- 500L

message("Reading ", dataset, " metadata with the native R H5AD reader")
sce <- readH5AD(input_path, use_hdf5 = TRUE, reader = "R", raw = FALSE)
obs <- as.data.frame(colData(sce)) %>%
  dplyr::transmute(
    cell_index = seq_len(n()),
    specimen_id = as.character(HTAN.Specimen.ID),
    donor_id = as.character(donor_id),
    polyp_type = as.character(Polyp_Type),
    sample_classification = as.character(Sample_Classification),
    cell_type = as.character(Cell_Type),
    route = unname(route_map[polyp_type])
  )
rm(sce)
gc(verbose = FALSE)

if (anyNA(obs$cell_index) || anyDuplicated(obs$cell_index)) {
  stop("Cell indices are missing or duplicated")
}
if (anyNA(obs$specimen_id) || anyNA(obs$donor_id)) {
  stop("Specimen or donor identifiers are missing")
}

specimen_metadata_audit <- obs %>%
  dplyr::distinct(specimen_id, donor_id, polyp_type, route) %>%
  dplyr::count(specimen_id, name = "n_distinct_metadata_rows")
if (any(specimen_metadata_audit$n_distinct_metadata_rows != 1L)) {
  stop("At least one specimen maps to multiple donor or route records")
}

state_strata <- obs %>%
  dplyr::filter(!is.na(route)) %>%
  dplyr::count(specimen_id, donor_id, route, cell_type, name = "n_cells")

state_support <- state_strata %>%
  dplyr::filter(cell_type %in% canonical_states, n_cells >= minimum_cells) %>%
  dplyr::distinct(cell_type, route, donor_id) %>%
  dplyr::count(cell_type, route, name = "n_donors") %>%
  tidyr::pivot_wider(
    names_from = route,
    values_from = n_donors,
    values_fill = 0L
  )

paired_support <- state_strata %>%
  dplyr::filter(cell_type %in% canonical_states, n_cells >= minimum_cells) %>%
  dplyr::distinct(cell_type, donor_id, route) %>%
  dplyr::count(cell_type, donor_id, name = "n_routes") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(n_paired_donors = sum(n_routes == 2L), .groups = "drop")

state_eligibility <- tibble::tibble(cell_type = canonical_states) %>%
  dplyr::left_join(state_support, by = "cell_type") %>%
  dplyr::left_join(paired_support, by = "cell_type") %>%
  dplyr::mutate(
    dplyr::across(c(normal, conventional_adenoma, n_paired_donors), ~ tidyr::replace_na(.x, 0L)),
    eligible = normal >= minimum_donors_per_route &
      conventional_adenoma >= minimum_donors_per_route &
      n_paired_donors >= minimum_paired_donors
  )

if (dataset == "discovery") {
  eligible_states <- state_eligibility %>%
    dplyr::filter(eligible) %>%
    dplyr::pull(cell_type)
  if (!identical(sort(eligible_states), sort(c("ABS", "GOB", "TAC")))) {
    stop(
      "Sampling-only state eligibility did not reproduce the frozen ABS/GOB/TAC set: ",
      paste(eligible_states, collapse = ",")
    )
  }
} else {
  eligible_states <- c("ABS", "GOB", "TAC")
}

eligible_strata <- state_strata %>%
  dplyr::filter(cell_type %in% eligible_states, n_cells >= minimum_cells) %>%
  dplyr::arrange(cell_type, donor_id, route, specimen_id) %>%
  dplyr::mutate(
    pseudobulk_id = paste(specimen_id, cell_type, sep = "__"),
    pseudobulk_index = row_number()
  )

if (anyDuplicated(eligible_strata$pseudobulk_id)) {
  stop("Pseudobulk identifiers are not unique")
}

obs_selected <- obs %>%
  dplyr::inner_join(
    eligible_strata %>% dplyr::select(specimen_id, cell_type, pseudobulk_index),
    by = c("specimen_id", "cell_type")
  )

if (nrow(obs_selected) != sum(eligible_strata$n_cells)) {
  stop("Selected-cell count does not equal the eligible-stratum inventory")
}

read_h5ad_categorical <- function(file, group) {
  codes <- as.integer(h5read(file, paste0(group, "/codes")))
  categories <- as.character(h5read(file, paste0(group, "/categories")))
  values <- rep(NA_character_, length(codes))
  valid <- codes >= 0L
  values[valid] <- categories[codes[valid] + 1L]
  values
}

feature_id <- as.character(h5read(input_path, "/raw/var/_index"))
feature_symbol <- read_h5ad_categorical(input_path, "/raw/var/feature_name")
if (length(feature_id) != length(feature_symbol)) {
  stop("Raw feature identifiers and symbols have different lengths")
}

valid_feature <- !is.na(feature_symbol) & nzchar(feature_symbol)
unique_symbols <- unique(feature_symbol[valid_feature])
symbol_index <- integer(length(feature_symbol))
symbol_index[valid_feature] <- match(feature_symbol[valid_feature], unique_symbols)
feature_map <- tibble::tibble(
  raw_feature_index = seq_along(feature_id),
  feature_id = feature_id,
  gene = feature_symbol,
  valid_symbol = valid_feature,
  merged_gene_index = symbol_index,
  symbol_multiplicity = ave(
    as.integer(valid_feature),
    feature_symbol,
    FUN = function(x) sum(x > 0L)
  )
)

raw_attributes <- h5readAttributes(input_path, "/raw/X")
raw_shape <- as.integer(raw_attributes$shape)
if (!identical(raw_shape, c(nrow(obs), length(feature_id)))) {
  stop(
    "Raw matrix shape does not match metadata/features: ",
    paste(raw_shape, collapse = "x")
  )
}
if (!identical(as.character(raw_attributes$`encoding-type`), "csr_matrix")) {
  stop("raw/X is not encoded as a CSR matrix")
}

indptr <- as.numeric(h5read(input_path, "/raw/X/indptr"))
if (length(indptr) != nrow(obs) + 1L || indptr[1L] != 0) {
  stop("raw/X indptr is inconsistent with the cell count")
}

cell_to_pseudobulk <- integer(nrow(obs))
cell_to_pseudobulk[obs_selected$cell_index] <- obs_selected$pseudobulk_index

counts <- Matrix(
  0,
  nrow = length(unique_symbols),
  ncol = nrow(eligible_strata),
  sparse = TRUE
)
maximum_fractional_deviation <- 0
processed_nonzero <- 0
retained_nonzero <- 0

message(
  "Aggregating raw CSR counts across ",
  nrow(obs),
  " cells in blocks of ",
  block_size
)
for (block_start in seq.int(1L, nrow(obs), by = block_size)) {
  block_end <- min(block_start + block_size - 1L, nrow(obs))
  start_pointer <- indptr[block_start]
  end_pointer <- indptr[block_end + 1L]
  block_nnz <- as.integer(end_pointer - start_pointer)
  if (block_nnz == 0L) {
    next
  }

  values <- h5read(
    input_path,
    "/raw/X/data",
    start = as.integer(start_pointer + 1),
    count = block_nnz
  )
  indices <- as.integer(h5read(
    input_path,
    "/raw/X/indices",
    start = as.integer(start_pointer + 1),
    count = block_nnz
  )) + 1L
  nonzero_per_cell <- diff(indptr[block_start:(block_end + 1L)])
  cells <- rep.int(seq.int(block_start, block_end), times = nonzero_per_cell)

  if (length(values) != length(indices) || length(values) != length(cells)) {
    stop("CSR block components have inconsistent lengths")
  }

  deviation <- max(abs(values - round(values)))
  if (is.finite(deviation)) {
    maximum_fractional_deviation <- max(maximum_fractional_deviation, deviation)
  }

  pseudobulk_index <- cell_to_pseudobulk[cells]
  gene_index <- symbol_index[indices]
  retain <- pseudobulk_index > 0L & gene_index > 0L
  processed_nonzero <- processed_nonzero + length(values)
  retained_nonzero <- retained_nonzero + sum(retain)

  if (any(retain)) {
    counts <- counts + sparseMatrix(
      i = gene_index[retain],
      j = pseudobulk_index[retain],
      x = as.numeric(values[retain]),
      dims = dim(counts),
      giveCsparse = TRUE
    )
  }
}

if (maximum_fractional_deviation > 1e-8) {
  stop("raw/X contains non-integer values beyond numerical tolerance")
}

counts <- drop0(counts)
rownames(counts) <- unique_symbols
colnames(counts) <- eligible_strata$pseudobulk_id

if (any(colSums(counts) <= 0)) {
  stop("At least one eligible pseudobulk has zero retained counts")
}
if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
  stop("Pseudobulk matrix dimnames are not unique")
}

pseudobulk_metadata <- eligible_strata %>%
  dplyr::transmute(
    pseudobulk_id,
    specimen_id,
    donor_id,
    route = factor(route, levels = c("normal", "conventional_adenoma")),
    cell_type = factor(cell_type, levels = c("ABS", "GOB", "TAC")),
    n_cells,
    library_size = as.numeric(colSums(counts))
  )

if (!identical(as.character(pseudobulk_metadata$pseudobulk_id), colnames(counts))) {
  stop("Pseudobulk metadata and count columns are misaligned")
}

state_inventory <- state_strata %>%
  dplyr::group_by(cell_type, route) %>%
  dplyr::summarise(
    n_cells = sum(n_cells),
    n_specimens = n_distinct(specimen_id),
    n_donors = n_distinct(donor_id),
    .groups = "drop"
  )

input_contract <- list(
  counts = counts,
  metadata = pseudobulk_metadata,
  state_eligibility = state_eligibility,
  minimum_cells = minimum_cells,
  minimum_donors_per_route = minimum_donors_per_route,
  minimum_paired_donors = minimum_paired_donors,
  input_sha256 = actual_input_sha256,
  contract_sha256 = digest::digest(contract_path, algo = "sha256", file = TRUE)
)
saveRDS(
  input_contract,
  file.path(out_dir, paste0(dataset, "_state_pseudobulk.rds")),
  compress = "xz"
)

write.table(
  pseudobulk_metadata,
  file.path(out_dir, "pseudobulk_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  state_inventory,
  file.path(out_dir, "state_inventory_all_routes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  state_eligibility,
  file.path(out_dir, "state_eligibility.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  feature_map,
  file.path(out_dir, "raw_feature_symbol_map.tsv.gz"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = paste0("state_aware_build_", dataset, "_pseudobulk_v1"),
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input = list(
    path = normalizePath(input_path),
    sha256 = actual_input_sha256,
    raw_shape_cells_by_features = raw_shape,
    raw_encoding = as.character(raw_attributes$`encoding-type`)
  ),
  contract = list(
    path = normalizePath(contract_path),
    sha256 = input_contract$contract_sha256
  ),
  eligibility = list(
    minimum_cells_per_specimen_state = minimum_cells,
    minimum_donors_per_route = minimum_donors_per_route,
    minimum_paired_donors = minimum_paired_donors,
    retained_states = eligible_states,
    state_rule = if (dataset == "discovery") {
      "sampling-only discovery eligibility"
    } else {
      "states frozen from discovery; validation support reported without reselection"
    }
  ),
  output = list(
    genes_after_symbol_merge = nrow(counts),
    pseudobulks = ncol(counts),
    selected_cells = nrow(obs_selected),
    processed_nonzero_entries = processed_nonzero,
    retained_nonzero_entries = retained_nonzero,
    maximum_fractional_deviation = maximum_fractional_deviation,
    duplicate_symbols_merged = sum(table(feature_symbol[valid_feature]) > 1L)
  ),
  random_seed = NULL,
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "pseudobulk_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Completed: ",
  nrow(counts),
  " genes x ",
  ncol(counts),
  " pseudobulks from ",
  nrow(obs_selected),
  " cells"
)
