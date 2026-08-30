#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
  library(SummarizedExperiment)
  library(zellkonverter)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
input_path <- file.path(
  root,
  "data_sources",
  "Chen_Cell_2021_CELLxGENE",
  "chen_discovery_epithelial.h5ad"
)
contract_path <- file.path(
  root,
  "analysis",
  "contracts",
  "state_aware_program_rederivation_v1_2026-08-29.md"
)
out_dir <- file.path(root, "results", "state_aware_program_v1", "purity_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_input_sha256 <-
  "2e2b5aa4909fd314ea0c6dba17f6ced16de95c1bb3317e69a3edfe2caff17440"
expected_contract_sha256 <-
  "0f39a03154408d14bb0bbe0c1e4f55d498e94d403bd6b45d4be5dd8918555f69"
if (!identical(
  digest::digest(input_path, algo = "sha256", file = TRUE),
  expected_input_sha256
)) {
  stop("Discovery H5AD hash changed before epithelial-purity audit")
}
if (!identical(
  digest::digest(contract_path, algo = "sha256", file = TRUE),
  expected_contract_sha256
)) {
  stop("Frozen parent contract changed before epithelial-purity audit")
}

marker_sets <- list(
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20"),
  immune = c("PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CD68", "LAPTM5"),
  stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1"),
  endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "ENG"),
  audit_targets = c("INPP5D", "TIMP1", "S100P", "LCN2")
)
all_markers <- unique(unlist(marker_sets, use.names = FALSE))

message("Reading metadata and targeted markers directly from raw/X")
sce <- readH5AD(input_path, use_hdf5 = TRUE, reader = "R", raw = FALSE)
full_metadata <- as.data.frame(colData(sce))
full_metadata$cell_index <- seq_len(nrow(full_metadata))
rm(sce)
gc(verbose = FALSE)
route_map <- c(
  NL = "normal",
  TA = "conventional_adenoma",
  TV = "conventional_adenoma",
  TVA = "conventional_adenoma"
)
full_metadata$route <- unname(route_map[as.character(full_metadata$Polyp_Type)])
full_metadata$cell_type <- as.character(full_metadata$Cell_Type)
full_metadata$specimen_id <- as.character(full_metadata$HTAN.Specimen.ID)
full_metadata$donor_id <- as.character(full_metadata$donor_id)

eligible_cells <- which(
  !is.na(full_metadata$route) &
    full_metadata$cell_type %in% c("ABS", "GOB", "TAC")
)
metadata <- full_metadata[eligible_cells, , drop = FALSE]

read_h5ad_categorical <- function(file, group) {
  codes <- as.integer(h5read(file, paste0(group, "/codes")))
  categories <- as.character(h5read(file, paste0(group, "/categories")))
  values <- rep(NA_character_, length(codes))
  valid <- codes >= 0L
  values[valid] <- categories[codes[valid] + 1L]
  values
}

feature_symbol <- read_h5ad_categorical(input_path, "/raw/var/feature_name")
available <- all_markers %in% feature_symbol
if (!all(available)) {
  warning(
    "Unavailable purity markers: ",
    paste(all_markers[!available], collapse = ", ")
  )
}
marker_names <- all_markers[available]
feature_to_marker <- match(feature_symbol, marker_names)
raw_attributes <- h5readAttributes(input_path, "/raw/X")
raw_shape <- as.integer(raw_attributes$shape)
if (!identical(raw_shape, c(nrow(full_metadata), length(feature_symbol)))) {
  stop("raw/X shape is inconsistent with metadata and raw features")
}
if (!identical(as.character(raw_attributes$`encoding-type`), "csr_matrix")) {
  stop("raw/X is not encoded as a CSR matrix")
}
indptr <- as.numeric(h5read(input_path, "/raw/X/indptr"))
if (length(indptr) != nrow(full_metadata) + 1L || indptr[1L] != 0) {
  stop("raw/X indptr is inconsistent with cell metadata")
}

cell_to_selected <- integer(nrow(full_metadata))
cell_to_selected[eligible_cells] <- seq_along(eligible_cells)
marker_counts <- Matrix(
  0,
  nrow = length(marker_names),
  ncol = length(eligible_cells),
  sparse = TRUE
)
maximum_fractional_deviation <- 0
block_size <- 1000L
for (block_start in seq.int(1L, nrow(full_metadata), by = block_size)) {
  block_end <- min(block_start + block_size - 1L, nrow(full_metadata))
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
  feature_indices <- as.integer(h5read(
    input_path,
    "/raw/X/indices",
    start = as.integer(start_pointer + 1),
    count = block_nnz
  )) + 1L
  nonzero_per_cell <- diff(indptr[block_start:(block_end + 1L)])
  cells <- rep.int(seq.int(block_start, block_end), times = nonzero_per_cell)
  deviation <- max(abs(values - round(values)))
  if (is.finite(deviation)) {
    maximum_fractional_deviation <- max(maximum_fractional_deviation, deviation)
  }
  marker_index <- feature_to_marker[feature_indices]
  selected_cell_index <- cell_to_selected[cells]
  retain <- !is.na(marker_index) & selected_cell_index > 0L
  if (any(retain)) {
    marker_counts <- marker_counts + sparseMatrix(
      i = marker_index[retain],
      j = selected_cell_index[retain],
      x = as.numeric(values[retain]),
      dims = dim(marker_counts),
      giveCsparse = TRUE
    )
  }
}
if (maximum_fractional_deviation > 1e-8) {
  stop("raw/X contains non-integer values beyond numerical tolerance")
}
rownames(marker_counts) <- marker_names
colnames(marker_counts) <- as.character(metadata$cell_index)

detected <- marker_counts > 0
module_positive <- function(markers) {
  markers <- intersect(markers, rownames(detected))
  if (!length(markers)) {
    return(rep(FALSE, ncol(detected)))
  }
  colSums(detected[markers, , drop = FALSE]) > 0
}

metadata$epithelial_marker_positive <- module_positive(marker_sets$epithelial)
metadata$immune_marker_positive <- module_positive(marker_sets$immune)
metadata$stromal_marker_positive <- module_positive(marker_sets$stromal)
metadata$endothelial_marker_positive <- module_positive(marker_sets$endothelial)
metadata$any_nonepithelial_marker_positive <-
  metadata$immune_marker_positive |
  metadata$stromal_marker_positive |
  metadata$endothelial_marker_positive

for (gene in intersect(marker_sets$audit_targets, rownames(marker_counts))) {
  metadata[[paste0(gene, "_count")]] <- marker_counts[gene, ]
  metadata[[paste0(gene, "_positive")]] <- marker_counts[gene, ] > 0
}

stratum_id <- interaction(
  metadata$specimen_id,
  metadata$cell_type,
  drop = TRUE,
  lex.order = TRUE,
  sep = "__"
)
split_indices <- split(seq_len(nrow(metadata)), stratum_id)
stratum_summary <- do.call(rbind, lapply(split_indices, function(index) {
  first <- index[1L]
  frame <- data.frame(
    specimen_id = metadata$specimen_id[first],
    donor_id = metadata$donor_id[first],
    route = metadata$route[first],
    cell_type = metadata$cell_type[first],
    n_cells = length(index),
    epithelial_positive_fraction = mean(metadata$epithelial_marker_positive[index]),
    immune_positive_fraction = mean(metadata$immune_marker_positive[index]),
    stromal_positive_fraction = mean(metadata$stromal_marker_positive[index]),
    endothelial_positive_fraction = mean(metadata$endothelial_marker_positive[index]),
    any_nonepithelial_positive_fraction = mean(
      metadata$any_nonepithelial_marker_positive[index]
    ),
    stringsAsFactors = FALSE
  )
  for (gene in intersect(marker_sets$audit_targets, rownames(marker_counts))) {
    count <- metadata[[paste0(gene, "_count")]][index]
    frame[[paste0(gene, "_positive_fraction")]] <- mean(count > 0)
    frame[[paste0(gene, "_mean_raw_count")]] <- mean(count)
  }
  frame
}))
rownames(stratum_summary) <- NULL

route_summary <- do.call(rbind, lapply(
  split(seq_len(nrow(metadata)), interaction(metadata$cell_type, metadata$route)),
  function(index) {
    data.frame(
      cell_type = metadata$cell_type[index[1L]],
      route = metadata$route[index[1L]],
      n_cells = length(index),
      n_donors = length(unique(metadata$donor_id[index])),
      epithelial_positive_fraction = mean(metadata$epithelial_marker_positive[index]),
      immune_positive_fraction = mean(metadata$immune_marker_positive[index]),
      stromal_positive_fraction = mean(metadata$stromal_marker_positive[index]),
      endothelial_positive_fraction = mean(metadata$endothelial_marker_positive[index]),
      any_nonepithelial_positive_fraction = mean(
        metadata$any_nonepithelial_marker_positive[index]
      ),
      stringsAsFactors = FALSE
    )
  }
))
rownames(route_summary) <- NULL

coexpression <- do.call(rbind, lapply(
  intersect(marker_sets$audit_targets, rownames(marker_counts)),
  function(gene) {
    positive <- metadata[[paste0(gene, "_positive")]]
    do.call(rbind, lapply(
      split(seq_len(nrow(metadata)), interaction(metadata$cell_type, metadata$route)),
      function(index) {
        gene_positive_index <- index[positive[index]]
        data.frame(
          gene = gene,
          cell_type = metadata$cell_type[index[1L]],
          route = metadata$route[index[1L]],
          n_cells = length(index),
          n_gene_positive = length(gene_positive_index),
          gene_positive_fraction = mean(positive[index]),
          fraction_gene_positive_with_epithelial_marker = if (
            length(gene_positive_index)
          ) {
            mean(metadata$epithelial_marker_positive[gene_positive_index])
          } else {
            NA_real_
          },
          fraction_gene_positive_with_immune_marker = if (
            length(gene_positive_index)
          ) {
            mean(metadata$immune_marker_positive[gene_positive_index])
          } else {
            NA_real_
          },
          fraction_gene_positive_with_any_nonepithelial_marker = if (
            length(gene_positive_index)
          ) {
            mean(metadata$any_nonepithelial_marker_positive[gene_positive_index])
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      }
    ))
  }
))
rownames(coexpression) <- NULL

marker_detection <- do.call(rbind, lapply(
  rownames(marker_counts),
  function(gene) {
    do.call(rbind, lapply(
      split(seq_len(nrow(metadata)), interaction(metadata$cell_type, metadata$route)),
      function(index) {
        counts <- marker_counts[gene, index]
        data.frame(
          gene = gene,
          cell_type = metadata$cell_type[index[1L]],
          route = metadata$route[index[1L]],
          n_cells = length(index),
          positive_fraction = mean(counts > 0),
          mean_raw_count = mean(counts),
          stringsAsFactors = FALSE
        )
      }
    ))
  }
))
rownames(marker_detection) <- NULL

write.table(
  stratum_summary,
  file.path(out_dir, "specimen_state_purity_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  route_summary,
  file.path(out_dir, "cell_state_route_purity_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  coexpression,
  file.path(out_dir, "audit_target_cellular_coexpression.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  marker_detection,
  file.path(out_dir, "marker_detection_by_state_route.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_audit_epithelial_purity_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_sha256 = expected_input_sha256,
  contract_sha256 = expected_contract_sha256,
  cells_audited = nrow(metadata),
  states = c("ABS", "GOB", "TAC"),
  routes = c("normal", "conventional_adenoma"),
  marker_sets = marker_sets,
  assay = "raw/X integer counts",
  scope = "targeted contamination and co-expression audit; not a decontamination algorithm",
  package_versions = as.list(vapply(
    c("Matrix", "SummarizedExperiment", "zellkonverter", "digest", "jsonlite"),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "purity_audit_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Purity audit complete for ",
  nrow(metadata),
  " eligible epithelial cells"
)
