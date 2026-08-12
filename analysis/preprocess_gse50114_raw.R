#!/usr/bin/env Rscript

# Label-free preprocessing of complete GSE50114 Agilent one-colour raw arrays.

suppressPackageStartupMessages({
  library(limma)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "analysis", "preprocess_gse50114_raw.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
SOURCE_TAR <- file.path(ROOT, "data_sources", "GEO_sporadic_adenoma_validation", "GSE50114_RAW.tar")
OUT_DIR <- file.path(ROOT, "results", "external_sporadic_adenoma_validation")
OUT_MATRIX <- file.path(OUT_DIR, "GSE50114_raw_limma_log2_quantile.tsv.gz")
OUT_QC <- file.path(OUT_DIR, "GSE50114_raw_preprocessing_qc.tsv")
OUT_MANIFEST <- file.path(OUT_DIR, "GSE50114_raw_preprocessing_manifest.tsv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(SOURCE_TAR))

work_dir <- tempfile("gse50114_raw_")
dir.create(work_dir, recursive = TRUE)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

utils::untar(SOURCE_TAR, exdir = work_dir)
files <- sort(list.files(work_dir, pattern = "^GSM.*\\.txt\\.gz$", full.names = TRUE))
stopifnot(length(files) == 46)

raw <- read.maimages(files, source = "agilent", green.only = TRUE, verbose = FALSE)
sample_ids <- sub("_.*$", "", basename(files))
stopifnot(length(unique(sample_ids)) == 46)
colnames(raw$E) <- sample_ids

corrected <- backgroundCorrect(raw, method = "normexp", offset = 50, verbose = FALSE)
normalized <- normalizeBetweenArrays(corrected, method = "quantile")

keep <- normalized$genes$ControlType == 0 &
  !is.na(normalized$genes$ProbeName) &
  nzchar(normalized$genes$ProbeName)
probe_ids <- as.character(normalized$genes$ProbeName[keep])
matrix <- avereps(normalized$E[keep, , drop = FALSE], ID = probe_ids)
probe_ids <- rownames(matrix)
stopifnot(!anyDuplicated(probe_ids))
stopifnot(all(is.finite(matrix)))

output <- data.frame(feature_id = rownames(matrix), matrix, check.names = FALSE)
con <- gzfile(OUT_MATRIX, open = "wt", encoding = "UTF-8")
write.table(output, con, sep = "\t", row.names = FALSE, quote = FALSE)
close(con)

sample_qc <- data.frame(
  sample_id = colnames(matrix),
  n_noncontrol_probes = nrow(matrix),
  minimum = apply(matrix, 2, min),
  q25 = apply(matrix, 2, quantile, probs = 0.25),
  median = apply(matrix, 2, median),
  q75 = apply(matrix, 2, quantile, probs = 0.75),
  maximum = apply(matrix, 2, max),
  stringsAsFactors = FALSE
)
write.table(sample_qc, OUT_QC, sep = "\t", row.names = FALSE, quote = FALSE)

manifest <- data.frame(
  source_file = "data_sources/GEO_sporadic_adenoma_validation/GSE50114_RAW.tar",
  n_raw_files = length(files),
  n_noncontrol_probes = nrow(matrix),
  n_samples = ncol(matrix),
  background_correction = "normexp",
  background_offset = 50,
  between_array_normalization = "quantile",
  group_labels_used = FALSE,
  output_file = "results/external_sporadic_adenoma_validation/GSE50114_raw_limma_log2_quantile.tsv.gz",
  stringsAsFactors = FALSE
)
write.table(manifest, OUT_MANIFEST, sep = "\t", row.names = FALSE, quote = FALSE)

message("GSE50114 raw preprocessing complete: ", nrow(matrix), " probes x ", ncol(matrix), " samples")
