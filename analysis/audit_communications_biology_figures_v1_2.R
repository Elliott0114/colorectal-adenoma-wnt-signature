#!/usr/bin/env Rscript

# Structural and text-level audit for the canonical Communications Biology
# figure package. Visual inspection is recorded separately in the audit report.

options(stringsAsFactors = FALSE)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  file.path(getwd(), "analysis", "audit_communications_biology_figures_v1_2.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
FIG_DIR <- file.path(ROOT, "figures", "communications_biology_v1.2")
SOURCE_DIR <- file.path(FIG_DIR, "source_data")

expected <- data.frame(
  figure = c(
    paste0("figure", 1:6),
    paste0("figureS", 1:8)
  ),
  stem = c(
    "figure1_discovery_core_and_objective_reduction",
    "figure2_independent_replication_and_ffpe",
    "figure3_rna_atac_regulatory_support",
    "figure4_crc_atlas_cross_sectional_recurrence",
    "figure5_empirical_and_virtual_perturbation_support",
    "figure6_spatial_and_protein_context",
    "figureS1_core_composition_and_portability",
    "figureS2_external_and_ffpe_sensitivity",
    "figureS3_signature_transparency_and_random_benchmark",
    "figureS4_rna_atac_robustness",
    "figureS5_crc_atlas_source_audit",
    "figureS6_perturbation_boundaries",
    "figureS7_virtual_knockout_robustness",
    "figureS8_spatial_and_protein_assayability"
  ),
  expected_panels = c(6, 6, 6, 5, 6, 7, 4, 4, 4, 4, 4, 4, 3, 4)
)

formats <- c("svg", "pdf", "tiff", "png")
format_rows <- do.call(rbind, lapply(seq_len(nrow(expected)), function(i) {
  paths <- file.path(FIG_DIR, paste0(expected$stem[i], ".", formats))
  info <- file.info(paths)
  data.frame(
    figure = expected$figure[i], stem = expected$stem[i], format = formats,
    exists = file.exists(paths), nonempty = file.exists(paths) & info$size > 1000,
    bytes = ifelse(file.exists(paths), info$size, NA_real_), path = paths
  )
}))

svg_dimensions_mm <- function(stem) {
  path <- file.path(FIG_DIR, paste0(stem, ".svg"))
  header <- paste(readLines(path, n = 2, warn = FALSE), collapse = " ")
  match <- regexec(
    "width='([0-9.]+)pt' height='([0-9.]+)pt'", header,
    perl = TRUE
  )
  values <- regmatches(header, match)[[1]]
  if (length(values) != 3) stop("Could not read SVG dimensions: ", path)
  c(width = as.numeric(values[2]), height = as.numeric(values[3])) * 25.4 / 72
}

# Replace any inherited export manifest with one describing the current
# Communications Biology package. This prevents stale paths or canvas sizes
# from surviving after targeted figure revisions.
manifest_rows <- do.call(rbind, lapply(seq_len(nrow(expected)), function(i) {
  dimensions <- svg_dimensions_mm(expected$stem[i])
  paths <- file.path(FIG_DIR, paste0(expected$stem[i], ".", formats))
  data.frame(
    figure = expected$stem[i],
    format = toupper(formats),
    file = file.path(
      "figures", "communications_biology_v1.2",
      paste0(expected$stem[i], ".", formats)
    ),
    width_mm = round(dimensions[["width"]], 3),
    height_mm = round(dimensions[["height"]], 3),
    resolution_dpi = c(NA, NA, 600, 300),
    file_size_bytes = file.info(paths)$size,
    sha256 = vapply(
      paths,
      function(path) strsplit(
        system2("sha256sum", path, stdout = TRUE), "[[:space:]]+"
      )[[1]][1],
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}))
write.table(
  manifest_rows, file.path(SOURCE_DIR, "figure_export_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

read_svg <- function(stem) {
  path <- file.path(FIG_DIR, paste0(stem, ".svg"))
  if (!file.exists(path)) return("")
  paste(readLines(path, warn = FALSE), collapse = "\n")
}
svg_text <- setNames(vapply(expected$stem, read_svg, character(1)), expected$figure)

resource_landscape <- read.delim(
  file.path(SOURCE_DIR, "figure1a_resource_landscape.tsv"),
  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)

panel_rows <- do.call(rbind, lapply(seq_len(nrow(expected)), function(i) {
  text <- svg_text[[expected$figure[i]]]
  tags <- vapply(letters[seq_len(expected$expected_panels[i])], function(tag) {
    length(gregexpr(paste0(">", tag, "<"), text, fixed = TRUE)[[1]][
      gregexpr(paste0(">", tag, "<"), text, fixed = TRUE)[[1]] > 0
    ])
  }, integer(1))
  data.frame(
    figure = expected$figure[i], expected_panels = expected$expected_panels[i],
    detected_panel_tags = sum(tags), each_expected_tag_once = all(tags == 1)
  )
}))

contains_all <- function(text, terms) all(vapply(terms, grepl, logical(1), x = text, fixed = TRUE))
contains_none <- function(text, terms) !any(vapply(terms, grepl, logical(1), x = text, fixed = TRUE))

text_checks <- data.frame(
  check = c(
    "figure1_complete_resource_overview",
    "figure1_six_resource_icons_defined",
    "figure1_complete_selection_counts",
    "figure1_threshold_language",
    "figure1_no_internal_audit_copy",
    "figure1_typography_spacing",
    "figure1_no_legacy_abstract_axis",
    "figure2_numeric_precision_matches_results",
    "figure4_readable_study_labels",
    "supplementary4_no_internal_variable_names",
    "supplementary6_no_internal_variable_names",
    "supplementary6_standardised_gene_case",
    "figure5_frozen_signature_wording",
    "all_panel_tags_present_once"
  ),
  pass = c(
    contains_all(svg_text[["figure1"]], c(
      "Public resources (n = 22)", "sc/snRNA-seq", "Bulk and FFPE RNA",
      "Matched RNA–ATAC", "Perturbation RNA", "Spatial RNA", "Proteomics",
      "60 specimens · 27 donors", "44 specimens · 24 donors",
      "5 cohorts · 203 samples · 161 patient clusters",
      "51 FFPE pairs (102 specimens)", "40 samples · 12 patients",
      "6 datasets · 8 contrasts", "CRC Atlas", "6 tissue sections",
      "4 datasets"
    )),
    nrow(resource_landscape) == 6 &&
      length(unique(resource_landscape$icon_key)) == 6 &&
      all(resource_landscape$icon_key %in% c(
        "single_cell", "bulk_ffpe", "multiome", "perturbation", "spatial",
        "proteomics"
      )),
    contains_all(svg_text[["figure1"]], c(
      "33,698", "6,127", "2,496", "1,504", "287-gene core",
      "62 portable genes", "12-gene signature"
    )),
    contains_all(svg_text[["figure1"]], c(
      "Same direction in ≥90%", "BH FDR ≤0.05",
      "Cross-platform portability", "Frozen"
    )),
    contains_none(svg_text[["figure1"]], c(
      "DISCOVERY DATA ONLY", "GENE SELECTION PERMITTED",
      "NO GENE SELECTION OR RETUNING AFTER FREEZE", "PORTABILITY AUDIT",
      "labels hidden", "criteria in b", "fidelity in c",
      "FROZEN BEFORE ALL TESTS", "no preset gene count"
    )),
    contains_all(svg_text[["figure1"]], c(
      "Public resources (n = 22)", "Analysis design", "Discovery (1 resource)",
      "Held-out test (1 partition)"
    )) && contains_none(svg_text[["figure1"]], c(
      "PUBLIC RESOURCES  ·", "textLength="
    )),
    contains_none(svg_text[["figure1"]], c(
      "Genes retained (log scale)", "Expressed / non-technical",
      "Directionally stable", "Bootstrap CI excludes zero"
    )),
    contains_all(svg_text[["figure2"]], c(
      "1.830", "1.614", "1.000", "0.997", "0.993", "0.694"
    )),
    contains_all(svg_text[["figure4"]], c(
      "Chen 2021", "Pelka 2021", "MUI Innsbruck", "Source study"
    )) && contains_none(svg_text[["figure4"]], c("Chen_2021", "Pelka_2021")),
    contains_none(svg_text[["figureS4"]], c(
      "atac_wnt_stem", "atac_wnt_minus_housekeeping",
      "atac_wnt_tcf_ascl2_axis"
    )) && contains_all(svg_text[["figureS4"]], c(
      "WNT/stem accessibility", "WNT/TCF/ASCL2 − housekeeping"
    )),
    contains_none(svg_text[["figureS6"]], c(
      "route_up", "route_down", "wnt_stem", "proliferation_control",
      "epithelial_control"
    )) && contains_all(svg_text[["figureS6"]], c(
      "Signature up arm", "Epithelial identity control"
    )),
    contains_all(svg_text[["figureS6"]], c(
      "APC restoration + KRAS", "ASCL2 knockout", "Conditional WNT silencing"
    )) && contains_none(svg_text[["figureS6"]], c("Apc restoration", "Ascl2 KO")),
    contains_all(svg_text[["figure5"]], c(
      "FROZEN 12-GENE SIGNATURE", "Up arm · 6 genes", "Down arm · 6 genes"
    )),
    all(panel_rows$each_expected_tag_once)
  )
)

tsv_paths <- list.files(SOURCE_DIR, pattern = "tsv$", full.names = TRUE)
tsv_rows <- do.call(rbind, lapply(tsv_paths, function(path) {
  fields <- count.fields(path, sep = "\t", quote = "", blank.lines.skip = FALSE)
  data.frame(
    file = basename(path), lines = length(fields),
    expected_fields = if (length(fields)) fields[1] else NA_integer_,
    min_fields = if (length(fields)) min(fields) else NA_integer_,
    max_fields = if (length(fields)) max(fields) else NA_integer_,
    rectangular = length(fields) > 0 && all(fields == fields[1])
  )
}))

summary_checks <- data.frame(
  check = c(
    "14_canonical_figures", "all_four_formats_present", "all_files_nonempty",
    "canonical_manifest_current", "panel_tag_contracts", "text_contracts",
    "source_tsvs_rectangular"
  ),
  pass = c(
    nrow(expected) == 14,
    all(format_rows$exists),
    all(format_rows$nonempty),
    nrow(manifest_rows) == nrow(expected) * length(formats) &&
      all(grepl("figures/communications_biology_v1.2/", manifest_rows$file,
                fixed = TRUE)) &&
      all(manifest_rows$width_mm == 170) &&
      all(manifest_rows$height_mm <= 225),
    all(panel_rows$each_expected_tag_once),
    all(text_checks$pass),
    all(tsv_rows$rectangular)
  )
)

write.table(format_rows, file.path(SOURCE_DIR, "figure_package_format_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(panel_rows, file.path(SOURCE_DIR, "figure_package_panel_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(text_checks, file.path(SOURCE_DIR, "figure_package_text_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(tsv_rows, file.path(SOURCE_DIR, "figure_package_source_data_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(summary_checks, file.path(SOURCE_DIR, "figure_package_audit_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

if (!all(summary_checks$pass)) {
  stop("Figure package audit failed: ",
       paste(summary_checks$check[!summary_checks$pass], collapse = ", "))
}
message("Figure package audit PASS: 14 figures, four formats, text and source-data contracts")
