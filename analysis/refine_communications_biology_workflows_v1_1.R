#!/usr/bin/env Rscript

# Communications Biology workflow refinement.
#
# The analytical panels and numerical results are inherited unchanged from the
# locked v2.8 figure package. This script replaces only Fig. 1a and Fig. 3a,
# then exports a complete journal-specific figure directory from R.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(grid)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  file.path(getwd(), "analysis", "refine_communications_biology_workflows_v1_1.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
LEGACY_DIR <- file.path(ROOT, "figures", "jtm_submission_v2.8")
OUT_DIR <- file.path(ROOT, "figures", "communications_biology_v1.1")
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("dplyr", "ggplot2", "patchwork", "scales", "svglite", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing locked-environment R packages: ", paste(missing_packages, collapse = ", "))
}

# Load the verified plots into an isolated environment. The v2.8 renderer keeps
# its unchanged v2.6 panels in a nested environment, so the required objects are
# collected explicitly after all locked numerical assertions have run.
env_v28 <- new.env(parent = globalenv())
sys.source(
  file.path(ROOT, "analysis", "plot_jtm_submission_figures_v2_8.R"),
  envir = env_v28
)
env_v26 <- env_v28$base
required_v28 <- c("COL", "p1b", "p1c", "p1e")
required_v26 <- c("p1d", "p1f", "p3b", "p3c", "p3d", "p3e", "p3f")
if (!all(vapply(required_v28, exists, logical(1), envir = env_v28,
                inherits = FALSE)) ||
    !all(vapply(required_v26, exists, logical(1), envir = env_v26,
                inherits = FALSE))) {
  stop("The verified legacy renderer did not expose all required panels")
}
base <- list(
  COL = env_v28$COL,
  p1b = env_v28$p1b,
  p1c = env_v28$p1c,
  p1d = env_v26$p1d,
  p1e = env_v28$p1e,
  p1f = env_v26$p1f,
  p3b = env_v26$p3b,
  p3c = env_v26$p3c,
  p3d = env_v26$p3d,
  p3e = env_v26$p3e,
  p3f = env_v26$p3f
)

# Start with an exact copy of the verified figure package, then overwrite only
# the two figures whose schematic panels change.
legacy_files <- list.files(LEGACY_DIR, recursive = TRUE, full.names = TRUE,
                           include.dirs = FALSE)
for (source in legacy_files) {
  relative <- substring(source, nchar(LEGACY_DIR) + 2L)
  destination <- file.path(OUT_DIR, relative)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, destination, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Could not copy legacy figure asset: ", relative)
  }
}

JOURNAL_FONT <- "Arial"
COL <- base$COL
clean_panel <- function(plot) {
  plot + labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL) +
    theme(
      text = element_text(family = JOURNAL_FONT, colour = COL[["ink"]]),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.tag = element_text(size = 8.2, face = "bold", family = JOURNAL_FONT)
    )
}
tag_theme <- theme(
  plot.tag.position = c(0, 1),
  plot.tag = element_text(size = 8.2, face = "bold", family = JOURNAL_FONT,
                          colour = COL[["ink"]])
)

card_grob <- function(fill, border, radius_mm = 1.8, lwd = 0.9) {
  roundrectGrob(
    r = unit(radius_mm, "mm"),
    gp = gpar(fill = fill, col = border, lwd = lwd, linejoin = "round")
  )
}

icon_grob <- function(icon, accent) {
  ink <- unname(COL[["ink"]])
  neutral <- unname(COL[["neutral"]])
  pale <- alpha(accent, 0.13)
  line_gp <- gpar(col = ink, lwd = 0.9, lineend = "round", linejoin = "round")
  accent_gp <- gpar(col = accent, fill = accent, lwd = 0.9)
  open_gp <- gpar(col = accent, fill = "white", lwd = 1.0)

  switch(
    icon,
    donors = grobTree(
      circleGrob(x = c(0.22, 0.50, 0.78), y = 0.69, r = 0.075, gp = accent_gp),
      segmentsGrob(x0 = c(0.22, 0.50, 0.78), x1 = c(0.22, 0.50, 0.78),
                   y0 = 0.57, y1 = 0.32, gp = line_gp),
      segmentsGrob(x0 = c(0.12, 0.40, 0.68), x1 = c(0.32, 0.60, 0.88),
                   y0 = 0.46, y1 = 0.46, gp = line_gp),
      segmentsGrob(x0 = 0.15, x1 = 0.85, y0 = 0.20, y1 = 0.20,
                   gp = gpar(col = neutral, lwd = 0.7))
    ),
    lock = grobTree(
      roundrectGrob(x = 0.34, y = 0.50, width = 0.48, height = 0.62,
                    r = unit(1.0, "mm"),
                    gp = gpar(col = ink, fill = "white", lwd = 0.9)),
      segmentsGrob(x0 = 0.18, x1 = 0.46, y0 = c(0.66, 0.50, 0.34),
                   y1 = c(0.66, 0.50, 0.34),
                   gp = gpar(col = accent, lwd = 1.2, lineend = "round")),
      circleGrob(x = 0.74, y = 0.62, r = 0.14, gp = open_gp),
      rectGrob(x = 0.74, y = 0.40, width = 0.32, height = 0.34,
               gp = gpar(col = accent, fill = pale, lwd = 1.0)),
      circleGrob(x = 0.74, y = 0.44, r = 0.03, gp = accent_gp),
      segmentsGrob(x0 = 0.74, x1 = 0.74, y0 = 0.41, y1 = 0.31,
                   gp = gpar(col = accent, lwd = 1.0))
    ),
    replication = grobTree(
      rectGrob(x = 0.34, y = 0.55, width = 0.40, height = 0.52,
               gp = gpar(col = neutral, fill = "white", lwd = 0.8)),
      rectGrob(x = 0.46, y = 0.45, width = 0.40, height = 0.52,
               gp = gpar(col = ink, fill = "white", lwd = 0.9)),
      segmentsGrob(x0 = 0.33, x1 = 0.56, y0 = c(0.56, 0.46, 0.36),
                   y1 = c(0.56, 0.46, 0.36),
                   gp = gpar(col = accent, lwd = 1.0, lineend = "round")),
      circleGrob(x = 0.78, y = 0.39, r = 0.16, gp = accent_gp),
      polylineGrob(x = c(0.69, 0.76, 0.88), y = c(0.39, 0.31, 0.49),
                   gp = gpar(col = "white", lwd = 1.6, lineend = "round",
                             linejoin = "round"))
    ),
    network = grobTree(
      segmentsGrob(x0 = c(0.50, 0.50, 0.50, 0.50),
                   y0 = c(0.50, 0.50, 0.50, 0.50),
                   x1 = c(0.20, 0.34, 0.68, 0.82),
                   y1 = c(0.30, 0.76, 0.76, 0.30),
                   gp = gpar(col = neutral, lwd = 0.85)),
      circleGrob(x = 0.50, y = 0.50, r = 0.12, gp = accent_gp),
      circleGrob(x = c(0.20, 0.34, 0.68, 0.82),
                 y = c(0.30, 0.76, 0.76, 0.30), r = 0.07, gp = open_gp)
    ),
    tissue = grobTree(
      rectGrob(x = 0.50, y = 0.50, width = 0.72, height = 0.58,
               gp = gpar(col = ink, fill = "white", lwd = 0.9)),
      segmentsGrob(x0 = c(0.38, 0.62), x1 = c(0.38, 0.62),
                   y0 = 0.22, y1 = 0.78,
                   gp = gpar(col = neutral, lwd = 0.65)),
      segmentsGrob(x0 = 0.14, x1 = 0.86, y0 = c(0.41, 0.59), y1 = c(0.41, 0.59),
                   gp = gpar(col = neutral, lwd = 0.65)),
      circleGrob(x = c(0.26, 0.50, 0.74), y = c(0.68, 0.50, 0.32),
                 r = 0.065, gp = accent_gp)
    ),
    rna = grobTree(
      rectGrob(x = 0.50, y = 0.50, width = 0.66, height = 0.58,
               gp = gpar(col = ink, fill = "white", lwd = 0.9)),
      segmentsGrob(x0 = 0.28, x1 = 0.72, y0 = c(0.66, 0.52, 0.38),
                   y1 = c(0.66, 0.52, 0.38),
                   gp = gpar(col = accent, lwd = 1.2, lineend = "round")),
      circleGrob(x = 0.76, y = 0.23, r = 0.11, gp = accent_gp)
    ),
    atac = grobTree(
      polylineGrob(x = c(0.25, 0.38, 0.62, 0.75),
                   y = c(0.73, 0.35, 0.65, 0.27),
                   gp = gpar(col = accent, lwd = 1.4, lineend = "round")),
      polylineGrob(x = c(0.25, 0.38, 0.62, 0.75),
                   y = c(0.27, 0.65, 0.35, 0.73),
                   gp = gpar(col = ink, lwd = 1.0, lineend = "round")),
      segmentsGrob(x0 = c(0.31, 0.44, 0.56, 0.69),
                   x1 = c(0.31, 0.44, 0.56, 0.69),
                   y0 = c(0.56, 0.47, 0.47, 0.56),
                   y1 = c(0.44, 0.53, 0.53, 0.44),
                   gp = gpar(col = neutral, lwd = 0.7))
    ),
    association = grobTree(
      segmentsGrob(x0 = 0.20, x1 = 0.78, y0 = 0.24, y1 = 0.24, gp = line_gp),
      segmentsGrob(x0 = 0.20, x1 = 0.20, y0 = 0.24, y1 = 0.78, gp = line_gp),
      segmentsGrob(x0 = 0.26, x1 = 0.74, y0 = 0.31, y1 = 0.70,
                   gp = gpar(col = neutral, lwd = 1.0)),
      circleGrob(x = c(0.30, 0.40, 0.52, 0.63, 0.72),
                 y = c(0.35, 0.43, 0.49, 0.63, 0.66), r = 0.045,
                 gp = accent_gp)
    )
  )
}

add_card <- function(plot, xmin, xmax, ymin, ymax, fill, border, icon,
                     title, primary, secondary, badge = NULL,
                     title_colour = border, icon_fraction = 0.31) {
  xmid <- (xmin + xmax) / 2
  width <- xmax - xmin
  plot <- plot +
    annotation_custom(card_grob(fill, border), xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax) +
    annotation_custom(
      card_grob("white", alpha(border, 0.36), radius_mm = 1.2, lwd = 0.7),
      xmin = xmid - width * 0.18, xmax = xmid + width * 0.18,
      ymin = ymax - (ymax - ymin) * 0.46,
      ymax = ymax - (ymax - ymin) * 0.12
    ) +
    annotation_custom(
      icon_grob(icon, border),
      xmin = xmid - width * icon_fraction / 2,
      xmax = xmid + width * icon_fraction / 2,
      ymin = ymax - (ymax - ymin) * 0.42,
      ymax = ymax - (ymax - ymin) * 0.15
    ) +
    annotate("segment", x = xmin + width * 0.13, xend = xmax - width * 0.13,
             y = ymax - 0.07, yend = ymax - 0.07, colour = border,
             linewidth = 1.05, lineend = "round") +
    annotate("text", x = xmid, y = ymin + (ymax - ymin) * 0.45,
             label = title, family = JOURNAL_FONT, fontface = "bold",
             size = 1.65, colour = title_colour) +
    annotate("text", x = xmid, y = ymin + (ymax - ymin) * 0.29,
             label = primary, family = JOURNAL_FONT, fontface = "bold",
             size = 1.47, colour = COL[["ink"]], lineheight = 0.94) +
    annotate("text", x = xmid, y = ymin + (ymax - ymin) * 0.14,
             label = secondary, family = JOURNAL_FONT, size = 1.30,
             colour = COL[["neutral"]], lineheight = 0.92)
  if (!is.null(badge)) {
    plot <- plot +
      annotate("label", x = xmin + width * 0.13, y = ymax - 0.17,
               label = badge, family = JOURNAL_FONT, fontface = "bold",
               size = 1.12, fill = COL[["ink"]], colour = "white",
               linewidth = 0, label.padding = unit(0.58, "mm"),
               label.r = unit(1.8, "mm"))
  }
  plot
}

# ---------------------------------------------------------------------------
# Figure 1a: one complete workflow with a branching evidence architecture
# ---------------------------------------------------------------------------

p1a <- ggplot() +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 3.25), clip = "off") +
  theme_void(base_family = JOURNAL_FONT) +
  theme(plot.margin = margin(1.8, 2.2, 1.0, 2.2, "mm"))

p1a <- add_card(
  p1a, 0.10, 2.35, 0.48, 2.93, "#FCF6F3", COL[["route"]], "donors",
  "DEFINE THE STATE", "27 discovery donors", "donor-aware FDR\n+ bootstrap stability"
)
p1a <- add_card(
  p1a, 2.75, 5.00, 0.48, 2.93, "#F5F8FA", COL[["wnt"]], "lock",
  "REDUCE AND FREEZE", "287 → 62 → 12 genes", "platform gate + grouped fit\nvalidation outcomes unseen"
)

p1a <- p1a +
  annotate("segment", x = 2.37, xend = 2.70, y = 1.70, yend = 1.70,
           colour = COL[["neutral"]], linewidth = 0.52,
           arrow = arrow(type = "closed", angle = 24, length = unit(1.35, "mm"))) +
  annotate("segment", x = 5.02, xend = 5.38, y = 1.70, yend = 1.70,
           colour = COL[["neutral"]], linewidth = 0.52,
           arrow = arrow(type = "closed", angle = 24, length = unit(1.35, "mm"))) +
  annotation_custom(card_grob("#FAFBFB", COL[["neutral_light"]], lwd = 0.85),
                    xmin = 5.45, xmax = 11.90, ymin = 0.28, ymax = 3.07) +
  annotate("segment", x = 5.75, xend = 11.60, y = 2.98, yend = 2.98,
           colour = COL[["ink"]], linewidth = 0.92, lineend = "round") +
  annotate("text", x = 8.68, y = 2.75, label = "TEST THE FROZEN SIGNATURE",
           family = JOURNAL_FONT, fontface = "bold", size = 1.78,
           colour = COL[["ink"]])

p1a <- add_card(
  p1a, 5.70, 7.55, 0.66, 2.48, "#FCF9F1", COL[["adenoma"]], "replication",
  "REPLICATION", "held-out + 5 cohorts", "51 paired FFPE specimens",
  icon_fraction = 0.38
)
p1a <- add_card(
  p1a, 7.75, 9.60, 0.66, 2.48, "#F7F5FA", COL[["context"]], "network",
  "PATHWAY SUPPORT", "RNA–ATAC", "empirical + virtual\nperturbation",
  icon_fraction = 0.38
)
p1a <- add_card(
  p1a, 9.80, 11.65, 0.66, 2.48, "#F3F8F7", COL[["crc"]], "tissue",
  "TISSUE CONTEXT", "CRC Atlas + spatial", "public proteomic anchors",
  icon_fraction = 0.38
)
p1a <- p1a +
  annotate("text", x = 8.68, y = 0.43,
           label = "FROZEN OUTPUT  ·  12-gene research signature",
           family = JOURNAL_FONT, fontface = "bold", size = 1.30,
           colour = COL[["neutral"]])

# ---------------------------------------------------------------------------
# Figure 3a: paired multiomic design, not a second study workflow
# ---------------------------------------------------------------------------

p3a <- ggplot() +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 2.35), clip = "off") +
  theme_void(base_family = JOURNAL_FONT) +
  theme(plot.margin = margin(1.2, 2.0, 0.8, 2.0, "mm"))

p3a <- add_card(
  p3a, 0.20, 2.45, 0.34, 2.05, "#F7F8F9", COL[["neutral"]], "donors",
  "MATCHED SPECIMENS", "40 samples", "12 patients · normal/polyp",
  icon_fraction = 0.34
)
p3a <- add_card(
  p3a, 3.25, 6.15, 1.24, 2.18, "#FCF6F3", COL[["route"]], "rna",
  "EPITHELIAL snRNA-seq", "frozen 12-gene RNA score", "aggregated by patient and state",
  icon_fraction = 0.22
)
p3a <- add_card(
  p3a, 3.25, 6.15, 0.10, 1.04, "#F5F8FA", COL[["wnt"]], "atac",
  "MATCHED ATAC-seq", "promoter accessibility", "WNT and TCF/ASCL2 axes",
  icon_fraction = 0.22
)
p3a <- add_card(
  p3a, 8.05, 11.75, 0.48, 1.86, "#F7F5FA", COL[["context"]], "association",
  "PATIENT-AWARE ASSOCIATION", "paired RNA–ATAC evidence", "adjusted for tissue state,\nproliferation and depth",
  icon_fraction = 0.20
)

p3a <- p3a +
  annotate("segment", x = 2.48, xend = 2.88, y = 1.20, yend = 1.20,
           colour = COL[["neutral"]], linewidth = 0.50) +
  annotate("segment", x = 2.88, xend = 3.18, y = 1.20, yend = 1.71,
           colour = COL[["neutral"]], linewidth = 0.50,
           arrow = arrow(type = "closed", angle = 24, length = unit(1.25, "mm"))) +
  annotate("segment", x = 2.88, xend = 3.18, y = 1.20, yend = 0.57,
           colour = COL[["neutral"]], linewidth = 0.50,
           arrow = arrow(type = "closed", angle = 24, length = unit(1.25, "mm"))) +
  annotate("segment", x = 6.18, xend = 7.18, y = 1.71, yend = 1.18,
           colour = COL[["neutral"]], linewidth = 0.50) +
  annotate("segment", x = 6.18, xend = 7.18, y = 0.57, yend = 1.18,
           colour = COL[["neutral"]], linewidth = 0.50) +
  annotate("segment", x = 7.18, xend = 7.98, y = 1.18, yend = 1.18,
           colour = COL[["neutral"]], linewidth = 0.50,
           arrow = arrow(type = "closed", angle = 24, length = unit(1.25, "mm")))

# Reassemble the two figures. Every quantitative panel remains unchanged.
fig1 <- clean_panel(p1a) /
  (clean_panel(base$p1b) | clean_panel(base$p1c)) /
  (clean_panel(base$p1d) | clean_panel(base$p1e) | clean_panel(base$p1f)) +
  plot_layout(heights = c(1.22, 1, 1), widths = c(1.05, 1, 0.82)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

fig3 <- clean_panel(p3a) /
  (clean_panel(base$p3b) | clean_panel(base$p3c)) /
  (clean_panel(base$p3d) | clean_panel(base$p3e) | clean_panel(base$p3f)) +
  plot_layout(heights = c(0.66, 1.05, 0.95), widths = c(1.02, 0.94, 0.94)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

export_figure <- function(plot, stem, width_mm, height_mm) {
  paths <- c(
    SVG = file.path(OUT_DIR, paste0(stem, ".svg")),
    PDF = file.path(OUT_DIR, paste0(stem, ".pdf")),
    TIFF = file.path(OUT_DIR, paste0(stem, ".tiff")),
    PNG = file.path(OUT_DIR, paste0(stem, ".png"))
  )
  ggsave(paths[["SVG"]], plot = plot, device = svglite::svglite,
         width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(paths[["PDF"]], plot = plot, device = grDevices::cairo_pdf,
         width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(paths[["TIFF"]], plot = plot, device = ragg::agg_tiff,
         width = width_mm, height = height_mm, units = "mm", dpi = 600,
         compression = "lzw", background = "white")
  ggsave(paths[["PNG"]], plot = plot, device = ragg::agg_png,
         width = width_mm, height = height_mm, units = "mm", dpi = 300,
         background = "white")
  data.frame(
    figure = stem, format = names(paths), file = unname(paths),
    width_mm = width_mm, height_mm = height_mm,
    resolution_dpi = c(NA, NA, 600, 300), stringsAsFactors = FALSE
  )
}

exports <- bind_rows(
  export_figure(fig1, "figure1_discovery_core_and_objective_reduction", 170, 192),
  export_figure(fig3, "figure3_rna_atac_regulatory_support", 170, 164)
)
exports$file_size_bytes <- file.info(exports$file)$size

workflow_source <- data.frame(
  phase = c("Define the state", "Reduce and freeze", "Replication",
            "Pathway support", "Tissue context", "Boundary"),
  detail = c(
    "27 discovery donors; donor-aware FDR and bootstrap stability; 287-gene core",
    "Label-blind platform gate and grouped reconstruction; 287 to 62 to 12 genes; validation outcomes unseen",
    "Held-out source data, five external cohorts and 51 paired FFPE specimens",
    "Matched RNA-ATAC plus empirical and virtual perturbation",
    "CRC Atlas, spatial transcriptomics and public proteomic anchors",
    "Frozen 12-gene research signature; clinical utility not tested"
  ),
  visual_role = c("discovery", "locked reduction", "independent replication",
                  "regulatory and perturbation support", "localisation and assay context",
                  "claim boundary")
)
write.table(workflow_source, file.path(SOURCE_DIR, "figure1a_workflow.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

multiomic_source <- data.frame(
  component = c("Matched specimens", "Epithelial snRNA-seq", "Matched ATAC-seq",
                "Patient-aware association"),
  detail = c("40 normal or polyp samples from 12 patients",
             "Frozen 12-gene RNA score aggregated by patient and state",
             "WNT and TCF/ASCL2 promoter-accessibility axes",
             "Tissue-state, proliferation and depth adjusted; patient-clustered inference")
)
write.table(multiomic_source, file.path(SOURCE_DIR, "figure3a_multiomic_design.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

svg_text <- vapply(
  exports$file[exports$format == "SVG"],
  function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
  character(1)
)
qa <- data.frame(
  check = c("four_formats_each", "all_files_nonempty", "editable_svg_text",
            "single_global_workflow", "figure3_is_paired_design", "no_legacy_100_gene_text"),
  pass = c(
    all(table(exports$figure) == 4L),
    all(exports$file_size_bytes > 0),
    all(grepl("<text", svg_text, fixed = TRUE)),
    grepl("TEST THE FROZEN SIGNATURE", svg_text[[1]], fixed = TRUE),
    grepl("MATCHED SPECIMENS", svg_text[[2]], fixed = TRUE),
    !any(grepl("100-gene|50 up + 50 down|10-gene", svg_text))
  )
)
write.table(qa, file.path(SOURCE_DIR, "workflow_figure_qa.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(qa$pass)) {
  stop("Workflow figure QA failed: ", paste(qa$check[!qa$pass], collapse = ", "))
}

session_lines <- capture.output(print(sessionInfo()))
session_lines <- gsub(ROOT, "<REPOSITORY_ROOT>", session_lines, fixed = TRUE)
writeLines(session_lines, file.path(SOURCE_DIR, "workflow_figure_sessionInfo.txt"),
           useBytes = TRUE)

message("Wrote Communications Biology workflow figures: ", OUT_DIR)
