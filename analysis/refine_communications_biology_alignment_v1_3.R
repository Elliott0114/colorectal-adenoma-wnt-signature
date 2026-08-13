#!/usr/bin/env Rscript

# Targeted final alignment pass for the Communications Biology figure package.
#
# Only panels with a confirmed alignment or schematic defect are redrawn.
# Numerical data, geoms, statistics and panel meanings are inherited unchanged
# from the audited v1.2 renderer or its frozen panel-level source data.

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
  file.path(getwd(), "analysis", "refine_communications_biology_alignment_v1_3.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT_DIR <- file.path(ROOT, "figures", "communications_biology_v1.2")
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
RESULT_DIR <- file.path(ROOT, "results", "objective_compact_panel_v2_7")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

# The audited renderer is standalone and exposes the verified Fig. 2, Fig. 4,
# Supplementary Fig. 4 and Supplementary Fig. 6 panel objects.
audit_env <- new.env(parent = globalenv())
sys.source(
  file.path(ROOT, "analysis", "refine_communications_biology_figure_audit_fixes_v1_2.R"),
  envir = audit_env
)

COL <- audit_env$COL
JOURNAL_FONT <- audit_env$JOURNAL_FONT
clean_panel <- audit_env$clean_panel
tag_theme <- audit_env$tag_theme
theme_journal <- audit_env$theme_journal

read_tsv <- function(path) {
  read.delim(
    path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
    quote = "", comment.char = "", fileEncoding = "UTF-8"
  )
}

read_source <- function(filename) read_tsv(file.path(SOURCE_DIR, filename))

no_guides <- function(plot) {
  plot + guides(
    colour = "none", fill = "none", shape = "none",
    linetype = "none", linewidth = "none", alpha = "none", size = "none"
  )
}

legend_band_theme <- theme(
  legend.position = "top",
  legend.direction = "horizontal",
  legend.box = "horizontal",
  legend.box.just = "center",
  legend.justification = "center",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.spacing.x = unit(1.2, "mm"),
  legend.spacing.y = unit(0, "mm"),
  legend.key.height = unit(2.2, "mm"),
  legend.key.width = unit(3.2, "mm")
)

export_figure <- function(plot, stem, width_mm, height_mm) {
  paths <- c(
    SVG = file.path(OUT_DIR, paste0(stem, ".svg")),
    PDF = file.path(OUT_DIR, paste0(stem, ".pdf")),
    TIFF = file.path(OUT_DIR, paste0(stem, ".tiff")),
    PNG = file.path(OUT_DIR, paste0(stem, ".png"))
  )
  ggsave(
    paths[["SVG"]], plot = plot,
    device = function(...) svglite::svglite(..., fix_text_size = FALSE),
    width = width_mm, height = height_mm, units = "mm", bg = "white"
  )
  ggsave(
    paths[["PDF"]], plot = plot, device = grDevices::cairo_pdf,
    width = width_mm, height = height_mm, units = "mm", bg = "white"
  )
  ggsave(
    paths[["TIFF"]], plot = plot, device = ragg::agg_tiff,
    width = width_mm, height = height_mm, units = "mm", dpi = 600,
    compression = "lzw", background = "white"
  )
  ggsave(
    paths[["PNG"]], plot = plot, device = ragg::agg_png,
    width = width_mm, height = height_mm, units = "mm", dpi = 300,
    background = "white"
  )
  data.frame(
    figure = stem, format = names(paths), file = unname(paths),
    width_mm = width_mm, height_mm = height_mm,
    resolution_dpi = c(NA, NA, 600, 300), stringsAsFactors = FALSE
  )
}

grid_2x2_with_guide <- function(p_a, p_b, p_c, p_d, widths = c(1, 1),
                                heights = c(1, 0.13, 1)) {
  design <- c(
    area(1, 1), area(1, 2),
    area(2, 1, 2, 2),
    area(3, 1), area(3, 2)
  )
  p_a + p_b + guide_area() + p_c + p_d +
    plot_layout(
      design = design, guides = "collect", heights = heights, widths = widths
    ) +
    plot_annotation(tag_levels = "a", theme = tag_theme) &
    legend_band_theme
}

# -----------------------------------------------------------------------------
# Figure 2: one shared arm legend between two strict three-column rows
# -----------------------------------------------------------------------------

fig2 <-
  no_guides(clean_panel(audit_env$p2a)) +
  no_guides(clean_panel(audit_env$p2b)) +
  no_guides(clean_panel(audit_env$p2c)) +
  guide_area() +
  no_guides(clean_panel(audit_env$p2d)) +
  clean_panel(audit_env$p2e) +
  no_guides(clean_panel(audit_env$p2f)) +
  plot_layout(
    design = c(
      area(1, 1), area(1, 2), area(1, 3),
      area(2, 1, 2, 3),
      area(3, 1), area(3, 2), area(3, 3)
    ),
    guides = "collect", heights = c(0.92, 0.13, 1.08),
    widths = c(1, 0.96, 0.91)
  ) +
  plot_annotation(tag_levels = "a", theme = tag_theme) &
  legend_band_theme

# Figure 4 already had reader-facing study labels; its source-study legend is
# moved into a dedicated band so panels d and e share the same plotting top.
fig4 <-
  no_guides(clean_panel(audit_env$p4a)) +
  no_guides(clean_panel(audit_env$p4b)) +
  no_guides(clean_panel(audit_env$p4c)) +
  guide_area() +
  clean_panel(audit_env$p4d) +
  no_guides(clean_panel(audit_env$p4e)) +
  plot_layout(
    design = c(
      area(1, 1, 1, 2),
      area(2, 1), area(2, 2),
      area(3, 1, 3, 2),
      area(4, 1), area(4, 2)
    ),
    guides = "collect", heights = c(0.82, 0.86, 0.13, 1.48),
    widths = c(1.02, 0.98)
  ) +
  plot_annotation(tag_levels = "a", theme = tag_theme) &
  legend_band_theme

# -----------------------------------------------------------------------------
# Supplementary Figs. 2, 4, 5 and 6: explicit two-column grids
# -----------------------------------------------------------------------------

loo_cohort <- read_tsv(file.path(
  RESULT_DIR, "validation_external_adjusted_leave_one_cohort_out.tsv"
)) %>%
  filter(panel_id == "objective_12") %>%
  mutate(excluded_cohort = factor(excluded_cohort, levels = rev(excluded_cohort)))
pS2a <- ggplot(loo_cohort, aes(adenoma_coef_sd, excluded_cohort)) +
  geom_vline(
    xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35
  ) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = excluded_cohort),
    linewidth = 0.70, colour = COL[["route"]]
  ) +
  geom_point(size = 1.9, colour = COL[["route"]]) +
  labs(x = "Adjusted effect after cohort omission", y = "Excluded cohort") +
  theme_journal(base_size = 6.0)

figS2 <- grid_2x2_with_guide(
  no_guides(clean_panel(pS2a)), clean_panel(audit_env$p2e),
  no_guides(clean_panel(audit_env$p2a)), no_guides(clean_panel(audit_env$p2d)),
  heights = c(1, 0.12, 1)
)

figS4 <- grid_2x2_with_guide(
  no_guides(clean_panel(audit_env$pS4a)),
  no_guides(clean_panel(audit_env$pS4b)),
  no_guides(clean_panel(audit_env$pS4c)),
  clean_panel(audit_env$pS4d),
  heights = c(1, 0.12, 1)
)

# Supplementary Fig. 5 now shares the same author-year labels as main Fig. 4.
figS5 <- grid_2x2_with_guide(
  no_guides(clean_panel(audit_env$p4b)),
  no_guides(clean_panel(audit_env$p4c)),
  clean_panel(audit_env$p4d),
  no_guides(clean_panel(audit_env$p4e)),
  widths = c(1.02, 0.98), heights = c(0.86, 0.13, 1.42)
)

figS6 <- grid_2x2_with_guide(
  clean_panel(audit_env$pS6a),
  no_guides(clean_panel(audit_env$pS6b)),
  no_guides(clean_panel(audit_env$pS6c)),
  clean_panel(audit_env$pS6d),
  heights = c(1, 0.16, 1)
)

# -----------------------------------------------------------------------------
# Supplementary Figure 7: four-stage workflow plus unchanged robustness panels
# -----------------------------------------------------------------------------

card_grob <- function(fill, border, radius_mm = 1.7, lwd = 0.85) {
  roundrectGrob(
    r = unit(radius_mm, "mm"),
    gp = gpar(fill = fill, col = border, lwd = lwd, linejoin = "round")
  )
}

workflow_icon <- function(icon, accent) {
  ink <- unname(COL[["ink"]])
  neutral <- unname(COL[["neutral"]])
  pale <- alpha(accent, 0.13)
  line_gp <- gpar(col = ink, lwd = 0.9, lineend = "round", linejoin = "round")
  accent_gp <- gpar(col = accent, fill = accent, lwd = 0.9)
  open_gp <- gpar(col = accent, fill = "white", lwd = 1.0)

  switch(
    icon,
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
    donors = grobTree(
      circleGrob(x = c(0.22, 0.50, 0.78), y = 0.69, r = 0.075,
                 gp = accent_gp),
      segmentsGrob(x0 = c(0.22, 0.50, 0.78), x1 = c(0.22, 0.50, 0.78),
                   y0 = 0.57, y1 = 0.32, gp = line_gp),
      segmentsGrob(x0 = c(0.12, 0.40, 0.68), x1 = c(0.32, 0.60, 0.88),
                   y0 = 0.46, y1 = 0.46, gp = line_gp),
      segmentsGrob(x0 = 0.15, x1 = 0.85, y0 = 0.20, y1 = 0.20,
                   gp = gpar(col = neutral, lwd = 0.7))
    ),
    network = grobTree(
      segmentsGrob(
        x0 = c(0.50, 0.50, 0.50, 0.50), y0 = c(0.50, 0.50, 0.50, 0.50),
        x1 = c(0.20, 0.34, 0.68, 0.82), y1 = c(0.30, 0.76, 0.76, 0.30),
        gp = gpar(col = neutral, lwd = 0.85)
      ),
      circleGrob(x = 0.50, y = 0.50, r = 0.12, gp = accent_gp),
      circleGrob(x = c(0.20, 0.34, 0.68, 0.82),
                 y = c(0.30, 0.76, 0.76, 0.30), r = 0.07, gp = open_gp),
      segmentsGrob(x0 = 0.34, x1 = 0.66, y0 = 0.66, y1 = 0.34,
                   gp = gpar(col = "white", lwd = 2.4)),
      segmentsGrob(x0 = 0.34, x1 = 0.66, y0 = 0.66, y1 = 0.34,
                   gp = gpar(col = accent, lwd = 1.0))
    ),
    validation = grobTree(
      rectGrob(x = 0.32, y = 0.52, width = 0.44, height = 0.56,
               gp = gpar(col = ink, fill = "white", lwd = 0.9)),
      segmentsGrob(x0 = 0.19, x1 = 0.43, y0 = c(0.65, 0.50, 0.35),
                   y1 = c(0.65, 0.50, 0.35),
                   gp = gpar(col = accent, lwd = 1.1, lineend = "round")),
      circleGrob(x = 0.72, y = 0.47, r = 0.18, gp = open_gp),
      polylineGrob(x = c(0.62, 0.69, 0.83), y = c(0.47, 0.38, 0.58),
                   gp = gpar(col = accent, lwd = 1.6, lineend = "round",
                             linejoin = "round"))
    )
  )
}

add_workflow_card <- function(plot, xmin, xmax, fill, border, icon, step,
                              title, primary, secondary) {
  ymin <- 0.18
  ymax <- 2.52
  xmid <- (xmin + xmax) / 2
  width <- xmax - xmin
  plot +
    annotation_custom(
      card_grob(fill, border), xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
    ) +
    annotation_custom(
      card_grob("white", alpha(border, 0.34), radius_mm = 1.1, lwd = 0.70),
      xmin = xmid - width * 0.18, xmax = xmid + width * 0.18,
      ymin = 1.62, ymax = 2.31
    ) +
    annotation_custom(
      workflow_icon(icon, border),
      xmin = xmid - width * 0.13, xmax = xmid + width * 0.13,
      ymin = 1.72, ymax = 2.21
    ) +
    annotate(
      "segment", x = xmin + width * 0.13, xend = xmax - width * 0.13,
      y = ymax - 0.07, yend = ymax - 0.07, colour = border,
      linewidth = 1.0, lineend = "round"
    ) +
    annotate(
      "point", x = xmin + width * 0.13, y = ymax - 0.18,
      shape = 21, size = 3.4, stroke = 0.65, fill = border, colour = "white"
    ) +
    annotate(
      "text", x = xmin + width * 0.13, y = ymax - 0.18,
      label = step, family = JOURNAL_FONT, fontface = "bold", size = 1.04,
      colour = "white"
    ) +
    annotate(
      "text", x = xmid, y = 1.34, label = title,
      family = JOURNAL_FONT, fontface = "bold", size = 1.70, colour = border
    ) +
    annotate(
      "text", x = xmid, y = 0.91, label = primary,
      family = JOURNAL_FONT, fontface = "bold", size = 1.47,
      colour = COL[["ink"]], lineheight = 0.94
    ) +
    annotate(
      "text", x = xmid, y = 0.45, label = secondary,
      family = JOURNAL_FONT, size = 1.27, colour = COL[["neutral"]],
      lineheight = 0.92
    )
}

workflow_s7 <- data.frame(
  step = paste0("0", 1:4),
  title = c("PRESPECIFY", "ASSEMBLE", "DELETE IN SILICO", "TEST AGAINST NULLS"),
  primary = c(
    "287-gene core\n12-gene signature",
    "13 held-out donors\n1,664 epithelial cells",
    "13 targets\ntwo initialisation seeds",
    "10,000 matched sets\nper endpoint"
  ),
  secondary = c(
    "targets and weights fixed", "128 cells per donor",
    "unsigned GenKI impact", "KL primary · EMD sensitivity"
  ),
  stringsAsFactors = FALSE
)

pS7a <- ggplot() +
  coord_cartesian(xlim = c(0, 10.4), ylim = c(0.02, 2.72), clip = "off") +
  theme_void(base_family = JOURNAL_FONT) +
  theme(plot.margin = margin(1.6, 2.2, 1.0, 2.2, "mm"))

card_bounds <- data.frame(
  xmin = c(0.10, 2.73, 5.36, 7.99),
  xmax = c(2.35, 4.98, 7.61, 10.24),
  fill = c("#FCF6F3", "#F5F8FA", "#F2F8F7", "#FBF7EF"),
  border = c(COL[["route"]], COL[["wnt"]], COL[["pass"]], COL[["adenoma"]]),
  icon = c("lock", "donors", "network", "validation"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(card_bounds))) {
  pS7a <- add_workflow_card(
    pS7a, card_bounds$xmin[i], card_bounds$xmax[i],
    card_bounds$fill[i], card_bounds$border[i], card_bounds$icon[i],
    workflow_s7$step[i], workflow_s7$title[i], workflow_s7$primary[i],
    workflow_s7$secondary[i]
  )
}
for (x in c(2.38, 5.01, 7.64)) {
  pS7a <- pS7a + annotate(
    "segment", x = x, xend = x + 0.31, y = 1.34, yend = 1.34,
    colour = COL[["neutral"]], linewidth = 0.52,
    arrow = arrow(type = "closed", angle = 24, length = unit(1.25, "mm"))
  )
}

robustness <- read_source("figureS8b_seed_metric_robustness.tsv") %>%
  mutate(measure = factor(measure, levels = c("Across seeds", "KL vs EMD")))
robustness_summary <- robustness %>%
  group_by(measure) %>%
  summarise(
    median_rho = median(rho), min_rho = min(rho),
    min_target = target[which.min(rho)], .groups = "drop"
  )
pS7b <- ggplot(robustness, aes(measure, rho, colour = measure)) +
  geom_hline(
    yintercept = c(0.5, 0.75), linewidth = 0.3,
    linetype = c("dotted", "22"), colour = COL[["neutral_light"]]
  ) +
  geom_boxplot(
    width = 0.38, outlier.shape = NA, linewidth = 0.42, fill = "white"
  ) +
  geom_point(
    size = 1.45, alpha = 0.82,
    position = position_jitter(width = 0.09, height = 0, seed = 20260810)
  ) +
  geom_point(
    data = robustness_summary, aes(y = median_rho), shape = 95, size = 6.2
  ) +
  geom_text(
    data = filter(robustness_summary, min_rho < 0.75),
    aes(y = min_rho, label = paste0(min_target, " ", sprintf("%.2f", min_rho))),
    nudge_y = -0.035, family = JOURNAL_FONT, size = 1.6,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Across seeds" = COL[["route"]], "KL vs EMD" = "#3C8D88"),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(0.43, 1.01), breaks = c(0.5, 0.75, 1.0)) +
  labs(x = NULL, y = "Spearman ρ") +
  theme_journal(base_size = 5.8) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

internal_panel <- read_source("figureS8c_internal_panel_boundary.tsv") %>%
  mutate(target = factor(target, levels = target[order(matched_z)]))
pS7c <- ggplot(internal_panel, aes(matched_z, target, colour = arm_label)) +
  geom_vline(
    xintercept = 0, linetype = "22", linewidth = 0.35,
    colour = COL[["neutral"]]
  ) +
  geom_segment(aes(x = 0, xend = matched_z, yend = target), linewidth = 0.58) +
  geom_point(aes(fill = supported), shape = 21, size = 1.9, stroke = 0.60) +
  geom_text(
    aes(x = text_x, label = q_label, hjust = text_hjust),
    size = 1.25, family = JOURNAL_FONT, colour = COL[["ink"]]
  ) +
  scale_colour_manual(
    values = c("Up arm" = COL[["route"]], "Down arm" = COL[["wnt"]]),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(`TRUE` = COL[["route"]], `FALSE` = "white"), guide = "none"
  ) +
  coord_cartesian(
    xlim = c(min(internal_panel$matched_z) - 0.6,
             max(internal_panel$matched_z) + 1.15), clip = "off"
  ) +
  labs(x = "Panel member → remaining signature (matched-null z)", y = NULL) +
  theme_journal(base_size = 5.55) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 4.0))

figS7 <-
  clean_panel(pS7a) +
  guide_area() +
  no_guides(clean_panel(pS7b)) +
  clean_panel(pS7c) +
  plot_layout(
    design = c(
      area(1, 1, 1, 2), area(2, 1, 2, 2),
      area(3, 1), area(3, 2)
    ),
    guides = "collect", heights = c(0.68, 0.13, 1.10),
    widths = c(0.72, 1.28)
  ) +
  plot_annotation(tag_levels = "a", theme = tag_theme) &
  legend_band_theme

# -----------------------------------------------------------------------------
# Export and QA
# -----------------------------------------------------------------------------

exports <- bind_rows(
  export_figure(fig2, "figure2_independent_replication_and_ffpe", 170, 148),
  export_figure(fig4, "figure4_crc_atlas_cross_sectional_recurrence", 170, 202),
  export_figure(figS2, "figureS2_external_and_ffpe_sensitivity", 170, 142),
  export_figure(figS4, "figureS4_rna_atac_robustness", 170, 142),
  export_figure(figS5, "figureS5_crc_atlas_source_audit", 170, 150),
  export_figure(figS6, "figureS6_perturbation_boundaries", 170, 148),
  export_figure(figS7, "figureS7_virtual_knockout_robustness", 170, 118)
)

workflow_s7_source <- workflow_s7
workflow_s7_source$primary <- gsub("\n", "; ", workflow_s7_source$primary, fixed = TRUE)
write.table(
  workflow_s7_source,
  file.path(SOURCE_DIR, "figureS7a_virtual_deletion_workflow.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)
exports$file_size_bytes <- file.info(exports$file)$size
write.table(
  exports, file.path(SOURCE_DIR, "figure_alignment_export_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

read_svg <- function(stem) {
  paste(readLines(file.path(OUT_DIR, paste0(stem, ".svg")), warn = FALSE),
        collapse = "\n")
}
expected_panels <- c(
  figure2_independent_replication_and_ffpe = 6,
  figure4_crc_atlas_cross_sectional_recurrence = 5,
  figureS2_external_and_ffpe_sensitivity = 4,
  figureS4_rna_atac_robustness = 4,
  figureS5_crc_atlas_source_audit = 4,
  figureS6_perturbation_boundaries = 4,
  figureS7_virtual_knockout_robustness = 3
)
svg_text <- setNames(vapply(names(expected_panels), read_svg, character(1)),
                     names(expected_panels))
panel_tags_ok <- vapply(names(expected_panels), function(stem) {
  tags <- letters[seq_len(expected_panels[[stem]])]
  counts <- vapply(tags, function(tag) {
    matches <- gregexpr(paste0(">", tag, "<"), svg_text[[stem]], fixed = TRUE)[[1]]
    sum(matches > 0)
  }, integer(1))
  all(counts == 1L)
}, logical(1))

qa <- data.frame(
  check = c(
    "7_revised_figures_four_formats", "all_exports_nonempty",
    "all_panel_tags_once", "supplementary7_uses_four_stage_workflow",
    "supplementary7_has_no_legacy_audit_sentence",
    "supplementary5_has_reader_facing_study_labels",
    "svg_text_remains_editable"
  ),
  pass = c(
    nrow(exports) == 7L * 4L && all(table(exports$figure) == 4L),
    all(exports$file_size_bytes > 1000), all(panel_tags_ok),
    all(vapply(
      c("PRESPECIFY", "ASSEMBLE", "DELETE IN SILICO", "TEST AGAINST NULLS"),
      grepl, logical(1), x = svg_text[["figureS7_virtual_knockout_robustness"]],
      fixed = TRUE
    )),
    !grepl(
      "Frozen-gene validation only: no replacement, reweighting or target expansion",
      svg_text[["figureS7_virtual_knockout_robustness"]], fixed = TRUE
    ),
    all(vapply(
      c("Chen 2021", "Pelka 2021", "MUI Innsbruck"), grepl, logical(1),
      x = svg_text[["figureS5_crc_atlas_source_audit"]], fixed = TRUE
    )) && !grepl(
      "Chen_2021", svg_text[["figureS5_crc_atlas_source_audit"]], fixed = TRUE
    ),
    all(grepl("<text", svg_text, fixed = TRUE))
  ),
  stringsAsFactors = FALSE
)
write.table(
  qa, file.path(SOURCE_DIR, "figure_alignment_qa.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)
if (!all(qa$pass)) {
  stop("Figure alignment QA failed: ", paste(qa$check[!qa$pass], collapse = ", "))
}

writeLines(
  capture.output(sessionInfo()),
  file.path(SOURCE_DIR, "figure_alignment_sessionInfo.txt"), useBytes = TRUE
)
message("Targeted figure alignment package written to: ", OUT_DIR)
