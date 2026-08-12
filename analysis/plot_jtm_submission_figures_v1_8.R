#!/usr/bin/env Rscript

# Biological-narrative refinement for the six-main-figure Journal of
# Translational Medicine submission. Figure 1a is the only study workflow;
# Figure 2 is fully quantitative and Figure 5a is a biological mechanism map.
# All illustrations are original R vectors. This script only reads existing
# result tables and does not refit models or alter locked estimates.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  file.path(getwd(), "analysis", "plot_jtm_submission_figures_v1_8.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT_DIR <- Sys.getenv(
  "JTM_FIGURE_DIR",
  unset = file.path(ROOT, "figures", "jtm_submission_v1.8")
)
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "dplyr", "ggplot2", "ggrepel", "patchwork", "scales", "tidyr",
  "svglite", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing locked-environment R packages: ", paste(missing_packages, collapse = ", "))
}

read_tsv <- function(path) {
  read.delim(
    file.path(ROOT, path), check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

write_tsv <- function(frame, filename) {
  write.table(
    frame, file.path(SOURCE_DIR, filename), sep = "\t", row.names = FALSE,
    quote = FALSE, na = ""
  )
}

JTM_FONT <- "Arial"
COL <- c(
  ink = "#25292D",
  route = "#C85E3D",
  route_light = "#E9B6A5",
  wnt = "#356F9D",
  wnt_light = "#B8D0E1",
  neutral = "#727B84",
  neutral_light = "#D8DDE1",
  neutral_pale = "#F3F5F6",
  adenoma = "#D9A640",
  crc = "#3F9A8B",
  uncertain = "#9A8877",
  context = "#7C6AA3"
)

theme_jtm <- function(base_size = 7.1) {
  theme_classic(base_size = base_size, base_family = JTM_FONT) +
    theme(
      text = element_text(family = JTM_FONT, colour = COL[["ink"]]),
      axis.line = element_line(linewidth = 0.32, colour = COL[["ink"]]),
      axis.ticks = element_line(linewidth = 0.28, colour = COL[["ink"]]),
      axis.ticks.length = grid::unit(0.9, "mm"),
      axis.title = element_text(size = base_size, colour = COL[["ink"]], margin = margin(t = 1.2, r = 1.2)),
      axis.text = element_text(size = base_size - 0.55, colour = COL[["ink"]]),
      panel.grid = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 0.35, face = "bold"),
      legend.text = element_text(size = base_size - 0.75),
      legend.key.height = grid::unit(2.4, "mm"),
      legend.key.width = grid::unit(3.2, "mm"),
      legend.spacing.x = grid::unit(0.8, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.15, face = "bold", colour = COL[["ink"]]),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(1.6, 1.8, 1.6, 1.8, unit = "mm"),
      plot.tag = element_text(size = 8.2, face = "bold", colour = COL[["ink"]])
    )
}

tag_theme <- theme(
  plot.tag.position = c(0, 1),
  plot.tag = element_text(size = 8.2, face = "bold", family = JTM_FONT,
                          colour = COL[["ink"]])
)

clean_panel <- function(plot) {
  plot + labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL) +
    theme(
      text = element_text(family = JTM_FONT, colour = COL[["ink"]]),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.tag = element_text(size = 8.2, face = "bold", family = JTM_FONT,
                              colour = COL[["ink"]])
    )
}

workflow_icon_grob <- function(icon, accent) {
  ink <- unname(COL[["ink"]])
  neutral <- unname(COL[["neutral"]])
  pale <- scales::alpha(accent, 0.16)
  thin <- 0.85
  arrow_small <- grid::arrow(
    length = grid::unit(1.05, "mm"), type = "closed", angle = 24
  )
  line_gp <- grid::gpar(col = ink, lwd = thin, lineend = "round", linejoin = "round")
  accent_gp <- grid::gpar(col = accent, fill = accent, lwd = thin)
  open_accent_gp <- grid::gpar(col = accent, fill = "white", lwd = 1.05)

  switch(
    icon,
    discover = grid::grobTree(
      grid::circleGrob(x = 0.14, y = 0.68, r = 0.065, gp = accent_gp),
      grid::circleGrob(x = 0.14, y = 0.45, r = 0.065, gp = accent_gp),
      grid::circleGrob(x = 0.14, y = 0.22, r = 0.065, gp = accent_gp),
      grid::segmentsGrob(
        x0 = c(0.22, 0.22, 0.22), y0 = c(0.68, 0.45, 0.22),
        x1 = c(0.43, 0.43, 0.43), y1 = c(0.68, 0.45, 0.22),
        gp = grid::gpar(col = neutral, lwd = thin)
      ),
      grid::segmentsGrob(x0 = 0.52, y0 = 0.15, x1 = 0.52, y1 = 0.82, gp = line_gp),
      grid::segmentsGrob(x0 = 0.52, y0 = 0.15, x1 = 0.92, y1 = 0.15, gp = line_gp),
      grid::segmentsGrob(
        x0 = c(0.59, 0.70, 0.81), y0 = 0.15,
        x1 = c(0.59, 0.70, 0.81), y1 = c(0.43, 0.68, 0.55),
        gp = grid::gpar(col = accent, lwd = 3.0, lineend = "butt")
      )
    ),
    lock = grid::grobTree(
      grid::roundrectGrob(
        x = 0.31, y = 0.50, width = 0.42, height = 0.70, r = grid::unit(1.1, "mm"),
        gp = grid::gpar(col = ink, fill = "white", lwd = thin)
      ),
      grid::segmentsGrob(
        x0 = c(0.18, 0.18, 0.18), y0 = c(0.67, 0.50, 0.33),
        x1 = c(0.43, 0.43, 0.43), y1 = c(0.67, 0.50, 0.33),
        gp = grid::gpar(col = accent, lwd = 1.25, lineend = "round")
      ),
      grid::circleGrob(x = 0.74, y = 0.63, r = 0.15, gp = open_accent_gp),
      grid::rectGrob(
        x = 0.74, y = 0.41, width = 0.33, height = 0.36,
        gp = grid::gpar(col = accent, fill = pale, lwd = 1.05)
      ),
      grid::circleGrob(x = 0.74, y = 0.45, r = 0.033, gp = accent_gp),
      grid::segmentsGrob(x0 = 0.74, y0 = 0.42, x1 = 0.74, y1 = 0.32,
                         gp = grid::gpar(col = accent, lwd = 1.0))
    ),
    validate = grid::grobTree(
      grid::segmentsGrob(x0 = 0.18, y0 = 0.57, x1 = 0.51, y1 = 0.57,
                         gp = grid::gpar(col = neutral, lwd = 1.0)),
      grid::circleGrob(x = 0.18, y = 0.57, r = 0.11, gp = open_accent_gp),
      grid::circleGrob(x = 0.51, y = 0.57, r = 0.11, gp = open_accent_gp),
      grid::segmentsGrob(x0 = 0.18, y0 = 0.46, x1 = 0.18, y1 = 0.27, gp = line_gp),
      grid::segmentsGrob(x0 = 0.51, y0 = 0.46, x1 = 0.51, y1 = 0.27, gp = line_gp),
      grid::circleGrob(x = 0.79, y = 0.48, r = 0.18,
                       gp = grid::gpar(col = accent, fill = accent, lwd = 1.0)),
      grid::polylineGrob(
        x = c(0.69, 0.77, 0.90), y = c(0.48, 0.37, 0.61),
        gp = grid::gpar(col = "white", lwd = 1.7, lineend = "round", linejoin = "round")
      )
    ),
    transport = grid::grobTree(
      grid::segmentsGrob(
        x0 = c(0.50, 0.50, 0.50, 0.50), y0 = c(0.50, 0.50, 0.50, 0.50),
        x1 = c(0.18, 0.36, 0.66, 0.82), y1 = c(0.25, 0.77, 0.77, 0.25),
        gp = grid::gpar(col = neutral, lwd = thin)
      ),
      grid::circleGrob(x = 0.50, y = 0.50, r = 0.12, gp = accent_gp),
      grid::circleGrob(x = 0.18, y = 0.25, r = 0.075, gp = open_accent_gp),
      grid::circleGrob(x = 0.36, y = 0.77, r = 0.075, gp = open_accent_gp),
      grid::circleGrob(x = 0.66, y = 0.77, r = 0.075, gp = open_accent_gp),
      grid::circleGrob(x = 0.82, y = 0.25, r = 0.075, gp = open_accent_gp)
    ),
    replicate = grid::grobTree(
      grid::segmentsGrob(
        x0 = c(0.12, 0.12, 0.12), y0 = c(0.74, 0.50, 0.26),
        x1 = c(0.60, 0.60, 0.60), y1 = c(0.74, 0.50, 0.26),
        gp = grid::gpar(col = neutral, lwd = 1.15, lineend = "round")
      ),
      grid::circleGrob(x = c(0.20, 0.34, 0.48), y = 0.74, r = 0.045, gp = accent_gp),
      grid::circleGrob(x = c(0.20, 0.34, 0.48), y = 0.50, r = 0.045, gp = accent_gp),
      grid::circleGrob(x = c(0.20, 0.34, 0.48), y = 0.26, r = 0.045, gp = accent_gp),
      grid::segmentsGrob(x0 = 0.65, y0 = 0.50, x1 = 0.88, y1 = 0.50,
                         arrow = arrow_small, gp = grid::gpar(col = accent, lwd = 1.2))
    ),
    pool = grid::grobTree(
      grid::segmentsGrob(x0 = 0.48, y0 = 0.12, x1 = 0.48, y1 = 0.86,
                         gp = grid::gpar(col = neutral, lwd = 0.7, lty = 2)),
      grid::segmentsGrob(
        x0 = c(0.16, 0.29, 0.37), y0 = c(0.72, 0.52, 0.32),
        x1 = c(0.63, 0.77, 0.86), y1 = c(0.72, 0.52, 0.32),
        gp = grid::gpar(col = ink, lwd = 1.15, lineend = "round")
      ),
      grid::circleGrob(x = c(0.40, 0.55, 0.63), y = c(0.72, 0.52, 0.32),
                       r = 0.045, gp = accent_gp),
      grid::polygonGrob(
        x = c(0.48, 0.58, 0.68, 0.58), y = c(0.13, 0.21, 0.13, 0.05),
        gp = grid::gpar(col = accent, fill = pale, lwd = 1.0)
      )
    ),
    ffpe = grid::grobTree(
      grid::polygonGrob(
        x = c(0.12, 0.63, 0.78, 0.27), y = c(0.28, 0.28, 0.73, 0.73),
        gp = grid::gpar(col = ink, fill = pale, lwd = 0.95)
      ),
      grid::circleGrob(x = 0.45, y = 0.50, r = 0.12,
                       gp = grid::gpar(col = accent, fill = scales::alpha(accent, 0.55), lwd = 0.8)),
      grid::rectGrob(x = 0.78, y = 0.34, width = 0.27, height = 0.52,
                     gp = grid::gpar(col = accent, fill = "white", lwd = 1.0)),
      grid::segmentsGrob(x0 = 0.70, y0 = 0.50, x1 = 0.86, y1 = 0.50,
                         gp = grid::gpar(col = accent, lwd = 1.25))
    ),
    expand = grid::grobTree(
      grid::segmentsGrob(x0 = 0.18, y0 = 0.50, x1 = 0.43, y1 = 0.50,
                         arrow = arrow_small, gp = grid::gpar(col = accent, lwd = 1.1)),
      grid::segmentsGrob(
        x0 = c(0.43, 0.43, 0.43), y0 = c(0.50, 0.50, 0.50),
        x1 = c(0.78, 0.78, 0.78), y1 = c(0.76, 0.50, 0.24),
        gp = grid::gpar(col = neutral, lwd = thin)
      ),
      grid::circleGrob(x = 0.18, y = 0.50, r = 0.085, gp = accent_gp),
      grid::circleGrob(x = 0.78, y = 0.76, r = 0.07, gp = open_accent_gp),
      grid::circleGrob(x = 0.78, y = 0.50, r = 0.07, gp = open_accent_gp),
      grid::circleGrob(x = 0.78, y = 0.24, r = 0.07, gp = open_accent_gp)
    ),
    apc = grid::grobTree(
      grid::circleGrob(x = 0.24, y = 0.50, r = 0.20, gp = open_accent_gp),
      grid::circleGrob(x = c(0.15, 0.24, 0.33, 0.20, 0.29),
                       y = c(0.50, 0.62, 0.49, 0.38, 0.37), r = 0.045, gp = accent_gp),
      grid::segmentsGrob(x0 = 0.46, y0 = 0.50, x1 = 0.70, y1 = 0.50,
                         arrow = arrow_small, gp = grid::gpar(col = neutral, lwd = 1.1)),
      grid::circleGrob(x = 0.82, y = 0.50, r = 0.12, gp = open_accent_gp),
      grid::segmentsGrob(x0 = 0.75, y0 = 0.43, x1 = 0.89, y1 = 0.57,
                         gp = grid::gpar(col = accent, lwd = 1.45))
    ),
    wnt = grid::grobTree(
      grid::segmentsGrob(x0 = 0.08, y0 = 0.28, x1 = 0.92, y1 = 0.28,
                         gp = grid::gpar(col = neutral, lwd = 1.1)),
      grid::rectGrob(x = 0.22, y = 0.36, width = 0.13, height = 0.24,
                     gp = grid::gpar(col = accent, fill = pale, lwd = 1.0)),
      grid::segmentsGrob(x0 = 0.22, y0 = 0.48, x1 = 0.54, y1 = 0.62,
                         arrow = arrow_small, gp = grid::gpar(col = accent, lwd = 1.1)),
      grid::circleGrob(x = 0.74, y = 0.63, r = 0.18,
                       gp = grid::gpar(col = ink, fill = "white", lwd = 0.95)),
      grid::circleGrob(x = 0.74, y = 0.63, r = 0.07, gp = accent_gp)
    ),
    drug = grid::grobTree(
      grid::roundrectGrob(x = 0.29, y = 0.58, width = 0.38, height = 0.20,
                          r = grid::unit(2.3, "mm"),
                          gp = grid::gpar(col = accent, fill = "white", lwd = 1.05)),
      grid::segmentsGrob(x0 = 0.29, y0 = 0.48, x1 = 0.29, y1 = 0.68,
                         gp = grid::gpar(col = accent, lwd = 1.0)),
      grid::rectGrob(x = 0.71, y = 0.45, width = 0.36, height = 0.50,
                     gp = grid::gpar(col = ink, fill = "white", lwd = thin)),
      grid::circleGrob(x = c(0.62, 0.72, 0.80, 0.62, 0.72, 0.80),
                       y = c(0.56, 0.56, 0.56, 0.35, 0.35, 0.35), r = 0.035,
                       gp = accent_gp)
    ),
    virtual = grid::grobTree(
      grid::segmentsGrob(
        x0 = c(0.20, 0.20, 0.50, 0.50, 0.50),
        y0 = c(0.50, 0.50, 0.50, 0.50, 0.50),
        x1 = c(0.50, 0.50, 0.80, 0.80, 0.80),
        y1 = c(0.72, 0.28, 0.74, 0.50, 0.26),
        gp = grid::gpar(col = neutral, lwd = thin)
      ),
      grid::circleGrob(x = 0.20, y = 0.50, r = 0.07, gp = open_accent_gp),
      grid::circleGrob(x = 0.50, y = 0.72, r = 0.07, gp = open_accent_gp),
      grid::circleGrob(x = 0.50, y = 0.28, r = 0.07, gp = open_accent_gp),
      grid::circleGrob(x = 0.80, y = 0.74, r = 0.07, gp = open_accent_gp),
      grid::circleGrob(x = 0.80, y = 0.50, r = 0.07, gp = open_accent_gp),
      grid::circleGrob(x = 0.80, y = 0.26, r = 0.07, gp = open_accent_gp),
      grid::segmentsGrob(x0 = 0.43, y0 = 0.21, x1 = 0.57, y1 = 0.35,
                         gp = grid::gpar(col = accent, lwd = 1.45)),
      grid::segmentsGrob(x0 = 0.43, y0 = 0.35, x1 = 0.57, y1 = 0.21,
                         gp = grid::gpar(col = accent, lwd = 1.45))
    ),
    grid::grobTree(
      grid::circleGrob(x = 0.5, y = 0.5, r = 0.20, gp = open_accent_gp)
    )
  )
}

workflow_panel <- function(data) {
  links <- data.frame(
    x = head(data$x, -1) + 0.42,
    xend = tail(data$x, -1) - 0.42,
    y = 0,
    yend = 0
  )
  plot <- ggplot(data, aes(x, 0)) +
    geom_segment(
      data = links,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, linewidth = 0.50, colour = COL[["neutral"]],
      arrow = arrow(length = grid::unit(1.35, "mm"), type = "closed", angle = 23)
    ) +
    scale_colour_identity() +
    coord_cartesian(
      xlim = range(data$x) + c(-0.55, 0.55), ylim = c(-0.48, 0.48), clip = "off"
    ) +
    theme_void(base_size = 7, base_family = JTM_FONT) +
    theme(
      text = element_text(family = JTM_FONT, colour = COL[["ink"]]),
      plot.margin = margin(1.5, 3.0, 0.9, 3.0, "mm"),
      plot.tag = element_text(size = 8.2, face = "bold", family = JTM_FONT)
    )

  for (i in seq_len(nrow(data))) {
    card <- grid::roundrectGrob(
      r = grid::unit(1.5, "mm"),
      gp = grid::gpar(
        col = COL[["neutral_light"]], fill = data$fill[i], lwd = 0.85,
        linejoin = "round"
      )
    )
    icon_well <- grid::roundrectGrob(
      r = grid::unit(1.3, "mm"),
      gp = grid::gpar(
        col = scales::alpha(data$accent[i], 0.34), fill = "white", lwd = 0.75
      )
    )
    plot <- plot +
      annotation_custom(card, xmin = data$x[i] - 0.42, xmax = data$x[i] + 0.42,
                        ymin = -0.43, ymax = 0.43) +
      annotation_custom(icon_well, xmin = data$x[i] - 0.20, xmax = data$x[i] + 0.20,
                        ymin = 0.025, ymax = 0.355) +
      annotation_custom(
        workflow_icon_grob(data$icon[i], data$accent[i]),
        xmin = data$x[i] - 0.165, xmax = data$x[i] + 0.165,
        ymin = 0.065, ymax = 0.325
      )
  }

  # Annotation grobs are added after the ggplot layers; redraw text and step
  # badges on top so all exported formats retain identical hierarchy.
  plot +
    geom_segment(
      aes(x = x - 0.30, xend = x + 0.30, y = 0.415, yend = 0.415, colour = accent),
      linewidth = 1.05, lineend = "round"
    ) +
    geom_point(
      aes(x = x - 0.355, y = 0.345), inherit.aes = FALSE,
      shape = 21, size = 2.45, stroke = 0.34, fill = COL[["ink"]], colour = "white"
    ) +
    geom_text(
      aes(x = x - 0.355, y = 0.345, label = step), size = 1.50,
      fontface = "bold", colour = "white", family = JTM_FONT
    ) +
    geom_text(
      aes(y = -0.105, label = stage), size = 1.82, fontface = "bold",
      colour = COL[["ink"]], family = JTM_FONT
    ) +
    geom_text(
      aes(y = -0.285, label = detail), size = 1.48, lineheight = 0.92,
      colour = COL[["ink"]], family = JTM_FONT
    )
}

study_overview_panel <- function(branches) {
  branch_links_left <- branches %>%
    transmute(x = 5.30, xend = 5.62, y = y, yend = y)
  branch_links_right <- branches %>%
    transmute(x = 9.58, xend = 9.92, y = y, yend = y)

  plot <- ggplot() +
    geom_segment(
      aes(x = 2.25, xend = 2.78, y = 3.50, yend = 3.50),
      linewidth = 0.52, colour = COL[["neutral"]],
      arrow = arrow(length = grid::unit(1.35, "mm"), type = "closed", angle = 23)
    ) +
    geom_segment(
      aes(x = 5.05, xend = 5.30, y = 3.50, yend = 3.50),
      linewidth = 0.52, colour = COL[["neutral"]]
    ) +
    geom_segment(
      aes(x = 5.30, xend = 5.30, y = min(branches$y), yend = max(branches$y)),
      linewidth = 0.42, colour = COL[["neutral"]]
    ) +
    geom_segment(
      data = branch_links_left,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.46, colour = COL[["neutral"]],
      arrow = arrow(length = grid::unit(1.15, "mm"), type = "closed", angle = 23)
    ) +
    geom_segment(
      data = branch_links_right,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.42, colour = COL[["neutral"]]
    ) +
    geom_segment(
      aes(x = 9.92, xend = 9.92, y = min(branches$y), yend = max(branches$y)),
      linewidth = 0.42, colour = COL[["neutral"]]
    ) +
    geom_segment(
      aes(x = 9.92, xend = 10.23, y = 3.50, yend = 3.50),
      linewidth = 0.52, colour = COL[["neutral"]],
      arrow = arrow(length = grid::unit(1.35, "mm"), type = "closed", angle = 23)
    ) +
    coord_cartesian(xlim = c(0.05, 12.75), ylim = c(0.25, 6.75), clip = "off") +
    theme_void(base_size = 7, base_family = JTM_FONT) +
    theme(
      text = element_text(family = JTM_FONT, colour = COL[["ink"]]),
      plot.margin = margin(1.8, 2.5, 1.2, 2.5, "mm"),
      plot.tag = element_text(size = 8.2, face = "bold", family = JTM_FONT)
    )

  major_cards <- list(
    list(xmin = 0.15, xmax = 2.25, ymin = 2.12, ymax = 4.88,
         fill = "#F5F8FA", border = COL[["wnt"]]),
    list(xmin = 2.78, xmax = 5.05, ymin = 2.12, ymax = 4.88,
         fill = "#FCF6F3", border = COL[["route"]]),
    list(xmin = 10.23, xmax = 12.70, ymin = 1.78, ymax = 5.22,
         fill = "#F4F8F7", border = COL[["crc"]])
  )
  for (card in major_cards) {
    plot <- plot + annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(1.8, "mm"),
        gp = grid::gpar(
          col = scales::alpha(card$border, 0.48), fill = card$fill, lwd = 1.0
        )
      ),
      xmin = card$xmin, xmax = card$xmax, ymin = card$ymin, ymax = card$ymax
    )
  }

  plot <- plot +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(1.4, "mm"),
        gp = grid::gpar(col = scales::alpha(COL[["wnt"]], 0.34), fill = "white", lwd = 0.8)
      ),
      xmin = 0.70, xmax = 1.70, ymin = 3.62, ymax = 4.55
    ) +
    annotation_custom(
      workflow_icon_grob("discover", COL[["wnt"]]),
      xmin = 0.84, xmax = 1.56, ymin = 3.78, ymax = 4.40
    ) +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(1.4, "mm"),
        gp = grid::gpar(col = scales::alpha(COL[["route"]], 0.34), fill = "white", lwd = 0.8)
      ),
      xmin = 3.42, xmax = 4.42, ymin = 3.62, ymax = 4.55
    ) +
    annotation_custom(
      workflow_icon_grob("lock", COL[["route"]]),
      xmin = 3.56, xmax = 4.28, ymin = 3.78, ymax = 4.40
    )

  for (i in seq_len(nrow(branches))) {
    card <- grid::roundrectGrob(
      r = grid::unit(1.35, "mm"),
      gp = grid::gpar(
        col = scales::alpha(branches$accent[i], 0.42),
        fill = branches$fill[i], lwd = 0.82
      )
    )
    icon_well <- grid::roundrectGrob(
      r = grid::unit(1.0, "mm"),
      gp = grid::gpar(
        col = scales::alpha(branches$accent[i], 0.30), fill = "white", lwd = 0.68
      )
    )
    plot <- plot +
      annotation_custom(card, xmin = 5.62, xmax = 9.58,
                        ymin = branches$y[i] - 0.46, ymax = branches$y[i] + 0.46) +
      annotation_custom(icon_well, xmin = 5.76, xmax = 6.40,
                        ymin = branches$y[i] - 0.31, ymax = branches$y[i] + 0.31) +
      annotation_custom(
        workflow_icon_grob(branches$icon[i], branches$accent[i]),
        xmin = 5.84, xmax = 6.32,
        ymin = branches$y[i] - 0.23, ymax = branches$y[i] + 0.23
      )
  }

  plot +
    annotate("segment", x = 0.55, xend = 1.85, y = 4.78, yend = 4.78,
             colour = COL[["wnt"]], linewidth = 1.05, lineend = "round") +
    annotate("segment", x = 3.15, xend = 4.68, y = 4.78, yend = 4.78,
             colour = COL[["route"]], linewidth = 1.05, lineend = "round") +
    annotate("segment", x = 10.62, xend = 12.31, y = 5.12, yend = 5.12,
             colour = COL[["crc"]], linewidth = 1.05, lineend = "round") +
    annotate("text", x = 1.20, y = 3.24, label = "DEFINE IN ADENOMA",
             size = 2.05, fontface = "bold", family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 1.20, y = 2.63,
             label = "Chen donor medians\n1,000 donor bootstraps",
             size = 1.62, lineheight = 0.94, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 3.92, y = 3.24, label = "FIXED 100-GENE STATE",
             size = 2.05, fontface = "bold", family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 3.92, y = 2.66,
             label = "50 up + 50 down\nno reselection · equal weights",
             size = 1.62, lineheight = 0.94, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("label", x = 3.92, y = 2.24, label = "validation data unseen",
             size = 1.35, family = JTM_FONT, colour = COL[["route"]], fill = "white",
             linewidth = 0.20, label.padding = grid::unit(0.9, "mm")) +
    geom_text(
      data = branches, aes(x = 6.56, y = y + 0.16, label = title, colour = accent),
      hjust = 0, size = 1.75, fontface = "bold", family = JTM_FONT,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = branches, aes(x = 6.56, y = y - 0.18, label = detail),
      hjust = 0, size = 1.42, family = JTM_FONT, colour = COL[["ink"]],
      inherit.aes = FALSE
    ) +
    geom_label(
      data = branches, aes(x = 9.20, y = y, label = figure, colour = accent),
      size = 1.25, fontface = "bold", family = JTM_FONT, fill = "white",
      linewidth = 0.18, label.padding = grid::unit(0.75, "mm"),
      inherit.aes = FALSE
    ) +
    scale_colour_identity() +
    annotate("text", x = 11.46, y = 4.53, label = "BIOLOGICAL CONCLUSION",
             size = 1.72, fontface = "bold", family = JTM_FONT, colour = COL[["crc"]]) +
    annotate("text", x = 11.46, y = 3.73,
             label = "Pre-invasive adenoma\nepithelial state",
             size = 1.75, fontface = "bold", lineheight = 0.94,
             family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("segment", x = 10.66, xend = 10.94, y = 3.06, yend = 3.06,
             colour = COL[["route"]], linewidth = 1.05, lineend = "round") +
    annotate("text", x = 11.02, y = 3.06, hjust = 0,
             label = "WNT / stem–progenitor ↑",
             size = 1.35, fontface = "bold", family = JTM_FONT, colour = COL[["route"]]) +
    annotate("segment", x = 10.66, xend = 10.94, y = 2.61, yend = 2.61,
             colour = COL[["wnt"]], linewidth = 1.05, lineend = "round") +
    annotate("text", x = 11.02, y = 2.61, hjust = 0,
             label = "Mature differentiation ↓",
             size = 1.35, fontface = "bold", family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 11.46, y = 2.13,
             label = "reproduced across independent evidence layers",
             size = 1.20, family = JTM_FONT, colour = COL[["neutral"]])
}

mechanism_state_icon_grob <- function(type, accent) {
  pale <- scales::alpha(accent, 0.18)
  nucleus_gp <- grid::gpar(col = accent, fill = "white", lwd = 0.62)
  cell_gp <- grid::gpar(col = accent, fill = pale, lwd = 0.78)
  if (type == "stem") {
    return(grid::grobTree(
      grid::circleGrob(x = 0.25, y = 0.35, r = 0.19, gp = cell_gp),
      grid::circleGrob(x = 0.51, y = 0.62, r = 0.19, gp = cell_gp),
      grid::circleGrob(x = 0.76, y = 0.35, r = 0.19, gp = cell_gp),
      grid::circleGrob(x = 0.25, y = 0.35, r = 0.055, gp = nucleus_gp),
      grid::circleGrob(x = 0.51, y = 0.62, r = 0.055, gp = nucleus_gp),
      grid::circleGrob(x = 0.76, y = 0.35, r = 0.055, gp = nucleus_gp)
    ))
  }
  grid::grobTree(
    grid::roundrectGrob(x = 0.23, y = 0.46, width = 0.24, height = 0.68,
                       r = grid::unit(0.8, "mm"), gp = cell_gp),
    grid::roundrectGrob(x = 0.51, y = 0.46, width = 0.24, height = 0.68,
                       r = grid::unit(0.8, "mm"), gp = cell_gp),
    grid::roundrectGrob(x = 0.79, y = 0.46, width = 0.24, height = 0.68,
                       r = grid::unit(0.8, "mm"), gp = cell_gp),
    grid::circleGrob(x = 0.23, y = 0.35, r = 0.050, gp = nucleus_gp),
    grid::circleGrob(x = 0.51, y = 0.35, r = 0.050, gp = nucleus_gp),
    grid::circleGrob(x = 0.79, y = 0.35, r = 0.050, gp = nucleus_gp),
    grid::segmentsGrob(x0 = c(0.16, 0.44, 0.72), y0 = 0.82,
                       x1 = c(0.30, 0.58, 0.86), y1 = 0.82,
                       gp = grid::gpar(col = accent, lwd = 0.85, lineend = "round"))
  )
}

canonical_receptor_grob <- function() {
  helix_y <- seq(0.47, 0.83, length.out = 7)
  loop_side <- rep(c(0.54, 0.29), length.out = 6)
  grid::grobTree(
    # Frizzled: extracellular cysteine-rich domain plus a seven-pass bundle.
    grid::circleGrob(
      x = 0.15, y = 0.59, r = 0.075,
      gp = grid::gpar(
        col = COL[["wnt"]], fill = scales::alpha(COL[["wnt"]], 0.10), lwd = 0.95
      )
    ),
    grid::segmentsGrob(
      x0 = 0.21, y0 = 0.62, x1 = 0.31, y1 = 0.83,
      gp = grid::gpar(col = COL[["wnt"]], lwd = 0.85, lineend = "round")
    ),
    grid::segmentsGrob(
      x0 = rep(0.31, 7), x1 = rep(0.52, 7),
      y0 = helix_y, y1 = helix_y,
      gp = grid::gpar(col = COL[["wnt"]], lwd = 1.20, lineend = "round")
    ),
    grid::segmentsGrob(
      x0 = loop_side, x1 = loop_side,
      y0 = helix_y[-7], y1 = helix_y[-1],
      gp = grid::gpar(col = COL[["wnt"]], lwd = 0.62, lineend = "round")
    ),
    grid::segmentsGrob(
      x0 = 0.52, y0 = 0.47, x1 = 0.64, y1 = 0.43,
      gp = grid::gpar(col = COL[["wnt"]], lwd = 0.80, lineend = "round")
    ),
    # LRP5/6: extracellular repeat domains, one transmembrane helix and a
    # phosphorylated cytoplasmic tail. It is part of the same receptor complex.
    grid::circleGrob(
      x = c(0.10, 0.18, 0.26), y = c(0.30, 0.30, 0.30),
      r = c(0.064, 0.071, 0.064),
      gp = grid::gpar(
        col = COL[["wnt"]], fill = scales::alpha(COL[["wnt_light"]], 0.30), lwd = 0.78
      )
    ),
    grid::segmentsGrob(
      x0 = 0.31, x1 = 0.52, y0 = 0.30, y1 = 0.30,
      gp = grid::gpar(col = COL[["wnt"]], lwd = 1.65, lineend = "round")
    ),
    grid::segmentsGrob(
      x0 = 0.52, x1 = 0.84, y0 = 0.30, y1 = 0.30,
      gp = grid::gpar(col = COL[["wnt"]], lwd = 0.86, lineend = "round")
    ),
    grid::circleGrob(
      x = c(0.64, 0.74, 0.84), y = rep(0.30, 3), r = 0.031,
      gp = grid::gpar(col = COL[["wnt"]], fill = "white", lwd = 0.67)
    ),
    grid::textGrob(
      "P", x = c(0.64, 0.74, 0.84), y = rep(0.30, 3),
      gp = grid::gpar(
        col = COL[["wnt"]], fontsize = 4.5, fontface = "bold", fontfamily = JTM_FONT
      )
    )
  )
}

canonical_proteasome_grob <- function() {
  grid::grobTree(
    grid::roundrectGrob(
      x = 0.50, y = 0.50, width = 0.66, height = 0.50,
      r = grid::unit(1.3, "mm"),
      gp = grid::gpar(
        col = COL[["neutral"]], fill = scales::alpha(COL[["neutral_light"]], 0.46), lwd = 0.85
      )
    ),
    grid::segmentsGrob(
      x0 = c(0.24, 0.36, 0.48, 0.60, 0.72), y0 = 0.30,
      x1 = c(0.24, 0.36, 0.48, 0.60, 0.72), y1 = 0.70,
      gp = grid::gpar(col = COL[["neutral"]], lwd = 0.50)
    ),
    grid::segmentsGrob(
      x0 = 0.20, x1 = 0.80, y0 = c(0.37, 0.63), y1 = c(0.37, 0.63),
      gp = grid::gpar(col = COL[["neutral"]], lwd = 0.65)
    ),
    grid::textGrob(
      "26S", x = 0.50, y = 0.50,
      gp = grid::gpar(col = COL[["ink"]], fontsize = 5.6, fontface = "bold", fontfamily = JTM_FONT)
    )
  )
}

mechanism_panel <- function(interventions, evidence_strip) {
  flow_arrow <- grid::arrow(
    length = grid::unit(1.25, "mm"), type = "closed", angle = 24
  )
  small_arrow <- grid::arrow(
    length = grid::unit(1.05, "mm"), type = "closed", angle = 24
  )
  plot <- ggplot() +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(2.1, "mm"),
        gp = grid::gpar(col = COL[["neutral_light"]], fill = "#FAFBFB", lwd = 0.90)
      ),
      xmin = 0.15, xmax = 13.05, ymin = 2.55, ymax = 8.00
    ) +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(2.1, "mm"),
        gp = grid::gpar(
          col = scales::alpha(COL[["route"]], 0.50), fill = "#FCF6F3", lwd = 1.00
        )
      ),
      xmin = 13.25, xmax = 15.90, ymin = 2.55, ymax = 8.00
    ) +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(1.45, "mm"),
        gp = grid::gpar(
          col = scales::alpha(COL[["route"]], 0.45), fill = "#FFF9F6", lwd = 0.82
        )
      ),
      xmin = 10.80, xmax = 12.88, ymin = 5.42, ymax = 6.54
    ) +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(1.45, "mm"),
        gp = grid::gpar(
          col = scales::alpha(COL[["wnt"]], 0.45), fill = "#F6F9FB", lwd = 0.82
        )
      ),
      xmin = 10.80, xmax = 12.88, ymin = 4.05, ymax = 5.17
    ) +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(1.6, "mm"),
        gp = grid::gpar(col = COL[["neutral_light"]], fill = "#F7F8F8", lwd = 0.82)
      ),
      xmin = 0.15, xmax = 15.90, ymin = 0.05, ymax = 1.12
    ) +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(10, "mm"),
        gp = grid::gpar(
          col = scales::alpha(COL[["context"]], 0.64), fill = "white", lwd = 0.98
        )
      ),
      xmin = 8.45, xmax = 10.62, ymin = 4.58, ymax = 7.15
    ) +
    annotation_custom(
      mechanism_state_icon_grob("stem", COL[["route"]]),
      xmin = 10.96, xmax = 11.47, ymin = 5.66, ymax = 6.30
    ) +
    annotation_custom(
      mechanism_state_icon_grob("differentiated", COL[["wnt"]]),
      xmin = 10.96, xmax = 11.47, ymin = 4.29, ymax = 4.93
    ) +
    annotation_custom(
      canonical_receptor_grob(), xmin = 1.82, xmax = 2.82, ymin = 5.03, ymax = 6.48
    ) +
    annotation_custom(
      canonical_proteasome_grob(), xmin = 8.88, xmax = 9.72, ymin = 3.35, ymax = 4.14
    ) +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(1.2, "mm"),
        gp = grid::gpar(
          col = scales::alpha(COL[["neutral"]], 0.64), fill = "white", lwd = 0.82
        )
      ),
      xmin = 4.05, xmax = 6.05, ymin = 3.22, ymax = 4.52
    ) +
    geom_segment(
      data = interventions,
      aes(x = target_x, xend = x, y = target_y, yend = anchor_y, colour = colour),
      linewidth = 0.40, linetype = "dotted", inherit.aes = FALSE
    ) +
    geom_segment(aes(x = 2.16, xend = 2.16, y = 2.95, yend = 7.38),
                 linewidth = 1.20, colour = COL[["neutral_light"]]) +
    geom_segment(aes(x = 2.25, xend = 2.25, y = 2.95, yend = 7.38),
                 linewidth = 0.58, colour = "white") +
    geom_segment(aes(x = 1.10, xend = 1.78, y = 5.80, yend = 5.80),
                 linewidth = 0.68, colour = COL[["wnt"]], arrow = flow_arrow) +
    geom_segment(aes(x = 2.82, xend = 3.24, y = 5.80, yend = 5.80),
                 linewidth = 0.68, colour = COL[["wnt"]], arrow = flow_arrow) +
    geom_segment(aes(x = 3.72, xend = 6.22, y = 5.80, yend = 5.80),
                 linewidth = 0.68, colour = COL[["wnt"]], arrow = flow_arrow) +
    geom_segment(aes(x = 7.52, xend = 8.38, y = 5.80, yend = 5.80),
                 linewidth = 0.68, colour = COL[["wnt"]], arrow = flow_arrow) +
    geom_segment(aes(x = 3.48, xend = 3.48, y = 5.49, yend = 4.83),
                 linewidth = 0.50, colour = COL[["wnt"]]) +
    geom_segment(aes(x = 3.48, xend = 5.05, y = 4.83, yend = 4.83),
                 linewidth = 0.50, colour = COL[["wnt"]]) +
    geom_segment(aes(x = 5.05, xend = 5.05, y = 4.83, yend = 4.63),
                 linewidth = 0.50, colour = COL[["wnt"]]) +
    geom_segment(aes(x = 4.78, xend = 5.32, y = 4.61, yend = 4.61),
                 linewidth = 1.05, colour = COL[["wnt"]], lineend = "round") +
    geom_segment(aes(x = 6.06, xend = 6.55, y = 3.82, yend = 3.82),
                 linewidth = 0.55, colour = COL[["neutral"]], arrow = small_arrow) +
    geom_segment(aes(x = 7.18, xend = 7.54, y = 3.82, yend = 3.82),
                 linewidth = 0.55, colour = COL[["neutral"]], arrow = small_arrow) +
    geom_segment(aes(x = 8.25, xend = 8.82, y = 3.82, yend = 3.82),
                 linewidth = 0.55, colour = COL[["neutral"]], arrow = small_arrow) +
    geom_segment(aes(x = 10.63, xend = 10.78, y = 5.97, yend = 5.97),
                 linewidth = 0.58, colour = COL[["route"]], arrow = small_arrow) +
    geom_segment(aes(x = 10.63, xend = 10.78, y = 4.60, yend = 4.60),
                 linewidth = 0.58, colour = COL[["wnt"]], arrow = small_arrow) +
    geom_segment(aes(x = 12.88, xend = 13.08, y = 5.97, yend = 5.97),
                 linewidth = 0.48, colour = COL[["route"]]) +
    geom_segment(aes(x = 12.88, xend = 13.08, y = 4.60, yend = 4.60),
                 linewidth = 0.48, colour = COL[["wnt"]]) +
    geom_segment(aes(x = 13.08, xend = 13.08, y = 4.60, yend = 5.97),
                 linewidth = 0.48, colour = COL[["neutral"]]) +
    geom_segment(aes(x = 13.08, xend = 13.24, y = 5.28, yend = 5.28),
                 linewidth = 0.65, colour = COL[["ink"]], arrow = flow_arrow) +
    geom_point(aes(x = 0.72, y = 5.80), shape = 21, size = 4.9, stroke = 0.72,
               fill = scales::alpha(COL[["wnt"]], 0.18), colour = COL[["wnt"]]) +
    geom_point(aes(x = 3.48, y = 5.80), shape = 23, size = 4.15, stroke = 0.76,
               fill = "white", colour = COL[["wnt"]]) +
    geom_point(
      data = data.frame(
        x = c(6.48, 6.86, 7.22),
        y = c(5.72, 5.95, 5.72)
      ),
      aes(x, y), inherit.aes = FALSE, shape = 21, size = 2.85, stroke = 0.60,
      fill = scales::alpha(COL[["wnt"]], 0.16), colour = COL[["wnt"]]
    ) +
    geom_point(aes(x = 6.86, y = 3.82), shape = 21, size = 3.55, stroke = 0.62,
               fill = "white", colour = COL[["neutral"]]) +
    geom_point(
      data = data.frame(x = c(7.70, 7.93, 8.16), y = c(3.70, 3.91, 3.70)),
      aes(x, y), inherit.aes = FALSE, shape = 21, size = 2.18, stroke = 0.50,
      fill = scales::alpha(COL[["neutral_light"]], 0.55), colour = COL[["neutral"]]
    ) +
    annotate("text", x = 0.72, y = 5.80, label = "WNT",
             size = 1.46, fontface = "bold", family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 1.80, y = 4.92, label = "FZD",
             size = 1.07, fontface = "bold", family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 2.66, y = 4.92, label = "LRP5/6",
             size = 1.07, fontface = "bold", family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 3.48, y = 5.80, label = "DVL",
             size = 1.28, fontface = "bold", family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 6.86, y = 6.40, label = "stabilised β-catenin",
             size = 1.63, fontface = "bold", family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 7.70, y = 6.08, label = "nuclear entry",
             size = 1.12, family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 4.25, y = 5.05, label = "inhibits complex",
             size = 1.10, family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 4.62, y = 6.14, label = "p-LRP5/6 recruits AXIN",
             size = 1.00, family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 5.05, y = 4.29, label = "β-catenin\ndestruction complex",
             size = 1.15, lineheight = 0.88, fontface = "bold",
             family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("label", x = 4.56, y = 3.83, label = "APC", size = 1.05,
             fontface = "bold", family = JTM_FONT, colour = COL[["route"]], fill = "#FFF8F4",
             linewidth = 0.18, label.padding = grid::unit(0.65, "mm")) +
    annotate("label", x = 5.47, y = 3.83, label = "AXIN", size = 1.05,
             fontface = "bold", family = JTM_FONT, colour = COL[["ink"]], fill = "white",
             linewidth = 0.18, label.padding = grid::unit(0.65, "mm")) +
    annotate("label", x = 4.56, y = 3.43, label = "GSK3β", size = 1.00,
             fontface = "bold", family = JTM_FONT, colour = COL[["ink"]], fill = "white",
             linewidth = 0.18, label.padding = grid::unit(0.58, "mm")) +
    annotate("label", x = 5.47, y = 3.43, label = "CK1", size = 1.05,
             fontface = "bold", family = JTM_FONT, colour = COL[["ink"]], fill = "white",
             linewidth = 0.18, label.padding = grid::unit(0.65, "mm")) +
    annotate("text", x = 6.42, y = 4.18, label = "CK1 / GSK3β\nphosphorylation",
             size = 0.98, lineheight = 0.88, family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 6.86, y = 3.82, label = "P",
             size = 1.05, fontface = "bold", family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 6.86, y = 3.24, label = "p–β-catenin",
             size = 1.08, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 7.93, y = 3.70, label = "Ub",
             size = 0.92, fontface = "bold", family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 7.93, y = 3.24, label = "β-TrCP / ubiquitin",
             size = 1.03, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 9.30, y = 3.15, label = "proteasomal\ndegradation",
             size = 1.02, lineheight = 0.90, family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 9.54, y = 6.86, label = "NUCLEUS",
             size = 1.18, fontface = "bold", family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 9.54, y = 6.14, label = "β-catenin + TCF7L2 / LEF1",
             size = 1.38, fontface = "bold", family = JTM_FONT, colour = COL[["context"]]) +
    annotate("segment", x = 8.90, xend = 10.18, y = 5.83, yend = 5.83,
             linewidth = 0.28, colour = COL[["neutral_light"]]) +
    annotate("text", x = 9.54, y = 5.42, label = "ASCL2 transcriptional programme",
             size = 1.25, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 11.65, y = 5.98, label = "Stem / progenitor\ngenes ↑",
             hjust = 0, size = 1.47, lineheight = 0.91, fontface = "bold",
             family = JTM_FONT, colour = COL[["route"]]) +
    annotate("text", x = 11.65, y = 4.61, label = "Differentiation\ngenes ↓",
             hjust = 0, size = 1.47, lineheight = 0.91, fontface = "bold",
             family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 14.58, y = 7.25, label = "LOCKED EPITHELIAL\nROUTE",
             size = 1.70, lineheight = 0.90, fontface = "bold", family = JTM_FONT,
             colour = COL[["route"]]) +
    annotate("segment", x = 13.67, xend = 15.49, y = 6.64, yend = 6.64,
             linewidth = 0.34, colour = COL[["route_light"]]) +
    annotate("segment", x = 13.62, xend = 14.10, y = 5.75, yend = 5.75,
             linewidth = 1.45, colour = COL[["route"]], lineend = "round") +
    annotate("text", x = 14.28, y = 5.75, label = "50-gene up arm", hjust = 0,
             size = 1.38, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("segment", x = 13.62, xend = 14.10, y = 5.00, yend = 5.00,
             linewidth = 1.45, colour = COL[["wnt"]], lineend = "round") +
    annotate("text", x = 14.28, y = 5.00, label = "50-gene down arm", hjust = 0,
             size = 1.38, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 14.58, y = 3.93, label = "integrated without refitting\nor context-specific reweighting",
             size = 1.18, lineheight = 0.95, family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 0.42, y = 7.68, label = "EXTRACELLULAR SIGNAL",
             hjust = 0,
             size = 1.28, fontface = "bold", family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 3.05, y = 7.68, label = "CYTOPLASMIC CONTROL",
             hjust = 0,
             size = 1.28, fontface = "bold", family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 8.55, y = 7.68, label = "TRANSCRIPTIONAL SWITCH",
             hjust = 0,
             size = 1.28, fontface = "bold", family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 2.47, y = 3.54, label = "cell membrane",
             angle = 90, size = 1.10, family = JTM_FONT, colour = COL[["neutral"]])

  plot <- plot +
    annotate("label", x = 0.92, y = 6.78, label = "WNT ON or APC loss",
             size = 1.08, fontface = "bold", family = JTM_FONT,
             colour = COL[["route"]], fill = "white", linewidth = 0.22,
             label.padding = grid::unit(0.72, "mm")) +
    annotate("label", x = 0.75, y = 3.84, label = "WNT OFF",
             size = 1.08, fontface = "bold", family = JTM_FONT,
             colour = COL[["wnt"]], fill = "white", linewidth = 0.22,
             label.padding = grid::unit(0.72, "mm")) +
    annotate("text", x = 0.75, y = 3.43, label = "destruction complex active",
             hjust = 0.5, size = 1.00, family = JTM_FONT, colour = COL[["neutral"]])

  for (i in seq_len(nrow(interventions))) {
    plot <- plot + annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(1.25, "mm"),
        gp = grid::gpar(
          col = scales::alpha(interventions$colour[i], 0.78), fill = "white", lwd = 0.82
        )
      ),
      xmin = interventions$xmin[i], xmax = interventions$xmax[i],
      ymin = interventions$ymin[i], ymax = interventions$ymax[i]
    )
  }

  plot +
    geom_point(
      data = interventions, aes(x = target_x, y = target_y, colour = colour),
      shape = 21, size = 1.45, stroke = 0.48, fill = "white", inherit.aes = FALSE
    ) +
    geom_text(
      data = interventions,
      aes(x = x, y = y + 0.12, label = label, colour = colour),
      size = 1.17, fontface = "bold", family = JTM_FONT, inherit.aes = FALSE
    ) +
    geom_text(
      data = interventions,
      aes(x = x, y = y - 0.14, label = effect),
      size = 1.03, family = JTM_FONT, colour = COL[["ink"]], inherit.aes = FALSE
    ) +
    geom_segment(
      data = evidence_strip,
      aes(x = accent_x, xend = accent_x, y = 0.27, yend = 0.90, colour = colour),
      linewidth = 1.10, lineend = "round", inherit.aes = FALSE
    ) +
    geom_text(
      data = evidence_strip,
      aes(x = text_x, y = 0.86, label = title, colour = colour),
      hjust = 0, size = 1.20, fontface = "bold", family = JTM_FONT,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = evidence_strip,
      aes(x = text_x, y = detail_y, label = detail),
      hjust = 0, size = 1.05, lineheight = 0.92, family = JTM_FONT,
      colour = COL[["ink"]], inherit.aes = FALSE
    ) +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0.00, 16.05), ylim = c(0.00, 8.15), clip = "off") +
    theme_void(base_size = 7, base_family = JTM_FONT) +
    theme(
      text = element_text(family = JTM_FONT, colour = COL[["ink"]]),
      plot.margin = margin(1.3, 2.0, 0.7, 2.0, "mm"),
      plot.tag = element_text(size = 8.2, face = "bold", family = JTM_FONT)
    )
}

mechanism_panel_minimal <- function(interventions) {
  flow_arrow <- grid::arrow(
    length = grid::unit(1.25, "mm"), type = "closed", angle = 24
  )
  small_arrow <- grid::arrow(
    length = grid::unit(1.00, "mm"), type = "closed", angle = 24
  )

  ggplot() +
    annotate(
      "rect", xmin = 0.15, xmax = 12.90, ymin = 4.05, ymax = 5.88,
      fill = scales::alpha(COL[["wnt_light"]], 0.12), colour = NA
    ) +
    annotate(
      "rect", xmin = 0.15, xmax = 9.72, ymin = 1.72, ymax = 3.35,
      fill = scales::alpha(COL[["neutral_light"]], 0.14), colour = NA
    ) +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(10, "mm"),
        gp = grid::gpar(
          col = scales::alpha(COL[["context"]], 0.64), fill = "white", lwd = 0.98
        )
      ),
      xmin = 8.30, xmax = 10.48, ymin = 3.62, ymax = 6.00
    ) +
    annotation_custom(
      canonical_proteasome_grob(), xmin = 8.88, xmax = 9.70, ymin = 2.12, ymax = 2.92
    ) +
    annotation_custom(
      mechanism_state_icon_grob("stem", COL[["route"]]),
      xmin = 10.72, xmax = 11.22, ymin = 4.87, ymax = 5.47
    ) +
    annotation_custom(
      mechanism_state_icon_grob("differentiated", COL[["wnt"]]),
      xmin = 10.72, xmax = 11.22, ymin = 3.87, ymax = 4.47
    ) +
    annotation_custom(
      grid::roundrectGrob(
        r = grid::unit(1.15, "mm"),
        gp = grid::gpar(
          col = scales::alpha(COL[["neutral"]], 0.68), fill = "white", lwd = 0.82
        )
      ),
      xmin = 4.05, xmax = 6.05, ymin = 1.84, ymax = 3.18
    ) +
    geom_segment(aes(x = 2.16, xend = 2.16, y = 1.55, yend = 6.02),
                 linewidth = 1.18, colour = COL[["neutral_light"]]) +
    geom_segment(aes(x = 2.25, xend = 2.25, y = 1.55, yend = 6.02),
                 linewidth = 0.56, colour = "white") +
    annotation_custom(
      canonical_receptor_grob(), xmin = 1.68, xmax = 3.02, ymin = 4.12, ymax = 5.62
    ) +
    geom_segment(aes(x = 1.10, xend = 1.78, y = 4.85, yend = 4.85),
                 linewidth = 0.68, colour = COL[["wnt"]], arrow = flow_arrow) +
    geom_segment(aes(x = 2.82, xend = 3.24, y = 4.85, yend = 4.85),
                 linewidth = 0.68, colour = COL[["wnt"]], arrow = flow_arrow) +
    geom_segment(aes(x = 3.72, xend = 6.20, y = 4.85, yend = 4.85),
                 linewidth = 0.68, colour = COL[["wnt"]], arrow = flow_arrow) +
    geom_segment(aes(x = 7.50, xend = 8.23, y = 4.85, yend = 4.85),
                 linewidth = 0.68, colour = COL[["wnt"]], arrow = flow_arrow) +
    geom_segment(aes(x = 3.48, xend = 3.48, y = 4.55, yend = 3.50),
                 linewidth = 0.50, colour = COL[["wnt"]]) +
    geom_segment(aes(x = 3.48, xend = 5.05, y = 3.50, yend = 3.50),
                 linewidth = 0.50, colour = COL[["wnt"]]) +
    geom_segment(aes(x = 5.05, xend = 5.05, y = 3.50, yend = 3.31),
                 linewidth = 0.50, colour = COL[["wnt"]]) +
    geom_segment(aes(x = 4.78, xend = 5.32, y = 3.29, yend = 3.29),
                 linewidth = 1.04, colour = COL[["wnt"]], lineend = "round") +
    geom_segment(aes(x = 4.56, xend = 4.56, y = 1.58, yend = 1.82),
                 linewidth = 0.52, colour = COL[["route"]]) +
    geom_segment(aes(x = 4.31, xend = 4.81, y = 1.82, yend = 1.82),
                 linewidth = 1.02, colour = COL[["route"]], lineend = "round") +
    geom_segment(aes(x = 6.06, xend = 6.55, y = 2.53, yend = 2.53),
                 linewidth = 0.55, colour = COL[["neutral"]], arrow = small_arrow) +
    geom_segment(aes(x = 7.18, xend = 7.54, y = 2.53, yend = 2.53),
                 linewidth = 0.55, colour = COL[["neutral"]], arrow = small_arrow) +
    geom_segment(aes(x = 8.25, xend = 8.82, y = 2.53, yend = 2.53),
                 linewidth = 0.55, colour = COL[["neutral"]], arrow = small_arrow) +
    geom_segment(aes(x = 10.48, xend = 10.68, y = 5.17, yend = 5.17),
                 linewidth = 0.58, colour = COL[["route"]], arrow = small_arrow) +
    geom_segment(aes(x = 10.48, xend = 10.68, y = 4.17, yend = 4.17),
                 linewidth = 0.58, colour = COL[["wnt"]], arrow = small_arrow) +
    geom_segment(aes(x = 12.72, xend = 12.98, y = 5.17, yend = 5.17),
                 linewidth = 0.46, colour = COL[["route"]]) +
    geom_segment(aes(x = 12.72, xend = 12.98, y = 4.17, yend = 4.17),
                 linewidth = 0.46, colour = COL[["wnt"]]) +
    geom_segment(aes(x = 12.98, xend = 12.98, y = 4.17, yend = 5.17),
                 linewidth = 0.46, colour = COL[["neutral"]]) +
    geom_segment(aes(x = 12.98, xend = 13.24, y = 4.67, yend = 4.67),
                 linewidth = 0.62, colour = COL[["ink"]], arrow = flow_arrow) +
    geom_point(aes(x = 0.72, y = 4.85), shape = 21, size = 4.8, stroke = 0.70,
               fill = scales::alpha(COL[["wnt"]], 0.17), colour = COL[["wnt"]]) +
    geom_point(aes(x = 3.48, y = 4.85), shape = 23, size = 4.05, stroke = 0.74,
               fill = "white", colour = COL[["wnt"]]) +
    geom_point(
      data = data.frame(
        x = c(6.48, 6.86, 7.22), y = c(4.77, 5.00, 4.77)
      ),
      aes(x, y), inherit.aes = FALSE, shape = 21, size = 2.78, stroke = 0.58,
      fill = scales::alpha(COL[["wnt"]], 0.16), colour = COL[["wnt"]]
    ) +
    geom_point(aes(x = 6.86, y = 2.53), shape = 21, size = 3.45, stroke = 0.60,
               fill = "white", colour = COL[["neutral"]]) +
    geom_point(
      data = data.frame(x = c(7.70, 7.93, 8.16), y = c(2.41, 2.62, 2.41)),
      aes(x, y), inherit.aes = FALSE, shape = 21, size = 2.10, stroke = 0.48,
      fill = scales::alpha(COL[["neutral_light"]], 0.55), colour = COL[["neutral"]]
    ) +
    annotate("text", x = 0.25, y = 5.60, label = "WNT PRESENT",
             hjust = 0, size = 1.30, fontface = "bold", family = JTM_FONT,
             colour = COL[["wnt"]]) +
    annotate("text", x = 0.25, y = 3.08, label = "WNT ABSENT\n(APC intact)",
             hjust = 0, size = 1.18, lineheight = 0.92, fontface = "bold", family = JTM_FONT,
             colour = COL[["neutral"]]) +
    annotate("text", x = 0.72, y = 4.85, label = "WNT",
             size = 1.42, fontface = "bold", family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 2.34, y = 4.01, label = "FZD–LRP5/6 complex",
             size = 1.20, fontface = "bold", family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 3.48, y = 4.85, label = "DVL",
             size = 1.34, fontface = "bold", family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 4.57, y = 5.19, label = "p-LRP5/6 recruits AXIN",
             size = 0.96, family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 4.20, y = 3.71, label = "inhibits",
             size = 1.00, family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 6.86, y = 5.46, label = "stabilised β-catenin",
             size = 1.65, fontface = "bold", family = JTM_FONT, colour = COL[["wnt"]]) +
    annotate("text", x = 7.74, y = 5.09, label = "nuclear entry",
             size = 1.15, family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 5.05, y = 2.97, label = "β-catenin\ndestruction complex",
             size = 1.20, lineheight = 0.88, fontface = "bold",
             family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("label", x = 4.56, y = 2.46, label = "APC", size = 1.00,
             fontface = "bold", family = JTM_FONT, colour = COL[["route"]], fill = "#FFF8F4",
             linewidth = 0.17, label.padding = grid::unit(0.60, "mm")) +
    annotate("text", x = 4.56, y = 1.40, label = "APC loss",
             size = 1.00, fontface = "bold", family = JTM_FONT,
             colour = COL[["route"]]) +
    annotate("label", x = 5.47, y = 2.46, label = "AXIN", size = 1.00,
             fontface = "bold", family = JTM_FONT, colour = COL[["ink"]], fill = "white",
             linewidth = 0.17, label.padding = grid::unit(0.60, "mm")) +
    annotate("label", x = 4.56, y = 2.06, label = "GSK3β", size = 0.96,
             fontface = "bold", family = JTM_FONT, colour = COL[["ink"]], fill = "white",
             linewidth = 0.17, label.padding = grid::unit(0.54, "mm")) +
    annotate("label", x = 5.47, y = 2.06, label = "CK1", size = 1.00,
             fontface = "bold", family = JTM_FONT, colour = COL[["ink"]], fill = "white",
             linewidth = 0.17, label.padding = grid::unit(0.60, "mm")) +
    annotate("text", x = 6.86, y = 2.53, label = "P", size = 1.00,
             fontface = "bold", family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 6.86, y = 1.52, label = "p–β-catenin",
             size = 1.02, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 7.93, y = 2.41, label = "Ub", size = 0.88,
             fontface = "bold", family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 7.93, y = 1.52, label = "β-TrCP / ubiquitin",
             size = 0.98, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 9.29, y = 1.52, label = "proteasomal degradation",
             size = 0.98, family = JTM_FONT, colour = COL[["neutral"]]) +
    annotate("text", x = 6.48, y = 2.96, label = "phosphorylation",
             size = 0.93, family = JTM_FONT,
             colour = COL[["neutral"]]) +
    annotate("text", x = 9.39, y = 5.70, label = "NUCLEUS",
             size = 1.12, fontface = "bold", family = JTM_FONT,
             colour = COL[["neutral"]]) +
    annotate("text", x = 9.39, y = 4.98, label = "β-catenin + TCF7L2 / LEF1",
             size = 1.38, fontface = "bold", family = JTM_FONT,
             colour = COL[["context"]]) +
    annotate("segment", x = 8.77, xend = 10.02, y = 4.69, yend = 4.69,
             linewidth = 0.27, colour = COL[["neutral_light"]]) +
    annotate("text", x = 9.39, y = 4.34, label = "ASCL2 target programme",
             size = 1.28, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 11.36, y = 5.17, label = "Stem / progenitor ↑",
             hjust = 0, size = 1.43, fontface = "bold", family = JTM_FONT,
             colour = COL[["route"]]) +
    annotate("text", x = 11.36, y = 4.17, label = "Differentiation ↓",
             hjust = 0, size = 1.43, fontface = "bold", family = JTM_FONT,
             colour = COL[["wnt"]]) +
    annotate("text", x = 14.26, y = 5.52, label = "FIXED EPITHELIAL PROGRAMME",
             size = 1.62, fontface = "bold", family = JTM_FONT,
             colour = COL[["route"]]) +
    annotate("segment", x = 13.45, xend = 15.07, y = 5.22, yend = 5.22,
             linewidth = 0.30, colour = COL[["route_light"]]) +
    annotate("segment", x = 13.45, xend = 13.91, y = 4.68, yend = 4.68,
             linewidth = 1.40, colour = COL[["route"]], lineend = "round") +
    annotate("text", x = 14.08, y = 4.68, label = "50-gene up arm", hjust = 0,
             size = 1.36, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("segment", x = 13.45, xend = 13.91, y = 4.08, yend = 4.08,
             linewidth = 1.40, colour = COL[["wnt"]], lineend = "round") +
    annotate("text", x = 14.08, y = 4.08, label = "50-gene down arm", hjust = 0,
             size = 1.36, family = JTM_FONT, colour = COL[["ink"]]) +
    annotate("text", x = 0.25, y = 6.31, label = "EXTRACELLULAR",
             hjust = 0, size = 1.32, fontface = "bold", family = JTM_FONT,
             colour = COL[["neutral"]]) +
    annotate("text", x = 2.82, y = 6.31, label = "CYTOPLASM",
             hjust = 0, size = 1.32, fontface = "bold", family = JTM_FONT,
             colour = COL[["neutral"]]) +
    annotate("text", x = 8.37, y = 6.31, label = "TRANSCRIPTION",
             hjust = 0, size = 1.32, fontface = "bold", family = JTM_FONT,
             colour = COL[["neutral"]]) +
    annotate("text", x = 2.47, y = 2.74, label = "cell membrane",
             angle = 90, size = 1.02, family = JTM_FONT,
             colour = COL[["neutral"]]) +
    annotate("segment", x = 0.18, xend = 12.85, y = 1.14, yend = 1.14,
             linewidth = 0.26, colour = COL[["neutral_light"]]) +
    geom_segment(
      data = interventions,
      aes(x = x - 0.23, xend = x + 0.23, y = 0.86, yend = 0.86, colour = colour),
      linewidth = 1.05, lineend = "round", inherit.aes = FALSE
    ) +
    geom_text(
      data = interventions,
      aes(x = x, y = 0.46, label = label, colour = colour),
      size = 1.22, fontface = "bold", family = JTM_FONT, inherit.aes = FALSE
    ) +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0.00, 15.55), ylim = c(0.00, 6.48), clip = "off") +
    theme_void(base_size = 7, base_family = JTM_FONT) +
    theme(
      text = element_text(family = JTM_FONT, colour = COL[["ink"]]),
      plot.margin = margin(1.2, 1.8, 0.8, 1.8, "mm"),
      plot.tag = element_text(size = 8.2, face = "bold", family = JTM_FONT)
    )
}

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
    figure = stem, format = names(paths), file = normalizePath(unname(paths)),
    width_mm = width_mm, height_mm = height_mm,
    resolution_dpi = c(NA, NA, 600, 300), stringsAsFactors = FALSE
  )
}

# Reuse verified panel-building code in isolated environments. Its legacy
# exports are written only to the R session's temporary directory.
old_figure_dir <- Sys.getenv("JTM_FIGURE_DIR", unset = NA_character_)
old_figure_stem <- Sys.getenv("JTM_FIGURE_STEM", unset = NA_character_)

env_main <- new.env(parent = globalenv())
Sys.setenv(JTM_FIGURE_DIR = file.path(tempdir(), "jtm_v16_mainline_components"))
sys.source(file.path(ROOT, "analysis", "plot_jtm_mainline_figures_v0_5.R"), envir = env_main)

env_gap <- new.env(parent = globalenv())
Sys.setenv(JTM_FIGURE_DIR = file.path(tempdir(), "jtm_v16_gap_components"))
sys.source(file.path(ROOT, "analysis", "plot_jtm_v06_new_supplementary_figures.R"), envir = env_gap)

env_closure <- new.env(parent = globalenv())
Sys.setenv(
  JTM_FIGURE_DIR = file.path(tempdir(), "jtm_v16_closure_components"),
  JTM_FIGURE_STEM = "closure_component"
)
sys.source(file.path(ROOT, "analysis", "plot_computational_closure_validation.R"), envir = env_closure)

if (is.na(old_figure_dir)) Sys.unsetenv("JTM_FIGURE_DIR") else Sys.setenv(JTM_FIGURE_DIR = old_figure_dir)
if (is.na(old_figure_stem)) Sys.unsetenv("JTM_FIGURE_STEM") else Sys.setenv(JTM_FIGURE_STEM = old_figure_stem)

# -----------------------------------------------------------------------------
# Figure 1: study workflow, biological identity and held-out validation (5 panels)
# -----------------------------------------------------------------------------

study_branches <- data.frame(
  y = c(6.05, 4.78, 3.50, 2.22, 0.95),
  title = c(
    "Independent replication", "Regulatory concordance", "Recurrence in CRC states",
    "Bidirectional perturbation", "Tissue-level readouts"
  ),
  detail = c(
    "Held-out · 5 cohorts · 51 FFPE pairs · 7 analyses",
    "Becker single-nucleus · matched RNA–ATAC",
    "CRC Atlas · 33 source-study omissions",
    "APC · WNT · ASCL2 · TCF7L2 · virtual KO",
    "6 Visium sections · 4 proteomic datasets"
  ),
  figure = c("Figs. 1–2", "Fig. 3", "Fig. 4", "Fig. 5", "Fig. 6"),
  icon = c("replicate", "wnt", "transport", "virtual", "ffpe"),
  accent = unname(c(
    COL[["route"]], COL[["wnt"]], COL[["crc"]], COL[["context"]], COL[["adenoma"]]
  )),
  fill = c("#FCF6F3", "#F5F8FA", "#F3F8F7", "#F7F5FA", "#FCF9F1"),
  stringsAsFactors = FALSE
)
study_overview_nodes <- bind_rows(
  data.frame(
    node = c("Adenoma discovery", "Fixed epithelial state", "Biological conclusion"),
    role = c(
      "Donor-level discovery and bootstrap selection",
      "Fixed 50-up/50-down programme",
      "WNT/stem-progenitor-high, differentiation-low state"
    ),
    figure = c("Figure 1", "Figures 1–6", "Figures 1–6"),
    stringsAsFactors = FALSE
  ),
  study_branches %>% transmute(node = title, role = detail, figure = figure)
)
p1a <- study_overview_panel(study_branches)

signature <- read_tsv("results/figure_data_locked/fig1_locked_signature_genes.tsv")
biological_anchor_genes <- c(
  "OLFM4", "ASCL2", "LGR5", "AXIN2", "RNF43", "ZNRF3", "NKD1", "EPHB3",
  "FABP1", "PHGR1", "PCK1", "LGALS4", "CA2", "AQP8", "HMGCS2", "GUCA2A"
)
signature_anchors <- signature %>%
  filter(gene %in% biological_anchor_genes) %>%
  mutate(
    direction = factor(
      ifelse(
        signature_direction == "adenoma_up",
        "WNT / stem–progenitor ↑",
        "Mature differentiation ↓"
      ),
      levels = c("WNT / stem–progenitor ↑", "Mature differentiation ↓")
    ),
    gene = reorder(gene, discovery_effect_adenoma_minus_normal)
  )
stopifnot(nrow(signature_anchors) == length(biological_anchor_genes))

p1c <- ggplot(signature_anchors,
              aes(discovery_effect_adenoma_minus_normal, gene, fill = direction)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = COL[["neutral"]]) +
  scale_fill_manual(values = c(
    "WNT / stem–progenitor ↑" = COL[["route"]],
    "Mature differentiation ↓" = COL[["wnt"]]
  )) +
  labs(x = "Discovery donor-median effect", y = NULL, fill = NULL) +
  theme_jtm() +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.text = element_text(size = 4.9),
    legend.key.width = grid::unit(2.4, "mm"),
    axis.text.y = element_text(size = 5.8)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

heldout_data <- env_main$chen_plot %>% filter(dataset == "Held-out validation")
heldout_ann <- env_main$chen_ann %>% filter(dataset == "Held-out validation")
heldout_panel <- ggplot(
  heldout_data, aes(tissue, score__ca_route_signature, colour = tissue)
) +
  geom_boxplot(
    width = 0.52, outlier.shape = NA, linewidth = 0.40,
    colour = COL[["ink"]], fill = "#FAFBFB"
  ) +
  geom_point(
    position = position_jitter(width = 0.10, seed = 20260710),
    alpha = 0.74, size = 0.92
  ) +
  geom_text(
    data = heldout_ann, aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE, size = 1.82, lineheight = 0.92, family = JTM_FONT
  ) +
  scale_colour_manual(values = c(
    "Normal" = COL[["neutral"]], "Conventional adenoma" = COL[["route"]]
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "Epithelial programme score") +
  theme_jtm() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 14, hjust = 1, size = 5.7)
  )

p1b <- clean_panel(env_main$p1b) +
  scale_colour_manual(
    values = c(
      "Not locked" = COL[["neutral_light"]],
      "Adenoma-down" = COL[["wnt"]],
      "Adenoma-up" = COL[["route"]]
    ),
    breaks = c("Adenoma-down", "Adenoma-up"), name = NULL
  ) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

p1e <- ggplot(env_main$chen_pairs, aes(dataset, delta, colour = dataset)) +
  geom_hline(
    yintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_point(
    position = position_jitter(width = 0.085, seed = 20260710),
    size = 1.08, alpha = 0.82
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 23, size = 2.5, stroke = 0.48,
    fill = "white", colour = COL[["ink"]]
  ) +
  geom_text(
    data = env_main$pair_ann,
    aes(x = dataset, y = y, label = label, hjust = hjust),
    inherit.aes = FALSE, size = 1.82, family = JTM_FONT
  ) +
  scale_colour_manual(values = c(Discovery = COL[["adenoma"]], `Held-out` = COL[["route"]])) +
  labs(x = NULL, y = "Paired score change\n(adenoma − normal)") +
  theme_jtm() +
  theme(legend.position = "none")

fig1_mid <- (clean_panel(p1b) | clean_panel(p1c)) +
  plot_layout(widths = c(1.12, 0.88))
fig1_bottom <- (clean_panel(heldout_panel) | clean_panel(p1e)) +
  plot_layout(widths = c(1.0, 1.0))
fig1 <- clean_panel(p1a) / fig1_mid / fig1_bottom +
  plot_layout(heights = c(1.38, 1.06, 0.96)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(study_overview_nodes, "figure1a_study_evidence_architecture.tsv")
write_tsv(signature_anchors, "figure1c_biological_anchor_genes.tsv")

# -----------------------------------------------------------------------------
# Figure 2: multicohort, FFPE and expanded transcriptomic replication (6 panels)
# -----------------------------------------------------------------------------

p2b <- ggplot(
  env_main$external_primary,
  aes(clustered_standardized_mean_difference, cohort_label, colour = evidence)
) +
  geom_vline(
    xintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_segment(
    aes(x = clustered_standardized_ci_low, xend = clustered_standardized_ci_high,
        yend = cohort_label),
    linewidth = 0.68, lineend = "round"
  ) +
  geom_point(size = 1.95) +
  scale_colour_manual(
    values = c(Clear = COL[["route"]], Imprecise = COL[["uncertain"]]), name = NULL
  ) +
  coord_cartesian(xlim = c(-0.80, 2.65)) +
  labs(x = "Standardised adenoma effect", y = NULL) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

p2c <- ggplot(env_main$robust, aes(estimate, model)) +
  geom_vline(
    xintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_segment(
    aes(x = loo_low, xend = loo_high, yend = model),
    linewidth = 2.4, colour = "#E1E5E8", lineend = "round"
  ) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = model),
    linewidth = 0.68, colour = COL[["ink"]], lineend = "round"
  ) +
  geom_point(size = 1.95, colour = COL[["route"]]) +
  geom_text(
    aes(x = ci_high + 0.08, label = sprintf("%.2f SD", estimate)),
    hjust = 0, size = 1.74, family = JTM_FONT, colour = COL[["ink"]]
  ) +
  coord_cartesian(xlim = c(0, 2.36), clip = "off") +
  labs(x = "Standardised adenoma effect", y = NULL) +
  theme_jtm()

specificity_effects <- bind_rows(
  env_main$external_grade %>%
    filter(comparison == "high_grade_vs_low_grade") %>%
    transmute(
      comparison = "High-grade vs low-grade",
      estimate = clustered_standardized_mean_difference,
      ci_low = clustered_standardized_ci_low,
      ci_high = clustered_standardized_ci_high,
      p_value = p_patient_clustered_standardized_ols,
      q_value = q_value_bh,
      evidence = "Supportive"
    ),
  env_main$external_tests %>%
    filter(
      cohort == "GSE40362", signature_size_per_direction == 50,
      comparison %in% c("hyperplastic_vs_normal", "adenoma_vs_hyperplastic")
    ) %>%
    transmute(
      comparison = ifelse(
        comparison == "hyperplastic_vs_normal",
        "Hyperplastic vs normal", "Adenoma vs hyperplastic"
      ),
      estimate = clustered_standardized_mean_difference,
      ci_low = clustered_standardized_ci_low,
      ci_high = clustered_standardized_ci_high,
      p_value = p_patient_clustered_standardized_ols,
      q_value = q_value_within_comparison_and_size,
      evidence = ifelse(comparison == "Hyperplastic vs normal", "Near null", "Supportive")
    )
) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = rev(c(
        "High-grade vs low-grade", "Adenoma vs hyperplastic", "Hyperplastic vs normal"
      ))
    ),
    label = sprintf("%.2f SD", estimate)
  )

p2_specificity <- ggplot(specificity_effects, aes(estimate, comparison, colour = evidence)) +
  geom_vline(
    xintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = comparison),
    linewidth = 0.68, lineend = "round"
  ) +
  geom_point(size = 1.95) +
  geom_text(
    aes(x = ci_high + 0.08, label = label), hjust = 0,
    size = 1.65, family = JTM_FONT, colour = COL[["ink"]]
  ) +
  scale_colour_manual(values = c("Supportive" = COL[["route"]], "Near null" = COL[["neutral"]])) +
  coord_cartesian(xlim = c(-1.18, 2.72), clip = "off") +
  labs(x = "Standardised programme difference", y = NULL) +
  theme_jtm(base_size = 6.7) +
  theme(
    legend.position = "none", axis.text.y = element_text(size = 5.45),
    plot.margin = margin(1.6, 3.0, 1.6, 1.8, "mm")
  )

p2d <- ggplot(env_gap$paired_long, aes(tissue_group, route_score, group = patient_id)) +
  geom_line(linewidth = 0.30, alpha = 0.22, colour = COL[["neutral"]]) +
  geom_point(aes(colour = tissue_group), size = 0.78, alpha = 0.76) +
  geom_point(
    data = env_gap$paired_medians, aes(tissue_group, route_score), inherit.aes = FALSE,
    shape = 23, size = 2.3, stroke = 0.50, fill = "white", colour = COL[["ink"]]
  ) +
  scale_colour_manual(values = c("Reference" = COL[["wnt"]], "Adenoma" = COL[["route"]])) +
  scale_x_discrete(labels = c("Reference" = "Ref.", "Adenoma" = "Aden.")) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = NULL, y = "Epithelial programme score") + theme_jtm() +
  theme(legend.position = "none")

p2e <- clean_panel(env_gap$pS5b) +
  scale_colour_manual(
    values = c(
      "Patient-pair bootstrap 95% CI" = COL[["route"]],
      "HC3 95% CI" = COL[["wnt"]]
    ),
    labels = c(
      "Patient-pair bootstrap 95% CI" = "Pair bootstrap",
      "HC3 95% CI" = "HC3"
    )
  ) +
  labs(x = "Adenoma − reference programme change") +
  theme_jtm(base_size = 6.8) +
  theme(legend.position = "top", legend.justification = "left")

expanded_v11 <- env_gap$expanded_plot %>%
  mutate(
    route_short = recode(
      as.character(route_label),
      "GSE164541 (paired)" = "GSE164541 · paired",
      "Microdissected crypt" = "Crypt · Warsaw",
      "Microdissected mucosa" = "Mucosa · Warsaw",
      "Macrodissected" = "Macro · Warsaw"
    ),
    route_display_text = case_when(
      as.character(cluster_label) == "GSE164541" ~ route_short,
      as.character(cluster_label) == "Warsaw" ~ route_short,
      TRUE ~ paste(route_short, as.character(cluster_label), sep = " · ")
    )
  )
expanded_v11_levels <- expanded_v11 %>%
  arrange(route_label) %>%
  pull(route_display_text)
expanded_v11 <- expanded_v11 %>%
  mutate(route_display = factor(route_display_text, levels = expanded_v11_levels))
p2f <- ggplot(expanded_v11, aes(clustered_mean_difference, route_display)) +
  geom_hline(
    yintercept = seq_along(levels(expanded_v11$route_display)),
    linewidth = 0.25, colour = COL[["neutral_pale"]]
  ) +
  geom_vline(
    xintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_segment(
    aes(x = clustered_ci_low, xend = clustered_ci_high, yend = route_display,
        colour = near_null),
    linewidth = 0.62, lineend = "round"
  ) +
  geom_point(aes(colour = near_null), size = 1.75) +
  scale_colour_manual(values = c(`FALSE` = COL[["route"]], `TRUE` = COL[["neutral"]])) +
  coord_cartesian(xlim = c(-0.55, 3.45), clip = "off") +
  labs(x = "Adenoma − normal programme effect", y = NULL) +
  theme_jtm(base_size = 6.6) +
  theme(
    legend.position = "none", axis.text.y = element_text(size = 4.9),
    plot.margin = margin(1.4, 1.6, 1.4, 1.6, "mm")
  )

fig2_row1 <- (clean_panel(p2b) | clean_panel(p2c)) +
  plot_layout(widths = c(1.08, 0.92))
fig2_row2 <- (clean_panel(p2_specificity) | clean_panel(p2d)) +
  plot_layout(widths = c(0.92, 1.08))
fig2_row3 <- (clean_panel(p2e) | clean_panel(p2f)) +
  plot_layout(widths = c(0.92, 1.08))
fig2 <- fig2_row1 / fig2_row2 / fig2_row3 +
  plot_layout(heights = c(0.94, 1.05, 0.98)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(specificity_effects, "figure2c_histological_specificity.tsv")
write_tsv(expanded_v11, "figure2f_expanded_route_display.tsv")

# -----------------------------------------------------------------------------
# Figure 3: Becker transfer and RNA-ATAC support (6 panels)
# -----------------------------------------------------------------------------

b3a <- clean_panel(env_main$p2a) +
  scale_colour_manual(values = c(
    Normal = COL[["neutral"]], Polyp = COL[["adenoma"]], CRC = COL[["crc"]]
  )) +
  labs(y = "Epithelial programme score") +
  theme_jtm() +
  theme(legend.position = "none")

b3b <- ggplot(env_main$becker_forest, aes(coef, comparison)) +
  geom_vline(
    xintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = comparison),
    linewidth = 0.68, lineend = "round", colour = COL[["ink"]]
  ) +
  geom_point(size = 1.95, colour = COL[["route"]]) +
  geom_text(
    aes(x = 2.48, label = p_label), hjust = 0, size = 1.72,
    family = JTM_FONT, colour = COL[["ink"]]
  ) +
  coord_cartesian(xlim = c(-0.20, 3.25), clip = "off") +
  labs(x = "Adjusted coefficient (95% CI)", y = NULL) +
  theme_jtm()

b3c <- clean_panel(env_main$p2c) +
  scale_colour_manual(values = c(Normal = COL[["neutral"]], Polyp = COL[["adenoma"]])) +
  labs(y = "Epithelial RNA programme") +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

b3d <- ggplot(env_main$rna_forest, aes(coef, label)) +
  geom_vline(
    xintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = label),
    linewidth = 0.66, lineend = "round", colour = COL[["ink"]]
  ) +
  geom_point(size = 1.9, colour = COL[["route"]]) +
  coord_cartesian(xlim = c(-0.05, 0.90)) +
  labs(x = "RNA-from-ATAC coefficient (95% CI)", y = NULL) +
  theme_jtm()

b3e <- clean_panel(env_main$pS3b) +
  scale_colour_manual(values = c(
    "Patient clustered" = COL[["route"]],
    "Patient fixed effect" = COL[["wnt"]]
  )) +
  theme_jtm(base_size = 6.9) +
  theme(legend.position = "top", legend.justification = "left")

b3f <- clean_panel(env_main$pS3c) +
  scale_colour_manual(values = c(
    "Locked route" = COL[["route"]],
    "WNT/stemness" = COL[["wnt"]],
    "Proliferation" = COL[["neutral"]]
  ), labels = c(
    "Locked route" = "Programme", "WNT/stemness" = "WNT/stem.",
    "Proliferation" = "Prolif."
  )) +
  facet_wrap(
    ~locus, nrow = 1,
    labeller = as_labeller(c(
      "WNT-route loci" = "WNT-programme loci",
      "TCF/ASCL2-axis loci" = "TCF/ASCL2-axis loci"
    ))
  ) +
  theme_jtm(base_size = 6.9) +
  theme(
    legend.position = "top", legend.justification = "left",
    axis.text.x = element_text(angle = 24, hjust = 1, size = 5.6)
  )

fig3_row1 <- (clean_panel(b3a) | clean_panel(b3b)) +
  plot_layout(widths = c(1.12, 0.88))
fig3_row2 <- (clean_panel(b3c) | clean_panel(b3d)) +
  plot_layout(widths = c(1.12, 0.88))
fig3_row3 <- (clean_panel(b3e) | clean_panel(b3f)) +
  plot_layout(widths = c(1.12, 0.88))
fig3 <- fig3_row1 / fig3_row2 / fig3_row3 +
  plot_layout(heights = c(0.82, 1.08, 0.94)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Figure 4: CRC Atlas state recapitulation and source audit (5 panels)
# -----------------------------------------------------------------------------

atlas_state_colours <- c(
  "Normal\nepithelium" = COL[["neutral"]],
  "Polyp\nepithelium" = COL[["adenoma"]],
  "Polyp\ncancer" = "#E6BF70",
  "Primary\nepithelium" = "#78A9C8",
  "Primary\ncancer" = COL[["crc"]],
  "Metastasis\nepithelium" = "#A393BC",
  "Metastasis\ncancer" = COL[["context"]]
)

p4a <- ggplot(
  env_main$atlas_plot,
  aes(state, score__ca_route_signature, colour = state)
) +
  geom_boxplot(
    width = 0.50, outlier.shape = NA, linewidth = 0.38,
    colour = COL[["ink"]], fill = "#FAFBFB"
  ) +
  geom_point(
    position = position_jitter(width = 0.12, seed = 20260710),
    alpha = 0.31, size = 0.52
  ) +
  geom_text(
    data = env_main$atlas_n,
    aes(x = state, y = -Inf, label = paste0("n = ", n)),
    inherit.aes = FALSE, vjust = -0.32, size = 1.65, family = JTM_FONT
  ) +
  scale_colour_manual(values = atlas_state_colours) +
  labs(x = NULL, y = "Epithelial programme score") +
  theme_jtm() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1, size = 5.7, lineheight = 0.88)
  )

p4b <- ggplot(env_main$carrier_models, aes(coef, state, colour = programme)) +
  geom_vline(
    xintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = state),
    linewidth = 0.60, lineend = "round", position = position_dodge(width = 0.44)
  ) +
  geom_point(size = 1.75, position = position_dodge(width = 0.44)) +
  scale_colour_manual(values = c(
    "Locked route" = COL[["route"]], "WNT/stemness" = COL[["wnt"]]
  ), labels = c("Locked route" = "Programme", "WNT/stemness" = "WNT/stemness")) +
  labs(x = "Coefficient (95% CI)", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.9) +
  theme(legend.position = "top", legend.justification = "left")

p4c <- ggplot(env_main$influence_plot, aes(full_coef, state_label, colour = programme)) +
  geom_vline(
    xintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_segment(
    aes(x = loo_min_coef, xend = loo_max_coef, yend = state_label, group = programme),
    linewidth = 2.15, alpha = 0.24, position = position_dodge(width = 0.46),
    lineend = "round"
  ) +
  geom_segment(
    aes(x = full_ci_low, xend = full_ci_high, yend = state_label, group = programme),
    linewidth = 0.56, position = position_dodge(width = 0.46), lineend = "round"
  ) +
  geom_point(size = 1.65, position = position_dodge(width = 0.46)) +
  scale_colour_manual(values = c(
    "Locked route" = COL[["route"]], "WNT/stemness" = COL[["wnt"]]
  ), labels = c("Locked route" = "Programme", "WNT/stemness" = "WNT/stemness")) +
  labs(x = "Adjusted coefficient", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.9) +
  theme(legend.position = "none")

compact_study_labels <- function(study_ids) {
  base <- sub("^([^_]+)_([0-9]{4}).*$", "\\1 \\2", study_ids)
  base[study_ids == "MUI_Innsbruck"] <- "MUI Innsbruck"
  duplicate_base <- duplicated(base) | duplicated(base, fromLast = TRUE)
  suffix <- case_when(
    grepl("Cancer_Lett", study_ids) ~ "Cancer Lett",
    grepl("PLoS_Genet", study_ids) ~ "PLoS Genet",
    TRUE ~ ""
  )
  ifelse(duplicate_base, paste0(base, " · ", suffix), base)
}

study_order <- sort(unique(env_main$support_grid$study_id))
study_label_map <- setNames(compact_study_labels(study_order), study_order)
state_compact <- c(
  "Normal epith.", "Polyp epith.", "Polyp cancer", "Primary epith.",
  "Primary cancer", "Metastatic epith.", "Metastatic cancer"
)
support_compact <- env_main$support_grid %>%
  mutate(
    study_short = factor(
      study_label_map[study_id], levels = rev(unname(study_label_map[study_order]))
    ),
    state_short = factor(state, levels = levels(state), labels = state_compact)
  )
p4d <- ggplot(support_compact, aes(state_short, study_short, fill = fill_value)) +
  geom_tile(colour = "white", linewidth = 0.30) +
  geom_text(
    aes(label = donor_label), size = 1.52, colour = COL[["ink"]], family = JTM_FONT
  ) +
  scale_fill_gradient(
    low = "#E7EFF4", high = COL[["wnt"]],
    breaks = log10(c(2, 6, 21, 81)), labels = c("1", "5", "20", "80"),
    na.value = COL[["neutral_pale"]], name = "Donors"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.6, base_family = JTM_FONT) +
  theme(
    panel.grid = element_blank(), axis.title = element_blank(),
    axis.text.x = element_text(angle = 43, hjust = 1, vjust = 1, size = 5.0),
    axis.text.y = element_text(size = 5.0),
    legend.position = "top", legend.title = element_text(size = 5.6, face = "bold"),
    legend.text = element_text(size = 5.2),
    legend.key.width = grid::unit(13, "mm"), legend.key.height = grid::unit(2.0, "mm"),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm")
  ) +
  guides(fill = guide_colourbar(title.position = "left", title.hjust = 0.5))

within_compact <- env_main$within_primary %>%
  mutate(
    study_short = factor(
      study_label_map[as.character(sub(" ", "_", study))],
      levels = rev(unique(study_label_map[as.character(sub(" ", "_", study))]))
    )
  )
if (any(is.na(within_compact$study_short))) {
  within_compact <- env_main$within_primary %>%
    mutate(
      study_id_recovered = gsub(" ", "_", as.character(study)),
      study_short = factor(
        compact_study_labels(study_id_recovered),
        levels = rev(unique(compact_study_labels(study_id_recovered)))
      )
    )
}
p4e <- ggplot(within_compact, aes(coef, study_short)) +
  geom_vline(
    xintercept = 0, linewidth = 0.32, linetype = "22", colour = COL[["neutral"]]
  ) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = study_short),
    linewidth = 0.58, lineend = "round", colour = COL[["ink"]]
  ) +
  geom_point(size = 1.72, colour = COL[["route"]]) +
  geom_text(
    aes(x = max(ci_high) + 0.13, label = label), hjust = 0,
    size = 1.55, family = JTM_FONT, colour = COL[["neutral"]]
  ) +
  coord_cartesian(
    xlim = c(min(within_compact$ci_low) - 0.05, max(within_compact$ci_high) + 0.70),
    clip = "off"
  ) +
  labs(x = "Primary cancer − normal coefficient", y = NULL) +
  theme_jtm(base_size = 6.7) +
  theme(axis.text.y = element_text(size = 5.2))

fig4_row2 <- (clean_panel(p4b) | clean_panel(p4c)) +
  plot_layout(widths = c(1, 1))
fig4_row3 <- (clean_panel(p4d) | clean_panel(p4e)) +
  plot_layout(widths = c(1.32, 0.68))
fig4 <- clean_panel(p4a) / fig4_row2 / fig4_row3 +
  plot_layout(heights = c(0.75, 0.91, 1.34)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(support_compact, "figure4d_compact_source_support.tsv")
write_tsv(within_compact, "figure4e_compact_within_study_effects.tsv")

# -----------------------------------------------------------------------------
# Figure 5: perturbation validation and virtual deletion (6 panels)
# -----------------------------------------------------------------------------

mechanism_interventions <- data.frame(
  id = c("wnt_withdrawal", "apc_perturbation", "tf_perturbation"),
  x = c(0.95, 5.02, 9.42),
  label = c(
    "WNT withdrawal · programme ↓", "APC loss ↑ · restoration ↓",
    "TF knockout / virtual TCF7L2 deletion · programme ↓"
  ),
  effect = c(
    "programme ↓", "loss ↑ · restoration ↓", "programme ↓"
  ),
  colour = unname(c(
    COL[["wnt"]], COL[["route"]], COL[["context"]]
  )),
  stringsAsFactors = FALSE
)
mechanism_nodes <- data.frame(
  node = c(
    "WNT ligand", "FZD–LRP5/6 receptor complex", "DVL",
    "APC–AXIN–GSK3β–CK1 destruction complex", "β-catenin phosphorylation",
    "β-TrCP ubiquitination", "26S proteasomal degradation",
    "β-catenin stabilisation and nuclear entry",
    "β-catenin–TCF7L2/LEF1–ASCL2 transcriptional state",
    "Stem/progenitor arm", "Differentiation arm", "Fixed epithelial programme"
  ),
  inferred_effect = c(
    "binds receptor", "activates DVL", "inhibits the destruction complex",
    "phosphorylates β-catenin when active", "permits β-TrCP recognition",
    "targets β-catenin for degradation", "reduces cytoplasmic β-catenin",
    "increases when WNT is on or APC is lost", "activates",
    "increases", "decreases", "integrates both arms"
  ),
  stringsAsFactors = FALSE
)
p5a <- mechanism_panel_minimal(mechanism_interventions)

apc_scores <- read_tsv(
  "results/perturbation_validation_locked_route/gse125472_sample_scores.tsv"
) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "APC")),
    condition = case_when(
      genotype == "WT" & wnt_rspo == "with" ~ "WT\n+WNT/RSPO",
      genotype == "APC" & wnt_rspo == "with" ~ "APC-KO\n+WNT/RSPO",
      genotype == "WT" & wnt_rspo == "without" ~ "WT\n−WNT/RSPO",
      genotype == "APC" & wnt_rspo == "without" ~ "APC-KO\n−WNT/RSPO"
    ),
    condition = factor(condition, levels = c(
      "WT\n+WNT/RSPO", "APC-KO\n+WNT/RSPO", "WT\n−WNT/RSPO", "APC-KO\n−WNT/RSPO"
    )),
    x = c(1.0, 2.0, 3.4, 4.4)[match(condition, levels(condition))]
  )
apc_segments <- apc_scores %>%
  select(donor_id, wnt_rspo, genotype, route_score, x) %>%
  pivot_wider(names_from = genotype, values_from = c(route_score, x))
donor_colours <- c(Donor1 = COL[["wnt"]], Donor2 = COL[["route"]], Donor3 = COL[["crc"]])
p5b <- ggplot(apc_scores, aes(genotype, route_score, colour = donor_id)) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = COL[["neutral_light"]]) +
  geom_segment(
    data = apc_segments,
    aes(x = 1, xend = 2, y = route_score_WT, yend = route_score_APC, colour = donor_id),
    inherit.aes = FALSE, linewidth = 0.45, alpha = 0.78
  ) +
  geom_point(size = 1.75, stroke = 0.42) +
  facet_wrap(~wnt_rspo, labeller = as_labeller(c(with = "+WNT/RSPO", without = "−WNT/RSPO"))) +
  scale_x_discrete(labels = c(WT = "WT", APC = "APC-KO")) +
  scale_colour_manual(values = donor_colours) +
  labs(x = NULL, y = "Epithelial programme score", colour = NULL) +
  theme_jtm() +
  theme(
    axis.text.x = element_text(size = 5.6), legend.position = "top",
    legend.justification = "left", strip.text = element_text(size = 5.9)
  )

perturbation_labels <- c(
  "GSE125472|APC_vs_WT_with_Wnt" = "GSE125472 APC-KO (+WNT)",
  "GSE125472|APC_vs_WT_without_Wnt" = "GSE125472 APC-KO (−WNT)",
  "GSE171910|conditional_wnt_silencing" = "GSE171910 WNT off",
  "GSE130822|ascl2_ko_vs_resting_wt" = "GSE130822 ASCL2-KO",
  "GSE135328_HT29|tcf7l2_ko_vs_wt" = "HT29 TCF7L2-KO",
  "GSE135328_HCT116|tcf7l2_ko_vs_wt" = "HCT116 TCF7L2-KO",
  "GSE114059|trametinib_vs_dmso" = "GSE114059 trametinib",
  "GSE114059|pri724_reversal_of_trametinib" = "GSE114059 PRI-724 reversal",
  "GSE67186|apc_restoration_shApc" = "GSE67186 Apc restored",
  "GSE67186|apc_restoration_shApc_Kras" = "GSE67186 Apc/Kras restored"
)
forest_summary <- env_closure$forest_summary %>%
  mutate(
    context = factor(
      key, levels = rev(names(perturbation_labels)), labels = rev(unname(perturbation_labels))
    ),
    status_label = factor(
      status,
      levels = c("supportive_specific", "supportive_direction_only", "discordant", "exploratory_low_coverage"),
      labels = c("Matched", "Direction", "Discordant", "Low coverage")
    )
  )
forest_units <- env_closure$forest_units %>%
  mutate(context = factor(
    key, levels = rev(names(perturbation_labels)), labels = rev(unname(perturbation_labels))
  ))
p5c <- ggplot(forest_summary, aes(y = context)) +
  geom_hline(
    yintercept = seq_along(levels(forest_summary$context)),
    linewidth = 0.26, colour = COL[["neutral_pale"]]
  ) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = COL[["neutral"]]) +
  geom_errorbar(aes(xmin = aligned_low, xmax = aligned_high), orientation = "y",
                width = 0.18, linewidth = 0.50, colour = COL[["ink"]]) +
  geom_point(data = forest_units, aes(x = aligned, y = context), inherit.aes = FALSE,
             position = position_jitter(height = 0.08, width = 0, seed = 20260808),
             size = 1.05, colour = COL[["neutral"]], alpha = 0.72) +
  geom_point(aes(x = aligned_mean, colour = status_label, shape = coverage),
             size = 1.95, stroke = 0.58) +
  scale_colour_manual(values = c(
    "Matched" = COL[["route"]], "Direction" = COL[["wnt"]],
    "Discordant" = COL[["uncertain"]], "Low coverage" = COL[["neutral"]]
  )) +
  scale_shape_manual(values = c("≥80% route coverage" = 18, "<80% route coverage" = 23)) +
  labs(x = "Direction-aligned programme effect", y = NULL, colour = NULL, shape = NULL) +
  theme_jtm(base_size = 6.5) +
  theme(
    axis.text.y = element_text(size = 5.05), legend.position = "top",
    legend.justification = "left", legend.text = element_text(size = 4.7),
    legend.key.width = grid::unit(2.4, "mm"),
    legend.spacing.x = grid::unit(0.25, "mm")
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE), shape = "none")

component_short_labels <- c(
  "Locked\nroute" = "Full\nprogramme",
  "Up\narm" = "Up\narm",
  "Differentiation\nloss" = "Diff.\nloss",
  "WNT/\nstemness" = "WNT/\nstem.",
  "Proliferation\nco-movement" = "Prolif."
)
component_v11 <- env_closure$component_heat %>%
  mutate(
    context_short = factor(
      gsub("\n", " ", as.character(context)),
      levels = gsub("\n", " ", levels(context))
    ),
    component_short = factor(
      component_short_labels[as.character(component)],
      levels = unname(component_short_labels)
    )
  )
p5d <- ggplot(component_v11, aes(component_short, context_short, fill = aligned)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = COL[["wnt"]], mid = "white", high = COL[["route"]], midpoint = 0,
    limits = c(-1.6, 1.6), oob = squish, name = "Aligned effect"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.4, base_family = JTM_FONT) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 4.8, angle = 0, hjust = 0.5, lineheight = 0.86),
    axis.text.y = element_text(size = 4.65, lineheight = 0.84),
    legend.position = "bottom", legend.title = element_text(size = 5.2),
    legend.text = element_text(size = 4.8),
    legend.key.width = grid::unit(14, "mm"), legend.key.height = grid::unit(1.8, "mm"),
    plot.margin = margin(1.4, 1.6, 1.4, 1.6, "mm")
  ) +
  guides(fill = guide_colourbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = grid::unit(22, "mm"), barheight = grid::unit(1.8, "mm")
  ))

pharm_v11 <- env_closure$pharm %>%
  mutate(
    row_short = case_when(
      comparison == "trametinib_vs_dmso" & feature == "route_score" ~ "Trametinib · programme",
      comparison == "trametinib_vs_dmso" & feature == "wnt_stem" ~ "Trametinib · WNT",
      comparison == "pri724_reversal_of_trametinib" & feature == "route_score" ~ "PRI-724 rescue · programme",
      comparison == "pri724_reversal_of_trametinib" & feature == "wnt_stem" ~ "PRI-724 rescue · WNT"
    ),
    row_short = factor(
      row_short,
      levels = rev(c(
        "Trametinib · programme", "Trametinib · WNT",
        "PRI-724 rescue · programme", "PRI-724 rescue · WNT"
      ))
    )
  )
pharm_v11_summary <- pharm_v11 %>%
  group_by(row_short, endpoint) %>%
  summarise(mean_aligned = mean(aligned), .groups = "drop")
p5e <- ggplot(pharm_v11, aes(aligned, row_short, colour = endpoint)) +
  geom_vline(xintercept = 0, linewidth = 0.32, colour = COL[["neutral"]]) +
  geom_point(
    size = 1.35, alpha = 0.82,
    position = position_jitter(height = 0.07, width = 0, seed = 20260808)
  ) +
  geom_point(
    data = pharm_v11_summary,
    aes(x = mean_aligned, y = row_short, fill = endpoint),
    inherit.aes = FALSE, shape = 23, size = 2.15, stroke = 0.40,
    colour = COL[["ink"]]
  ) +
  scale_colour_manual(values = c(
    "Locked route" = COL[["route"]], "WNT/stemness" = COL[["wnt"]]
  ), labels = c("Locked route" = "Programme", "WNT/stemness" = "WNT/stemness")) +
  scale_fill_manual(values = c(
    "Locked route" = COL[["route"]], "WNT/stemness" = COL[["wnt"]]
  ), labels = c("Locked route" = "Programme", "WNT/stemness" = "WNT/stemness")) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0)) +
  labs(x = "Direction-aligned effect", y = NULL) +
  theme_jtm(base_size = 6.4) +
  theme(
    axis.text.y = element_text(size = 4.75), legend.position = "none",
    plot.margin = margin(1.4, 1.6, 1.4, 1.6, "mm")
  )
p5f <- ggplot(env_closure$tcf_grid, aes(method, context, fill = scaled)) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = label), size = 1.50, colour = COL[["ink"]], family = JTM_FONT) +
  scale_fill_gradient2(
    low = COL[["wnt"]], mid = "white", high = COL[["route"]], midpoint = 0,
    limits = c(-1, 1), na.value = COL[["neutral_pale"]], guide = "none"
  ) +
  scale_x_discrete(labels = c(
    "Locked route" = "Programme", "ULM activity" = "ULM",
    "Signed-mean activity" = "Signed", "Edge deletion" = "Graph"
  )) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.4, base_family = JTM_FONT) +
  theme(
    panel.grid = element_blank(), legend.position = "none",
    axis.text.x = element_text(size = 4.55, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 4.9),
    plot.margin = margin(1.4, 1.6, 1.4, 1.6, "mm")
  )

fig5_row2 <- (clean_panel(p5b) | clean_panel(p5c)) +
  plot_layout(widths = c(0.82, 1.18))
fig5_row3 <- (clean_panel(p5d) | clean_panel(p5e) | clean_panel(p5f)) +
  plot_layout(widths = c(1.14, 1.03, 0.83))
fig5 <- clean_panel(p5a) / fig5_row2 / fig5_row3 +
  plot_layout(heights = c(2.05, 1.28, 1.16)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(mechanism_nodes, "figure5a_apc_wnt_mechanism_nodes.tsv")
write_tsv(mechanism_interventions, "figure5a_perturbation_readouts.tsv")
write_tsv(apc_scores, "figure5b_apc_organoid_scores.tsv")

# -----------------------------------------------------------------------------
# Figure 6: spatial localisation and protein triangulation (7 panels)
# -----------------------------------------------------------------------------

spatial_spots <- read_tsv(
  "results/computational_closure_validation/spatial_locked_route_spot_scores.tsv.gz"
) %>% filter(in_tissue == 1)
eligible_maps <- spatial_spots %>%
  filter(pathology_group %in% c("tumor", "non_neoplastic_epithelium")) %>%
  count(sample_id, pathology_group) %>%
  count(sample_id, name = "eligible_groups") %>%
  filter(eligible_groups == 2)
representative_sample <- spatial_spots %>%
  filter(sample_id %in% eligible_maps$sample_id,
         pathology_group %in% c("tumor", "non_neoplastic_epithelium")) %>%
  count(sample_id, name = "eligible_spots") %>%
  arrange(desc(eligible_spots), sample_id) %>%
  slice(1) %>% pull(sample_id)
map_data <- spatial_spots %>% filter(sample_id == representative_sample)

pathology_levels <- c(
  "tumor", "non_neoplastic_epithelium", "tumor_stroma", "stroma",
  "submucosa", "exclude"
)
pathology_labels <- c(
  tumor = "Tumour", non_neoplastic_epithelium = "Non-neoplastic epithelium",
  tumor_stroma = "Tumour stroma", stroma = "Stroma",
  submucosa = "Submucosa", exclude = "Excluded/other"
)
pathology_colours <- c(
  tumor = COL[["route"]], non_neoplastic_epithelium = COL[["neutral"]],
  tumor_stroma = COL[["context"]], stroma = COL[["crc"]],
  submucosa = COL[["adenoma"]], exclude = COL[["neutral_light"]]
)
map_data <- map_data %>%
  mutate(pathology_display = factor(
    pathology_group, levels = pathology_levels, labels = pathology_labels[pathology_levels]
  ))

map_theme <- theme_void(base_size = 7, base_family = JTM_FONT) +
  theme(
    text = element_text(family = JTM_FONT, colour = COL[["ink"]]),
    legend.position = "bottom", legend.title = element_text(size = 6.2, face = "bold"),
    legend.text = element_text(size = 5.25), legend.key.height = grid::unit(2.0, "mm"),
    legend.key.width = grid::unit(3.0, "mm"),
    legend.spacing.x = grid::unit(0.55, "mm"),
    plot.margin = margin(1.2, 1.6, 1.2, 1.6, "mm"),
    plot.tag = element_text(size = 8.2, face = "bold", family = JTM_FONT)
  )

p6a <- ggplot(map_data, aes(pxl_col, -pxl_row, colour = pathology_display)) +
  geom_point(size = 0.70, alpha = 0.92) +
  scale_colour_manual(values = setNames(pathology_colours, pathology_labels[names(pathology_colours)]),
                      drop = FALSE, name = NULL) +
  coord_equal() + map_theme +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1.65)))

score_limit <- quantile(
  abs(c(map_data$route_score, map_data$route_residual_prolif_epithelial)),
  probs = 0.98, na.rm = TRUE
)
make_spatial_score <- function(feature, legend_title) {
  ggplot(map_data, aes(pxl_col, -pxl_row, colour = .data[[feature]])) +
    geom_point(size = 0.74, alpha = 0.94, na.rm = TRUE) +
    scale_colour_gradient2(
      low = COL[["wnt"]], mid = "white", high = COL[["route"]], midpoint = 0,
      limits = c(-score_limit, score_limit), oob = squish, name = legend_title
    ) +
    coord_equal() + map_theme +
    guides(colour = guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(22, "mm"), barheight = grid::unit(2.2, "mm")
    ))
}
p6b <- make_spatial_score("route_score", "Raw programme")
p6c <- make_spatial_score("route_residual_prolif_epithelial", "Adjusted programme")

spatial_main <- env_closure$spatial %>%
  mutate(
    feature_short = factor(
      as.character(feature_label),
      levels = c("Raw route\nP = 0.0313", "Adjusted route\nP = 0.0313"),
      labels = c("Raw", "Adjusted")
    )
  )
spatial_main_median <- spatial_main %>%
  group_by(feature_short) %>% summarise(median_difference = median(difference), .groups = "drop")
p6d <- ggplot(spatial_main, aes(feature_short, difference, group = sample_id)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = COL[["neutral"]]) +
  geom_line(colour = COL[["neutral_light"]], linewidth = 0.48) +
  geom_point(aes(colour = feature_short), size = 1.45) +
  geom_point(data = spatial_main_median, aes(feature_short, median_difference),
             inherit.aes = FALSE, shape = 18, size = 2.5, colour = COL[["ink"]]) +
  scale_colour_manual(values = c(Raw = COL[["route_light"]], Adjusted = COL[["route"]])) +
  labs(x = NULL, y = "Tumour − non-neoplastic\nepithelium") + theme_jtm() +
  theme(axis.text.x = element_text(size = 5.8), legend.position = "none")

pxd2137_main <- env_gap$pxd2137_plot %>%
  mutate(gene_display = factor(gene, levels = rev(c("OLFM4", "FABP1", "SOX9", "CTNNB1", "ETHE1"))))
p6e <- ggplot(pxd2137_main, aes(y = gene_display, colour = role)) +
  geom_hline(
    yintercept = seq_along(levels(pxd2137_main$gene_display)),
    linewidth = 0.25, colour = COL[["neutral_pale"]]
  ) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(
    data = filter(pxd2137_main, is.finite(age_sex_adjusted_log2_effect)),
    aes(x = age_sex_adjusted_ci_low, xend = age_sex_adjusted_ci_high, yend = gene_display),
    linewidth = 0.65
  ) +
  geom_point(
    data = filter(pxd2137_main, is.finite(age_sex_adjusted_log2_effect)),
    aes(x = age_sex_adjusted_log2_effect), size = 1.9
  ) +
  geom_text(
    data = filter(pxd2137_main, !is.finite(age_sex_adjusted_log2_effect)),
    aes(x = 0.15, label = "not estimable"), hjust = 0, size = 1.55, colour = COL[["neutral"]]
  ) +
  scale_colour_manual(values = c(
    "Direct locked arm" = COL[["route"]], "Context marker" = COL[["context"]],
    "Reserve locked marker" = COL[["wnt"]]
  ), labels = c("Direct locked arm" = "Direct", "Context marker" = "Context",
                "Reserve locked marker" = "Reserve")) +
  coord_cartesian(xlim = c(-2.0, 5.8), clip = "off") +
  labs(x = "Adjusted log2 effect", y = NULL, colour = NULL) + theme_jtm(base_size = 6.4) +
  theme(legend.position = "none", axis.text.y = element_text(size = 5.6))

pxd445_main <- env_gap$pxd445_plot %>%
  mutate(gene_display = factor(gene, levels = rev(c("OLFM4", "FABP1", "SORD", "ETHE1"))))
p6f <- ggplot(pxd445_main,
              aes(median_adenoma_minus_normal_log2, gene_display)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_hline(
    yintercept = seq_along(levels(pxd445_main$gene_display)),
    linewidth = 0.25, colour = COL[["neutral_pale"]]
  ) +
  geom_point(
    aes(size = direction_concordant_fraction, fill = fdr_supported),
    shape = 21, stroke = 0.38, colour = COL[["ink"]], alpha = 0.96
  ) +
  scale_fill_manual(values = c(`TRUE` = COL[["route"]], `FALSE` = COL[["neutral"]]),
                    labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "Direction only")) +
  scale_size_continuous(range = c(1.5, 2.9), limits = c(0.5, 1),
                        labels = percent_format(accuracy = 1), name = "Concordant pairs") +
  coord_cartesian(xlim = c(-1.15, 1.10)) +
  labs(x = "Median paired log2 effect", y = NULL, fill = NULL) + theme_jtm(base_size = 6.4) +
  theme(legend.position = "none", axis.text.y = element_text(size = 5.6)) +
  guides(fill = guide_legend(nrow = 1), size = guide_legend(nrow = 1))

protein_matrix_v11 <- env_gap$matrix_plot %>%
  mutate(
    marker = gsub("locked up", "programme-up", marker, fixed = TRUE),
    marker = gsub("locked down", "programme-down", marker, fixed = TRUE),
    marker = gsub("reserve down", "reserve", marker, fixed = TRUE)
  )
protein_matrix_colours <- c(
  "FDR-supported difference" = "#E8B09E",
  "Direction only / imprecise" = "#F3D8CF",
  "Detectability only" = "#D7E5EE",
  "Not detected" = "#ECEFF1",
  "Not estimable" = "#F4F5F6",
  "Not evaluated" = "#F4F5F6",
  "Direct locked arm" = "#F1D39B",
  "Context marker" = "#DDD5E8",
  "Reserve locked marker" = "#D2E1EA"
)
p6g <- ggplot(protein_matrix_v11, aes(evidence_source, marker, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.42) +
  geom_text(
    aes(label = display), size = 1.68, lineheight = 0.88,
    colour = COL[["ink"]], family = JTM_FONT
  ) +
  scale_fill_manual(values = protein_matrix_colours, guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.4, base_family = JTM_FONT) +
  theme(
    panel.grid = element_blank(),
    axis.text.x.top = element_text(size = 5.05, lineheight = 0.88, colour = COL[["ink"]]),
    axis.text.y = element_text(size = 5.4, colour = COL[["ink"]]),
    axis.title = element_blank(),
    plot.margin = margin(1.2, 1.7, 1.2, 1.7, "mm")
  )

fig6_row2 <- (clean_panel(p6d) | clean_panel(p6e) | clean_panel(p6f)) +
  plot_layout(widths = c(0.76, 1.12, 1.12))
fig6 <- (clean_panel(p6a) | clean_panel(p6b) | clean_panel(p6c)) /
  fig6_row2 /
  clean_panel(p6g) +
  plot_layout(heights = c(1.03, 1.02, 0.84), widths = c(0.98, 1.01, 1.01)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(map_data, "figure6a-c_representative_spatial_spots.tsv")
write_tsv(protein_matrix_v11, "figure6g_protein_evidence_matrix.tsv")

# -----------------------------------------------------------------------------
# Supplementary figures: supporting sensitivity and audit detail
# -----------------------------------------------------------------------------

figS1 <- clean_panel(env_main$figS1) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
supp_s2c <- clean_panel(env_main$pS2c) +
  labs(y = "Epithelial programme score")
supp_s2d <- clean_panel(env_main$pS2d) +
  labs(y = "Epithelial programme score")
figS2 <- (clean_panel(env_main$pS2a) | clean_panel(env_main$pS2b)) /
  (supp_s2c | supp_s2d) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

supp_s3a <- clean_panel(env_main$pS3a) +
  labs(y = "Patient-median epithelial programme")
supp_s3c <- clean_panel(env_main$pS3c) +
  scale_colour_manual(
    values = c(
      "Locked route" = COL[["route"]],
      "WNT/stemness" = COL[["wnt"]],
      "Proliferation" = COL[["neutral"]]
    ),
    labels = c(
      "Locked route" = "Programme",
      "WNT/stemness" = "WNT/stemness",
      "Proliferation" = "Proliferation"
    )
  ) +
  facet_wrap(
    ~locus, nrow = 1,
    labeller = as_labeller(c(
      "WNT-route loci" = "WNT-programme loci",
      "TCF/ASCL2-axis loci" = "TCF/ASCL2-axis loci"
    ))
  )
figS3 <- (supp_s3a | clean_panel(env_main$pS3b)) / supp_s3c +
  plot_annotation(tag_levels = "a", theme = tag_theme)

atlas_influence <- read_tsv("results/figure_data_locked/figs2_atlas_study_influence.tsv")
state_levels <- c(
  "normal_epithelial", "polyp_epithelial", "polyp_cancer", "tumor_epithelial",
  "tumor_cancer", "metastasis_epithelial", "metastasis_cancer"
)
state_short <- c(
  "Polyp epithelium", "Polyp cancer", "Primary epithelium",
  "Primary cancer", "Metastatic epithelium", "Metastatic cancer"
)
atlas_influence_plot <- atlas_influence %>%
  mutate(
    state_label = factor(state, levels = rev(state_levels[-1]), labels = rev(state_short)),
    programme = ifelse(outcome == "score__ca_route_signature", "Epithelial programme", "WNT/stemness"),
    min_label = paste0("min n=", minimum_target_donors_after_omission)
  )
make_atlas_influence <- function(programme_name, colour) {
  part <- atlas_influence_plot %>% filter(programme == programme_name)
  right_edge <- max(part$loo_max_coef) + 0.14
  ggplot(part, aes(full_coef, state_label)) +
    geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
    geom_segment(aes(x = loo_min_coef, xend = loo_max_coef, yend = state_label),
                 linewidth = 1.5, colour = COL[["neutral_light"]], lineend = "round") +
    geom_segment(aes(x = full_ci_low, xend = full_ci_high, yend = state_label),
                 linewidth = 0.62, colour = colour) +
    geom_point(size = 1.8, colour = colour) +
    geom_text(aes(x = right_edge, label = min_label), hjust = 1, size = 1.7,
              colour = COL[["neutral"]]) +
    coord_cartesian(xlim = c(-0.02, right_edge + 0.02), clip = "off") +
    labs(x = paste0(programme_name, " adjusted coefficient"), y = NULL) +
    theme_jtm()
}
pS4b <- make_atlas_influence("Epithelial programme", COL[["route"]])
pS4c <- make_atlas_influence("WNT/stemness", COL[["wnt"]])
figS4 <- clean_panel(env_main$pS4a) | (clean_panel(pS4b) / clean_panel(pS4c)) +
  plot_layout(widths = c(1.08, 0.92)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

feature_coverage <- read_tsv("results/computational_closure_validation/feature_coverage.tsv") %>%
  mutate(
    dataset = factor(dataset, levels = rev(unique(dataset))),
    feature = factor(feature, levels = c(
      "route_up", "route_down", "wnt_stem", "proliferation_control", "epithelial_control"
    ), labels = c("Up arm", "Down arm", "WNT/stem", "Proliferation", "Epithelial")),
    label = percent(coverage_fraction, accuracy = 1)
  )
pS5a <- ggplot(feature_coverage, aes(feature, dataset, fill = coverage_fraction)) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = label), size = 1.65) +
  scale_fill_gradient(low = "#F3E4DF", high = COL[["route"]], limits = c(0, 1),
                      name = "Coverage") +
  labs(x = NULL, y = NULL) + theme_jtm() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 28, hjust = 1), legend.position = "top")

evidence <- read_tsv("results/computational_closure_validation/evidence_closure_matrix.tsv") %>%
  filter(layer == "empirical_perturbation", expected_direction != 0) %>%
  mutate(
    context = factor(paste(dataset, comparison, sep = "\n"),
                     levels = rev(paste(dataset, comparison, sep = "\n"))),
    status = factor(status, levels = c(
      "supportive_specific", "supportive_direction_only", "discordant", "exploratory_low_coverage"
    ), labels = c("Expression-matched", "Direction only", "Discordant", "Low coverage")),
    minus_log10_p = -log10(specificity_p)
  )
pS5b <- ggplot(evidence, aes(minus_log10_p, context, colour = status)) +
  geom_vline(xintercept = -log10(0.05), linetype = "22", linewidth = 0.35,
             colour = COL[["neutral"]]) +
  geom_segment(aes(x = 0, xend = minus_log10_p, yend = context), linewidth = 0.55) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(
    "Expression-matched" = COL[["route"]], "Direction only" = COL[["wnt"]],
    "Discordant" = COL[["uncertain"]], "Low coverage" = COL[["neutral"]]
  )) +
  labs(x = expression(-log[10](P[matched])), y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.5) +
  theme(axis.text.y = element_text(size = 5.0, lineheight = 0.86), legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

virtual_panel <- read_tsv("results/computational_closure_validation/virtual_tf_knockout_panel.tsv") %>%
  mutate(
    n_paths = n_direct_paths + n_two_hop_paths,
    tf = factor(tf, levels = tf[order(combined_virtual_ko_route_impact)])
  )
pS5c <- ggplot(virtual_panel, aes(combined_virtual_ko_route_impact, tf)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = COL[["neutral"]]) +
  geom_segment(aes(x = 0, xend = combined_virtual_ko_route_impact, yend = tf),
               linewidth = 0.55, colour = COL[["neutral"]]) +
  geom_point(aes(size = n_paths, fill = tf == "TCF7L2"), shape = 21,
             colour = COL[["ink"]], stroke = 0.35) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = COL[["route"]]), guide = "none") +
  scale_size_continuous(range = c(1.5, 3.5), name = "Signed paths") +
  labs(x = "Virtual-KO programme impact", y = NULL) + theme_jtm() +
  theme(legend.position = "top")

ulm_coverage <- read_tsv("results/computational_closure_validation/collectri_ulm_target_coverage.tsv") %>%
  mutate(dataset = factor(dataset, levels = rev(unique(dataset))))
pS5d <- ggplot(ulm_coverage, aes(tf, dataset, fill = log10(n_targets + 1))) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = n_targets), size = 1.65) +
  scale_fill_gradient(low = "#E7EEF3", high = COL[["wnt"]], name = "log10 targets") +
  labs(x = NULL, y = NULL) + theme_jtm() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")

figS5 <- (clean_panel(pS5a) | clean_panel(pS5b)) /
  (clean_panel(pS5c) | clean_panel(pS5d)) +
  plot_layout(heights = c(1.1, 0.9)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(feature_coverage, "figureS5a_perturbation_feature_coverage.tsv")
write_tsv(evidence, "figureS5b_expression_matched_support.tsv")
write_tsv(virtual_panel, "figureS5c_virtual_tf_deletion.tsv")
write_tsv(ulm_coverage, "figureS5d_regulon_target_coverage.tsv")

spatial_effects <- read_tsv(
  "results/computational_closure_validation/spatial_locked_route_section_effects.tsv"
) %>%
  filter(comparison == "tumor_vs_non_neoplastic_epithelium") %>%
  mutate(
    feature_label = recode(
      feature, route_score = "Raw programme", route_residual_prolif_epithelial = "Adjusted programme",
      wnt_stem = "WNT/stemness", proliferation_control = "Proliferation",
      epithelial_control = "Epithelial control"
    ),
    feature_label = factor(feature_label, levels = c(
      "Raw programme", "Adjusted programme", "WNT/stemness", "Proliferation", "Epithelial control"
    ))
  )
pS6a <- ggplot(spatial_effects, aes(difference, sample_id, colour = feature_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = COL[["neutral"]]) +
  geom_point(size = 1.45, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(
    "Raw programme" = COL[["route"]], "Adjusted programme" = "#8F3D2A",
    "WNT/stemness" = COL[["wnt"]], "Proliferation" = COL[["crc"]],
    "Epithelial control" = COL[["context"]]
  ), drop = FALSE) +
  labs(x = "Tumour − non-neoplastic epithelium", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.4) +
  theme(axis.text.y = element_text(size = 5.2), legend.position = "top")

eligible_sample_labels <- setNames(
  paste("Section", seq_along(unique(spatial_effects$sample_id))),
  unique(spatial_effects$sample_id)
)
pathology_summary <- read_tsv(
  "results/computational_closure_validation/spatial_locked_route_pathology_summary.tsv"
) %>%
  select(sample_id, pathology_group, route_score, route_residual_prolif_epithelial,
         wnt_stem, proliferation_control) %>%
  filter(sample_id %in% unique(spatial_effects$sample_id)) %>%
  pivot_longer(c(route_score, route_residual_prolif_epithelial, wnt_stem, proliferation_control),
               names_to = "feature", values_to = "median_score") %>%
  mutate(
    feature = factor(feature, levels = c(
      "route_score", "route_residual_prolif_epithelial", "wnt_stem", "proliferation_control"
    ), labels = c("Raw programme", "Adjusted programme", "WNT/stemness", "Proliferation")),
    pathology_group = factor(pathology_group),
    sample_label = factor(
      eligible_sample_labels[sample_id], levels = unname(eligible_sample_labels)
    )
  )
pS6b <- ggplot(pathology_summary, aes(feature, pathology_group, fill = median_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_wrap(~sample_label, ncol = 3) +
  scale_fill_gradient2(low = COL[["wnt"]], mid = "white", high = COL[["route"]], midpoint = 0,
                       oob = squish, name = "Median score") +
  labs(x = NULL, y = NULL) + theme_jtm(base_size = 5.8) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 4.7),
        axis.text.y = element_text(size = 4.7), strip.text = element_text(size = 4.8),
        legend.position = "top")

protein_matrix_s6 <- env_gap$matrix_plot %>%
  mutate(
    marker = recode(
      as.character(marker),
      "OLFM4 · locked up" = "OLFM4 · programme-up",
      "FABP1 · locked down" = "FABP1 · programme-down",
      "ETHE1 · reserve down" = "ETHE1 · reserve"
    ),
    marker = factor(
      marker,
      levels = rev(c(
        "OLFM4 · programme-up", "FABP1 · programme-down", "SOX9 · context",
        "CTNNB1 · context", "ETHE1 · reserve"
      ))
    )
  )
pS6c_v18 <- ggplot(protein_matrix_s6, aes(evidence_source, marker, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = display), size = 1.80, colour = COL[["ink"]]) +
  scale_fill_manual(values = protein_matrix_colours, guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_jtm(base_size = 6.1) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x.top = element_text(size = 5.0, lineheight = 0.90),
    axis.text.y = element_text(size = 5.5),
    panel.background = element_rect(fill = COL[["neutral_pale"]], colour = NA)
  )
figS6 <- clean_panel(pS6a) / (clean_panel(pS6b) | clean_panel(pS6c_v18)) +
  plot_layout(heights = c(0.72, 1.28), widths = c(1.08, 0.92)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(spatial_effects, "figureS6a_spatial_section_controls.tsv")
write_tsv(pathology_summary, "figureS6b_spatial_pathology_summary.tsv")

# -----------------------------------------------------------------------------
# Export, traceability and numerical/format QA
# -----------------------------------------------------------------------------

exports <- bind_rows(
  export_figure(fig1, "figure1_discovery_lock_and_heldout_validation", 170, 180),
  export_figure(fig2, "figure2_multicohort_and_ffpe_replication", 170, 178),
  export_figure(fig3, "figure3_cross_modal_rna_atac_support", 170, 142),
  export_figure(fig4, "figure4_crc_atlas_recapitulation", 170, 160),
  export_figure(fig5, "figure5_perturbation_validation", 170, 210),
  export_figure(fig6, "figure6_spatial_and_protein_triangulation", 170, 165),
  export_figure(figS1, "figureS1_programme_and_platform_coverage", 170, 112),
  export_figure(figS2, "figureS2_external_validation_sensitivity", 170, 132),
  export_figure(figS3, "figureS3_rna_atac_robustness", 170, 124),
  export_figure(figS4, "figureS4_atlas_source_audit", 170, 140),
  export_figure(figS5, "figureS5_perturbation_boundary_and_coverage", 170, 140),
  export_figure(figS6, "figureS6_spatial_and_protein_assayability", 170, 145)
)
exports$file_size_bytes <- file.info(exports$file)$size
exports$sha256 <- vapply(
  exports$file,
  function(path) strsplit(system2("sha256sum", path, stdout = TRUE), "[[:space:]]+")[[1]][1],
  character(1)
)
write_tsv(exports, "figure_export_manifest.tsv")

panel_trace <- data.frame(
  figure = rep(c("Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Figure 6"),
               c(5, 6, 6, 5, 6, 7)),
  panel = unlist(lapply(c(5, 6, 6, 5, 6, 7), function(n) letters[seq_len(n)])),
  claim = c(
    "Study workflow and biological conclusion", "Discovery stability", "Biological anchor genes",
    "Held-out discrimination", "Donor-paired confirmation",
    "Five independent cohorts", "Pooled and adjusted robustness", "Histological specificity",
    "Patient-paired FFPE validation", "FFPE robustness estimands", "Expanded validation analyses",
    "Becker tissue-state transfer", "Patient-clustered Becker model", "Matched RNA-ATAC coupling",
    "Adjusted RNA-from-ATAC models", "Clustered and fixed-effect sensitivity", "Regulatory-distance sensitivity",
    "CRC Atlas state distributions", "Adjusted carrier-state effects", "Study-omission ranges",
    "Source-study support", "Within-study primary-cancer contrasts",
    "APC–WNT–TCF/ASCL2 mechanistic evidence map", "APC-KO donor-paired organoids", "Cross-perturbation programme effects",
    "Component concordance", "PDO pharmacological boundary", "TCF7L2 regulatory consensus",
    "Representative pathology map", "Raw spatial programme", "Adjusted spatial programme",
    "Six-section paired effects", "PXD002137 abundance", "PXD000445 paired sensitivity",
    "Protein evidence roles"
  ),
  independent_unit = c(
    "evidence layer", "gene", "gene", "specimen", "donor pair",
    "sample or patient cluster", "sample or patient cluster", "patient cluster",
    "patient pair", "patient pair", "analysis within recruitment cluster",
    "sample", "patient cluster", "matched sample and patient cluster",
    "matched sample and patient cluster", "matched sample and patient", "matched sample",
    "donor-carrier", "donor cluster", "donor cluster", "donor-carrier", "donor-carrier",
    "biological mechanism node", "organoid donor", "experimental context", "experimental context",
    "patient-derived organoid", "clone, cell line or regulatory graph",
    "Visium spot", "Visium spot", "Visium spot", "Visium section",
    "specimen", "patient pair", "protein marker by resource"
  ),
  source_group = c(
    rep("results/figure_data_locked", 5), rep("transcriptomic validation result tables", 6),
    rep("results/figure_data_locked", 6), rep("CRC Atlas locked result tables", 5),
    rep("results/computational_closure_validation", 6),
    rep("spatial and public proteomic result tables", 7)
  ),
  stringsAsFactors = FALSE
)
write_tsv(panel_trace, "panel_source_trace.tsv")

qa <- data.frame(
  check = c(
    "main_figure_count", "main_panel_count", "supplementary_figure_count",
    "formats_per_figure", "all_widths_170mm", "all_heights_within_225mm",
    "all_exports_nonempty", "all_submission_files_under_10mb",
    "locked_signature_genes", "apc_organoid_donors", "ffpe_patient_pairs",
    "atlas_eligible_omissions", "spatial_positive_sections", "representative_spatial_sample"
  ),
  observed = c(
    n_distinct(exports$figure[grepl("^figure[1-6]_", exports$figure)]),
    nrow(panel_trace),
    n_distinct(exports$figure[grepl("^figureS", exports$figure)]),
    min(table(exports$figure)), all(exports$width_mm == 170), all(exports$height_mm <= 225),
    all(exports$file_size_bytes > 0),
    all(exports$file_size_bytes[exports$format %in% c("PDF", "TIFF")] < 10 * 1024^2),
    nrow(signature), n_distinct(apc_scores$donor_id),
    n_distinct(env_gap$paired_long$patient_id), unique(atlas_influence$n_eligible_omissions),
    sum(env_closure$spatial$difference[env_closure$spatial$feature_label == "Raw route\nP = 0.0313"] > 0),
    representative_sample
  ),
  expected = c(6, 35, 6, 4, TRUE, TRUE, TRUE, TRUE, 100, 3, 51, 33, 6, representative_sample),
  stringsAsFactors = FALSE
)
qa$pass <- mapply(function(observed, expected) {
  if (expected %in% c("TRUE", "FALSE")) return(as.character(observed) == expected)
  suppressWarnings({
    observed_num <- as.numeric(observed)
    expected_num <- as.numeric(expected)
  })
  if (is.finite(observed_num) && is.finite(expected_num)) return(abs(observed_num - expected_num) < 1e-8)
  as.character(observed) == as.character(expected)
}, qa$observed, qa$expected)
write_tsv(qa, "figure_build_qc.tsv")
if (!all(qa$pass)) stop("One or more v1.8 figure QA checks failed.")

message("Wrote six submission-refined main and six supplementary JTM figures to: ", OUT_DIR)
