#!/usr/bin/env Rscript

# Communications Biology Figure 1 v1.9.
#
# Panels a and b use an original, literature-informed vector design: a compact
# discovery-freeze-validation study spine and a true vertical gene-attrition
# waterfall. Quantitative results are read from locked source-data exports and
# are not recomputed here.

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
  file.path(getwd(), "analysis", "revise_communications_biology_figure1_v1_2.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
SOURCE_PACKAGE <- file.path(ROOT, "figures", "communications_biology_v1.1")
OUT_DIR <- file.path(ROOT, "figures", "communications_biology_v1.2")
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

JOURNAL_FONT <- "Arial"
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
  context = "#7C6AA3",
  pass = "#4D8C75"
)

read_source <- function(filename) {
  read.delim(
    file.path(SOURCE_PACKAGE, "source_data", filename), sep = "\t",
    stringsAsFactors = FALSE, check.names = FALSE, quote = "",
    comment.char = ""
  )
}

write_source <- function(frame, filename) {
  frame[] <- lapply(frame, function(column) {
    if (is.character(column)) {
      column <- gsub("[\r\n]+", "; ", column)
      column <- gsub("[[:space:]]+", " ", column)
      trimws(column)
    } else {
      column
    }
  })
  write.table(
    frame, file.path(SOURCE_DIR, filename), sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )
}

theme_journal <- function(base_size = 7.0) {
  theme_classic(base_size = base_size, base_family = JOURNAL_FONT) +
    theme(
      text = element_text(family = JOURNAL_FONT, colour = COL[["ink"]]),
      axis.line = element_line(linewidth = 0.32, colour = COL[["ink"]]),
      axis.ticks = element_line(linewidth = 0.28, colour = COL[["ink"]]),
      axis.ticks.length = unit(0.9, "mm"),
      axis.title = element_text(size = base_size, colour = COL[["ink"]]),
      axis.text = element_text(size = base_size - 0.55, colour = COL[["ink"]]),
      panel.grid = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 0.35, face = "bold"),
      legend.text = element_text(size = base_size - 0.75),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.15, face = "bold"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(1.6, 1.8, 1.6, 1.8, unit = "mm"),
      plot.tag = element_text(size = 8.2, face = "bold", colour = COL[["ink"]])
    )
}

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
  plot.tag = element_text(
    size = 8.2, face = "bold", family = JOURNAL_FONT, colour = COL[["ink"]]
  )
)

round_box <- function(fill, border, radius_mm = 1.6, lwd = 0.85) {
  roundrectGrob(
    r = unit(radius_mm, "mm"),
    gp = gpar(fill = fill, col = border, lwd = lwd, linejoin = "round")
  )
}

lock_icon <- function(colour = COL[["ink"]]) {
  grobTree(
    xsplineGrob(
      x = c(0.29, 0.30, 0.50, 0.70, 0.71),
      y = c(0.54, 0.78, 0.88, 0.78, 0.54),
      shape = c(0, 1, 1, 1, 0), open = TRUE,
      gp = gpar(col = colour, lwd = 1.25)
    ),
    roundrectGrob(
      x = 0.5, y = 0.36, width = 0.58, height = 0.43,
      r = unit(0.08, "snpc"),
      gp = gpar(fill = "white", col = colour, lwd = 1.1)
    ),
    circleGrob(0.5, 0.39, 0.045, gp = gpar(fill = colour, col = NA)),
    segmentsGrob(0.5, 0.5, 0.35, 0.27, gp = gpar(col = colour, lwd = 1.1))
  )
}

resource_icon_grob <- function(type, colour = "#54758B", fill = "#E7F0F4") {
  line_gp <- gpar(
    col = colour, fill = NA, lwd = 1.05, linecap = "round",
    linejoin = "round"
  )
  fill_gp <- gpar(col = colour, fill = fill, lwd = 1.00)
  nucleus_gp <- gpar(col = NA, fill = alpha(colour, 0.68))

  if (type == "single_cell") {
    return(grobTree(
      circleGrob(0.34, 0.54, 0.20, gp = fill_gp),
      circleGrob(0.34, 0.54, 0.055, gp = nucleus_gp),
      circleGrob(0.65, 0.63, 0.155, gp = fill_gp),
      circleGrob(0.65, 0.63, 0.043, gp = nucleus_gp),
      circleGrob(0.62, 0.34, 0.125, gp = fill_gp),
      circleGrob(0.62, 0.34, 0.036, gp = nucleus_gp)
    ))
  }

  if (type == "bulk_ffpe") {
    return(grobTree(
      roundrectGrob(
        x = 0.50, y = 0.46, width = 0.37, height = 0.58,
        r = unit(0.08, "snpc"), gp = fill_gp
      ),
      rectGrob(
        x = 0.50, y = 0.80, width = 0.48, height = 0.12,
        gp = gpar(col = colour, fill = "white", lwd = 1.00)
      ),
      segmentsGrob(0.34, 0.42, 0.66, 0.42, gp = line_gp),
      xsplineGrob(
        x = c(0.39, 0.47, 0.54, 0.62),
        y = c(0.34, 0.50, 0.31, 0.51),
        shape = c(0, 1, 1, 0), open = TRUE, gp = line_gp
      )
    ))
  }

  if (type == "multiome") {
    return(grobTree(
      circleGrob(0.50, 0.51, 0.34, gp = fill_gp),
      segmentsGrob(0.50, 0.22, 0.50, 0.80, gp = line_gp),
      xsplineGrob(
        x = c(0.25, 0.34, 0.27, 0.39),
        y = c(0.32, 0.43, 0.57, 0.69),
        shape = c(0, 1, 1, 0), open = TRUE, gp = line_gp
      ),
      segmentsGrob(
        x0 = c(0.61, 0.73, 0.61, 0.61, 0.61),
        y0 = c(0.31, 0.31, 0.40, 0.50, 0.60),
        x1 = c(0.61, 0.73, 0.73, 0.73, 0.73),
        y1 = c(0.70, 0.70, 0.40, 0.50, 0.60),
        gp = line_gp
      )
    ))
  }

  if (type == "perturbation") {
    return(grobTree(
      segmentsGrob(0.16, 0.50, 0.84, 0.50, gp = line_gp),
      rectGrob(
        x = 0.54, y = 0.50, width = 0.15, height = 0.17,
        gp = gpar(col = NA, fill = "white")
      ),
      circleGrob(0.28, 0.29, 0.095, gp = fill_gp),
      circleGrob(0.28, 0.71, 0.095, gp = fill_gp),
      segmentsGrob(
        x0 = c(0.36, 0.36), y0 = c(0.34, 0.66),
        x1 = c(0.68, 0.61), y1 = c(0.62, 0.39), gp = line_gp
      ),
      circleGrob(0.68, 0.62, 0.026, gp = nucleus_gp),
      circleGrob(0.61, 0.39, 0.026, gp = nucleus_gp)
    ))
  }

  if (type == "spatial") {
    grid_x <- rep(seq(0.27, 0.73, length.out = 4), 4)
    grid_y <- rep(seq(0.28, 0.72, length.out = 4), each = 4)
    return(grobTree(
      circleGrob(
        x = grid_x, y = grid_y, r = 0.022,
        gp = gpar(col = NA, fill = alpha(colour, 0.55))
      ),
      polygonGrob(
        x = c(0.23, 0.35, 0.58, 0.77, 0.70, 0.48, 0.29),
        y = c(0.43, 0.73, 0.78, 0.59, 0.31, 0.22, 0.28),
        gp = gpar(col = colour, fill = alpha(fill, 0.66), lwd = 1.00)
      )
    ))
  }

  if (type == "proteomics") {
    return(grobTree(
      segmentsGrob(0.17, 0.25, 0.84, 0.25, gp = line_gp),
      segmentsGrob(
        x0 = c(0.25, 0.35, 0.45, 0.56, 0.67, 0.77),
        y0 = rep(0.25, 6),
        x1 = c(0.25, 0.35, 0.45, 0.56, 0.67, 0.77),
        y1 = c(0.44, 0.69, 0.52, 0.80, 0.60, 0.38),
        gp = line_gp
      ),
      circleGrob(
        x = c(0.31, 0.43, 0.55, 0.67),
        y = c(0.79, 0.70, 0.77, 0.69), r = 0.035,
        gp = gpar(col = colour, fill = "white", lwd = 0.85)
      ),
      segmentsGrob(
        x0 = c(0.345, 0.465, 0.585),
        y0 = c(0.765, 0.72, 0.745),
        x1 = c(0.395, 0.515, 0.635),
        y1 = c(0.725, 0.755, 0.71), gp = line_gp
      )
    ))
  }

  stop("Unknown resource icon type: ", type)
}

# -----------------------------------------------------------------------------
# Panel a: resource landscape, sample counts and fixed evidence roles
# -----------------------------------------------------------------------------

dataset_overview <- data.frame(
  order = 1:4,
  evidence_role = c(
    "Discovery and reduction", "Internal held-out test",
    "Independent replication", "Orthogonal support and tissue context"
  ),
  public_resources = c(1, 1, 7, 13),
  data_types = c(
    "epithelial scRNA-seq", "epithelial scRNA-seq",
    "bulk RNA, FFPE RNA and snRNA",
    "RNA–ATAC, perturbation, CRC atlas, spatial RNA and proteomics"
  ),
  independent_units = c(
    "60 specimens; 27 donors; 37 donor–tissue profiles",
    "44 specimens; 24 donors; 7 paired donors",
    "203 samples; 161 patient clusters; 51 FFPE pairs (102 specimens)",
    "40 RNA–ATAC samples/12 patients; 398 CRC Atlas donors; 6 spatial sections"
  ),
  purpose = c(
    "define the 287-gene core and select the 12-gene signature",
    "test the frozen core and signature within the source study",
    "test transfer across five cohorts, archival tissue and independent nuclei",
    "assess regulatory, perturbational, anatomical and protein-level support"
  ),
  signature_status = c(
    "gene selection permitted", "frozen before scoring",
    "frozen before scoring", "frozen before scoring"
  ),
  stringsAsFactors = FALSE
)
stopifnot(sum(dataset_overview$public_resources) == 22)

resource_landscape <- data.frame(
  order = 1:6,
  icon_key = c(
    "single_cell", "bulk_ffpe", "multiome", "perturbation", "spatial",
    "proteomics"
  ),
  resource_type = c(
    "sc/snRNA-seq", "Bulk and FFPE RNA", "Matched RNA–ATAC",
    "Perturbation RNA", "Spatial RNA", "Proteomics"
  ),
  displayed_scope = c(
    "Chen · Becker · CRC Atlas",
    "5 cohorts · 51 FFPE pairs",
    "40 samples · 12 patients",
    "6 datasets · 8 contrasts",
    "6 tissue sections",
    "4 datasets"
  ),
  stringsAsFactors = FALSE
)

resource_boundaries <- seq(0.25, 15.75, length.out = 7)
resource_card_centre <- head(resource_boundaries, -1) + diff(resource_boundaries) / 2

p1a <- ggplot() +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 7.28), clip = "off") +
  theme_void(base_family = JOURNAL_FONT) +
  theme(plot.margin = margin(1.0, 2.3, 0.8, 2.3, "mm")) +
  annotate(
    "text", x = 0.25, y = 7.03,
    label = "Public resources (n = 22)",
    hjust = 0, family = JOURNAL_FONT, fontface = "bold",
    size = 1.78, colour = COL[["ink"]]
  ) +
  annotate(
    "rect", xmin = 0.25, xmax = 15.75, ymin = 5.70, ymax = 6.78,
    fill = "#F7F8F9", colour = alpha(COL[["neutral"]], 0.52),
    linewidth = 0.34
  ) +
  annotate(
    "segment", x = resource_boundaries[2:6], xend = resource_boundaries[2:6],
    y = 5.70, yend = 6.78, colour = alpha(COL[["neutral"]], 0.34),
    linewidth = 0.28
  ) +
  annotate(
    "text", x = resource_card_centre, y = 6.13,
    label = resource_landscape$resource_type,
    family = JOURNAL_FONT, fontface = "bold", size = 1.55,
    colour = COL[["ink"]]
  ) +
  annotate(
    "text", x = resource_card_centre, y = 5.86,
    label = resource_landscape$displayed_scope,
    family = JOURNAL_FONT, size = 1.27,
    colour = COL[["neutral"]]
  ) +
  annotate(
    "text", x = 0.25, y = 5.28, label = "Analysis design",
    hjust = 0, family = JOURNAL_FONT, fontface = "bold",
    size = 1.64, colour = COL[["neutral"]]
  ) +
  annotate(
    "segment", x = 1.78, xend = 15.75, y = 5.28, yend = 5.28,
    colour = alpha(COL[["neutral"]], 0.34), linewidth = 0.30
  ) +
  annotation_custom(
    round_box("white", alpha(COL[["route"]], 0.82), 1.2, 0.90),
    xmin = 0.25, xmax = 3.25, ymin = 0.72, ymax = 4.96
  ) +
  annotate(
    "segment", x = 0.45, xend = 3.05, y = 4.83, yend = 4.83,
    colour = COL[["route"]], linewidth = 1.00, lineend = "round"
  ) +
  annotate(
    "text", x = 0.53, y = 4.48, label = "Discovery (1 resource)",
    hjust = 0, family = JOURNAL_FONT, fontface = "bold",
    size = 1.64, colour = COL[["route"]]
  ) +
  annotate(
    "text", x = 0.53, y = 3.83, label = "Epithelial scRNA-seq",
    hjust = 0, family = JOURNAL_FONT, fontface = "bold",
    size = 2.08, colour = COL[["ink"]]
  ) +
  annotate(
    "text", x = 0.53, y = 2.75,
    label = "60 specimens · 27 donors\n37 donor–tissue profiles",
    hjust = 0, vjust = 0.5, family = JOURNAL_FONT,
    size = 1.80, lineheight = 1.06, colour = COL[["ink"]]
  ) +
  annotate(
    "segment", x = 3.30, xend = 3.67, y = 2.78, yend = 2.78,
    colour = COL[["neutral"]], linewidth = 0.62,
    arrow = arrow(type = "closed", angle = 23, length = unit(1.30, "mm"))
  ) +
  annotation_custom(
    round_box("white", alpha(COL[["neutral"]], 0.66), 1.2, 0.82),
    xmin = 3.78, xmax = 6.68, ymin = 0.72, ymax = 4.96
  ) +
  annotate(
    "text", x = 4.06, y = 4.48, label = "Derivation",
    hjust = 0, family = JOURNAL_FONT, fontface = "bold",
    size = 1.64, colour = COL[["neutral"]]
  ) +
  annotate(
    "rect", xmin = 4.08, xmax = 6.38, ymin = 3.28, ymax = 3.90,
    fill = "#FBECE7", colour = alpha(COL[["route"]], 0.84), linewidth = 0.48
  ) +
  annotate(
    "text", x = 5.23, y = 3.59, label = "287-gene core",
    family = JOURNAL_FONT, fontface = "bold", size = 1.70,
    colour = COL[["route"]]
  ) +
  annotate(
    "segment", x = 5.23, xend = 5.23, y = 3.17, yend = 2.84,
    colour = COL[["neutral"]], linewidth = 0.52,
    arrow = arrow(type = "closed", angle = 23, length = unit(1.00, "mm"))
  ) +
  annotate(
    "rect", xmin = 4.08, xmax = 6.38, ymin = 2.12, ymax = 2.74,
    fill = "#FBF3DE", colour = alpha(COL[["adenoma"]], 0.88), linewidth = 0.48
  ) +
  annotate(
    "text", x = 5.23, y = 2.43, label = "12-gene signature",
    family = JOURNAL_FONT, fontface = "bold", size = 1.70,
    colour = COL[["adenoma"]]
  ) +
  annotate(
    "segment", x = 6.73, xend = 7.08, y = 2.78, yend = 2.78,
    colour = COL[["neutral"]], linewidth = 0.62,
    arrow = arrow(type = "closed", angle = 23, length = unit(1.30, "mm"))
  ) +
  annotate(
    "segment", x = 7.38, xend = 7.38, y = 0.68, yend = 4.96,
    colour = COL[["ink"]], linewidth = 0.56, linetype = "22"
  ) +
  annotation_custom(
    lock_icon(COL[["ink"]]), xmin = 7.15, xmax = 7.61,
    ymin = 4.66, ymax = 5.12
  ) +
  annotate(
    "text", x = 7.38, y = 0.42, label = "Frozen",
    family = JOURNAL_FONT, fontface = "bold", size = 1.37,
    colour = COL[["ink"]]
  ) +
  annotate(
    "segment", x = 7.52, xend = 8.18, y = 2.78, yend = 2.78,
    colour = COL[["neutral"]], linewidth = 0.62,
    arrow = arrow(type = "closed", angle = 23, length = unit(1.30, "mm"))
  ) +
  annotate(
    "segment", x = 8.27, xend = 8.27, y = 1.33, yend = 4.43,
    colour = alpha(COL[["neutral"]], 0.78), linewidth = 0.55
  )

for (i in seq_len(nrow(resource_landscape))) {
  p1a <- p1a + annotation_custom(
    resource_icon_grob(resource_landscape$icon_key[i]),
    xmin = resource_card_centre[i] - 0.28,
    xmax = resource_card_centre[i] + 0.28,
    ymin = 6.32,
    ymax = 6.72
  )
}

branch_specs <- data.frame(
  ymin = c(3.90, 2.23, 0.53),
  ymax = c(4.96, 3.75, 2.10),
  centre = c(4.43, 2.99, 1.32),
  fill = c("#F7FAFC", "#FCFAF5", "#FAF9FC"),
  border = c(COL[["wnt"]], COL[["adenoma"]], COL[["context"]]),
  heading = c(
    "Held-out test (1 partition)",
    "Independent RNA replication (7 resources)",
    "Multi-omic support (13 resources)"
  ),
  body = c(
    "44 specimens · 24 donors · 7 paired donors",
    "5 cohorts · 203 samples · 161 patient clusters\n51 FFPE pairs (102 specimens) · snRNA transfer",
    "Regulatory, perturbational and tissue-context support\nRNA–ATAC · perturbation · CRC Atlas · spatial RNA · proteomics"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(branch_specs))) {
  spec <- branch_specs[i, ]
  p1a <- p1a +
    annotate(
      "segment", x = 8.27, xend = 8.55, y = spec$centre, yend = spec$centre,
      colour = spec$border, linewidth = 0.58,
      arrow = arrow(type = "closed", angle = 23, length = unit(1.20, "mm"))
    ) +
    annotation_custom(
      round_box(spec$fill, alpha(spec$border, 0.70), 1.0, 0.78),
      xmin = 8.63, xmax = 15.75, ymin = spec$ymin, ymax = spec$ymax
    ) +
    annotate(
      "segment", x = 8.82, xend = 8.82,
      y = spec$ymin + 0.16, yend = spec$ymax - 0.16,
      colour = spec$border, linewidth = 0.95, lineend = "round"
    ) +
    annotate(
      "text", x = 9.08, y = spec$ymax - 0.27, label = spec$heading,
      hjust = 0, family = JOURNAL_FONT, fontface = "bold",
      size = if (i == 2) 1.55 else 1.62, colour = spec$border
    ) +
    annotate(
      "text", x = 9.08,
      y = spec$ymax - if (i == 3) 0.53 else 0.57,
      label = spec$body,
      hjust = 0, vjust = 1, family = JOURNAL_FONT,
      size = if (i == 3) 1.42 else 1.49, lineheight = 0.95,
      colour = COL[["ink"]]
    )
}

# -----------------------------------------------------------------------------
# Panel b: explicit, non-arbitrary selection path
# -----------------------------------------------------------------------------

selection_flow <- data.frame(
  order = 1:7,
  count = c(33698, 6127, 2496, 1504, 287, 62, 12),
  stage = c(
    "Assayed epithelial features", "Expressed and technically eligible genes",
    "Same-direction genes", "Non-zero bootstrap effects",
    "Threshold-defined 287-gene biological core",
    "Label-blind cross-platform portable genes",
    "Frozen balanced 12-gene signature"
  ),
  criterion_note = c(
    "starting universe",
    "mean expression >0.001 plus technical exclusions",
    "same direction in at least 90% of whole-donor bootstraps",
    "whole-donor bootstrap 95% interval excludes zero",
    "BH FDR at most 0.05 plus agreeing model directions",
    "protein coding and detected on every validation platform",
    "balanced leave-one-donor-out reconstruction plus Kneedle"
  ),
  removed_at_transition = c(NA, 27571, 3631, 992, 1217, 225, 50),
  phase = c(rep("Threshold-defined biological core", 5),
            "Label-blind portability audit", "Discovery-only compact readout"),
  stringsAsFactors = FALSE
)
stopifnot(
  identical(selection_flow$count, c(33698, 6127, 2496, 1504, 287, 62, 12)),
  all(diff(selection_flow$count) < 0)
)

row_y <- c(7.42, 6.40, 5.38, 4.36, 3.34, 1.92, 0.54)
bar_left <- 2.78
bar_width <- 1.55 + 3.20 *
  (log10(selection_flow$count) - log10(min(selection_flow$count))) /
  (log10(max(selection_flow$count)) - log10(min(selection_flow$count)))
bar_right <- bar_left + bar_width
bar_height <- 0.60
node_count <- format(selection_flow$count, big.mark = ",", scientific = FALSE)
node_label <- c(
  "Assayed epithelial\nfeatures",
  "Expressed and technically\neligible genes",
  "Same-direction\ngenes",
  "Non-zero bootstrap\neffects",
  "287-gene core\n89 up · 198 down",
  "62 portable genes\n30 up · 32 down",
  "12-gene signature\n6 up · 6 down"
)
node_fill <- c(rep("#FFFFFF", 4), "#FBECE7", "#EAF2F7", "#FBF3DE")
node_border <- c(rep("#AEB6BD", 4), COL[["route"]], COL[["wnt"]],
                 COL[["adenoma"]])
criteria <- c(
  "Expression threshold >0.001\nand technical filters",
  "Same direction in ≥90%\nof donor bootstraps",
  "95% donor-bootstrap interval\nexcludes zero",
  "BH FDR ≤0.05 and\nconcordant model effects",
  "Protein coding and detected\nacross validation platforms",
  "Donor hold-out fidelity\nand Kneedle"
)
removed <- c("−27,571", "−3,631", "−992", "−1,217", "−225", "−50")

p1b <- ggplot() +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 8.52), clip = "off") +
  theme_void(base_family = JOURNAL_FONT) +
  theme(plot.margin = margin(1.0, 2.3, 0.8, 2.3, "mm")) +
  annotate(
    "rect", xmin = 0.08, xmax = 15.92, ymin = 2.78, ymax = 8.34,
    fill = alpha(COL[["route"]], 0.045), colour = alpha(COL[["route"]], 0.22),
    linewidth = 0.32
  ) +
  annotate(
    "segment", x = 0.14, xend = 0.14, y = 2.91, yend = 8.21,
    colour = COL[["route"]], linewidth = 1.08, lineend = "round"
  ) +
  annotate(
    "rect", xmin = 0.08, xmax = 15.92, ymin = 1.46, ymax = 2.58,
    fill = alpha(COL[["wnt"]], 0.055), colour = alpha(COL[["wnt"]], 0.25),
    linewidth = 0.32
  ) +
  annotate(
    "segment", x = 0.14, xend = 0.14, y = 1.56, yend = 2.48,
    colour = COL[["wnt"]], linewidth = 1.08, lineend = "round"
  ) +
  annotate(
    "rect", xmin = 0.08, xmax = 15.92, ymin = 0.05, ymax = 1.28,
    fill = alpha(COL[["adenoma"]], 0.075), colour = alpha(COL[["adenoma"]], 0.28),
    linewidth = 0.32
  ) +
  annotate(
    "segment", x = 0.14, xend = 0.14, y = 0.15, yend = 1.18,
    colour = COL[["adenoma"]], linewidth = 1.08, lineend = "round"
  ) +
  annotate(
    "text", x = 0.36, y = 8.07, label = "Biological core",
    hjust = 0, family = JOURNAL_FONT, fontface = "bold", size = 1.72,
    colour = COL[["route"]]
  ) +
  annotate(
    "text", x = 15.66, y = 8.07, label = "Excluded",
    hjust = 1, family = JOURNAL_FONT, fontface = "bold", size = 1.42,
    colour = COL[["neutral"]]
  ) +
  annotate(
    "text", x = 0.36, y = 2.34, label = "Cross-platform portability",
    hjust = 0, family = JOURNAL_FONT, fontface = "bold", size = 1.60,
    colour = COL[["wnt"]]
  ) +
  annotate(
    "text", x = 0.36, y = 1.06, label = "Compact signature",
    hjust = 0, family = JOURNAL_FONT, fontface = "bold", size = 1.65,
    colour = COL[["adenoma"]]
  )

for (i in seq_along(row_y)) {
  p1b <- p1b +
    annotation_custom(
      round_box(node_fill[i], node_border[i], 1.65, 0.88),
      xmin = bar_left, xmax = bar_right[i],
      ymin = row_y[i] - bar_height / 2, ymax = row_y[i] + bar_height / 2
    ) +
    annotate(
      "text", x = bar_left + bar_width[i] / 2, y = row_y[i],
      label = node_count[i], family = JOURNAL_FONT, fontface = "bold",
      size = if (i >= 5) 2.00 else 1.91,
      colour = if (i <= 4) COL[["ink"]] else node_border[i]
    ) +
    annotate(
      "text", x = 0.36, y = row_y[i],
      label = node_label[i], hjust = 0, vjust = 0.5,
      family = JOURNAL_FONT, fontface = if (i >= 5) "bold" else "plain",
      size = if (i >= 5) 1.62 else 1.58, lineheight = 0.96,
      colour = COL[["ink"]]
    )
}

for (i in seq_len(length(row_y) - 1)) {
  midpoint_y <- (row_y[i] + row_y[i + 1]) / 2
  transition_colour <- if (i == 4) COL[["route"]] else if (i == 5) {
    COL[["wnt"]]
  } else if (i == 6) {
    COL[["adenoma"]]
  } else {
    COL[["neutral"]]
  }
  p1b <- p1b +
    annotate(
      "segment", x = bar_right[i], xend = bar_right[i],
      y = row_y[i] - bar_height / 2, yend = midpoint_y,
      colour = transition_colour, linewidth = 0.52
    ) +
    annotate(
      "segment", x = bar_right[i], xend = bar_right[i + 1],
      y = midpoint_y, yend = midpoint_y,
      colour = transition_colour, linewidth = 0.52
    ) +
    annotate(
      "segment", x = bar_right[i + 1], xend = bar_right[i + 1],
      y = midpoint_y, yend = row_y[i + 1] + bar_height / 2,
      colour = transition_colour, linewidth = 0.52,
      arrow = arrow(type = "closed", angle = 23, length = unit(1.14, "mm"))
    ) +
    annotate(
      "segment", x = max(bar_right[i], bar_right[i + 1]) + 0.17,
      xend = 8.02, y = midpoint_y, yend = midpoint_y,
      colour = alpha(COL[["neutral"]], 0.48), linewidth = 0.34
    ) +
    annotate(
      "point", x = 8.12, y = midpoint_y, shape = 21, size = 1.08,
      fill = "white", colour = transition_colour, stroke = 0.46
    ) +
    annotate(
      "text", x = 8.34, y = midpoint_y, label = criteria[i],
      hjust = 0, vjust = 0.5, family = JOURNAL_FONT,
      fontface = "bold", size = 1.58, lineheight = 0.96,
      colour = COL[["ink"]]
    ) +
    annotate(
      "text", x = 15.66, y = midpoint_y, label = removed[i],
      hjust = 1, family = JOURNAL_FONT, fontface = "bold", size = 1.43,
      colour = alpha(transition_colour, 0.96)
    )
}

p1b <- p1b +
  annotation_custom(
    lock_icon(COL[["ink"]]), xmin = 4.49, xmax = 4.87,
    ymin = 0.35, ymax = 0.73
  ) +
  annotate(
    "text", x = 5.00, y = 0.54, label = "Frozen",
    hjust = 0, family = JOURNAL_FONT, fontface = "bold", size = 1.42,
    colour = COL[["ink"]]
  )

# -----------------------------------------------------------------------------
# Panels c–f: locked quantitative evidence, re-rendered from source-data tables
# -----------------------------------------------------------------------------

kneedle <- read_source("figure1c_kneedle.tsv")
p1c <- ggplot(kneedle, aes(total_genes, oof_spearman)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = alpha(COL[["route"]], 0.13)) +
  geom_line(colour = COL[["route"]], linewidth = 0.75) +
  geom_point(
    data = filter(kneedle, total_genes == 12), size = 2.35,
    colour = COL[["adenoma"]]
  ) +
  geom_vline(
    xintercept = 12, linetype = "22", colour = COL[["adenoma"]], linewidth = 0.45
  ) +
  annotate(
    "label", x = 17, y = 0.80,
    label = "Kneedle: 12 genes\n6 up + 6 down\nOOF ρ = 0.929",
    size = 1.43, family = JOURNAL_FONT, fill = "white", linewidth = 0.20,
    label.padding = unit(0.65, "mm")
  ) +
  coord_cartesian(xlim = c(2, 60), ylim = c(0.42, 1.01)) +
  labs(
    x = "Balanced signature size (genes)",
    y = "Leave-one-donor-out fidelity to 287-gene core"
  ) +
  theme_journal(base_size = 6.3)

hallmark <- read_source("figure1d_hallmark.tsv")
hallmark <- hallmark %>%
  mutate(
    label = factor(label, levels = label[order(rank_biserial)]),
    q_label = ifelse(q_bh < 0.001,
                     paste0("q=", formatC(q_bh, format = "e", digits = 1)),
                     paste0("q=", formatC(q_bh, format = "f", digits = 3)))
  )
p1d <- ggplot(hallmark, aes(rank_biserial, label, colour = direction)) +
  geom_vline(xintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = rank_biserial, yend = label), linewidth = 0.75) +
  geom_point(size = 2.0) +
  geom_text(
    aes(label = q_label),
    nudge_x = ifelse(hallmark$rank_biserial > 0, 0.055, -0.055),
    hjust = ifelse(hallmark$rank_biserial > 0, 0, 1),
    size = 1.27, family = JOURNAL_FONT, colour = COL[["ink"]]
  ) +
  scale_colour_manual(
    values = c("Adenoma-up" = COL[["route"]], "Adenoma-down" = COL[["wnt"]]),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-0.93, 0.79), clip = "off") +
  labs(
    x = "Rank-biserial enrichment across 6,127 eligible genes",
    y = NULL
  ) +
  theme_journal(base_size = 6.3) +
  theme(axis.text.y = element_text(size = 5.05), plot.margin = margin(1.6, 4, 1.6, 2, "mm"))

chen_scores <- read_source("figure1e_chen_scores.tsv") %>%
  mutate(
    panel = factor(panel, levels = c("287-gene core", "12-gene panel")),
    tissue = factor(tissue, levels = c("Normal", "Adenoma"))
  )
p1e <- ggplot(chen_scores, aes(tissue, panel_score, colour = tissue)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.43, alpha = 0.16) +
  geom_jitter(width = 0.12, size = 0.78, alpha = 0.68) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(Normal = COL[["wnt"]], Adenoma = COL[["route"]])) +
  labs(x = NULL, y = "Within-dataset score") +
  theme_journal(base_size = 6.3) +
  theme(
    legend.position = "none", strip.text = element_text(size = 5.35),
    axis.text.x = element_text(size = 4.80)
  )

chen_pair <- read_source("figure1f_chen_pairs.tsv") %>%
  mutate(tissue = factor(tissue, levels = c("Normal", "Adenoma")))
p1f <- ggplot(chen_pair, aes(tissue, score, group = donor_id)) +
  geom_line(colour = COL[["neutral"]], alpha = 0.40, linewidth = 0.42) +
  geom_point(aes(colour = tissue), size = 1.30) +
  stat_summary(
    aes(group = 1), fun = median, geom = "point", shape = 23,
    fill = "white", size = 2.35, stroke = 0.48, colour = COL[["ink"]]
  ) +
  annotate(
    "text", x = 1.5, y = max(chen_pair$score) + 0.35,
    label = "7/7 increased\npaired P = 0.0156", size = 1.48,
    family = JOURNAL_FONT, lineheight = 0.90
  ) +
  scale_colour_manual(values = c(Normal = COL[["wnt"]], Adenoma = COL[["route"]])) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(x = NULL, y = "12-gene score") +
  theme_journal(base_size = 6.3) +
  theme(legend.position = "none", axis.text.x = element_text(size = 4.80))

row_cd <- clean_panel(p1c) | clean_panel(p1d)
row_cd <- row_cd + plot_layout(widths = c(0.90, 1.10))
row_ef <- clean_panel(p1e) | clean_panel(p1f)
row_ef <- row_ef + plot_layout(widths = c(1.28, 0.72))

fig1 <- clean_panel(p1a) / clean_panel(p1b) / row_cd / row_ef +
  plot_layout(heights = c(1.25, 1.32, 0.94, 0.90)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

export_figure <- function(plot, stem, width_mm, height_mm) {
  paths <- c(
    SVG = file.path(OUT_DIR, paste0(stem, ".svg")),
    PDF = file.path(OUT_DIR, paste0(stem, ".pdf")),
    TIFF = file.path(OUT_DIR, paste0(stem, ".tiff")),
    PNG = file.path(OUT_DIR, paste0(stem, ".png"))
  )
  ggsave(
    paths[["SVG"]], plot = plot, device = svglite::svglite,
    width = width_mm, height = height_mm, units = "mm", bg = "white",
    fix_text_size = FALSE
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

# Confirm that the remaining audited figures are present. This Figure 1-only
# renderer must not overwrite targeted fixes already applied to Figs. 2–6 or
# the supplementary figures.
canonical_stems <- c(
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
)
formats <- c(".svg", ".pdf", ".tiff", ".png")
for (stem in canonical_stems) {
  for (extension in formats) {
    destination <- file.path(OUT_DIR, paste0(stem, extension))
    source <- file.path(SOURCE_PACKAGE, paste0(stem, extension))
    if (!file.exists(destination) && file.exists(source)) {
      if (!file.copy(source, destination, overwrite = FALSE, copy.mode = TRUE)) {
        stop("Could not restore missing audited figure asset: ", destination)
      }
    }
    if (!file.exists(destination)) stop("Missing audited figure asset: ", destination)
  }
}

# Replace only the panel-a and panel-b source records.
write_source(dataset_overview, "figure1a_workflow.tsv")
write_source(resource_landscape, "figure1a_resource_landscape.tsv")
write_source(selection_flow, "figure1b_selection_flow.tsv")

export_row <- export_figure(
  fig1, "figure1_discovery_core_and_objective_reduction", 170, 218
)
write_source(export_row, "figure1_v1_9_export_manifest.tsv")

# Exact-text QA catches regressions back to the abstract y-axis vocabulary.
svg_text <- paste(readLines(export_row$file[export_row$format == "SVG"], warn = FALSE),
                  collapse = "\n")
qa <- data.frame(
  check = c(
    "resource_total_is_22", "six_resource_types_are_explicit",
    "six_resource_icons_are_defined",
    "selection_counts_are_monotone",
    "no_abstract_assayed_features_axis", "threshold_language_present",
    "key_validation_counts_present", "frozen_boundary_present",
    "internal_audit_copy_removed", "typography_spacing_normalised",
    "all_four_formats_exported"
  ),
  pass = c(
    sum(dataset_overview$public_resources) == 22 &&
      grepl("Public resources (n = 22)", svg_text, fixed = TRUE),
    all(vapply(resource_landscape$resource_type, grepl, logical(1),
               x = svg_text, fixed = TRUE)),
    length(unique(resource_landscape$icon_key)) == 6 &&
      all(resource_landscape$icon_key %in% c(
        "single_cell", "bulk_ffpe", "multiome", "perturbation", "spatial",
        "proteomics"
      )),
    all(diff(selection_flow$count) < 0),
    !grepl("Genes retained \\(log scale\\)", svg_text),
    grepl("BH FDR", svg_text) && grepl("90%", svg_text),
    grepl("44 specimens · 24 donors", svg_text) &&
      grepl("5 cohorts · 203 samples · 161 patient clusters", svg_text) &&
      grepl("51 FFPE pairs \\(102 specimens\\)", svg_text) &&
      grepl("40 samples · 12 patients", svg_text) &&
      grepl("6 datasets · 8 contrasts", svg_text) &&
      grepl("CRC Atlas", svg_text) &&
      grepl("6 tissue sections", svg_text) &&
      grepl("4 datasets", svg_text),
    grepl("Frozen", svg_text),
    !grepl(
      paste(c(
        "DISCOVERY DATA ONLY", "GENE SELECTION PERMITTED",
        "NO GENE SELECTION", "PORTABILITY AUDIT", "labels hidden",
        "criteria in b", "fidelity in c", "FROZEN BEFORE ALL TESTS",
        "no preset gene count"
      ), collapse = "|"),
      svg_text
    ),
    grepl("Public resources (n = 22)", svg_text, fixed = TRUE) &&
      grepl("Analysis design", svg_text) &&
      !grepl("PUBLIC RESOURCES  ·", svg_text, fixed = TRUE) &&
      !grepl("textLength=", svg_text, fixed = TRUE),
    all(file.exists(export_row$file))
  ),
  stringsAsFactors = FALSE
)
write_source(qa, "figure1_v1_9_qa.tsv")
if (!all(qa$pass)) {
  stop("Figure 1 v1.9 QA failed: ", paste(qa$check[!qa$pass], collapse = ", "))
}

writeLines(capture.output(sessionInfo()), file.path(SOURCE_DIR, "figure1_v1_9_sessionInfo.txt"))
message("Figure 1 v1.9 exported to: ", OUT_DIR)
