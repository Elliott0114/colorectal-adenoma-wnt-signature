#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(grid)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
out_dir <- file.path(root, "figures", "communications_biology_v3.0")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

font <- "Arial"
colour <- c(
  ink = "#24292E", muted = "#6B747C", line = "#C7CDD2",
  normal = "#3274A1", adenoma = "#D45D3F", down = "#3274A1",
  up = "#D45D3F", gold = "#D6A23A", purple = "#7669A8",
  green = "#3F8A70", pale_blue = "#EDF4F8", pale_orange = "#FBF0EC",
  pale_gold = "#FBF6E9", pale_purple = "#F2F0F7", pale = "#F5F6F7"
)

read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, check.names = FALSE, quote = "", comment.char = "")
}

write_source <- function(x, name) {
  write.table(
    x, file.path(source_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )
}

theme_journal <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = font) +
    theme(
      text = element_text(family = font, colour = colour[["ink"]]),
      axis.line = element_line(linewidth = 0.35, colour = colour[["ink"]]),
      axis.ticks = element_line(linewidth = 0.3, colour = colour[["ink"]]),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.6),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(2, 2, 2, 2, "mm"),
      plot.tag = element_text(size = 9, face = "bold")
    )
}

export_figure <- function(plot, stem, width_mm = 178, height_mm = 225) {
  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf, bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(stem, ".svg")), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = function(...) svglite::svglite(..., fix_text_size = FALSE),
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(stem, ".tiff")), plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 600,
    device = ragg::agg_tiff, compression = "lzw", background = "white"
  )
  ggsave(
    file.path(out_dir, paste0(stem, ".png")), plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 300,
    device = ragg::agg_png, background = "white"
  )
}

common <- read_tsv(file.path(
  result_root, "common_effects", "cross_state_common_effects.tsv.gz"
))
leaveout <- read_tsv(file.path(
  result_root, "donor_leaveout_stability", "donor_leaveout_rank_stability.tsv"
))
# a. Data landscape and one discovery-freeze-validation path.
resources <- data.frame(
  order = 1:6,
  abbreviation = c("sc", "RNA", "ATAC", "PERT", "ST", "MS"),
  type = c(
    "sc/snRNA-seq", "Bulk and FFPE\nRNA", "Matched\nRNA–ATAC",
    "Genetic and drug\nperturbation", "Spatial\nRNA", "Proteomics"
  ),
  scope = c(
    "Chen · Becker\nCRC Atlas", "5 cohorts\n51 FFPE pairs",
    "40 samples\n12 patients", "6 datasets",
    "6 sections", "4 datasets"
  ),
  x = 1:6,
  stringsAsFactors = FALSE
)

p1a <- ggplot() +
  annotate(
    "text", x = 0.3, y = 4.88,
    label = "22 public resources across six data modalities",
    hjust = 0, family = font, fontface = "bold", size = 3.0,
    colour = colour[["ink"]]
  ) +
  geom_rect(
    data = resources,
    aes(xmin = x - 0.43, xmax = x + 0.43, ymin = 2.92, ymax = 4.48),
    fill = "white", colour = colour[["line"]], linewidth = 0.45
  ) +
  geom_point(
    data = resources, aes(x = x, y = 4.10), shape = 21, size = 5.7,
    stroke = 0.55, fill = colour[["pale_blue"]], colour = colour[["normal"]]
  ) +
  geom_text(
    data = resources, aes(x = x, y = 4.10, label = abbreviation),
    family = font, fontface = "bold", size = 2.35,
    colour = colour[["normal"]]
  ) +
  geom_text(
    data = resources, aes(x = x, y = 3.60, label = type),
    family = font, fontface = "bold", size = 1.75, lineheight = 0.92,
    colour = colour[["ink"]]
  ) +
  geom_text(
    data = resources, aes(x = x, y = 3.15, label = scope),
    family = font, size = 1.55, lineheight = 0.92,
    colour = colour[["muted"]]
  ) +
  annotate(
    "rect", xmin = 0.55, xmax = 2.05, ymin = 0.55, ymax = 2.40,
    fill = colour[["pale_orange"]], colour = colour[["adenoma"]],
    linewidth = 0.55
  ) +
  annotate(
    "text", x = 0.73, y = 2.13, label = "Discovery",
    hjust = 0, family = font, fontface = "bold", size = 2.6,
    colour = colour[["adenoma"]]
  ) +
  annotate(
    "text", x = 0.73, y = 1.66,
    label = "Chen epithelial scRNA-seq\n27 donors\nABS · GOB · TAC",
    hjust = 0, family = font, size = 1.95, lineheight = 0.98
  ) +
  annotate(
    "rect", xmin = 2.55, xmax = 3.75, ymin = 0.55, ymax = 2.40,
    fill = colour[["pale_gold"]], colour = colour[["gold"]], linewidth = 0.55
  ) +
  annotate(
    "text", x = 2.72, y = 2.13, label = "Freeze",
    hjust = 0, family = font, fontface = "bold", size = 2.6,
    colour = colour[["gold"]]
  ) +
  annotate(
    "text", x = 2.72, y = 1.66,
    label = "1,843 genes\nFixed directions\nand score",
    hjust = 0, family = font, size = 1.95, lineheight = 0.98
  ) +
  annotate(
    "rect", xmin = 4.25, xmax = 6.45, ymin = 0.55, ymax = 2.40,
    fill = colour[["pale_blue"]], colour = colour[["normal"]], linewidth = 0.55
  ) +
  annotate(
    "text", x = 4.43, y = 2.13, label = "Validation",
    hjust = 0, family = font, fontface = "bold", size = 2.6,
    colour = colour[["normal"]]
  ) +
  annotate(
    "text", x = 4.43, y = 1.61,
    label = paste(
      "23 donor-disjoint subjects",
      "5 transcriptomic cohorts",
      "51 FFPE pairs",
      "Regulatory and perturbation support",
      sep = "\n"
    ),
    hjust = 0, family = font, size = 1.85, lineheight = 0.98
  ) +
  annotate(
    "segment", x = 2.08, xend = 2.47, y = 1.47, yend = 1.47,
    linewidth = 0.6, colour = colour[["muted"]],
    arrow = arrow(type = "closed", length = unit(1.6, "mm"))
  ) +
  annotate(
    "segment", x = 3.78, xend = 4.17, y = 1.47, yend = 1.47,
    linewidth = 0.6, colour = colour[["muted"]],
    arrow = arrow(type = "closed", length = unit(1.6, "mm"))
  ) +
  coord_cartesian(xlim = c(0.25, 6.75), ylim = c(0.40, 5.02), clip = "off") +
  theme_void(base_family = font) +
  theme(plot.margin = margin(1, 2, 1, 2, "mm"))

# b. Objective reduction from the assayed feature space.
selection <- data.frame(
  step = factor(
    c("Assayed", "Jointly testable", "State-shared"),
    levels = c("Assayed", "Jointly testable", "State-shared")
  ),
  count = c(33703, 8221, 1843),
  rule = c(
    "features in discovery object",
    "expressed in ABS, GOB and TAC",
    "FDR ≤ 0.05\nsame sign in ABS/GOB/TAC\nlfsr ≤ 0.05"
  ),
  width = c(3.0, 2.25, 1.55),
  y = c(3.2, 2.1, 1.0),
  half_height = c(0.40, 0.40, 0.56),
  count_offset = c(0.10, 0.10, 0.20),
  rule_offset = c(-0.17, -0.17, -0.28),
  rule_size = c(1.65, 1.65, 1.42),
  fill = c(colour[["pale"]], colour[["pale_blue"]], colour[["pale_orange"]])
)

p1b <- ggplot() +
  geom_rect(
    data = selection,
    aes(xmin = -width / 2, xmax = width / 2, ymin = y - half_height, ymax = y + half_height),
    fill = selection$fill, colour = c(
      colour[["line"]], colour[["normal"]], colour[["adenoma"]]
    ), linewidth = 0.5
  ) +
  geom_text(
    data = selection,
    aes(x = 0, y = y + count_offset, label = comma(count)),
    family = font, fontface = "bold", size = 3.0
  ) +
  geom_text(
    data = selection, aes(x = 0, y = y + rule_offset, label = rule, size = rule_size),
    family = font, lineheight = 0.82,
    colour = colour[["muted"]]
  ) +
  scale_size_identity() +
  annotate(
    "text", x = 0, y = 0.20, label = "884 adenoma-up · 959 adenoma-down",
    family = font, fontface = "bold", size = 2.2,
    colour = colour[["adenoma"]]
  ) +
  annotate(
    "segment", x = 0, xend = 0, y = 2.78, yend = 2.52,
    colour = colour[["muted"]], linewidth = 0.45,
    arrow = arrow(type = "closed", length = unit(1.25, "mm"))
  ) +
  annotate(
    "segment", x = 0, xend = 0, y = 1.68, yend = 1.58,
    colour = colour[["muted"]], linewidth = 0.45,
    arrow = arrow(type = "closed", length = unit(1.25, "mm"))
  ) +
  coord_cartesian(xlim = c(-1.7, 1.7), ylim = c(0.02, 3.68), clip = "off") +
  theme_void(base_family = font) +
  theme(plot.margin = margin(2, 3, 2, 3, "mm"))

# c. Complete common-effect ranking.
ranking <- common %>%
  arrange(common_z) %>%
  mutate(
    rank = row_number(),
    class = case_when(
      strict_state_shared & shared_direction == "up" ~ "Adenoma-up",
      strict_state_shared & shared_direction == "down" ~ "Adenoma-down",
      TRUE ~ "Other testable genes"
    )
  )

p1c <- ggplot(ranking, aes(rank, common_z, colour = class)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = colour[["line"]]) +
  geom_point(size = 0.65, alpha = 0.68) +
  scale_colour_manual(values = c(
    "Adenoma-down" = colour[["down"]],
    "Other testable genes" = "#C9CED2",
    "Adenoma-up" = colour[["up"]]
  )) +
  scale_x_continuous(labels = comma, breaks = c(1, 4000, 8221)) +
  annotate("text", x = 500, y = -8.7, label = "adenoma-down", hjust = 0,
           family = font, size = 1.8, colour = colour[["down"]]) +
  annotate("text", x = 7900, y = 11.2, label = "adenoma-up", hjust = 1,
           family = font, size = 1.8, colour = colour[["up"]]) +
  guides(colour = "none") +
  labs(x = "Common-effect gene rank", y = "Signed common effect (z)") +
  theme_journal() +
  theme(plot.margin = margin(3, 2, 2, 2, "mm"))

# d. Whole-donor stability of the complete common ranking.
leaveout_common <- leaveout %>%
  filter(scope == "common_GLS") %>%
  arrange(effect_spearman) %>%
  mutate(index = row_number())

p1f <- ggplot(leaveout_common, aes(index, effect_spearman)) +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = 0.97, ymax = 1,
    fill = colour[["pale_blue"]], colour = NA
  ) +
  geom_hline(yintercept = 0.99, linewidth = 0.35, linetype = 2,
             colour = colour[["normal"]]) +
  geom_point(size = 1.8, colour = colour[["normal"]]) +
  scale_y_continuous(limits = c(0.968, 1.0005), breaks = c(0.97, 0.98, 0.99, 1.00)) +
  scale_x_continuous(breaks = c(1, 9, 18, 27)) +
  labs(x = "Discovery donor left out", y = "Rank correlation with complete fit") +
  theme_journal()

figure <- p1a /
  (p1b | p1c) /
  p1f +
  plot_layout(heights = c(1.30, 1, 0.70)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(family = font, size = 9, face = "bold"),
    plot.tag.location = "panel",
    plot.tag.position = c(0.026, 0.988)
  )

write_source(resources, "figure1a_resources.tsv")
write_source(selection, "figure1b_selection.tsv")
write_source(ranking, "figure1c_ranking.tsv")
write_source(leaveout_common, "figure1d_leaveout.tsv")

export_figure(
  figure,
  "figure1_study_design_and_programme_derivation",
  width_mm = 178,
  height_mm = 170
)

cat("Figure 1 exported to ", out_dir, "\n", sep = "")
