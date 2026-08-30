#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")

root <- normalizePath(".", mustWork = TRUE)
external_root <- file.path(root, "results", "state_aware_program_v1", "external_validation")
extended_root <- file.path(root, "results", "state_aware_program_v1", "extended_validation", "perturbation_spatial")
out_dir <- file.path(root, "figures", "communications_biology_v2.0")
source_dir <- file.path(out_dir, "source_data")

sample_scores <- read_tsv(file.path(external_root, "perturbation_sample_scores.tsv"))
apc_effects <- read_tsv(file.path(external_root, "apc_organoid_effects.tsv"))
tcf_clones <- read_tsv(file.path(external_root, "tcf7l2_clone_effects.tsv")) %>%
  mutate(direction_matches_expected = tolower(as.character(direction_matches_expected)) == "true")
genetic_summary <- read_tsv(file.path(extended_root, "perturbation_effect_summary.tsv"))
genetic_units <- read_tsv(file.path(extended_root, "perturbation_unit_effects.tsv"))

arrow_closed <- grid::arrow(type = "closed", length = grid::unit(1.7, "mm"))
ellipse_data <- function(cx, cy, rx, ry, n = 180L) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + rx * cos(theta), y = cy + ry * sin(theta))
}
nucleus <- ellipse_data(8.95, 2.72, 0.78, 0.93)

# a. Canonical WNT pathway and the experimentally interrogated nodes.
p5a <- ggplot() +
  annotate("rect", xmin = 0.15, xmax = 13.85, ymin = 1.40, ymax = 4.15,
           fill = figure_colours[["pale_blue"]], colour = NA) +
  annotate("rect", xmin = 3.10, xmax = 7.15, ymin = 0.22, ymax = 1.15,
           fill = figure_colours[["pale"]], colour = NA) +
  annotate("text", x = 0.35, y = 3.92, label = "WNT signalling or APC loss",
           hjust = 0, family = figure_font, fontface = "bold", size = 2.55,
           colour = figure_colours[["normal"]]) +
  annotate("point", x = 0.95, y = 2.72, shape = 21, size = 5.0, stroke = 0.7,
           fill = "white", colour = figure_colours[["normal"]]) +
  annotate("text", x = 0.95, y = 2.72, label = "WNT", family = figure_font,
           fontface = "bold", size = 2.2, colour = figure_colours[["normal"]]) +
  annotate("segment", x = 1.35, xend = 2.25, y = 2.72, yend = 2.72,
           colour = figure_colours[["normal"]], linewidth = 0.75, arrow = arrow_closed) +
  annotate("rect", xmin = 2.35, xmax = 2.52, ymin = 2.20, ymax = 3.24,
           fill = "white", colour = figure_colours[["normal"]], linewidth = 0.6) +
  annotate("rect", xmin = 2.63, xmax = 2.80, ymin = 2.20, ymax = 3.24,
           fill = "white", colour = figure_colours[["normal"]], linewidth = 0.6) +
  annotate("text", x = 2.58, y = 1.92, label = "FZD–LRP5/6",
           family = figure_font, fontface = "bold", size = 2.05) +
  annotate("segment", x = 2.92, xend = 4.05, y = 2.72, yend = 2.72,
           colour = figure_colours[["normal"]], linewidth = 0.75, arrow = arrow_closed) +
  annotate("label", x = 4.75, y = 2.72, label = "APC · AXIN\nGSK3β",
           family = figure_font, size = 2.0, lineheight = 0.92,
           fill = "white", colour = figure_colours[["grey"]], linewidth = 0.3,
           label.padding = grid::unit(0.8, "mm")) +
  annotate("segment", x = 5.45, xend = 6.15, y = 2.72, yend = 2.72,
           colour = figure_colours[["normal"]], linewidth = 0.75, arrow = arrow_closed) +
  annotate("point", x = c(6.38, 6.66, 6.94), y = c(2.58, 2.84, 2.58),
           shape = 21, size = c(3.0, 3.6, 3.0), stroke = 0.55,
           fill = "white", colour = figure_colours[["normal"]]) +
  annotate("text", x = 6.66, y = 3.35, label = "stabilised β-catenin",
           family = figure_font, fontface = "bold", size = 2.2,
           colour = figure_colours[["normal"]]) +
  annotate("segment", x = 7.28, xend = 8.05, y = 2.72, yend = 2.72,
           colour = figure_colours[["normal"]], linewidth = 0.75, arrow = arrow_closed) +
  geom_polygon(data = nucleus, aes(x, y), inherit.aes = FALSE,
               fill = "white", colour = figure_colours[["purple"]], linewidth = 0.55) +
  annotate("text", x = 8.95, y = 2.78, label = "β-catenin\nTCF/LEF · ASCL2",
           family = figure_font, fontface = "bold", size = 2.0,
           colour = figure_colours[["purple"]], lineheight = 0.95) +
  annotate("text", x = 8.95, y = 3.82, label = "Nucleus", family = figure_font,
           size = 1.75, colour = figure_colours[["muted"]]) +
  annotate("segment", x = 9.78, xend = 10.55, y = 2.72, yend = 2.72,
           colour = figure_colours[["purple"]], linewidth = 0.75, arrow = arrow_closed) +
  annotate("text", x = 10.72, y = 3.12, label = "WNT/stem regulation ↑",
           hjust = 0, family = figure_font, fontface = "bold", size = 2.05,
           colour = figure_colours[["adenoma"]]) +
  annotate("text", x = 10.72, y = 2.35, label = "Mature epithelial functions ↓",
           hjust = 0, family = figure_font, fontface = "bold", size = 2.05,
           colour = figure_colours[["normal"]]) +
  annotate("segment", x = 12.60, xend = 13.17, y = 2.72, yend = 2.72,
           colour = figure_colours[["ink"]], linewidth = 0.7, arrow = arrow_closed) +
  annotate("text", x = 13.25, y = 2.72, label = "Epithelial\nidentity programme",
           hjust = 0, family = figure_font, fontface = "bold", size = 2.0,
           lineheight = 0.95, colour = figure_colours[["ink"]]) +
  annotate("label", x = 4.60, y = 0.75, label = "APC · AXIN · GSK3β active",
           family = figure_font, size = 1.8, fill = "white",
           colour = figure_colours[["grey"]], linewidth = 0.25) +
  annotate("segment", x = 5.78, xend = 6.48, y = 0.75, yend = 0.75,
           colour = figure_colours[["grey"]], linewidth = 0.58, arrow = arrow_closed) +
  annotate("text", x = 6.65, y = 0.75, label = "β-catenin degradation",
           hjust = 0, family = figure_font, size = 1.8,
           colour = figure_colours[["grey"]]) +
  annotate("segment", x = 0.55, xend = 2.15, y = 0.08, yend = 0.08,
           colour = figure_colours[["normal"]], linewidth = 1.0, lineend = "round") +
  annotate("text", x = 1.35, y = -0.12, label = "WNT/RSPO withdrawal",
           family = figure_font, fontface = "bold", size = 1.8,
           colour = figure_colours[["normal"]]) +
  annotate("segment", x = 3.85, xend = 5.25, y = 0.08, yend = 0.08,
           colour = figure_colours[["adenoma"]], linewidth = 1.0, lineend = "round") +
  annotate("text", x = 4.55, y = -0.12, label = "APC loss / restoration",
           family = figure_font, fontface = "bold", size = 1.8,
           colour = figure_colours[["adenoma"]]) +
  annotate("segment", x = 8.15, xend = 9.75, y = 0.08, yend = 0.08,
           colour = figure_colours[["purple"]], linewidth = 1.0, lineend = "round") +
  annotate("text", x = 8.95, y = -0.12, label = "TCF7L2 / ASCL2 editing",
           family = figure_font, fontface = "bold", size = 1.8,
           colour = figure_colours[["purple"]]) +
  coord_cartesian(xlim = c(0, 14.2), ylim = c(-0.28, 4.3), clip = "off") +
  theme_void(base_family = figure_font) +
  theme(plot.margin = margin(2, 2, 1, 2, "mm"))

# b. Donor-level APC-by-WNT experiment for the full programme.
apc_samples <- sample_scores %>%
  filter(dataset == "GSE125472", signature_id == "state_shared_1843") %>%
  mutate(
    condition = factor(
      paste(genotype, wnt_rspo),
      levels = c("WT with", "APC with", "WT without", "APC without"),
      labels = c("WT\nWNT+", "APC-KO\nWNT+", "WT\nWNT−", "APC-KO\nWNT−")
    )
  )

p5b <- ggplot(apc_samples, aes(condition, route_score, group = donor_id, colour = donor_id)) +
  geom_line(linewidth = 0.6, alpha = 0.72) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c(Donor1 = figure_colours[["adenoma"]], Donor2 = figure_colours[["normal"]], Donor3 = figure_colours[["purple"]])) +
  labs(x = NULL, y = "Full-programme score") +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

# c. Prespecified APC/WNT contrasts for full and reduced readouts.
apc_plot <- apc_effects %>%
  filter(comparison %in% c("WT_withdrawal", "APC_vs_WT_without_Wnt", "genotype_by_Wnt_interaction")) %>%
  mutate(
    comparison_label = factor(
      recode(
        comparison,
        WT_withdrawal = "WNT withdrawal in WT",
        APC_vs_WT_without_Wnt = "APC-KO vs WT, WNT−",
        genotype_by_Wnt_interaction = "Genotype × WNT interaction"
      ),
      levels = rev(c("WNT withdrawal in WT", "APC-KO vs WT, WNT−", "Genotype × WNT interaction"))
    ),
    readout = factor(signature_id, levels = c("state_shared_1843", "compact_8"), labels = c("Full programme", "Eight-gene candidate"))
  )

p5c <- ggplot(apc_plot, aes(mean_difference, comparison_label, colour = readout, shape = readout)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = min_difference, xmax = max_difference), height = 0, linewidth = 0.62, position = position_dodge(width = 0.36)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.36)) +
  scale_colour_manual(values = c("Full programme" = figure_colours[["grey"]], "Eight-gene candidate" = figure_colours[["adenoma"]])) +
  scale_shape_manual(values = c("Full programme" = 16, "Eight-gene candidate" = 18)) +
  labs(x = "Mean donor difference (range; n = 3)", y = NULL) +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

# d. Independent genetic interventions using the reduced readout.
independent_labels <- c(
  apc_restoration_shApc = "Apc restoration",
  apc_restoration_shApc_Kras = "Apc restoration + Kras",
  ascl2_ko_vs_resting_wt = "Ascl2 knockout",
  conditional_wnt_silencing = "Conditional WNT silencing"
)
independent <- genetic_summary %>%
  filter(feature == "route_score", comparison %in% names(independent_labels)) %>%
  mutate(
    label = factor(independent_labels[comparison], levels = rev(unname(independent_labels)))
  )
independent_units <- genetic_units %>%
  filter(feature == "route_score", comparison %in% names(independent_labels)) %>%
  mutate(label = factor(independent_labels[comparison], levels = levels(independent$label)))

p5d <- ggplot(independent, aes(mean_difference, label)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = bootstrap_mean_ci_low, xmax = bootstrap_mean_ci_high), height = 0, linewidth = 0.6, colour = figure_colours[["normal"]]) +
  geom_point(data = independent_units, aes(difference, label), inherit.aes = FALSE,
             size = 1.3, alpha = 0.55, colour = figure_colours[["normal"]],
             position = position_jitter(height = 0.08, width = 0)) +
  geom_point(size = 2.25, colour = figure_colours[["normal"]]) +
  labs(x = "Change in eight-gene score", y = NULL) +
  theme_cb()

# e. TCF7L2-edited clones, with both readout scales shown for each clone.
tcf_plot <- tcf_clones %>%
  mutate(
    readout = factor(signature_id, levels = c("state_shared_1843", "compact_8"), labels = c("Full programme", "Eight-gene candidate")),
    clone = paste(cell_line, clone_id, sep = " · ")
  )

p5e <- ggplot(tcf_plot, aes(readout, difference_vs_WT, group = clone, colour = cell_line)) +
  geom_hline(yintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_line(linewidth = 0.45, alpha = 0.55) +
  geom_point(aes(shape = genotype), size = 2.0) +
  scale_colour_manual(values = c(HCT116 = figure_colours[["purple"]], HT29 = figure_colours[["green"]])) +
  scale_shape_manual(values = c(KO = 16, Het = 17)) +
  labs(x = NULL, y = "TCF7L2-edited clone − WT") +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

middle <- p5b | p5c
bottom <- p5d | p5e
figure <- p5a / middle / bottom +
  plot_layout(heights = c(0.72, 1, 0.92)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

pathway_source <- data.frame(
  node = c("WNT/FZD-LRP", "APC-AXIN-GSK3B", "beta-catenin", "TCF_LEF_ASCL2", "identity_programme"),
  role = c("ligand_receptor", "destruction_complex", "signal_effector", "nuclear_regulation", "measured_outcome")
)
write_source(pathway_source, source_dir, "figure5a_pathway_nodes.tsv")
write_source(apc_samples, source_dir, "figure5b_apc_wnt_donor_scores.tsv")
write_source(apc_plot, source_dir, "figure5c_apc_wnt_contrasts.tsv")
write_source(independent_units, source_dir, "figure5d_independent_genetic_units.tsv")
write_source(tcf_plot, source_dir, "figure5e_tcf7l2_clones.tsv")

export_cb_figure(
  figure, out_dir,
  "figure5_genetic_perturbation_support",
  width_mm = 178, height_mm = 205
)

cat("Figure 5 exported to ", out_dir, "\n", sep = "")
