#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")
suppressPackageStartupMessages({
  library(ggrepel)
})

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
out_dir <- file.path(root, "figures", "communications_biology_v5.0")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

external_root <- file.path(result_root, "external_validation")
two_axis_root <- file.path(result_root, "identity_reversal_target_prioritization_v1")
module_root <- file.path(
  result_root, "functional_architecture_exploratory_v2_1"
)

sample_scores <- read_tsv(file.path(external_root, "perturbation_sample_scores.tsv"))
apc_effects <- read_tsv(file.path(external_root, "apc_organoid_effects.tsv"))
two_axis <- read_tsv(file.path(two_axis_root, "perturbation_two_component_summary.tsv"))
module_effects <- read_tsv(file.path(
  module_root, "module_perturbation_protein", "module_perturbation_summary.tsv"
))
module_correlations <- read_tsv(file.path(
  module_root, "integration_summary", "module_perturbation_profile_correlations.tsv"
))

arrow_closed <- grid::arrow(type = "closed", length = grid::unit(1.45, "mm"))

# a. Minimal WNT-pathway context leading to two separately measured responses.
p6a <- ggplot() +
  annotate(
    "rect", xmin = 0.20, xmax = 10.15, ymin = 0.38, ymax = 2.72,
    fill = figure_colours[["pale_blue"]], colour = NA
  ) +
  annotate(
    "text", x = 0.40, y = 2.47, label = "APC–WNT regulatory context",
    hjust = 0, family = figure_font, fontface = "bold", size = 2.45,
    colour = figure_colours[["normal"]]
  ) +
  annotate(
    "point", x = 0.85, y = 1.55, shape = 21, size = 4.9,
    stroke = 0.6, fill = "white", colour = figure_colours[["normal"]]
  ) +
  annotate(
    "text", x = 0.85, y = 1.55, label = "WNT", family = figure_font,
    fontface = "bold", size = 2.0, colour = figure_colours[["normal"]]
  ) +
  annotate(
    "segment", x = 1.25, xend = 1.92, y = 1.55, yend = 1.55,
    linewidth = 0.62, colour = figure_colours[["normal"]], arrow = arrow_closed
  ) +
  annotate(
    "label", x = 2.55, y = 1.55, label = "FZD–LRP5/6",
    family = figure_font, size = 1.85, fill = "white",
    colour = figure_colours[["normal"]], linewidth = 0.25
  ) +
  annotate(
    "segment", x = 3.18, xend = 3.78, y = 1.55, yend = 1.55,
    linewidth = 0.62, colour = figure_colours[["normal"]], arrow = arrow_closed
  ) +
  annotate(
    "label", x = 4.55, y = 1.55,
    label = "APC loss / destruction\ncomplex inhibition",
    family = figure_font, size = 1.68, lineheight = 0.90, fill = "white",
    colour = figure_colours[["grey"]], linewidth = 0.25
  ) +
  annotate(
    "segment", x = 5.35, xend = 5.96, y = 1.55, yend = 1.55,
    linewidth = 0.62, colour = figure_colours[["purple"]], arrow = arrow_closed
  ) +
  annotate(
    "label", x = 6.75, y = 1.55, label = "β-catenin–TCF/LEF\nASCL2",
    family = figure_font, size = 1.80, lineheight = 0.9, fill = "white",
    colour = figure_colours[["purple"]], linewidth = 0.25
  ) +
  annotate(
    "segment", x = 7.55, xend = 8.14, y = 1.55, yend = 1.55,
    linewidth = 0.62, colour = figure_colours[["purple"]], arrow = arrow_closed
  ) +
  annotate(
    "text", x = 8.30, y = 1.93, label = "WNT/stem response ↑",
    hjust = 0, family = figure_font, fontface = "bold", size = 1.95,
    colour = figure_colours[["adenoma"]]
  ) +
  annotate(
    "text", x = 8.30, y = 1.18, label = "Mature epithelial functions ↓",
    hjust = 0, family = figure_font, fontface = "bold", size = 1.95,
    colour = figure_colours[["normal"]]
  ) +
  annotate(
    "segment", x = 10.35, xend = 11.05, y = 1.55, yend = 1.55,
    linewidth = 0.58, colour = figure_colours[["ink"]], arrow = arrow_closed
  ) +
  annotate(
    "rect", xmin = 11.15, xmax = 14.70, ymin = 0.38, ymax = 2.72,
    fill = "white", colour = figure_colours[["line"]], linewidth = 0.45
  ) +
  annotate(
    "text", x = 11.40, y = 2.43, label = "Measure independently",
    hjust = 0, family = figure_font, fontface = "bold", size = 2.25
  ) +
  annotate(
    "text", x = 11.40, y = 1.86, label = "1  WNT/stem suppression →",
    hjust = 0, family = figure_font, size = 1.90,
    colour = figure_colours[["adenoma"]]
  ) +
  annotate(
    "text", x = 11.40, y = 1.18, label = "2  Mature-function restoration →",
    hjust = 0, family = figure_font, size = 1.90,
    colour = figure_colours[["normal"]]
  ) +
  coord_cartesian(xlim = c(0, 14.9), ylim = c(0.22, 2.90), clip = "off") +
  theme_void(base_family = figure_font) +
  theme(plot.margin = margin(2, 2, 1, 2, "mm"))

# b. Donor-matched APC-by-WNT organoid scores.
apc_samples <- sample_scores %>%
  filter(dataset == "GSE125472", signature_id == "state_shared_1843") %>%
  mutate(
    condition = factor(
      paste(genotype, wnt_rspo),
      levels = c("WT with", "APC with", "WT without", "APC without"),
      labels = c("WT\nWNT+", "APC-KO\nWNT+", "WT\nWNT−", "APC-KO\nWNT−")
    )
  )

p6b <- ggplot(
  apc_samples,
  aes(condition, route_score, group = donor_id, colour = donor_id)
) +
  geom_line(linewidth = 0.62, alpha = 0.75) +
  geom_point(size = 2.05) +
  scale_colour_manual(values = c(
    Donor1 = figure_colours[["adenoma"]],
    Donor2 = figure_colours[["normal"]],
    Donor3 = figure_colours[["purple"]]
  )) +
  labs(x = NULL, y = "State-shared remodelling score") +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

# c. Prespecified APC/WNT contrasts.
apc_plot <- apc_effects %>%
  filter(
    signature_id == "state_shared_1843",
    comparison %in% c(
      "WT_withdrawal", "APC_vs_WT_without_Wnt", "genotype_by_Wnt_interaction"
    )
  ) %>%
  mutate(
    comparison_label = factor(
      recode(
        comparison,
        WT_withdrawal = "WNT withdrawal in WT",
        APC_vs_WT_without_Wnt = "APC-KO vs WT, WNT−",
        genotype_by_Wnt_interaction = "Genotype × WNT interaction"
      ),
      levels = rev(c(
        "WNT withdrawal in WT", "APC-KO vs WT, WNT−",
        "Genotype × WNT interaction"
      ))
    )
  )

p6c <- ggplot(apc_plot, aes(mean_difference, comparison_label)) +
  geom_vline(
    xintercept = 0, colour = figure_colours[["muted"]],
    linetype = 2, linewidth = 0.35
  ) +
  geom_errorbarh(
    aes(xmin = min_difference, xmax = max_difference),
    height = 0, linewidth = 0.65, colour = figure_colours[["adenoma"]]
  ) +
  geom_point(size = 2.25, colour = figure_colours[["adenoma"]]) +
  labs(x = "Mean donor difference (range; n = 3)", y = NULL) +
  theme_cb()

# d. Two-coordinate perturbation map.
label_map <- c(
  apc_restoration_shApc_Kras = "APC restoration\n+ Kras",
  wnt_rspo_withdrawal_in_WT = "WNT withdrawal\nWT organoids",
  tcf7l2_KO_vs_WT = "TCF7L2 KO",
  conditional_wnt_silencing = "Conditional WNT\nsilencing",
  wnt_rspo_withdrawal_in_APC_KO = "WNT withdrawal\nAPC-KO",
  apc_restoration_shApc = "APC restoration",
  pri724_reversal_of_trametinib = "PRI-724 after\ntrametinib",
  doxycycline_control_shRenilla = "Doxycycline control",
  tcf7l2_heterozygous_vs_WT = "TCF7L2 heterozygous",
  trametinib_vs_dmso = "Trametinib",
  ascl2_ko_vs_resting_wt = "ASCL2 KO"
)
two_axis <- two_axis %>%
  mutate(label = unname(label_map[comparison]))
quadrant_colours <- c(
  "both axes favourable" = figure_colours[["green"]],
  "WNT/stem only" = figure_colours[["adenoma"]],
  "mature-function only" = figure_colours[["gold"]],
  "neither axis favourable" = figure_colours[["grey"]]
)

p6d <- ggplot(
  two_axis,
  aes(mean_wnt_stem_suppression, mean_mature_function_restoration, colour = quadrant)
) +
  annotate(
    "rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
    fill = "#ECF7F2", colour = NA
  ) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = figure_colours[["muted"]]) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = figure_colours[["muted"]]) +
  geom_errorbar(
    aes(ymin = mature_function_ci_low, ymax = mature_function_ci_high),
    width = 0, linewidth = 0.35, alpha = 0.72
  ) +
  geom_errorbar(
    aes(xmin = wnt_stem_ci_low, xmax = wnt_stem_ci_high),
    width = 0, orientation = "y", linewidth = 0.35, alpha = 0.72
  ) +
  geom_point(aes(shape = n_units > 1), size = 2.15, stroke = 0.4) +
  geom_text_repel(
    aes(label = label), family = figure_font, size = 1.55,
    lineheight = 0.88, seed = 20260831, box.padding = 0.28,
    point.padding = 0.14, segment.size = 0.22,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_colour_manual(values = quadrant_colours) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16)) +
  coord_cartesian(xlim = c(-1.05, 2.02), ylim = c(-0.92, 0.78), clip = "off") +
  labs(
    x = "WNT/stem suppression (favourable →)",
    y = "Mature-function restoration (favourable →)"
  ) +
  guides(colour = "none", shape = "none") +
  theme_cb()

# e. Direction-oriented module responses.
module_order <- c("M02", "M03", "M06", "M04", "M05", "M09", "M10")
comparison_labels <- c(
  pri724_reversal_of_trametinib = "PRI-724",
  trametinib_vs_dmso = "Trametinib",
  wnt_rspo_withdrawal_in_APC_KO = "WNT− / APC-KO",
  wnt_rspo_withdrawal_in_WT = "WNT− / WT",
  ascl2_ko_vs_resting_wt = "ASCL2 KO",
  tcf7l2_KO_vs_WT = "TCF7L2 KO",
  conditional_wnt_silencing = "WNT off",
  apc_restoration_shApc = "APC restore",
  apc_restoration_shApc_Kras = "APC restore + Kras"
)
module_plot <- module_effects %>%
  filter(
    module %in% module_order,
    interpretation_role == "causal or pathway perturbation",
    comparison %in% names(comparison_labels)
  ) %>%
  mutate(
    comparison = factor(
      comparison, levels = names(comparison_labels), labels = comparison_labels
    ),
    module = factor(module, levels = rev(module_order)),
    consistent_n3 = n_units >= 3 &
      tolower(as.character(all_units_positive_reversal)) == "true"
  )

p6e <- ggplot(module_plot, aes(comparison, module, fill = mean_module_reversal)) +
  geom_tile(colour = "white", linewidth = 0.30) +
  geom_point(
    data = module_plot %>% filter(consistent_n3),
    shape = 21, size = 1.35, stroke = 0.28,
    fill = "white", colour = figure_colours[["ink"]]
  ) +
  scale_fill_gradient2(
    low = figure_colours[["normal"]], mid = "white",
    high = figure_colours[["green"]], midpoint = 0,
    limits = c(-1.4, 1.4), oob = scales::squish,
    name = "Module reversal"
  ) +
  labs(x = NULL, y = NULL) +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 4.6, angle = 42, hjust = 1),
    legend.position = "top", legend.justification = "left",
    panel.border = element_rect(
      fill = NA, colour = figure_colours[["line"]], linewidth = 0.3
    )
  )

# f. Module response-profile correlation matrix.
correlation_long <- module_correlations %>%
  filter(module %in% module_order) %>%
  tidyr::pivot_longer(
    all_of(module_order), names_to = "module_2", values_to = "correlation"
  ) %>%
  mutate(
    module = factor(module, levels = rev(module_order)),
    module_2 = factor(module_2, levels = module_order)
  )

p6f <- ggplot(correlation_long, aes(module_2, module, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 0.32) +
  geom_text(
    aes(label = sprintf("%.2f", correlation)),
    family = figure_font, size = 1.40,
    colour = ifelse(abs(correlation_long$correlation) > 0.72, "white", figure_colours[["ink"]])
  ) +
  scale_fill_gradient2(
    low = figure_colours[["normal"]], mid = "white",
    high = figure_colours[["adenoma"]], midpoint = 0,
    limits = c(-1, 1), name = "Profile correlation"
  ) +
  labs(x = NULL, y = NULL) +
  coord_equal() +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 5.2, angle = 36, hjust = 1),
    axis.text.y = element_text(size = 5.2),
    legend.position = "top", legend.justification = "left"
  )

figure <- p6a / (p6b | p6c) / (p6d | (p6e / p6f)) +
  plot_layout(heights = c(0.55, 0.92, 1.55), widths = c(1.02, 0.98)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

pathway_source <- data.frame(
  node = c(
    "WNT", "FZD_LRP5_6", "APC_AXIN_GSK3B", "beta_catenin_TCF_LEF_ASCL2",
    "WNT_stem_response", "mature_epithelial_functions"
  ),
  role = c(
    "ligand", "receptor", "destruction_complex", "nuclear_regulation",
    "response_coordinate_1", "response_coordinate_2"
  )
)
write_source(pathway_source, source_dir, "figure6a_pathway_and_response_coordinates.tsv")
write_source(apc_samples, source_dir, "figure6b_apc_wnt_donor_scores.tsv")
write_source(apc_plot, source_dir, "figure6c_apc_wnt_contrasts.tsv")
write_source(two_axis, source_dir, "figure6d_two_coordinate_perturbations.tsv")
write_source(module_plot, source_dir, "figure6e_module_perturbation_responses.tsv")
write_source(correlation_long, source_dir, "figure6f_module_profile_correlations.tsv")

export_cb_figure(
  figure, out_dir, "figure6_separable_perturbation_responses",
  width_mm = 178, height_mm = 218
)

message("Communications Biology v5 Figure 6 exported")
