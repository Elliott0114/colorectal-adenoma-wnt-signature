#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")
suppressPackageStartupMessages(library(tidyr))

root <- normalizePath(".", mustWork = TRUE)
state_root <- file.path(root, "results", "state_aware_program_v1")
module_root <- file.path(
  state_root, "functional_architecture_exploratory_v2_1"
)
two_axis_root <- file.path(
  state_root, "identity_reversal_target_prioritization_v1"
)
out_dir <- file.path(root, "figures", "communications_biology_v5.0")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

module_order <- c("M02", "M03", "M06", "M04", "M05", "M09", "M10")
comparison_order <- c(
  "pri724_reversal_of_trametinib", "trametinib_vs_dmso",
  "wnt_rspo_withdrawal_in_APC_KO", "wnt_rspo_withdrawal_in_WT",
  "ascl2_ko_vs_resting_wt", "tcf7l2_KO_vs_WT",
  "conditional_wnt_silencing", "apc_restoration_shApc",
  "apc_restoration_shApc_Kras"
)
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
all_labels <- c(
  comparison_labels,
  doxycycline_control_shRenilla = "Doxycycline control",
  tcf7l2_heterozygous_vs_WT = "TCF7L2 heterozygous"
)

# a. Outcome-independent module coverage in every perturbation dataset.
coverage <- read_tsv(file.path(
  module_root, "module_perturbation_protein", "module_perturbation_coverage.tsv"
)) %>%
  filter(module %in% module_order) %>%
  mutate(
    dataset = recode(
      dataset,
      GSE135328_HCT116 = "GSE135328 HCT116",
      GSE135328_HT29 = "GSE135328 HT29"
    ),
    dataset = factor(dataset, levels = rev(unique(dataset))),
    module = factor(module, levels = module_order)
  )

p8a <- ggplot(coverage, aes(module, dataset, fill = coverage_fraction)) +
  geom_tile(colour = "white", linewidth = 0.32) +
  geom_text(
    aes(label = percent(coverage_fraction, accuracy = 1)),
    family = figure_font, size = 1.45
  ) +
  scale_fill_gradient(
    low = figure_colours[["pale"]], high = figure_colours[["adenoma"]],
    limits = c(0, 1), name = "Coverage"
  ) +
  labs(x = "Co-expression module", y = NULL) +
  theme_minimal(base_family = figure_font, base_size = 6.4) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top", legend.justification = "left",
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

# b. Donor-matched APC-by-WNT contrasts.
apc <- read_tsv(file.path(
  state_root, "external_validation", "apc_organoid_effects.tsv"
)) %>%
  filter(
    signature_id == "state_shared_1843",
    comparison %in% c(
      "WT_withdrawal", "APC_vs_WT_without_Wnt", "genotype_by_Wnt_interaction"
    )
  ) %>%
  mutate(
    comparison = factor(
      recode(
        comparison,
        WT_withdrawal = "WNT withdrawal in WT",
        APC_vs_WT_without_Wnt = "APC-KO vs WT, WNT−",
        genotype_by_Wnt_interaction = "Genotype × WNT"
      ),
      levels = rev(c(
        "WNT withdrawal in WT", "APC-KO vs WT, WNT−", "Genotype × WNT"
      ))
    )
  )

p8b <- ggplot(apc, aes(mean_difference, comparison)) +
  geom_vline(
    xintercept = 0, linetype = 2, linewidth = 0.35,
    colour = figure_colours[["muted"]]
  ) +
  geom_errorbar(
    aes(xmin = min_difference, xmax = max_difference),
    width = 0, orientation = "y", linewidth = 0.65,
    colour = figure_colours[["adenoma"]]
  ) +
  geom_point(size = 2.2, colour = figure_colours[["adenoma"]]) +
  labs(x = "Mean donor difference (range; n = 3)", y = NULL) +
  theme_cb()

# c. Complete two-coordinate effects and intervals.
two_axis <- read_tsv(file.path(
  two_axis_root, "perturbation_two_component_summary.tsv"
)) %>%
  mutate(
    display = unname(all_labels[comparison]),
    display = factor(display, levels = rev(unname(all_labels)))
  )

two_axis_long <- bind_rows(
  two_axis %>% transmute(
    dataset, comparison, display, n_units,
    axis = "WNT/stem suppression", estimate = mean_wnt_stem_suppression,
    ci_low = wnt_stem_ci_low, ci_high = wnt_stem_ci_high
  ),
  two_axis %>% transmute(
    dataset, comparison, display, n_units,
    axis = "Mature-function restoration", estimate = mean_mature_function_restoration,
    ci_low = mature_function_ci_low, ci_high = mature_function_ci_high
  )
) %>%
  mutate(axis = factor(
    axis, levels = c("WNT/stem suppression", "Mature-function restoration")
  ))

p8c <- ggplot(two_axis_long, aes(estimate, display, colour = axis, shape = axis)) +
  geom_vline(
    xintercept = 0, linetype = 2, linewidth = 0.35,
    colour = figure_colours[["muted"]]
  ) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y",
    linewidth = 0.52, position = position_dodge(width = 0.38)
  ) +
  geom_point(size = 1.95, position = position_dodge(width = 0.38)) +
  scale_colour_manual(values = c(
    "WNT/stem suppression" = figure_colours[["adenoma"]],
    "Mature-function restoration" = figure_colours[["normal"]]
  )) +
  scale_shape_manual(values = c(
    "WNT/stem suppression" = 16, "Mature-function restoration" = 18
  )) +
  labs(x = "Direction-oriented perturbation effect (favourable →)", y = NULL) +
  theme_cb(base_size = 6.4) +
  theme(legend.position = "top", legend.justification = "left")

# d. Complete module-by-perturbation response matrix.
module_plot <- read_tsv(file.path(
  module_root, "module_perturbation_protein", "module_perturbation_summary.tsv"
)) %>%
  filter(
    module %in% module_order,
    interpretation_role == "causal or pathway perturbation",
    comparison %in% comparison_order
  ) %>%
  mutate(
    comparison = factor(
      comparison, levels = comparison_order, labels = comparison_labels[comparison_order]
    ),
    module = factor(module, levels = rev(module_order)),
    consistent_n3 = n_units >= 3 &
      tolower(as.character(all_units_positive_reversal)) == "true"
  )

p8d <- ggplot(module_plot, aes(comparison, module, fill = mean_module_reversal)) +
  geom_tile(colour = "white", linewidth = 0.30) +
  geom_point(
    data = module_plot %>% filter(consistent_n3),
    shape = 21, size = 1.45, stroke = 0.28,
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
    axis.text.x = element_text(size = 5.0, angle = 38, hjust = 1),
    legend.position = "top", legend.justification = "left"
  )

# e. Unit-level directional agreement for comparisons with replicated units.
unit_effects <- read_tsv(file.path(
  two_axis_root, "perturbation_two_component_unit_effects.tsv"
)) %>%
  filter(comparison %in% two_axis$comparison) %>%
  mutate(display = unname(all_labels[comparison]))

direction_counts <- bind_rows(
  unit_effects %>% transmute(
    comparison, display, unit_id, axis = "WNT/stem suppression",
    favourable = wnt_stem_suppression > 0
  ),
  unit_effects %>% transmute(
    comparison, display, unit_id, axis = "Mature-function restoration",
    favourable = mature_function_restoration > 0
  )
) %>%
  group_by(comparison, display, axis) %>%
  summarise(
    n_units = n(), n_favourable = sum(favourable, na.rm = TRUE),
    fraction_favourable = n_favourable / n_units, .groups = "drop"
  ) %>%
  filter(n_units > 1) %>%
  mutate(
    display = factor(display, levels = rev(unname(all_labels))),
    axis = factor(
      axis, levels = c("WNT/stem suppression", "Mature-function restoration")
    )
  )

p8e <- ggplot(
  direction_counts,
  aes(fraction_favourable, display, colour = axis, shape = axis)
) +
  geom_segment(
    aes(x = 0, xend = fraction_favourable, yend = display),
    position = position_dodge(width = 0.42),
    linewidth = 0.52, colour = figure_colours[["line"]]
  ) +
  geom_point(size = 2.0, position = position_dodge(width = 0.42)) +
  geom_text(
    aes(label = paste0(n_favourable, "/", n_units)),
    position = position_dodge(width = 0.42), nudge_x = 0.055,
    hjust = 0, family = figure_font, size = 1.55,
    colour = figure_colours[["muted"]], show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "WNT/stem suppression" = figure_colours[["adenoma"]],
    "Mature-function restoration" = figure_colours[["normal"]]
  )) +
  scale_shape_manual(values = c(
    "WNT/stem suppression" = 16, "Mature-function restoration" = 18
  )) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1), limits = c(0, 1.16),
    breaks = c(0, 0.5, 1)
  ) +
  labs(x = "Independent units moving favourably", y = NULL) +
  theme_cb(base_size = 6.4) +
  theme(legend.position = "none")

figure <- (p8a | p8b) / p8c / (p8d | p8e) +
  plot_layout(heights = c(0.82, 1.15, 1.18), widths = c(1.06, 0.94)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

write_source(coverage, source_dir, "figureS8a_module_coverage.tsv")
write_source(apc, source_dir, "figureS8b_apc_organoid_contrasts.tsv")
write_source(two_axis_long, source_dir, "figureS8c_two_coordinate_effects.tsv")
write_source(module_plot, source_dir, "figureS8d_module_perturbation_matrix.tsv")
write_source(direction_counts, source_dir, "figureS8e_unit_direction_counts.tsv")

export_cb_figure(
  figure, out_dir, "figureS8_complete_perturbation_context",
  width_mm = 178, height_mm = 235
)

message("Communications Biology v5 Supplementary Figure 8 exported")
