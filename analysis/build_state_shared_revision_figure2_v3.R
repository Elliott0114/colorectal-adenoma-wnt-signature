#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")

root <- normalizePath(".", mustWork = TRUE)
donor_root <- file.path(root, "results", "state_shared_revision_v2", "donor_site")
fine_root <- file.path(root, "results", "state_shared_revision_v2", "fine_state_models")
out_dir <- file.path(root, "figures", "communications_biology_v2.1")
source_dir <- file.path(out_dir, "source_data")

gene_effects <- read_tsv(file.path(donor_root, "donor_disjoint_common_gene_effects.tsv.gz")) %>%
  mutate(strict_state_shared = tolower(as.character(strict_state_shared)) == "true") %>%
  filter(strict_state_shared, is.finite(common_effect), is.finite(discovery_common_effect))
broad_effects <- read_tsv(file.path(donor_root, "donor_disjoint_score_effects.tsv")) %>%
  filter(score == "full_programme_score")
fine_effects <- read_tsv(file.path(fine_root, "fine_state_adjusted_route_effects.tsv")) %>%
  mutate(singular_fit = tolower(as.character(singular_fit)) == "true")
inventory <- read_tsv(file.path(fine_root, "fine_state_inventory.tsv"))
decomposition <- read_tsv(file.path(fine_root, "programme_composition_decomposition.tsv"))
replication_summary <- read_tsv(file.path(donor_root, "donor_disjoint_replication_summary.tsv"))

# a. Gene-level transfer into the donor-disjoint partition.
p2a <- ggplot(
  gene_effects,
  aes(discovery_common_effect, common_effect, colour = shared_direction)
) +
  geom_hline(yintercept = 0, colour = figure_colours[["line"]], linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = figure_colours[["line"]], linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, colour = figure_colours[["grey"]], linetype = 2, linewidth = 0.4) +
  geom_point(size = 0.8, alpha = 0.55) +
  scale_colour_manual(values = c(up = figure_colours[["adenoma"]], down = figure_colours[["normal"]])) +
  annotate(
    "text", x = -1.85, y = 2.9,
    label = sprintf(
      "ρ = %.3f\n%.1f%% common direction",
      replication_summary$discovery_validation_effect_spearman,
      100 * replication_summary$common_direction_match_fraction
    ),
    hjust = 0, vjust = 1, family = figure_font, size = 2.0,
    colour = figure_colours[["muted"]]
  ) +
  labs(x = "Discovery common effect", y = "Donor-disjoint common effect") +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

# b. Broad-state effects before and after anatomical and demographic adjustment.
broad_plot <- broad_effects %>%
  mutate(
    scope = factor(scope, levels = rev(c("all", "ABS", "GOB", "TAC")), labels = rev(c("All states", "ABS", "GOB", "TAC"))),
    model_label = recode(model, unadjusted = "Unadjusted", site_age_sex_adjusted = "Adjusted")
  )

p2b <- ggplot(broad_plot, aes(estimate, scope, colour = model_label, shape = model_label)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.6, position = position_dodge(width = 0.32)) +
  geom_point(size = 2.15, position = position_dodge(width = 0.32)) +
  scale_colour_manual(values = c(Unadjusted = figure_colours[["grey"]], Adjusted = figure_colours[["adenoma"]])) +
  scale_shape_manual(values = c(Unadjusted = 16, Adjusted = 18)) +
  labs(x = "Adenoma effect on full-programme score", y = NULL) +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

# c. Programme-blind fine-state proportions in donor-disjoint data.
fine_proportions <- inventory %>%
  filter(partition == "validation", k == 4) %>%
  group_by(broad_state, route) %>%
  mutate(
    proportion = n_cells / sum(n_cells),
    fine_label = sub("^[A-Z]+_", "", fine_state),
    route_label = factor(route, levels = c("normal", "conventional_adenoma"), labels = c("Normal", "Adenoma"))
  ) %>%
  ungroup()

p2c <- ggplot(fine_proportions, aes(route_label, proportion, fill = fine_label)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  facet_grid(. ~ broad_state) +
  scale_fill_manual(values = c(F1 = "#315E73", F2 = "#6F9BB1", F3 = "#D6A23A", F4 = "#D45D3F")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Cells assigned to each fine state") +
  theme_cb() +
  theme(
    legend.position = "top", legend.justification = "left",
    panel.spacing.x = unit(6, "mm")
  )

# d. Fine-state-adjusted effects across clustering resolutions.
fine_effect_plot <- fine_effects %>%
  filter(model == "unweighted") %>%
  mutate(
    scope = factor(scope, levels = rev(c("all", "ABS", "GOB", "TAC")), labels = rev(c("All states", "ABS", "GOB", "TAC"))),
    resolution = factor(paste0("k = ", k), levels = c("k = 3", "k = 4", "k = 5"))
  )

p2d <- ggplot(fine_effect_plot, aes(estimate, scope, colour = resolution)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.58, position = position_dodge(width = 0.42)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.42)) +
  scale_colour_manual(values = c("k = 3" = figure_colours[["normal"]], "k = 4" = figure_colours[["adenoma"]], "k = 5" = figure_colours[["purple"]])) +
  labs(x = "Fine-state-adjusted route effect", y = NULL) +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

# e. Primary exact decomposition in donor-disjoint validation.
primary_decomposition <- decomposition %>%
  filter(
    partition == "validation", k == 4, scope == "all",
    decomposition_type == "fine_state",
    component %in% c("total", "composition", "within")
  ) %>%
  mutate(
    component_label = factor(
      recode(component, total = "Observed total", composition = "Fine-state composition", within = "Within fine states"),
      levels = rev(c("Observed total", "Fine-state composition", "Within fine states"))
    ),
    colour_group = ifelse(component == "within", "Within", "Other")
  )

p2e <- ggplot(primary_decomposition, aes(estimate, component_label, colour = colour_group)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.68) +
  geom_point(size = 2.25) +
  scale_colour_manual(values = c(Within = figure_colours[["adenoma"]], Other = figure_colours[["grey"]])) +
  guides(colour = "none") +
  labs(x = "Contribution to adenoma–normal difference", y = NULL) +
  theme_cb()

# f. Fraction of the total difference attributed to within-state change.
within_fraction <- decomposition %>%
  filter(
    decomposition_type == "fine_state", scope == "all",
    component %in% c("total", "within")
  ) %>%
  select(partition, k, component, estimate) %>%
  tidyr::pivot_wider(names_from = component, values_from = estimate) %>%
  mutate(
    fraction = within / total,
    partition_label = factor(partition, levels = c("discovery", "validation"), labels = c("Discovery", "Validation"))
  )

p2f <- ggplot(within_fraction, aes(k, fraction, colour = partition_label, group = partition_label)) +
  geom_hline(yintercept = 0.5, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.2) +
  geom_text(
    data = within_fraction %>% filter(k == 5),
    aes(label = partition_label), hjust = 0, nudge_x = 0.08,
    family = figure_font, size = 1.8, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(Discovery = figure_colours[["grey"]], Validation = figure_colours[["adenoma"]])) +
  scale_x_continuous(breaks = c(3, 4, 5), limits = c(3, 5.72)) +
  scale_y_continuous(limits = c(0.45, 0.9), labels = percent_format(accuracy = 1)) +
  labs(x = "Fine states per broad state", y = "Within-state fraction") +
  guides(colour = "none") +
  theme_cb() +
  theme(plot.margin = margin(2, 4, 2, 2, "mm"))

top <- p2a | p2b
bottom <- p2d | p2e | p2f
figure <- top / p2c / bottom +
  plot_layout(heights = c(1, 0.72, 1)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

write_source(gene_effects, source_dir, "figure2a_donor_disjoint_gene_effects.tsv")
write_source(broad_plot, source_dir, "figure2b_broad_state_effects.tsv")
write_source(fine_proportions, source_dir, "figure2c_fine_state_proportions.tsv")
write_source(fine_effect_plot, source_dir, "figure2d_fine_state_adjusted_effects.tsv")
write_source(primary_decomposition, source_dir, "figure2e_primary_decomposition.tsv")
write_source(within_fraction, source_dir, "figure2f_within_fraction.tsv")

export_cb_figure(
  figure, out_dir,
  "figure2_donor_disjoint_identity_remodelling",
  width_mm = 178, height_mm = 200
)

cat("Figure 2 exported to ", out_dir, "\n", sep = "")
