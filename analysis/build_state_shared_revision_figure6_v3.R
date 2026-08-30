#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")

root <- normalizePath(".", mustWork = TRUE)
panel_root <- file.path(root, "results", "state_aware_program_v1", "panel_derivation")
heldout_root <- file.path(root, "results", "state_aware_program_v1", "heldout_validation")
compact_root <- file.path(root, "results", "state_shared_revision_v2", "compact_rank")
external_rank_root <- file.path(root, "results", "state_shared_revision_v2", "external_rank")
out_dir <- file.path(root, "figures", "communications_biology_v2.1")
source_dir <- file.path(out_dir, "source_data")

panel <- read_tsv(file.path(panel_root, "compact_state_shared_panel_frozen.tsv"))
oof_curve <- read_tsv(file.path(panel_root, "discovery_grouped_oof_fidelity_curve.tsv"))
heldout_genes <- read_tsv(file.path(heldout_root, "heldout_compact_panel_gene_validation.tsv"))
heldout_fidelity <- read_tsv(file.path(compact_root, "heldout_single_sample_rank_fidelity.tsv"))
random_internal <- read_tsv(file.path(compact_root, "random_eight_gene_benchmark.tsv"))
random_internal_summary <- read_tsv(file.path(compact_root, "random_eight_gene_benchmark_summary.tsv"))
random_external <- read_tsv(file.path(external_rank_root, "random_eight_gene_benchmark.tsv"))
external_random_summary <- read_tsv(file.path(external_rank_root, "random_eight_gene_benchmark_summary.tsv"))

# a. Prespecified reduction from the frozen programme to one tractable candidate.
p6a <- ggplot() +
  annotate("rect", xmin = 0.12, xmax = 1.75, ymin = 2.55, ymax = 3.45,
           fill = figure_colours[["pale_orange"]], colour = figure_colours[["adenoma"]], linewidth = 0.5) +
  annotate("text", x = 0.94, y = 3.10, label = "1,843 genes", family = figure_font,
           fontface = "bold", size = 2.7, colour = figure_colours[["adenoma"]]) +
  annotate("text", x = 0.94, y = 2.78, label = "frozen biological programme",
           family = figure_font, size = 1.65, colour = figure_colours[["muted"]]) +
  annotate("segment", x = 0.94, xend = 0.94, y = 2.50, yend = 2.10,
           colour = figure_colours[["grey"]], linewidth = 0.55,
           arrow = grid::arrow(type = "closed", length = grid::unit(1.5, "mm"))) +
  annotate("rect", xmin = 0.28, xmax = 1.60, ymin = 1.18, ymax = 2.05,
           fill = figure_colours[["pale_blue"]], colour = figure_colours[["normal"]], linewidth = 0.5) +
  annotate("text", x = 0.94, y = 1.72, label = "53 genes", family = figure_font,
           fontface = "bold", size = 2.55, colour = figure_colours[["normal"]]) +
  annotate("text", x = 0.94, y = 1.40, label = "protein coding · six platforms\n15 up · 38 down",
           family = figure_font, size = 1.65, lineheight = 0.92,
           colour = figure_colours[["muted"]]) +
  annotate("segment", x = 0.94, xend = 0.94, y = 1.13, yend = 0.75,
           colour = figure_colours[["grey"]], linewidth = 0.55,
           arrow = grid::arrow(type = "closed", length = grid::unit(1.5, "mm"))) +
  annotate("rect", xmin = 0.42, xmax = 1.46, ymin = 0.02, ymax = 0.70,
           fill = figure_colours[["pale_gold"]], colour = figure_colours[["gold"]], linewidth = 0.5) +
  annotate("text", x = 0.94, y = 0.44, label = "8-gene candidate", family = figure_font,
           fontface = "bold", size = 2.55, colour = figure_colours[["gold"]]) +
  annotate("text", x = 0.94, y = 0.19, label = "4 balanced pairs",
           family = figure_font, size = 1.65, colour = figure_colours[["muted"]]) +
  annotate("text", x = 1.86, y = 1.60,
           label = "Donor-held-out\nreconstruction\ndata-defined knee",
           hjust = 0, family = figure_font, size = 1.75, lineheight = 0.96,
           colour = figure_colours[["ink"]]) +
  coord_cartesian(xlim = c(0, 2.75), ylim = c(-0.05, 3.55), clip = "off") +
  theme_void(base_family = figure_font) +
  theme(plot.margin = margin(2, 2, 2, 2, "mm"))

# b. Selected genes, discovery/validation effects and bootstrap stability.
gene_plot <- heldout_genes %>%
  mutate(
    arm = factor(arm, levels = c("up", "down")),
    gene = factor(gene, levels = rev(panel$gene[order(panel$pair_step, panel$arm)]))
  ) %>%
  select(gene, arm, common_effect, validation_common_effect, donor_bootstrap_selection_frequency) %>%
  tidyr::pivot_longer(c(common_effect, validation_common_effect), names_to = "partition", values_to = "effect") %>%
  mutate(partition = factor(partition, levels = c("common_effect", "validation_common_effect"), labels = c("Discovery", "Donor-disjoint validation")))

p6b <- ggplot(gene_plot, aes(effect, gene, colour = partition, shape = partition)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_point(size = 2.0, position = position_dodge(width = 0.38), na.rm = TRUE) +
  scale_colour_manual(values = c(Discovery = figure_colours[["grey"]], "Donor-disjoint validation" = figure_colours[["adenoma"]])) +
  scale_shape_manual(values = c(Discovery = 16, "Donor-disjoint validation" = 18)) +
  geom_text(
    data = gene_plot %>% filter(partition == "Discovery"),
    aes(x = 7.0, y = gene, label = sprintf("%.0f%%", 100 * donor_bootstrap_selection_frequency)),
    inherit.aes = FALSE, hjust = 1, family = figure_font, size = 1.65,
    colour = figure_colours[["muted"]]
  ) +
  scale_x_continuous(limits = c(-1.6, 7.2)) +
  labs(x = "Common effect", y = NULL) +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

# c. Panel-size fidelity curve in discovery donors.
p6c <- ggplot(oof_curve, aes(total_genes, monotone_oof_spearman)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = figure_colours[["pale_orange"]], colour = NA) +
  geom_line(colour = figure_colours[["adenoma"]], linewidth = 0.75) +
  geom_point(size = 1.55, colour = figure_colours[["adenoma"]]) +
  geom_vline(xintercept = 8, colour = figure_colours[["gold"]], linetype = 2, linewidth = 0.5) +
  annotate("point", x = 8, y = oof_curve$monotone_oof_spearman[oof_curve$total_genes == 8],
           size = 2.5, colour = figure_colours[["gold"]]) +
  annotate("text", x = 8.7, y = 0.88, label = "8 genes", hjust = 0,
           family = figure_font, size = 1.8, colour = figure_colours[["gold"]]) +
  scale_x_continuous(breaks = c(2, 8, 16, 24, 30)) +
  labs(x = "Genes in balanced candidate", y = "Donor-held-out fidelity") +
  theme_cb()

# d. Label-independent single-sample fidelity in donor-disjoint profiles.
heldout_plot <- heldout_fidelity %>%
  mutate(
    scope = factor(scope, levels = rev(c("all", "ABS", "GOB", "TAC")), labels = rev(c("All states", "ABS", "GOB", "TAC"))),
    target_label = factor(target, levels = c("full_programme_score", "single_sample_full_rank_score"), labels = c("Discovery-standardised full", "Single-sample full rank"))
  )

p6d <- ggplot(heldout_plot, aes(spearman, scope, colour = target_label, shape = target_label)) +
  geom_errorbarh(aes(xmin = donor_bootstrap_ci_low, xmax = donor_bootstrap_ci_high), height = 0, linewidth = 0.58, position = position_dodge(width = 0.36)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.36)) +
  scale_colour_manual(values = c("Discovery-standardised full" = figure_colours[["grey"]], "Single-sample full rank" = figure_colours[["normal"]])) +
  scale_shape_manual(values = c("Discovery-standardised full" = 16, "Single-sample full rank" = 18)) +
  scale_x_continuous(limits = c(0.4, 0.96), breaks = c(0.5, 0.7, 0.9)) +
  labs(x = "Spearman fidelity", y = NULL) +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left", legend.text = element_text(size = 5.5))

# e. Internal and external direction-balanced random-panel benchmarks.
benchmark_values <- bind_rows(
  random_internal %>%
    transmute(context = "Donor-disjoint validation", fidelity = spearman_with_full_programme),
  random_external %>%
    transmute(context = "Median across external cohorts", fidelity = median_cohort_spearman)
) %>%
  mutate(context = factor(context, levels = c("Donor-disjoint validation", "Median across external cohorts")))

benchmark_reference <- data.frame(
  context = factor(
    c("Donor-disjoint validation", "Median across external cohorts"),
    levels = levels(benchmark_values$context)
  ),
  observed = c(
    random_internal_summary$observed_compact_spearman,
    external_random_summary$observed_median_cohort_spearman
  ),
  random_q95 = c(random_internal_summary$random_q95, external_random_summary$random_q95),
  empirical_p = c(
    random_internal_summary$empirical_upper_tail_p,
    external_random_summary$empirical_upper_tail_p
  )
)

p6e <- ggplot(benchmark_values, aes(fidelity)) +
  geom_histogram(bins = 32, fill = figure_colours[["line"]], colour = "white", linewidth = 0.22) +
  geom_vline(data = benchmark_reference, aes(xintercept = random_q95), colour = figure_colours[["normal"]], linetype = 2, linewidth = 0.6) +
  geom_vline(data = benchmark_reference, aes(xintercept = observed), colour = figure_colours[["adenoma"]], linewidth = 0.7) +
  geom_text(
    data = benchmark_reference,
    aes(x = observed, y = Inf, label = sprintf("Observed\nP = %.3f", empirical_p)),
    hjust = -0.08, vjust = 1.2, family = figure_font, size = 1.6,
    colour = figure_colours[["adenoma"]], inherit.aes = FALSE
  ) +
  facet_wrap(~context, ncol = 1, scales = "free_y") +
  labs(x = "Direction-balanced random-panel fidelity", y = "Panels") +
  theme_cb() +
  theme(strip.text = element_text(size = 5.8, face = "bold"))

top <- p6a | p6b
bottom <- p6c | p6d | p6e
figure <- top / bottom +
  plot_layout(heights = c(0.95, 1.05)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

selection_source <- data.frame(
  stage = c("Frozen programme", "Measurable candidate universe", "Candidate reduced readout"),
  n_genes = c(1843, 53, 8),
  rule = c("confidence-defined programme", "protein coding and present on six platforms", "donor-held-out data-defined knee")
)
write_source(selection_source, source_dir, "figure6a_reduction_path.tsv")
write_source(gene_plot, source_dir, "figure6b_candidate_genes.tsv")
write_source(oof_curve, source_dir, "figure6c_oof_fidelity_curve.tsv")
write_source(heldout_plot, source_dir, "figure6d_heldout_single_sample_fidelity.tsv")
write_source(benchmark_values, source_dir, "figure6e_random_panel_benchmark_values.tsv")
write_source(benchmark_reference, source_dir, "figure6e_random_panel_benchmark_reference.tsv")

export_cb_figure(
  figure, out_dir,
  "figure6_reduced_measurement_candidate",
  width_mm = 178, height_mm = 165
)

cat("Figure 6 exported to ", out_dir, "\n", sep = "")
