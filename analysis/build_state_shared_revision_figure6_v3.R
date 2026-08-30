#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")

root <- normalizePath(".", mustWork = TRUE)
panel_root <- file.path(root, "results", "state_aware_program_v1", "panel_derivation")
heldout_root <- file.path(root, "results", "state_aware_program_v1", "heldout_validation")
external_root <- file.path(root, "results", "state_aware_program_v1", "external_validation")
compact_root <- file.path(root, "results", "state_shared_revision_v2", "compact_rank")
external_rank_root <- file.path(root, "results", "state_shared_revision_v2", "external_rank")
out_dir <- file.path(root, "figures", "communications_biology_v2.0")
source_dir <- file.path(out_dir, "source_data")

panel <- read_tsv(file.path(panel_root, "compact_state_shared_panel_frozen.tsv"))
oof_curve <- read_tsv(file.path(panel_root, "discovery_grouped_oof_fidelity_curve.tsv"))
heldout_genes <- read_tsv(file.path(heldout_root, "heldout_compact_panel_gene_validation.tsv"))
heldout_fidelity <- read_tsv(file.path(compact_root, "heldout_single_sample_rank_fidelity.tsv"))
random_internal <- read_tsv(file.path(compact_root, "random_eight_gene_benchmark.tsv"))
random_internal_summary <- read_tsv(file.path(compact_root, "random_eight_gene_benchmark_summary.tsv"))
external_fidelity <- read_tsv(file.path(external_rank_root, "external_rank_fidelity.tsv"))
external_random_summary <- read_tsv(file.path(external_rank_root, "random_eight_gene_benchmark_summary.tsv"))
ffpe <- read_tsv(file.path(external_root, "ffpe_sample_scores.tsv.gz")) %>%
  filter(tissue_group %in% c("adenoma", "normal"), signature_id == "compact_8") %>%
  group_by(patient_id, tissue_group) %>%
  summarise(score = mean(programme_score), .groups = "drop")

# a. Objective reduction from the frozen programme to a portable candidate.
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
  annotate("text", x = 0.94, y = 0.44, label = "8 genes", family = figure_font,
           fontface = "bold", size = 2.55, colour = figure_colours[["gold"]]) +
  annotate("text", x = 0.94, y = 0.19, label = "4 balanced pairs",
           family = figure_font, size = 1.65, colour = figure_colours[["muted"]]) +
  annotate("text", x = 1.86, y = 1.60,
           label = "Donor-held-out\nreconstruction\nKneedle + one-SE",
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

# e. Internal random-panel benchmark.
p6e <- ggplot(random_internal, aes(spearman_with_full_programme)) +
  geom_histogram(bins = 34, fill = figure_colours[["line"]], colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = random_internal_summary$random_q95, colour = figure_colours[["normal"]], linetype = 2, linewidth = 0.65) +
  geom_vline(xintercept = random_internal_summary$observed_compact_spearman, colour = figure_colours[["adenoma"]], linewidth = 0.75) +
  annotate("text", x = random_internal_summary$observed_compact_spearman, y = Inf,
           label = "observed", hjust = -0.08, vjust = 1.25, family = figure_font,
           size = 1.75, colour = figure_colours[["adenoma"]]) +
  annotate("text", x = random_internal_summary$random_q95, y = Inf,
           label = "random q95", hjust = 1.05, vjust = 2.45, family = figure_font,
           size = 1.65, colour = figure_colours[["normal"]]) +
  labs(x = "Random eight-gene panel fidelity", y = "Panels") +
  theme_cb()

# f. External single-sample fidelity with the external random benchmark shown.
external_plot <- external_fidelity %>%
  mutate(cohort = factor(cohort, levels = rev(c("GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820"))))

p6f <- ggplot(external_plot, aes(spearman_compact_rank_vs_reference_full, cohort)) +
  geom_vline(xintercept = external_random_summary$random_q95, colour = figure_colours[["normal"]], linetype = 2, linewidth = 0.55) +
  geom_vline(xintercept = external_random_summary$observed_median_cohort_spearman, colour = figure_colours[["adenoma"]], linewidth = 0.55) +
  geom_segment(aes(x = 0.55, xend = spearman_compact_rank_vs_reference_full, yend = cohort), colour = figure_colours[["line"]], linewidth = 0.65) +
  geom_point(size = 2.2, colour = figure_colours[["adenoma"]]) +
  scale_x_continuous(limits = c(0.52, 0.91), breaks = c(0.6, 0.7, 0.8, 0.9)) +
  labs(x = "External compact-to-full fidelity", y = NULL) +
  theme_cb()

# g. Reduced readout in paired FFPE tissue.
ffpe_pairs <- ffpe %>%
  tidyr::pivot_wider(names_from = tissue_group, values_from = score) %>%
  filter(!is.na(normal), !is.na(adenoma)) %>%
  tidyr::pivot_longer(c(normal, adenoma), names_to = "tissue_group", values_to = "score") %>%
  mutate(tissue_group = factor(tissue_group, levels = c("normal", "adenoma"), labels = c("Adjacent mucosa", "Adenoma")))
ffpe_summary <- ffpe_pairs %>%
  group_by(tissue_group) %>%
  summarise(median = median(score), .groups = "drop")

p6g <- ggplot(ffpe_pairs, aes(tissue_group, score, group = patient_id)) +
  geom_line(colour = figure_colours[["line"]], linewidth = 0.35, alpha = 0.75) +
  geom_point(aes(colour = tissue_group), size = 1.15, alpha = 0.85) +
  geom_point(data = ffpe_summary, aes(tissue_group, median, group = 1), shape = 23,
             size = 3.0, fill = "white", colour = figure_colours[["ink"]], stroke = 0.55) +
  scale_colour_manual(values = c("Adjacent mucosa" = figure_colours[["normal"]], Adenoma = figure_colours[["adenoma"]])) +
  guides(colour = "none") +
  annotate("text", x = 1.5, y = max(ffpe_pairs$score) + 0.12,
           label = "47/51 pairs increased", family = figure_font, size = 1.9,
           colour = figure_colours[["muted"]]) +
  labs(x = NULL, y = "Eight-gene FFPE score") +
  theme_cb()

top <- p6a | p6b
middle <- p6c | p6d | p6e
bottom <- p6f | p6g
figure <- top / middle / bottom +
  plot_layout(heights = c(0.92, 1, 0.9)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

selection_source <- data.frame(
  stage = c("Frozen programme", "Portable candidate universe", "Balanced reduced candidate"),
  n_genes = c(1843, 53, 8),
  rule = c("confidence-defined programme", "protein coding and present on six platforms", "donor-held-out Kneedle and one-SE")
)
write_source(selection_source, source_dir, "figure6a_reduction_path.tsv")
write_source(gene_plot, source_dir, "figure6b_candidate_genes.tsv")
write_source(oof_curve, source_dir, "figure6c_oof_fidelity_curve.tsv")
write_source(heldout_plot, source_dir, "figure6d_heldout_single_sample_fidelity.tsv")
write_source(random_internal, source_dir, "figure6e_internal_random_benchmark.tsv")
write_source(external_plot, source_dir, "figure6f_external_fidelity.tsv")
write_source(ffpe_pairs, source_dir, "figure6g_ffpe_pairs.tsv")

export_cb_figure(
  figure, out_dir,
  "figure6_reduced_measurement_candidate",
  width_mm = 178, height_mm = 205
)

cat("Figure 6 exported to ", out_dir, "\n", sep = "")
