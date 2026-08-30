#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")
suppressPackageStartupMessages(library(tidyr))

root <- normalizePath(".", mustWork = TRUE)
state_root <- file.path(root, "results", "state_aware_program_v1")
revision_root <- file.path(root, "results", "state_shared_revision_v2")
out_dir <- file.path(root, "figures", "communications_biology_v2.0")
source_dir <- file.path(out_dir, "source_data")

save_supplement <- function(plot, number, stem, height_mm = 190) {
  export_cb_figure(
    tagged(plot), out_dir,
    paste0("figureS", number, "_", stem),
    width_mm = 178, height_mm = height_mm
  )
}

state_labels <- c(ABS = "Absorptive", GOB = "Goblet", TAC = "Transit-amplifying")
route_labels <- c(conventional_adenoma = "Adenoma", normal = "Normal")

# Supplementary Figure 1: sampling, overlap exclusion and donor stability.
eligibility <- read_tsv(file.path(state_root, "discovery_pseudobulk", "state_eligibility.tsv")) %>%
  filter(cell_type %in% names(state_labels)) %>%
  select(cell_type, conventional_adenoma, normal) %>%
  pivot_longer(c(conventional_adenoma, normal), names_to = "route", values_to = "n_donors") %>%
  mutate(partition = "Discovery")
validation_scores <- read_tsv(file.path(revision_root, "donor_site", "donor_disjoint_programme_scores.tsv.gz"))
validation_support <- validation_scores %>%
  distinct(donor_id, route, cell_type) %>%
  count(cell_type, route, name = "n_donors") %>%
  mutate(partition = "Donor-disjoint validation")
support <- bind_rows(eligibility, validation_support) %>%
  mutate(
    cell_type = factor(state_labels[cell_type], levels = unname(state_labels)),
    route = factor(route_labels[route], levels = c("Normal", "Adenoma")),
    partition = factor(partition, levels = c("Discovery", "Donor-disjoint validation"))
  )
p1a <- ggplot(support, aes(cell_type, n_donors, fill = route)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  facet_wrap(~partition) +
  scale_fill_manual(values = c(Normal = figure_colours[["normal"]], Adenoma = figure_colours[["adenoma"]])) +
  labs(x = NULL, y = "Donors represented") +
  theme_cb() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")

overlap <- read_tsv(file.path(revision_root, "donor_site", "donor_overlap_exclusion_audit.tsv"))
overlap_plot <- data.frame(
  measure = rep(c("Donors", "Pseudobulk profiles"), each = 2),
  stage = rep(c("Before overlap removal", "After overlap removal"), 2),
  value = c(
    overlap$validation_donors_before, overlap$validation_donors_after,
    overlap$validation_profiles_before, overlap$validation_profiles_after
  )
) %>%
  mutate(stage = factor(stage, levels = c("Before overlap removal", "After overlap removal")))
p1b <- ggplot(overlap_plot, aes(stage, value, fill = stage)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = value), vjust = -0.35, family = figure_font, size = 2.1) +
  facet_wrap(~measure, scales = "free_y") +
  scale_fill_manual(values = c("Before overlap removal" = figure_colours[["grey"]], "After overlap removal" = figure_colours[["adenoma"]])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = NULL) +
  guides(fill = "none") +
  theme_cb() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

pseudobulk <- read_tsv(file.path(state_root, "discovery_pseudobulk", "pseudobulk_metadata.tsv")) %>%
  filter(cell_type %in% names(state_labels)) %>%
  mutate(
    cell_type = factor(state_labels[cell_type], levels = unname(state_labels)),
    route = factor(route_labels[route], levels = c("Normal", "Adenoma"))
  )
p1c <- ggplot(pseudobulk, aes(cell_type, library_size, colour = route)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.45) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.62), size = 0.8, alpha = 0.55) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_colour_manual(values = c(Normal = figure_colours[["normal"]], Adenoma = figure_colours[["adenoma"]])) +
  labs(x = NULL, y = "Discovery pseudobulk library size") +
  theme_cb() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")

primary <- read_tsv(file.path(state_root, "discovery_models", "state_specific_primary_effects.tsv.gz")) %>%
  select(gene, cell_type, primary = logFC)
paired <- read_tsv(file.path(state_root, "discovery_models", "state_specific_paired_sensitivity.tsv.gz")) %>%
  select(gene, cell_type, paired = logFC)
paired_comparison <- inner_join(primary, paired, by = c("gene", "cell_type"))
paired_rho <- cor(paired_comparison$primary, paired_comparison$paired, method = "spearman", use = "complete.obs")
p1d <- ggplot(paired_comparison, aes(primary, paired, colour = cell_type)) +
  geom_hline(yintercept = 0, colour = figure_colours[["line"]], linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = figure_colours[["line"]], linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = figure_colours[["grey"]], linewidth = 0.4) +
  geom_point(size = 0.45, alpha = 0.22) +
  annotate("text", x = -Inf, y = Inf, label = sprintf("ρ = %.3f", paired_rho), hjust = -0.08, vjust = 1.15,
           family = figure_font, size = 2.1, colour = figure_colours[["muted"]]) +
  scale_colour_manual(values = c(ABS = figure_colours[["adenoma"]], GOB = figure_colours[["normal"]], TAC = figure_colours[["purple"]])) +
  labs(x = "Primary mixed-model effect", y = "Paired-donor sensitivity effect") +
  theme_cb() +
  theme(legend.position = "top")

leaveout <- read_tsv(file.path(state_root, "donor_leaveout_stability", "donor_leaveout_rank_stability.tsv")) %>%
  filter(scope == "common_GLS") %>%
  arrange(effect_spearman) %>%
  mutate(index = row_number())
p1e <- ggplot(leaveout, aes(index, effect_spearman)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.97, ymax = 1, fill = figure_colours[["pale_blue"]], colour = NA) +
  geom_hline(yintercept = median(leaveout$effect_spearman), linetype = 2, colour = figure_colours[["normal"]], linewidth = 0.45) +
  geom_point(size = 1.9, colour = figure_colours[["normal"]]) +
  scale_y_continuous(limits = c(0.968, 1.0005), breaks = c(0.97, 0.98, 0.99, 1.00)) +
  labs(x = "Discovery donor left out", y = "Rank correlation") +
  theme_cb()

same_site <- read_tsv(file.path(revision_root, "donor_site", "same_site_paired_score_effects.tsv")) %>%
  filter(score == "full_programme_score") %>%
  mutate(cell_type = factor(state_labels[cell_type], levels = rev(unname(state_labels))))
p1f <- ggplot(same_site, aes(mean_difference, cell_type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65, colour = figure_colours[["adenoma"]]) +
  geom_point(size = 2.2, colour = figure_colours[["adenoma"]]) +
  geom_text(aes(label = paste0("n = ", n_donors)), nudge_y = -0.22, family = figure_font, size = 1.8, colour = figure_colours[["muted"]]) +
  labs(x = "Same-donor, same-site adenoma − normal", y = NULL) +
  theme_cb()

fig_s1 <- (p1a | p1b) / (p1c | p1d) / (p1e | p1f) + plot_annotation(tag_levels = "a")
save_supplement(fig_s1, 1, "sampling_and_donor_stability", 205)
write_source(support, source_dir, "figureS1a_donor_support.tsv")
write_source(overlap_plot, source_dir, "figureS1b_overlap_exclusion.tsv")
write_source(pseudobulk, source_dir, "figureS1c_library_sizes.tsv")
write_source(paired_comparison, source_dir, "figureS1d_paired_sensitivity.tsv")
write_source(leaveout, source_dir, "figureS1e_leaveout.tsv")
write_source(same_site, source_dir, "figureS1f_same_site.tsv")

# Supplementary Figure 2: transparent audit of the historical 287- and 12-gene objects.
legacy_coverage <- read_tsv(file.path(state_root, "legacy_audit", "legacy_signature_coverage_summary.tsv"))
coverage_plot <- legacy_coverage %>%
  filter(signature %in% c("original_287", "existing_12")) %>%
  select(signature, frozen_genes, testable_genes, strict_state_shared) %>%
  pivot_longer(c(frozen_genes, testable_genes, strict_state_shared), names_to = "stage", values_to = "genes") %>%
  mutate(
    signature = recode(signature, original_287 = "Historical 287", existing_12 = "Historical 12"),
    stage = factor(stage, levels = c("frozen_genes", "testable_genes", "strict_state_shared"),
                   labels = c("Frozen", "Testable in all states", "State-shared"))
  )
p2a <- ggplot(coverage_plot, aes(stage, genes, fill = signature)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = genes), position = position_dodge(width = 0.72), vjust = -0.35, family = figure_font, size = 1.8) +
  scale_fill_manual(values = c("Historical 287" = figure_colours[["grey"]], "Historical 12" = figure_colours[["gold"]])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(x = NULL, y = "Genes") +
  theme_cb() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")

legacy_gene <- read_tsv(file.path(state_root, "legacy_audit", "original_287_gene_level_audit.tsv.gz")) %>%
  filter(in_common_universe) %>%
  mutate(expected_arm = factor(expected_arm, levels = c("down", "up"), labels = c("Adenoma-down", "Adenoma-up")))
p2b <- ggplot(legacy_gene, aes(expected_arm, common_z, fill = expected_arm)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = figure_colours[["muted"]]) +
  geom_violin(width = 0.72, alpha = 0.28, colour = NA) +
  geom_boxplot(width = 0.28, outlier.shape = NA, fill = "white", linewidth = 0.45) +
  geom_jitter(width = 0.12, size = 0.65, alpha = 0.35) +
  scale_fill_manual(values = c("Adenoma-down" = figure_colours[["normal"]], "Adenoma-up" = figure_colours[["adenoma"]])) +
  guides(fill = "none") +
  labs(x = NULL, y = "State-aware common effect (z)") +
  theme_cb()

legacy_enrichment <- read_tsv(file.path(state_root, "legacy_audit", "original_287_competitive_enrichment.tsv")) %>%
  mutate(
    label = recode(gene_set, original_287_up = "Historical up arm", original_287_down = "Historical down arm"),
    evidence = -log10(FDR),
    label = factor(label, levels = rev(c("Historical up arm", "Historical down arm")))
  )
p2c <- ggplot(legacy_enrichment, aes(evidence, label, colour = Direction)) +
  geom_segment(aes(x = 0, xend = evidence, yend = label), linewidth = 0.8) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(Down = figure_colours[["normal"]], Up = figure_colours[["adenoma"]])) +
  labs(x = "−log10 competitive-enrichment FDR", y = NULL) +
  guides(colour = "none") +
  theme_cb()

legacy_null <- read_tsv(file.path(state_root, "legacy_audit", "original_287_matched_null_summary.tsv")) %>%
  mutate(
    label = recode(metric,
                   up_signed_mean_z = "Up arm", down_signed_mean_z = "Down arm", joint_signed_mean_z = "Joint"),
    label = factor(label, levels = rev(c("Up arm", "Down arm", "Joint")))
  )
p2d <- ggplot(legacy_null, aes(observed, label)) +
  geom_errorbarh(aes(xmin = null_ci_low, xmax = null_ci_high), height = 0, linewidth = 1.0, colour = figure_colours[["line"]]) +
  geom_point(size = 2.4, colour = figure_colours[["adenoma"]]) +
  geom_point(aes(x = null_mean), shape = 4, size = 2.0, stroke = 0.7, colour = figure_colours[["grey"]]) +
  labs(x = "Observed signed evidence versus matched null", y = NULL) +
  theme_cb()

historical_12 <- read_tsv(file.path(state_root, "legacy_audit", "existing_12_gene_level_audit.tsv")) %>%
  select(gene, expected_arm, posterior_mean_ABS, posterior_mean_GOB, posterior_mean_TAC) %>%
  pivot_longer(starts_with("posterior_mean_"), names_to = "state", values_to = "effect") %>%
  mutate(
    state = sub("posterior_mean_", "", state),
    gene = factor(gene, levels = rev(unique(gene))),
    state = factor(state, levels = c("ABS", "GOB", "TAC"))
  )
limit12 <- max(abs(historical_12$effect), na.rm = TRUE)
p2e <- ggplot(historical_12, aes(state, gene, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = ifelse(is.na(effect), "NA", sprintf("%.1f", effect))), family = figure_font, size = 1.8) +
  scale_fill_gradient2(low = figure_colours[["normal"]], mid = "white", high = figure_colours[["adenoma"]], midpoint = 0,
                       limits = c(-limit12, limit12), na.value = figure_colours[["pale"]], name = "Posterior effect") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = figure_font, base_size = 7) +
  theme(panel.grid = element_blank(), legend.position = "right", plot.margin = margin(2, 2, 2, 2, "mm"))

legacy_membership <- legacy_coverage %>%
  filter(signature %in% c("original_287_up", "original_287_down", "existing_12")) %>%
  mutate(
    label = recode(signature, original_287_up = "287 up", original_287_down = "287 down", existing_12 = "Historical 12"),
    outside = frozen_genes - testable_genes,
    testable_not_shared = testable_genes - strict_state_shared
  ) %>%
  select(label, strict_state_shared, testable_not_shared, outside) %>%
  pivot_longer(-label, names_to = "class", values_to = "genes") %>%
  mutate(class = factor(class, levels = c("strict_state_shared", "testable_not_shared", "outside"),
                        labels = c("State-shared", "Other testable", "Not testable")))
p2f <- ggplot(legacy_membership, aes(label, genes, fill = class)) +
  geom_col(width = 0.66) +
  scale_fill_manual(values = c("State-shared" = figure_colours[["adenoma"]],
                               "Other testable" = figure_colours[["gold"]],
                               "Not testable" = figure_colours[["line"]])) +
  labs(x = NULL, y = "Genes") +
  theme_cb() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top", legend.text = element_text(size = 5.8))

fig_s2 <- (p2a | p2b) / (p2c | p2d) / (p2e | p2f) + plot_annotation(tag_levels = "a")
save_supplement(fig_s2, 2, "historical_gene_set_audit", 210)
write_source(coverage_plot, source_dir, "figureS2a_legacy_coverage.tsv")
write_source(legacy_gene, source_dir, "figureS2b_legacy_gene_effects.tsv")
write_source(legacy_enrichment, source_dir, "figureS2c_legacy_enrichment.tsv")
write_source(legacy_null, source_dir, "figureS2d_legacy_matched_null.tsv")
write_source(historical_12, source_dir, "figureS2e_legacy_12.tsv")
write_source(legacy_membership, source_dir, "figureS2f_legacy_membership.tsv")

# Supplementary Figure 3: fine-state resolution, model and composition sensitivities.
p3a <- ggplot() +
  annotate("rect", xmin = 0.1, xmax = 1.55, ymin = 0.35, ymax = 1.65, fill = figure_colours[["pale_blue"]], colour = figure_colours[["normal"]], linewidth = 0.5) +
  annotate("text", x = 0.825, y = 1.30, label = "Normal discovery cells", family = figure_font, fontface = "bold", size = 2.2) +
  annotate("text", x = 0.825, y = 0.88, label = "Programme and nuisance\ngenes excluded", family = figure_font, size = 1.7, lineheight = 0.95, colour = figure_colours[["muted"]]) +
  annotate("segment", x = 1.62, xend = 2.05, y = 1, yend = 1, arrow = grid::arrow(type = "closed", length = grid::unit(1.5, "mm")), colour = figure_colours[["grey"]], linewidth = 0.55) +
  annotate("rect", xmin = 2.12, xmax = 3.55, ymin = 0.35, ymax = 1.65, fill = figure_colours[["pale_gold"]], colour = figure_colours[["gold"]], linewidth = 0.5) +
  annotate("text", x = 2.835, y = 1.28, label = "Frozen reference", family = figure_font, fontface = "bold", size = 2.2) +
  annotate("text", x = 2.835, y = 0.88, label = "1,500 genes · 20 PCs\nk = 3, 4 or 5", family = figure_font, size = 1.7, lineheight = 0.95, colour = figure_colours[["muted"]]) +
  annotate("segment", x = 3.62, xend = 4.05, y = 1, yend = 1, arrow = grid::arrow(type = "closed", length = grid::unit(1.5, "mm")), colour = figure_colours[["grey"]], linewidth = 0.55) +
  annotate("rect", xmin = 4.12, xmax = 5.75, ymin = 0.35, ymax = 1.65, fill = figure_colours[["pale_orange"]], colour = figure_colours[["adenoma"]], linewidth = 0.5) +
  annotate("text", x = 4.935, y = 1.28, label = "Projection and test", family = figure_font, fontface = "bold", size = 2.2) +
  annotate("text", x = 4.935, y = 0.88, label = "All cells projected\nwithout refitting", family = figure_font, size = 1.7, lineheight = 0.95, colour = figure_colours[["muted"]]) +
  coord_cartesian(xlim = c(0, 5.9), ylim = c(0.2, 1.8), clip = "off") +
  theme_void(base_family = figure_font)

fine_diagnostics <- read_tsv(file.path(revision_root, "fine_states", "fine_state_model_diagnostics.tsv")) %>%
  mutate(broad_state = factor(broad_state, levels = c("ABS", "GOB", "TAC")))
p3b <- ggplot(fine_diagnostics, aes(k, silhouette, colour = broad_state)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = c(ABS = figure_colours[["adenoma"]], GOB = figure_colours[["normal"]], TAC = figure_colours[["purple"]])) +
  scale_x_continuous(breaks = c(3, 4, 5)) +
  labs(x = "Fine states per broad state", y = "Normal-reference silhouette") +
  theme_cb() +
  theme(legend.position = "top")

fine_effects <- read_tsv(file.path(revision_root, "fine_state_models", "fine_state_adjusted_route_effects.tsv")) %>%
  filter(k == 4) %>%
  mutate(
    scope = factor(scope, levels = rev(c("all", "ABS", "GOB", "TAC")), labels = rev(c("All states", "ABS", "GOB", "TAC"))),
    model = factor(model, levels = c("unweighted", "cell_count_weighted"), labels = c("Equal strata", "Cell-count weighted"))
  )
p3c <- ggplot(fine_effects, aes(estimate, scope, colour = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.6, position = position_dodge(width = 0.34)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.34)) +
  scale_colour_manual(values = c("Equal strata" = figure_colours[["adenoma"]], "Cell-count weighted" = figure_colours[["normal"]])) +
  scale_shape_manual(values = c("Equal strata" = 16, "Cell-count weighted" = 18)) +
  labs(x = "Fine-state-adjusted adenoma effect", y = NULL) +
  theme_cb() +
  theme(legend.position = "top")

decomposition <- read_tsv(file.path(revision_root, "fine_state_models", "programme_composition_decomposition.tsv"))
fraction <- decomposition %>%
  filter(decomposition_type == "fine_state", scope == "all", component %in% c("total", "within")) %>%
  select(partition, k, component, estimate) %>%
  pivot_wider(names_from = component, values_from = estimate) %>%
  mutate(
    within_fraction = within / total,
    partition = factor(partition, levels = c("discovery", "validation"), labels = c("Discovery", "Donor-disjoint validation"))
  )
p3d <- ggplot(fraction, aes(k, within_fraction, colour = partition, group = partition)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(Discovery = figure_colours[["grey"]], "Donor-disjoint validation" = figure_colours[["adenoma"]])) +
  scale_x_continuous(breaks = c(3, 4, 5)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.45, 0.9)) +
  labs(x = "Fine states per broad state", y = "Within-state fraction") +
  theme_cb() +
  theme(legend.position = "top")

decomp_primary <- decomposition %>%
  filter(k == 4, scope == "all", decomposition_type == "fine_state", component %in% c("composition", "within")) %>%
  mutate(
    component = factor(component, levels = c("composition", "within"), labels = c("Fine-state composition", "Within fine states")),
    partition = factor(partition, levels = c("discovery", "validation"), labels = c("Discovery", "Donor-disjoint validation"))
  )
p3e <- ggplot(decomp_primary, aes(estimate, component, colour = partition, shape = partition)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.6, position = position_dodge(width = 0.34)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.34)) +
  scale_colour_manual(values = c(Discovery = figure_colours[["grey"]], "Donor-disjoint validation" = figure_colours[["adenoma"]])) +
  scale_shape_manual(values = c(Discovery = 16, "Donor-disjoint validation" = 18)) +
  labs(x = "Contribution to adenoma–normal difference", y = NULL) +
  theme_cb() +
  theme(legend.position = "top")

purity <- read_tsv(file.path(state_root, "purity_adjusted_sensitivity", "primary_vs_purity_adjusted.tsv.gz"))
p3f <- ggplot(purity, aes(logFC_primary, logFC_purity_adjusted, colour = cell_type)) +
  geom_hline(yintercept = 0, colour = figure_colours[["line"]], linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = figure_colours[["line"]], linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = figure_colours[["grey"]], linewidth = 0.4) +
  geom_point(size = 0.55, alpha = 0.30) +
  facet_wrap(~cell_type) +
  scale_colour_manual(values = c(ABS = figure_colours[["adenoma"]], GOB = figure_colours[["normal"]], TAC = figure_colours[["purple"]])) +
  guides(colour = "none") +
  labs(x = "Primary gene effect", y = "Purity-adjusted gene effect") +
  theme_cb()

dslab <- read_tsv(file.path(state_root, "dslab_cnv_validation", "dslab_state_shared_cnv_composition_decomposition.tsv")) %>%
  filter(population == "conventional_adenomas", score == "state_shared_1843")
p3g <- ggplot(dslab, aes(observed_mean, common_composition_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = figure_colours[["grey"]], linewidth = 0.45) +
  geom_point(size = 2.3, colour = figure_colours[["adenoma"]]) +
  geom_text(aes(label = patient_token), nudge_y = 0.0012, family = figure_font, size = 1.55, colour = figure_colours[["muted"]]) +
  labs(x = "Observed patient score", y = "Common-composition score") +
  theme_cb()

dslab_stats <- read_tsv(file.path(state_root, "dslab_cnv_validation", "dslab_state_shared_cnv_validation_statistics.tsv")) %>%
  filter(population == "conventional_adenomas", score == "state_shared_1843",
         metric == "composition_only_range_over_observed_range")
p3h <- ggplot(dslab_stats, aes(estimate, "Full programme")) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.75, colour = figure_colours[["adenoma"]]) +
  geom_point(size = 2.4, colour = figure_colours[["adenoma"]]) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.16)) +
  labs(x = "Score range generated by CNV composition", y = NULL) +
  theme_cb()

fig_s3 <- p3a / (p3b | p3c) / (p3d | p3e) / (p3f | p3g | p3h) +
  plot_layout(heights = c(0.48, 1, 1, 1)) + plot_annotation(tag_levels = "a")
save_supplement(fig_s3, 3, "fine_state_and_composition_sensitivities", 240)
write_source(fine_diagnostics, source_dir, "figureS3b_fine_state_diagnostics.tsv")
write_source(fine_effects, source_dir, "figureS3c_weighted_models.tsv")
write_source(fraction, source_dir, "figureS3d_within_fraction.tsv")
write_source(decomp_primary, source_dir, "figureS3e_decomposition.tsv")
write_source(purity, source_dir, "figureS3f_purity_adjustment.tsv")
write_source(dslab, source_dir, "figureS3g_dslab_common_composition.tsv")
write_source(dslab_stats, source_dir, "figureS3h_dslab_range_fraction.tsv")

# Supplementary Figure 4: ranked pathway, regulator and anchor-gene structure.
pathways <- read_tsv(file.path(state_root, "interpretation", "pathway_competitive_enrichment_all.tsv.gz"))
clean_label <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}
hallmark <- pathways %>%
  filter(grepl("HALLMARK", gene_set), FDR <= 0.10) %>%
  arrange(FDR) %>%
  slice_head(n = 12) %>%
  mutate(
    label = factor(clean_label(gene_set), levels = rev(clean_label(gene_set))),
    evidence = ifelse(Direction == "Up", 1, -1) * pmin(-log10(FDR), 20)
  )
p4a <- ggplot(hallmark, aes(evidence, label, colour = Direction)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = evidence, yend = label), linewidth = 0.75) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = c(Down = figure_colours[["normal"]], Up = figure_colours[["adenoma"]])) +
  labs(x = "Down ← signed −log10 FDR → Up", y = NULL) +
  guides(colour = "none") +
  theme_cb()

functional <- pathways %>%
  filter(!grepl("HALLMARK", gene_set), FDR <= 0.05, interpretive_family != "other") %>%
  group_by(interpretive_family) %>%
  slice_min(FDR, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(FDR) %>%
  slice_head(n = 12) %>%
  mutate(
    label_text = stringr::str_wrap(clean_label(gene_set), width = 30),
    label = factor(label_text, levels = rev(label_text)),
    evidence = ifelse(Direction == "Up", 1, -1) * pmin(-log10(FDR), 20)
  )
p4b <- ggplot(functional, aes(evidence, label, colour = Direction)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = evidence, yend = label), linewidth = 0.75) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = c(Down = figure_colours[["normal"]], Up = figure_colours[["adenoma"]])) +
  labs(x = "Down ← signed −log10 FDR → Up", y = NULL) +
  guides(colour = "none") +
  theme_cb()

regulators <- read_tsv(file.path(state_root, "interpretation", "collectri_regulator_activity.tsv")) %>%
  arrange(q_value) %>%
  mutate(source = factor(source, levels = rev(source)))
p4c <- ggplot(regulators, aes(score, source)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = score, yend = source, colour = score > 0), linewidth = 0.75) +
  geom_point(aes(colour = score > 0), size = 2.2) +
  scale_colour_manual(values = c(`FALSE` = figure_colours[["normal"]], `TRUE` = figure_colours[["adenoma"]])) +
  guides(colour = "none") +
  labs(x = "CollecTRI activity statistic", y = NULL) +
  theme_cb()

common <- read_tsv(file.path(state_root, "common_effects", "cross_state_common_effects.tsv.gz"))
anchor_order <- c("ASCL2", "AXIN2", "RNF43", "ZNRF3", "EPHB2", "OLFM4", "CA2", "FABP1", "COX6C", "ACAA2")
anchors <- common %>%
  filter(gene %in% anchor_order) %>%
  select(gene, posterior_mean_ABS, posterior_mean_GOB, posterior_mean_TAC) %>%
  pivot_longer(-gene, names_to = "state", values_to = "effect") %>%
  mutate(
    state = sub("posterior_mean_", "", state),
    gene = factor(gene, levels = rev(anchor_order)),
    state = factor(state, levels = c("ABS", "GOB", "TAC"))
  )
anchor_limit <- max(abs(anchors$effect), na.rm = TRUE)
p4d <- ggplot(anchors, aes(state, gene, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.1f", effect)), family = figure_font, size = 2.0) +
  scale_fill_gradient2(low = figure_colours[["normal"]], mid = "white", high = figure_colours[["adenoma"]],
                       midpoint = 0, limits = c(-anchor_limit, anchor_limit), name = "Posterior effect") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = figure_font, base_size = 7) +
  theme(panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2, "mm"))

fig_s4 <- (p4a | p4b) / (p4c | p4d) + plot_annotation(tag_levels = "a")
save_supplement(fig_s4, 4, "functional_and_regulatory_structure", 175)
write_source(hallmark, source_dir, "figureS4a_hallmark.tsv")
write_source(functional, source_dir, "figureS4b_functional_families.tsv")
write_source(regulators, source_dir, "figureS4c_collectri.tsv")
write_source(anchors, source_dir, "figureS4d_anchor_genes.tsv")

cat("Supplementary Figures S1–S4 exported to ", out_dir, "\n", sep = "")
