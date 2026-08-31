#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")
suppressPackageStartupMessages(library(tidyr))

root <- normalizePath(".", mustWork = TRUE)
state_root <- file.path(root, "results", "state_aware_program_v1")
revision_root <- file.path(root, "results", "state_shared_revision_v2")
full_extended_root <- file.path(state_root, "extended_validation_full_programme")
out_dir <- Sys.getenv(
  "CB_OUT_DIR",
  unset = file.path(root, "figures", "communications_biology_v2.1")
)
source_dir <- file.path(out_dir, "source_data")
supplement_number_offset <- as.integer(Sys.getenv("CB_SUPPLEMENT_NUMBER_OFFSET", unset = "0"))

base_write_source <- write_source
write_source <- function(data, out_dir, filename) {
  match <- regexec("^figureS([0-9]+)(.*)$", filename)
  parts <- regmatches(filename, match)[[1]]
  if (length(parts) == 3 && supplement_number_offset != 0) {
    filename <- paste0(
      "figureS", as.integer(parts[[2]]) + supplement_number_offset, parts[[3]]
    )
  }
  base_write_source(data, out_dir, filename)
}

stem_for_version <- function(stem) {
  if (supplement_number_offset == 1) {
    return(switch(
      stem,
      compact_derivation_and_benchmarks = "reduced_readout_audit",
      empirical_perturbation_context = "complete_perturbation_context",
      stem
    ))
  }
  stem
}

save_supplement <- function(plot, number, stem, height_mm = 190) {
  export_cb_figure(
    tagged(plot), out_dir,
    paste0(
      "figureS", number + supplement_number_offset, "_", stem_for_version(stem)
    ),
    width_mm = 178, height_mm = height_mm
  )
}

signature_labels <- c(state_shared_1843 = "State-shared response", compact_8 = "Eight-gene candidate")
signature_colours <- c("State-shared response" = figure_colours[["grey"]], "Eight-gene candidate" = figure_colours[["adenoma"]])

# Supplementary Figure 5: external models, platform transfer and FFPE gene-level audit.
external_tests <- read_tsv(file.path(state_root, "external_validation", "external_cohort_tests.tsv")) %>%
  filter(signature_id %in% names(signature_labels)) %>%
  mutate(
    readout = factor(signature_labels[signature_id], levels = unname(signature_labels)),
    cohort = factor(cohort, levels = rev(c("GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820")))
  )
p5a <- ggplot(external_tests, aes(clustered_standardized_mean_difference, cohort, colour = readout, shape = readout)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = clustered_standardized_ci_low, xmax = clustered_standardized_ci_high),
                 height = 0, linewidth = 0.6, position = position_dodge(width = 0.35)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = signature_colours) +
  scale_shape_manual(values = c("State-shared response" = 16, "Eight-gene candidate" = 18)) +
  labs(x = "Adenoma effect (standardised mean difference)", y = NULL) +
  theme_cb() +
  theme(legend.position = "top")

pooled <- read_tsv(file.path(state_root, "external_validation", "external_pooled_models.tsv")) %>%
  filter(excluded_cohort == "__NONE__", signature_id %in% names(signature_labels)) %>%
  transmute(signature_id, model = "One-stage", estimate = adenoma_coef_sd, ci_low, ci_high)
adjusted <- read_tsv(file.path(state_root, "external_validation", "external_proliferation_adjusted_models.tsv")) %>%
  filter(excluded_cohort == "__NONE__", signature_id %in% names(signature_labels)) %>%
  transmute(signature_id, model = "Proliferation-adjusted", estimate = adenoma_coef_sd, ci_low, ci_high)
pooled_plot <- bind_rows(pooled, adjusted) %>%
  mutate(
    readout = factor(signature_labels[signature_id], levels = rev(unname(signature_labels))),
    model = factor(model, levels = c("One-stage", "Proliferation-adjusted"))
  )
p5b <- ggplot(pooled_plot, aes(estimate, readout, colour = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65, position = position_dodge(width = 0.35)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = c("One-stage" = figure_colours[["grey"]], "Proliferation-adjusted" = figure_colours[["normal"]])) +
  scale_shape_manual(values = c("One-stage" = 16, "Proliferation-adjusted" = 18)) +
  labs(x = "Pooled adenoma coefficient (SD)", y = NULL) +
  theme_cb() +
  theme(legend.position = "top")

leaveout <- read_tsv(file.path(state_root, "external_validation", "external_adjusted_leave_one_cohort_out.tsv")) %>%
  filter(signature_id %in% names(signature_labels)) %>%
  mutate(
    readout = factor(signature_labels[signature_id], levels = unname(signature_labels)),
    excluded_cohort = factor(excluded_cohort, levels = rev(c("GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820")))
  )
p5c <- ggplot(leaveout, aes(adenoma_coef_sd, excluded_cohort, colour = readout, shape = readout)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.58, position = position_dodge(width = 0.34)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.34)) +
  scale_colour_manual(values = signature_colours) +
  scale_shape_manual(values = c("State-shared response" = 16, "Eight-gene candidate" = 18)) +
  labs(x = "Proliferation-adjusted effect", y = "Cohort omitted") +
  theme_cb() +
  theme(legend.position = "top")

ffpe_scores <- read_tsv(file.path(state_root, "external_validation", "ffpe_sample_scores.tsv.gz")) %>%
  filter(tissue_group %in% c("normal", "adenoma"), signature_id %in% names(signature_labels)) %>%
  group_by(patient_id, signature_id, tissue_group) %>%
  summarise(score = mean(programme_score), .groups = "drop") %>%
  pivot_wider(names_from = c(signature_id, tissue_group), values_from = score, names_glue = "{signature_id}_{tissue_group}") %>%
  filter(complete.cases(state_shared_1843_adenoma, state_shared_1843_normal,
                        compact_8_adenoma, compact_8_normal)) %>%
  mutate(
    full_delta = state_shared_1843_adenoma - state_shared_1843_normal,
    compact_delta = compact_8_adenoma - compact_8_normal
  )
ffpe_delta_rho <- cor(ffpe_scores$full_delta, ffpe_scores$compact_delta, method = "spearman", use = "complete.obs")
p5d <- ggplot(ffpe_scores, aes(full_delta, compact_delta)) +
  geom_hline(yintercept = 0, colour = figure_colours[["line"]], linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = figure_colours[["line"]], linewidth = 0.3) +
  geom_smooth(method = "lm", se = TRUE, colour = figure_colours[["grey"]], fill = figure_colours[["pale"]], linewidth = 0.55) +
  geom_point(size = 1.7, colour = figure_colours[["adenoma"]], alpha = 0.75) +
  annotate("text", x = -Inf, y = Inf, label = sprintf("ρ = %.3f", ffpe_delta_rho), hjust = -0.08, vjust = 1.15,
           family = figure_font, size = 2.0, colour = figure_colours[["muted"]]) +
  labs(x = "State-shared response paired difference", y = "Eight-gene paired difference") +
  theme_cb()

ffpe_genes <- read_tsv(file.path(state_root, "external_validation", "ffpe_compact_gene_tests.tsv")) %>%
  mutate(
    arm = ifelse(expected_direction_x > 0, "Adenoma-up", "Adenoma-down"),
    gene = factor(gene, levels = rev(gene[order(mean_paired_delta)]))
  )
p5e <- ggplot(ffpe_genes, aes(mean_paired_delta, gene, colour = arm)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = mean_delta_ci_low, xmax = mean_delta_ci_high), height = 0, linewidth = 0.65) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = c("Adenoma-up" = figure_colours[["adenoma"]], "Adenoma-down" = figure_colours[["normal"]])) +
  labs(x = "Paired FFPE gene difference", y = NULL) +
  theme_cb() +
  theme(legend.position = "top")

external_fidelity <- read_tsv(file.path(state_root, "external_validation", "external_compact_full_fidelity.tsv")) %>%
  select(cohort, n_samples, spearman_compact_vs_full) %>%
  bind_rows(
    read_tsv(file.path(state_root, "external_validation", "ffpe_compact_full_fidelity.tsv")) %>%
      transmute(cohort = "GSE117606 FFPE", n_samples, spearman_compact_vs_full)
  ) %>%
  mutate(cohort = factor(cohort, levels = rev(c("GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820", "GSE117606 FFPE"))))
p5f <- ggplot(external_fidelity, aes(spearman_compact_vs_full, cohort)) +
  geom_segment(aes(x = 0, xend = spearman_compact_vs_full, yend = cohort), colour = figure_colours[["line"]], linewidth = 0.7) +
  geom_point(size = 2.2, colour = figure_colours[["adenoma"]]) +
  geom_text(aes(label = paste0("n=", n_samples)), nudge_y = -0.22, family = figure_font, size = 1.7, colour = figure_colours[["muted"]]) +
  scale_x_continuous(limits = c(0, 1.03), breaks = c(0, 0.5, 1.0)) +
  labs(x = "Candidate-to-full Spearman correlation", y = NULL) +
  theme_cb()

fig_s5 <- (p5a | p5b) / (p5c | p5d) / (p5e | p5f) + plot_annotation(tag_levels = "a")
save_supplement(fig_s5, 4, "external_and_ffpe_sensitivities", 210)
write_source(external_tests, source_dir, "figureS4a_external_readouts.tsv")
write_source(pooled_plot, source_dir, "figureS4b_pooled_models.tsv")
write_source(leaveout, source_dir, "figureS4c_adjusted_leaveout.tsv")
write_source(ffpe_scores, source_dir, "figureS4d_ffpe_score_differences.tsv")
write_source(ffpe_genes, source_dir, "figureS4e_ffpe_genes.tsv")
write_source(external_fidelity, source_dir, "figureS4f_compact_full_fidelity.tsv")

# Supplementary Figure 6: complete compact-readout derivation and benchmark audit.
candidate_summary <- read_tsv(file.path(state_root, "panel_derivation", "candidate_universe_summary.tsv")) %>%
  pivot_longer(c(strict_state_shared_genes, portable_protein_coding), names_to = "stage", values_to = "genes") %>%
  mutate(
    arm = factor(arm, levels = c("down", "up"), labels = c("Adenoma-down", "Adenoma-up")),
    stage = factor(stage, levels = c("strict_state_shared_genes", "portable_protein_coding"),
                   labels = c("State-shared", "Platform-measurable protein-coding"))
  )
p6a <- ggplot(candidate_summary, aes(stage, genes, fill = arm)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = genes), position = position_dodge(width = 0.72), vjust = -0.35, family = figure_font, size = 1.9) +
  scale_fill_manual(values = c("Adenoma-down" = figure_colours[["normal"]], "Adenoma-up" = figure_colours[["adenoma"]])) +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Genes (log scale)") +
  theme_cb() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "top")

state_fidelity <- read_tsv(file.path(state_root, "panel_derivation", "discovery_state_specific_fidelity_curve.tsv")) %>%
  mutate(cell_type = factor(cell_type, levels = c("ABS", "GOB", "TAC")))
p6b <- ggplot(state_fidelity, aes(total_genes, oof_spearman, colour = cell_type)) +
  geom_vline(xintercept = 8, linetype = 2, colour = figure_colours[["gold"]], linewidth = 0.5) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.45) +
  scale_colour_manual(values = c(ABS = figure_colours[["adenoma"]], GOB = figure_colours[["normal"]], TAC = figure_colours[["purple"]])) +
  labs(x = "Genes in balanced candidate", y = "State-specific out-of-fold fidelity") +
  theme_cb() +
  theme(legend.position = "top")

selected_genes <- c("EPHB2", "REG1A", "LTBP1", "RNF43", "CALM2", "COX6C", "B2M", "ACAA2")
selection_stability <- read_tsv(file.path(state_root, "panel_derivation", "candidate_selection_stability.tsv")) %>%
  mutate(
    selection_frequency_donor_bootstrap = replace_na(selection_frequency_donor_bootstrap, 0),
    selected = gene %in% selected_genes
  ) %>%
  arrange(desc(selection_frequency_donor_bootstrap)) %>%
  slice_head(n = 20) %>%
  bind_rows(
    read_tsv(file.path(state_root, "panel_derivation", "candidate_selection_stability.tsv")) %>%
      mutate(selection_frequency_donor_bootstrap = replace_na(selection_frequency_donor_bootstrap, 0), selected = gene %in% selected_genes) %>%
      filter(selected)
  ) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(selection_frequency_donor_bootstrap) %>%
  mutate(gene = factor(gene, levels = gene))
p6c <- ggplot(selection_stability, aes(selection_frequency_donor_bootstrap, gene, colour = selected)) +
  geom_segment(aes(x = 0, xend = selection_frequency_donor_bootstrap, yend = gene), colour = figure_colours[["line"]], linewidth = 0.6) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c(`FALSE` = figure_colours[["grey"]], `TRUE` = figure_colours[["adenoma"]])) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.03)) +
  guides(colour = "none") +
  labs(x = "Whole-donor bootstrap selection frequency", y = NULL) +
  theme_cb()

heldout_genes <- read_tsv(file.path(state_root, "heldout_validation", "heldout_compact_panel_gene_validation.tsv")) %>%
  select(gene, arm, logFC_ABS, logFC_GOB, logFC_TAC) %>%
  pivot_longer(starts_with("logFC_"), names_to = "state", values_to = "effect") %>%
  mutate(
    state = sub("logFC_", "", state),
    gene = factor(gene, levels = rev(selected_genes)),
    state = factor(state, levels = c("ABS", "GOB", "TAC"))
  )
heldout_limit <- max(abs(heldout_genes$effect), na.rm = TRUE)
p6d <- ggplot(heldout_genes, aes(state, gene, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = ifelse(is.na(effect), "NA", sprintf("%.1f", effect))), family = figure_font, size = 2.0) +
  scale_fill_gradient2(low = figure_colours[["normal"]], mid = "white", high = figure_colours[["adenoma"]], midpoint = 0,
                       limits = c(-heldout_limit, heldout_limit), na.value = figure_colours[["pale"]], name = "Validation effect") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = figure_font, base_size = 7) +
  theme(panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2, "mm"))

rank_effects <- read_tsv(file.path(revision_root, "compact_rank", "heldout_single_sample_rank_route_effects.tsv")) %>%
  mutate(scope = factor(scope, levels = rev(c("all", "ABS", "GOB", "TAC")), labels = rev(c("All states", "ABS", "GOB", "TAC"))))
p6e <- ggplot(rank_effects, aes(estimate, scope)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65, colour = figure_colours[["adenoma"]]) +
  geom_point(size = 2.2, colour = figure_colours[["adenoma"]]) +
  labs(x = "Single-sample eight-gene adenoma effect", y = NULL) +
  theme_cb()

internal_random <- read_tsv(file.path(revision_root, "compact_rank", "random_eight_gene_benchmark.tsv")) %>%
  transmute(partition = "Donor-disjoint", fidelity = spearman_with_full_programme)
external_random <- read_tsv(file.path(revision_root, "external_rank", "random_eight_gene_benchmark.tsv")) %>%
  transmute(partition = "External-cohort median", fidelity = median_cohort_spearman)
random_benchmark <- bind_rows(internal_random, external_random) %>%
  mutate(partition = factor(partition, levels = c("Donor-disjoint", "External-cohort median")))
internal_summary <- read_tsv(file.path(revision_root, "compact_rank", "random_eight_gene_benchmark_summary.tsv")) %>%
  transmute(partition = "Donor-disjoint", observed = observed_compact_spearman, q95 = random_q95)
external_summary <- read_tsv(file.path(revision_root, "external_rank", "random_eight_gene_benchmark_summary.tsv")) %>%
  transmute(partition = "External-cohort median", observed = observed_median_cohort_spearman, q95 = random_q95)
benchmark_summary <- bind_rows(internal_summary, external_summary) %>%
  mutate(partition = factor(partition, levels = levels(random_benchmark$partition)))
p6f <- ggplot(random_benchmark, aes(fidelity)) +
  geom_histogram(bins = 32, fill = figure_colours[["line"]], colour = "white", linewidth = 0.2) +
  geom_vline(data = benchmark_summary, aes(xintercept = observed), colour = figure_colours[["adenoma"]], linewidth = 0.65) +
  geom_vline(data = benchmark_summary, aes(xintercept = q95), colour = figure_colours[["normal"]], linetype = 2, linewidth = 0.65) +
  facet_wrap(~partition, scales = "free_y") +
  labs(x = "Direction-balanced random-panel fidelity", y = "Panels") +
  theme_cb()

fig_s6 <- (p6a | p6b) / (p6c | p6d) / (p6e | p6f) + plot_annotation(tag_levels = "a")
save_supplement(fig_s6, 5, "compact_derivation_and_benchmarks", 210)
write_source(candidate_summary, source_dir, "figureS5a_candidate_summary.tsv")
write_source(state_fidelity, source_dir, "figureS5b_state_fidelity.tsv")
write_source(selection_stability, source_dir, "figureS5c_selection_stability.tsv")
write_source(heldout_genes, source_dir, "figureS5d_heldout_gene_effects.tsv")
write_source(rank_effects, source_dir, "figureS5e_rank_effects.tsv")
write_source(random_benchmark, source_dir, "figureS5f_random_benchmarks.tsv")
write_source(benchmark_summary, source_dir, "figureS5f_benchmark_summary.tsv")

# Supplementary Figure 7: patient-aware multi-omic, atlas, spatial and protein context.
becker_adjusted <- read_tsv(file.path(full_extended_root, "becker", "becker_full_programme_adjusted_models.tsv")) %>%
  filter(term == "disease_stage_group_polyp", outcome %in% c("all__ca_route_signature", "epi__ca_route_signature")) %>%
  mutate(
    outcome_label = factor(recode(outcome, all__ca_route_signature = "All nuclei", epi__ca_route_signature = "Epithelial-marker positive"),
                           levels = rev(c("All nuclei", "Epithelial-marker positive"))),
    ci_low = coef - 1.96 * se_hc1,
    ci_high = coef + 1.96 * se_hc1
  )
p7a <- ggplot(becker_adjusted, aes(coef, outcome_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65, colour = figure_colours[["adenoma"]]) +
  geom_point(size = 2.2, colour = figure_colours[["adenoma"]]) +
  labs(x = "Adjusted polyp coefficient", y = NULL) +
  theme_cb()

rna_atac <- read_tsv(file.path(full_extended_root, "becker_rna_atac", "becker_locked_rna_atac_patient_cluster_models.tsv")) %>%
  mutate(
    label = factor(
      recode(analysis_id,
             locked_route__wnt_tss = "WNT/stem TSS",
             locked_route__wnt_tss_minus_housekeeping = "WNT/stem TSS − housekeeping",
             locked_route__wnt_tcf_ascl2_axis = "WNT/TCF/ASCL2 axis",
             locked_route__wnt_tcf_ascl2_axis_minus_housekeeping = "WNT/TCF/ASCL2 − housekeeping"),
      levels = rev(c("WNT/stem TSS", "WNT/stem TSS − housekeeping", "WNT/TCF/ASCL2 axis", "WNT/TCF/ASCL2 − housekeeping"))
    )
  )
p7b <- ggplot(rna_atac, aes(coef, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65, colour = figure_colours[["purple"]]) +
  geom_point(size = 2.2, colour = figure_colours[["purple"]]) +
  labs(x = "Patient-clustered adjusted coefficient", y = NULL) +
  theme_cb()

window <- read_tsv(file.path(full_extended_root, "becker_rna_atac", "becker_full_programme_regulatory_window_rna_correlations.tsv")) %>%
  filter(rna_feature == "rna_epi__ca_route_signature") %>%
  mutate(
    locus = factor(recode(locus_set, wnt_route_loci = "WNT/stem loci", wnt_tcf_ascl2_axis_loci = "WNT/TCF/ASCL2 loci")),
    distance = factor(distance_bin,
                      levels = c("tss_core_1kb", "promoter_proximal_2_5kb", "proximal_regulatory_10kb", "distal_flank_20kb"),
                      labels = c("TSS 1 kb", "2–5 kb", "10 kb", "20 kb"))
  )
p7c <- ggplot(window, aes(distance, spearman_rho, colour = locus, group = locus)) +
  geom_hline(yintercept = 0, colour = figure_colours[["line"]], linewidth = 0.35) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = c("WNT/stem loci" = figure_colours[["adenoma"]], "WNT/TCF/ASCL2 loci" = figure_colours[["normal"]])) +
  labs(x = "Regulatory distance window", y = "Spearman correlation with RNA") +
  theme_cb() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")

atlas_support <- read_tsv(file.path(full_extended_root, "crc_atlas", "atlas_locked_state_study_support.tsv")) %>%
  filter(carrier_group %in% c("polyp_epithelial", "polyp_cancer", "tumor_epithelial", "tumor_cancer")) %>%
  mutate(
    state = recode(carrier_group, polyp_epithelial = "Polyp epi.", polyp_cancer = "Polyp cancer-like",
                   tumor_epithelial = "Primary epi.", tumor_cancer = "Primary cancer-like"),
    study = gsub("_", " ", sub("_[A-Za-z]+$", "", study_id))
  )
p7d <- ggplot(atlas_support, aes(state, study, fill = n_donors)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = figure_colours[["pale_orange"]], high = figure_colours[["adenoma"]], name = "Donors") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = figure_font, base_size = 6.3) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), plot.margin = margin(2, 2, 2, 2, "mm"))

within_study <- read_tsv(file.path(full_extended_root, "crc_atlas", "atlas_locked_within_study_contrasts.tsv")) %>%
  filter(outcome == "score__ca_route_signature", state %in% c("polyp_epithelial", "polyp_cancer", "tumor_epithelial", "tumor_cancer")) %>%
  arrange(desc(n_target_donors)) %>%
  slice_head(n = 14) %>%
  mutate(
    label_text = paste(gsub("_", " ", sub("_[A-Za-z]+$", "", study_id)),
                       recode(state, polyp_epithelial = "polyp epi.", polyp_cancer = "polyp cancer-like",
                              tumor_epithelial = "primary epi.", tumor_cancer = "primary cancer-like"), sep = " · "),
    label = factor(label_text, levels = rev(label_text))
  )
p7e <- ggplot(within_study, aes(coef, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.58, colour = figure_colours[["purple"]]) +
  geom_point(size = 1.9, colour = figure_colours[["purple"]]) +
  labs(x = "Within-study difference from normal epithelium", y = NULL) +
  theme_cb(base_size = 6.4)

spatial <- read_tsv(file.path(full_extended_root, "perturbation_spatial", "spatial_full_programme_section_effects.tsv")) %>%
  filter(comparison == "tumor_vs_non_neoplastic_epithelium", feature %in% c("route_score", "route_residual_prolif_epithelial")) %>%
  mutate(feature = factor(feature, levels = c("route_score", "route_residual_prolif_epithelial"), labels = c("Raw", "Adjusted")))
p7f <- ggplot(spatial, aes(feature, difference, group = sample_id)) +
  geom_hline(yintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_line(colour = figure_colours[["line"]], linewidth = 0.5) +
  geom_point(aes(colour = feature), size = 2.0) +
  scale_colour_manual(values = c(Raw = figure_colours[["adenoma"]], Adjusted = figure_colours[["normal"]])) +
  guides(colour = "none") +
  labs(x = NULL, y = "Tumour − non-neoplastic epithelium") +
  theme_cb()

protein <- read_tsv(file.path(root, "figures", "communications_biology_v1.2", "source_data", "figure6g_protein_evidence_matrix.tsv")) %>%
  mutate(
    gene = factor(gene, levels = rev(c("OLFM4", "CA2", "FABP1"))),
    resource = factor(resource, levels = c("PXD002137 / Differential", "PXD000445 / Paired", "PXD017269 / FFPE", "PXD046999 / DVP"))
  )
p7g <- ggplot(protein, aes(resource, gene, fill = role)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = display), family = figure_font, size = 1.75, lineheight = 0.9) +
  scale_fill_manual(values = c("Directional difference" = figure_colours[["pale_orange"]],
                               "Direction only" = figure_colours[["pale_gold"]],
                               "Detectability" = figure_colours[["pale_blue"]],
                               "Not detected" = figure_colours[["pale"]])) +
  guides(fill = "none") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = figure_font, base_size = 6.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), plot.margin = margin(2, 2, 2, 2, "mm"))

fig_s7 <- (p7a | p7b) / (p7c | p7f) / (p7d | p7e) / p7g +
  plot_layout(heights = c(0.9, 0.9, 1.25, 0.7)) + plot_annotation(tag_levels = "a")
save_supplement(fig_s7, 6, "multiomic_atlas_spatial_protein", 245)
write_source(becker_adjusted, source_dir, "figureS6a_becker_adjusted.tsv")
write_source(rna_atac, source_dir, "figureS6b_rna_atac_models.tsv")
write_source(window, source_dir, "figureS6c_regulatory_windows.tsv")
write_source(atlas_support, source_dir, "figureS6d_atlas_support.tsv")
write_source(within_study, source_dir, "figureS6e_within_study.tsv")
write_source(spatial, source_dir, "figureS6f_spatial_sections.tsv")
write_source(protein, source_dir, "figureS6g_protein_context.tsv")

# Supplementary Figure 7: complete empirical perturbation context.
feature_coverage <- read_tsv(file.path(full_extended_root, "perturbation_spatial", "perturbation_full_programme_feature_coverage.tsv")) %>%
  filter(feature %in% c("route_up", "route_down", "wnt_stem", "proliferation_control", "epithelial_control")) %>%
  mutate(
    feature = factor(feature, levels = c("route_up", "route_down", "wnt_stem", "proliferation_control", "epithelial_control"),
                     labels = c("Response up", "Response down", "WNT/stem", "Proliferation", "Epithelial control")),
    dataset = factor(dataset, levels = rev(unique(dataset)))
  )
p8a <- ggplot(feature_coverage, aes(feature, dataset, fill = coverage_fraction)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = percent(coverage_fraction, accuracy = 1)), family = figure_font, size = 1.7) +
  scale_fill_gradient(low = figure_colours[["pale"]], high = figure_colours[["adenoma"]], limits = c(0, 1), name = "Coverage") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = figure_font, base_size = 6.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 28, hjust = 1), plot.margin = margin(2, 2, 2, 2, "mm"))

perturbation_summary <- read_tsv(file.path(full_extended_root, "perturbation_spatial", "perturbation_full_programme_effect_summary.tsv")) %>%
  filter(feature == "route_score") %>%
  mutate(
    label_text = paste(dataset, gsub("_", " ", comparison), sep = " · "),
    label = factor(label_text, levels = rev(label_text)),
    match = n_expected_direction == n_units
  )
perturbation_units <- read_tsv(file.path(full_extended_root, "perturbation_spatial", "perturbation_full_programme_unit_effects.tsv")) %>%
  filter(feature == "route_score") %>%
  mutate(
    label_text = paste(dataset, gsub("_", " ", comparison), sep = " · "),
    label = factor(label_text, levels = levels(perturbation_summary$label))
  )
p8b <- ggplot(perturbation_summary, aes(mean_difference, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = bootstrap_mean_ci_low, xmax = bootstrap_mean_ci_high), height = 0, linewidth = 0.58, colour = figure_colours[["grey"]]) +
  geom_point(data = perturbation_units, aes(difference, label), inherit.aes = FALSE,
             size = 1.1, alpha = 0.45, colour = figure_colours[["normal"]], position = position_jitter(height = 0.08, width = 0)) +
  geom_point(aes(colour = match), size = 2.0) +
  scale_colour_manual(values = c(`FALSE` = figure_colours[["gold"]], `TRUE` = figure_colours[["adenoma"]]), na.value = figure_colours[["grey"]]) +
  guides(colour = "none") +
  labs(x = "Change in state-shared response score", y = NULL) +
  theme_cb(base_size = 6.2)

direction_audit <- perturbation_summary %>%
  filter(is.finite(n_expected_direction), n_units > 0) %>%
  mutate(fraction_expected = n_expected_direction / n_units)
p8c <- ggplot(direction_audit, aes(fraction_expected, label)) +
  geom_segment(aes(x = 0, xend = fraction_expected, yend = label), colour = figure_colours[["line"]], linewidth = 0.65) +
  geom_point(size = 2.0, colour = figure_colours[["adenoma"]]) +
  geom_text(aes(label = paste0(n_expected_direction, "/", n_units)), nudge_y = -0.21, family = figure_font, size = 1.55, colour = figure_colours[["muted"]]) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.03)) +
  labs(x = "Independent units in prespecified direction", y = NULL) +
  theme_cb(base_size = 6.2)

apc <- read_tsv(file.path(state_root, "external_validation", "apc_organoid_effects.tsv")) %>%
  filter(
    signature_id == "state_shared_1843",
    comparison %in% c("WT_withdrawal", "APC_vs_WT_without_Wnt", "genotype_by_Wnt_interaction")
  ) %>%
  mutate(
    comparison = factor(recode(comparison,
                               WT_withdrawal = "WNT withdrawal in WT",
                               APC_vs_WT_without_Wnt = "APC-KO vs WT, WNT−",
                               genotype_by_Wnt_interaction = "Genotype × WNT"),
                        levels = rev(c("WNT withdrawal in WT", "APC-KO vs WT, WNT−", "Genotype × WNT")))
  )
p8d <- ggplot(apc, aes(mean_difference, comparison)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = min_difference, xmax = max_difference), height = 0, linewidth = 0.6, colour = figure_colours[["adenoma"]]) +
  geom_point(size = 2.0, colour = figure_colours[["adenoma"]]) +
  labs(x = "Mean donor difference (range; n = 3)", y = NULL) +
  theme_cb()

fig_s8 <- (p8a | p8d) / (p8b | p8c) + plot_annotation(tag_levels = "a")
save_supplement(fig_s8, 7, "empirical_perturbation_context", 175)
write_source(feature_coverage, source_dir, "figureS7a_perturbation_coverage.tsv")
write_source(apc, source_dir, "figureS7b_apc_organoids.tsv")
write_source(perturbation_summary, source_dir, "figureS7c_perturbation_summary.tsv")
write_source(perturbation_units, source_dir, "figureS7c_perturbation_units.tsv")
write_source(direction_audit, source_dir, "figureS7d_direction_audit.tsv")

cat("Supplementary Figures S4–S7 exported to ", out_dir, "\n", sep = "")
