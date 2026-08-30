#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")
suppressPackageStartupMessages(library(tidyr))

root <- normalizePath(".", mustWork = TRUE)
state_root <- file.path(root, "results", "state_aware_program_v1")
revision_root <- file.path(root, "results", "state_shared_revision_v2")
external_root <- file.path(state_root, "external_validation")
meta_root <- file.path(revision_root, "external_meta")
out_dir <- file.path(root, "figures", "communications_biology_v3.0")
source_dir <- file.path(out_dir, "source_data")

external_order <- c("GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820")

# a. Cohort-specific effects and the prespecified random-effects synthesis.
cohort_effects <- read_tsv(file.path(external_root, "external_cohort_tests.tsv")) %>%
  filter(signature_id == "state_shared_1843") %>%
  transmute(
    estimate_id = cohort,
    estimate = clustered_standardized_mean_difference,
    ci_low = clustered_standardized_ci_low,
    ci_high = clustered_standardized_ci_high,
    estimate_type = "Cohort"
  )

meta_complete <- read_tsv(file.path(meta_root, "random_effects_meta_summary.tsv")) %>%
  filter(signature_id == "state_shared_1843", excluded_cohort == "__NONE__") %>%
  transmute(
    estimate_id = "REML pooled",
    estimate = pooled_standardized_effect,
    ci_low,
    ci_high,
    estimate_type = "Pooled"
  )

forest_data <- bind_rows(cohort_effects, meta_complete) %>%
  mutate(
    estimate_id = factor(estimate_id, levels = rev(c(external_order, "REML pooled"))),
    estimate_type = factor(estimate_type, levels = c("Cohort", "Pooled"))
  )

p4a <- ggplot(forest_data, aes(estimate, estimate_id, colour = estimate_type, shape = estimate_type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65) +
  geom_point(size = 2.15) +
  scale_colour_manual(values = c("Cohort" = figure_colours[["grey"]], "Pooled" = figure_colours[["adenoma"]])) +
  scale_shape_manual(values = c("Cohort" = 16, "Pooled" = 18)) +
  guides(colour = "none", shape = "none") +
  labs(x = "Adenoma effect (standardised mean difference)", y = NULL) +
  theme_cb()

# b. Complete-data model sensitivity without introducing the reduced readout.
one_stage <- read_tsv(file.path(external_root, "external_pooled_models.tsv")) %>%
  filter(signature_id == "state_shared_1843", excluded_cohort == "__NONE__") %>%
  transmute(model = "One-stage cohort-adjusted", estimate = adenoma_coef_sd, ci_low, ci_high)

proliferation_adjusted <- read_tsv(file.path(external_root, "external_proliferation_adjusted_models.tsv")) %>%
  filter(signature_id == "state_shared_1843", excluded_cohort == "__NONE__") %>%
  transmute(model = "One-stage + proliferation", estimate = adenoma_coef_sd, ci_low, ci_high)

pooled_sensitivity <- bind_rows(
  meta_complete %>% transmute(model = "REML random effects", estimate, ci_low, ci_high),
  one_stage,
  proliferation_adjusted
) %>%
  mutate(model = factor(model, levels = rev(c(
    "REML random effects", "One-stage cohort-adjusted", "One-stage + proliferation"
  ))))

p4b <- ggplot(pooled_sensitivity, aes(estimate, model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.68,
                 colour = figure_colours[["normal"]]) +
  geom_point(size = 2.2, colour = figure_colours[["normal"]]) +
  labs(x = "Pooled adenoma effect (SD)", y = NULL) +
  theme_cb()

# c. Two complementary leave-one-cohort-out audits.
meta_leaveout <- read_tsv(file.path(meta_root, "random_effects_leave_one_cohort_out.tsv")) %>%
  filter(signature_id == "state_shared_1843") %>%
  transmute(
    excluded_cohort,
    model = "Random effects",
    estimate = pooled_standardized_effect,
    ci_low,
    ci_high
  )

adjusted_leaveout <- read_tsv(file.path(external_root, "external_adjusted_leave_one_cohort_out.tsv")) %>%
  filter(signature_id == "state_shared_1843") %>%
  transmute(
    excluded_cohort,
    model = "Proliferation-adjusted",
    estimate = adenoma_coef_sd,
    ci_low,
    ci_high
  )

leaveout_data <- bind_rows(meta_leaveout, adjusted_leaveout) %>%
  mutate(
    excluded_cohort = factor(excluded_cohort, levels = rev(external_order)),
    model = factor(model, levels = c("Random effects", "Proliferation-adjusted"))
  )

p4c <- ggplot(leaveout_data, aes(estimate, excluded_cohort, colour = model, shape = model)) +
  geom_vline(xintercept = meta_complete$estimate, linetype = 3,
             colour = figure_colours[["line"]], linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.58,
                 position = position_dodge(width = 0.38)) +
  geom_point(size = 1.95, position = position_dodge(width = 0.38)) +
  scale_colour_manual(values = c(
    "Random effects" = figure_colours[["grey"]],
    "Proliferation-adjusted" = figure_colours[["normal"]]
  )) +
  scale_shape_manual(values = c("Random effects" = 16, "Proliferation-adjusted" = 18)) +
  labs(x = "Adenoma effect after cohort omission", y = "Cohort omitted") +
  theme_cb() +
  theme(legend.position = "top")

# d. Patient-level FFPE paired differences for the complete programme.
ffpe_differences <- read_tsv(file.path(external_root, "ffpe_sample_scores.tsv.gz")) %>%
  filter(signature_id == "state_shared_1843", tissue_group %in% c("normal", "adenoma")) %>%
  group_by(patient_id, tissue_group) %>%
  summarise(score = mean(programme_score), .groups = "drop") %>%
  pivot_wider(names_from = tissue_group, values_from = score) %>%
  filter(complete.cases(normal, adenoma)) %>%
  mutate(delta = adenoma - normal) %>%
  arrange(delta) %>%
  mutate(
    patient_rank = row_number(),
    direction = factor(ifelse(delta > 0, "Adenoma higher", "Adenoma not higher"),
                       levels = c("Adenoma not higher", "Adenoma higher"))
  )

positive_pairs <- sum(ffpe_differences$delta > 0)
total_pairs <- nrow(ffpe_differences)

p4d <- ggplot(ffpe_differences, aes(patient_rank, delta, colour = direction)) +
  geom_hline(yintercept = 0, linetype = 2, colour = figure_colours[["muted"]], linewidth = 0.35) +
  geom_segment(aes(xend = patient_rank, y = 0, yend = delta), linewidth = 0.55) +
  geom_point(size = 1.55) +
  annotate(
    "text", x = Inf, y = Inf,
    label = sprintf("%d/%d pairs above zero", positive_pairs, total_pairs),
    hjust = 1.04, vjust = 1.2, family = figure_font, size = 2.0,
    colour = figure_colours[["muted"]]
  ) +
  scale_colour_manual(values = c(
    "Adenoma not higher" = figure_colours[["normal"]],
    "Adenoma higher" = figure_colours[["adenoma"]]
  )) +
  guides(colour = "none") +
  labs(x = "FFPE pair ranked by difference", y = "Adenoma − adjacent-mucosa score") +
  theme_cb()

# e. Arm-specific platform coverage, including archival FFPE tissue.
platform_coverage <- bind_rows(
  read_tsv(file.path(external_root, "external_gene_coverage.tsv")),
  read_tsv(file.path(external_root, "ffpe_gene_coverage.tsv"))
) %>%
  filter(signature_id == "state_shared_1843") %>%
  select(cohort, up_coverage, down_coverage) %>%
  pivot_longer(c(up_coverage, down_coverage), names_to = "arm", values_to = "coverage") %>%
  mutate(
    arm = factor(arm, levels = c("up_coverage", "down_coverage"),
                 labels = c("Adenoma-up", "Adenoma-down")),
    cohort = factor(cohort, levels = rev(c(external_order, "GSE117606"))),
    coverage_label = percent(coverage, accuracy = 1)
  )

p4e <- ggplot(platform_coverage, aes(arm, cohort, fill = coverage)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = coverage_label), family = figure_font, size = 2.1,
            colour = figure_colours[["ink"]]) +
  scale_fill_gradientn(
    colours = c(figure_colours[["pale"]], figure_colours[["pale_blue"]], figure_colours[["normal"]]),
    values = rescale(c(0, 0.5, 1)), limits = c(0, 1), labels = percent,
    name = "Coverage"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = figure_font, base_size = 7) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6.2),
    axis.text = element_text(colour = figure_colours[["ink"]]),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

figure_s4 <- ((p4a | p4b) / (p4c | p4d) / p4e) +
  plot_layout(heights = c(1, 1, 0.82)) +
  plot_annotation(tag_levels = "a")

export_cb_figure(
  tagged(figure_s4), out_dir,
  "figureS4_external_and_ffpe_sensitivities",
  width_mm = 178, height_mm = 205
)

write_source(forest_data, source_dir, "figureS4a_external_forest.tsv")
write_source(pooled_sensitivity, source_dir, "figureS4b_model_sensitivity.tsv")
write_source(leaveout_data, source_dir, "figureS4c_leave_one_cohort_out.tsv")
write_source(ffpe_differences, source_dir, "figureS4d_ffpe_paired_differences.tsv")
write_source(platform_coverage, source_dir, "figureS4e_platform_coverage.tsv")

stopifnot(total_pairs == 51L, positive_pairs == 41L)
message("Supplementary Figure 4 written with complete-programme evidence only.")
