#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1", "external_validation")
revision_root <- file.path(root, "results", "state_shared_revision_v2", "external_meta")
out_dir <- file.path(root, "figures", "communications_biology_v2.1")
source_dir <- file.path(out_dir, "source_data")

cohort_tests <- read_tsv(file.path(result_root, "external_cohort_tests.tsv")) %>%
  filter(signature_id == "state_shared_1843") %>%
  mutate(
    cohort = factor(cohort, levels = rev(c("GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820")))
  )
meta <- read_tsv(file.path(revision_root, "random_effects_meta_summary.tsv")) %>%
  filter(signature_id == "state_shared_1843")
leaveout <- read_tsv(file.path(revision_root, "random_effects_leave_one_cohort_out.tsv")) %>%
  filter(signature_id == "state_shared_1843") %>%
  mutate(excluded_cohort = factor(excluded_cohort, levels = levels(cohort_tests$cohort)))
boundary <- read_tsv(file.path(revision_root, "histology_and_grade_boundary_analysis.tsv"))
ffpe <- read_tsv(file.path(result_root, "ffpe_sample_scores.tsv.gz")) %>%
  filter(tissue_group %in% c("adenoma", "normal")) %>%
  group_by(patient_id, tissue_group, signature_id) %>%
  summarise(score = mean(programme_score), .groups = "drop")
coverage <- read_tsv(file.path(result_root, "external_gene_coverage.tsv")) %>%
  filter(signature_id == "state_shared_1843") %>%
  mutate(
    up_percent = 100 * up_coverage,
    down_percent = 100 * down_coverage,
    cohort = factor(cohort, levels = c("GSE8671", "GSE40362", "GSE72820", "GSE50114", "GSE41657"))
  )

# a. Cohort-specific estimates and prespecified random-effects synthesis.
forest <- cohort_tests %>%
  transmute(
    label = as.character(cohort),
    order = as.numeric(cohort) + 1,
    estimate = clustered_standardized_mean_difference,
    ci_low = clustered_standardized_ci_low,
    ci_high = clustered_standardized_ci_high,
    type = "Cohort"
  ) %>%
  bind_rows(data.frame(
    label = "Random-effects pooled",
    order = 1,
    estimate = meta$pooled_standardized_effect,
    ci_low = meta$ci_low,
    ci_high = meta$ci_high,
    type = "Pooled"
  )) %>%
  mutate(label = factor(label, levels = label[order(order)]))

p3a <- ggplot(forest, aes(estimate, label)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, colour = type), height = 0, linewidth = 0.65) +
  geom_point(aes(colour = type, shape = type), size = 2.25) +
  scale_colour_manual(values = c(Cohort = figure_colours[["adenoma"]], Pooled = figure_colours[["ink"]])) +
  scale_shape_manual(values = c(Cohort = 16, Pooled = 18)) +
  guides(colour = "none", shape = "none") +
  labs(x = "Adenoma effect (standardised mean difference)", y = NULL) +
  theme_cb()

# b. Leave-one-cohort-out random-effects estimates.
p3b <- ggplot(leaveout, aes(pooled_standardized_effect, excluded_cohort)) +
  geom_vline(xintercept = meta$pooled_standardized_effect, colour = figure_colours[["line"]], linewidth = 0.45) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65, colour = figure_colours[["normal"]]) +
  geom_point(size = 2.25, colour = figure_colours[["normal"]]) +
  labs(x = "Pooled effect after cohort omission", y = "Cohort omitted") +
  theme_cb()

# c. One-cohort histology boundary and grade contrasts.
selected_boundary <- boundary %>%
  filter(signature_id == "state_shared_1843") %>%
  filter(
    (cohort == "GSE40362" & contrast %in% c("adenoma_vs_normal", "hyperplastic_vs_normal", "adenoma_vs_hyperplastic")) |
      (cohort == "GSE41657" & contrast %in% c("low_grade_adenoma_vs_normal", "high_grade_adenoma_vs_normal", "high_vs_low_grade_adenoma"))
  ) %>%
  mutate(
    label = recode(
      contrast,
      adenoma_vs_normal = "Adenoma − normal",
      hyperplastic_vs_normal = "Hyperplastic − normal",
      adenoma_vs_hyperplastic = "Adenoma − hyperplastic",
      low_grade_adenoma_vs_normal = "Low grade − normal",
      high_grade_adenoma_vs_normal = "High grade − normal",
      high_vs_low_grade_adenoma = "High grade − low grade"
    ),
    cohort_label = recode(cohort, GSE40362 = "Histology", GSE41657 = "Grade"),
    label = factor(label, levels = rev(c(
      "Adenoma − normal", "Hyperplastic − normal", "Adenoma − hyperplastic",
      "Low grade − normal", "High grade − normal", "High grade − low grade"
    )))
  )

p3c <- ggplot(selected_boundary, aes(mean_difference, label, colour = cohort_label)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(Histology = figure_colours[["purple"]], Grade = figure_colours[["green"]])) +
  labs(x = "Full-programme score difference", y = NULL) +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

# d. Paired archival FFPE transfer for the full programme.
ffpe_full <- ffpe %>%
  filter(signature_id == "state_shared_1843") %>%
  tidyr::pivot_wider(names_from = tissue_group, values_from = score) %>%
  filter(!is.na(normal), !is.na(adenoma)) %>%
  mutate(patient_index = rank(adenoma - normal, ties.method = "first")) %>%
  tidyr::pivot_longer(c(normal, adenoma), names_to = "tissue_group", values_to = "score") %>%
  mutate(tissue_group = factor(tissue_group, levels = c("normal", "adenoma"), labels = c("Adjacent mucosa", "Adenoma")))

ffpe_summary <- ffpe_full %>%
  group_by(tissue_group) %>%
  summarise(median = median(score), .groups = "drop")

p3d <- ggplot(ffpe_full, aes(tissue_group, score, group = patient_id)) +
  geom_line(colour = figure_colours[["line"]], linewidth = 0.35, alpha = 0.75) +
  geom_point(aes(colour = tissue_group), size = 1.25, alpha = 0.85) +
  geom_point(data = ffpe_summary, aes(tissue_group, median, group = 1), shape = 23, size = 3.1, fill = "white", colour = figure_colours[["ink"]], stroke = 0.55) +
  scale_colour_manual(values = c("Adjacent mucosa" = figure_colours[["normal"]], Adenoma = figure_colours[["adenoma"]])) +
  guides(colour = "none") +
  annotate("text", x = 1.5, y = max(ffpe_full$score) + 0.12,
           label = "41/51 pairs increased", family = figure_font, size = 2.0,
           colour = figure_colours[["muted"]]) +
  labs(x = NULL, y = "Full-programme FFPE score") +
  theme_cb()

# e. Platform coverage of both full-programme arms.
coverage_long <- coverage %>%
  select(cohort, up_percent, down_percent) %>%
  tidyr::pivot_longer(c(up_percent, down_percent), names_to = "arm", values_to = "coverage") %>%
  mutate(arm = recode(arm, up_percent = "Adenoma-up", down_percent = "Adenoma-down"))

p3e <- ggplot(coverage_long, aes(cohort, coverage, fill = arm)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  scale_fill_manual(values = c("Adenoma-up" = figure_colours[["adenoma"]], "Adenoma-down" = figure_colours[["normal"]])) +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 50, 100), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Programme genes measured") +
  theme_cb() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "inside",
    legend.position.inside = c(0.04, 0.46),
    legend.justification = c(0, 0.5),
    legend.direction = "vertical",
    legend.key.height = unit(2.2, "mm"),
    legend.key.width = unit(3.0, "mm"),
    legend.background = element_rect(fill = scales::alpha("white", 0.86), colour = NA)
  )

figure <- (p3a | p3b) /
  (p3c | p3d | p3e) +
  plot_layout(heights = c(1, 1.05), widths = c(1.05, 1, 0.78)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

write_source(forest, source_dir, "figure3a_external_meta.tsv")
write_source(leaveout, source_dir, "figure3b_leave_one_cohort_out.tsv")
write_source(selected_boundary, source_dir, "figure3c_histology_grade.tsv")
write_source(ffpe_full, source_dir, "figure3d_ffpe_pairs.tsv")
write_source(coverage_long, source_dir, "figure3e_platform_coverage.tsv")

export_cb_figure(
  figure, out_dir,
  "figure3_external_recurrence_and_archival_transfer",
  width_mm = 178, height_mm = 155
)

cat("Figure 3 exported to ", out_dir, "\n", sep = "")
