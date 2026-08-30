#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")

root <- normalizePath(".", mustWork = TRUE)
validation_root <- file.path(root, "results", "state_aware_program_v1", "extended_validation")
rna_atac_root <- file.path(validation_root, "becker_rna_atac")
atlas_root <- file.path(validation_root, "crc_atlas")
out_dir <- file.path(root, "figures", "communications_biology_v2.0")
source_dir <- file.path(out_dir, "source_data")

paired <- read_tsv(file.path(rna_atac_root, "becker_rna_atac_paired_scores.tsv")) %>%
  filter(analysis_set == "normal_polyp") %>%
  mutate(tissue = ifelse(is_polyp == 1, "Polyp", "Normal / unaffected"))
cluster_models <- read_tsv(file.path(rna_atac_root, "becker_locked_rna_atac_patient_cluster_models.tsv"))
patient_correlations <- read_tsv(file.path(rna_atac_root, "becker_locked_rna_atac_patient_median_correlations.tsv"))
atlas <- read_tsv(file.path(atlas_root, "atlas_compact_8_donor_scores.tsv"))
atlas_models <- read_tsv(file.path(atlas_root, "atlas_locked_leave_one_study_out.tsv")) %>%
  mutate(estimable = tolower(as.character(estimable)) == "true")

# a. Matched sample-level RNA and WNT-related accessibility.
p4a <- ggplot(
  paired,
  aes(atac_tss__wnt_stem, rna_epi__ca_route_signature, colour = tissue)
) +
  geom_smooth(method = "lm", se = TRUE, colour = figure_colours[["grey"]], fill = figure_colours[["pale"]], linewidth = 0.55) +
  geom_point(size = 2.0, alpha = 0.9) +
  scale_colour_manual(values = c("Normal / unaffected" = figure_colours[["normal"]], Polyp = figure_colours[["adenoma"]])) +
  labs(x = "WNT/stem TSS accessibility", y = "Epithelial eight-gene RNA score") +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

# b. Patient-clustered, covariate-adjusted RNA–accessibility associations.
model_labels <- c(
  locked_route__wnt_tss = "WNT/stem TSS",
  locked_route__wnt_tss_minus_housekeeping = "WNT/stem TSS − housekeeping",
  locked_route__wnt_tcf_ascl2_axis = "WNT/TCF/ASCL2 axis",
  locked_route__wnt_tcf_ascl2_axis_minus_housekeeping = "WNT/TCF/ASCL2 − housekeeping"
)
model_plot <- cluster_models %>%
  filter(analysis_id %in% names(model_labels)) %>%
  mutate(
    label = factor(model_labels[analysis_id], levels = rev(unname(model_labels)))
  )

p4b <- ggplot(model_plot, aes(coef, label)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65, colour = figure_colours[["purple"]]) +
  geom_point(size = 2.25, colour = figure_colours[["purple"]]) +
  labs(x = "Adjusted coefficient (95% CI)", y = NULL) +
  theme_cb()

# c. Patient-median association, which preserves the patient as the displayed unit.
patient_median <- paired %>%
  group_by(patient_id) %>%
  summarise(
    atac_wnt_tss = median(atac_tss__wnt_stem),
    rna_score = median(rna_epi__ca_route_signature),
    .groups = "drop"
  )
median_result <- patient_correlations %>%
  filter(analysis_id == "locked_route__wnt_tss")

p4c <- ggplot(patient_median, aes(atac_wnt_tss, rna_score)) +
  geom_smooth(method = "lm", se = TRUE, colour = figure_colours[["grey"]], fill = figure_colours[["pale"]], linewidth = 0.55) +
  geom_point(size = 2.2, colour = figure_colours[["adenoma"]]) +
  annotate(
    "text",
    x = min(patient_median$atac_wnt_tss, na.rm = TRUE) + 0.05,
    y = max(patient_median$rna_score, na.rm = TRUE) - 0.02,
    label = sprintf("Spearman ρ = %.2f\npermutation P = %.3f", median_result$spearman_rho, median_result$p_value_patient_permutation),
    hjust = 0, vjust = 1, family = figure_font, size = 2.0,
    colour = figure_colours[["muted"]]
  ) +
  labs(x = "Patient-median WNT/stem accessibility", y = "Patient-median RNA score") +
  theme_cb()

# d. Independent epithelial-atlas distribution across pathological contexts.
carrier_labels <- c(
  normal_epithelial = "Normal\nepithelium",
  polyp_epithelial = "Polyp\nepithelium",
  polyp_cancer = "Polyp\ncancer-like",
  tumor_epithelial = "Primary\nepithelium",
  tumor_cancer = "Primary\ncancer-like",
  metastasis_epithelial = "Metastasis\nepithelium",
  metastasis_cancer = "Metastasis\ncancer-like"
)
atlas_plot <- atlas %>%
  filter(carrier_group %in% names(carrier_labels)) %>%
  mutate(
    carrier = factor(carrier_labels[carrier_group], levels = unname(carrier_labels))
  )
atlas_palette <- c(
  "Normal\nepithelium" = figure_colours[["normal"]],
  "Polyp\nepithelium" = figure_colours[["gold"]],
  "Polyp\ncancer-like" = figure_colours[["adenoma"]],
  "Primary\nepithelium" = figure_colours[["purple"]],
  "Primary\ncancer-like" = figure_colours[["green"]],
  "Metastasis\nepithelium" = "#517C73",
  "Metastasis\ncancer-like" = "#315E55"
)

p4d <- ggplot(atlas_plot, aes(carrier, score__ca_route_signature, colour = carrier)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.45, fill = "white") +
  geom_jitter(width = 0.14, height = 0, size = 0.65, alpha = 0.32) +
  scale_colour_manual(values = atlas_palette) +
  guides(colour = "none") +
  labs(x = NULL, y = "Donor–study carrier score") +
  theme_cb() +
  theme(axis.text.x = element_text(size = 5.2, angle = 35, hjust = 1))

# e. Atlas effects adjusted for source study.
atlas_base <- atlas_models %>%
  filter(
    outcome == "score__ca_route_signature",
    omitted_study == "__NONE__", estimable
  ) %>%
  mutate(
    state_label = recode(
      state,
      polyp_epithelial = "Polyp epithelium",
      polyp_cancer = "Polyp cancer-like",
      tumor_epithelial = "Primary epithelium",
      tumor_cancer = "Primary cancer-like",
      metastasis_epithelial = "Metastasis epithelium",
      metastasis_cancer = "Metastasis cancer-like"
    ),
    state_label = factor(state_label, levels = rev(c(
      "Polyp epithelium", "Polyp cancer-like", "Primary epithelium",
      "Primary cancer-like", "Metastasis epithelium", "Metastasis cancer-like"
    )))
  )

p4e <- ggplot(atlas_base, aes(coef, state_label)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.65, colour = figure_colours[["adenoma"]]) +
  geom_point(size = 2.2, colour = figure_colours[["adenoma"]]) +
  labs(x = "Study-adjusted difference from normal epithelium", y = NULL) +
  theme_cb()

# f. Leave-one-study-out coefficient ranges for the two premalignant carriers.
atlas_omission <- atlas_models %>%
  filter(
    outcome == "score__ca_route_signature",
    omitted_study != "__NONE__", estimable,
    state %in% c("polyp_epithelial", "polyp_cancer")
  ) %>%
  group_by(state) %>%
  summarise(
    minimum = min(coef), maximum = max(coef), median = median(coef),
    n_omissions = n(), n_positive = sum(coef > 0), .groups = "drop"
  ) %>%
  mutate(
    state_label = factor(
      recode(state, polyp_epithelial = "Polyp epithelium", polyp_cancer = "Polyp cancer-like"),
      levels = c("Polyp cancer-like", "Polyp epithelium")
    )
  )

p4f <- ggplot(atlas_omission, aes(median, state_label)) +
  geom_vline(xintercept = 0, colour = figure_colours[["muted"]], linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = minimum, xmax = maximum), height = 0, linewidth = 0.75, colour = figure_colours[["purple"]]) +
  geom_point(size = 2.3, colour = figure_colours[["purple"]]) +
  geom_text(aes(label = paste0(n_positive, "/", n_omissions, " positive")), nudge_y = -0.23, family = figure_font, size = 1.8, colour = figure_colours[["muted"]]) +
  scale_x_continuous(limits = c(0, 1.0)) +
  labs(x = "Coefficient range after omission", y = NULL) +
  theme_cb()

top <- p4a | p4b | p4c
bottom <- (p4d | p4e | p4f) + plot_layout(widths = c(1.25, 1, 0.9))
figure <- top / bottom +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

write_source(paired, source_dir, "figure4a_matched_rna_atac.tsv")
write_source(model_plot, source_dir, "figure4b_adjusted_rna_atac_models.tsv")
write_source(patient_median, source_dir, "figure4c_patient_medians.tsv")
write_source(atlas_plot, source_dir, "figure4d_atlas_carriers.tsv")
write_source(atlas_base, source_dir, "figure4e_atlas_adjusted_effects.tsv")
write_source(atlas_omission, source_dir, "figure4f_atlas_study_omission.tsv")

export_cb_figure(
  figure, out_dir,
  "figure4_regulatory_and_epithelial_context",
  width_mm = 178, height_mm = 170
)

cat("Figure 4 exported to ", out_dir, "\n", sep = "")
