#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")

root <- normalizePath(".", mustWork = TRUE)

# Reuse the audited RNA–ATAC and CRC Atlas panels.
source("analysis/build_state_shared_revision_figure4_v3.R")

out_dir <- file.path(root, "figures", "communications_biology_v5.0")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

p5a <- p4a +
  labs(
    x = "WNT/stem TSS accessibility",
    y = "State-shared remodelling score"
  )

p5b <- p4b +
  labs(x = "Patient-clustered coefficient (95% CI)", y = NULL)

p5c <- p4e +
  labs(x = "Study-adjusted difference from normal epithelium", y = NULL)

p5d <- p4f +
  labs(x = "Coefficient range after study omission", y = NULL)

spatial_path <- file.path(
  root, "results", "state_aware_program_v1",
  "extended_validation_full_programme", "perturbation_spatial",
  "spatial_full_programme_section_effects.tsv"
)
spatial <- read_tsv(spatial_path) %>%
  filter(
    comparison == "tumor_vs_non_neoplastic_epithelium",
    feature %in% c("route_score", "route_residual_prolif_epithelial")
  ) %>%
  mutate(
    analysis = factor(
      feature,
      levels = c("route_score", "route_residual_prolif_epithelial"),
      labels = c("Raw", "Adjusted")
    )
  )

p5e <- ggplot(spatial, aes(analysis, difference, group = sample_id)) +
  geom_hline(
    yintercept = 0, colour = figure_colours[["muted"]],
    linetype = 2, linewidth = 0.35
  ) +
  geom_line(colour = figure_colours[["line"]], linewidth = 0.45) +
  geom_point(
    aes(fill = analysis), shape = 21, size = 2.2,
    colour = "white", stroke = 0.35
  ) +
  scale_fill_manual(values = c(
    Raw = figure_colours[["grey"]],
    Adjusted = figure_colours[["adenoma"]]
  )) +
  annotate(
    "text", x = 1, y = max(spatial$difference) + 0.018,
    label = "5/6 positive", family = figure_font, size = 1.8,
    colour = figure_colours[["muted"]]
  ) +
  annotate(
    "text", x = 2, y = max(spatial$difference) + 0.018,
    label = "6/6 positive", family = figure_font, size = 1.8,
    colour = figure_colours[["adenoma"]]
  ) +
  guides(fill = "none") +
  labs(
    x = NULL,
    y = "Pathology-region score difference"
  ) +
  theme_cb()

protein_path <- file.path(
  root, "results", "public_adenoma_protein_triangulation",
  "candidate_public_protein_evidence_matrix.tsv"
)
protein <- read_tsv(protein_path) %>%
  filter(gene %in% c("OLFM4", "EPHB3", "CA2", "PRDX6")) %>%
  mutate(
    direction = factor(
      ifelse(expected_direction > 0, "Adenoma-up", "Adenoma-down"),
      levels = c("Adenoma-up", "Adenoma-down")
    ),
    gene = factor(gene, levels = rev(c("OLFM4", "EPHB3", "CA2", "PRDX6")))
  )

p5f <- ggplot(protein, aes(age_sex_adjusted_log2_effect, gene, colour = direction)) +
  geom_vline(
    xintercept = 0, colour = figure_colours[["muted"]],
    linetype = 2, linewidth = 0.35
  ) +
  geom_errorbarh(
    aes(xmin = age_sex_adjusted_ci_low, xmax = age_sex_adjusted_ci_high),
    height = 0, linewidth = 0.68
  ) +
  geom_point(size = 2.35) +
  scale_colour_manual(values = c(
    "Adenoma-up" = figure_colours[["adenoma"]],
    "Adenoma-down" = figure_colours[["normal"]]
  )) +
  labs(x = "Adjusted adenoma effect (log2 protein abundance)", y = NULL) +
  theme_cb() +
  theme(legend.position = "top", legend.justification = "left")

figure <- (p5a | p5b) / (p5c | p5d) / (p5e | p5f) +
  plot_layout(heights = c(1.04, 0.92, 0.90), widths = c(1.05, 0.95)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

write_source(paired, source_dir, "figure5a_matched_rna_atac.tsv")
write_source(model_plot, source_dir, "figure5b_patient_clustered_rna_atac.tsv")
write_source(atlas_base, source_dir, "figure5c_atlas_adjusted_effects.tsv")
write_source(atlas_omission, source_dir, "figure5d_atlas_leave_one_study_out.tsv")
write_source(spatial, source_dir, "figure5e_spatial_section_effects.tsv")
write_source(protein, source_dir, "figure5f_reciprocal_protein_effects.tsv")

export_cb_figure(
  figure, out_dir, "figure5_multimodal_tissue_context",
  width_mm = 178, height_mm = 188
)

message("Communications Biology v5 Figure 5 exported")
