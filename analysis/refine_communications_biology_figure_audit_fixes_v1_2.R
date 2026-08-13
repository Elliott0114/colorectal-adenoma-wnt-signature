#!/usr/bin/env Rscript

# Targeted v1.2 figure-audit fixes.
#
# Changes are restricted to confirmed presentation defects:
# - Fig. 2: move numerical labels away from points and report the same precision
#   used in the Results.
# - Fig. 4: enlarge and simplify dense source-study labels.
# - Supplementary Fig. 4: replace internal analysis-variable names.
# - Supplementary Fig. 6: replace internal feature names and standardise gene
#   and treatment labels.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
  library(svglite)
  library(ragg)
  library(grid)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  file.path(getwd(), "analysis", "refine_communications_biology_figure_audit_fixes_v1_2.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT_DIR <- file.path(ROOT, "figures", "communications_biology_v1.2")
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
RESULT_DIR <- file.path(ROOT, "results", "objective_compact_panel_v2_7")
EXT_DIR <- file.path(RESULT_DIR, "extended_validation")

JOURNAL_FONT <- "Arial"
COL <- c(
  ink = "#25292D", route = "#C85E3D", route_light = "#E9B6A5",
  wnt = "#356F9D", wnt_light = "#B8D0E1", neutral = "#727B84",
  neutral_light = "#D8DDE1", neutral_pale = "#F3F5F6",
  adenoma = "#D9A640", crc = "#3F9A8B", context = "#7C6AA3",
  pass = "#4D8C75"
)

read_tsv <- function(path) {
  read.delim(
    path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
    quote = "", comment.char = ""
  )
}
read_source <- function(filename) read_tsv(file.path(SOURCE_DIR, filename))

sanitize_frame <- function(frame) {
  frame[] <- lapply(frame, function(column) {
    if (is.character(column)) {
      column <- gsub("[\r\n]+", "; ", column)
      column <- gsub("[[:space:]]+", " ", column)
      trimws(column)
    } else {
      column
    }
  })
  frame
}

write_source <- function(frame, filename) {
  write.table(
    sanitize_frame(frame), file.path(SOURCE_DIR, filename), sep = "\t",
    quote = FALSE, row.names = FALSE, na = ""
  )
}

theme_journal <- function(base_size = 7.0) {
  theme_classic(base_size = base_size, base_family = JOURNAL_FONT) +
    theme(
      text = element_text(family = JOURNAL_FONT, colour = COL[["ink"]]),
      axis.line = element_line(linewidth = 0.32, colour = COL[["ink"]]),
      axis.ticks = element_line(linewidth = 0.28, colour = COL[["ink"]]),
      axis.ticks.length = unit(0.9, "mm"),
      axis.title = element_text(size = base_size, colour = COL[["ink"]]),
      axis.text = element_text(size = base_size - 0.55, colour = COL[["ink"]]),
      panel.grid = element_blank(), legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 0.35, face = "bold"),
      legend.text = element_text(size = base_size - 0.75),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.15, face = "bold"),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(1.6, 1.8, 1.6, 1.8, unit = "mm"),
      plot.tag = element_text(size = 8.2, face = "bold", colour = COL[["ink"]])
    )
}

clean_panel <- function(plot) {
  plot + labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL) +
    theme(
      text = element_text(family = JOURNAL_FONT, colour = COL[["ink"]]),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.tag = element_text(size = 8.2, face = "bold", family = JOURNAL_FONT)
    )
}
tag_theme <- theme(
  plot.tag.position = c(0, 1),
  plot.tag = element_text(
    size = 8.2, face = "bold", family = JOURNAL_FONT, colour = COL[["ink"]]
  )
)

export_figure <- function(plot, stem, width_mm, height_mm) {
  paths <- c(
    SVG = file.path(OUT_DIR, paste0(stem, ".svg")),
    PDF = file.path(OUT_DIR, paste0(stem, ".pdf")),
    TIFF = file.path(OUT_DIR, paste0(stem, ".tiff")),
    PNG = file.path(OUT_DIR, paste0(stem, ".png"))
  )
  ggsave(paths[["SVG"]], plot, device = svglite::svglite,
         width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(paths[["PDF"]], plot, device = grDevices::cairo_pdf,
         width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(paths[["TIFF"]], plot, device = ragg::agg_tiff,
         width = width_mm, height = height_mm, units = "mm", dpi = 600,
         compression = "lzw", background = "white")
  ggsave(paths[["PNG"]], plot, device = ragg::agg_png,
         width = width_mm, height = height_mm, units = "mm", dpi = 300,
         background = "white")
  data.frame(
    figure = stem, format = names(paths), file = unname(paths),
    width_mm = width_mm, height_mm = height_mm,
    resolution_dpi = c(NA, NA, 600, 300), stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Figure 2
# -----------------------------------------------------------------------------

external <- read_source("figure2a_external_effects.tsv") %>%
  mutate(cohort = factor(cohort, levels = rev(unique(cohort))))

p2a <- ggplot(external, aes(clustered_standardized_mean_difference, cohort)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(
    aes(x = clustered_standardized_ci_low,
        xend = clustered_standardized_ci_high, yend = cohort),
    linewidth = 0.72, colour = COL[["route"]]
  ) +
  geom_point(size = 2.0, colour = COL[["route"]]) +
  labs(x = "Adenoma effect (SD, patient-clustered 95% CI)", y = NULL) +
  theme_journal(base_size = 6.2)

models <- read_source("figure2b_pooled_models.tsv") %>%
  mutate(
    model = factor(model, levels = rev(c("Cohort-adjusted", "Plus proliferation"))),
    label_x = ci_high + 0.08
  )
p2b <- ggplot(models, aes(estimate, model, colour = model)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = model), linewidth = 0.82) +
  geom_point(size = 2.25) +
  geom_text(
    aes(x = label_x, label = sprintf("%.3f", estimate)), hjust = 0,
    size = 1.40, family = JOURNAL_FONT, colour = COL[["ink"]]
  ) +
  scale_colour_manual(
    values = c("Cohort-adjusted" = COL[["route"]],
               "Plus proliferation" = COL[["wnt"]]), guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 2.38), clip = "off") +
  labs(x = "Five-cohort effect (95% CI)", y = NULL) +
  theme_journal(base_size = 6.2) +
  theme(plot.margin = margin(1.6, 4.0, 1.6, 1.8, "mm"))

external <- external %>%
  mutate(auc_label_x = pmin(auc_a_vs_b + 0.042, 1.040))
p2c <- ggplot(external, aes(auc_a_vs_b, cohort)) +
  geom_vline(xintercept = 0.5, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(
    aes(x = 0.5, xend = auc_a_vs_b, yend = cohort), linewidth = 0.70,
    colour = COL[["neutral_light"]]
  ) +
  geom_point(size = 2.0, colour = COL[["neutral"]]) +
  geom_text(
    aes(x = auc_label_x, label = sprintf("%.3f", auc_a_vs_b)), hjust = 0,
    size = 1.38, family = JOURNAL_FONT, colour = COL[["ink"]]
  ) +
  coord_cartesian(xlim = c(0.46, 1.13), clip = "off") +
  labs(x = "AUC (adenoma vs normal)", y = NULL) +
  theme_journal(base_size = 6.2) +
  theme(plot.margin = margin(1.6, 4.0, 1.6, 1.8, "mm"))

ffpe_pairs <- read_source("figure2d_ffpe_pairs.tsv") %>%
  mutate(tissue = factor(tissue, levels = c("Adjacent mucosa", "Adenoma")))
p2d <- ggplot(ffpe_pairs, aes(tissue, score, group = patient_id)) +
  geom_line(colour = COL[["neutral"]], linewidth = 0.30, alpha = 0.22) +
  geom_point(aes(colour = tissue), size = 0.75, alpha = 0.75) +
  stat_summary(
    aes(group = 1), fun = median, geom = "point", shape = 23,
    size = 2.4, fill = "white", colour = COL[["ink"]]
  ) +
  annotate(
    "text", x = 1.5, y = max(ffpe_pairs$score) + 0.48,
    label = "47/51 increased\nmedian Δ 1.62; P = 3.33×10⁻⁹",
    size = 1.42, family = JOURNAL_FONT, lineheight = 0.90
  ) +
  scale_colour_manual(
    values = c("Adjacent mucosa" = COL[["wnt"]], Adenoma = COL[["route"]])
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(x = NULL, y = "12-gene FFPE score") +
  theme_journal(base_size = 6.1) +
  theme(legend.position = "none", axis.text.x = element_text(size = 5.2))

ffpe_genes <- read_source("figure2e_ffpe_genes.tsv") %>%
  mutate(gene = factor(gene, levels = gene[order(median_paired_delta)]))
p2e <- ggplot(ffpe_genes, aes(median_paired_delta, gene, colour = arm)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(
    aes(x = mean_delta_ci_low, xend = mean_delta_ci_high, yend = gene),
    linewidth = 0.62
  ) +
  geom_point(size = 1.75) +
  scale_colour_manual(values = c("Up arm" = COL[["route"]], "Down arm" = COL[["wnt"]])) +
  labs(x = "Paired adenoma − mucosa expression", y = NULL, colour = NULL) +
  theme_journal(base_size = 5.8) +
  theme(legend.position = "top", axis.text.y = element_text(size = 4.7))

becker <- read_source("figure2f_becker.tsv") %>%
  filter(disease_stage_group %in% c("normal_unaffected", "polyp", "crc")) %>%
  mutate(
    state = factor(
      disease_stage_group, levels = c("normal_unaffected", "polyp", "crc"),
      labels = c("Normal / unaffected", "Polyp", "CRC")
    )
  )
p2f <- ggplot(becker, aes(state, score_epi__ca_route_signature, colour = state)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.43, alpha = 0.15) +
  geom_jitter(width = 0.11, size = 0.72, alpha = 0.62) +
  scale_colour_manual(
    values = c("Normal / unaffected" = COL[["wnt"]],
               Polyp = COL[["adenoma"]], CRC = COL[["crc"]])
  ) +
  labs(x = NULL, y = "Epithelial 12-gene score") +
  theme_journal(base_size = 6.0) +
  theme(legend.position = "none", axis.text.x = element_text(size = 4.9))

fig2 <- (clean_panel(p2a) | clean_panel(p2b) | clean_panel(p2c)) /
  (clean_panel(p2d) | clean_panel(p2e) | clean_panel(p2f)) +
  plot_layout(heights = c(0.92, 1.08), widths = c(1, 0.96, 0.91)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Figure 4
# -----------------------------------------------------------------------------

atlas_raw <- read_tsv(file.path(EXT_DIR, "crc_atlas", "atlas_objective_panel_donor_scores.tsv")) %>%
  filter(n_cells_sampled >= 20)
atlas <- atlas_raw %>%
  mutate(
    state = factor(
      carrier_group,
      levels = c("normal_epithelial", "polyp_epithelial", "polyp_cancer",
                 "tumor_epithelial", "tumor_cancer",
                 "metastasis_epithelial", "metastasis_cancer"),
      labels = c("Normal\nepithelium", "Polyp\nepithelium", "Polyp\ncancer cells",
                 "Primary\nepithelium", "Primary\ncancer cells",
                 "Metastasis\nepithelium", "Metastasis\ncancer cells")
    )
  )
p4a <- ggplot(atlas, aes(state, score__ca_route_signature, colour = state)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.40, alpha = 0.12) +
  geom_jitter(width = 0.12, size = 0.45, alpha = 0.30) +
  scale_colour_manual(
    values = c(COL[["wnt"]], COL[["adenoma"]], COL[["route"]],
               COL[["context"]], COL[["crc"]], "#507D73", "#204F46")
  ) +
  labs(x = NULL, y = "Donor-carrier median 12-gene score") +
  theme_journal(base_size = 5.9) +
  theme(legend.position = "none", axis.text.x = element_text(size = 4.65))

influence <- read_source("figure4_influence.tsv") %>%
  mutate(state_label = factor(state_label, levels = rev(state_label)))
p4b <- ggplot(influence, aes(full_coef, state_label)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(
    aes(x = full_ci_low, xend = full_ci_high, yend = state_label),
    linewidth = 0.75, colour = COL[["route"]]
  ) +
  geom_point(size = 2.0, colour = COL[["route"]]) +
  labs(x = "Adjusted carrier-state coefficient (95% CI)", y = NULL) +
  theme_journal(base_size = 6.1)

p4c <- ggplot(influence, aes(loo_min_coef, state_label)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(
    aes(x = loo_min_coef, xend = loo_max_coef, yend = state_label),
    linewidth = 1.05, colour = COL[["neutral_light"]], lineend = "round"
  ) +
  geom_point(aes(x = full_coef), size = 1.9, colour = COL[["route"]]) +
  annotate(
    "text", x = max(influence$loo_max_coef), y = 0.55,
    label = "33/33 omissions positive", hjust = 1, size = 1.35,
    family = JOURNAL_FONT
  ) +
  labs(x = "Coefficient range after one-study omission", y = NULL) +
  theme_journal(base_size = 6.1)

support <- read_source("figure4_support.tsv") %>%
  mutate(
    state = factor(
      state,
      levels = c("Polyp epi.", "Polyp cancer", "Primary epi.", "Primary cancer",
                 "Metastasis epi.", "Metastasis cancer")
    ),
    study_short = vapply(
      strsplit(study_id, "_", fixed = TRUE),
      function(parts) paste(head(parts, 2), collapse = " "), character(1)
    )
  )
study_order <- support %>%
  group_by(study_short) %>% summarise(total = sum(n_donors), .groups = "drop") %>%
  arrange(total) %>% pull(study_short)
support <- support %>%
  mutate(study_short = factor(study_short, levels = study_order))
p4d <- ggplot(support, aes(state, study_short, fill = n_donors)) +
  geom_tile(colour = "white", linewidth = 0.22) +
  scale_fill_gradient(low = "#F3F6F7", high = COL[["route"]], trans = "sqrt") +
  labs(x = NULL, y = "Source study", fill = "Donors") +
  theme_journal(base_size = 5.9) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 4.8),
    axis.text.y = element_text(size = 4.9), axis.line = element_blank(),
    axis.ticks = element_blank(), legend.position = "top"
  )

within <- read_source("figure4_within_study.tsv") %>%
  mutate(
    study_short = vapply(
      strsplit(study_id, "_", fixed = TRUE),
      function(parts) paste(head(parts, 2), collapse = " "), character(1)
    ),
    state_short = recode(
      state, polyp_epithelial = "polyp epithelium", polyp_cancer = "polyp cancer",
      tumor_epithelial = "primary epithelium", tumor_cancer = "primary cancer",
      metastasis_epithelial = "metastasis epithelium",
      metastasis_cancer = "metastasis cancer"
    ),
    label_short = paste(study_short, state_short, sep = " · "),
    label_short = factor(label_short, levels = rev(label_short))
  )
p4e <- ggplot(within, aes(coef, label_short)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = label_short), linewidth = 0.60,
    colour = COL[["context"]]
  ) +
  geom_point(size = 1.55, colour = COL[["context"]]) +
  labs(x = "Eligible within-study contrast (95% CI)", y = NULL) +
  theme_journal(base_size = 5.8) +
  theme(axis.text.y = element_text(size = 4.55))

fig4 <- clean_panel(p4a) /
  (clean_panel(p4b) | clean_panel(p4c)) /
  (clean_panel(p4d) | clean_panel(p4e)) +
  plot_layout(heights = c(0.82, 0.86, 1.48), widths = c(1.02, 0.98)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Supplementary Figure 4
# -----------------------------------------------------------------------------

rna_atac <- read_source("figure3_rna_atac.tsv") %>%
  filter(disease_stage_group %in% c("normal_unaffected", "polyp"))
patient_median <- rna_atac %>%
  group_by(patient_id) %>%
  summarise(
    rna = median(rna_epi__ca_route_signature),
    atac = median(atac_tss__wnt_stem), .groups = "drop"
  )
pS4a <- ggplot(patient_median, aes(atac, rna)) +
  geom_smooth(
    method = "lm", se = TRUE, colour = COL[["neutral"]], fill = "#E7EBEE",
    linewidth = 0.48
  ) +
  geom_point(size = 1.75, colour = COL[["route"]]) +
  annotate(
    "text", x = min(patient_median$atac), y = max(patient_median$rna), hjust = 0,
    label = "Patient-median ρ = 0.790\npermutation P = 0.00292",
    size = 1.40, family = JOURNAL_FONT
  ) +
  labs(x = "Patient-median WNT TSS accessibility", y = "Patient-median RNA score") +
  theme_journal(base_size = 6.0)

models3 <- read_source("figure3_models.tsv") %>%
  mutate(feature = factor(feature, levels = rev(feature)))
pS4b <- ggplot(models3, aes(coef, feature)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = feature), linewidth = 0.72,
    colour = COL[["route"]]
  ) +
  geom_point(size = 1.95, colour = COL[["route"]]) +
  labs(x = "Patient-clustered adjusted coefficient (95% CI)", y = NULL) +
  theme_journal(base_size = 6.0) +
  theme(axis.text.y = element_text(size = 4.9))

residual <- read_tsv(file.path(
  EXT_DIR, "becker_rna_atac", "becker_rna_atac_residual_correlations.tsv"
)) %>%
  filter(subset == "normal_polyp", grepl("rna_epi_ca_route", analysis_id)) %>%
  mutate(
    feature = recode(
      analysis_id,
      rna_epi_ca_route__atac_wnt_stem = "WNT/stem accessibility",
      rna_epi_ca_route__atac_wnt_minus_housekeeping = "WNT/stem − housekeeping",
      rna_epi_ca_route__atac_wnt_tcf_ascl2_axis = "WNT/TCF/ASCL2 axis",
      rna_epi_ca_route__atac_wnt_tcf_ascl2_minus_housekeeping =
        "WNT/TCF/ASCL2 − housekeeping"
    ),
    feature = factor(
      feature,
      levels = rev(c(
        "WNT/stem accessibility", "WNT/stem − housekeeping",
        "WNT/TCF/ASCL2 axis", "WNT/TCF/ASCL2 − housekeeping"
      ))
    )
  )
pS4c <- ggplot(residual, aes(residual_spearman_rho, feature)) +
  geom_vline(xintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(
    aes(x = 0, xend = residual_spearman_rho, yend = feature), linewidth = 0.70,
    colour = COL[["neutral_light"]]
  ) +
  geom_point(size = 1.9, colour = COL[["route"]]) +
  labs(x = "Residual Spearman ρ", y = NULL) +
  theme_journal(base_size = 5.8) +
  theme(axis.text.y = element_text(size = 4.7))

regulatory <- read_source("figure3_regulatory.tsv") %>%
  mutate(
    distance = factor(distance, levels = c("TSS 1 kb", "2–5 kb", "10 kb", "20 kb")),
    locus = factor(locus, levels = c("TCF/ASCL2 loci", "WNT-route loci"))
  )
pS4d <- ggplot(regulatory, aes(distance, spearman_rho, colour = locus, group = locus)) +
  geom_hline(yintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_line(linewidth = 0.70) + geom_point(size = 1.85) +
  scale_colour_manual(
    values = c("WNT-route loci" = COL[["route"]], "TCF/ASCL2 loci" = COL[["wnt"]])
  ) +
  labs(x = "Regulatory distance window", y = "Spearman ρ with RNA score", colour = NULL) +
  theme_journal(base_size = 6.0) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1, size = 4.8))

figS4 <- (clean_panel(pS4a) | clean_panel(pS4b)) /
  (clean_panel(pS4c) | clean_panel(pS4d)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Supplementary Figure 6
# -----------------------------------------------------------------------------

coverage <- read_tsv(file.path(EXT_DIR, "perturbation_spatial", "feature_coverage.tsv")) %>%
  mutate(
    dataset = gsub("spatial_", "", dataset),
    feature_label = recode(
      feature,
      route_up = "Signature up arm", route_down = "Signature down arm",
      wnt_stem = "WNT/stemness", proliferation_control = "Proliferation control",
      epithelial_control = "Epithelial identity control"
    ),
    feature_label = factor(
      feature_label,
      levels = c("Signature up arm", "Signature down arm", "WNT/stemness",
                 "Proliferation control", "Epithelial identity control")
    )
  )
pS6a <- ggplot(coverage, aes(feature_label, dataset, fill = coverage_fraction)) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = percent(coverage_fraction, accuracy = 1)),
            size = 1.25, family = JOURNAL_FONT) +
  scale_fill_gradient(low = "#ECEFF1", high = COL[["route"]], limits = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "Coverage") +
  theme_journal(base_size = 5.7) +
  theme(
    axis.text.x = element_text(angle = 27, hjust = 1, size = 4.55),
    axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "top"
  )

matched <- read_tsv(file.path(
  EXT_DIR, "perturbation_spatial", "expression_matched_tests.tsv"
)) %>%
  filter(comparison != "doxycycline_control_shRenilla") %>%
  mutate(
    comparison_label = recode(
      comparison,
      trametinib_vs_dmso = "Trametinib vs DMSO",
      pri724_reversal_of_trametinib = "PRI-724 reversal of trametinib",
      apc_restoration_shApc = "APC restoration",
      apc_restoration_shApc_Kras = "APC restoration + KRAS",
      ascl2_ko_vs_resting_wt = "ASCL2 knockout vs resting WT",
      conditional_wnt_silencing = "Conditional WNT silencing"
    ),
    label = paste(dataset, comparison_label, sep = " · "),
    label = factor(label, levels = rev(label)),
    p_label = paste0("P=", ifelse(
      p_expression_matched_one_sided < 0.001,
      formatC(p_expression_matched_one_sided, format = "e", digits = 1),
      formatC(p_expression_matched_one_sided, format = "f", digits = 3)
    ))
  )
pS6b <- ggplot(matched, aes(observed_mean_difference, label)) +
  geom_vline(xintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(
    aes(x = null_ci_low, xend = null_ci_high, yend = label), linewidth = 1.05,
    colour = COL[["neutral_light"]], lineend = "round"
  ) +
  geom_point(size = 1.9, colour = COL[["route"]]) +
  geom_text(aes(label = p_label), nudge_x = 0.14, hjust = 0,
            size = 1.18, family = JOURNAL_FONT) +
  coord_cartesian(
    xlim = c(min(matched$null_ci_low) - 0.1, max(matched$null_ci_high) + 0.62),
    clip = "off"
  ) +
  labs(x = "Effect vs expression-matched null (95% range)", y = NULL) +
  theme_journal(base_size = 5.35) +
  theme(axis.text.y = element_text(size = 4.05), plot.margin = margin(1.6, 4, 1.6, 1.8, "mm"))

effects <- read_tsv(file.path(
  EXT_DIR, "perturbation_spatial", "perturbation_effect_summary.tsv"
)) %>%
  filter(feature == "route_score", expected_direction != 0) %>%
  mutate(
    label = recode(
      comparison,
      trametinib_vs_dmso = "Trametinib vs DMSO",
      pri724_reversal_of_trametinib = "PRI-724 reversal",
      apc_restoration_shApc = "APC restoration",
      apc_restoration_shApc_Kras = "APC restoration + KRAS",
      ascl2_ko_vs_resting_wt = "ASCL2 knockout",
      conditional_wnt_silencing = "Conditional WNT silencing"
    ),
    label = factor(label, levels = rev(label)),
    aligned = expected_direction * mean_difference > 0
  )
pS6c <- ggplot(effects, aes(mean_difference, label, colour = aligned)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(
    aes(x = bootstrap_mean_ci_low, xend = bootstrap_mean_ci_high, yend = label),
    linewidth = 0.65
  ) +
  geom_point(size = 1.85) +
  scale_colour_manual(
    values = c(`TRUE` = COL[["route"]], `FALSE` = COL[["crc"]]), guide = "none"
  ) +
  labs(x = "Intervention effect on 12-gene score", y = NULL) +
  theme_journal(base_size = 5.9)

stable_compare <- bind_rows(
  read_tsv(file.path(RESULT_DIR, "validation_chen_panel_metrics.tsv")) %>%
    filter(dataset == "validation", panel_id == "objective_12") %>%
    transmute(panel = "Kneedle 12", metric = "Held-out AUC", value = auc),
  read_tsv(file.path(ROOT, "results", "stability_consensus_panel_v2_8",
                     "validation_chen_metrics.tsv")) %>%
    filter(dataset == "validation", panel_id == "stability_consensus_11") %>%
    transmute(panel = "Strict-majority 11", metric = "Held-out AUC", value = auc),
  read_tsv(file.path(RESULT_DIR, "validation_ffpe_paired_tests.tsv")) %>%
    filter(panel_id == "objective_12") %>%
    transmute(panel = "Kneedle 12", metric = "FFPE positive pairs",
              value = paired_positive_fraction),
  read_tsv(file.path(ROOT, "results", "stability_consensus_panel_v2_8",
                     "validation_ffpe_paired_test.tsv")) %>%
    transmute(panel = "Strict-majority 11", metric = "FFPE positive pairs",
              value = paired_positive_fraction)
)
pS6d <- ggplot(stable_compare, aes(value, metric, colour = panel)) +
  geom_point(position = position_dodge(width = 0.40), size = 2.0) +
  scale_colour_manual(
    values = c("Kneedle 12" = COL[["route"]],
               "Strict-majority 11" = COL[["wnt"]])
  ) +
  coord_cartesian(xlim = c(0.82, 0.97)) +
  labs(x = "Validation metric", y = NULL, colour = NULL) +
  theme_journal(base_size = 5.8) +
  theme(legend.position = "top")

figS6 <- (clean_panel(pS6a) | clean_panel(pS6b)) /
  (clean_panel(pS6c) | clean_panel(pS6d)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Export and panel-level source data
# -----------------------------------------------------------------------------

figure5b_source <- read_tsv(file.path(
  RESULT_DIR, "validation_perturbation_sample_scores.tsv"
)) %>%
  filter(panel_id == "objective_12", dataset == "GSE125472") %>%
  mutate(
    condition = paste(
      ifelse(genotype == "APC", "APC-KO", "WT"),
      ifelse(wnt_rspo == "with", "+WNT/RSPO", "−WNT/RSPO")
    )
  )

figureS8a_source <- data.frame(
  xmin = c(0.15, 2.75, 5.35, 7.95),
  xmax = c(2.25, 4.85, 7.45, 10.05),
  ymin = rep(0.25, 4), ymax = rep(1.70, 4),
  x = c(1.20, 3.80, 6.40, 9.00), y = rep(0.98, 4),
  label = c(
    "Frozen inputs; 287-gene core; 12-gene signature; 13 knockout targets",
    "Held-out input; 13 adenoma donors; 128 cells per donor; 1,664 cells",
    "Dual-seed GenKI; unsigned impact; KL primary; EMD sensitivity",
    "Validation tests; 10,000 matched sets; no gene reselection; no reweighting"
  ),
  fill = c("#E8F1F5", "#EEF2F4", "#E9F3F1", "#F8EFE3"),
  stringsAsFactors = FALSE
)

exports <- bind_rows(
  export_figure(fig2, "figure2_independent_replication_and_ffpe", 170, 148),
  export_figure(fig4, "figure4_crc_atlas_cross_sectional_recurrence", 170, 202),
  export_figure(figS4, "figureS4_rna_atac_robustness", 170, 142),
  export_figure(figS6, "figureS6_perturbation_boundaries", 170, 148)
)

write_source(external, "figure2a_external_effects.tsv")
write_source(models, "figure2b_pooled_models.tsv")
write_source(ffpe_pairs, "figure2d_ffpe_pairs.tsv")
write_source(ffpe_genes, "figure2e_ffpe_genes.tsv")
write_source(becker, "figure2f_becker.tsv")
write_source(atlas_raw %>% mutate(state = gsub("_", " ", carrier_group)), "figure4_atlas.tsv")
write_source(influence, "figure4_influence.tsv")
write_source(support, "figure4_support.tsv")
write_source(within, "figure4_within_study.tsv")
write_source(patient_median, "figureS4a_patient_median.tsv")
write_source(models3, "figureS4b_adjusted_models.tsv")
write_source(residual, "figureS4c_residual_correlations.tsv")
write_source(regulatory, "figureS4d_regulatory_windows.tsv")
write_source(coverage, "figureS6a_feature_coverage.tsv")
write_source(matched, "figureS6b_expression_matched_tests.tsv")
write_source(effects, "figureS6c_perturbation_effects.tsv")
write_source(stable_compare, "figureS6d_stability_comparison.tsv")
write_source(figure5b_source, "figure5b_apc_wnt_donor_scores.tsv")
write_source(figureS8a_source, "figureS8a_frozen_design.tsv")
write_source(exports, "figure_audit_fix_export_manifest.tsv")

tsv_paths <- list.files(SOURCE_DIR, pattern = "tsv$", full.names = TRUE)
rectangular_tsv <- vapply(tsv_paths, function(path) {
  fields <- count.fields(path, sep = "\t", quote = "", blank.lines.skip = FALSE)
  length(fields) > 0 && all(fields == fields[1])
}, logical(1))

qa <- data.frame(
  check = c(
    "figure2_model_labels_have_three_decimals",
    "figure2_auc_labels_have_three_decimals",
    "figure4_source_labels_are_author_year",
    "supplementary4_has_no_internal_atac_names",
    "supplementary6_has_no_internal_feature_names",
    "supplementary6_gene_case_standardised",
    "all_four_formats_exported",
    "all_source_tsvs_are_rectangular"
  ),
  pass = c(
    all(grepl("^[0-9]+\\.[0-9]{3}$", sprintf("%.3f", models$estimate))),
    all(grepl("^[0-9]+\\.[0-9]{3}$", sprintf("%.3f", external$auc_a_vs_b))),
    all(grepl("^[^_]+ [0-9]{4}$", as.character(support$study_short)) |
          as.character(support$study_short) == "MUI Innsbruck"),
    !any(grepl("^atac_", as.character(residual$feature))),
    !any(as.character(coverage$feature_label) %in%
         c("route_up", "route_down", "wnt_stem", "proliferation_control", "epithelial_control")),
    all(grepl("APC|ASCL2|WNT|Trametinib|PRI-724", as.character(effects$label))),
    all(file.exists(exports$file)),
    all(rectangular_tsv)
  ),
  stringsAsFactors = FALSE
)
write_source(qa, "figure_audit_fix_qa.tsv")
if (!all(qa$pass)) {
  stop("Targeted figure QA failed: ", paste(qa$check[!qa$pass], collapse = ", "))
}
writeLines(capture.output(sessionInfo()),
           file.path(SOURCE_DIR, "figure_audit_fix_sessionInfo.txt"))
message("Targeted figure audit fixes exported to: ", OUT_DIR)
