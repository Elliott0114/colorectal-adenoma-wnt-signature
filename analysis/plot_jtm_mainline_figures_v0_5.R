#!/usr/bin/env Rscript

# Recompose the locked computational evidence into a single-mainline figure set
# for Journal of Translational Medicine. This script reads immutable source
# tables and does not refit or alter any statistical model.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(tidyr)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[1])
} else {
  file.path(getwd(), "analysis", "plot_jtm_mainline_figures_v0_5.R")
}

ROOT <- normalizePath(file.path(dirname(script_path), ".."))
DATA_DIR <- file.path(ROOT, "results", "figure_data_locked")
OUT_DIR <- Sys.getenv(
  "JTM_FIGURE_DIR",
  unset = file.path(ROOT, "figures", "jtm_mainline_v0.5")
)
DERIVED_DIR <- file.path(OUT_DIR, "source_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)

read_source <- function(filename) {
  read.delim(
    file.path(DATA_DIR, filename),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

write_source <- function(frame, filename) {
  write.table(
    frame,
    file.path(DERIVED_DIR, filename),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

superscript_integer <- function(value) {
  chars <- strsplit(as.character(value), "", fixed = TRUE)[[1]]
  paste0(chartr("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹", chars), collapse = "")
}

format_p <- function(p) {
  vapply(p, function(value) {
    if (!is.finite(value)) return("not estimable")
    if (value < 0.001) {
      exponent <- floor(log10(value))
      mantissa <- value / (10 ^ exponent)
      return(sprintf("%.2f × 10%s", mantissa, superscript_integer(exponent)))
    }
    sprintf("%.3f", value)
  }, character(1))
}

JTM_FONT <- "DejaVu Sans"

COL <- c(
  ink = "#202428",
  neutral = "#6E7781",
  neutral_light = "#D9DEE3",
  neutral_pale = "#F2F4F5",
  route = "#C45A3E",
  route_light = "#E9B6A8",
  wnt = "#2F6F9F",
  wnt_light = "#B8D0E1",
  adenoma = "#D9A441",
  crc = "#2F8F83",
  metastasis = "#8064A2",
  uncertain = "#8C7A6B"
)

theme_jtm <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = JTM_FONT) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = COL[["ink"]]),
      axis.ticks = element_line(linewidth = 0.30, colour = COL[["ink"]]),
      axis.ticks.length = unit(1.1, "mm"),
      axis.title = element_text(size = base_size, colour = COL[["ink"]]),
      axis.text = element_text(size = base_size - 0.6, colour = COL[["ink"]]),
      panel.grid = element_blank(),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.7),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width = unit(4.0, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_text(size = base_size - 1.0, colour = COL[["neutral"]], hjust = 0),
      plot.margin = margin(2.4, 2.4, 2.4, 2.4, unit = "mm"),
      plot.tag = element_text(size = 9, face = "bold", colour = COL[["ink"]])
    )
}

tag_theme <- theme(
  plot.tag.position = c(0, 1),
  plot.tag = element_text(size = 9, face = "bold", colour = COL[["ink"]])
)

export_figure <- function(plot, stem, width_mm = 170, height_mm = 120) {
  paths <- c(
    SVG = file.path(OUT_DIR, paste0(stem, ".svg")),
    PDF = file.path(OUT_DIR, paste0(stem, ".pdf")),
    TIFF = file.path(OUT_DIR, paste0(stem, ".tiff")),
    PNG = file.path(OUT_DIR, paste0(stem, ".png"))
  )

  ggsave(paths[["SVG"]], plot = plot, device = svglite::svglite,
         width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(paths[["PDF"]], plot = plot, device = grDevices::cairo_pdf,
         width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(paths[["TIFF"]], plot = plot, device = ragg::agg_tiff,
         width = width_mm, height = height_mm, units = "mm", dpi = 600,
         compression = "lzw", background = "white")
  ggsave(paths[["PNG"]], plot = plot, device = ragg::agg_png,
         width = width_mm, height = height_mm, units = "mm", dpi = 300,
         background = "white")

  data.frame(
    figure = stem,
    format = names(paths),
    file = unname(paths),
    width_mm = width_mm,
    height_mm = height_mm,
    resolution_dpi = c(NA, NA, 600, 300),
    stringsAsFactors = FALSE
  )
}

forest_panel <- function(data, estimate, low, high, y, colour, title, subtitle, xlab,
                         xlim = NULL, point_size = 1.9) {
  p <- ggplot(data, aes(x = {{ estimate }}, y = {{ y }})) +
    geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
    geom_segment(aes(x = {{ low }}, xend = {{ high }}, yend = {{ y }}),
                 linewidth = 0.65, colour = COL[["ink"]]) +
    geom_point(size = point_size, colour = colour) +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    theme_jtm()
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim, clip = "off")
  p
}

# -----------------------------------------------------------------------------
# Shared locked inputs and numerical gates
# -----------------------------------------------------------------------------

audit <- read_source("fig1_discovery_stability_audit.tsv")
signature <- read_source("fig1_locked_signature_genes.tsv")
chen <- read_source("fig1_chen_locked_scores.tsv")
chen_discrimination <- read_source("fig1_chen_locked_discrimination.tsv")
chen_paired_tests <- read_source("fig1_chen_locked_paired_tests.tsv")
becker <- read_source("fig1_becker_locked_scores.tsv")
becker_models <- read_source("fig1_becker_patient_cluster_models.tsv")

external_tests <- read_source("figs3_external_cohort_tests.tsv")
external_models <- read_source("figs3_external_one_stage_models.tsv")
external_loo <- read_source("figs3_external_leave_one_cohort_out.tsv")
external_prolif <- read_source("figs3_external_proliferation_adjusted_model.tsv")
external_prolif_loo <- read_source("figs3_external_proliferation_adjusted_leave_one_cohort_out.tsv")
external_scores <- read_source("figs3_external_sample_scores.tsv")
external_grade <- read_source("figs3_external_grade_sensitivity.tsv")
external_coverage <- read_source("figs3_external_signature_coverage.tsv")

paired <- read_source("fig2_locked_rna_atac_paired_scores.tsv") %>%
  filter(disease_stage_group %in% c("normal_unaffected", "polyp"))
rna_atac_cor <- read_source("fig2_locked_rna_atac_correlations.tsv")
rna_atac_models <- read_source("fig2_locked_rna_atac_patient_cluster_models.tsv")
rna_atac_fixed <- read_source("fig2_locked_rna_atac_patient_fixed_effect_models.tsv")
rna_atac_medians <- read_source("fig2_locked_rna_atac_patient_median_correlations.tsv")
rna_atac_boot <- read_source("fig2_locked_rna_atac_patient_cluster_bootstrap.tsv")
window_cor <- read_source("fig2_locked_regulatory_window_correlations.tsv")

atlas_scores <- read_source("fig3_atlas_locked_donor_scores.tsv")
atlas_models <- read_source("fig3_atlas_donor_cluster_models.tsv")
atlas_influence <- read_source("figs2_atlas_study_influence.tsv")
atlas_support <- read_source("figs2_atlas_state_study_support.tsv")
atlas_within <- read_source("figs2_atlas_within_study_contrasts.tsv")
platform_coverage <- read_source("supp_locked_signature_platform_coverage.tsv")

stopifnot(nrow(signature) == 100)
stopifnot(all(tolower(as.character(signature$validation_used_for_selection)) == "false"))
stopifnot(abs(chen_discrimination$auc_adenoma_vs_normal[chen_discrimination$dataset == "validation"] - 0.9285714) < 1e-6)
stopifnot(nrow(paired) == 40, length(unique(paired$patient_id)) == 12)
stopifnot(all(atlas_influence$n_eligible_omissions == 33))
stopifnot(all(atlas_influence$loo_positive_fraction == 1))

# -----------------------------------------------------------------------------
# Main Figure 1: lock, held-out validation, and independent sporadic cohorts
# -----------------------------------------------------------------------------

workflow <- data.frame(
  x = 1:4,
  label = c(
    "DISCOVER\nChen donor medians\n1,000 bootstraps",
    "LOCK\n50 up + 50 down\nvalidation unseen",
    "VALIDATE\nHeld-out Chen\npaired donors",
    "TRANSPORT\n5 sporadic cohorts\n203 samples / 161 patients"
  ),
  stage = factor(c("Discovery", "Lock", "Validation", "Transport"),
                 levels = c("Discovery", "Lock", "Validation", "Transport"))
)

p1a <- ggplot(workflow, aes(x, 1)) +
  geom_segment(
    data = data.frame(x = 1:3, xend = 2:4, y = 1, yend = 1),
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, linewidth = 0.55, colour = COL[["neutral"]],
    arrow = arrow(length = unit(1.6, "mm"), type = "closed")
  ) +
  geom_label(
    aes(label = label, fill = stage), size = 2.05, linewidth = 0.28,
    label.padding = unit(1.5, "mm"), lineheight = 0.92, colour = COL[["ink"]]
  ) +
  scale_fill_manual(values = c(
    Discovery = "#E7EEF3", Lock = "#F0D6CE",
    Validation = "#F5E7C7", Transport = "#DDECE8"
  )) +
  coord_cartesian(xlim = c(0.55, 4.45), ylim = c(0.72, 1.28), clip = "off") +
  labs(title = "Leakage-resistant study design") +
  theme_void(base_size = 7, base_family = JTM_FONT) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.margin = margin(4, 3, 3, 3, unit = "mm")
  )

audit <- audit %>%
  mutate(
    selected_flag = tolower(as.character(selected)) == "true",
    direction = case_when(
      selected_flag & discovery_effect_adenoma_minus_normal > 0 ~ "Adenoma-up",
      selected_flag & discovery_effect_adenoma_minus_normal < 0 ~ "Adenoma-down",
      TRUE ~ "Not locked"
    ),
    direction = factor(direction, levels = c("Not locked", "Adenoma-down", "Adenoma-up"))
  )
audit_labels <- audit %>%
  filter(FALSE)

p1b <- ggplot(audit, aes(discovery_effect_adenoma_minus_normal, direction_stability)) +
  geom_point(data = audit %>% filter(!selected_flag),
             colour = COL[["neutral_light"]], alpha = 0.20, size = 0.42, stroke = 0) +
  geom_point(data = audit %>% filter(selected_flag), aes(colour = direction),
             alpha = 0.90, size = 1.05, stroke = 0) +
  geom_hline(yintercept = 0.90, linetype = "22", linewidth = 0.35, colour = COL[["neutral"]]) +
  ggrepel::geom_text_repel(
    data = audit_labels, aes(label = gene, colour = direction),
    size = 1.9, box.padding = 0.25, point.padding = 0.12,
    min.segment.length = 0, segment.size = 0.25, max.overlaps = Inf,
    seed = 20260710, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Not locked" = COL[["neutral_light"]],
    "Adenoma-down" = COL[["wnt"]],
    "Adenoma-up" = COL[["route"]]
  ), breaks = c("Adenoma-down", "Adenoma-up"), name = NULL) +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(
    title = "Discovery-only stability screen",
    subtitle = "33,698 genes; dashed line, prespecified 0.90 stability threshold",
    x = "Discovery donor-median effect (adenoma - normal)",
    y = "Directional stability"
  ) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

chen_plot <- chen %>%
  filter(route_group %in% c("normal", "conventional_adenoma")) %>%
  mutate(
    dataset = factor(dataset, levels = c("discovery", "validation"),
                     labels = c("Discovery", "Held-out validation")),
    tissue = factor(route_group, levels = c("normal", "conventional_adenoma"),
                    labels = c("Normal", "Conventional adenoma"))
  )
chen_ann <- chen_discrimination %>%
  mutate(
    dataset = factor(dataset, levels = c("discovery", "validation"),
                     labels = c("Discovery", "Held-out validation")),
    label = sprintf("AUC %.3f\nP = %s", auc_adenoma_vs_normal, format_p(p_mannwhitney))
  ) %>%
  left_join(
    chen_plot %>% group_by(dataset) %>%
      summarise(y = max(score__ca_route_signature, na.rm = TRUE) + 0.55, .groups = "drop"),
    by = "dataset"
  )

p1c <- ggplot(chen_plot %>% filter(dataset == "Held-out validation"),
               aes(tissue, score__ca_route_signature, colour = tissue)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.43,
               colour = COL[["ink"]], fill = "white") +
  geom_point(position = position_jitter(width = 0.11, seed = 20260710),
             alpha = 0.72, size = 0.90) +
  geom_text(data = chen_ann %>% filter(dataset == "Held-out validation"),
            aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 1.95, lineheight = 0.95) +
  scale_colour_manual(values = c("Normal" = COL[["neutral"]],
                                 "Conventional adenoma" = COL[["route"]])) +
  labs(title = "Held-out Chen validation", x = NULL, y = "Locked route score") +
  theme_jtm() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 18, hjust = 1, size = 5.5))

external_primary <- external_tests %>%
  filter(comparison == "adenoma_vs_normal", signature_size_per_direction == 50) %>%
  mutate(
    cohort_label = factor(
      sprintf("%s  (%d/%d)", cohort, n_a, n_b),
      levels = sprintf(
        "%s  (%d/%d)",
        cohort[match(c("GSE72820", "GSE40362", "GSE41657", "GSE50114", "GSE8671"), cohort)],
        n_a[match(c("GSE72820", "GSE40362", "GSE41657", "GSE50114", "GSE8671"), cohort)],
        n_b[match(c("GSE72820", "GSE40362", "GSE41657", "GSE50114", "GSE8671"), cohort)]
      )
    ),
    evidence = ifelse(clustered_standardized_ci_low > 0, "Clear", "Imprecise")
  )
stopifnot(nrow(external_primary) == 5)

p1d <- ggplot(external_primary,
              aes(clustered_standardized_mean_difference, cohort_label, colour = evidence)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(aes(x = clustered_standardized_ci_low,
                   xend = clustered_standardized_ci_high, yend = cohort_label), linewidth = 0.68) +
  geom_point(size = 2.05) +
  scale_colour_manual(values = c(Clear = COL[["route"]], Imprecise = COL[["uncertain"]]), name = NULL) +
  coord_cartesian(xlim = c(-0.80, 2.65)) +
  labs(
    title = "External validation",
    subtitle = "Clustered 95% CIs",
    x = "Standardised adenoma effect", y = NULL
  ) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

chen_pairs <- chen %>%
  filter(route_group %in% c("normal", "conventional_adenoma")) %>%
  group_by(dataset, donor_id, route_group) %>%
  summarise(score = median(score__ca_route_signature, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = route_group, values_from = score) %>%
  filter(is.finite(normal), is.finite(conventional_adenoma)) %>%
  mutate(
    delta = conventional_adenoma - normal,
    dataset = factor(dataset, levels = c("discovery", "validation"),
                     labels = c("Discovery", "Held-out"))
  )
pair_ann <- chen_paired_tests %>%
  filter(score == "ca_route_signature", dataset %in% c("discovery", "validation")) %>%
  mutate(
    dataset = factor(dataset, levels = c("discovery", "validation"),
                     labels = c("Discovery", "Held-out")),
    label = sprintf("n = %d; P = %s", n_pairs, format_p(p_wilcoxon))
  ) %>%
  left_join(
    chen_pairs %>% group_by(dataset) %>% summarise(y = max(delta) + 0.40, .groups = "drop"),
    by = "dataset"
  ) %>%
  mutate(hjust = ifelse(dataset == "Discovery", 0, 1))

p1e <- ggplot(chen_pairs, aes(dataset, delta, colour = dataset)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_point(position = position_jitter(width = 0.10, seed = 20260710), size = 1.2, alpha = 0.82) +
  stat_summary(fun = median, geom = "crossbar", width = 0.48, linewidth = 0.58, colour = COL[["ink"]]) +
  geom_text(data = pair_ann, aes(x = dataset, y = y, label = label, hjust = hjust),
            inherit.aes = FALSE, size = 1.95) +
  scale_colour_manual(values = c(Discovery = COL[["adenoma"]], `Held-out` = COL[["route"]])) +
  labs(title = "Within-donor confirmation", x = NULL,
       y = "Paired score change\n(adenoma - normal)") +
  theme_jtm() +
  theme(legend.position = "none")

primary_model <- external_models %>%
  filter(signature_size_per_direction == 50, excluded_cohort == "__NONE__")
primary_prolif <- external_prolif %>%
  filter(signature_size_per_direction == 50, excluded_cohort == "__NONE__")
robust <- bind_rows(
  data.frame(
    model = "Primary cohort-fixed",
    estimate = primary_model$adenoma_coef_sd,
    ci_low = primary_model$ci_low,
    ci_high = primary_model$ci_high,
    loo_low = min(external_loo$adenoma_coef_sd[external_loo$signature_size_per_direction == 50]),
    loo_high = max(external_loo$adenoma_coef_sd[external_loo$signature_size_per_direction == 50])
  ),
  data.frame(
    model = "+ proliferation control",
    estimate = primary_prolif$adenoma_coef_sd,
    ci_low = primary_prolif$ci_low,
    ci_high = primary_prolif$ci_high,
    loo_low = min(external_prolif_loo$adenoma_coef_sd),
    loo_high = max(external_prolif_loo$adenoma_coef_sd)
  )
) %>%
  mutate(model = factor(model, levels = rev(c("Primary cohort-fixed", "+ proliferation control"))))

p1f <- ggplot(robust, aes(estimate, model)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(aes(x = loo_low, xend = loo_high, yend = model),
               linewidth = 2.0, colour = COL[["neutral_light"]], lineend = "round") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = model),
               linewidth = 0.70, colour = COL[["ink"]]) +
  geom_point(size = 2.05, colour = COL[["route"]]) +
  geom_text(aes(x = ci_high + 0.08, label = sprintf("%.2f SD", estimate)),
            hjust = 0, size = 1.85) +
  coord_cartesian(xlim = c(0, 2.35), clip = "off") +
  labs(
    title = "Pooled effect and cohort-exclusion stability",
    subtitle = "Black, 95% CI; pale bar, range after each cohort exclusion",
    x = "Standardised adenoma effect", y = NULL
  ) +
  theme_jtm()

fig1_design <- "
AAAAAA
BBCCDD
EEFFFF
"
fig1 <- p1a + p1b + p1c + p1d + p1e + p1f +
  plot_layout(design = fig1_design, heights = c(0.62, 1.10, 0.90)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_source(external_primary, "fig1_external_cohort_effects.tsv")
write_source(robust, "fig1_pooled_robustness.tsv")
write_source(chen_pairs, "fig1_chen_paired_changes.tsv")

# -----------------------------------------------------------------------------
# Main Figure 2: Becker transfer and matched RNA-ATAC coupling
# -----------------------------------------------------------------------------

becker_plot <- becker %>%
  mutate(group = factor(disease_stage_group,
                        levels = c("normal_unaffected", "polyp", "crc"),
                        labels = c("Normal", "Polyp", "CRC")))

p2a <- ggplot(becker_plot, aes(group, score_epi__ca_route_signature, colour = group)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.43,
               colour = COL[["ink"]], fill = "white") +
  geom_point(position = position_jitter(width = 0.12, seed = 20260710),
             size = 0.92, alpha = 0.72) +
  scale_colour_manual(values = c(Normal = COL[["neutral"]], Polyp = COL[["adenoma"]], CRC = COL[["crc"]])) +
  labs(title = "Becker epithelial transfer", subtitle = "Points are samples; inference is patient clustered",
       x = NULL, y = "Locked route score") +
  theme_jtm() +
  theme(legend.position = "none")

becker_forest <- becker_models %>%
  filter(outcome == "epi__ca_route_signature",
         term %in% c("disease_stage_group_polyp", "disease_stage_group_crc")) %>%
  mutate(
    comparison = factor(term,
                        levels = c("disease_stage_group_crc", "disease_stage_group_polyp"),
                        labels = c("CRC vs normal", "Polyp vs normal")),
    p_label = paste0("P = ", format_p(p_value))
  )

p2b <- forest_panel(
  becker_forest, coef, ci_low, ci_high, comparison, COL[["route"]],
  "Patient-clustered transfer model", "72 observations from 17 patients",
  "Adjusted coefficient (95% CI)", c(-0.20, 3.25)
) +
  geom_text(aes(x = 2.50, label = p_label), hjust = 0, size = 1.85)

primary_cor <- rna_atac_cor %>%
  filter(subset == "normal_polyp", analysis_id == "rna_epi_ca_route__atac_wnt_stem")
primary_boot <- rna_atac_boot %>% filter(analysis_id == "locked_route__wnt_tss")
paired_plot <- paired %>%
  mutate(group = factor(disease_stage_group, levels = c("normal_unaffected", "polyp"),
                        labels = c("Normal", "Polyp")))

p2c <- ggplot(paired_plot, aes(atac_tss__wnt_stem, rna_epi__ca_route_signature, colour = group)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = COL[["neutral"]], fill = COL[["neutral_light"]], linewidth = 0.55) +
  geom_point(size = 1.35, alpha = 0.86) +
  scale_colour_manual(values = c(Normal = COL[["neutral"]], Polyp = COL[["adenoma"]])) +
  labs(
    title = "Matched RNA-accessibility coupling",
    subtitle = sprintf("40 samples / 12 patients; rho = %.3f; cluster-bootstrap 95%% CI %.3f to %.3f",
                       primary_cor$spearman_rho, primary_boot$bootstrap_ci_low,
                       primary_boot$bootstrap_ci_high),
    x = "WNT/stemness TSS accessibility", y = "Locked epithelial RNA route", colour = NULL
  ) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

model_labels <- c(
  locked_route__wnt_tss = "WNT TSS",
  locked_route__wnt_tss_minus_housekeeping = "WNT − housekeeping",
  locked_route__wnt_tcf_ascl2_axis = "TCF/ASCL2 axis",
  locked_route__wnt_tcf_ascl2_axis_minus_housekeeping = "TCF/ASCL2 − housekeeping"
)
rna_forest <- rna_atac_models %>%
  mutate(label = factor(model_labels[analysis_id], levels = rev(unname(model_labels))))

p2d <- forest_panel(
  rna_forest, coef, ci_low, ci_high, label, COL[["route"]],
  "Patient-aware adjusted models", "Lesion, proliferation and TSS-depth adjusted",
  "RNA-from-ATAC coefficient (95% CI)", c(-0.05, 0.90)
)

fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_layout(widths = c(1.03, 0.97), heights = c(0.90, 1.10)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_source(becker_forest, "fig2_becker_patient_cluster_effects.tsv")
write_source(paired_plot %>% select(multiome_geo_accession, scrna_geo_accession, patient_id, group,
                                    rna_epi__ca_route_signature, atac_tss__wnt_stem),
             "fig2_matched_rna_atac_samples.tsv")
write_source(rna_forest, "fig2_rna_atac_adjusted_effects.tsv")

# -----------------------------------------------------------------------------
# Main Figure 3: cross-sectional recapitulation across CRC Atlas states
# -----------------------------------------------------------------------------

state_levels <- c(
  "normal_epithelial", "polyp_epithelial", "polyp_cancer",
  "tumor_epithelial", "tumor_cancer", "metastasis_epithelial", "metastasis_cancer"
)
state_labels <- c(
  "Normal\nepithelium", "Polyp\nepithelium", "Polyp\ncancer",
  "Primary\nepithelium", "Primary\ncancer", "Metastasis\nepithelium", "Metastasis\ncancer"
)
state_short <- c(
  "Polyp\nepithelium", "Polyp\ncancer", "Primary-tumour\nepithelium",
  "Primary\ncancer", "Metastatic\nepithelium", "Metastatic\ncancer"
)
state_colours <- setNames(c(COL[["neutral"]], COL[["adenoma"]], "#E4BC6A",
                            "#75A9C9", COL[["crc"]], "#A18DB8", COL[["metastasis"]]),
                          state_labels)

atlas_plot <- atlas_scores %>%
  filter(n_cells_sampled >= 20, carrier_group %in% state_levels) %>%
  mutate(state = factor(carrier_group, levels = state_levels, labels = state_labels))
atlas_n <- atlas_plot %>% count(state, name = "n")

p3a <- ggplot(atlas_plot, aes(state, score__ca_route_signature, colour = state)) +
  geom_boxplot(width = 0.54, outlier.shape = NA, linewidth = 0.42,
               colour = COL[["ink"]], fill = "white") +
  geom_point(position = position_jitter(width = 0.13, seed = 20260710),
             alpha = 0.34, size = 0.55) +
  geom_text(data = atlas_n, aes(x = state, y = -Inf, label = paste0("n = ", n)),
            inherit.aes = FALSE, vjust = -0.35, size = 1.75) +
  scale_colour_manual(values = state_colours) +
  labs(
    title = "Locked-route distributions across donor-carrier states",
    subtitle = "Displayed units contain at least 20 sampled epithelial/cancer cells",
    x = NULL, y = "Locked route score"
  ) +
  theme_jtm() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 22, hjust = 1, size = 5.8))

carrier_terms <- paste0("carrier_group_", state_levels[-1])
carrier_models <- atlas_models %>%
  filter(term %in% carrier_terms, outcome %in% c("ca_route_signature", "wnt_stem")) %>%
  mutate(
    state = factor(term, levels = rev(carrier_terms), labels = rev(state_short)),
    programme = factor(outcome, levels = c("ca_route_signature", "wnt_stem"),
                       labels = c("Locked route", "WNT/stemness"))
  )
stopifnot(all(carrier_models$n_observations == 670), all(carrier_models$n_clusters == 398))

p3b <- ggplot(carrier_models, aes(coef, state, colour = programme)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = state), linewidth = 0.62,
               position = position_dodge(width = 0.44)) +
  geom_point(size = 1.80, position = position_dodge(width = 0.44)) +
  scale_colour_manual(values = c("Locked route" = COL[["route"]], "WNT/stemness" = COL[["wnt"]])) +
  labs(
    title = "Adjusted state effects",
    subtitle = "Study, proliferation and cell count adjusted\nDonor-clustered SE",
    x = "Coefficient (95% CI)", y = NULL, colour = NULL
  ) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

influence_plot <- atlas_influence %>%
  mutate(
    state_label = factor(state, levels = rev(state_levels[-1]), labels = rev(state_short)),
    programme = factor(outcome,
                       levels = c("score__ca_route_signature", "score__wnt_stem"),
                       labels = c("Locked route", "WNT/stemness")),
    min_label = paste0("min n=", minimum_target_donors_after_omission)
  )

p3c <- ggplot(influence_plot, aes(full_coef, state_label, colour = programme)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(aes(x = loo_min_coef, xend = loo_max_coef, yend = state_label,
                   group = programme), linewidth = 1.65, alpha = 0.28,
               position = position_dodge(width = 0.46), lineend = "round") +
  geom_segment(aes(x = full_ci_low, xend = full_ci_high, yend = state_label,
                   group = programme), linewidth = 0.58,
               position = position_dodge(width = 0.46)) +
  geom_point(size = 1.65, position = position_dodge(width = 0.46)) +
  scale_colour_manual(values = c("Locked route" = COL[["route"]], "WNT/stemness" = COL[["wnt"]])) +
  labs(
    title = "Study-omission robustness",
    subtitle = "Pale: range across 33 omissions\nAll coefficients remained positive",
    x = "Adjusted coefficient", y = NULL, colour = NULL
  ) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

fig3 <- p3a / (p3b | p3c) +
  plot_layout(heights = c(0.96, 1.04), widths = c(1, 1)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_source(carrier_models, "fig3_adjusted_carrier_effects.tsv")
write_source(influence_plot, "fig3_leave_one_study_out_ranges.tsv")

# -----------------------------------------------------------------------------
# Supplementary Figure S1: signature composition and platform transportability
# -----------------------------------------------------------------------------

signature_display <- signature %>%
  group_by(signature_direction) %>%
  slice_min(rank_within_direction, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    direction = ifelse(signature_direction == "adenoma_up", "Adenoma-up", "Adenoma-down"),
    gene = reorder(gene, discovery_effect_adenoma_minus_normal)
  )

pS1a <- ggplot(signature_display,
               aes(discovery_effect_adenoma_minus_normal, gene, fill = direction)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = COL[["neutral"]]) +
  scale_fill_manual(values = c("Adenoma-up" = COL[["route"]], "Adenoma-down" = COL[["wnt"]])) +
  labs(
    title = "Highest-ranked locked genes",
    subtitle = "Top 12 per direction; ranks were fixed in Chen discovery donors",
    x = "Discovery donor-median effect", y = NULL, fill = NULL
  ) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

coverage_core <- platform_coverage %>%
  filter(dataset %in% c("Chen discovery", "Chen held-out", "Becker", "CRC Atlas")) %>%
  select(dataset, coverage_up, coverage_down)
coverage_external <- external_coverage %>%
  filter(signature_size_per_direction == 50) %>%
  transmute(dataset = cohort, coverage_up, coverage_down)
coverage_plot <- bind_rows(coverage_core, coverage_external) %>%
  pivot_longer(c(coverage_up, coverage_down), names_to = "component", values_to = "coverage") %>%
  mutate(
    component = factor(component, levels = c("coverage_up", "coverage_down"),
                       labels = c("Up", "Down")),
    dataset = factor(dataset, levels = rev(c("Chen discovery", "Chen held-out", "Becker", "CRC Atlas",
                                             "GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820"))),
    label = percent(coverage, accuracy = 1)
  )

pS1b <- ggplot(coverage_plot, aes(component, dataset, fill = coverage)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = label), size = 2.0, colour = COL[["ink"]]) +
  scale_fill_gradient(low = "#F3E4DF", high = COL[["route"]], limits = c(0.80, 1.00),
                      breaks = c(0.80, 0.90, 1.00),
                      labels = percent_format(accuracy = 1), name = "Coverage") +
  labs(
    title = "Platform coverage",
    subtitle = "All components met the 80% threshold",
    x = NULL, y = NULL
  ) +
  theme_jtm() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "top")

figS1 <- (pS1a | pS1b) +
  plot_layout(widths = c(1.08, 0.92)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_source(signature_display, "figureS1_top_locked_genes.tsv")
write_source(coverage_plot, "figureS1_platform_coverage.tsv")

# -----------------------------------------------------------------------------
# Supplementary Figure S2: external-cohort sensitivity and specificity
# -----------------------------------------------------------------------------

nested <- external_models %>%
  filter(excluded_cohort == "__NONE__") %>%
  mutate(size = factor(paste0(signature_size_per_direction, " + ", signature_size_per_direction),
                       levels = paste0(c(10, 20, 30, 50), " + ", c(10, 20, 30, 50))))

pS2a <- forest_panel(
  nested, adenoma_coef_sd, ci_low, ci_high, size, COL[["route"]],
  "Nested signature-size sensitivity", "203 samples from 161 patient clusters",
  "Standardised adenoma effect", c(1.55, 2.20)
)

loo_compare <- bind_rows(
  external_loo %>% transmute(excluded_cohort, estimate = adenoma_coef_sd,
                             ci_low, ci_high, model = "Primary"),
  external_prolif_loo %>% transmute(excluded_cohort, estimate = adenoma_coef_sd,
                                    ci_low, ci_high, model = "+ proliferation")
) %>%
  mutate(
    cohort = factor(excluded_cohort, levels = rev(c("GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820"))),
    model = factor(model, levels = c("Primary", "+ proliferation"))
  )

pS2b <- ggplot(loo_compare, aes(estimate, cohort, colour = model)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = cohort, group = model),
               position = position_dodge(width = 0.46), linewidth = 0.60) +
  geom_point(position = position_dodge(width = 0.46), size = 1.70) +
  scale_colour_manual(values = c(Primary = COL[["route"]], `+ proliferation` = COL[["wnt"]])) +
  labs(
    title = "Leave-one-cohort-out sensitivity",
    subtitle = "Labels identify the omitted cohort",
    x = "Standardised adenoma effect", y = NULL, colour = NULL
  ) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

grade_scores <- external_scores %>%
  filter(cohort == "GSE41657", grade_group %in% c("normal", "low_grade", "high_grade")) %>%
  mutate(grade = factor(grade_group, levels = c("normal", "low_grade", "high_grade"),
                        labels = c("Normal\n(n=12)", "Low grade\n(n=21)", "High grade\n(n=30)")))
grade_test <- external_grade %>% filter(comparison == "high_grade_vs_low_grade")

pS2c <- ggplot(grade_scores, aes(grade, route_score_k50, colour = grade)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, colour = COL[["ink"]], fill = "white", linewidth = 0.42) +
  geom_point(position = position_jitter(width = 0.10, seed = 20260710), size = 0.90, alpha = 0.72) +
  scale_colour_manual(values = c("Normal\n(n=12)" = COL[["neutral"]],
                                 "Low grade\n(n=21)" = COL[["adenoma"]],
                                 "High grade\n(n=30)" = COL[["route"]])) +
  labs(
    title = "Dysplasia-grade sensitivity",
    subtitle = sprintf("High vs low: %.2f SD; 95%% CI %.2f to %.2f; BH q = %s",
                       grade_test$clustered_standardized_mean_difference,
                       grade_test$clustered_standardized_ci_low,
                       grade_test$clustered_standardized_ci_high,
                       format_p(grade_test$q_value_bh)),
    x = NULL, y = "Locked route score"
  ) +
  theme_jtm() +
  theme(legend.position = "none")

hist_scores <- external_scores %>%
  filter(cohort == "GSE40362", tissue_group %in% c("normal", "hyperplastic", "adenoma")) %>%
  mutate(histology = factor(tissue_group, levels = c("normal", "hyperplastic", "adenoma"),
                            labels = c("Normal\n(n=8)", "Hyperplastic\n(n=8)", "Adenoma\n(n=8)")))
hyper_normal <- external_tests %>%
  filter(cohort == "GSE40362", comparison == "hyperplastic_vs_normal", signature_size_per_direction == 50)
adenoma_hyper <- external_tests %>%
  filter(cohort == "GSE40362", comparison == "adenoma_vs_hyperplastic", signature_size_per_direction == 50)

pS2d <- ggplot(hist_scores, aes(histology, route_score_k50, colour = histology)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, colour = COL[["ink"]], fill = "white", linewidth = 0.42) +
  geom_point(position = position_jitter(width = 0.10, seed = 20260710), size = 1.0, alpha = 0.76) +
  scale_colour_manual(values = c("Normal\n(n=8)" = COL[["neutral"]],
                                 "Hyperplastic\n(n=8)" = COL[["wnt_light"]],
                                 "Adenoma\n(n=8)" = COL[["route"]])) +
  labs(
    title = "Histology specificity",
    subtitle = sprintf("Hyperplastic vs normal P = %s\nAdenoma vs hyperplastic P = %s",
                       format_p(hyper_normal$primary_p_value),
                       format_p(adenoma_hyper$primary_p_value)),
    x = NULL, y = "Locked route score"
  ) +
  theme_jtm() +
  theme(legend.position = "none")

figS2 <- (pS2a | pS2b) / (pS2c | pS2d) +
  plot_layout(widths = c(0.92, 1.08), heights = c(0.92, 1.08)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_source(nested, "figureS2_nested_signature_effects.tsv")
write_source(loo_compare, "figureS2_leave_one_cohort_out.tsv")
write_source(grade_scores, "figureS2_grade_scores.tsv")
write_source(hist_scores, "figureS2_histology_scores.tsv")

# -----------------------------------------------------------------------------
# Supplementary Figure S3: patient-aware RNA-ATAC sensitivities
# -----------------------------------------------------------------------------

patient_medians <- paired_plot %>%
  group_by(patient_id) %>%
  summarise(
    rna_route = median(rna_epi__ca_route_signature, na.rm = TRUE),
    atac_wnt = median(atac_tss__wnt_stem, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )
primary_median <- rna_atac_medians %>% filter(analysis_id == "locked_route__wnt_tss")

pS3a <- ggplot(patient_medians, aes(atac_wnt, rna_route)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = COL[["neutral"]], fill = COL[["neutral_light"]], linewidth = 0.55) +
  geom_point(aes(size = n_samples), colour = COL[["route"]], alpha = 0.88) +
  scale_size_continuous(range = c(1.2, 2.5)) +
  labs(
    title = "Patient-median coupling",
    subtitle = sprintf("12 patients; rho = %.3f; 100,000-permutation P = %s",
                       primary_median$spearman_rho,
                       format_p(primary_median$p_value_patient_permutation)),
    x = "Patient-median WNT/stemness TSS accessibility",
    y = "Patient-median locked RNA route", size = "Samples"
  ) +
  theme_jtm() +
  theme(legend.position = "top")

model_sensitivity <- bind_rows(
  rna_atac_models %>% transmute(analysis_id, coef, ci_low, ci_high, inference = "Patient clustered"),
  rna_atac_fixed %>% transmute(analysis_id, coef, ci_low, ci_high, inference = "Patient fixed effect")
) %>%
  mutate(
    label = factor(model_labels[analysis_id], levels = rev(unname(model_labels))),
    inference = factor(inference, levels = c("Patient clustered", "Patient fixed effect"))
  )

pS3b <- ggplot(model_sensitivity, aes(coef, label, colour = inference)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label, group = inference),
               linewidth = 0.58, position = position_dodge(width = 0.46)) +
  geom_point(size = 1.65, position = position_dodge(width = 0.46)) +
  scale_colour_manual(values = c("Patient clustered" = COL[["route"]],
                                 "Patient fixed effect" = COL[["wnt"]])) +
  labs(
    title = "Patient-model sensitivity",
    subtitle = "Clustered versus fixed-effect estimates",
    x = "RNA-from-ATAC coefficient", y = NULL, colour = NULL
  ) +
  theme_jtm() +
  theme(legend.position = "top", legend.justification = "left")

window_plot <- window_cor %>%
  filter(subset == "normal_polyp") %>%
  mutate(
    distance = factor(distance_bin,
                      levels = c("tss_core_1kb", "promoter_proximal_2_5kb",
                                 "proximal_regulatory_10kb", "distal_flank_20kb"),
                      labels = c("TSS 1 kb", "2-5 kb", "10 kb", "20 kb")),
    feature = factor(rna_feature,
                     levels = c("rna_epi__ca_route_signature", "rna_epi__wnt_stem",
                                "rna_epi__proliferation_control"),
                     labels = c("Locked route", "WNT/stemness", "Proliferation")),
    locus = factor(locus_set,
                   levels = c("wnt_route_loci", "wnt_tcf_ascl2_axis_loci"),
                   labels = c("WNT-route loci", "TCF/ASCL2-axis loci"))
  )

pS3c <- ggplot(window_plot, aes(distance, spearman_rho, colour = feature, group = feature)) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = COL[["neutral"]]) +
  geom_line(linewidth = 0.58) +
  geom_point(size = 1.35) +
  facet_wrap(~locus, nrow = 1) +
  scale_colour_manual(values = c("Locked route" = COL[["route"]],
                                 "WNT/stemness" = COL[["wnt"]],
                                 "Proliferation" = COL[["neutral"]])) +
  coord_cartesian(ylim = c(-0.32, 0.90)) +
  labs(
    title = "Regulatory-distance sensitivity",
    subtitle = "Sample-level Spearman correlations are descriptive",
    x = "Accessibility window", y = "Spearman rho", colour = NULL
  ) +
  theme_jtm() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 22, hjust = 1))

figS3 <- (pS3a | pS3b) / pS3c +
  plot_layout(heights = c(1.03, 0.97)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_source(patient_medians, "figureS3_patient_medians.tsv")
write_source(model_sensitivity, "figureS3_patient_model_sensitivity.tsv")
write_source(window_plot, "figureS3_regulatory_window_correlations.tsv")

# -----------------------------------------------------------------------------
# Supplementary Figure S4: CRC Atlas source-study support
# -----------------------------------------------------------------------------

informative_studies <- atlas_support %>%
  filter(carrier_group != "normal_epithelial", n_donors > 0) %>%
  distinct(study_id) %>%
  pull(study_id)
support_grid <- expand_grid(study_id = sort(informative_studies), carrier_group = state_levels) %>%
  left_join(atlas_support %>% select(study_id, carrier_group, n_donors),
            by = c("study_id", "carrier_group")) %>%
  mutate(
    n_donors = replace_na(n_donors, 0),
    state = factor(
      carrier_group,
      levels = state_levels,
      labels = c(
        "Normal\nepithelium", "Polyp\nepithelium", "Polyp\ncancer",
        "Primary-tumour\nepithelium", "Primary\ncancer",
        "Metastatic\nepithelium", "Metastatic\ncancer"
      )
    ),
    study_label = factor(gsub("_", " ", study_id), levels = rev(gsub("_", " ", sort(informative_studies)))),
    fill_value = ifelse(n_donors == 0, NA_real_, log10(n_donors + 1)),
    donor_label = ifelse(n_donors > 0, as.character(n_donors), "")
  )

pS4a <- ggplot(support_grid, aes(state, study_label, fill = fill_value)) +
  geom_tile(colour = "white", linewidth = 0.32) +
  geom_text(aes(label = donor_label), size = 1.75, colour = COL[["ink"]]) +
  scale_fill_gradient(low = "#E4EEF4", high = COL[["wnt"]],
                      breaks = log10(c(2, 6, 21, 81)), labels = c("1", "5", "20", "80"),
                      na.value = COL[["neutral_pale"]], name = "Donors") +
  labs(title = "Atlas source-study support", subtitle = "Donor-carrier units; grey = none",
       x = NULL, y = NULL) +
  theme_jtm(base_size = 7.0) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 5.2),
    axis.text.y = element_text(size = 5.2), legend.position = "top",
    panel.background = element_rect(fill = COL[["neutral_pale"]], colour = NA)
  )

within_primary <- atlas_within %>%
  filter(outcome == "score__ca_route_signature", state == "tumor_cancer") %>%
  mutate(
    study = factor(gsub("_", " ", study_id), levels = rev(gsub("_", " ", study_id))),
    label = paste0("n = ", n_target_donors, "/", n_reference_donors)
  )

pS4b <- ggplot(within_primary, aes(coef, study)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = study), linewidth = 0.58, colour = COL[["ink"]]) +
  geom_point(size = 1.75, colour = COL[["route"]]) +
  geom_text(aes(x = max(ci_high) + 0.14, label = label), hjust = 0, size = 1.75, colour = COL[["neutral"]]) +
  coord_cartesian(xlim = c(min(within_primary$ci_low) - 0.05, max(within_primary$ci_high) + 0.72), clip = "off") +
  labs(
    title = "Within-study effects",
    subtitle = "Primary cancer; target/reference donors",
    x = "Primary cancer vs normal coefficient", y = NULL
  ) +
  theme_jtm(base_size = 7.0)

figS4 <- (pS4a | pS4b) +
  plot_layout(widths = c(1.18, 0.82)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_source(support_grid, "figureS4_study_support.tsv")
write_source(within_primary, "figureS4_within_study_primary_cancer.tsv")

# -----------------------------------------------------------------------------
# Provisional Figure 4 blueprint: no simulated data or outcome claims
# -----------------------------------------------------------------------------

image_slots <- expand_grid(
  tissue = factor(c("Normal / adjacent", "Conventional adenoma"),
                  levels = c("Normal / adjacent", "Conventional adenoma")),
  marker = factor(c("OLFM4", "SOX9", "beta-catenin", "FABP1", "Ki-67"),
                  levels = c("OLFM4", "SOX9", "beta-catenin", "FABP1", "Ki-67"))
) %>%
  mutate(x = as.numeric(marker), y = 3 - as.numeric(tissue))

p4a <- ggplot(image_slots, aes(x, y)) +
  geom_tile(width = 0.92, height = 0.82, fill = "#ECEFF1", colour = "#A7AFB6", linewidth = 0.35) +
  geom_text(aes(label = "Observed\nimage"), size = 1.8, colour = COL[["neutral"]], lineheight = 0.9) +
  scale_x_continuous(breaks = 1:5, labels = levels(image_slots$marker), position = "top") +
  scale_y_continuous(breaks = 1:2, labels = rev(levels(image_slots$tissue))) +
  coord_fixed(ratio = 0.62, clip = "off") +
  labs(title = "Representative serial-section IHC", subtitle = "Insert calibrated, blinded image crops with scale bars",
       x = NULL, y = NULL) +
  theme_void(base_size = 7, base_family = JTM_FONT) +
  theme(
    axis.text.x = element_text(size = 5.7, face = "bold", margin = margin(b = 1.2, unit = "mm")),
    axis.text.y = element_text(size = 5.7, face = "bold", hjust = 1),
    plot.title = element_text(size = 7.5, face = "bold"),
    plot.subtitle = element_text(size = 6.4, colour = "#4C5258"),
    plot.margin = margin(3, 3, 3, 3, unit = "mm")
  )

placeholder_panel <- function(title, subtitle, xlab, ylab = NULL) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "#F6F7F8", colour = "#B8BEC4", linewidth = 0.45) +
    annotate("text", x = 0.5, y = 0.55, label = "AWAITING\nOBSERVED IHC DATA",
             size = 2.25, fontface = "bold", colour = COL[["neutral"]], lineheight = 0.92) +
    annotate("text", x = 0.5, y = 0.28, label = subtitle,
             size = 1.72, colour = "#596169", lineheight = 0.95) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    labs(title = title, x = xlab, y = ylab) +
    theme_jtm() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
}

p4b <- placeholder_panel("Prespecified epithelial composite",
                         "Normal/adjacent vs adenoma\npatient-level points and 95% CI",
                         "Tissue group", "Composite score")
p4c <- placeholder_panel("Ki-67- and batch-adjusted effects",
                         "Composite plus four BH-adjusted\nmarker-specific estimates",
                         "Adjusted coefficient", NULL)
p4d <- placeholder_panel("Composite versus Ki-67",
                         "Show observed specimens\nwithout endpoint retuning",
                         "Epithelial Ki-67", "Composite score")
p4e <- placeholder_panel("Scoring reliability",
                         "Absolute-agreement ICC\nwith bootstrap 95% CI",
                         "ICC", NULL)

fig4_blueprint <- p4a / (p4b | p4c) / (p4d | p4e) +
  plot_layout(heights = c(1.18, 0.86, 0.86)) +
  plot_annotation(
    title = "PROVISIONAL FIGURE 4 BLUEPRINT — NOT FOR SUBMISSION",
    subtitle = "No simulated measurements or outcomes are shown.",
    tag_levels = "a",
    theme = tag_theme + theme(
      plot.title = element_text(size = 9, face = "bold", colour = COL[["route"]]),
      plot.subtitle = element_text(size = 7, colour = COL[["neutral"]])
    )
  )

# -----------------------------------------------------------------------------
# Export, traceability, and QA
# -----------------------------------------------------------------------------

export_parts <- list(
  export_figure(fig1, "figure1_mainline_discovery_external_validation", 170, 164),
  export_figure(fig2, "figure2_mainline_becker_rna_atac", 170, 122),
  export_figure(fig3, "figure3_mainline_crc_atlas_recapitulation", 170, 126),
  export_figure(figS1, "figureS1_signature_and_platform_coverage", 170, 112),
  export_figure(figS2, "figureS2_external_validation_sensitivity", 170, 132),
  export_figure(figS3, "figureS3_rna_atac_sensitivity", 170, 124),
  export_figure(figS4, "figureS4_atlas_source_support", 170, 140)
)
if (identical(Sys.getenv("JTM_INCLUDE_PROVISIONAL", unset = "0"), "1")) {
  export_parts <- append(
    export_parts,
    list(export_figure(
      fig4_blueprint,
      "figure4_PROVISIONAL_ihc_blueprint_NOT_FOR_SUBMISSION",
      170,
      174
    ))
  )
}
export_manifest <- bind_rows(export_parts)
export_manifest$file <- normalizePath(export_manifest$file, mustWork = TRUE)
export_manifest$file_size_bytes <- file.info(export_manifest$file)$size
export_manifest$sha256 <- vapply(
  export_manifest$file,
  function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), "\\s+")[[1]][1],
  character(1)
)
write.table(export_manifest, file.path(OUT_DIR, "figure_export_manifest.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

panel_trace <- tribble(
  ~figure, ~panel, ~claim, ~source_table, ~independent_unit, ~inference,
  "Figure 1", "a", "Validation labels were unavailable during discovery locking", "locked_figure_source_manifest.tsv", "workflow stage", "design schematic",
  "Figure 1", "b", "One hundred genes met the frozen stability rule", "fig1_discovery_stability_audit.tsv; fig1_locked_signature_genes.tsv", "gene", "1,000 donor-cluster bootstraps",
  "Figure 1", "c", "The frozen score separates held-out adenoma from normal tissue", "fig1_chen_locked_scores.tsv; fig1_chen_locked_discrimination.tsv", "specimen", "Mann-Whitney; rank AUC",
  "Figure 1", "d", "The score transports across five independent sporadic cohorts", "figs3_external_cohort_tests.tsv", "sample; patient cluster", "paired Wilcoxon or patient-clustered OLS",
  "Figure 1", "e", "The held-out direction is present within donors", "fig1_chen_locked_scores.tsv; fig1_chen_locked_paired_tests.tsv", "donor pair", "paired Wilcoxon",
  "Figure 1", "f", "The pooled effect survives proliferation adjustment and cohort exclusion", "figs3_external_one_stage_models.tsv; figs3_external_leave_one_cohort_out.tsv; figs3_external_proliferation_adjusted_model.tsv; figs3_external_proliferation_adjusted_leave_one_cohort_out.tsv", "sample; patient cluster", "cohort-fixed OLS with patient-clustered SE",
  "Figure 2", "a", "The locked route transfers across Becker tissue states", "fig1_becker_locked_scores.tsv", "sample", "descriptive distribution",
  "Figure 2", "b", "Becker transfer remains positive under patient-clustered inference", "fig1_becker_patient_cluster_models.tsv", "sample; patient cluster", "adjusted OLS with patient-clustered SE",
  "Figure 2", "c", "Matched RNA route activity covaries with WNT/stemness TSS accessibility", "fig2_locked_rna_atac_paired_scores.tsv; fig2_locked_rna_atac_patient_cluster_bootstrap.tsv", "matched sample; patient cluster", "Spearman effect with patient-cluster bootstrap CI",
  "Figure 2", "d", "RNA-accessibility effects survive prespecified covariate adjustment", "fig2_locked_rna_atac_patient_cluster_models.tsv", "matched sample; patient cluster", "adjusted OLS with patient-clustered SE",
  "Figure 3", "a", "The locked route is detectable across CRC Atlas states", "fig3_atlas_locked_donor_scores.tsv", "donor-carrier", "descriptive distribution",
  "Figure 3", "b", "Locked-route and curated WNT/stemness effects are positive after adjustment", "fig3_atlas_donor_cluster_models.tsv", "donor-carrier; donor cluster", "study-adjusted OLS with donor-clustered SE",
  "Figure 3", "c", "All state coefficients remain positive across 33 study omissions", "figs2_atlas_study_influence.tsv", "donor-carrier; donor cluster", "leave-one-study-out adjusted OLS",
  "Figure S1", "a-b", "Signature composition and platform coverage", "fig1_locked_signature_genes.tsv; supp_locked_signature_platform_coverage.tsv; figs3_external_signature_coverage.tsv", "gene and platform", "descriptive audit",
  "Figure S2", "a-d", "External-cohort size, exclusion, grade and histology sensitivities", "figs3_external_one_stage_models.tsv; figs3_external_leave_one_cohort_out.tsv; figs3_external_sample_scores.tsv; figs3_external_grade_sensitivity.tsv; figs3_external_cohort_tests.tsv", "sample; patient cluster", "patient-clustered or paired inference as specified",
  "Figure S3", "a-c", "Patient aggregation, fixed effects and distance-window sensitivities", "fig2_locked_rna_atac_patient_median_correlations.tsv; fig2_locked_rna_atac_patient_fixed_effect_models.tsv; fig2_locked_regulatory_window_correlations.tsv", "patient or matched sample", "permutation, adjusted model, descriptive Spearman",
  "Figure S4", "a-b", "Atlas source composition and within-study support", "figs2_atlas_state_study_support.tsv; figs2_atlas_within_study_contrasts.tsv", "donor-carrier", "descriptive support and donor-clustered within-study contrasts"
)
write.table(panel_trace, file.path(OUT_DIR, "panel_source_trace.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

qa <- data.frame(
  check = c(
    "locked_signature_gene_count",
    "validation_unseen_during_selection",
    "chen_heldout_auc",
    "external_cohort_count",
    "external_primary_samples",
    "external_patient_clusters",
    "external_one_stage_effect",
    "external_proliferation_adjusted_effect",
    "becker_observations",
    "becker_patients",
    "rna_atac_samples",
    "rna_atac_patients",
    "atlas_observations",
    "atlas_donors",
    "atlas_loo_all_positive",
    "export_count",
    "all_exports_nonempty",
    "all_exports_under_10mb"
  ),
  observed = c(
    nrow(signature),
    all(tolower(as.character(signature$validation_used_for_selection)) == "false"),
    chen_discrimination$auc_adenoma_vs_normal[chen_discrimination$dataset == "validation"],
    nrow(external_primary),
    primary_model$n_samples,
    primary_model$n_patient_clusters,
    primary_model$adenoma_coef_sd,
    primary_prolif$adenoma_coef_sd,
    unique(becker_forest$n_observations),
    unique(becker_forest$n_clusters),
    nrow(paired),
    length(unique(paired$patient_id)),
    unique(carrier_models$n_observations),
    unique(carrier_models$n_clusters),
    all(atlas_influence$loo_positive_fraction == 1),
    nrow(export_manifest),
    all(export_manifest$file_size_bytes > 0),
    all(export_manifest$file_size_bytes < 10 * 1024^2)
  ),
  expected = c(100, TRUE, 0.9285714, 5, 203, 161, 1.924521013305084,
               1.7410914795380723, 72, 17, 40, 12, 670, 398, TRUE, 28, TRUE, TRUE),
  stringsAsFactors = FALSE
)
qa$pass <- mapply(function(observed, expected) {
  if (expected %in% c("TRUE", "FALSE")) return(observed == expected)
  abs(as.numeric(observed) - as.numeric(expected)) < 1e-6
}, qa$observed, qa$expected)
write.table(qa, file.path(OUT_DIR, "figure_build_qc.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
if (!all(qa$pass)) stop("One or more JTM mainline figure QA checks failed.")

message("JTM mainline figure set written to: ", OUT_DIR)
