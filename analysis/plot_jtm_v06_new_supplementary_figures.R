#!/usr/bin/env Rscript

# Render the two evidence-gap figures identified in the v0.6 mock review.
# This script reads existing result tables and does not refit source analyses.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[1])
} else {
  file.path(getwd(), "analysis", "plot_jtm_v06_new_supplementary_figures.R")
}

ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT_DIR <- Sys.getenv(
  "JTM_FIGURE_DIR",
  unset = file.path(ROOT, "figures", "jtm_deep_v0.6")
)
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(
    file.path(ROOT, path),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

write_tsv <- function(frame, filename) {
  write.table(
    frame,
    file.path(SOURCE_DIR, filename),
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
  teal = "#2F8F83",
  yellow = "#D9A441",
  purple = "#8064A2"
)

theme_jtm <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = JTM_FONT) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = COL[["ink"]]),
      axis.ticks = element_line(linewidth = 0.30, colour = COL[["ink"]]),
      axis.ticks.length = grid::unit(1.1, "mm"),
      axis.title = element_text(size = base_size, colour = COL[["ink"]]),
      axis.text = element_text(size = base_size - 0.5, colour = COL[["ink"]]),
      panel.grid = element_blank(),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.7),
      legend.key.height = grid::unit(3, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.4, face = "bold"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(2.2, 2.4, 2.2, 2.4, unit = "mm"),
      plot.tag = element_text(size = 9, face = "bold", colour = COL[["ink"]])
    )
}

tag_theme <- theme(
  plot.tag.position = c(0, 1),
  plot.tag = element_text(size = 9, face = "bold", colour = COL[["ink"]])
)

export_figure <- function(plot, stem, width_mm, height_mm) {
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
    file = normalizePath(unname(paths), mustWork = TRUE),
    width_mm = width_mm,
    height_mm = height_mm,
    resolution_dpi = c(NA, NA, 600, 300),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Supplementary Figure S5: paired FFPE and expanded transcriptomic validation
# -----------------------------------------------------------------------------

gse_scores <- read_tsv(
  "results/gse117606_paired_route_validation/sample_scores_all_cohort_scaling.tsv"
)
gse_tests <- read_tsv(
  "results/gse117606_paired_route_validation/paired_route_tests.tsv"
)
gse_boot <- read_tsv(
  "results/gse117606_paired_route_validation/paired_route_bootstrap.tsv"
)
gse_adjusted <- read_tsv(
  "results/gse117606_paired_route_validation/paired_proliferation_adjusted_model.tsv"
)
expanded_tests <- read_tsv(
  "results/expanded_public_adenoma_validation/route_tests.tsv"
)

paired_wide <- gse_scores %>%
  filter(tissue_group %in% c("normal", "adenoma")) %>%
  group_by(patient_id, tissue_group) %>%
  summarise(route_score = median(route_score_k50, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = tissue_group, values_from = route_score) %>%
  filter(is.finite(normal), is.finite(adenoma)) %>%
  mutate(delta = adenoma - normal)

paired_long <- paired_wide %>%
  select(patient_id, normal, adenoma) %>%
  pivot_longer(c(normal, adenoma), names_to = "tissue_group", values_to = "route_score") %>%
  mutate(
    tissue_group = factor(
      tissue_group,
      levels = c("normal", "adenoma"),
    labels = c("Reference", "Adenoma")
    )
  )

gse_primary <- gse_tests %>%
  filter(signature_size_per_direction == 50, normalization_scope == "all_cohort_samples")
stopifnot(nrow(paired_wide) == 51, nrow(gse_primary) == 1)

paired_medians <- paired_long %>%
  group_by(tissue_group) %>%
  summarise(route_score = median(route_score), .groups = "drop")

pS5a <- ggplot(paired_long, aes(tissue_group, route_score, group = patient_id)) +
  geom_line(linewidth = 0.30, alpha = 0.24, colour = COL[["neutral"]]) +
  geom_point(aes(colour = tissue_group), size = 0.85, alpha = 0.72) +
  geom_point(
    data = paired_medians,
    aes(tissue_group, route_score),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.5,
    stroke = 0.55,
    fill = "white",
    colour = COL[["ink"]]
  ) +
  annotate(
    "label",
    x = 1.5,
    y = max(paired_long$route_score) + 0.35,
    label = sprintf(
      "%d/%d increased\nmedian difference %.3f; P = %s",
      sum(paired_wide$delta > 0),
      nrow(paired_wide),
      gse_primary$median_paired_difference,
      format_p(gse_primary$p_paired_wilcoxon)
    ),
    size = 2.05,
    linewidth = 0.22,
    fill = "white",
    colour = COL[["ink"]]
  ) +
  scale_colour_manual(values = c(
    "Reference" = COL[["wnt"]],
    "Adenoma" = COL[["route"]]
  )) +
  scale_x_discrete(expand = expansion(mult = c(0.10, 0.10))) +
  coord_cartesian(clip = "off") +
  labs(
    tag = "a",
    title = "Patient-paired FFPE validation",
    subtitle = "GSE117606; 51 patient pairs",
    x = NULL,
    y = "Frozen route score"
  ) +
  theme_jtm() +
  theme(legend.position = "none")

robustness <- bind_rows(
  data.frame(
    model = "Unadjusted\nmedian paired difference",
    estimate = gse_boot$median_delta_bootstrap_median,
    ci_low = gse_boot$median_delta_ci_low,
    ci_high = gse_boot$median_delta_ci_high,
    interval = "Patient-pair bootstrap 95% CI"
  ),
  data.frame(
    model = "Proliferation-adjusted\nmean difference",
    estimate = gse_adjusted$adjusted_mean_adenoma_minus_normal,
    ci_low = gse_adjusted$hc3_ci_low,
    ci_high = gse_adjusted$hc3_ci_high,
    interval = "HC3 95% CI"
  )
) %>%
  mutate(model = factor(model, levels = rev(model)))

pS5b <- ggplot(robustness, aes(estimate, model, colour = interval)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = model), linewidth = 0.65) +
  geom_point(size = 2.0) +
  geom_text(
    aes(label = sprintf("%.2f [%.2f, %.2f]", estimate, ci_low, ci_high)),
    nudge_y = 0.18,
    size = 1.85,
    colour = COL[["ink"]]
  ) +
  scale_colour_manual(values = c(
    "Patient-pair bootstrap 95% CI" = COL[["route"]],
    "HC3 95% CI" = COL[["wnt"]]
  ), labels = c(
    "Patient-pair bootstrap 95% CI" = "Pair-bootstrap CI",
    "HC3 95% CI" = "HC3 CI"
  )) +
  coord_cartesian(xlim = c(-0.10, 2.35), clip = "off") +
  labs(
    tag = "b",
    title = "Proliferation robustness",
    subtitle = "Distinct estimands and interval methods",
    x = "Adenoma-minus-reference route change",
    y = NULL,
    colour = NULL
  ) +
  theme_jtm(base_size = 6.6) +
  theme(legend.position = "top", legend.justification = "left")

route_labels <- c(
  GSE100179 = "GSE100179",
  GSE37364 = "GSE37364",
  GSE164541_paired = "GSE164541 (paired)",
  GSE20916_micro_crypt = "Microdissected crypt",
  GSE20916_micro_mucosa = "Microdissected mucosa",
  GSE20916_macro = "Macrodissected",
  GSE226739 = "GSE226739"
)
cluster_labels <- c(
  SOTE_Budapest = "SOTE",
  GSE164541_paired_set = "GSE164541",
  GSE20916_Warsaw = "Warsaw",
  Longhua_Shanghai = "Longhua"
)

expanded_plot <- expanded_tests %>%
  filter(tolower(as.character(primary_inference_eligible)) == "true") %>%
  mutate(
    route_label = unname(route_labels[cohort]),
    cluster_label = factor(
      unname(cluster_labels[recruitment_cluster]),
      levels = c("SOTE", "GSE164541", "Warsaw", "Longhua")
    ),
    route_label = factor(route_label, levels = rev(unname(route_labels))),
    near_null = cohort == "GSE226739"
  )
stopifnot(nrow(expanded_plot) == 7, n_distinct(expanded_plot$cluster_label) == 4)

pS5c <- ggplot(expanded_plot, aes(clustered_mean_difference, route_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_segment(
    aes(x = clustered_ci_low, xend = clustered_ci_high, yend = route_label, colour = near_null),
    linewidth = 0.62
  ) +
  geom_point(aes(colour = near_null), size = 1.8) +
  facet_grid(
    rows = vars(cluster_label),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_colour_manual(values = c(`FALSE` = COL[["route"]], `TRUE` = COL[["neutral"]])) +
  coord_cartesian(xlim = c(-0.55, 3.45), clip = "off") +
  labs(
    tag = "c",
    title = "Expanded cohort sensitivity",
    subtitle = "7 routes / 4 clusters; clustered 95% CIs",
    x = "Adenoma-minus-normal route effect",
    y = NULL
  ) +
  theme_jtm(base_size = 6.2) +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 5.3),
    axis.text.y = element_text(size = 5.6),
    panel.spacing.y = grid::unit(1.1, "mm")
  )

figS5 <- pS5a | (pS5b / pS5c) +
  plot_layout(widths = c(0.92, 1.08), heights = c(0.38, 0.62)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_tsv(paired_long, "figureS5a_gse117606_patient_tissue_scores.tsv")
write_tsv(robustness, "figureS5b_gse117606_robustness_estimates.tsv")
write_tsv(expanded_plot, "figureS5c_expanded_route_effects.tsv")

# -----------------------------------------------------------------------------
# Supplementary Figure S6: public protein triangulation
# -----------------------------------------------------------------------------

pxd2137 <- read_tsv(
  "results/public_adenoma_protein_triangulation/pxd002137_candidate_tests.tsv"
)
pxd445 <- read_tsv(
  "results/pxd000445_candidate_reanalysis/psm_candidate_paired_tests.tsv"
)
pxd17269 <- read_tsv(
  "results/public_adenoma_protein_triangulation/pxd017269_ffpe_detectability.tsv"
)
pxd46999 <- read_tsv(
  "results/public_adenoma_protein_triangulation/pxd046999_dvp_presence.tsv"
)

marker_order <- c("OLFM4", "FABP1", "SOX9", "CTNNB1", "ETHE1")
marker_labels <- c(
  OLFM4 = "OLFM4 · locked up",
  FABP1 = "FABP1 · locked down",
  SOX9 = "SOX9 · context",
  CTNNB1 = "CTNNB1 · context",
  ETHE1 = "ETHE1 · reserve down"
)

pxd2137_plot <- pxd2137 %>%
  filter(gene %in% marker_order) %>%
  mutate(
    marker = factor(
      unname(marker_labels[gene]),
      levels = rev(unname(marker_labels[marker_order]))
    ),
    role = case_when(
      gene %in% c("OLFM4", "FABP1") ~ "Direct locked arm",
      gene %in% c("SOX9", "CTNNB1") ~ "Context marker",
      TRUE ~ "Reserve locked marker"
    ),
    q_label = ifelse(
      is.finite(age_sex_adjusted_q_bh_candidates),
      paste0("q = ", format_p(age_sex_adjusted_q_bh_candidates)),
      "adjusted estimate\nnot estimable"
    )
  )
stopifnot(nrow(pxd2137_plot) == 5)

pS6a <- ggplot(pxd2137_plot, aes(y = marker, colour = role)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_blank(aes(x = 0)) +
  geom_segment(
    data = filter(pxd2137_plot, is.finite(age_sex_adjusted_log2_effect)),
    aes(
      x = age_sex_adjusted_ci_low,
      xend = age_sex_adjusted_ci_high,
      yend = marker
    ),
    linewidth = 0.68
  ) +
  geom_point(
    data = filter(pxd2137_plot, is.finite(age_sex_adjusted_log2_effect)),
    aes(x = age_sex_adjusted_log2_effect),
    size = 2.0
  ) +
  geom_label(
    aes(x = 6.2, label = q_label),
    hjust = 1,
    size = 1.85,
    colour = COL[["ink"]],
    fill = "white",
    linewidth = 0,
    label.padding = grid::unit(0.35, "mm")
  ) +
  scale_colour_manual(
    values = c(
      "Direct locked arm" = COL[["route"]],
      "Context marker" = COL[["purple"]],
      "Reserve locked marker" = COL[["wnt"]]
    ),
    labels = c(
      "Direct locked arm" = "Direct",
      "Context marker" = "Context",
      "Reserve locked marker" = "Reserve"
    )
  ) +
  coord_cartesian(xlim = c(-2.0, 6.3), clip = "off") +
  labs(
    tag = "a",
    title = "PXD002137 differential abundance",
    subtitle = "HC3 95% CIs; 8 normal, 16 adenoma",
    x = "Adjusted adenoma-minus-normal log2 effect",
    y = NULL,
    colour = NULL
  ) +
  theme_jtm(base_size = 6.4) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = grid::unit(3.2, "mm")
  )

paired_marker_order <- c("OLFM4", "FABP1", "SORD", "ETHE1")
paired_marker_labels <- c(
  OLFM4 = "OLFM4 · locked up",
  FABP1 = "FABP1 · locked down",
  SORD = "SORD · positive control",
  ETHE1 = "ETHE1 · reserve down"
)

pxd445_plot <- pxd445 %>%
  filter(analysis_set == "author_qc_21_pairs", gene %in% paired_marker_order) %>%
  mutate(
    marker = factor(
      unname(paired_marker_labels[gene]),
      levels = rev(unname(paired_marker_labels[paired_marker_order]))
    ),
    direction_concordant_fraction = ifelse(
      expected_direction > 0,
      paired_positive_fraction,
      1 - paired_positive_fraction
    ),
    label = sprintf(
      "n = %d; %.0f%% concordant; q = %s",
      n_pairs,
      100 * direction_concordant_fraction,
      format_p(q_value_bh_within_analysis_set)
    ),
    fdr_supported = q_value_bh_within_analysis_set < 0.05
  )
stopifnot(nrow(pxd445_plot) == 4)

pS6b <- ggplot(
  pxd445_plot,
  aes(median_adenoma_minus_normal_log2, marker, colour = fdr_supported)
) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = COL[["neutral"]]) +
  geom_point(aes(size = direction_concordant_fraction), alpha = 0.92) +
  geom_text(aes(x = 1.70, label = label), hjust = 1, nudge_y = 0.18,
            size = 1.75, colour = COL[["ink"]]) +
  scale_colour_manual(values = c(`TRUE` = COL[["route"]], `FALSE` = COL[["neutral"]])) +
  scale_size_continuous(range = c(1.4, 2.6), limits = c(0.50, 1.00)) +
  coord_cartesian(xlim = c(-1.15, 1.75), clip = "off") +
  labs(
    tag = "b",
    title = "PXD000445 paired sensitivity",
    subtitle = "Author-QC set; no confidence interval was invented",
    x = "Median paired adenoma-minus-normal log2 effect",
    y = NULL
  ) +
  theme_jtm(base_size = 6.2) +
  theme(legend.position = "none", axis.text.y = element_text(size = 5.7))

detect17269 <- pxd17269 %>%
  filter(gene %in% marker_order) %>%
  select(gene, pxd017269_detection_fraction)
detect46999 <- pxd46999 %>%
  filter(gene %in% marker_order) %>%
  select(gene, detected_in_nine_patient_dvp)

matrix_base <- data.frame(gene = marker_order, stringsAsFactors = FALSE) %>%
  left_join(
    pxd2137_plot %>% select(gene, age_sex_adjusted_q_bh_candidates, age_sex_adjusted_log2_effect),
    by = "gene"
  ) %>%
  left_join(
    pxd445_plot %>% select(gene, q_value_bh_within_analysis_set, direction_matches_prespecified),
    by = "gene"
  ) %>%
  left_join(detect17269, by = "gene") %>%
  left_join(detect46999, by = "gene")

matrix_plot <- bind_rows(
  matrix_base %>% transmute(
    gene,
    evidence_source = "PXD\n002137\ncontrast",
    status = case_when(
      is.finite(age_sex_adjusted_q_bh_candidates) & age_sex_adjusted_q_bh_candidates < 0.05 ~ "FDR-supported difference",
      is.finite(age_sex_adjusted_log2_effect) ~ "Direction only / imprecise",
      TRUE ~ "Not estimable"
    ),
    display = case_when(
      is.finite(age_sex_adjusted_q_bh_candidates) & age_sex_adjusted_q_bh_candidates < 0.05 ~ "FDR",
      is.finite(age_sex_adjusted_log2_effect) ~ "direction\nonly",
      TRUE ~ "not\nestimable"
    )
  ),
  matrix_base %>% transmute(
    gene,
    evidence_source = "PXD\n000445\npaired",
    status = case_when(
      is.finite(q_value_bh_within_analysis_set) & q_value_bh_within_analysis_set < 0.05 ~ "FDR-supported difference",
      tolower(as.character(direction_matches_prespecified)) == "true" ~ "Direction only / imprecise",
      TRUE ~ "Not evaluated"
    ),
    display = case_when(
      is.finite(q_value_bh_within_analysis_set) & q_value_bh_within_analysis_set < 0.05 ~ "FDR",
      tolower(as.character(direction_matches_prespecified)) == "true" ~ "direction\nonly",
      TRUE ~ "not\navailable"
    )
  ),
  matrix_base %>% transmute(
    gene,
    evidence_source = "PXD\n017269\nFFPE",
    status = ifelse(is.finite(pxd017269_detection_fraction), "Detectability only", "Not detected"),
    display = ifelse(
      is.finite(pxd017269_detection_fraction),
      percent(pxd017269_detection_fraction, accuracy = 1),
      "–"
    )
  ),
  matrix_base %>% transmute(
    gene,
    evidence_source = "PXD\n046999\nspatial",
    status = ifelse(
      tolower(as.character(detected_in_nine_patient_dvp)) == "true",
      "Detectability only",
      "Not detected"
    ),
    display = ifelse(
      tolower(as.character(detected_in_nine_patient_dvp)) == "true",
      "yes",
      "no"
    )
  ),
  matrix_base %>% transmute(
    gene,
    evidence_source = "Tissue\npanel\nrole",
    status = case_when(
      gene %in% c("OLFM4", "FABP1") ~ "Direct locked arm",
      gene %in% c("SOX9", "CTNNB1") ~ "Context marker",
      TRUE ~ "Reserve locked marker"
    ),
    display = case_when(
      gene %in% c("OLFM4", "FABP1") ~ "direct",
      gene %in% c("SOX9", "CTNNB1") ~ "context",
      TRUE ~ "reserve"
    )
  )
) %>%
  mutate(
    marker = factor(
      unname(marker_labels[gene]),
      levels = rev(unname(marker_labels[marker_order]))
    ),
    evidence_source = factor(
      evidence_source,
      levels = c(
        "PXD\n002137\ncontrast",
        "PXD\n000445\npaired",
        "PXD\n017269\nFFPE",
        "PXD\n046999\nspatial",
        "Tissue\npanel\nrole"
      )
    )
  )

matrix_colours <- c(
  "FDR-supported difference" = COL[["route"]],
  "Direction only / imprecise" = COL[["route_light"]],
  "Detectability only" = COL[["wnt_light"]],
  "Not detected" = COL[["neutral_light"]],
  "Not estimable" = COL[["neutral_pale"]],
  "Not evaluated" = COL[["neutral_pale"]],
  "Direct locked arm" = "#F1C98A",
  "Context marker" = "#D8CDE5",
  "Reserve locked marker" = "#C8DCE8"
)

pS6c <- ggplot(matrix_plot, aes(evidence_source, marker, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = display), size = 1.80, colour = COL[["ink"]]) +
  scale_fill_manual(values = matrix_colours, guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    tag = "c",
    title = "Evidence roles are not interchangeable",
    subtitle = "Detection-only resources lack a normal comparator",
    x = NULL,
    y = NULL
  ) +
  theme_jtm(base_size = 6.1) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x.top = element_text(size = 5.0, lineheight = 0.90),
    axis.text.y = element_text(size = 5.5),
    panel.background = element_rect(fill = COL[["neutral_pale"]], colour = NA)
  )

figS6 <- pS6a | (pS6b / pS6c) +
  plot_layout(widths = c(0.96, 1.14), heights = c(0.46, 0.54)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_tsv(pxd2137_plot, "figureS6a_pxd002137_candidate_effects.tsv")
write_tsv(pxd445_plot, "figureS6b_pxd000445_paired_candidates.tsv")
write_tsv(matrix_plot, "figureS6c_public_protein_evidence_roles.tsv")

export_manifest <- bind_rows(
  export_figure(figS5, "figureS5_paired_ffpe_and_expanded_transcriptomics", 170, 132),
  export_figure(figS6, "figureS6_public_protein_triangulation", 170, 126)
)
export_manifest$file_size_bytes <- file.info(export_manifest$file)$size
export_manifest$sha256 <- vapply(
  export_manifest$file,
  function(path) strsplit(system2("sha256sum", path, stdout = TRUE), "[[:space:]]+")[[1]][1],
  character(1)
)
write.table(
  export_manifest,
  file.path(OUT_DIR, "figure_export_manifest_new_supplementary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

source_manifest <- data.frame(
  figure = c("Supplementary Fig. S5", "Supplementary Fig. S6"),
  panels = c("a-c", "a-c"),
  independent_unit = c(
    "patient pair for GSE117606; route nested within four recruitment clusters for expanded sensitivity",
    "specimen/patient pair as specified; detectability resources are not differential comparisons"
  ),
  inference = c(
    "paired Wilcoxon, patient-pair bootstrap, HC3 adjusted model, patient-clustered OLS",
    "age/sex-adjusted HC3 OLS, paired Wilcoxon with BH correction, categorical detectability audit"
  ),
  stringsAsFactors = FALSE
)
write.table(
  source_manifest,
  file.path(OUT_DIR, "figure_source_manifest_new_supplementary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

stopifnot(
  all(export_manifest$file_size_bytes > 0),
  nrow(paired_wide) == 51,
  sum(paired_wide$delta > 0) == 49,
  nrow(expanded_plot) == 7,
  n_distinct(expanded_plot$cluster_label) == 4,
  any(expanded_plot$cohort == "GSE226739"),
  all(c("OLFM4", "FABP1", "SOX9", "CTNNB1", "ETHE1") %in% matrix_plot$gene)
)

message("Exported Supplementary Figs. S5-S6 to ", OUT_DIR)
