#!/usr/bin/env Rscript

# JTM v2.2 figure package. The 100-gene programme remains the reference; the
# fixed 10-gene panel is explicitly post hoc and exploratory. All drawing, assembly, export
# and preview generation are performed in R.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  file.path(getwd(), "analysis", "plot_jtm_submission_figures_v2_2.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT_DIR <- Sys.getenv(
  "JTM_FIGURE_DIR",
  unset = file.path(ROOT, "figures", "jtm_submission_v2.2")
)
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "dplyr", "ggplot2", "ggrepel", "patchwork", "scales", "tidyr",
  "svglite", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing locked-environment R packages: ", paste(missing_packages, collapse = ", "))
}

read_tsv <- function(path) {
  read.delim(
    file.path(ROOT, path), check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

write_tsv <- function(frame, filename) {
  write.table(
    frame, file.path(SOURCE_DIR, filename), sep = "\t", row.names = FALSE,
    quote = FALSE, na = ""
  )
}

JTM_FONT <- "Arial"
COL <- c(
  ink = "#25292D",
  route = "#C85E3D",
  route_light = "#E9B6A5",
  wnt = "#356F9D",
  wnt_light = "#B8D0E1",
  neutral = "#727B84",
  neutral_light = "#D8DDE1",
  neutral_pale = "#F3F5F6",
  adenoma = "#D9A640",
  crc = "#3F9A8B",
  uncertain = "#9A8877",
  context = "#7C6AA3",
  pass = "#4D8C75"
)

theme_jtm <- function(base_size = 7.1) {
  theme_classic(base_size = base_size, base_family = JTM_FONT) +
    theme(
      text = element_text(family = JTM_FONT, colour = COL[["ink"]]),
      axis.line = element_line(linewidth = 0.32, colour = COL[["ink"]]),
      axis.ticks = element_line(linewidth = 0.28, colour = COL[["ink"]]),
      axis.ticks.length = grid::unit(0.9, "mm"),
      axis.title = element_text(size = base_size, colour = COL[["ink"]]),
      axis.text = element_text(size = base_size - 0.55, colour = COL[["ink"]]),
      panel.grid = element_blank(),
      legend.background = element_blank(), legend.key = element_blank(),
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
      text = element_text(family = JTM_FONT, colour = COL[["ink"]]),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.tag = element_text(size = 8.2, face = "bold", family = JTM_FONT)
    )
}

tag_theme <- theme(
  plot.tag.position = c(0, 1),
  plot.tag = element_text(size = 8.2, face = "bold", family = JTM_FONT,
                          colour = COL[["ink"]])
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
    figure = stem, format = names(paths), file = normalizePath(unname(paths)),
    width_mm = width_mm, height_mm = height_mm,
    resolution_dpi = c(NA, NA, 600, 300), stringsAsFactors = FALSE
  )
}

# Source the verified v1.8 panel builders into an isolated environment. Their
# legacy exports are confined to the R temporary directory.
old_figure_dir <- Sys.getenv("JTM_FIGURE_DIR", unset = NA_character_)
env_v18 <- new.env(parent = globalenv())
Sys.setenv(JTM_FIGURE_DIR = file.path(tempdir(), "jtm_v22_v18_components"))
sys.source(file.path(ROOT, "analysis", "plot_jtm_submission_figures_v1_8.R"), envir = env_v18)
if (is.na(old_figure_dir)) Sys.unsetenv("JTM_FIGURE_DIR") else Sys.setenv(JTM_FIGURE_DIR = old_figure_dir)

# -----------------------------------------------------------------------------
# Figure 1: one complete study workflow, programme definition and held-out test
# -----------------------------------------------------------------------------

study_branches <- data.frame(
  y = c(6.05, 4.78, 3.50, 2.22, 0.95),
  title = c(
    "Independent replication", "Perturbation responsiveness",
    "Regulatory concordance", "Compact panel",
    "Spatial and protein readouts"
  ),
  detail = c(
    "Held-out · 5 cohorts · 51 FFPE pairs",
    "APC · WNT · ASCL2 · TCF7L2 · virtual KO",
    "Matched RNA–ATAC · CRC Atlas source audit",
    "Fixed 10 genes · no swapping · 7/7 gates",
    "6 Visium sections · OLFM4 gain · CA2 loss"
  ),
  figure = c("Figs. 1–2", "Fig. 3", "Fig. 4 / S5", "Fig. 5", "Fig. 6"),
  icon = c("replicate", "virtual", "wnt", "lock", "ffpe"),
  accent = unname(c(
    COL[["route"]], COL[["context"]], COL[["wnt"]],
    COL[["adenoma"]], COL[["crc"]]
  )),
  fill = c("#FCF6F3", "#F7F5FA", "#F5F8FA", "#FCF9F1", "#F3F8F7"),
  stringsAsFactors = FALSE
)

p1a <- env_v18$study_overview_panel(study_branches) +
  annotation_custom(
    grid::roundrectGrob(
      r = grid::unit(1.8, "mm"),
      gp = grid::gpar(
        col = scales::alpha(COL[["context"]], 0.48), fill = "#F7F5FA", lwd = 1.0
      )
    ),
    xmin = 0.15, xmax = 2.25, ymin = 2.12, ymax = 4.88
  ) +
  annotation_custom(
    grid::roundrectGrob(
      r = grid::unit(1.4, "mm"),
      gp = grid::gpar(
        col = scales::alpha(COL[["context"]], 0.34), fill = "white", lwd = 0.8
      )
    ),
    xmin = 0.70, xmax = 1.70, ymin = 3.62, ymax = 4.55
  ) +
  annotation_custom(
    env_v18$workflow_icon_grob("drug", COL[["context"]]),
    xmin = 0.84, xmax = 1.56, ymin = 3.78, ymax = 4.40
  ) +
  annotate("segment", x = 0.55, xend = 1.85, y = 4.78, yend = 4.78,
           colour = COL[["context"]], linewidth = 1.05, lineend = "round") +
  annotate("text", x = 1.20, y = 3.24, label = "CLINICAL UNMET NEED",
           size = 1.84, fontface = "bold", family = JTM_FONT, colour = COL[["ink"]]) +
  annotate("text", x = 1.20, y = 2.70,
           label = "Few reproducible tissue measures\nof early epithelial change",
           size = 1.48, lineheight = 0.94, family = JTM_FONT, colour = COL[["ink"]]) +
  annotate("label", x = 1.20, y = 2.24, label = "adenoma · archival tissue",
           size = 1.25, family = JTM_FONT, colour = COL[["context"]], fill = "white",
           linewidth = 0.20, label.padding = grid::unit(0.82, "mm")) +
  annotation_custom(
    grid::roundrectGrob(
      r = grid::unit(1.8, "mm"),
      gp = grid::gpar(
        col = scales::alpha(COL[["route"]], 0.48), fill = "#FCF6F3", lwd = 1.0
      )
    ),
    xmin = 2.78, xmax = 5.05, ymin = 2.12, ymax = 4.88
  ) +
  annotation_custom(
    grid::roundrectGrob(
      r = grid::unit(1.4, "mm"),
      gp = grid::gpar(
        col = scales::alpha(COL[["route"]], 0.34), fill = "white", lwd = 0.8
      )
    ),
    xmin = 3.42, xmax = 4.42, ymin = 3.62, ymax = 4.55
  ) +
  annotation_custom(
    env_v18$workflow_icon_grob("lock", COL[["route"]]),
    xmin = 3.56, xmax = 4.28, ymin = 3.78, ymax = 4.40
  ) +
  annotate("segment", x = 3.15, xend = 4.68, y = 4.78, yend = 4.78,
           colour = COL[["route"]], linewidth = 1.05, lineend = "round") +
  annotate("text", x = 3.92, y = 3.24, label = "DEFINE & LOCK STATE",
           size = 1.88, fontface = "bold", family = JTM_FONT, colour = COL[["ink"]]) +
  annotate("text", x = 3.92, y = 2.70,
           label = "Adenoma donor medians\n100 genes fixed before validation",
           size = 1.48, lineheight = 0.94, family = JTM_FONT, colour = COL[["ink"]]) +
  annotate("label", x = 3.92, y = 2.24, label = "no reselection · equal weights",
           size = 1.25, family = JTM_FONT, colour = COL[["route"]], fill = "white",
           linewidth = 0.20, label.padding = grid::unit(0.82, "mm")) +
  annotation_custom(
    grid::roundrectGrob(
      r = grid::unit(1.8, "mm"),
      gp = grid::gpar(col = scales::alpha(COL[["crc"]], 0.48), fill = "#F4F8F7", lwd = 1.0)
    ),
    xmin = 10.23, xmax = 12.70, ymin = 1.78, ymax = 5.22
  ) +
  annotate("segment", x = 10.62, xend = 12.31, y = 5.12, yend = 5.12,
           colour = COL[["crc"]], linewidth = 1.05, lineend = "round") +
  annotate("text", x = 11.46, y = 4.55, label = "TRANSLATIONAL OUTPUT",
           size = 1.72, fontface = "bold", family = JTM_FONT, colour = COL[["crc"]]) +
  annotate("text", x = 11.46, y = 3.75,
           label = "Compact candidate\ntissue panel",
           size = 1.72, fontface = "bold", lineheight = 0.94,
           family = JTM_FONT, colour = COL[["ink"]]) +
  annotate("segment", x = 10.66, xend = 10.94, y = 3.00, yend = 3.00,
           colour = COL[["route"]], linewidth = 1.05, lineend = "round") +
  annotate("text", x = 11.02, y = 3.00, hjust = 0, label = "100-gene reference programme",
           size = 1.30, fontface = "bold", family = JTM_FONT, colour = COL[["route"]]) +
  annotate("segment", x = 10.66, xend = 10.94, y = 2.55, yend = 2.55,
           colour = COL[["adenoma"]], linewidth = 1.05, lineend = "round") +
  annotate("text", x = 11.02, y = 2.55, hjust = 0, label = "10-gene compact representation",
           size = 1.30, fontface = "bold", family = JTM_FONT, colour = COL[["adenoma"]]) +
  annotate("text", x = 11.46, y = 2.08,
           label = "archival FFPE · prevention research",
           size = 1.13, family = JTM_FONT, colour = COL[["neutral"]])

fig1 <- clean_panel(p1a) / env_v18$fig1_mid / env_v18$fig1_bottom +
  plot_layout(heights = c(1.38, 1.06, 0.96)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
workflow_stages <- bind_rows(
  data.frame(
    stage_order = 1L, stage_type = "clinical_need", title = "Clinical unmet need",
    detail = "Few reproducible tissue measures of early epithelial change",
    figure = "Fig. 1a", stringsAsFactors = FALSE
  ),
  data.frame(
    stage_order = 2L, stage_type = "programme_lock", title = "Define and lock state",
    detail = "Adenoma donor medians; 100 genes fixed before validation",
    figure = "Fig. 1", stringsAsFactors = FALSE
  ),
  study_branches %>%
    transmute(
      stage_order = row_number() + 2L, stage_type = "evidence_layer",
      title = title, detail = detail, figure = figure
    ),
  data.frame(
    stage_order = 8L, stage_type = "translational_output", title = "Translational output",
    detail = "Compact candidate tissue panel for prospective analytical validation",
    figure = "Fig. 5 / Table 2", stringsAsFactors = FALSE
  )
)
write_tsv(workflow_stages, "figure1a_study_workflow.tsv")

# Existing quantitative pages are retained but reordered to follow the new
# evidence ladder: replication -> responsiveness -> regulatory concordance.
fig2 <- env_v18$fig2
fig3 <- env_v18$fig5
fig4 <- env_v18$fig3

# -----------------------------------------------------------------------------
# Figure 5: fixed exploratory 10-gene translational reduction
# -----------------------------------------------------------------------------

reduced <- read_tsv("results/translation_reduced_panel_v2_0/reduced_panel_definition.tsv") %>%
  mutate(
    arm = factor(
      panel_arm,
      levels = c("WNT_stem_progenitor_up", "mature_differentiation_down"),
      labels = c("WNT / stem–progenitor arm", "Mature-differentiation arm")
    ),
    gene = factor(gene, levels = gene[order(discovery_effect_adenoma_minus_normal)])
  )

p5a <- ggplot(reduced, aes(discovery_effect_adenoma_minus_normal, gene, colour = arm)) +
  geom_vline(xintercept = 0, colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = discovery_effect_adenoma_minus_normal, yend = gene),
               linewidth = 0.72, lineend = "round") +
  geom_point(size = 2.05) +
  scale_colour_manual(values = c(
    "WNT / stem–progenitor arm" = COL[["route"]],
    "Mature-differentiation arm" = COL[["wnt"]]
  )) +
  labs(x = "Discovery donor-median effect", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.8) +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 5.2), axis.text.y = element_text(size = 5.7))

concordance <- read_tsv(
  "results/translation_reduced_panel_v2_0/reduced_vs_100_gene_concordance.tsv"
) %>%
  filter(!(scope == "Chen" & cohort == "discovery")) %>%
  mutate(
    cohort_label = recode(cohort, validation = "Held-out Chen"),
    cohort_label = factor(cohort_label, levels = rev(cohort_label)),
    scope_label = recode(scope, Chen = "Held-out", external_sporadic = "External", FFPE = "FFPE")
  )
p5b <- ggplot(concordance, aes(spearman_rho_reduced_vs_100_gene, cohort_label,
                               colour = scope_label)) +
  geom_vline(xintercept = 0.80, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = 0.80, xend = spearman_rho_reduced_vs_100_gene, yend = cohort_label),
               linewidth = 0.62, lineend = "round") +
  geom_point(size = 1.95) +
  geom_text(aes(x = 0.974, label = sprintf("%.2f", spearman_rho_reduced_vs_100_gene)),
            hjust = 1, size = 1.60, family = JTM_FONT, colour = COL[["ink"]]) +
  scale_colour_manual(values = c("Held-out" = COL[["adenoma"]], External = COL[["route"]],
                                 FFPE = COL[["wnt"]])) +
  coord_cartesian(xlim = c(0.78, 0.982)) +
  labs(x = "Spearman correlation with 100-gene score", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.6) +
  theme(legend.position = "top", legend.justification = "left",
        axis.text.y = element_text(size = 5.4))

chen_scores <- read_tsv(
  "results/translation_reduced_panel_v2_0/chen_reduced_sample_scores.tsv"
) %>%
  filter(dataset == "validation", route_group %in% c("normal", "conventional_adenoma")) %>%
  mutate(tissue = factor(route_group,
                         levels = c("normal", "conventional_adenoma"),
                         labels = c("Normal", "Adenoma")))
chen_perf <- read_tsv(
  "results/translation_reduced_panel_v2_0/chen_reduced_performance.tsv"
) %>% filter(dataset == "validation")
p5c <- ggplot(chen_scores, aes(tissue, reduced_score, colour = tissue)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, colour = COL[["ink"]],
               fill = "#FAFBFB", linewidth = 0.40) +
  geom_point(position = position_jitter(width = 0.10, seed = 20260710),
             alpha = 0.76, size = 0.92) +
  annotate("text", x = 1.5, y = max(chen_scores$reduced_score) + 0.55,
           label = sprintf("AUC %.3f\nP = 4.25×10⁻⁶", chen_perf$auc),
           size = 1.76, lineheight = 0.92, family = JTM_FONT) +
  scale_colour_manual(values = c(Normal = COL[["neutral"]], Adenoma = COL[["route"]])) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "10-gene programme score") + theme_jtm() +
  theme(legend.position = "none")

external_tests <- read_tsv(
  "results/translation_reduced_panel_v2_0/external_reduced_cohort_tests.tsv"
) %>%
  transmute(
    label = cohort,
    estimate = clustered_standardized_mean_difference,
    ci_low = clustered_standardized_ci_low,
    ci_high = clustered_standardized_ci_high,
    type = "Cohort"
  )
pooled <- read_tsv(
  "results/translation_reduced_panel_v2_0/external_reduced_pooled_model.tsv"
) %>%
  transmute(label = "Pooled", estimate = adenoma_coef_sd,
            ci_low = ci_low, ci_high = ci_high, type = "Pooled")
external_display <- bind_rows(external_tests, pooled) %>%
  mutate(label = factor(label, levels = rev(c(external_tests$label, "Pooled"))))
p5d <- ggplot(external_display, aes(estimate, label, colour = type)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label), linewidth = 0.68,
               lineend = "round") +
  geom_point(aes(shape = type), size = 1.95) +
  scale_colour_manual(values = c(Cohort = COL[["route"]], Pooled = COL[["ink"]])) +
  scale_shape_manual(values = c(Cohort = 16, Pooled = 18)) +
  coord_cartesian(xlim = c(-0.85, 2.85)) +
  labs(x = "Standardised adenoma effect", y = NULL, colour = NULL, shape = NULL) +
  theme_jtm() + theme(legend.position = "none")

ffpe_scores <- read_tsv(
  "results/translation_reduced_panel_v2_0/ffpe_reduced_sample_scores.tsv"
)
complete_pairs <- ffpe_scores %>%
  filter(tissue_group %in% c("normal", "adenoma")) %>%
  distinct(patient_id, tissue_group) %>% count(patient_id) %>% filter(n == 2) %>% pull(patient_id)
ffpe_pair_long <- ffpe_scores %>%
  filter(patient_id %in% complete_pairs, tissue_group %in% c("normal", "adenoma")) %>%
  group_by(patient_id, tissue_group) %>%
  summarise(route_score_k5 = median(route_score_k5), .groups = "drop") %>%
  mutate(tissue = factor(tissue_group, levels = c("normal", "adenoma"),
                         labels = c("Reference", "Adenoma")))
ffpe_test <- read_tsv(
  "results/translation_reduced_panel_v2_0/ffpe_reduced_paired_test.tsv"
)
p5e <- ggplot(ffpe_pair_long, aes(tissue, route_score_k5, group = patient_id)) +
  geom_line(colour = COL[["neutral"]], linewidth = 0.30, alpha = 0.22) +
  geom_point(aes(colour = tissue), size = 0.78, alpha = 0.76) +
  stat_summary(aes(group = 1), fun = median, geom = "point", shape = 23, size = 2.35,
               stroke = 0.50, fill = "white", colour = COL[["ink"]]) +
  annotate("text", x = 1.5, y = max(ffpe_pair_long$route_score_k5) + 0.50,
           label = sprintf("49/51 increased\nmedian Δ %.2f; P = 1.25×10⁻⁹",
                           ffpe_test$median_paired_difference),
           size = 1.62, lineheight = 0.92, family = JTM_FONT) +
  scale_colour_manual(values = c(Reference = COL[["wnt"]], Adenoma = COL[["route"]])) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.15))) +
  labs(x = NULL, y = "10-gene FFPE score") + theme_jtm() +
  theme(legend.position = "none")

loo <- read_tsv(
  "results/translation_reduced_panel_v2_0/reduced_panel_leave_one_gene_out.tsv"
) %>%
  mutate(
    arm = factor(
      omitted_arm,
      levels = c("WNT_stem_progenitor_up", "mature_differentiation_down"),
      labels = c("Up arm", "Down arm")
    ),
    omitted_gene = factor(omitted_gene,
                          levels = omitted_gene[order(external_retained_fraction_vs_10_gene)])
  )
p5_loo <- ggplot(loo, aes(external_retained_fraction_vs_10_gene, omitted_gene, colour = arm)) +
  geom_vline(xintercept = 0.75, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_vline(xintercept = 1.00, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(aes(x = 0.75, xend = external_retained_fraction_vs_10_gene,
                   yend = omitted_gene), linewidth = 0.62, lineend = "round") +
  geom_point(size = 1.95) +
  scale_colour_manual(values = c("Up arm" = COL[["route"]], "Down arm" = COL[["wnt"]])) +
  coord_cartesian(xlim = c(0.74, 1.04)) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Pooled effect retained after gene omission", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.6) +
  theme(legend.position = "top", legend.justification = "left",
        axis.text.y = element_text(size = 5.5))

gates <- read_tsv(
  "results/translation_reduced_panel_v2_0/analysis_gate_summary.tsv"
) %>%
  mutate(
    observed_num = suppressWarnings(as.numeric(observed)),
    gate_label = factor(
      gate,
      levels = rev(gate),
      labels = rev(c("10/10 coverage", "Score concordance", "Held-out AUC",
                     "External effect", "FFPE pairs", "Perturbation direction",
                     "Leave-one-gene-out"))
    ),
    status = ifelse(pass, "Pass", "Fail"),
    observed_display = case_when(
      gate == "five_cohort_complete_10_of_10_coverage" ~ "5/5 cohorts; 10/10 each",
      gate == "score_concordance" ~ sprintf("min ρ = %.3f", observed_num),
      gate == "heldout_auc_retention" ~ sprintf("ΔAUC = %+.3f", observed_num),
      gate == "external_effect_retention" ~ sprintf("%.1f%% retained", 100 * observed_num),
      gate == "ffpe_pair_direction" ~ "49/51; P = 1.25×10⁻⁹",
      gate == "prespecified_perturbation_direction" ~ "APC 4/4; TCF7L2 KO 4/4",
      gate == "leave_one_gene_out" ~ sprintf("min %.1f%% retained", 100 * observed_num),
      TRUE ~ observed
    )
  )
p5f <- ggplot(gates, aes(1, gate_label, fill = status)) +
  geom_tile(width = 0.70, colour = "white", linewidth = 0.42) +
  geom_text(aes(label = ifelse(pass, "✓", "×")), size = 2.5,
            fontface = "bold", family = JTM_FONT, colour = "white") +
  geom_text(aes(x = 1.48, label = observed_display), hjust = 0, size = 1.45,
            family = JTM_FONT, colour = COL[["ink"]]) +
  scale_fill_manual(values = c(Pass = COL[["pass"]], Fail = "#B65B4B"), guide = "none") +
  coord_cartesian(xlim = c(0.55, 3.10), clip = "off") +
  labs(x = NULL, y = NULL) + theme_void(base_size = 6.2, base_family = JTM_FONT) +
  theme(axis.text.y = element_text(size = 5.2, colour = COL[["ink"]]),
        plot.margin = margin(1.4, 3.0, 1.4, 1.4, "mm"))

fig5 <- (clean_panel(p5a) | clean_panel(p5b)) /
  (clean_panel(p5c) | clean_panel(p5d)) /
  (clean_panel(p5e) | clean_panel(p5f)) +
  plot_layout(heights = c(1.05, 0.92, 1.02), widths = c(1, 1)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_tsv(reduced, "figure5a_reduced_panel_definition.tsv")
write_tsv(concordance, "figure5b_score_concordance.tsv")
write_tsv(chen_scores, "figure5c_heldout_scores.tsv")
write_tsv(external_display, "figure5d_external_effects.tsv")
write_tsv(ffpe_pair_long, "figure5e_ffpe_pairs.tsv")
write_tsv(gates, "figure5f_fidelity_gates.tsv")

# -----------------------------------------------------------------------------
# Figure 6: spatial localisation and protein anchors (OLFM4 / CA2 primary)
# -----------------------------------------------------------------------------

protein <- read_tsv(
  "results/translation_reduced_panel_v2_0/protein_anchor_evidence.tsv"
) %>%
  mutate(gene = factor(gene, levels = rev(c("OLFM4", "CA2", "FABP1"))))

p6e <- ggplot(protein, aes(age_sex_adjusted_log2_effect, gene,
                           colour = primary_tissue_anchor)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = age_sex_adjusted_ci_low, xend = age_sex_adjusted_ci_high,
                   yend = gene), linewidth = 0.70, lineend = "round") +
  geom_point(size = 2.05) +
  geom_text(aes(x = 5.85,
                label = sprintf("q = %.2g", age_sex_adjusted_q_bh_candidates)),
            hjust = 1,
            size = 1.55, family = JTM_FONT, colour = COL[["ink"]]) +
  scale_colour_manual(values = c(`TRUE` = COL[["route"]], `FALSE` = COL[["wnt"]])) +
  coord_cartesian(xlim = c(-4.65, 6.05), clip = "off") +
  labs(x = "PXD002137 adjusted log2 effect", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.5) + theme(legend.position = "none")

paired_protein <- protein %>%
  mutate(
    paired_status = case_when(
      is.na(n_pairs) ~ "Not detected",
      q_value_bh_within_analysis_set < 0.05 ~ "BH q < 0.05",
      TRUE ~ "Direction only"
    )
  )
p6f <- ggplot(paired_protein, aes(median_adenoma_minus_normal_log2, gene)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_point(data = filter(paired_protein, is.finite(median_adenoma_minus_normal_log2)),
             aes(fill = paired_status, size = paired_positive_fraction),
             shape = 21, colour = COL[["ink"]], stroke = 0.40) +
  geom_text(data = filter(paired_protein, !is.finite(median_adenoma_minus_normal_log2)),
            aes(x = 0, label = "not detected"), size = 1.55,
            family = JTM_FONT, colour = COL[["neutral"]]) +
  scale_fill_manual(values = c("BH q < 0.05" = COL[["route"]],
                               "Direction only" = COL[["neutral"]])) +
  scale_size_continuous(range = c(1.5, 2.9), limits = c(0.25, 1), guide = "none") +
  coord_cartesian(xlim = c(-1.25, 1.25)) +
  labs(x = "PXD000445 median paired log2 effect", y = NULL, fill = NULL) +
  theme_jtm(base_size = 6.5) + theme(legend.position = "none")

protein_matrix <- bind_rows(
  protein %>% transmute(
    gene, resource = "PXD002137\ndifferential",
    status = ifelse(age_sex_adjusted_q_bh_candidates < 0.05,
                    "FDR-supported difference", "Direction only / imprecise"),
    display = sprintf("%+.2f\nq %.2g", age_sex_adjusted_log2_effect,
                      age_sex_adjusted_q_bh_candidates)
  ),
  protein %>% transmute(
    gene, resource = "PXD000445\npaired",
    status = case_when(
      is.na(n_pairs) ~ "Not detected",
      q_value_bh_within_analysis_set < 0.05 ~ "FDR-supported difference",
      TRUE ~ "Direction only / imprecise"
    ),
    display = ifelse(is.na(n_pairs), "not\ndetected",
                     sprintf("%+.2f\nq %.2g", median_adenoma_minus_normal_log2,
                             q_value_bh_within_analysis_set))
  ),
  protein %>% transmute(
    gene, resource = "PXD017269\nFFPE detection",
    status = "Detectability only",
    display = percent(pxd017269_detection_fraction, accuracy = 1)
  ),
  protein %>% transmute(
    gene, resource = "PXD046999\nDVP detection",
    status = ifelse(detected_in_nine_patient_dvp, "Detectability only", "Not detected"),
    display = ifelse(detected_in_nine_patient_dvp, "detected", "not\ndetected")
  )
) %>%
  mutate(
    gene = factor(as.character(gene), levels = rev(c("OLFM4", "CA2", "FABP1"))),
    resource = factor(resource, levels = c(
      "PXD002137\ndifferential", "PXD000445\npaired",
      "PXD017269\nFFPE detection", "PXD046999\nDVP detection"
    ))
  )
protein_matrix_colours <- c(
  "FDR-supported difference" = "#E8B09E",
  "Direction only / imprecise" = "#F3D8CF",
  "Detectability only" = "#D7E5EE",
  "Not detected" = "#ECEFF1"
)
p6g <- ggplot(protein_matrix, aes(resource, gene, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = display), size = 1.75, lineheight = 0.88,
            family = JTM_FONT, colour = COL[["ink"]]) +
  scale_fill_manual(values = protein_matrix_colours, guide = "none") +
  scale_x_discrete(position = "top") + labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.3, base_family = JTM_FONT) +
  theme(panel.grid = element_blank(), axis.text.x.top = element_text(size = 5.1),
        axis.text.y = element_text(size = 5.7), axis.title = element_blank(),
        plot.margin = margin(1.2, 1.7, 1.2, 1.7, "mm"))

fig6_row2 <- (clean_panel(env_v18$p6d) | clean_panel(p6e) | clean_panel(p6f)) +
  plot_layout(widths = c(0.76, 1.12, 1.12))
fig6 <- (clean_panel(env_v18$p6a) | clean_panel(env_v18$p6b) | clean_panel(env_v18$p6c)) /
  fig6_row2 /
  clean_panel(p6g) +
  plot_layout(heights = c(1.03, 1.02, 0.84), widths = c(0.98, 1.01, 1.01)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(protein, "figure6e_f_protein_anchor_effects.tsv")
write_tsv(protein_matrix, "figure6g_protein_evidence_roles.tsv")

# -----------------------------------------------------------------------------
# Supplementary figures
# -----------------------------------------------------------------------------

figS1 <- env_v18$figS1
figS2 <- env_v18$figS2
pS3a <- ggplot(env_v18$feature_coverage,
               aes(feature, dataset, fill = coverage_fraction)) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = label), size = 1.65, family = JTM_FONT) +
  scale_fill_gradient(
    low = "#F3E4DF", high = COL[["route"]], limits = c(0, 1),
    breaks = c(0, 0.5, 1), labels = percent_format(accuracy = 1),
    name = "Coverage",
    guide = guide_colourbar(barwidth = grid::unit(28, "mm"),
                            barheight = grid::unit(2.8, "mm"),
                            title.position = "left")
  ) +
  labs(x = NULL, y = NULL) + theme_jtm() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 28, hjust = 1),
        legend.position = "top")
figS3 <- (clean_panel(pS3a) | clean_panel(env_v18$pS5b)) /
  (clean_panel(env_v18$pS5c) | clean_panel(env_v18$pS5d)) +
  plot_layout(heights = c(1.1, 0.9)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
figS4 <- env_v18$figS3
figS5 <- env_v18$fig4
figS7 <- clean_panel(env_v18$pS6a) /
  (clean_panel(env_v18$pS6b) | clean_panel(p6g)) +
  plot_layout(heights = c(0.72, 1.28), widths = c(1.08, 0.92)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(env_v18$feature_coverage, "figureS3a_perturbation_feature_coverage.tsv")
write_tsv(env_v18$evidence, "figureS3b_expression_matched_support.tsv")
write_tsv(env_v18$virtual_panel, "figureS3c_virtual_tf_deletion.tsv")
write_tsv(env_v18$ulm_coverage, "figureS3d_regulon_target_coverage.tsv")
write_tsv(env_v18$spatial_effects, "figureS7a_spatial_section_controls.tsv")
write_tsv(env_v18$pathology_summary, "figureS7b_spatial_pathology_summary.tsv")
write_tsv(protein_matrix, "figureS7c_protein_anchor_evidence_roles.tsv")

coverage <- read_tsv(
  "results/translation_reduced_panel_v2_0/reduced_panel_platform_coverage.tsv"
) %>%
  mutate(
    gene = factor(gene, levels = c("OLFM4", "ASCL2", "RNF43", "NKD1", "AXIN2",
                                   "FABP1", "CA2", "PCK1", "LGALS4", "AQP8")),
    cohort = factor(cohort, levels = rev(unique(cohort))),
    arm = ifelse(panel_arm == "WNT_stem_progenitor_up", "Up arm", "Down arm")
  )
pS6a <- ggplot(coverage, aes(gene, cohort, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = ifelse(present, "✓", "×")), size = 2.0, family = JTM_FONT) +
  scale_fill_manual(values = c(`TRUE` = "#DCECE5", `FALSE` = "#F0D7D1"), guide = "none") +
  labs(x = NULL, y = NULL) + theme_jtm(base_size = 6.3) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 32, hjust = 1, size = 5.2),
        axis.text.y = element_text(size = 5.2))

external_scores <- read_tsv(
  "results/translation_reduced_panel_v2_0/external_reduced_sample_scores.tsv"
)
external_corr <- concordance %>% filter(scope_label == "External") %>%
  transmute(cohort = as.character(cohort), label = sprintf("ρ = %.2f", spearman_rho_reduced_vs_100_gene))
pS6b <- ggplot(external_scores, aes(route_score_k50, route_score_k5, colour = tissue_group)) +
  geom_hline(yintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.25) +
  geom_point(size = 0.62, alpha = 0.60) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.45, colour = COL[["ink"]]) +
  geom_text(data = external_corr, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.10, vjust = 1.25, size = 1.45,
            family = JTM_FONT) +
  facet_wrap(~cohort, nrow = 1, scales = "free") +
  scale_colour_manual(values = c(normal = COL[["neutral"]], adenoma = COL[["route"]],
                                 crc = COL[["crc"]], hyperplastic = COL[["wnt"]]),
                      na.value = COL[["neutral_light"]]) +
  labs(x = "100-gene score", y = "10-gene score", colour = NULL) +
  theme_jtm(base_size = 5.8) +
  theme(legend.position = "top", strip.text = element_text(size = 5.1),
        axis.text = element_text(size = 4.5))

ffpe_corr <- concordance %>% filter(scope_label == "FFPE")
pS6c <- ggplot(ffpe_scores, aes(route_score_k50, route_score_k5, colour = tissue_group)) +
  geom_point(size = 0.72, alpha = 0.62) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.50, colour = COL[["ink"]]) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.20,
           label = sprintf("ρ = %.2f; n = %d",
                           ffpe_corr$spearman_rho_reduced_vs_100_gene, ffpe_corr$n),
           size = 1.65, family = JTM_FONT) +
  scale_colour_manual(values = c(normal = COL[["neutral"]], adenoma = COL[["route"]],
                                 crc = COL[["crc"]], ssa = COL[["wnt"]])) +
  labs(x = "100-gene FFPE score", y = "10-gene FFPE score", colour = NULL) +
  theme_jtm(base_size = 6.4) + theme(legend.position = "top")

apc_reduced <- read_tsv(
  "results/translation_reduced_panel_v2_0/reduced_apc_organoid_effects.tsv"
) %>% transmute(comparison, programme = "10-gene", mean_difference,
                ci_low = bootstrap_mean_ci_low, ci_high = bootstrap_mean_ci_high,
                expected_direction)
apc_full <- read_tsv(
  "results/perturbation_validation_locked_route/gse125472_contrast_summary.tsv"
) %>% filter(feature == "route_score") %>%
  transmute(comparison, programme = "100-gene", mean_difference,
            ci_low = bootstrap_mean_ci_low, ci_high = bootstrap_mean_ci_high,
            expected_direction)
apc_compare <- bind_rows(apc_full, apc_reduced) %>%
  mutate(
    comparison = factor(comparison, levels = c(
      "APC_vs_WT_with_Wnt", "APC_vs_WT_without_Wnt", "WT_withdrawal",
      "APC_withdrawal", "genotype_by_Wnt_interaction"
    ), labels = c("APC−WT\n+WNT", "APC−WT\n−WNT", "WT WNT\nwithdrawal",
                  "APC WNT\nwithdrawal", "Genotype ×\nWNT"))
  )
pS6d <- ggplot(apc_compare, aes(mean_difference, comparison, colour = programme)) +
  geom_vline(xintercept = 0, colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = comparison),
               position = position_dodge(width = 0.45), linewidth = 0.60) +
  geom_point(position = position_dodge(width = 0.45), size = 1.70) +
  scale_colour_manual(values = c("100-gene" = COL[["ink"]], "10-gene" = COL[["route"]])) +
  labs(x = "Programme-score difference", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.1) +
  theme(legend.position = "top", axis.text.y = element_text(size = 5.0))

tcf_reduced <- read_tsv(
  "results/translation_reduced_panel_v2_0/reduced_tcf7l2_clone_effects.tsv"
) %>% transmute(cell_line, clone_id, genotype, programme = "10-gene", difference_vs_WT)
tcf_full <- read_tsv(
  "results/perturbation_validation_locked_route/gse135328_clone_contrasts.tsv"
) %>% filter(feature == "route_score") %>%
  transmute(cell_line, clone_id, genotype, programme = "100-gene", difference_vs_WT)
tcf_compare <- bind_rows(tcf_full, tcf_reduced) %>%
  mutate(clone_label = factor(paste(cell_line, clone_id, genotype, sep = " · "),
                              levels = unique(paste(cell_line, clone_id, genotype, sep = " · "))))
pS6e <- ggplot(tcf_compare, aes(difference_vs_WT, clone_label, colour = programme)) +
  geom_vline(xintercept = 0, colour = COL[["neutral"]], linewidth = 0.35) +
  geom_point(position = position_dodge(width = 0.45), size = 1.75) +
  scale_colour_manual(values = c("100-gene" = COL[["ink"]], "10-gene" = COL[["wnt"]])) +
  labs(x = "TCF7L2-edited clone − WT", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.1) +
  theme(legend.position = "top", axis.text.y = element_text(size = 4.8))

pS6f <- p5_loo

figS6 <- clean_panel(pS6a) / clean_panel(pS6b) /
  (clean_panel(pS6c) | clean_panel(pS6d)) /
  (clean_panel(pS6e) | clean_panel(pS6f)) +
  plot_layout(heights = c(0.72, 0.88, 1.0, 1.0)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
write_tsv(coverage, "figureS6a_platform_coverage.tsv")
write_tsv(apc_compare, "figureS6d_apc_effect_comparison.tsv")
write_tsv(tcf_compare, "figureS6e_tcf7l2_effect_comparison.tsv")
write_tsv(loo, "figureS6f_leave_one_gene_out.tsv")

# -----------------------------------------------------------------------------
# Export, panel trace and QA
# -----------------------------------------------------------------------------

exports <- bind_rows(
  export_figure(fig1, "figure1_study_workflow_and_programme_definition", 170, 180),
  export_figure(fig2, "figure2_independent_replication_and_ffpe", 170, 178),
  export_figure(fig3, "figure3_perturbation_responsiveness", 170, 210),
  export_figure(fig4, "figure4_rna_atac_regulatory_support", 170, 142),
  export_figure(fig5, "figure5_reduced_translational_programme", 170, 190),
  export_figure(fig6, "figure6_spatial_and_protein_readouts", 170, 165),
  export_figure(figS1, "figureS1_programme_and_platform_coverage", 170, 112),
  export_figure(figS2, "figureS2_external_and_ffpe_sensitivity", 170, 132),
  export_figure(figS3, "figureS3_perturbation_boundaries", 170, 140),
  export_figure(figS4, "figureS4_rna_atac_robustness", 170, 124),
  export_figure(figS5, "figureS5_crc_atlas_recurrence_and_source_audit", 170, 160),
  export_figure(figS6, "figureS6_reduced_programme_sensitivity", 170, 170),
  export_figure(figS7, "figureS7_spatial_and_protein_assayability", 170, 145)
)
exports$file_size_bytes <- file.info(exports$file)$size
exports$sha256 <- vapply(
  exports$file,
  function(path) strsplit(system2("sha256sum", path, stdout = TRUE), "[[:space:]]+")[[1]][1],
  character(1)
)
write_tsv(exports, "figure_export_manifest.tsv")

main_counts <- c(5, 6, 6, 6, 6, 7)
panel_claims <- c(
  "Clinical unmet need, complete study workflow and translational output", "Discovery stability",
  "Biological anchors", "Held-out discrimination", "Donor-paired confirmation",
  "Five independent cohorts", "Pooled and proliferation-adjusted robustness",
  "Histological specificity", "Paired FFPE validation", "FFPE sensitivity estimands",
  "Expanded transcriptomic validation",
  "Canonical APC-WNT mechanism map", "APC-knockout donor effects",
  "Cross-perturbation effects", "Component responses", "Pharmacological boundary",
  "TCF7L2 regulatory consensus",
  "Becker tissue-state transfer", "Patient-clustered Becker model",
  "Matched RNA-ATAC correlation", "Adjusted RNA-from-ATAC model",
  "Patient fixed-effect sensitivity", "Regulatory-distance gradient",
  "Fixed 10-gene architecture", "Score concordance with 100-gene programme",
  "Held-out reduced-score performance", "Five-cohort reduced-score replication",
  "Paired FFPE reduced-score validation", "Seven descriptive fidelity gates",
  "Representative pathology map", "Raw spatial programme", "Adjusted spatial programme",
  "Six-section paired effects", "OLFM4 and CA2 differential protein anchors",
  "Paired proteomic direction", "Protein evidence roles"
)
stopifnot(length(panel_claims) == sum(main_counts))
panel_trace <- data.frame(
  figure = rep(paste("Figure", 1:6), main_counts),
  panel = unlist(lapply(main_counts, function(n) letters[seq_len(n)])),
  claim = panel_claims,
  primary_unit = c(
    rep("gene, specimen or evidence layer", 5),
    rep("patient or patient cluster", 6),
    rep("donor, model, clone or regulatory graph", 6),
    rep("sample with patient-aware inference", 6),
    rep("gene, specimen, patient pair or cohort", 6),
    rep("spot, section, specimen or patient pair", 7)
  ),
  stringsAsFactors = FALSE
)
write_tsv(panel_trace, "panel_source_trace.tsv")

svg_files <- exports$file[exports$format == "SVG"]
svg_has_text <- vapply(
  svg_files,
  function(path) any(grepl("<text", readLines(path, warn = FALSE, encoding = "UTF-8"), fixed = TRUE)),
  logical(1)
)
gate_table <- read_tsv("results/translation_reduced_panel_v2_0/analysis_gate_summary.tsv")
qa <- data.frame(
  check = c(
    "main_figure_count", "main_panel_count", "supplementary_figure_count",
    "formats_per_figure", "all_widths_170mm", "all_heights_within_225mm",
    "all_exports_nonempty", "pdf_tiff_under_10mb", "svg_text_editable",
    "reduced_panel_genes", "reduced_panel_gates", "primary_protein_anchors",
    "workflow_clinical_need"
  ),
  observed = c(
    length(unique(exports$figure[grepl("^figure[1-6]_", exports$figure)])),
    nrow(panel_trace),
    length(unique(exports$figure[grepl("^figureS", exports$figure)])),
    min(table(exports$figure)), all(exports$width_mm == 170),
    all(exports$height_mm <= 225), all(exports$file_size_bytes > 0),
    all(exports$file_size_bytes[exports$format %in% c("PDF", "TIFF")] < 10 * 1024^2),
    all(svg_has_text), nrow(reduced),
    sum(tolower(as.character(gate_table$pass)) == "true"),
    paste(sort(as.character(protein$gene[
      tolower(as.character(protein$primary_tissue_anchor)) == "true"
    ])), collapse = "+"),
    any(workflow_stages$stage_type == "clinical_need")
  ),
  expected = c(6, 36, 7, 4, TRUE, TRUE, TRUE, TRUE, TRUE, 10, 7, "CA2+OLFM4", TRUE),
  stringsAsFactors = FALSE
)
qa$pass <- mapply(function(observed, expected) {
  if (expected %in% c("TRUE", "FALSE")) return(as.character(observed) == expected)
  suppressWarnings({
    observed_num <- as.numeric(observed)
    expected_num <- as.numeric(expected)
  })
  if (is.finite(observed_num) && is.finite(expected_num)) {
    return(abs(observed_num - expected_num) < 1e-8)
  }
  as.character(observed) == as.character(expected)
}, qa$observed, qa$expected)
write_tsv(qa, "figure_build_qc.tsv")
if (!all(qa$pass)) stop("One or more v2.2 figure QA checks failed")

message("Wrote six main and seven supplementary JTM v2.2 figures to: ", OUT_DIR)
