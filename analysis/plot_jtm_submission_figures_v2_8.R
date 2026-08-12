#!/usr/bin/env Rscript

# JTM v2.8 figure package.
#
# Figure 5 integrates signed empirical perturbations and unsigned GenKI
# virtual-deletion results. Supplementary Figure S8 retains the frozen design,
# model robustness and the negative within-panel coherence boundary.

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
  file.path(getwd(), "analysis", "plot_jtm_submission_figures_v2_8.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT_DIR <- file.path(ROOT, "figures", "jtm_submission_v2.8")
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
RESULT_DIR <- file.path(ROOT, "results", "virtual_knockout_validation_v2_9")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

# Reuse the verified v2.6 visual system and all unchanged analytical panels.
# The legacy renderer is isolated in an environment; its own exports remain the
# verified v2.6 package and are not used as v2.8 deliverables.
base <- new.env(parent = globalenv())
sys.source(file.path(ROOT, "analysis", "plot_jtm_submission_figures_v2_6.R"), envir = base)

theme_jtm <- base$theme_jtm
tag_theme <- base$tag_theme
clean_panel <- base$clean_panel
COL <- base$COL
JTM_FONT <- base$JTM_FONT

read_tsv_path <- function(path) {
  read.delim(
    path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE, quote = "", comment.char = "", fileEncoding = "UTF-8"
  )
}
read_vk <- function(filename) read_tsv_path(file.path(RESULT_DIR, filename))
write_tsv <- function(frame, filename) {
  write.table(
    frame, file.path(SOURCE_DIR, filename), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )
}
fmt_p <- function(x) {
  ifelse(
    is.na(x), "",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
           ifelse(x < 0.10, formatC(x, format = "f", digits = 3),
                  formatC(x, format = "f", digits = 2)))
  )
}

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

replace_plot_labels <- function(plot, replacements) {
  for (i in seq_along(plot$layers)) {
    old <- plot$layers[[i]]$aes_params$label
    if (!is.null(old)) {
      hit <- as.character(old) %in% names(replacements)
      if (any(hit)) old[hit] <- unname(replacements[as.character(old[hit])])
      plot$layers[[i]]$aes_params$label <- old
    }
    layer_data <- plot$layers[[i]]$data
    if (is.data.frame(layer_data) && "label" %in% names(layer_data)) {
      labels <- as.character(layer_data$label)
      hit <- labels %in% names(replacements)
      if (any(hit)) labels[hit] <- unname(replacements[labels[hit]])
      layer_data$label <- labels
      plot$layers[[i]]$data <- layer_data
    }
  }
  plot
}

# -----------------------------------------------------------------------------
# Figure 1: one complete study workflow
# -----------------------------------------------------------------------------

workflow <- base$workflow
workflow$title[workflow$title == "Mechanistic support"] <- "Perturbation support"
workflow$detail[workflow$title == "Perturbation support"] <-
  "APC · WNT · ASCL2 · TCF7L2\nempirical + virtual"
p1a <- base$p1a
p1a$data <- workflow
p1a <- replace_plot_labels(p1a, c(
  "Error-controlled discovery core  →  fixed 12-gene representative panel  →  multi-layer validation" =
    "Threshold-defined 287-gene core  →  frozen 12-gene signature  →  multi-layer validation"
))
p1b <- base$p1b
flow <- base$flow
levels(flow$stage)[levels(flow$stage) == "Kneedle panel"] <- "Kneedle signature"
p1b$data <- flow
p1c <- base$p1c + labs(x = "Balanced signature size")
chen_scores <- base$chen_scores
levels(chen_scores$panel)[levels(chen_scores$panel) == "12-gene panel"] <- "12-gene signature"
p1e <- base$p1e
p1e$data <- chen_scores
fig1 <- clean_panel(p1a) /
  (clean_panel(p1b) | clean_panel(p1c)) /
  (clean_panel(base$p1d) | clean_panel(p1e) | clean_panel(base$p1f)) +
  plot_layout(heights = c(1.12, 1, 1), widths = c(1.05, 1, 0.82)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Figure 5: empirical direction plus unsigned virtual-deletion coupling
# -----------------------------------------------------------------------------

pathway_replacements <- c(
  "WNT withdrawal · programme ↓" = "WNT/RSPO withdrawal\nscore ↓",
  "APC loss ↑ · restoration ↓" = "APC loss ↑ · restoration ↓",
  "TF knockout / virtual TCF7L2 deletion · programme ↓" =
    "Experimental TF perturbation\nscore ↓",
  "ASCL2 target programme" = "ASCL2-linked transcription",
  "Stem / progenitor ↑" = "Stem/progenitor state ↑",
  "Differentiation ↓" = "Mature epithelial\nfunction ↓",
  "FIXED EPITHELIAL PROGRAMME" = "FROZEN 12-GENE SIGNATURE",
  "OBJECTIVE 12-GENE PANEL" = "FROZEN 12-GENE SIGNATURE",
  "50-gene up arm" = "Up arm · 6 genes",
  "50-gene down arm" = "Down arm · 6 genes"
)
p5a <- replace_plot_labels(base$p5a, pathway_replacements)

p5b <- base$p5b
p5c <- base$p5c

calibration <- read_vk("empirical_direction_calibration.tsv") %>%
  mutate(
    direction_supported = tolower(as.character(direction_supported)) %in%
      c("true", "t", "1", "yes"),
    display = case_when(
      comparison == "TCF7L2 KO/Het versus matched WT" ~ "TCF7L2 KO/Het vs WT",
      comparison == "APC_vs_WT_without_Wnt" ~ "APC loss vs WT (WNT−)",
      comparison == "WT_withdrawal" ~ "WNT/RSPO withdrawal",
      comparison == "genotype_by_Wnt_interaction" ~ "APC × WNT interaction",
      comparison == "ascl2_ko_vs_resting_wt" ~ "ASCL2 knockout",
      comparison == "conditional_wnt_silencing" ~ "Conditional WNT silencing",
      comparison == "apc_restoration_shApc" ~ "APC restoration",
      comparison == "apc_restoration_shApc_Kras" ~ "APC restoration + KRAS",
      TRUE ~ comparison
    ),
    display = factor(
      display,
      levels = rev(c(
        "APC loss vs WT (WNT−)", "APC × WNT interaction",
        "WNT/RSPO withdrawal", "TCF7L2 KO/Het vs WT", "ASCL2 knockout",
        "Conditional WNT silencing", "APC restoration",
        "APC restoration + KRAS"
      ))
    ),
    direction_label = recode(
      expected_direction, increase = "Expected increase",
      decrease = "Expected decrease"
    ),
    unit_class = ifelse(n_units == 1, "Descriptive n=1", "Multi-unit"),
    n_label = paste0("n=", n_units),
    text_x = mean_route_score_change + ifelse(mean_route_score_change >= 0, 0.11, -0.11),
    text_hjust = ifelse(mean_route_score_change >= 0, 0, 1)
  )
stopifnot(nrow(calibration) == 8L, all(calibration$direction_supported))

p5d <- ggplot(calibration, aes(y = display)) +
  geom_vline(xintercept = 0, linewidth = 0.38, colour = COL[["ink"]]) +
  geom_segment(
    aes(x = 0, xend = mean_route_score_change, yend = display,
        colour = direction_label), linewidth = 0.70, lineend = "round"
  ) +
  geom_point(
    aes(x = mean_route_score_change, fill = direction_label, shape = unit_class),
    size = 2.0, stroke = 0.45, colour = "white"
  ) +
  geom_text(
    aes(x = text_x, label = n_label, hjust = text_hjust),
    family = JTM_FONT, size = 1.30, colour = COL[["ink"]]
  ) +
  scale_colour_manual(values = c(
    "Expected increase" = "#C86B2B", "Expected decrease" = "#247BA0"
  ), guide = "none") +
  scale_fill_manual(values = c(
    "Expected increase" = "#C86B2B", "Expected decrease" = "#247BA0"
  ), guide = "none") +
  scale_shape_manual(values = c("Multi-unit" = 21, "Descriptive n=1" = 23), guide = "none") +
  scale_x_continuous(limits = c(-3.55, 3.55), breaks = -3:3) +
  labs(x = "Mean change in frozen 12-gene score", y = NULL) +
  theme_jtm(base_size = 5.25) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 3.65))

aggregate <- read_vk("genki_aggregate_validation_endpoints.tsv") %>%
  mutate(
    endpoint_label = recode(
      endpoint,
      upstream_context_to_core = "Upstream → core",
      upstream_context_to_panel = "Upstream → signature",
      within_panel_virtual_knockout_coherence = "Member → signature"
    ),
    endpoint_label = factor(
      endpoint_label,
      levels = rev(c("Upstream → core", "Upstream → signature", "Member → signature"))
    ),
    metric_label = factor(distance_metric, levels = c("KL", "EMD")),
    statistic = ifelse(
      distance_metric == "KL",
      paste0("q=", fmt_p(q_matched_primary_endpoints)),
      paste0("P=", fmt_p(p_matched_empirical_one_sided))
    ),
    supported = ifelse(
      distance_metric == "KL",
      !is.na(q_matched_primary_endpoints) & q_matched_primary_endpoints < 0.05,
      p_matched_empirical_one_sided < 0.05
    )
  )

p5e <- ggplot(aggregate, aes(x = matched_z, y = endpoint_label, colour = metric_label)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.35,
             colour = COL[["neutral"]]) +
  geom_segment(aes(x = 0, xend = matched_z, yend = endpoint_label),
               position = position_dodge(width = 0.45), linewidth = 0.52,
               colour = COL[["neutral_light"]]) +
  geom_point(aes(fill = supported, shape = metric_label),
             position = position_dodge(width = 0.45), size = 2.2, stroke = 0.65) +
  scale_colour_manual(values = c(KL = COL[["route"]], EMD = "#3C8D88"), name = NULL) +
  scale_fill_manual(values = c(`TRUE` = COL[["route"]], `FALSE` = "white"), guide = "none") +
  scale_shape_manual(values = c(KL = 21, EMD = 23), name = NULL) +
  coord_cartesian(xlim = c(-0.25, 5.35), clip = "off") +
  labs(x = "Enrichment vs matched sets (z)", y = NULL) +
  theme_jtm(base_size = 5.20) +
  theme(legend.position = "top", legend.key.width = unit(2.8, "mm"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 3.75))

upstream_tests <- read_vk("genki_fixed_gene_set_tests.tsv") %>%
  filter(
    distance_metric == "KL", target %in% c("TCF7L2", "ASCL2", "SOX4"),
    gene_set %in% c("measurable_fixed_287_core", "leave_target_out_fixed_12_panel")
  ) %>%
  mutate(
    target = factor(target, levels = rev(c("TCF7L2", "ASCL2", "SOX4"))),
    set_label = recode(
      gene_set, measurable_fixed_287_core = "287-gene core",
      leave_target_out_fixed_12_panel = "12-gene signature"
    ),
    set_label = factor(set_label, levels = c("287-gene core", "12-gene signature")),
    q_value = q_matched_primary_upstream_family,
    supported = !is.na(q_value) & q_value < 0.05,
    q_label = paste0("q=", fmt_p(q_value)),
    q_x = ifelse(matched_z > 3, matched_z - 0.18, matched_z + 0.18),
    q_hjust = ifelse(matched_z > 3, 1, 0)
  )
stopifnot(nrow(upstream_tests) == 6L)

p5f <- ggplot(upstream_tests, aes(x = matched_z, y = target)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.35,
             colour = COL[["neutral"]]) +
  geom_segment(aes(x = 0, xend = matched_z, yend = target),
               linewidth = 0.52, colour = COL[["neutral_light"]]) +
  geom_point(aes(fill = supported, colour = set_label), shape = 21,
             size = 2.0, stroke = 0.65) +
  geom_text(aes(x = q_x, label = q_label, hjust = q_hjust),
            position = position_nudge(y = 0.20), size = 1.23,
            family = JTM_FONT, colour = COL[["ink"]]) +
  facet_wrap(~set_label, nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = COL[["route"]], `FALSE` = "white"), guide = "none") +
  scale_colour_manual(values = c(
    "287-gene core" = COL[["route"]], "12-gene signature" = "#D4862A"
  ), guide = "none") +
  scale_x_continuous(limits = c(-1.9, 6.1), breaks = c(-1, 0, 2, 4, 6)) +
  labs(x = "Target-level matched-null z", y = NULL) +
  theme_jtm(base_size = 5.05) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 3.75), strip.text = element_text(size = 4.0),
        panel.spacing.x = unit(2.5, "mm"))

fig5_design <- "
AAAAAAAAAAAA
BBBBBCCCCCCC
DDDDDEEEFFFF
"
fig5 <- p5a + clean_panel(p5b) + clean_panel(p5c) + clean_panel(p5d) +
  clean_panel(p5e) + clean_panel(p5f) +
  plot_layout(design = fig5_design, heights = c(0.92, 0.95, 1.08)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Supplementary Figure S8: design, robustness and internal-panel boundary
# -----------------------------------------------------------------------------

workflow_s8 <- data.frame(
  xmin = c(0.15, 2.75, 5.35, 7.95), xmax = c(2.25, 4.85, 7.45, 10.05),
  ymin = 0.25, ymax = 1.70,
  x = c(1.20, 3.80, 6.40, 9.00), y = 0.98,
  label = c(
    "Frozen inputs\n287-gene core\n12-gene signature\n13 KO targets",
    "Held-out input\n13 adenoma donors\n128 cells per donor\n1,664 cells",
    "Dual-seed GenKI\nunsigned impact\nKL primary\nEMD sensitivity",
    "Validation tests\n10,000 matched sets\nno gene reselection\nno reweighting"
  ),
  fill = c("#E8F1F5", "#EEF2F4", "#E9F3F1", "#F8EFE3")
)
arrows_s8 <- data.frame(
  x = c(2.30, 4.90, 7.50), xend = c(2.70, 5.30, 7.90),
  y = 0.98, yend = 0.98
)
pS8a <- ggplot() +
  geom_rect(data = workflow_s8,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            colour = "white", linewidth = 0.55, show.legend = FALSE) +
  scale_fill_identity() +
  geom_segment(data = arrows_s8, aes(x = x, xend = xend, y = y, yend = yend),
               colour = COL[["neutral"]], linewidth = 0.48,
               arrow = arrow(length = unit(1.8, "mm"), type = "closed")) +
  geom_text(data = workflow_s8, aes(x = x, y = y, label = label),
            family = JTM_FONT, size = 2.15, lineheight = 0.95,
            colour = COL[["ink"]]) +
  annotate("text", x = 5.1, y = 2.04,
           label = "Frozen-gene validation only: no replacement, reweighting or target expansion",
           family = JTM_FONT, fontface = "bold", size = 2.25,
           colour = COL[["ink"]]) +
  coord_cartesian(xlim = c(0, 10.2), ylim = c(0.05, 2.25), clip = "off") +
  theme_void(base_family = JTM_FONT, base_size = 6.5) +
  theme(plot.margin = margin(3, 4, 3, 4, "mm"))

seed <- read_vk("genki_seed_stability.tsv") %>%
  filter(distance_metric == "KL") %>%
  transmute(target, measure = "Across seeds", rho = spearman_rho_between_seeds)
metric <- read_vk("genki_distance_metric_sensitivity.tsv") %>%
  transmute(target, measure = "KL vs EMD", rho = spearman_rho_kl_vs_emd)
robustness <- bind_rows(seed, metric) %>%
  mutate(measure = factor(measure, levels = c("Across seeds", "KL vs EMD")))
robustness_summary <- robustness %>%
  group_by(measure) %>%
  summarise(median_rho = median(rho), min_rho = min(rho),
            min_target = target[which.min(rho)], .groups = "drop")

pS8b <- ggplot(robustness, aes(x = measure, y = rho, colour = measure)) +
  geom_hline(yintercept = c(0.5, 0.75), linewidth = 0.3,
             linetype = c("dotted", "22"), colour = COL[["neutral_light"]]) +
  geom_boxplot(width = 0.38, outlier.shape = NA, linewidth = 0.42, fill = "white") +
  geom_point(size = 1.45, alpha = 0.82,
             position = position_jitter(width = 0.09, height = 0, seed = 20260810)) +
  geom_point(data = robustness_summary, aes(y = median_rho), shape = 95, size = 6.2) +
  geom_text(data = filter(robustness_summary, min_rho < 0.75),
            aes(y = min_rho, label = paste0(min_target, " ", sprintf("%.2f", min_rho))),
            nudge_y = -0.035, family = JTM_FONT, size = 1.6, show.legend = FALSE) +
  scale_colour_manual(values = c("Across seeds" = COL[["route"]], "KL vs EMD" = "#3C8D88"),
                      guide = "none") +
  scale_y_continuous(limits = c(0.43, 1.01), breaks = c(0.5, 0.75, 1.0)) +
  labs(x = NULL, y = "Spearman ρ") + theme_jtm(base_size = 5.8) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

panel_definition <- read_tsv_path(file.path(
  ROOT, "results", "objective_compact_panel_v2_7", "objective_compact_panel_frozen.tsv"
)) %>% select(gene, arm)
internal_panel <- read_vk("genki_fixed_gene_set_tests.tsv") %>%
  filter(
    distance_metric == "KL",
    target_role %in% c("panel_member", "upstream_context_and_panel_member"),
    gene_set == "leave_target_out_fixed_12_panel"
  ) %>%
  left_join(panel_definition, by = c("target" = "gene")) %>%
  mutate(
    q_value = q_matched_panel_connectivity_family,
    supported = !is.na(q_value) & q_value < 0.05,
    arm_label = ifelse(arm == "up", "Up arm", "Down arm"),
    target = factor(target, levels = target[order(matched_z)]),
    q_label = paste0("q=", fmt_p(q_value)),
    text_x = matched_z + ifelse(matched_z >= 0, 0.17, -0.17),
    text_hjust = ifelse(matched_z >= 0, 0, 1)
  )
stopifnot(nrow(internal_panel) == 12L, !any(internal_panel$supported))

pS8c <- ggplot(internal_panel, aes(x = matched_z, y = target, colour = arm_label)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.35,
             colour = COL[["neutral"]]) +
  geom_segment(aes(x = 0, xend = matched_z, yend = target), linewidth = 0.58) +
  geom_point(aes(fill = supported), shape = 21, size = 1.9, stroke = 0.60) +
  geom_text(aes(x = text_x, label = q_label, hjust = text_hjust),
            size = 1.25, family = JTM_FONT, colour = COL[["ink"]]) +
  scale_colour_manual(values = c("Up arm" = COL[["route"]], "Down arm" = "#D4862A"),
                      name = NULL) +
  scale_fill_manual(values = c(`TRUE` = COL[["route"]], `FALSE` = "white"), guide = "none") +
  coord_cartesian(xlim = c(min(internal_panel$matched_z) - 0.6,
                           max(internal_panel$matched_z) + 1.15), clip = "off") +
  labs(x = "Panel member → remaining signature (matched-null z)", y = NULL) +
  theme_jtm(base_size = 5.55) +
  theme(legend.position = "top", axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_text(size = 4.0))

figS8 <- pS8a / (clean_panel(pS8b) | clean_panel(pS8c)) +
  plot_layout(heights = c(0.62, 1.10), widths = c(0.72, 1.28)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Source data and exports
# -----------------------------------------------------------------------------

# Preserve the unchanged source tables in the new package, then overwrite the
# panels whose design or evidence role changed.
legacy_source <- file.path(ROOT, "figures", "jtm_submission_v2.6", "source_data")
legacy_files <- list.files(legacy_source, full.names = TRUE)
if (length(legacy_files)) file.copy(legacy_files, SOURCE_DIR, overwrite = TRUE)

write_tsv(workflow, "figure1a_workflow.tsv")
write_tsv(base$perturb_scores, "figure5b_apc_wnt_donor_scores.tsv")
write_tsv(base$apc, "figure5c_apc_wnt_contrasts.tsv")
write_tsv(calibration, "figure5d_empirical_direction_calibration.tsv")
write_tsv(aggregate, "figure5e_virtual_aggregate_endpoints.tsv")
write_tsv(upstream_tests, "figure5f_virtual_upstream_target_tests.tsv")
write_tsv(workflow_s8, "figureS8a_frozen_design.tsv")
write_tsv(robustness, "figureS8b_seed_metric_robustness.tsv")
write_tsv(internal_panel, "figureS8c_internal_panel_boundary.tsv")

exports <- bind_rows(
  export_figure(fig1, "figure1_discovery_core_and_objective_reduction", 170, 190),
  export_figure(base$fig2, "figure2_independent_replication_and_ffpe", 170, 148),
  export_figure(base$fig3, "figure3_rna_atac_regulatory_support", 170, 168),
  export_figure(base$fig4, "figure4_crc_atlas_cross_sectional_recurrence", 170, 186),
  export_figure(fig5, "figure5_empirical_and_virtual_perturbation_support", 170, 198),
  export_figure(base$fig6, "figure6_spatial_and_protein_context", 170, 182),
  export_figure(base$figS1, "figureS1_core_composition_and_portability", 170, 142),
  export_figure(base$figS2, "figureS2_external_and_ffpe_sensitivity", 170, 142),
  export_figure(base$figS6, "figureS3_signature_transparency_and_random_benchmark", 170, 160),
  export_figure(base$figS3, "figureS4_rna_atac_robustness", 170, 142),
  export_figure(base$figS4, "figureS5_crc_atlas_source_audit", 170, 150),
  export_figure(base$figS5, "figureS6_perturbation_boundaries", 170, 148),
  export_figure(figS8, "figureS7_virtual_knockout_robustness", 170, 118),
  export_figure(base$figS7, "figureS8_spatial_and_protein_assayability", 170, 145)
)
exports$file_size_bytes <- file.info(exports$file)$size
exports$sha256 <- vapply(
  exports$file,
  function(path) strsplit(system2("sha256sum", path, stdout = TRUE), "[[:space:]]+")[[1]][1],
  character(1)
)

svg_files <- exports$file[exports$format == "SVG"]
forbidden <- c(
  "TF knockout / virtual TCF7L2 deletion", "FIXED EPITHELIAL PROGRAMME",
  "OBJECTIVE 12-GENE PANEL", "Differentiation ↓", "Mechanistic support"
)
svg_text <- vapply(svg_files, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1))
exports$file <- sub(paste0(ROOT, "/"), "", exports$file, fixed = TRUE)
write_tsv(exports, "figure_export_manifest.tsv")
qa <- data.frame(
  check = c(
    "all_four_formats", "all_nonempty", "six_main_figures",
    "eight_supplementary_figures", "figure5_within_height_limit",
    "all_empirical_directions_supported", "negative_internal_coherence_retained",
    "forbidden_legacy_labels_absent"
  ),
  pass = c(
    all(table(exports$figure) == 4), all(exports$file_size_bytes > 0),
    length(unique(exports$figure[grepl("^figure[1-6]_", exports$figure)])) == 6,
    length(unique(exports$figure[grepl("^figureS", exports$figure)])) == 8,
    all(exports$height_mm[exports$figure == "figure5_empirical_and_virtual_perturbation_support"] <= 198),
    all(calibration$direction_supported), !any(internal_panel$supported),
    !any(vapply(forbidden, function(term) any(grepl(term, svg_text, fixed = TRUE)), logical(1)))
  )
)
write_tsv(qa, "figure_qa.tsv")
if (!all(qa$pass)) {
  stop("Figure v2.8 QA failed: ", paste(qa$check[!qa$pass], collapse = ", "))
}

session_lines <- capture.output(print(sessionInfo()))
session_lines <- gsub(ROOT, "<REPOSITORY_ROOT>", session_lines, fixed = TRUE)
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
if (nzchar(conda_prefix)) {
  session_lines <- gsub(conda_prefix, "<CONDA_ENV>", session_lines, fixed = TRUE)
}
writeLines(session_lines, file.path(SOURCE_DIR, "figure_sessionInfo.txt"), useBytes = TRUE)
message("Wrote JTM v2.8 R-only figure package: ", OUT_DIR)
