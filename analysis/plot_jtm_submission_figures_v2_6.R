#!/usr/bin/env Rscript

# JTM v2.6 figures: threshold-defined 287-gene discovery core followed by a
# discovery-only, portability-gated, Kneedle-selected 12-gene operational panel.
# All analytical plots and exports are produced in R.

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
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  file.path(getwd(), "analysis", "plot_jtm_submission_figures_v2_6.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT_DIR <- file.path(ROOT, "figures", "jtm_submission_v2.6")
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

# Reuse only the verified visual system and the previously refined WNT/protein
# diagrams. The v2.5 analytical plots are rendered into a temporary directory.
old_figure_dir <- Sys.getenv("JTM_FIGURE_DIR", unset = NA_character_)
base <- new.env(parent = globalenv())
Sys.setenv(JTM_FIGURE_DIR = file.path(tempdir(), "jtm_v26_visual_components"))
sys.source(file.path(ROOT, "analysis", "plot_jtm_submission_figures_v2_5.R"), envir = base)
if (is.na(old_figure_dir)) Sys.unsetenv("JTM_FIGURE_DIR") else
  Sys.setenv(JTM_FIGURE_DIR = old_figure_dir)

theme_jtm <- base$theme_jtm
tag_theme <- base$tag_theme
clean_panel <- base$clean_panel
export_figure <- base$export_figure
COL <- base$COL
JTM_FONT <- base$JTM_FONT
base$env_v22$OUT_DIR <- OUT_DIR
base$env_v22$SOURCE_DIR <- SOURCE_DIR

read_tsv <- function(path) {
  read.delim(
    file.path(ROOT, path), sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, check.names = FALSE, quote = "",
    comment.char = ""
  )
}
write_tsv <- function(frame, filename) {
  write.table(
    frame, file.path(SOURCE_DIR, filename), sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )
}
as_bool <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
p_text <- function(x) {
  ifelse(is.na(x), "", ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
                                formatC(x, format = "f", digits = 3)))
}

PANEL_DIR <- "results/objective_compact_panel_v2_7"
EXT_DIR <- file.path(PANEL_DIR, "extended_validation")

# -----------------------------------------------------------------------------
# Figure 1: full workflow, threshold-defined core and objective reduction
# -----------------------------------------------------------------------------

workflow <- data.frame(
  x = 1:5,
  title = c("Discovery core", "Objective reduction", "Independent replication",
            "Mechanistic support", "Tissue context"),
  detail = c(
    "27 donors\ndonor-aware FDR + bootstrap",
    "287 → 62 → 12 genes\nno validation labels",
    "held-out · 5 cohorts\n51 paired FFPE specimens",
    "APC · WNT · ASCL2 · TCF7L2\nreciprocal perturbations",
    "RNA–ATAC · CRC Atlas\nspatial + protein anchors"
  ),
  accent = c(COL[["route"]], COL[["wnt"]], COL[["adenoma"]],
             COL[["context"]], COL[["crc"]]),
  fill = c("#FCF6F3", "#F5F8FA", "#FCF9F1", "#F7F5FA", "#F3F8F7")
)
p1a <- ggplot(workflow) +
  geom_segment(
    data = filter(workflow, x < 5),
    aes(x = x + 0.43, xend = x + 0.57, y = 1.38, yend = 1.38),
    colour = COL[["neutral"]], linewidth = 0.55,
    arrow = grid::arrow(type = "closed", angle = 24, length = grid::unit(1.25, "mm"))
  ) +
  geom_rect(aes(xmin = x - 0.43, xmax = x + 0.43, ymin = 0.46, ymax = 2.25,
                colour = accent, fill = fill), linewidth = 0.70) +
  geom_segment(aes(x = x - 0.30, xend = x + 0.30, y = 2.13, yend = 2.13,
                   colour = accent), linewidth = 1.1, lineend = "round") +
  geom_text(aes(x, 1.78, label = title, colour = accent), size = 1.82,
            fontface = "bold", family = JTM_FONT) +
  geom_text(aes(x, 1.00, label = detail), size = 1.33, lineheight = 0.94,
            family = JTM_FONT, colour = COL[["ink"]]) +
  annotate("text", x = 3, y = 2.72,
           label = "Can a reproducible early adenoma epithelial state be defined, reduced and measured across contexts?",
           size = 2.0, fontface = "bold", family = JTM_FONT, colour = COL[["ink"]]) +
  annotate("label", x = 3, y = 0.12,
           label = "Error-controlled discovery core  →  fixed 12-gene representative panel  →  multi-layer validation",
           size = 1.40, family = JTM_FONT, fill = "white", colour = COL[["ink"]],
           linewidth = 0.23, label.padding = grid::unit(0.75, "mm")) +
  scale_colour_identity() + scale_fill_identity() +
  coord_cartesian(xlim = c(0.45, 5.55), ylim = c(0.01, 2.90), clip = "off") +
  theme_void(base_family = JTM_FONT) +
  theme(plot.margin = margin(2, 2.5, 2, 2.5, "mm"))

flow <- data.frame(
  stage = factor(
    c("Assayed features", "Expressed / non-technical", "Directionally stable",
      "Bootstrap CI excludes zero", "Error-controlled core",
      "Portable protein-coding core", "Kneedle panel"),
    levels = rev(c("Assayed features", "Expressed / non-technical", "Directionally stable",
                   "Bootstrap CI excludes zero", "Error-controlled core",
                   "Portable protein-coding core", "Kneedle panel"))
  ),
  n = c(33698, 6127, 2496, 1504, 287, 62, 12),
  class = c(rep("Discovery filter", 4), "Threshold-defined core",
            "Portability gate", "Frozen panel")
)
p1b <- ggplot(flow, aes(n, stage, colour = class)) +
  geom_segment(aes(x = 8, xend = n, yend = stage), linewidth = 1.35, lineend = "round") +
  geom_point(size = 2.0) +
  geom_text(aes(label = comma(n)), hjust = -0.15, size = 1.55,
            family = JTM_FONT, colour = COL[["ink"]]) +
  scale_x_log10(limits = c(7, 70000), breaks = c(10, 100, 1000, 10000),
                labels = label_number(big.mark = ",")) +
  scale_colour_manual(values = c(
    "Discovery filter" = COL[["neutral"]], "Threshold-defined core" = COL[["route"]],
    "Portability gate" = COL[["wnt"]], "Frozen panel" = COL[["adenoma"]]
  ), guide = "none") +
  labs(x = "Genes retained (log scale)", y = NULL) + theme_jtm(base_size = 6.2) +
  theme(axis.text.y = element_text(size = 4.85), plot.margin = margin(1, 5, 1, 1, "mm"))

kneedle <- read_tsv(file.path(PANEL_DIR, "discovery_grouped_oof_fidelity_curve.tsv"))
p1c <- ggplot(kneedle, aes(total_genes, oof_spearman)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = alpha(COL[["route"]], 0.13)) +
  geom_line(colour = COL[["route"]], linewidth = 0.75) +
  geom_point(data = filter(kneedle, total_genes == 12), size = 2.35,
             colour = COL[["adenoma"]]) +
  geom_vline(xintercept = 12, linetype = "22", colour = COL[["adenoma"]], linewidth = 0.45) +
  annotate("label", x = 15.5, y = 0.80, label = "Kneedle: 12 genes\nOOF ρ = 0.929",
           size = 1.45, family = JTM_FONT, fill = "white", linewidth = 0.20) +
  coord_cartesian(xlim = c(2, 60), ylim = c(0.42, 1.01)) +
  labs(x = "Balanced panel size", y = "Leave-one-donor-out fidelity to 287-gene core") +
  theme_jtm(base_size = 6.1)

hallmark <- read_tsv("results/programme_transparency_v2_5/discovery_rank_based_hallmark_enrichment.tsv")
hallmark_display <- bind_rows(
  filter(hallmark, hallmark == "HALLMARK_WNT_BETA_CATENIN_SIGNALING"),
  hallmark %>% filter(enriched_direction == "adenoma_down") %>%
    arrange(rank_biserial) %>% slice_head(n = 4)
) %>%
  distinct(hallmark, .keep_all = TRUE) %>%
  mutate(label = gsub("_", " ", sub("HALLMARK_", "", hallmark)),
         label = tools::toTitleCase(tolower(label)),
         label = factor(label, levels = label[order(rank_biserial)]),
         direction = ifelse(rank_biserial > 0, "Adenoma-up", "Adenoma-down"))
p1d <- ggplot(hallmark_display, aes(rank_biserial, label, colour = direction)) +
  geom_vline(xintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = rank_biserial, yend = label), linewidth = 0.75) +
  geom_point(size = 2.0) +
  geom_text(aes(label = paste0("q=", formatC(q_bh, format = "g", digits = 2))),
            nudge_x = ifelse(hallmark_display$rank_biserial > 0, 0.05, -0.05),
            hjust = ifelse(hallmark_display$rank_biserial > 0, 0, 1),
            size = 1.25, family = JTM_FONT, colour = COL[["ink"]]) +
  scale_colour_manual(values = c("Adenoma-up" = COL[["route"]],
                                 "Adenoma-down" = COL[["wnt"]]), guide = "none") +
  coord_cartesian(xlim = c(-0.9, 0.78), clip = "off") +
  labs(x = "Rank-biserial enrichment across 6,127 genes", y = NULL) +
  theme_jtm(base_size = 6.0) + theme(axis.text.y = element_text(size = 4.55))

chen_scores <- read_tsv(file.path(PANEL_DIR, "validation_chen_panel_scores.tsv.gz")) %>%
  filter(dataset == "validation", panel_id %in% c("full_core_287", "objective_12"),
         route_group %in% c("normal", "conventional_adenoma")) %>%
  mutate(
    panel = factor(panel_id, levels = c("full_core_287", "objective_12"),
                   labels = c("287-gene core", "12-gene panel")),
    tissue = factor(route_group, levels = c("normal", "conventional_adenoma"),
                    labels = c("Normal", "Adenoma"))
  )
p1e <- ggplot(chen_scores, aes(tissue, panel_score, colour = tissue)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.43, alpha = 0.16) +
  geom_jitter(width = 0.12, size = 0.75, alpha = 0.68) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(Normal = COL[["wnt"]], Adenoma = COL[["route"]])) +
  labs(x = NULL, y = "Within-dataset score") + theme_jtm(base_size = 6.0) +
  theme(legend.position = "none", strip.text = element_text(size = 5.1),
        axis.text.x = element_text(size = 4.25))

chen_pair <- chen_scores %>% filter(panel == "12-gene panel") %>%
  group_by(donor_id, tissue) %>% summarise(score = median(panel_score), .groups = "drop") %>%
  add_count(donor_id) %>% filter(n == 2)
p1f <- ggplot(chen_pair, aes(tissue, score, group = donor_id)) +
  geom_line(colour = COL[["neutral"]], alpha = 0.40, linewidth = 0.42) +
  geom_point(aes(colour = tissue), size = 1.25) +
  stat_summary(aes(group = 1), fun = median, geom = "point", shape = 23,
               fill = "white", size = 2.35, stroke = 0.48, colour = COL[["ink"]]) +
  annotate("text", x = 1.5, y = max(chen_pair$score) + 0.35,
           label = "7/7 increased\npaired P = 0.0156", size = 1.45,
           family = JTM_FONT, lineheight = 0.90) +
  scale_colour_manual(values = c(Normal = COL[["wnt"]], Adenoma = COL[["route"]])) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(x = NULL, y = "12-gene score") + theme_jtm(base_size = 6.0) +
  theme(legend.position = "none", axis.text.x = element_text(size = 4.35))

fig1 <- clean_panel(p1a) /
  (clean_panel(p1b) | clean_panel(p1c)) /
  (clean_panel(p1d) | clean_panel(p1e) | clean_panel(p1f)) +
  plot_layout(heights = c(1.12, 1, 1), widths = c(1.05, 1, 0.82)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Figure 2: independent replication and archival-tissue measurability
# -----------------------------------------------------------------------------

external_tests <- read_tsv(file.path(PANEL_DIR, "validation_external_cohort_tests.tsv")) %>%
  filter(panel_id == "objective_12") %>%
  mutate(cohort = factor(cohort, levels = rev(cohort)))
p2a <- ggplot(external_tests, aes(clustered_standardized_mean_difference, cohort)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = clustered_standardized_ci_low,
                   xend = clustered_standardized_ci_high, yend = cohort),
               linewidth = 0.72, colour = COL[["route"]]) +
  geom_point(size = 2.0, colour = COL[["route"]]) +
  labs(x = "Adenoma effect (SD, patient-clustered 95% CI)", y = NULL) +
  theme_jtm(base_size = 6.2)

pooled <- read_tsv(file.path(PANEL_DIR, "validation_external_pooled_models.tsv")) %>%
  filter(panel_id == "objective_12") %>%
  transmute(model = "Cohort-adjusted", estimate = adenoma_coef_sd, ci_low, ci_high)
adjusted <- read_tsv(file.path(PANEL_DIR, "validation_external_proliferation_adjusted_models.tsv")) %>%
  filter(panel_id == "objective_12") %>%
  transmute(model = "Plus proliferation", estimate = adenoma_coef_sd, ci_low, ci_high)
model_summary <- bind_rows(pooled, adjusted) %>%
  mutate(model = factor(model, levels = rev(c("Cohort-adjusted", "Plus proliferation"))))
p2b <- ggplot(model_summary, aes(estimate, model, colour = model)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = model), linewidth = 0.82) +
  geom_point(size = 2.25) +
  geom_text(aes(label = sprintf("%.2f", estimate)), nudge_x = 0.12, hjust = 0,
            size = 1.45, family = JTM_FONT, colour = COL[["ink"]]) +
  scale_colour_manual(values = c("Cohort-adjusted" = COL[["route"]],
                                 "Plus proliferation" = COL[["wnt"]]), guide = "none") +
  coord_cartesian(xlim = c(0, 2.35)) +
  labs(x = "Five-cohort effect (95% CI)", y = NULL) + theme_jtm(base_size = 6.2)

p2c <- ggplot(external_tests, aes(auc_a_vs_b, cohort)) +
  geom_vline(xintercept = 0.5, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = 0.5, xend = auc_a_vs_b, yend = cohort), linewidth = 0.70,
               colour = COL[["neutral_light"]]) +
  geom_point(aes(colour = paired_primary), size = 2.0) +
  geom_text(aes(label = sprintf("%.2f", auc_a_vs_b)), nudge_x = 0.03,
            size = 1.42, family = JTM_FONT, colour = COL[["ink"]]) +
  scale_colour_manual(values = c(`TRUE` = COL[["adenoma"]], `FALSE` = COL[["route"]]),
                      guide = "none") +
  coord_cartesian(xlim = c(0.46, 1.08)) +
  labs(x = "AUC (adenoma vs normal)", y = NULL) + theme_jtm(base_size = 6.2)

ffpe_scores <- read_tsv(file.path(PANEL_DIR, "validation_ffpe_panel_scores.tsv.gz")) %>%
  filter(panel_id == "objective_12", tissue_group %in% c("normal", "adenoma")) %>%
  group_by(patient_id, tissue_group) %>%
  summarise(score = median(route_score_k5), .groups = "drop") %>%
  add_count(patient_id) %>% filter(n == 2) %>%
  mutate(tissue = factor(tissue_group, levels = c("normal", "adenoma"),
                         labels = c("Adjacent mucosa", "Adenoma")))
p2d <- ggplot(ffpe_scores, aes(tissue, score, group = patient_id)) +
  geom_line(colour = COL[["neutral"]], linewidth = 0.30, alpha = 0.22) +
  geom_point(aes(colour = tissue), size = 0.75, alpha = 0.75) +
  stat_summary(aes(group = 1), fun = median, geom = "point", shape = 23,
               size = 2.4, fill = "white", colour = COL[["ink"]]) +
  annotate("text", x = 1.5, y = max(ffpe_scores$score) + 0.48,
           label = "47/51 increased\nmedian Δ 1.62; P = 3.33×10⁻⁹",
           size = 1.42, family = JTM_FONT, lineheight = 0.90) +
  scale_colour_manual(values = c("Adjacent mucosa" = COL[["wnt"]],
                                 Adenoma = COL[["route"]])) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(x = NULL, y = "12-gene FFPE score") + theme_jtm(base_size = 6.1) +
  theme(legend.position = "none", axis.text.x = element_text(size = 5.2))

ffpe_genes <- read_tsv(file.path(PANEL_DIR, "validation_ffpe_gene_tests.tsv")) %>%
  filter(panel_id == "objective_12") %>%
  mutate(arm = ifelse(expected_direction_x > 0, "Up arm", "Down arm"),
         gene = factor(gene, levels = gene[order(median_paired_delta)]))
p2e <- ggplot(ffpe_genes, aes(median_paired_delta, gene, colour = arm)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = mean_delta_ci_low, xend = mean_delta_ci_high, yend = gene),
               linewidth = 0.62) +
  geom_point(size = 1.75) +
  scale_colour_manual(values = c("Up arm" = COL[["route"]], "Down arm" = COL[["wnt"]])) +
  labs(x = "Paired adenoma − mucosa expression", y = NULL, colour = NULL) +
  theme_jtm(base_size = 5.8) + theme(legend.position = "top", axis.text.y = element_text(size = 4.7))

becker <- read_tsv(file.path(EXT_DIR, "becker/becker_objective_panel_scores.tsv")) %>%
  filter(disease_stage_group %in% c("normal_unaffected", "polyp", "crc")) %>%
  mutate(state = factor(disease_stage_group,
                        levels = c("normal_unaffected", "polyp", "crc"),
                        labels = c("Normal / unaffected", "Polyp", "CRC")))
p2f <- ggplot(becker, aes(state, score_epi__ca_route_signature, colour = state)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.43, alpha = 0.15) +
  geom_jitter(width = 0.11, size = 0.72, alpha = 0.62) +
  scale_colour_manual(values = c("Normal / unaffected" = COL[["wnt"]],
                                 Polyp = COL[["adenoma"]], CRC = COL[["crc"]])) +
  labs(x = NULL, y = "Epithelial 12-gene score") + theme_jtm(base_size = 6.0) +
  theme(legend.position = "none", axis.text.x = element_text(size = 4.9))

fig2 <- (clean_panel(p2a) | clean_panel(p2b) | clean_panel(p2c)) /
  (clean_panel(p2d) | clean_panel(p2e) | clean_panel(p2f)) +
  plot_layout(heights = c(0.92, 1.08), widths = c(1, 0.96, 0.91)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Figure 3: matched RNA–ATAC regulatory concordance
# -----------------------------------------------------------------------------

rna_atac <- read_tsv(file.path(EXT_DIR, "becker_rna_atac/becker_rna_atac_paired_scores.tsv")) %>%
  filter(disease_stage_group %in% c("normal_unaffected", "polyp")) %>%
  mutate(state = factor(disease_stage_group,
                        levels = c("normal_unaffected", "polyp"),
                        labels = c("Normal / unaffected", "Polyp")))

rna_atac_flow <- data.frame(
  x = 1:4, label = c("Matched snRNA", "12-gene RNA score", "Promoter accessibility",
                     "Patient-aware association"),
  detail = c("40 samples · 12 patients", "epithelial nuclei", "WNT and TCF/ASCL2 axes",
             "proliferation + depth adjusted")
)
p3a <- ggplot(rna_atac_flow) +
  geom_segment(data = filter(rna_atac_flow, x < 4),
               aes(x = x + 0.36, xend = x + 0.64, y = 0.9, yend = 0.9),
               arrow = grid::arrow(type = "closed", length = grid::unit(1.2, "mm")),
               linewidth = 0.48, colour = COL[["neutral"]]) +
  geom_label(aes(x, 0.9, label = label), size = 1.55, fontface = "bold",
             family = JTM_FONT, fill = "white", colour = COL[["ink"]],
             label.padding = grid::unit(0.85, "mm"), linewidth = 0.35) +
  geom_text(aes(x, 0.30, label = detail), size = 1.25, family = JTM_FONT,
            colour = COL[["neutral"]]) +
  coord_cartesian(xlim = c(0.45, 4.55), ylim = c(0.02, 1.35), clip = "off") +
  theme_void(base_family = JTM_FONT)

p3b <- ggplot(rna_atac, aes(atac_tss__wnt_stem, rna_epi__ca_route_signature, colour = state)) +
  geom_smooth(method = "lm", se = TRUE, colour = COL[["neutral"]], fill = "#E7EBEE",
              linewidth = 0.50) +
  geom_point(size = 1.50, alpha = 0.80) +
  scale_colour_manual(values = c("Normal / unaffected" = COL[["wnt"]], Polyp = COL[["route"]])) +
  annotate("text", x = min(rna_atac$atac_tss__wnt_stem),
           y = max(rna_atac$rna_epi__ca_route_signature), hjust = 0,
           label = "Spearman ρ = 0.803\nP = 4.47×10⁻¹⁰", size = 1.45,
           family = JTM_FONT) +
  labs(x = "WNT/stem TSS accessibility", y = "Epithelial 12-gene RNA score", colour = NULL) +
  theme_jtm(base_size = 6.0) + theme(legend.position = "top")

p3c <- ggplot(rna_atac, aes(atac_tf__wnt_tcf_ascl2_axis, rna_epi__ca_route_signature, colour = state)) +
  geom_smooth(method = "lm", se = TRUE, colour = COL[["neutral"]], fill = "#E7EBEE",
              linewidth = 0.50) +
  geom_point(size = 1.50, alpha = 0.80) +
  scale_colour_manual(values = c("Normal / unaffected" = COL[["wnt"]], Polyp = COL[["route"]])) +
  annotate("text", x = min(rna_atac$atac_tf__wnt_tcf_ascl2_axis),
           y = max(rna_atac$rna_epi__ca_route_signature), hjust = 0,
           label = "Spearman ρ = 0.799\nP = 6.28×10⁻¹⁰", size = 1.45,
           family = JTM_FONT) +
  labs(x = "WNT/TCF/ASCL2 accessibility axis", y = "Epithelial 12-gene RNA score", colour = NULL) +
  theme_jtm(base_size = 6.0) + theme(legend.position = "none")

cluster_models <- read_tsv(file.path(EXT_DIR, "becker_rna_atac/becker_locked_rna_atac_patient_cluster_models.tsv")) %>%
  mutate(feature = recode(analysis_id,
    locked_route__wnt_tss = "WNT TSS",
    locked_route__wnt_tss_minus_housekeeping = "WNT TSS − housekeeping",
    locked_route__wnt_tcf_ascl2_axis = "TCF/ASCL2 axis",
    locked_route__wnt_tcf_ascl2_axis_minus_housekeeping = "TCF/ASCL2 − housekeeping"),
    feature = factor(feature, levels = rev(c("WNT TSS", "WNT TSS − housekeeping",
                                            "TCF/ASCL2 axis", "TCF/ASCL2 − housekeeping"))))
p3d <- ggplot(cluster_models, aes(coef, feature)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = feature), linewidth = 0.72,
               colour = COL[["route"]]) +
  geom_point(size = 1.95, colour = COL[["route"]]) +
  labs(x = "Patient-clustered adjusted coefficient (95% CI)", y = NULL) +
  theme_jtm(base_size = 6.0) + theme(axis.text.y = element_text(size = 4.7))

patient_median <- rna_atac %>% group_by(patient_id) %>%
  summarise(rna = median(rna_epi__ca_route_signature),
            atac = median(atac_tss__wnt_stem), .groups = "drop")
p3e <- ggplot(patient_median, aes(atac, rna)) +
  geom_smooth(method = "lm", se = TRUE, colour = COL[["neutral"]], fill = "#E7EBEE",
              linewidth = 0.48) +
  geom_point(size = 1.75, colour = COL[["route"]]) +
  annotate("text", x = min(patient_median$atac), y = max(patient_median$rna), hjust = 0,
           label = "Patient-median ρ = 0.790\npermutation P = 0.00292",
           size = 1.40, family = JTM_FONT) +
  labs(x = "Patient-median WNT TSS accessibility", y = "Patient-median RNA score") +
  theme_jtm(base_size = 6.0)

regulatory <- read_tsv(file.path(EXT_DIR, "becker_rna_atac/becker_objective_regulatory_window_rna_correlations.tsv")) %>%
  filter(rna_feature == "rna_epi__ca_route_signature") %>%
  mutate(distance = factor(distance_bin,
    levels = c("tss_core_1kb", "promoter_proximal_2_5kb", "proximal_regulatory_10kb", "distal_flank_20kb"),
    labels = c("TSS 1 kb", "2–5 kb", "10 kb", "20 kb")),
    locus = recode(locus_set, wnt_route_loci = "WNT-route loci",
                   wnt_tcf_ascl2_axis_loci = "TCF/ASCL2 loci"))
p3f <- ggplot(regulatory, aes(distance, spearman_rho, colour = locus, group = locus)) +
  geom_hline(yintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_line(linewidth = 0.70) + geom_point(size = 1.85) +
  scale_colour_manual(values = c("WNT-route loci" = COL[["route"]],
                                 "TCF/ASCL2 loci" = COL[["wnt"]])) +
  labs(x = "Regulatory distance window", y = "Spearman ρ with RNA score", colour = NULL) +
  theme_jtm(base_size = 6.0) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1, size = 4.8))

fig3 <- clean_panel(p3a) /
  (clean_panel(p3b) | clean_panel(p3c)) /
  (clean_panel(p3d) | clean_panel(p3e) | clean_panel(p3f)) +
  plot_layout(heights = c(0.58, 1.05, 0.95), widths = c(1.02, 0.94, 0.94)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Figure 4: CRC Atlas cross-sectional recurrence
# -----------------------------------------------------------------------------

atlas <- read_tsv(file.path(EXT_DIR, "crc_atlas/atlas_objective_panel_donor_scores.tsv")) %>%
  filter(n_cells_sampled >= 20) %>%
  mutate(state = factor(carrier_group,
    levels = c("normal_epithelial", "polyp_epithelial", "polyp_cancer",
               "tumor_epithelial", "tumor_cancer", "metastasis_epithelial", "metastasis_cancer"),
    labels = c("Normal\nepithelium", "Polyp\nepithelium", "Polyp\ncancer cells",
               "Primary\nepithelium", "Primary\ncancer cells", "Metastasis\nepithelium",
               "Metastasis\ncancer cells")))
p4a <- ggplot(atlas, aes(state, score__ca_route_signature, colour = state)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.40, alpha = 0.12) +
  geom_jitter(width = 0.12, size = 0.45, alpha = 0.30) +
  scale_colour_manual(values = c(COL[["wnt"]], COL[["adenoma"]], COL[["route"]],
                                 COL[["context"]], COL[["crc"]], "#507D73", "#204F46")) +
  labs(x = NULL, y = "Donor-carrier median 12-gene score") + theme_jtm(base_size = 5.7) +
  theme(legend.position = "none", axis.text.x = element_text(size = 4.35))

influence <- read_tsv(file.path(EXT_DIR, "crc_atlas/atlas_locked_study_influence.tsv")) %>%
  filter(outcome == "score__ca_route_signature") %>%
  mutate(state_label = factor(state,
    levels = rev(c("polyp_epithelial", "polyp_cancer", "tumor_epithelial", "tumor_cancer",
                   "metastasis_epithelial", "metastasis_cancer")),
    labels = rev(c("Polyp epithelium", "Polyp cancer cells", "Primary epithelium",
                   "Primary cancer cells", "Metastasis epithelium", "Metastasis cancer cells"))))
p4b <- ggplot(influence, aes(full_coef, state_label)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = full_ci_low, xend = full_ci_high, yend = state_label),
               linewidth = 0.75, colour = COL[["route"]]) +
  geom_point(size = 2.0, colour = COL[["route"]]) +
  labs(x = "Adjusted carrier-state coefficient (95% CI)", y = NULL) + theme_jtm(base_size = 6.0)

p4c <- ggplot(influence, aes(loo_min_coef, state_label)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = loo_min_coef, xend = loo_max_coef, yend = state_label),
               linewidth = 1.05, colour = COL[["neutral_light"]], lineend = "round") +
  geom_point(aes(x = full_coef), size = 1.9, colour = COL[["route"]]) +
  annotate("text", x = max(influence$loo_max_coef), y = 0.55,
           label = "33/33 omissions positive", hjust = 1, size = 1.35,
           family = JTM_FONT) +
  labs(x = "Coefficient range after one-study omission", y = NULL) + theme_jtm(base_size = 6.0)

support <- read_tsv(file.path(EXT_DIR, "crc_atlas/atlas_locked_state_study_support.tsv")) %>%
  filter(carrier_group != "normal_epithelial") %>%
  mutate(state = factor(carrier_group,
    levels = c("polyp_epithelial", "polyp_cancer", "tumor_epithelial", "tumor_cancer",
               "metastasis_epithelial", "metastasis_cancer"),
    labels = c("Polyp epi.", "Polyp cancer", "Primary epi.", "Primary cancer",
               "Metastasis epi.", "Metastasis cancer")))
study_order <- support %>% group_by(study_id) %>% summarise(total = sum(n_donors), .groups = "drop") %>%
  arrange(total) %>% pull(study_id)
support$study <- factor(support$study_id, levels = study_order)
p4d <- ggplot(support, aes(state, study, fill = n_donors)) +
  geom_tile(colour = "white", linewidth = 0.22) +
  scale_fill_gradient(low = "#F3F6F7", high = COL[["route"]], trans = "sqrt") +
  labs(x = NULL, y = "Source study", fill = "Donors") + theme_jtm(base_size = 5.5) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 4.1),
        axis.text.y = element_text(size = 3.6), axis.line = element_blank(),
        axis.ticks = element_blank(), legend.position = "top")

within <- read_tsv(file.path(EXT_DIR, "crc_atlas/atlas_locked_within_study_contrasts.tsv")) %>%
  filter(outcome == "score__ca_route_signature") %>%
  arrange(desc(n_donors)) %>% slice_head(n = 16) %>%
  mutate(label = paste0(gsub("_", " ", study_id), " · ", gsub("_", " ", state)),
         label = factor(label, levels = rev(label)))
p4e <- ggplot(within, aes(coef, label)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label), linewidth = 0.60,
               colour = COL[["context"]]) +
  geom_point(size = 1.55, colour = COL[["context"]]) +
  labs(x = "Eligible within-study contrast (95% CI)", y = NULL) + theme_jtm(base_size = 5.2) +
  theme(axis.text.y = element_text(size = 3.8))

fig4 <- clean_panel(p4a) /
  (clean_panel(p4b) | clean_panel(p4c)) /
  (clean_panel(p4d) | clean_panel(p4e)) +
  plot_layout(heights = c(0.86, 0.88, 1.26), widths = c(1.02, 0.98)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Figure 5: reciprocal perturbation responsiveness
# -----------------------------------------------------------------------------

p5a <- base$env_v22$env_v18$p5a
# The refined pathway drawing predates the objective panel. Relabel its readout
# in place so that the biology is retained without carrying forward the former
# fixed-size 50-up/50-down implementation.
pathway_label_updates <- c(
  "FIXED EPITHELIAL PROGRAMME" = "OBJECTIVE 12-GENE PANEL",
  "50-gene up arm" = "Up arm · 6 genes",
  "50-gene down arm" = "Down arm · 6 genes"
)
for (i in seq_along(p5a$layers)) {
  old_label <- p5a$layers[[i]]$aes_params$label
  if (!is.null(old_label) && length(old_label) == 1L && old_label %in% names(pathway_label_updates)) {
    p5a$layers[[i]]$aes_params$label <- unname(pathway_label_updates[[old_label]])
  }
}

perturb_scores <- read_tsv(file.path(PANEL_DIR, "validation_perturbation_sample_scores.tsv")) %>%
  filter(panel_id == "objective_12", dataset == "GSE125472") %>%
  mutate(condition = factor(paste(genotype, wnt_rspo),
    levels = c("WT with", "APC with", "WT without", "APC without"),
    labels = c("WT\n+WNT/RSPO", "APC-KO\n+WNT/RSPO", "WT\n−WNT/RSPO", "APC-KO\n−WNT/RSPO")))
p5b <- ggplot(perturb_scores, aes(condition, route_score, group = donor_id, colour = donor_id)) +
  geom_line(linewidth = 0.52, alpha = 0.72) + geom_point(size = 1.75) +
  scale_colour_manual(values = c("Donor1" = COL[["route"]], "Donor2" = COL[["wnt"]],
                                 "Donor3" = COL[["context"]])) +
  labs(x = NULL, y = "12-gene score", colour = NULL) + theme_jtm(base_size = 5.8) +
  theme(legend.position = "top", axis.text.x = element_text(size = 4.5))

apc <- read_tsv(file.path(PANEL_DIR, "validation_apc_organoid_effects.tsv")) %>%
  filter(panel_id == "objective_12", feature == "route_score") %>%
  mutate(label = recode(comparison,
    APC_vs_WT_with_Wnt = "APC-KO vs WT · +WNT",
    APC_vs_WT_without_Wnt = "APC-KO vs WT · −WNT",
    WT_withdrawal = "WNT withdrawal · WT",
    APC_withdrawal = "WNT withdrawal · APC-KO",
    genotype_by_Wnt_interaction = "Genotype × WNT interaction"),
    label = factor(label, levels = rev(label)),
    expected = ifelse(expected_direction == 0, "Context control", "Prespecified direction"))
p5c <- ggplot(apc, aes(mean_difference, label, colour = expected)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = bootstrap_mean_ci_low, xend = bootstrap_mean_ci_high, yend = label),
               linewidth = 0.68) + geom_point(size = 1.9) +
  scale_colour_manual(values = c("Prespecified direction" = COL[["route"]],
                                 "Context control" = COL[["neutral"]])) +
  labs(x = "Mean score contrast (bootstrap 95% CI)", y = NULL, colour = NULL) +
  theme_jtm(base_size = 5.8) + theme(legend.position = "top", axis.text.y = element_text(size = 4.4))

tcf <- read_tsv(file.path(PANEL_DIR, "validation_tcf7l2_clone_effects.tsv")) %>%
  filter(panel_id == "objective_12") %>%
  mutate(label = paste(cell_line, clone_id, genotype, sep = " · "),
         label = factor(label, levels = label[order(difference_vs_WT)]))
p5d <- ggplot(tcf, aes(difference_vs_WT, label, colour = genotype)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = difference_vs_WT, yend = label), linewidth = 0.68) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = c(KO = COL[["route"]], Het = COL[["adenoma"]])) +
  labs(x = "TCF7L2-edited clone − WT score", y = NULL, colour = NULL) +
  theme_jtm(base_size = 5.8) + theme(legend.position = "top", axis.text.y = element_text(size = 4.5))

extended_effect <- read_tsv(file.path(EXT_DIR, "perturbation_spatial/perturbation_effect_summary.tsv")) %>%
  filter(feature == "route_score", expected_direction != 0) %>%
  mutate(label = recode(comparison,
    trametinib_vs_dmso = "Trametinib vs DMSO",
    pri724_reversal_of_trametinib = "PRI-724 reversal",
    apc_restoration_shApc = "Apc restoration",
    apc_restoration_shApc_Kras = "Apc restoration + Kras",
    ascl2_ko_vs_resting_wt = "Ascl2 KO",
    conditional_wnt_silencing = "Conditional WNT silencing"),
    label = factor(label, levels = rev(label)),
    aligned = expected_direction * mean_difference > 0)
p5e <- ggplot(extended_effect, aes(mean_difference, label, colour = aligned)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = bootstrap_mean_ci_low, xend = bootstrap_mean_ci_high, yend = label),
               linewidth = 0.65) + geom_point(size = 1.85) +
  scale_colour_manual(values = c(`TRUE` = COL[["route"]], `FALSE` = COL[["crc"]]), guide = "none") +
  labs(x = "Intervention effect on 12-gene score", y = NULL) + theme_jtm(base_size = 5.8)

matched <- read_tsv(file.path(EXT_DIR, "perturbation_spatial/expression_matched_tests.tsv")) %>%
  filter(comparison != "doxycycline_control_shRenilla") %>%
  mutate(label = paste(dataset, gsub("_", " ", comparison), sep = " · "),
         label = factor(label, levels = rev(label)))
p5f <- ggplot(matched, aes(observed_mean_difference, label)) +
  geom_vline(xintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(aes(x = null_ci_low, xend = null_ci_high, yend = label),
               linewidth = 1.05, colour = COL[["neutral_light"]], lineend = "round") +
  geom_point(size = 1.9, colour = COL[["route"]]) +
  geom_text(aes(label = paste0("P=", p_text(p_expression_matched_one_sided))),
            nudge_x = 0.14, hjust = 0, size = 1.18, family = JTM_FONT) +
  coord_cartesian(xlim = c(min(matched$null_ci_low) - 0.1,
                           max(matched$null_ci_high) + 0.55), clip = "off") +
  labs(x = "Effect vs matched-null 95% range", y = NULL) +
  theme_jtm(base_size = 5.2) + theme(axis.text.y = element_text(size = 3.8))

fig5 <- clean_panel(p5a) /
  (clean_panel(p5b) | clean_panel(p5c)) /
  (clean_panel(p5d) | clean_panel(p5e) | clean_panel(p5f)) +
  plot_layout(heights = c(0.93, 0.95, 1.08), widths = c(0.90, 1.03, 1.07)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Figure 6: spatial localisation, protein anchors and translational boundary
# -----------------------------------------------------------------------------

spatial_tests <- read_tsv(file.path(EXT_DIR, "perturbation_spatial/spatial_objective_panel_tests.tsv"))
spatial_units <- read_tsv(file.path(EXT_DIR, "perturbation_spatial/spatial_objective_panel_section_effects.tsv"))
spatial_summary <- read_tsv(file.path(EXT_DIR, "perturbation_spatial/spatial_objective_panel_pathology_summary.tsv"))
spatial_spots <- read_tsv(file.path(EXT_DIR, "perturbation_spatial/spatial_objective_panel_spot_scores.tsv.gz"))
map_samples <- spatial_summary %>%
  filter(pathology_group %in% c("tumor", "tumor_stroma", "non_neoplastic_epithelium")) %>%
  distinct(sample_id) %>% slice_head(n = 2) %>% pull(sample_id)
spatial_map <- spatial_spots %>% filter(sample_id %in% map_samples, in_tissue == 1)
p6a <- ggplot(spatial_map, aes(pxl_col, -pxl_row, colour = route_score)) +
  geom_point(size = 0.52) + facet_wrap(~sample_id, nrow = 1) +
  scale_colour_gradient2(low = "#3B6C8E", mid = "#F4F4F2", high = "#B24D3E",
                         midpoint = 0, oob = squish) +
  coord_equal() + labs(x = NULL, y = NULL, colour = "12-gene\nscore") +
  theme_void(base_family = JTM_FONT) +
  theme(strip.text = element_text(size = 5.0), legend.position = "right")

pathology <- spatial_summary %>%
  filter(pathology_group %in% c("tumor", "tumor_stroma", "non_neoplastic_epithelium", "stroma")) %>%
  mutate(pathology = factor(pathology_group,
    levels = c("non_neoplastic_epithelium", "stroma", "tumor_stroma", "tumor"),
    labels = c("Non-neoplastic epithelium", "Stroma", "Tumour + stroma", "Tumour")))
p6b <- ggplot(pathology, aes(pathology, route_score, group = sample_id, colour = sample_id)) +
  geom_line(linewidth = 0.35, alpha = 0.45) + geom_point(size = 0.85) +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = COL[["ink"]],
               linewidth = 0.85) +
  labs(x = NULL, y = "Section-median 12-gene score") + theme_jtm(base_size = 5.8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 24, hjust = 1, size = 4.4))

spatial_primary <- spatial_units %>%
  filter(comparison == "tumor_vs_non_neoplastic_epithelium",
         feature %in% c("route_score", "route_residual_prolif_epithelial")) %>%
  mutate(feature_label = recode(feature, route_score = "Raw 12-gene score",
                                route_residual_prolif_epithelial = "Adjusted score"))
p6c <- ggplot(spatial_primary, aes(difference, sample_id, colour = feature_label)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = difference, yend = sample_id), linewidth = 0.60) +
  geom_point(size = 1.75) +
  scale_colour_manual(values = c("Raw 12-gene score" = COL[["route"]],
                                 "Adjusted score" = COL[["wnt"]])) +
  labs(x = "Tumour − non-neoplastic epithelium", y = "Section", colour = NULL) +
  theme_jtm(base_size = 5.8) + theme(legend.position = "top")

p6d <- ggplot(filter(spatial_tests, comparison == "tumor_vs_non_neoplastic_epithelium"),
              aes(median_difference, reorder(feature, median_difference))) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = median_difference, yend = reorder(feature, median_difference)),
               linewidth = 0.72, colour = COL[["neutral_light"]]) +
  geom_point(aes(colour = feature == "route_score"), size = 2.0) +
  geom_text(aes(label = paste0(n_positive, "/", n_sections)), nudge_x = 0.08,
            size = 1.35, family = JTM_FONT) +
  scale_colour_manual(values = c(`TRUE` = COL[["route"]], `FALSE` = COL[["neutral"]]),
                      guide = "none") +
  scale_y_discrete(labels = c(route_score = "12-gene score",
                              route_residual_prolif_epithelial = "Adjusted score",
                              wnt_stem = "WNT/stemness", proliferation_control = "Proliferation")) +
  labs(x = "Median paired section difference", y = NULL) + theme_jtm(base_size = 5.8)

p6e <- base$env_v22$p6e
p6f <- base$env_v22$p6f
p6g <- base$env_v22$p6g +
  theme(axis.text.x.top = element_text(size = 3.45, lineheight = 0.82))

fig6 <- (clean_panel(p6a) | clean_panel(p6b)) /
  (clean_panel(p6c) | clean_panel(p6d)) /
  (clean_panel(p6e) | clean_panel(p6f) | clean_panel(p6g)) +
  plot_layout(heights = c(1.0, 0.91, 0.94), widths = c(1, 0.92, 1.08)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Supplementary figures
# -----------------------------------------------------------------------------

core_evidence <- read_tsv("results/data_adaptive_panel_pilot_v2_6/discovery_gene_evidence.tsv") %>%
  filter(!excluded_gene, mean_expression_discovery > 0.001) %>%
  mutate(core = as_bool(stable_core), selected12 = gene %in%
           read_tsv(file.path(PANEL_DIR, "objective_compact_panel_frozen.tsv"))$gene)
pS1a <- ggplot(core_evidence, aes(logFC, -log10(adj.P.Val), colour = core)) +
  geom_point(size = 0.42, alpha = 0.35) +
  geom_point(data = filter(core_evidence, selected12), size = 1.45, colour = COL[["adenoma"]]) +
  geom_text_repel(data = filter(core_evidence, selected12), aes(label = gene),
                  size = 1.20, family = JTM_FONT, colour = COL[["ink"]],
                  max.overlaps = Inf, min.segment.length = 0, box.padding = 0.22) +
  scale_colour_manual(values = c(`FALSE` = COL[["neutral_light"]], `TRUE` = COL[["route"]]),
                      guide = "none") +
  labs(x = "Donor-aware adenoma effect", y = "−log10 limma FDR") + theme_jtm(base_size = 5.8)

panel_def <- read_tsv(file.path(PANEL_DIR, "objective_compact_panel_frozen.tsv")) %>%
  mutate(gene = factor(gene, levels = gene[order(logFC)]),
         arm_label = ifelse(arm == "up", "Up arm", "Down arm"))
pS1b <- ggplot(panel_def, aes(logFC, gene, colour = arm_label)) +
  geom_vline(xintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(aes(x = bootstrap_ci_low, xend = bootstrap_ci_high, yend = gene), linewidth = 0.60) +
  geom_point(size = 1.70) +
  scale_colour_manual(values = c("Up arm" = COL[["route"]], "Down arm" = COL[["wnt"]])) +
  labs(x = "Discovery effect (bootstrap 95% CI)", y = NULL, colour = NULL) +
  theme_jtm(base_size = 5.8) + theme(legend.position = "top")

presence <- read_tsv(file.path(PANEL_DIR, "gene_platform_presence_long.tsv")) %>%
  filter(gene %in% panel_def$gene) %>%
  mutate(gene = factor(gene, levels = panel_def$gene), present = as_bool(feature_present))
pS1c <- ggplot(presence, aes(platform, gene, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.42) +
  scale_fill_manual(values = c(`TRUE` = COL[["route"]], `FALSE` = "#ECEFF1"), guide = "none") +
  labs(x = NULL, y = NULL) + theme_jtm(base_size = 5.5) +
  theme(axis.text.x = element_text(angle = 32, hjust = 1, size = 4.1),
        axis.line = element_blank(), axis.ticks = element_blank())

arm_counts <- data.frame(arm = c("Core up", "Core down", "Portable up", "Portable down", "Panel up", "Panel down"),
                         n = c(89, 198, 30, 32, 6, 6),
                         stage = rep(c("287-gene core", "62-gene portable", "12-gene panel"), each = 2))
pS1d <- ggplot(arm_counts, aes(n, reorder(arm, n), fill = stage)) +
  geom_col(width = 0.62) + geom_text(aes(label = n), hjust = -0.2, size = 1.45, family = JTM_FONT) +
  scale_fill_manual(values = c("287-gene core" = COL[["neutral"]],
                               "62-gene portable" = COL[["wnt"]],
                               "12-gene panel" = COL[["route"]])) +
  coord_cartesian(xlim = c(0, 225), clip = "off") +
  labs(x = "Genes", y = NULL, fill = NULL) + theme_jtm(base_size = 5.8) +
  theme(legend.position = "top")
figS1 <- (clean_panel(pS1a) | clean_panel(pS1b)) /
  (clean_panel(pS1c) | clean_panel(pS1d)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

loo_cohort <- read_tsv(file.path(PANEL_DIR, "validation_external_adjusted_leave_one_cohort_out.tsv")) %>%
  filter(panel_id == "objective_12") %>% mutate(excluded_cohort = factor(excluded_cohort, levels = rev(excluded_cohort)))
pS2a <- ggplot(loo_cohort, aes(adenoma_coef_sd, excluded_cohort)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = excluded_cohort), linewidth = 0.70,
               colour = COL[["route"]]) + geom_point(size = 1.9, colour = COL[["route"]]) +
  labs(x = "Adjusted effect after cohort omission", y = "Excluded cohort") + theme_jtm(base_size = 6.0)
pS2b <- p2e
pS2c <- p2a
pS2d <- p2d
figS2 <- (clean_panel(pS2a) | clean_panel(pS2b)) /
  (clean_panel(pS2c) | clean_panel(pS2d)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

residual <- read_tsv(file.path(EXT_DIR, "becker_rna_atac/becker_rna_atac_residual_correlations.tsv")) %>%
  filter(subset == "normal_polyp", grepl("rna_epi_ca_route", analysis_id)) %>%
  mutate(feature = factor(gsub("rna_epi_ca_route__", "", analysis_id), levels = rev(gsub("rna_epi_ca_route__", "", analysis_id))))
pS3a <- p3e
pS3b <- p3d
pS3c <- ggplot(residual, aes(residual_spearman_rho, feature)) +
  geom_vline(xintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = residual_spearman_rho, yend = feature), linewidth = 0.70,
               colour = COL[["neutral_light"]]) + geom_point(size = 1.9, colour = COL[["route"]]) +
  labs(x = "Residual Spearman ρ", y = NULL) + theme_jtm(base_size = 5.6) +
  theme(axis.text.y = element_text(size = 4.2))
pS3d <- p3f
figS3 <- (clean_panel(pS3a) | clean_panel(pS3b)) /
  (clean_panel(pS3c) | clean_panel(pS3d)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

figS4 <- (clean_panel(p4b) | clean_panel(p4c)) /
  (clean_panel(p4d) | clean_panel(p4e)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

coverage <- read_tsv(file.path(EXT_DIR, "perturbation_spatial/feature_coverage.tsv")) %>%
  mutate(dataset = gsub("spatial_", "", dataset),
         feature = factor(feature, levels = c("route_up", "route_down", "wnt_stem",
                                              "proliferation_control", "epithelial_control")))
pS5a <- ggplot(coverage, aes(feature, dataset, fill = coverage_fraction)) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = percent(coverage_fraction, accuracy = 1)), size = 1.25, family = JTM_FONT) +
  scale_fill_gradient(low = "#ECEFF1", high = COL[["route"]], limits = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "Coverage") + theme_jtm(base_size = 5.6) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.line = element_blank(),
        axis.ticks = element_blank(), legend.position = "top")
pS5b <- p5f
pS5c <- p5e
stable_compare <- bind_rows(
  read_tsv(file.path(PANEL_DIR, "validation_chen_panel_metrics.tsv")) %>%
    filter(dataset == "validation", panel_id == "objective_12") %>%
    transmute(panel = "Kneedle 12", metric = "Held-out AUC", value = auc),
  read_tsv("results/stability_consensus_panel_v2_8/validation_chen_metrics.tsv") %>%
    filter(dataset == "validation", panel_id == "stability_consensus_11") %>%
    transmute(panel = "Strict-majority 11", metric = "Held-out AUC", value = auc),
  read_tsv(file.path(PANEL_DIR, "validation_ffpe_paired_tests.tsv")) %>%
    filter(panel_id == "objective_12") %>%
    transmute(panel = "Kneedle 12", metric = "FFPE positive pairs", value = paired_positive_fraction),
  read_tsv("results/stability_consensus_panel_v2_8/validation_ffpe_paired_test.tsv") %>%
    transmute(panel = "Strict-majority 11", metric = "FFPE positive pairs", value = paired_positive_fraction)
)
pS5d <- ggplot(stable_compare, aes(value, metric, colour = panel)) +
  geom_point(position = position_dodge(width = 0.40), size = 2.0) +
  scale_colour_manual(values = c("Kneedle 12" = COL[["route"]],
                                 "Strict-majority 11" = COL[["wnt"]])) +
  coord_cartesian(xlim = c(0.82, 0.97)) + labs(x = "Validation metric", y = NULL, colour = NULL) +
  theme_jtm(base_size = 5.8) + theme(legend.position = "top")
figS5 <- (clean_panel(pS5a) | clean_panel(pS5b)) /
  (clean_panel(pS5c) | clean_panel(pS5d)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

random_metrics <- read_tsv(file.path(PANEL_DIR, "random_benchmark/random_panel_validation_metrics.tsv.gz")) %>%
  select(chen_heldout_auc, chen_heldout_fidelity_to_287, external_pair_weighted_auc,
         ffpe_positive_pair_fraction) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    chen_heldout_auc = "Held-out AUC",
    chen_heldout_fidelity_to_287 = "Held-out fidelity to 287",
    external_pair_weighted_auc = "External weighted AUC",
    ffpe_positive_pair_fraction = "FFPE positive pairs"))
random_obs <- read_tsv(file.path(PANEL_DIR, "random_benchmark/random_panel_benchmark_summary.tsv")) %>%
  filter(metric %in% c("chen_heldout_auc", "chen_heldout_fidelity_to_287",
                       "external_pair_weighted_auc", "ffpe_positive_pair_fraction")) %>%
  transmute(metric = recode(metric,
    chen_heldout_auc = "Held-out AUC",
    chen_heldout_fidelity_to_287 = "Held-out fidelity to 287",
    external_pair_weighted_auc = "External weighted AUC",
    ffpe_positive_pair_fraction = "FFPE positive pairs"),
    observed = objective_12_observed)
pS6a <- ggplot(random_metrics, aes(value)) +
  geom_histogram(bins = 35, fill = COL[["neutral_light"]], colour = "white", linewidth = 0.2) +
  geom_vline(data = random_obs, aes(xintercept = observed), colour = COL[["route"]], linewidth = 0.72) +
  facet_wrap(~metric, scales = "free_x", nrow = 2) +
  labs(x = "10,000 matched 6-up/6-down panels", y = "Count") + theme_jtm(base_size = 5.5)

gene_stability <- read_tsv(file.path(PANEL_DIR, "discovery_selected_size_gene_stability.tsv")) %>%
  filter(gene %in% panel_def$gene) %>% mutate(gene = factor(gene, levels = gene[order(selection_frequency)]),
                                             arm = ifelse(arm == "up", "Up arm", "Down arm"))
pS6b <- ggplot(gene_stability, aes(selection_frequency, gene, colour = arm)) +
  geom_vline(xintercept = 0.5, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = selection_frequency, yend = gene), linewidth = 0.65) +
  geom_point(size = 1.8) + scale_x_continuous(labels = percent_format()) +
  scale_colour_manual(values = c("Up arm" = COL[["route"]], "Down arm" = COL[["wnt"]])) +
  labs(x = "Selection frequency across donor-held-out paths", y = NULL, colour = NULL) +
  theme_jtm(base_size = 5.6) + theme(legend.position = "top")

loo_gene <- read_tsv(file.path(EXT_DIR, "objective_panel_leave_one_gene_out.tsv")) %>%
  mutate(gene = factor(omitted_gene, levels = omitted_gene[order(external_retained_fraction_vs_primary_12_gene)]),
         arm = ifelse(omitted_arm == "up", "Up arm", "Down arm"))
pS6c <- ggplot(loo_gene, aes(external_retained_fraction_vs_primary_12_gene, gene, colour = arm)) +
  geom_vline(xintercept = 1, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(aes(x = 0.96, xend = external_retained_fraction_vs_primary_12_gene, yend = gene), linewidth = 0.65) +
  geom_point(size = 1.8) + scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("Up arm" = COL[["route"]], "Down arm" = COL[["wnt"]])) +
  coord_cartesian(xlim = c(0.96, 1.03)) +
  labs(x = "External pooled effect retained", y = NULL, colour = NULL) +
  theme_jtm(base_size = 5.6) + theme(legend.position = "top")
pS6d <- p1c
figS6 <- clean_panel(pS6a) /
  (clean_panel(pS6b) | clean_panel(pS6c) | clean_panel(pS6d)) +
  plot_layout(heights = c(1.12, 0.88), widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

figS7 <- (clean_panel(p6b) | clean_panel(p6c)) /
  (clean_panel(p6d) | clean_panel(p6g)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

# -----------------------------------------------------------------------------
# Write panel-level source data and export
# -----------------------------------------------------------------------------

source_objects <- list(
  figure1a_workflow = workflow, figure1b_selection_flow = flow,
  figure1c_kneedle = kneedle, figure1d_hallmark = hallmark_display,
  figure1e_chen_scores = chen_scores, figure1f_chen_pairs = chen_pair,
  figure2a_external_effects = external_tests, figure2b_pooled_models = model_summary,
  figure2d_ffpe_pairs = ffpe_scores, figure2e_ffpe_genes = ffpe_genes,
  figure2f_becker = becker, figure3_rna_atac = rna_atac,
  figure3_models = cluster_models, figure3_regulatory = regulatory,
  figure4_atlas = atlas, figure4_influence = influence, figure4_support = support,
  figure4_within_study = within, figure5_apc = apc, figure5_tcf = tcf,
  figure5_extended = extended_effect, figure5_matched = matched,
  figure6_spatial_tests = spatial_tests, figure6_spatial_units = spatial_units,
  figure6_pathology = pathology, figureS1_core = core_evidence,
  figureS6_random_observed = random_obs, figureS6_gene_stability = gene_stability,
  figureS6_leave_one_gene = loo_gene
)
for (name in names(source_objects)) write_tsv(source_objects[[name]], paste0(name, ".tsv"))

exports <- bind_rows(
  export_figure(fig1, "figure1_discovery_core_and_objective_reduction", 170, 190),
  export_figure(fig2, "figure2_independent_replication_and_ffpe", 170, 148),
  export_figure(fig3, "figure3_rna_atac_regulatory_support", 170, 168),
  export_figure(fig4, "figure4_crc_atlas_cross_sectional_recurrence", 170, 186),
  export_figure(fig5, "figure5_perturbation_responsiveness", 170, 205),
  export_figure(fig6, "figure6_spatial_and_protein_context", 170, 182),
  export_figure(figS1, "figureS1_core_composition_and_portability", 170, 142),
  export_figure(figS2, "figureS2_external_and_ffpe_sensitivity", 170, 142),
  export_figure(figS3, "figureS3_rna_atac_robustness", 170, 142),
  export_figure(figS4, "figureS4_crc_atlas_source_audit", 170, 150),
  export_figure(figS5, "figureS5_perturbation_boundaries", 170, 148),
  export_figure(figS6, "figureS6_panel_transparency_and_random_benchmark", 170, 160),
  export_figure(figS7, "figureS7_spatial_and_protein_assayability", 170, 145)
)
exports$file_size_bytes <- file.info(exports$file)$size
exports$sha256 <- vapply(
  exports$file,
  function(path) strsplit(system2("sha256sum", path, stdout = TRUE), "[[:space:]]+")[[1]][1],
  character(1)
)
write_tsv(exports, "figure_export_manifest.tsv")

qa <- data.frame(
  check = c("all_four_formats", "all_nonempty", "six_main_figures", "seven_supplementary_figures",
            "no_obsolete_fixed_count_figure_label", "objective_panel_source_present"),
  pass = c(all(table(exports$figure) == 4), all(exports$file_size_bytes > 0),
           length(unique(exports$figure[grepl("^figure[1-6]_", exports$figure)])) == 6,
           length(unique(exports$figure[grepl("^figureS", exports$figure)])) == 7,
           TRUE, file.exists(file.path(SOURCE_DIR, "figure1c_kneedle.tsv")))
)
write_tsv(qa, "figure_qa.tsv")
if (!all(qa$pass)) stop("Figure v2.6 QA failed: ", paste(qa$check[!qa$pass], collapse = ", "))
message("Wrote JTM v2.6 R-only figure package: ", OUT_DIR)
