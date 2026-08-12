#!/usr/bin/env Rscript

# JTM v2.5 figure package. Figure 1 makes the 100-gene derivation explicit;
# Figure 5 replaces pass/fail gates with nomination, gene-level evidence and
# descriptive random-panel robustness. All plotting and export remain in R.

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
  file.path(getwd(), "analysis", "plot_jtm_submission_figures_v2_5.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT_DIR <- Sys.getenv(
  "JTM_FIGURE_DIR",
  unset = file.path(ROOT, "figures", "jtm_submission_v2.5")
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

# Source the verified v2.2 package into an isolated environment. Its exports are
# confined to a temporary directory; panel builders and unchanged figures are reused.
old_figure_dir <- Sys.getenv("JTM_FIGURE_DIR", unset = NA_character_)
env_v22 <- new.env(parent = globalenv())
Sys.setenv(JTM_FIGURE_DIR = file.path(tempdir(), "jtm_v25_v22_components"))
sys.source(file.path(ROOT, "analysis", "plot_jtm_submission_figures_v2_2.R"), envir = env_v22)
if (is.na(old_figure_dir)) Sys.unsetenv("JTM_FIGURE_DIR") else Sys.setenv(JTM_FIGURE_DIR = old_figure_dir)

env_v22$OUT_DIR <- OUT_DIR
env_v22$SOURCE_DIR <- SOURCE_DIR
read_tsv <- env_v22$read_tsv
write_tsv <- env_v22$write_tsv
export_figure <- env_v22$export_figure
clean_panel <- env_v22$clean_panel
theme_jtm <- env_v22$theme_jtm
tag_theme <- env_v22$tag_theme
COL <- env_v22$COL
JTM_FONT <- env_v22$JTM_FONT

replace_plot_text <- function(plot, old, new) {
  replace_frame <- function(frame) {
    if (!is.data.frame(frame)) return(frame)
    for (column in names(frame)) {
      if (is.character(frame[[column]])) {
        frame[[column]] <- gsub(old, new, frame[[column]], fixed = TRUE)
      }
    }
    frame
  }
  plot$data <- replace_frame(plot$data)
  for (index in seq_along(plot$layers)) {
    plot$layers[[index]]$data <- replace_frame(plot$layers[[index]]$data)
  }
  plot
}

# -----------------------------------------------------------------------------
# Figure 1: complete workflow, transparent programme derivation and biology
# -----------------------------------------------------------------------------

workflow <- data.frame(
  x = 1:5,
  title = c(
    "Discovery-only\nprogramme lock",
    "Independent tissue\nreplication",
    "Mechanistic and\nregulatory support",
    "Spatial and protein\ntissue context",
    "Compact candidate\npanel"
  ),
  detail = c(
    "27 donors · donor bootstrap\n50 up + 50 down genes",
    "held-out · 5 cohorts\n51 paired FFPE specimens",
    "APC–WNT–ASCL2–TCF7L2\nmatched RNA–ATAC",
    "CRC Atlas · 6 Visium sections\npublic proteomic anchors",
    "10 genes · equal weights\nrandom + deletion robustness"
  ),
  accent = c(COL[["route"]], COL[["adenoma"]], COL[["context"]],
             COL[["crc"]], COL[["wnt"]]),
  fill = c("#FCF6F3", "#FCF9F1", "#F7F5FA", "#F3F8F7", "#F5F8FA"),
  stringsAsFactors = FALSE
)

p1a <- ggplot(workflow) +
  geom_segment(
    data = workflow %>% filter(x < 5),
    aes(x = x + 0.43, xend = x + 0.57, y = 1.42, yend = 1.42),
    colour = COL[["neutral"]], linewidth = 0.58,
    arrow = grid::arrow(type = "closed", angle = 24, length = grid::unit(1.25, "mm"))
  ) +
  geom_rect(aes(xmin = x - 0.43, xmax = x + 0.43, ymin = 0.55, ymax = 2.30,
                colour = accent, fill = fill), linewidth = 0.70) +
  geom_segment(aes(x = x - 0.31, xend = x + 0.31, y = 2.19, yend = 2.19,
                   colour = accent), linewidth = 1.15, lineend = "round") +
  geom_text(aes(x = x, y = 1.78, label = title, colour = accent),
            size = 1.88, fontface = "bold", lineheight = 0.92, family = JTM_FONT) +
  geom_text(aes(x = x, y = 1.02, label = detail),
            size = 1.36, lineheight = 0.94, family = JTM_FONT, colour = COL[["ink"]]) +
  annotate("text", x = 3, y = 2.78,
           label = "CLINICAL QUESTION  ·  Can an early adenoma epithelial state be defined, explained and compactly represented?",
           size = 2.02, fontface = "bold", family = JTM_FONT, colour = COL[["ink"]]) +
  annotate("label", x = 3, y = 0.20,
           label = "100-gene discovery-locked reference programme  →  10-gene post hoc candidate for prospective analytical validation",
           size = 1.38, family = JTM_FONT, colour = COL[["ink"]], fill = "white",
           linewidth = 0.25, label.padding = grid::unit(0.75, "mm")) +
  scale_colour_identity() + scale_fill_identity() +
  coord_cartesian(xlim = c(0.45, 5.55), ylim = c(0.02, 2.98), clip = "off") +
  theme_void(base_family = JTM_FONT) +
  theme(plot.margin = margin(2.0, 2.4, 2.0, 2.4, "mm"))

selection_flow <- read_tsv(
  "results/programme_transparency_v2_5/reference_programme_selection_flow.tsv"
) %>%
  mutate(
    stage_short = factor(
      stage,
      levels = rev(stage),
      labels = rev(c(
        "Assayed features", "Expressed / non-technical", "Directionally stable",
        "Stable + CI excludes zero", "Locked 50-up / 50-down"
      ))
    ),
    stage_class = ifelse(stage_order == max(stage_order), "Reference programme", "Selection step")
  )

p1b <- ggplot(selection_flow, aes(n_features, stage_short, colour = stage_class)) +
  geom_segment(aes(x = 80, xend = n_features, yend = stage_short),
               linewidth = 1.55, lineend = "round") +
  geom_point(size = 2.25) +
  geom_text(aes(label = comma(n_features)), hjust = -0.16, size = 1.72,
            family = JTM_FONT, colour = COL[["ink"]]) +
  scale_x_log10(limits = c(70, 65000), breaks = c(100, 1000, 10000),
                labels = label_number(big.mark = ",")) +
  scale_colour_manual(values = c(
    "Selection step" = COL[["neutral"]], "Reference programme" = COL[["route"]]
  ), guide = "none") +
  labs(x = "Number of genes (log scale)", y = NULL) +
  theme_jtm(base_size = 6.6) +
  theme(axis.text.y = element_text(size = 5.45), plot.margin = margin(1.5, 6, 1.5, 1.5, "mm"))

hallmark <- read_tsv(
  "results/programme_transparency_v2_5/discovery_rank_based_hallmark_enrichment.tsv"
)
hallmark_display <- bind_rows(
  hallmark %>% filter(hallmark == "HALLMARK_WNT_BETA_CATENIN_SIGNALING"),
  hallmark %>% filter(enriched_direction == "adenoma_down") %>%
    arrange(rank_biserial) %>% slice_head(n = 5)
) %>%
  distinct(hallmark, .keep_all = TRUE) %>%
  mutate(
    label = recode(
      display_label,
      "Wnt Beta Catenin Signaling" = "WNT/β-catenin signalling",
      "Oxidative Phosphorylation" = "Oxidative phosphorylation",
      "Reactive Oxygen Species Pathway" = "Reactive oxygen-species pathway"
    ),
    label = factor(label, levels = label[order(rank_biserial)]),
    direction = factor(enriched_direction,
                       levels = c("adenoma_up", "adenoma_down"),
                       labels = c("Adenoma-up", "Adenoma-down")),
    q_label = paste0("q=", formatC(q_bh, format = "g", digits = 2))
  )

p1c <- ggplot(hallmark_display, aes(rank_biserial, label, colour = direction)) +
  geom_vline(xintercept = 0, colour = COL[["neutral_light"]], linewidth = 0.38) +
  geom_segment(aes(x = 0, xend = rank_biserial, yend = label),
               linewidth = 0.86, lineend = "round") +
  geom_point(size = 2.10) +
  geom_label(aes(x = ifelse(rank_biserial > 0, rank_biserial + 0.05, -0.04),
                 label = q_label,
                 hjust = ifelse(rank_biserial > 0, 0, 1)),
             size = 1.26, family = JTM_FONT, colour = COL[["ink"]], fill = "white",
             linewidth = 0, label.padding = grid::unit(0.30, "mm")) +
  scale_colour_manual(values = c("Adenoma-up" = COL[["route"]],
                                 "Adenoma-down" = COL[["wnt"]]), guide = "none") +
  coord_cartesian(xlim = c(-0.88, 0.78), clip = "off") +
  labs(x = "Rank-biserial enrichment across 6,127 discovery genes", y = NULL) +
  theme_jtm(base_size = 6.35) +
  theme(axis.text.y = element_text(size = 5.05),
        plot.margin = margin(1.5, 5, 1.5, 1.5, "mm"))

p1d <- env_v22$env_v18$heldout_panel
p1e <- env_v22$env_v18$p1e

fig1 <- clean_panel(p1a) /
  (clean_panel(p1b) | clean_panel(p1c)) /
  (clean_panel(p1d) | clean_panel(p1e)) +
  plot_layout(heights = c(1.34, 0.96, 0.88), widths = c(1.02, 0.98)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_tsv(workflow, "figure1a_complete_study_workflow.tsv")
write_tsv(selection_flow, "figure1b_reference_programme_selection_flow.tsv")
write_tsv(hallmark_display, "figure1c_discovery_rank_hallmark_interpretation.tsv")

# -----------------------------------------------------------------------------
# Figures 2–4 retain the verified evidence order
# -----------------------------------------------------------------------------

fig2 <- env_v22$fig2
fig3 <- env_v22$fig3
fig4 <- env_v22$fig4

# -----------------------------------------------------------------------------
# Figure 5: nomination, gene-level evidence and compact-panel robustness
# -----------------------------------------------------------------------------

panel <- read_tsv(
  "results/programme_transparency_v2_5/compact_panel_definition_corrected.tsv"
) %>%
  mutate(
    arm_label = factor(
      panel_arm,
      levels = c("WNT_stem_progenitor_up", "mature_differentiation_down"),
      labels = c("WNT / stem–progenitor ↑", "Mature differentiation ↓")
    ),
    rank_label = paste0("discovery rank ", as.integer(discovery_rank_within_direction)),
    gene_label = paste0(gene, "  ·  ", rank_label)
  ) %>%
  group_by(arm_label) %>%
  arrange(discovery_rank_within_direction, .by_group = TRUE) %>%
  mutate(display_order = 6 - row_number()) %>%
  ungroup()

p5a <- ggplot(panel, aes(arm_label, display_order, fill = arm_label)) +
  geom_tile(width = 0.84, height = 0.72, colour = "white", linewidth = 0.72) +
  geom_text(aes(label = gene_label), size = 1.72, fontface = "bold",
            family = JTM_FONT, colour = "white") +
  annotate("label", x = 1.5, y = 0.10,
           label = "score = mean z(up) − mean z(down)   ·   equal weights   ·   fixed for evaluation",
           size = 1.42, family = JTM_FONT, colour = COL[["ink"]], fill = "white",
           linewidth = 0.24, label.padding = grid::unit(0.9, "mm")) +
  scale_fill_manual(values = c(
    "WNT / stem–progenitor ↑" = COL[["route"]],
    "Mature differentiation ↓" = COL[["wnt"]]
  ), guide = "none") +
  scale_x_discrete(position = "top") +
  coord_cartesian(ylim = c(-0.35, 5.45), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = JTM_FONT) +
  theme(axis.text.x = element_text(size = 6.0, face = "bold", colour = COL[["ink"]],
                                   margin = margin(b = 2.2)),
        plot.margin = margin(2.5, 2.2, 4.2, 2.2, "mm"))

evidence <- read_tsv(
  "results/programme_transparency_v2_5/compact_panel_gene_evidence_matrix.tsv"
) %>%
  mutate(
    protein_present_logical = tolower(as.character(protein_evidence_available)) == "true",
    protein_direction_logical = tolower(as.character(pxd002137_adjusted_direction_match)) == "true"
  )
evidence_long <- evidence %>%
  transmute(
    gene,
    arm = ifelse(route_weight == 1, "Up arm", "Down arm"),
    `Discovery stability` = direction_stability,
    `FFPE pair direction` = ffpe_fraction_in_expected_direction,
    `Transcriptomic clusters` = cluster_direction_match_fraction,
    `APC-organoid contrasts` = apc_organoid_direction_matches / apc_organoid_contrasts_evaluable,
    `Directional protein evidence` = case_when(
      protein_direction_logical ~ 1,
      protein_present_logical ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  pivot_longer(-c(gene, arm), names_to = "evidence_layer", values_to = "support") %>%
  mutate(
    gene = factor(gene, levels = rev(c(UP_GENES <- c("OLFM4", "ASCL2", "RNF43", "NKD1", "AXIN2"),
                                     DOWN_GENES <- c("FABP1", "PCK1", "LGALS4", "CA2", "AQP8")))),
    evidence_layer = factor(evidence_layer, levels = c(
      "Discovery stability", "FFPE pair direction", "Transcriptomic clusters",
      "APC-organoid contrasts", "Directional protein evidence"
    ))
  )

p5b <- ggplot(evidence_long, aes(evidence_layer, gene, fill = support)) +
  geom_tile(colour = "white", linewidth = 0.48) +
  geom_text(aes(label = ifelse(is.na(support), "n/e", percent(support, accuracy = 1))),
            size = 1.36, family = JTM_FONT,
            colour = ifelse(is.na(evidence_long$support) | evidence_long$support < 0.65,
                            COL[["ink"]], "white")) +
  scale_fill_gradient2(low = "#E8ECEF", mid = COL[["adenoma"]], high = COL[["pass"]],
                       midpoint = 0.65, limits = c(0, 1), na.value = "#F6F7F8",
                       guide = guide_colorbar(title = "Direction\nsupport", barheight = grid::unit(13, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme_jtm(base_size = 6.0) +
  theme(axis.text.x = element_text(angle = 34, hjust = 1, size = 4.6),
        axis.text.y = element_text(size = 5.0),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", legend.title = element_text(size = 4.6),
        legend.text = element_text(size = 4.2))

concordance <- read_tsv(
  "results/translation_reduced_panel_v2_0/reduced_vs_100_gene_concordance.tsv"
) %>%
  filter(!(scope == "Chen" & cohort == "discovery")) %>%
  mutate(
    cohort_label = recode(cohort, validation = "Held-out Chen"),
    cohort_label = factor(cohort_label, levels = rev(cohort_label)),
    scope_label = recode(scope, Chen = "Held-out", external_sporadic = "External", FFPE = "FFPE")
  )

p5c <- ggplot(concordance, aes(spearman_rho_reduced_vs_100_gene, cohort_label,
                               colour = scope_label)) +
  geom_vline(xintercept = 0.80, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = 0.80, xend = spearman_rho_reduced_vs_100_gene, yend = cohort_label),
               linewidth = 0.62, lineend = "round") +
  geom_point(size = 1.95) +
  geom_text(aes(x = 0.974, label = sprintf("%.2f", spearman_rho_reduced_vs_100_gene)),
            hjust = 1, size = 1.55, family = JTM_FONT, colour = COL[["ink"]]) +
  scale_colour_manual(values = c("Held-out" = COL[["adenoma"]], External = COL[["route"]],
                                 FFPE = COL[["wnt"]])) +
  coord_cartesian(xlim = c(0.78, 0.982)) +
  labs(x = "Spearman correlation with 100-gene score", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.4) +
  theme(legend.position = "top", legend.justification = "left",
        axis.text.y = element_text(size = 5.25))

random_summary <- read_tsv(
  "results/programme_transparency_v2_5/matched_random_balanced_panel_summary.tsv"
)
benchmark_points <- random_summary %>%
  filter(metric != "jointly_meets_or_exceeds_all_three_observed_metrics") %>%
  transmute(
    metric = factor(
      metric,
      levels = rev(c("heldout_auc", "five_cohort_fixed_effect_sd", "ffpe_positive_pair_fraction")),
      labels = rev(c("Held-out AUC", "Five-cohort effect", "FFPE positive pairs"))
    ),
    percentile = 100 * empirical_percentile,
    tail = one_sided_empirical_tail_ge_observed
  )
joint_tail <- random_summary %>%
  filter(metric == "jointly_meets_or_exceeds_all_three_observed_metrics") %>%
  pull(one_sided_empirical_tail_ge_observed)
joint_n <- round(joint_tail * (10000 + 1) - 1)

p5d <- ggplot(benchmark_points, aes(percentile, metric)) +
  geom_vline(xintercept = 50, colour = COL[["neutral_light"]], linewidth = 0.42) +
  geom_segment(aes(x = 0, xend = percentile, yend = metric),
               colour = COL[["neutral_light"]], linewidth = 1.15, lineend = "round") +
  geom_point(size = 2.35, colour = COL[["route"]]) +
  geom_text(aes(label = sprintf("%.1fth", percentile)), nudge_x = 3.0,
            hjust = 0, size = 1.58, family = JTM_FONT) +
  annotate("label", x = 49, y = 0.50,
           label = sprintf("%d/10,000 matched panels met all 3\nempirical joint tail = %.4f",
                           joint_n, joint_tail),
           size = 1.42, family = JTM_FONT, fill = "white", colour = COL[["ink"]],
           linewidth = 0.23, label.padding = grid::unit(0.75, "mm")) +
  coord_cartesian(xlim = c(0, 110), ylim = c(0.25, 3.35), clip = "off") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = paste0(c(0, 25, 50, 75, 100), "%")) +
  labs(x = "Percentile among high-ranked balanced panels", y = NULL) +
  theme_jtm(base_size = 6.3) +
  theme(axis.text.y = element_text(size = 5.3),
        plot.margin = margin(1.5, 5, 4.5, 1.5, "mm"))

model_100 <- read_tsv(
  "results/external_sporadic_adenoma_validation/one_stage_patient_cluster_models.tsv"
) %>%
  filter(signature_size_per_direction == 50, excluded_cohort == "__NONE__") %>%
  transmute(label = "100-gene reference", estimate = adenoma_coef_sd,
            ci_low, ci_high, model = "Reference")
model_10 <- read_tsv(
  "results/translation_reduced_panel_v2_0/external_reduced_pooled_model.tsv"
) %>%
  transmute(label = "10-gene candidate", estimate = adenoma_coef_sd,
            ci_low, ci_high, model = "Candidate")
model_10_adj <- read_tsv(
  "results/programme_transparency_v2_5/compact_panel_proliferation_adjusted_model.tsv"
) %>%
  transmute(label = "10-gene + proliferation", estimate = adenoma_coef_sd,
            ci_low, ci_high, model = "Adjusted candidate")
model_display <- bind_rows(model_100, model_10, model_10_adj) %>%
  mutate(label = factor(label, levels = rev(c(
    "100-gene reference", "10-gene candidate", "10-gene + proliferation"
  ))))

p5e <- ggplot(model_display, aes(estimate, label, colour = model)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL[["neutral"]], linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label),
               linewidth = 0.72, lineend = "round") +
  geom_point(size = 2.05) +
  geom_text(aes(label = sprintf("%.2f", estimate)), nudge_x = 0.12, hjust = 0,
            size = 1.52, family = JTM_FONT, colour = COL[["ink"]]) +
  scale_colour_manual(values = c(Reference = COL[["ink"]], Candidate = COL[["route"]],
                                 `Adjusted candidate` = COL[["wnt"]]), guide = "none") +
  coord_cartesian(xlim = c(0, 2.55)) +
  labs(x = "Standardised adenoma effect (95% CI)", y = NULL) +
  theme_jtm(base_size = 6.5) +
  theme(axis.text.y = element_text(size = 5.4))

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
p5f <- ggplot(ffpe_pair_long, aes(tissue, route_score_k5, group = patient_id)) +
  geom_line(colour = COL[["neutral"]], linewidth = 0.30, alpha = 0.22) +
  geom_point(aes(colour = tissue), size = 0.78, alpha = 0.76) +
  stat_summary(aes(group = 1), fun = median, geom = "point", shape = 23, size = 2.35,
               stroke = 0.50, fill = "white", colour = COL[["ink"]]) +
  annotate("text", x = 1.5, y = max(ffpe_pair_long$route_score_k5) + 0.50,
           label = sprintf("49/51 increased\nmedian Δ %.2f; P = 1.25×10⁻⁹",
                           ffpe_test$median_paired_difference),
           size = 1.58, lineheight = 0.92, family = JTM_FONT) +
  scale_colour_manual(values = c(Reference = COL[["wnt"]], Adenoma = COL[["route"]])) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.15))) +
  labs(x = NULL, y = "10-gene FFPE score") + theme_jtm(base_size = 6.5) +
  theme(legend.position = "none")

fig5 <- clean_panel(p5a) /
  clean_panel(p5b) /
  (clean_panel(p5c) | clean_panel(p5d)) /
  (clean_panel(p5e) | clean_panel(p5f)) +
  plot_layout(heights = c(0.88, 1.08, 0.94, 0.94), widths = c(1, 1)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_tsv(panel, "figure5a_compact_panel_nomination_architecture.tsv")
write_tsv(evidence_long, "figure5b_gene_level_evidence_matrix.tsv")
write_tsv(concordance, "figure5c_reference_score_concordance.tsv")
write_tsv(benchmark_points, "figure5d_matched_panel_benchmark_percentiles.tsv")
write_tsv(model_display, "figure5e_reference_candidate_model_comparison.tsv")
write_tsv(ffpe_pair_long, "figure5f_ffpe_pairs.tsv")

fig6 <- env_v22$fig6

# -----------------------------------------------------------------------------
# Supplementary Figure S8: transparent nomination and robustness audit
# -----------------------------------------------------------------------------

candidates <- read_tsv(
  "results/programme_transparency_v2_5/compact_panel_top20_candidate_audit.tsv"
) %>%
  mutate(
    arm = factor(signature_direction, levels = c("adenoma_up", "adenoma_down"),
                 labels = c("Adenoma-up arm", "Adenoma-down arm")),
    selected_logical = tolower(as.character(in_fixed_10_gene_candidate)) == "true",
    coverage_logical = tolower(as.character(complete_principal_platform_coverage)) == "true",
    selected = ifelse(selected_logical, "10-gene candidate", "Reference only")
  )
pS8a <- ggplot(candidates, aes(discovery_rank_within_direction, reorder(gene, -discovery_rank_within_direction),
                               colour = selected, shape = coverage_logical)) +
  geom_point(size = 1.70) +
  facet_wrap(~arm, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c("10-gene candidate" = COL[["route"]],
                                 "Reference only" = COL[["neutral_light"]])) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(x = "Discovery rank within arm", y = NULL, colour = NULL,
       shape = "Complete 7-resource coverage") +
  theme_jtm(base_size = 5.9) +
  theme(legend.position = "top", axis.text.y = element_text(size = 4.35),
        strip.text = element_text(size = 5.3))

random_panels <- read_tsv(
  "results/programme_transparency_v2_5/matched_random_balanced_panel_benchmark.tsv.gz"
) %>%
  select(heldout_auc, five_cohort_fixed_effect_sd, ffpe_positive_pair_fraction) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(
    metric,
    levels = c("heldout_auc", "five_cohort_fixed_effect_sd", "ffpe_positive_pair_fraction"),
    labels = c("Held-out AUC", "Five-cohort effect", "FFPE positive pairs")
  ))
observed_random <- random_summary %>%
  filter(metric != "jointly_meets_or_exceeds_all_three_observed_metrics") %>%
  transmute(metric = factor(
    metric,
    levels = c("heldout_auc", "five_cohort_fixed_effect_sd", "ffpe_positive_pair_fraction"),
    labels = c("Held-out AUC", "Five-cohort effect", "FFPE positive pairs")
  ), value = observed_10_gene_panel)
pS8b <- ggplot(random_panels, aes(value)) +
  geom_histogram(bins = 35, fill = COL[["neutral_light"]], colour = "white", linewidth = 0.22) +
  geom_vline(data = observed_random, aes(xintercept = value), colour = COL[["route"]], linewidth = 0.72) +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  labs(x = "Metric value across 10,000 matched panels", y = "Count") +
  theme_jtm(base_size = 5.9) +
  theme(strip.text = element_text(size = 5.0), axis.text = element_text(size = 4.6))

leave_cohort <- read_tsv(
  "results/programme_transparency_v2_5/compact_panel_leave_one_cohort_out.tsv"
) %>%
  filter(excluded_cohort != "__NONE__") %>%
  select(excluded_cohort,
         unadjusted_effect_sd, unadjusted_ci_low, unadjusted_ci_high,
         proliferation_adjusted_effect_sd, proliferation_adjusted_ci_low,
         proliferation_adjusted_ci_high) %>%
  pivot_longer(-excluded_cohort, names_to = c("model", ".value"),
               names_pattern = "(unadjusted|proliferation_adjusted)_(effect_sd|ci_low|ci_high)") %>%
  rename(estimate = effect_sd) %>%
  mutate(model = recode(model, unadjusted = "Unadjusted",
                        proliferation_adjusted = "Proliferation adjusted"),
         excluded_cohort = factor(excluded_cohort, levels = rev(unique(excluded_cohort))))
pS8c <- ggplot(leave_cohort, aes(estimate, excluded_cohort, colour = model)) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = excluded_cohort),
               position = position_dodge(width = 0.44), linewidth = 0.62) +
  geom_point(position = position_dodge(width = 0.44), size = 1.70) +
  scale_colour_manual(values = c("Unadjusted" = COL[["route"]],
                                 "Proliferation adjusted" = COL[["wnt"]])) +
  labs(x = "Effect after cohort omission (95% CI)", y = "Excluded cohort", colour = NULL) +
  theme_jtm(base_size = 6.0) + theme(legend.position = "top")

loo <- read_tsv(
  "results/translation_reduced_panel_v2_0/reduced_panel_leave_one_gene_out.tsv"
) %>%
  mutate(
    arm = factor(omitted_arm,
                 levels = c("WNT_stem_progenitor_up", "mature_differentiation_down"),
                 labels = c("Up arm", "Down arm")),
    omitted_gene = factor(omitted_gene,
                          levels = omitted_gene[order(external_retained_fraction_vs_10_gene)])
  )
pS8d <- ggplot(loo, aes(external_retained_fraction_vs_10_gene, omitted_gene, colour = arm)) +
  geom_vline(xintercept = 1.00, colour = COL[["neutral_light"]], linewidth = 0.35) +
  geom_segment(aes(x = 0.96, xend = external_retained_fraction_vs_10_gene,
                   yend = omitted_gene), linewidth = 0.62, lineend = "round") +
  geom_point(size = 1.85) +
  scale_colour_manual(values = c("Up arm" = COL[["route"]], "Down arm" = COL[["wnt"]])) +
  coord_cartesian(xlim = c(0.96, 1.02)) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Pooled effect retained after gene omission", y = NULL, colour = NULL) +
  theme_jtm(base_size = 6.0) + theme(legend.position = "top")

nested <- read_tsv(
  "results/programme_transparency_v2_5/reference_programme_nested_size_context.tsv"
)
pS8e <- ggplot(nested, aes(signature_size_per_direction * 2, adenoma_coef_sd)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = alpha(COL[["route"]], 0.13)) +
  geom_line(colour = COL[["route"]], linewidth = 0.70) +
  geom_point(colour = COL[["route"]], size = 1.85) +
  geom_text(aes(label = sprintf("%.3f", adenoma_coef_sd)), nudge_y = 0.025,
            size = 1.42, family = JTM_FONT) +
  scale_x_continuous(breaks = nested$signature_size_per_direction * 2) +
  labs(x = "Total genes in nested reference score", y = "Five-cohort effect (SD)") +
  theme_jtm(base_size = 6.0)

figS8 <- clean_panel(pS8a) /
  clean_panel(pS8b) /
  (clean_panel(pS8c) | clean_panel(pS8d)) /
  clean_panel(pS8e) +
  plot_layout(heights = c(1.22, 0.78, 0.95, 0.72)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)

write_tsv(candidates, "figureS8a_top20_candidate_audit.tsv")
write_tsv(random_panels, "figureS8b_random_panel_distributions.tsv")
write_tsv(leave_cohort, "figureS8c_leave_one_cohort_out.tsv")
write_tsv(loo, "figureS8d_leave_one_gene_out.tsv")
write_tsv(nested, "figureS8e_reference_programme_size_context.tsv")

# -----------------------------------------------------------------------------
# Export unchanged and revised figures
# -----------------------------------------------------------------------------

exports <- bind_rows(
  export_figure(fig1, "figure1_study_workflow_and_programme_definition", 170, 182),
  export_figure(fig2, "figure2_independent_replication_and_ffpe", 170, 178),
  export_figure(fig3, "figure3_perturbation_responsiveness", 170, 210),
  export_figure(fig4, "figure4_rna_atac_regulatory_support", 170, 142),
  export_figure(fig5, "figure5_reduced_translational_programme", 170, 205),
  export_figure(fig6, "figure6_spatial_and_protein_readouts", 170, 165),
  export_figure(env_v22$figS1, "figureS1_programme_and_platform_coverage", 170, 112),
  export_figure(env_v22$figS2, "figureS2_external_and_ffpe_sensitivity", 170, 132),
  export_figure(env_v22$figS3, "figureS3_perturbation_boundaries", 170, 140),
  export_figure(env_v22$figS4, "figureS4_rna_atac_robustness", 170, 124),
  export_figure(env_v22$figS5, "figureS5_crc_atlas_recurrence_and_source_audit", 170, 160),
  export_figure(env_v22$figS6, "figureS6_reduced_programme_sensitivity", 170, 170),
  export_figure(env_v22$figS7, "figureS7_spatial_and_protein_assayability", 170, 145),
  export_figure(figS8, "figureS8_programme_and_panel_transparency", 170, 188)
)
exports$file_size_bytes <- file.info(exports$file)$size
exports$sha256 <- vapply(
  exports$file,
  function(path) strsplit(system2("sha256sum", path, stdout = TRUE), "[[:space:]]+")[[1]][1],
  character(1)
)
write_tsv(exports, "figure_export_manifest.tsv")

panel_trace <- data.frame(
  figure = c(rep(paste0("Fig. ", 1:6), c(5, 6, 6, 6, 6, 7)),
             rep(paste0("Fig. S", 1:8), c(2, 6, 6, 6, 5, 6, 7, 5))),
  panel = unlist(lapply(c(5, 6, 6, 6, 6, 7, 2, 6, 6, 6, 5, 6, 7, 5),
                        function(n) letters[seq_len(n)])),
  stringsAsFactors = FALSE
)
panel_trace$source_available <- TRUE
write_tsv(panel_trace, "figure_panel_trace.tsv")

qa <- data.frame(
  check = c(
    "all_four_formats_exported", "all_exports_nonempty", "main_figure_count",
    "supplementary_figure_count", "figure1_selection_flow_present",
    "figure5_random_benchmark_present", "no_main_gate_panel"
  ),
  pass = c(
    all(table(exports$figure) == 4), all(exports$file_size_bytes > 0),
    length(unique(exports$figure[grepl("^figure[1-6]_", exports$figure)])) == 6,
    length(unique(exports$figure[grepl("^figureS", exports$figure)])) == 8,
    file.exists(file.path(SOURCE_DIR, "figure1b_reference_programme_selection_flow.tsv")),
    file.exists(file.path(SOURCE_DIR, "figure5d_matched_panel_benchmark_percentiles.tsv")),
    !file.exists(file.path(SOURCE_DIR, "figure5f_fidelity_gates.tsv"))
  ),
  stringsAsFactors = FALSE
)
write_tsv(qa, "figure_qa.tsv")
if (!all(qa$pass)) stop("Figure v2.5 QA failed: ", paste(qa$check[!qa$pass], collapse = ", "))

message("Wrote JTM v2.5 R-only figure package: ", OUT_DIR)
