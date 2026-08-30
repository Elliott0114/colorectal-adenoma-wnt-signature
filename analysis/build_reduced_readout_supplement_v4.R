#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")
suppressPackageStartupMessages(library(tidyr))

root <- normalizePath(".", mustWork = TRUE)
state_root <- file.path(root, "results", "state_aware_program_v1")
revision_root <- file.path(root, "results", "state_shared_revision_v2")
out_dir <- file.path(root, "figures", "communications_biology_v3.0")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

candidate_summary <- read_tsv(file.path(
  state_root, "panel_derivation", "candidate_universe_summary.tsv"
)) %>%
  pivot_longer(
    c(strict_state_shared_genes, portable_protein_coding),
    names_to = "stage", values_to = "genes"
  ) %>%
  mutate(
    arm = factor(
      arm, levels = c("down", "up"),
      labels = c("Adenoma-down", "Adenoma-up")
    ),
    stage = factor(
      stage,
      levels = c("strict_state_shared_genes", "portable_protein_coding"),
      labels = c("State-shared", "Platform-measurable\nprotein-coding")
    )
  )

p5a <- ggplot(candidate_summary, aes(stage, genes, fill = arm)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(
    aes(label = genes), position = position_dodge(width = 0.72),
    vjust = -0.35, family = figure_font, size = 1.9
  ) +
  scale_fill_manual(values = c(
    "Adenoma-down" = figure_colours[["normal"]],
    "Adenoma-up" = figure_colours[["adenoma"]]
  )) +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Genes (log scale)") +
  theme_cb() +
  theme(legend.position = "top", axis.text.x = element_text(size = 5.4))

state_fidelity <- read_tsv(file.path(
  state_root, "panel_derivation", "discovery_state_specific_fidelity_curve.tsv"
)) %>%
  mutate(cell_type = factor(cell_type, levels = c("ABS", "GOB", "TAC")))

p5b <- ggplot(
  state_fidelity, aes(total_genes, oof_spearman, colour = cell_type)
) +
  geom_vline(
    xintercept = 8, linetype = 2, colour = figure_colours[["gold"]],
    linewidth = 0.5
  ) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.45) +
  scale_colour_manual(values = c(
    ABS = figure_colours[["adenoma"]],
    GOB = figure_colours[["normal"]],
    TAC = figure_colours[["purple"]]
  )) +
  labs(x = "Genes in balanced readout", y = "Donor-held-out fidelity") +
  theme_cb() +
  theme(legend.position = "top")

selected_genes <- c(
  "EPHB2", "REG1A", "LTBP1", "RNF43",
  "CALM2", "COX6C", "B2M", "ACAA2"
)
heldout_genes <- read_tsv(file.path(
  state_root, "heldout_validation", "heldout_compact_panel_gene_validation.tsv"
)) %>%
  select(gene, arm, logFC_ABS, logFC_GOB, logFC_TAC) %>%
  pivot_longer(starts_with("logFC_"), names_to = "state", values_to = "effect") %>%
  mutate(
    state = factor(
      sub("logFC_", "", state), levels = c("ABS", "GOB", "TAC")
    ),
    gene = factor(gene, levels = rev(selected_genes))
  )
heldout_limit <- max(abs(heldout_genes$effect), na.rm = TRUE)

p5c <- ggplot(heldout_genes, aes(state, gene, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(
    aes(label = ifelse(is.na(effect), "NA", sprintf("%.1f", effect))),
    family = figure_font, size = 2.0
  ) +
  scale_fill_gradient2(
    low = figure_colours[["normal"]], mid = "white",
    high = figure_colours[["adenoma"]], midpoint = 0,
    limits = c(-heldout_limit, heldout_limit),
    na.value = figure_colours[["pale"]], name = "Held-out effect"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = figure_font, base_size = 7) +
  theme(
    panel.grid = element_blank(), legend.position = "top",
    axis.text = element_text(size = 5.7),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

internal_random <- read_tsv(file.path(
  revision_root, "compact_rank", "random_eight_gene_benchmark.tsv"
)) %>%
  transmute(
    partition = "Donor-disjoint",
    fidelity = spearman_with_full_programme
  )
external_random <- read_tsv(file.path(
  revision_root, "external_rank", "random_eight_gene_benchmark.tsv"
)) %>%
  transmute(
    partition = "External-cohort median",
    fidelity = median_cohort_spearman
  )
random_benchmark <- bind_rows(internal_random, external_random) %>%
  mutate(partition = factor(
    partition, levels = c("Donor-disjoint", "External-cohort median")
  ))
internal_summary <- read_tsv(file.path(
  revision_root, "compact_rank", "random_eight_gene_benchmark_summary.tsv"
)) %>%
  transmute(
    partition = "Donor-disjoint",
    observed = observed_compact_spearman, q95 = random_q95
  )
external_summary <- read_tsv(file.path(
  revision_root, "external_rank", "random_eight_gene_benchmark_summary.tsv"
)) %>%
  transmute(
    partition = "External-cohort median",
    observed = observed_median_cohort_spearman, q95 = random_q95
  )
benchmark_summary <- bind_rows(internal_summary, external_summary) %>%
  mutate(partition = factor(
    partition, levels = levels(random_benchmark$partition)
  ))

p5d <- ggplot(random_benchmark, aes(fidelity)) +
  geom_histogram(
    bins = 32, fill = figure_colours[["line"]], colour = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    data = benchmark_summary, aes(xintercept = observed),
    colour = figure_colours[["adenoma"]], linewidth = 0.65
  ) +
  geom_vline(
    data = benchmark_summary, aes(xintercept = q95),
    colour = figure_colours[["normal"]], linetype = 2, linewidth = 0.65
  ) +
  facet_wrap(~partition, ncol = 1, scales = "free_y") +
  labs(x = "Direction-balanced random-panel fidelity", y = "Count") +
  theme_cb()

figure <- (p5a | p5b) / (p5c | p5d) +
  plot_annotation(tag_levels = "a")
export_cb_figure(
  tagged(figure), out_dir, "figureS5_reduced_readout_audit",
  width_mm = 178, height_mm = 150
)

write_source(candidate_summary, source_dir, "figureS5a_candidate_summary.tsv")
write_source(state_fidelity, source_dir, "figureS5b_state_fidelity.tsv")
write_source(heldout_genes, source_dir, "figureS5c_heldout_gene_effects.tsv")
write_source(random_benchmark, source_dir, "figureS5d_random_benchmarks.tsv")
write_source(benchmark_summary, source_dir, "figureS5d_benchmark_summary.tsv")

message("Simplified reduced-readout Supplementary Figure 5 exported")
