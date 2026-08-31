#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")

root <- normalizePath(".", mustWork = TRUE)
out_dir <- file.path(root, "figures", "communications_biology_v5.0")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

# These scripts expose the audited pathway and exploratory WGCNA panels.
source("analysis/build_functional_architecture_figure4_v1.R")
source("analysis/build_exploratory_wgcna_integration_v2_1.R")

# Both source scripts retain their historical output directories; restore the
# v5 targets before composing the submission figures.
out_dir <- file.path(root, "figures", "communications_biology_v5.0")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

p_b_clean <- p_b +
  theme(axis.text.x = element_text(size = 4.8, angle = 38, hjust = 1))
p_c_clean <- p_c +
  guides(size = "none") +
  theme(axis.text.x = element_text(size = 5.0))

main_figure <- (p4a | p4b) / (p_b_clean | p_c_clean) +
  plot_layout(heights = c(1.03, 1.0), widths = c(1.12, 0.88)) +
  plot_annotation(tag_levels = "a")
main_figure <- tagged(main_figure)
export_cb_figure(
  main_figure, out_dir, "figure4_functional_module_architecture",
  width_mm = 178, height_mm = 150
)

supplement_figure <- p_a / p_d +
  plot_layout(heights = c(0.82, 1.0)) +
  plot_annotation(tag_levels = "a")
supplement_figure <- tagged(supplement_figure)
export_cb_figure(
  supplement_figure, out_dir, "figureS4_wgcna_structure_and_context",
  width_mm = 178, height_mm = 128
)

write_source(
  heatmap_data, source_dir,
  "figure4a_pathway_replication_heatmap.tsv"
)
write_source(
  curve_data, source_dir,
  "figure4b_running_enrichment.tsv"
)
write_source(
  association_plot_data, source_dir,
  "figure4c_module_rank_enrichment.tsv"
)
write_source(
  community_plot_data, source_dir,
  "figure4d_module_pathway_overlap.tsv"
)
write_source(
  dendrogram_data$segments, source_dir,
  "figureS4a_consensus_dendrogram_segments.tsv"
)
write_source(
  module_strip, source_dir,
  "figureS4a_consensus_module_strip.tsv"
)
write_source(
  evidence_plot_data, source_dir,
  "figureS4b_module_cross_context_evidence.tsv"
)

message("Communications Biology v5 Figure 4 and Supplementary Figure 4 exported")
