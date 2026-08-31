#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")
suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
  library(tidyr)
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
set.seed(20260831)

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
source_root <- file.path(result_root, "functional_architecture_v1")
run_root_env <- Sys.getenv("STATE_AWARE_MODULE_RUN_ROOT", unset = "")
run_root <- if (nzchar(run_root_env)) {
  normalizePath(run_root_env, mustWork = TRUE)
} else {
  file.path(result_root, "functional_architecture_exploratory_v2")
}
wgcna_root <- file.path(source_root, "consensus_wgcna")
integration_root <- file.path(run_root, "integration_summary")
run_label <- basename(run_root)
figure_stem <- if (grepl("v2_1$", run_label)) {
  "exploratory_wgcna_module_integration_v2_1"
} else {
  "exploratory_wgcna_module_integration_v2"
}
figure_root <- file.path(
  root, "figures", "communications_biology_v3.0", run_label
)
source_data_root <- file.path(figure_root, "source_data")
dir.create(integration_root, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_root, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  audit_object = file.path(
    wgcna_root, "audit", "consensus_network_audit_object.rds"
  ),
  modules = file.path(source_root, "consensus_modules.tsv"),
  rank_enrichment = file.path(wgcna_root, "module_rank_enrichment.tsv"),
  preservation = file.path(wgcna_root, "module_preservation.tsv"),
  community_overlap = file.path(
    wgcna_root, "module_pathway_community_overlap.tsv"
  ),
  validation = file.path(run_root, "module_validation.tsv"),
  external = file.path(
    run_root, "module_external_validation", "module_external_cohort_effects.tsv"
  ),
  ffpe = file.path(
    run_root, "module_external_validation", "module_ffpe_paired_effects.tsv"
  ),
  orthogonal = file.path(
    run_root, "module_orthogonal_context",
    "module_orthogonal_context_summary.tsv"
  ),
  perturbation = file.path(
    run_root, "module_perturbation_protein",
    "module_perturbation_summary.tsv"
  ),
  sentinels = file.path(run_root, "protein_priorities.tsv"),
  routing_contract = normalizePath(Sys.getenv(
    "STATE_AWARE_MODULE_ROUTING_ADDENDUM_PATH",
    unset = file.path(
      root, "analysis", "contracts",
      "state_aware_functional_architecture_exploratory_routing_v2_2026-08-31.md"
    )
  ), mustWork = TRUE)
)
if (!all(file.exists(unlist(paths)))) {
  stop("At least one exploratory WGCNA integration input is missing")
}

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)
as_bool <- function(value) {
  if (is.logical(value)) value else tolower(as.character(value)) == "true"
}
validation <- read_tsv(paths$validation)
selected_modules <- validation$module[as_bool(validation$analysis_route_pass)]
expected_modules <- c("M02", "M03", "M04", "M05", "M06", "M09", "M10")
if (!identical(sort(selected_modules), sort(expected_modules))) {
  stop("The exploratory module route changed")
}
module_order <- c("M02", "M03", "M06", "M04", "M05", "M09", "M10")

module_palette <- c(
  M01 = "#4E79A7", M02 = "#E15759", M03 = "#59A14F",
  M04 = "#B07AA1", M05 = "#F28E2B", M06 = "#76B7B2",
  M07 = "#9C755F", M08 = "#EDC948", M09 = "#79706E",
  M10 = "#AF7AA1", M11 = "#6B747C", grey = "#D4D9DD"
)

clean_pathway <- function(value) {
  value <- sub("^(HALLMARK_|REACTOME_|GOBP_)", "", value)
  value <- gsub("_", " ", value)
  value <- tools::toTitleCase(tolower(value))
  value <- gsub("Atp", "ATP", value, fixed = TRUE)
  value <- gsub("Tca", "TCA", value, fixed = TRUE)
  value <- gsub("Mtorc1", "mTORC1", value, fixed = TRUE)
  value <- dplyr::recode(
    value,
    "mTORC1 Signaling" = "mTORC1 signaling",
    "Oxidative Phosphorylation" = "Oxidative phosphorylation",
    "Asparagine N Linked Glycosylation" = "Asparagine N-linked glycosylation",
    "Proton Motive Force Driven ATP Synthesis" =
      "Proton-motive-force-driven ATP synthesis",
    "Protein Localization" = "Protein localization",
    "Pkmts Methylate Histone Lysines" = "Histone lysine methylation"
  )
  value
}

extract_dendrogram_segments <- function(node, positions) {
  if (is.leaf(node)) {
    label <- attr(node, "label")
    return(list(
      x = unname(positions[[label]]),
      height = 0,
      segments = data.frame(
        x = numeric(), y = numeric(), xend = numeric(), yend = numeric()
      )
    ))
  }
  children <- lapply(node, extract_dendrogram_segments, positions = positions)
  height <- attr(node, "height")
  child_x <- vapply(children, `[[`, numeric(1), "x")
  child_height <- vapply(children, `[[`, numeric(1), "height")
  vertical <- data.frame(
    x = child_x, y = child_height, xend = child_x, yend = height
  )
  horizontal <- data.frame(
    x = min(child_x), y = height, xend = max(child_x), yend = height
  )
  list(
    x = mean(range(child_x)), height = height,
    segments = bind_rows(
      lapply(children, `[[`, "segments"), vertical, horizontal
    )
  )
}

audit <- readRDS(paths$audit_object)
tree <- audit$primary$dendrograms[[1L]]
block_genes <- audit$primary$blockGenes[[1L]]
network_gene_names <- names(audit$primary$colors)[block_genes]
tree$labels <- network_gene_names
dendrogram <- as.dendrogram(tree)
ordered_genes <- labels(dendrogram)
positions <- setNames(seq_along(ordered_genes), ordered_genes)
dendrogram_data <- extract_dendrogram_segments(dendrogram, positions)

modules <- read_tsv(paths$modules)
gene_module <- setNames(modules$module, modules$gene)
module_strip <- data.frame(
  gene = ordered_genes,
  x = seq_along(ordered_genes),
  module = unname(gene_module[ordered_genes]),
  stringsAsFactors = FALSE
)
module_strip$module[is.na(module_strip$module)] <- "grey"
maximum_height <- max(dendrogram_data$segments$yend)
strip_y <- -0.026 * maximum_height
p_a <- ggplot() +
  geom_segment(
    data = dendrogram_data$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.10, colour = alpha(figure_colours[["ink"]], 0.36)
  ) +
  geom_tile(
    data = module_strip,
    aes(x = x, y = strip_y, fill = module),
    width = 1.1, height = 0.024 * maximum_height
  ) +
  scale_fill_manual(
    values = module_palette,
    breaks = sprintf("M%02d", 1:11),
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(-0.055 * maximum_height, maximum_height),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = NULL, y = "Consensus dissimilarity") +
  theme_cb() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "bottom", legend.justification = "left",
    legend.margin = margin(0, 0, 0, 0, "mm"),
    panel.border = element_rect(
      fill = NA, colour = figure_colours[["line"]], linewidth = 0.3
    )
  )

rank_enrichment <- read_tsv(paths$rank_enrichment)
all_module_order <- sprintf("M%02d", 1:11)
context_order <- c(
  "discovery_common", "discovery_ABS", "discovery_GOB", "discovery_TAC",
  "heldout_common", "heldout_ABS", "heldout_GOB", "heldout_TAC"
)
context_labels <- c(
  discovery_common = "Disc.\ncommon", discovery_ABS = "Disc.\nABS",
  discovery_GOB = "Disc.\nGOB", discovery_TAC = "Disc.\nTAC",
  heldout_common = "Held\ncommon", heldout_ABS = "Held\nABS",
  heldout_GOB = "Held\nGOB", heldout_TAC = "Held\nTAC"
)
association_plot_data <- rank_enrichment %>%
  filter(module %in% all_module_order, context %in% context_order) %>%
  mutate(
    signed_evidence = direction_sign * pmin(
      -log10(pmax(fdr, 1e-12)), 12
    ),
    context = factor(context, levels = context_order, labels = context_labels),
    module = factor(module, levels = rev(all_module_order))
  )
p_b <- ggplot(
  association_plot_data, aes(context, module, fill = signed_evidence)
) +
  geom_tile(colour = "white", linewidth = 0.30) +
  scale_fill_gradient2(
    low = figure_colours[["normal"]], mid = "white",
    high = figure_colours[["adenoma"]], midpoint = 0,
    limits = c(-12, 12), breaks = c(-12, 0, 12),
    name = expression("Signed " * -log[10] * " FDR")
  ) +
  labs(x = NULL, y = NULL) +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 5.2),
    legend.position = "top", legend.justification = "left",
    panel.border = element_rect(
      fill = NA, colour = figure_colours[["line"]], linewidth = 0.3
    )
  )

community_overlap <- read_tsv(paths$community_overlap) %>%
  filter(module %in% selected_modules, fdr <= 0.05) %>%
  group_by(module) %>%
  arrange(fdr, desc(overlap), .by_group = TRUE) %>%
  slice_head(n = 2L) %>%
  ungroup()
selected_communities <- unique(community_overlap$community_id)
community_plot_data <- read_tsv(paths$community_overlap) %>%
  filter(
    module %in% selected_modules,
    community_id %in% selected_communities,
    fdr <= 0.05
  ) %>%
  mutate(
    evidence = pmin(-log10(pmax(fdr, 1e-12)), 12),
    pathway = clean_pathway(community_label),
    pathway = factor(pathway, levels = rev(unique(pathway))),
    module = factor(module, levels = module_order)
  )
p_c <- ggplot(
  community_plot_data,
  aes(module, pathway, size = overlap, fill = evidence)
) +
  geom_point(shape = 21, colour = "white", stroke = 0.3) +
  scale_x_discrete(drop = FALSE) +
  scale_size_continuous(range = c(1.2, 4.2), name = "Genes") +
  scale_fill_gradient(
    low = figure_colours[["pale_gold"]], high = figure_colours[["gold"]],
    name = expression(-log[10] * " overlap FDR")
  ) +
  labs(x = NULL, y = NULL) +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 5.4),
    axis.text.y = element_text(size = 5.0),
    legend.position = "top", legend.justification = "left"
  )

external <- read_tsv(paths$external) %>%
  transmute(
    module,
    context = cohort,
    estimate = oriented_adenoma_effect_sd,
    p_value
  )
ffpe <- read_tsv(paths$ffpe) %>%
  transmute(
    module, context = "Paired FFPE",
    estimate = median_oriented_paired_difference,
    p_value = p_paired_wilcoxon
  )
orthogonal <- read_tsv(paths$orthogonal)
orthogonal_long <- bind_rows(
  orthogonal %>% transmute(
    module, context = "Becker snRNA-seq", estimate = becker_effect_sd,
    p_value = becker_p_value
  ),
  orthogonal %>% transmute(
    module, context = "RNA–ATAC", estimate = rna_atac_patient_median_spearman_rho,
    p_value = rna_atac_patient_median_spearman_p
  ),
  orthogonal %>% transmute(
    module, context = "CRC Atlas", estimate = atlas_effect_sd,
    p_value = atlas_p_value
  ),
  orthogonal %>% transmute(
    module, context = "Spatial", estimate = spatial_median_direction_oriented_effect,
    p_value = spatial_p_paired_wilcoxon
  )
)
evidence_context_order <- c(
  "GSE8671", "GSE50114", "GSE41657", "GSE40362", "GSE72820",
  "Paired FFPE",
  "Becker snRNA-seq", "RNA–ATAC", "CRC Atlas", "Spatial"
)
evidence_plot_data <- bind_rows(external, ffpe, orthogonal_long) %>%
  filter(module %in% selected_modules) %>%
  mutate(
    context = factor(context, levels = evidence_context_order),
    module = factor(module, levels = rev(module_order)),
    direction = factor(
      ifelse(estimate >= 0, "Expected direction", "Opposite direction"),
      levels = c("Expected direction", "Opposite direction")
    ),
    evidence_strength = pmin(-log10(pmax(p_value, 1e-8)), 8)
  )
p_d <- ggplot(
  evidence_plot_data,
  aes(context, module, size = evidence_strength, fill = direction)
) +
  geom_point(shape = 21, colour = figure_colours[["ink"]], stroke = 0.25) +
  scale_fill_manual(values = c(
    "Expected direction" = figure_colours[["adenoma"]],
    "Opposite direction" = figure_colours[["normal"]]
  )) +
  scale_size_continuous(
    range = c(0.65, 4.4), breaks = c(1, 2, 4, 6, 8),
    name = expression(-log[10] * " P")
  ) +
  labs(x = NULL, y = NULL) +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 5.5, angle = 32, hjust = 1),
    legend.position = "top", legend.justification = "left",
    panel.grid.major = element_line(
      colour = figure_colours[["pale"]], linewidth = 0.25
    )
  )

comparison_labels <- c(
  pri724_reversal_of_trametinib = "PRI-724 after\ntrametinib",
  trametinib_vs_dmso = "Trametinib",
  wnt_rspo_withdrawal_in_APC_KO = "WNT/RSPO withdrawal\nAPC-KO",
  wnt_rspo_withdrawal_in_WT = "WNT/RSPO withdrawal\nWT",
  ascl2_ko_vs_resting_wt = "ASCL2 KO",
  tcf7l2_KO_vs_WT = "TCF7L2 KO",
  conditional_wnt_silencing = "Conditional WNT\nsilencing",
  apc_restoration_shApc = "APC restoration",
  apc_restoration_shApc_Kras = "APC restoration\n+ KRAS"
)
perturbation <- read_tsv(paths$perturbation) %>%
  filter(
    module %in% selected_modules,
    interpretation_role == "causal or pathway perturbation",
    comparison %in% names(comparison_labels)
  ) %>%
  mutate(
    comparison = factor(
      comparison, levels = names(comparison_labels), labels = comparison_labels
    ),
    module = factor(module, levels = rev(module_order)),
    consistent_n3 = n_units >= 3L & as_bool(all_units_positive_reversal)
  )
p_e <- ggplot(
  perturbation, aes(comparison, module, fill = mean_module_reversal)
) +
  geom_tile(colour = "white", linewidth = 0.30) +
  geom_point(
    data = perturbation %>% filter(consistent_n3),
    shape = 21, size = 1.45, stroke = 0.3,
    fill = "white", colour = figure_colours[["ink"]]
  ) +
  scale_fill_gradient2(
    low = figure_colours[["normal"]], mid = "white",
    high = figure_colours[["green"]], midpoint = 0,
    limits = c(-1.4, 1.4), oob = squish,
    name = "Module reversal"
  ) +
  labs(x = NULL, y = NULL) +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 5.2, angle = 32, hjust = 1),
    legend.position = "top", legend.justification = "left",
    panel.border = element_rect(
      fill = NA, colour = figure_colours[["line"]], linewidth = 0.3
    )
  )

write_source(
  dendrogram_data$segments, source_data_root,
  "exploratory_wgcna_v2a_dendrogram_segments.tsv"
)
write_source(
  module_strip, source_data_root,
  "exploratory_wgcna_v2a_module_strip.tsv"
)
write_source(
  association_plot_data, source_data_root,
  "exploratory_wgcna_v2b_module_association.tsv"
)
write_source(
  community_plot_data, source_data_root,
  "exploratory_wgcna_v2c_module_pathway_overlap.tsv"
)
write_source(
  evidence_plot_data, source_data_root,
  "exploratory_wgcna_v2d_independent_context_evidence.tsv"
)
write_source(
  perturbation, source_data_root,
  "exploratory_wgcna_v2e_perturbation_response.tsv"
)

perturbation_matrix <- perturbation %>%
  select(module, comparison, mean_module_reversal) %>%
  pivot_wider(names_from = comparison, values_from = mean_module_reversal)
numeric_matrix <- as.matrix(perturbation_matrix[, -1L, drop = FALSE])
rownames(numeric_matrix) <- as.character(perturbation_matrix$module)
profile_correlation <- cor(t(numeric_matrix), use = "pairwise.complete.obs")
write.table(
  data.frame(module = rownames(profile_correlation), profile_correlation),
  file.path(integration_root, "module_perturbation_profile_correlations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

integrated_summary <- validation %>%
  filter(module %in% selected_modules) %>%
  select(
    module, module_size, strict_programme_overlap, heldout_direction,
    heldout_common_fdr, minimum_preservation_z, pooled_effect, ci_low, ci_high,
    p_value, i2_percent, external_gate_pass, technical_gate_pass,
    becker_effect_sd, becker_p_value,
    rna_atac_patient_median_spearman_rho,
    rna_atac_patient_median_spearman_p,
    atlas_effect_sd, atlas_p_value,
    spatial_median_direction_oriented_effect,
    spatial_positive_section_fraction, spatial_p_paired_wilcoxon,
    n_independent_datasets_positive_reversal,
    one_model_n3_all_positive, perturbation_gate_pass,
    n_sentinel_proteins, n_regulatory_nodes
  ) %>%
  arrange(match(module, module_order))
write.table(
  integrated_summary,
  file.path(integration_root, "module_integrated_evidence.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

figure <- p_a / (p_b | p_c) / p_d / p_e +
  plot_layout(heights = c(0.78, 1.12, 0.82, 0.90)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)
export_cb_figure(
  figure, figure_root, figure_stem,
  width_mm = 178, height_mm = 220
)

manifest <- list(
  analysis = "build_exploratory_wgcna_integration_v2",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = 20260831,
  status = "post-result exploratory presentation; not confirmatory validation",
  selected_modules = selected_modules,
  input_sha256 = as.list(vapply(paths, sha256, character(1))),
  output_sha256 = list(
    integrated_summary = sha256(file.path(
      integration_root, "module_integrated_evidence.tsv"
    )),
    perturbation_correlations = sha256(file.path(
      integration_root, "module_perturbation_profile_correlations.tsv"
    )),
    figure_pdf = sha256(file.path(
      figure_root, paste0(figure_stem, ".pdf")
    )),
    figure_svg = sha256(file.path(
      figure_root, paste0(figure_stem, ".svg")
    ))
  ),
  package_versions = list(
    R = as.character(getRversion()),
    WGCNA = as.character(packageVersion("WGCNA")),
    ggplot2 = as.character(packageVersion("ggplot2")),
    dplyr = as.character(packageVersion("dplyr")),
    patchwork = as.character(packageVersion("patchwork"))
  )
)
write_json(
  manifest,
  file.path(integration_root, "exploratory_wgcna_integration_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)
message("Exploratory WGCNA integration figure and source data exported")
