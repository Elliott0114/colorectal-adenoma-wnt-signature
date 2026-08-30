#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")
suppressPackageStartupMessages({
  library(igraph)
  library(msigdbr)
})

set.seed(20260830)
root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(
  root, "results", "state_aware_program_v1", "functional_architecture_v1"
)
pathway_root <- file.path(result_root, "pathway_replication")
wgcna_root <- file.path(result_root, "consensus_wgcna")
out_dir <- file.path(root, "figures", "communications_biology_v3.0")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

pathway_long <- read_tsv(file.path(result_root, "pathway_replication.tsv"))
pathway_summary <- read_tsv(file.path(
  pathway_root, "pathway_replication_summary.tsv"
))
representatives <- read_tsv(file.path(
  pathway_root, "replicated_pathway_representatives.tsv"
))
edges <- read_tsv(file.path(
  pathway_root, "replicated_pathway_network_edges.tsv"
))
communities <- read_tsv(file.path(
  pathway_root, "replicated_pathway_communities.tsv"
))
ranking <- read_tsv(file.path(
  root, "results", "state_aware_program_v1", "common_effects",
  "cross_state_common_effects.tsv.gz"
))
heldout <- read_tsv(file.path(
  root, "results", "state_aware_program_v1", "heldout_validation",
  "heldout_cross_state_common_effects.tsv.gz"
))

clean_pathway <- function(value) {
  value <- sub("^(HALLMARK_|REACTOME_|GOBP_)", "", value)
  value <- gsub("_", " ", value)
  value <- tools::toTitleCase(tolower(value))
  value <- gsub("Atp", "ATP", value, fixed = TRUE)
  value <- gsub("Tca", "TCA", value, fixed = TRUE)
  value <- gsub("Nadh", "NADH", value, fixed = TRUE)
  value <- gsub("Jak Stat", "JAK–STAT", value, fixed = TRUE)
  value
}

selected_representatives <- pathway_summary %>%
  filter(
    replicated,
    collection %in% c("Hallmark 2026.1.Hs", "Reactome 2026.1.Hs")
  ) %>%
  mutate(replication_strength = pmin(
    -log10(pmax(discovery_common_fdr, .Machine$double.xmin)),
    -log10(pmax(heldout_common_fdr, .Machine$double.xmin))
  )) %>%
  group_by(collection) %>%
  arrange(desc(replication_strength), gene_set, .by_group = TRUE) %>%
  slice_head(n = 6L) %>%
  ungroup() %>%
  arrange(collection, desc(replication_strength), gene_set) %>%
  mutate(
    base_label = clean_pathway(gene_set),
    collection_short = recode(
      collection,
      "Hallmark 2026.1.Hs" = "Hallmark",
      "Reactome 2026.1.Hs" = "Reactome",
      "GO:BP 2026.1.Hs" = "GO BP"
    )
  ) %>%
  group_by(base_label) %>%
  mutate(
    display_label = if (n() > 1L) {
      paste0(base_label, " (", collection_short, ")")
    } else {
      base_label
    }
  ) %>%
  ungroup()
selected_key <- paste(
  selected_representatives$collection,
  selected_representatives$gene_set,
  sep = "__"
)

context_labels <- c(
  discovery_common = "Disc.\ncommon",
  discovery_ABS = "Disc.\nABS",
  discovery_GOB = "Disc.\nGOB",
  discovery_TAC = "Disc.\nTAC",
  heldout_common = "Held-out\ncommon",
  heldout_ABS = "Held-out\nABS",
  heldout_GOB = "Held-out\nGOB",
  heldout_TAC = "Held-out\nTAC"
)
heatmap_data <- pathway_long %>%
  mutate(key = paste(collection, gene_set, sep = "__")) %>%
  filter(key %in% selected_key) %>%
  mutate(
    signed_evidence = direction_sign * pmin(-log10(pmax(fdr, 1e-15)), 15),
    display_label = selected_representatives$display_label[
      match(key, selected_key)
    ],
    display_label = factor(
      display_label,
      levels = rev(selected_representatives$display_label)
    ),
    partition = factor(
      ifelse(grepl("^discovery", context), "Discovery", "Held-out"),
      levels = c("Discovery", "Held-out")
    ),
    state_axis = factor(
      sub("^(discovery|heldout)_", "", context),
      levels = c("common", "ABS", "GOB", "TAC"),
      labels = c("Common", "ABS", "GOB", "TAC")
    ),
    context_label = factor(
      unname(context_labels[context]),
      levels = unname(context_labels)
    )
  )

p4a <- ggplot(heatmap_data, aes(state_axis, display_label, fill = signed_evidence)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_point(
    data = heatmap_data %>% filter(fdr <= 0.10),
    shape = 21, size = 0.65, stroke = 0.2, fill = "white", colour = "white"
  ) +
  scale_fill_gradient2(
    low = figure_colours[["normal"]], mid = "white",
    high = figure_colours[["adenoma"]], midpoint = 0,
    limits = c(-15, 15), breaks = c(-15, -5, 0, 5, 15),
    name = expression("Signed " * -log[10] * " FDR")
  ) +
  facet_grid(cols = vars(partition), scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL) +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 5.2),
    axis.text.y = element_text(size = 5.3),
    strip.text.x = element_text(size = 5.4, face = "bold"),
    strip.background = element_rect(fill = "white", colour = NA),
    panel.spacing.x = grid::unit(1.2, "mm"),
    legend.position = "top",
    legend.justification = "left",
    panel.border = element_rect(fill = NA, colour = figure_colours[["line"]], linewidth = 0.3)
  )

hallmark_frame <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- lapply(
  split(hallmark_frame$gene_symbol, hallmark_frame$gs_name),
  unique
)
curve_candidates <- pathway_summary %>%
  filter(replicated, collection == "Hallmark 2026.1.Hs") %>%
  arrange(discovery_common_fdr) %>%
  slice_head(n = 3L)
ranked <- sort(setNames(ranking$common_z, ranking$gene), decreasing = TRUE)

running_enrichment <- function(statistic, genes, pathway_name) {
  hit <- names(statistic) %in% genes
  n_hit <- sum(hit)
  if (!n_hit || n_hit == length(hit)) {
    return(data.frame())
  }
  hit_weight <- abs(statistic)
  hit_weight[!hit] <- 0
  hit_weight <- hit_weight / sum(hit_weight)
  miss_weight <- ifelse(hit, 0, 1 / sum(!hit))
  score <- cumsum(hit_weight - miss_weight)
  data.frame(
    rank = seq_along(score), enrichment_score = score,
    hit = hit, pathway = clean_pathway(pathway_name),
    stringsAsFactors = FALSE
  )
}
curve_data <- do.call(rbind, lapply(curve_candidates$gene_set, function(pathway) {
  running_enrichment(ranked, hallmark_sets[[pathway]], pathway)
}))
curve_hits <- curve_data %>%
  filter(hit) %>%
  group_by(pathway) %>%
  group_modify(function(frame, key) {
    index <- unique(round(seq(
      1, nrow(frame), length.out = min(80L, nrow(frame))
    )))
    frame[index, , drop = FALSE]
  }) %>%
  ungroup()

p4b <- ggplot(curve_data, aes(rank, enrichment_score, colour = pathway)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = figure_colours[["line"]]) +
  geom_line(linewidth = 0.7) +
  geom_rug(
    data = curve_hits, aes(x = rank), sides = "b", inherit.aes = FALSE,
    linewidth = 0.18, alpha = 0.32, colour = figure_colours[["muted"]]
  ) +
  scale_colour_manual(values = c(
    figure_colours[["normal"]], figure_colours[["green"]],
    figure_colours[["purple"]]
  )) +
  labs(x = "8,221-gene common-effect rank", y = "Running enrichment score") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = NULL)) +
  theme_cb() +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.text = element_text(size = 5.2),
    legend.key.width = grid::unit(6, "mm")
  )

node_data <- communities %>%
  left_join(
    pathway_summary %>% select(
      collection, gene_set, discovery_common_direction,
      discovery_common_fdr, heldout_common_fdr
    ),
    by = c("collection", "gene_set")
  ) %>%
  mutate(
    community_size = ave(node_id, community_id, FUN = length),
    direction = discovery_common_direction,
    evidence = pmin(
      -log10(pmax(discovery_common_fdr, .Machine$double.xmin)),
      -log10(pmax(heldout_common_fdr, .Machine$double.xmin))
    )
  )
graph <- graph_from_data_frame(
  edges, directed = FALSE,
  vertices = data.frame(name = node_data$node_id, stringsAsFactors = FALSE)
)
layout <- layout_with_fr(graph, weights = if (ecount(graph)) E(graph)$jaccard else NULL,
                         niter = 1500, grid = "nogrid")
positions <- data.frame(
  node_id = V(graph)$name, x = layout[, 1], y = layout[, 2],
  stringsAsFactors = FALSE
) %>% left_join(node_data, by = "node_id")
edge_plot <- if (nrow(edges)) {
  edges %>%
    left_join(positions %>% select(node_id, x, y), by = c("from" = "node_id")) %>%
    rename(x_from = x, y_from = y) %>%
    left_join(positions %>% select(node_id, x, y), by = c("to" = "node_id")) %>%
    rename(x_to = x, y_to = y)
} else {
  data.frame()
}
label_communities <- node_data %>%
  filter(
    direction == "Down",
    collection == "Hallmark 2026.1.Hs"
  ) %>%
  group_by(community_id) %>%
  arrange(discovery_common_fdr, gene_set, .by_group = TRUE) %>%
  slice_head(n = 1L) %>%
  ungroup() %>%
  transmute(
    community_id,
    community_size,
    best_evidence = evidence,
    direction,
    label_gene_set = gene_set
  ) %>%
  arrange(desc(community_size), desc(best_evidence)) %>%
  slice_head(n = 6L)
label_positions <- positions %>%
  filter(community_id %in% label_communities$community_id) %>%
  group_by(community_id) %>%
  summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
  left_join(label_communities, by = "community_id") %>%
  mutate(label = vapply(
    clean_pathway(label_gene_set),
    function(value) paste(strwrap(value, width = 27), collapse = "\n"),
    character(1)
  ))

p4c <- ggplot() +
  {if (nrow(edge_plot)) geom_segment(
    data = edge_plot,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linewidth = jaccard),
    colour = figure_colours[["line"]], alpha = 0.65
  )} +
  geom_point(
    data = positions,
    aes(x, y, fill = direction, size = evidence, shape = collection),
    colour = "white", stroke = 0.35, alpha = 0.95
  ) +
  ggrepel::geom_label_repel(
    data = label_positions,
    aes(x, y, label = label),
    family = figure_font, size = 1.75, lineheight = 0.9,
    fill = alpha("white", 0.88), colour = figure_colours[["ink"]],
    label.size = 0, label.padding = grid::unit(0.6, "mm"),
    box.padding = grid::unit(0.8, "mm"), point.padding = 0,
    min.segment.length = 0, segment.size = 0.25,
    segment.colour = figure_colours[["line"]], max.overlaps = Inf,
    seed = 20260830
  ) +
  scale_fill_manual(values = c(Down = figure_colours[["normal"]], Up = figure_colours[["adenoma"]])) +
  scale_shape_manual(values = c(
    "Hallmark 2026.1.Hs" = 21,
    "Reactome 2026.1.Hs" = 22,
    "GO:BP 2026.1.Hs" = 24
  )) +
  scale_size_continuous(range = c(1.5, 4.2), guide = "none") +
  scale_linewidth_continuous(range = c(0.25, 0.85), guide = "none") +
  guides(
    fill = guide_legend(title = "Direction", nrow = 1, byrow = TRUE),
    shape = "none"
  ) +
  coord_equal(clip = "off") +
  labs(fill = "Direction", shape = "Collection") +
  theme_void(base_family = figure_font) +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.title = element_text(size = 5.2),
    legend.text = element_text(size = 5.2),
    legend.key.width = grid::unit(4.5, "mm"),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

leading <- read_tsv(file.path(
  pathway_root, "replicated_pathway_leading_edges.tsv"
)) %>%
  filter(gene_set %in% curve_candidates$gene_set)
leading_genes <- unique(unlist(strsplit(
  leading$leading_edge[!is.na(leading$leading_edge)], ";", fixed = TRUE
)))
gene_evidence <- ranking %>%
  filter(gene %in% leading_genes) %>%
  arrange(desc(abs(common_z))) %>%
  slice_head(n = 24L)
gene_context <- bind_rows(
  data.frame(
    gene = ranking$gene, partition = "Discovery",
    ABS = ranking$logFC_ABS / ranking$raw_se_ABS,
    GOB = ranking$logFC_GOB / ranking$raw_se_GOB,
    TAC = ranking$logFC_TAC / ranking$raw_se_TAC
  ),
  data.frame(
    gene = heldout$gene, partition = "Held-out",
    ABS = heldout$logFC_ABS / heldout$raw_se_ABS,
    GOB = heldout$logFC_GOB / heldout$raw_se_GOB,
    TAC = heldout$logFC_TAC / heldout$raw_se_TAC
  )
) %>%
  filter(gene %in% gene_evidence$gene) %>%
  tidyr::pivot_longer(c(ABS, GOB, TAC), names_to = "state", values_to = "z") %>%
  mutate(
    context = factor(
      paste(partition, state, sep = "\n"),
      levels = c(
        "Discovery\nABS", "Discovery\nGOB", "Discovery\nTAC",
        "Held-out\nABS", "Held-out\nGOB", "Held-out\nTAC"
      )
    ),
    gene = factor(gene, levels = rev(gene_evidence$gene))
  )

p4d <- ggplot(gene_context, aes(context, gene, fill = pmax(pmin(z, 8), -8))) +
  geom_tile(colour = "white", linewidth = 0.28) +
  scale_y_discrete(expand = expansion(add = c(0.15, 1.0))) +
  scale_fill_gradient2(
    low = figure_colours[["normal"]], mid = "white",
    high = figure_colours[["adenoma"]], midpoint = 0,
    limits = c(-8, 8), breaks = c(-8, 0, 8), name = "Gene z"
  ) +
  labs(x = NULL, y = NULL) +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 5.4), axis.text.y = element_text(size = 5.1),
    legend.position = "top", legend.justification = "left",
    panel.border = element_rect(fill = NA, colour = figure_colours[["line"]], linewidth = 0.3)
  )

module_summary_path <- file.path(wgcna_root, "module_internal_gate_summary.tsv")
has_modules <- file.exists(module_summary_path) &&
  any(read_tsv(module_summary_path)$internal_gate_pass)

write_source(heatmap_data, source_dir, "figure4a_pathway_replication_heatmap.tsv")
write_source(curve_data, source_dir, "figure4b_running_enrichment.tsv")
write_source(positions, source_dir, "figure4c_pathway_network_nodes.tsv")
write_source(edge_plot, source_dir, "figure4c_pathway_network_edges.tsv")
write_source(gene_context, source_dir, "figure4d_leading_edge_gene_heatmap.tsv")

if (has_modules) {
  module_summary <- read_tsv(module_summary_path) %>% filter(internal_gate_pass)
  preservation <- read_tsv(file.path(wgcna_root, "module_preservation.tsv")) %>%
    filter(module %in% module_summary$module) %>%
    mutate(
      state = factor(state, levels = c("ABS", "GOB", "TAC")),
      module = factor(module, levels = rev(module_summary$module))
    )
  p4e <- ggplot(preservation, aes(state, module, fill = pmin(zsummary, 10))) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(
      aes(label = sprintf("%.1f", zsummary), colour = zsummary >= 6),
      family = figure_font, size = 1.65
    ) +
    scale_fill_gradient(
      low = figure_colours[["pale_blue"]], high = figure_colours[["normal"]],
      limits = c(0, 10), name = expression(Z[summary])
    ) +
    scale_colour_manual(
      values = c(`TRUE` = "white", `FALSE` = figure_colours[["ink"]]),
      guide = "none"
    ) +
    labs(x = NULL, y = NULL) +
    theme_cb() +
    theme(legend.position = "top", legend.justification = "left")

  overlap <- read_tsv(file.path(
    wgcna_root, "module_pathway_community_overlap.tsv"
  )) %>%
    filter(module %in% module_summary$module, fdr <= 0.10) %>%
    group_by(module, community_id, community_label) %>%
    slice_min(fdr, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(fdr) %>%
    group_by(module) %>%
    slice_head(n = 6L) %>%
    ungroup() %>%
    mutate(
      evidence = pmin(-log10(pmax(fdr, 1e-12)), 12),
      pathway = clean_pathway(community_label),
      pathway = factor(pathway, levels = rev(unique(pathway))),
      module = factor(module, levels = module_summary$module)
    )
  p4f <- ggplot(overlap, aes(module, pathway, size = overlap, fill = evidence)) +
    geom_point(shape = 21, colour = "white", stroke = 0.3) +
    scale_size_continuous(range = c(1.4, 4.8), name = "Genes") +
    scale_fill_gradient(
      low = figure_colours[["pale_gold"]], high = figure_colours[["gold"]],
      name = expression(-log[10] * " overlap FDR")
    ) +
    labs(x = NULL, y = NULL) +
    theme_cb() +
    theme(
      axis.text.y = element_text(size = 5.1),
      legend.position = "top", legend.justification = "left"
    )

  write_source(preservation, source_dir, "figure4e_module_preservation.tsv")
  write_source(overlap, source_dir, "figure4f_module_function_overlap.tsv")
  figure <- (p4a | p4b) / (p4c | p4d) / (p4e | p4f) +
    plot_layout(heights = c(1.05, 1, 0.82)) +
    plot_annotation(tag_levels = "a")
  height_mm <- 205
} else {
  figure <- (p4a | p4b) / (p4c | p4d) +
    plot_layout(heights = c(1.05, 1)) +
    plot_annotation(tag_levels = "a")
  height_mm <- 160
}

figure <- tagged(figure)
export_cb_figure(
  figure, out_dir, "figure4_functional_architecture",
  width_mm = 178, height_mm = height_mm
)
message("Functional-architecture Figure 4 exported to ", out_dir)
