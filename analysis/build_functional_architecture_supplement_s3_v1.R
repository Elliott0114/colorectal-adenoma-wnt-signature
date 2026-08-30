#!/usr/bin/env Rscript

source("analysis/state_shared_revision_figure_utils_v3.R")

set.seed(20260830)
root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(
  root, "results", "state_aware_program_v1", "functional_architecture_v1"
)
pathway_root <- file.path(result_root, "pathway_replication")
out_dir <- file.path(root, "figures", "communications_biology_v3.0")
source_dir <- file.path(out_dir, "source_data")

pathway_long <- read_tsv(file.path(result_root, "pathway_replication.tsv"))
pathway_summary <- read_tsv(file.path(
  pathway_root, "pathway_replication_summary.tsv"
))
ranking <- read_tsv(file.path(
  root, "results", "state_aware_program_v1", "common_effects",
  "cross_state_common_effects.tsv.gz"
))
heldout <- read_tsv(file.path(
  root, "results", "state_aware_program_v1", "heldout_validation",
  "heldout_cross_state_common_effects.tsv.gz"
))

collection_labels <- c(
  "Hallmark 2026.1.Hs" = "Hallmark",
  "Reactome 2026.1.Hs" = "Reactome",
  "GO:BP 2026.1.Hs" = "GO Biological Process"
)
collection_palette <- c(
  Hallmark = figure_colours[["adenoma"]],
  Reactome = figure_colours[["normal"]],
  `GO Biological Process` = figure_colours[["purple"]]
)

clean_pathway <- function(value) {
  value <- sub("^(HALLMARK_|REACTOME_|GOBP_)", "", value)
  value <- gsub("_", " ", value)
  value <- tools::toTitleCase(tolower(value))
  value <- gsub("Atp", "ATP", value, fixed = TRUE)
  value <- gsub("Mtorc1", "mTORC1", value, fixed = TRUE)
  value
}

# a. Prespecified pathway attrition by collection.
stage_counts <- pathway_summary %>%
  mutate(collection_label = unname(collection_labels[collection])) %>%
  group_by(collection_label) %>%
  summarise(
    Tested = n(),
    `Discovery FDR ≤ 0.05` = sum(
      is.finite(discovery_common_fdr) & discovery_common_fdr <= 0.05
    ),
    `Replicated in held-out` = sum(replicated, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    -collection_label, names_to = "stage", values_to = "n_pathways"
  ) %>%
  mutate(
    stage = factor(
      stage,
      levels = c("Tested", "Discovery FDR ≤ 0.05", "Replicated in held-out")
    )
  )

p3a <- ggplot(
  stage_counts,
  aes(stage, n_pathways, colour = collection_label, group = collection_label)
) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.1) +
  geom_text(
    aes(label = scales::comma(n_pathways)),
    nudge_y = 0.08, family = figure_font, size = 1.65,
    show.legend = FALSE
  ) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 3000),
    labels = scales::comma
  ) +
  scale_colour_manual(values = collection_palette) +
  labs(x = NULL, y = "Pathways (log scale)") +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 5.2),
    legend.position = "top", legend.justification = "left"
  )

# b. Discovery versus held-out evidence for all replicated pathways.
replicated <- pathway_summary %>%
  filter(replicated) %>%
  mutate(
    collection_label = unname(collection_labels[collection]),
    discovery_signed_evidence = ifelse(
      discovery_common_direction == "Up", 1, -1
    ) * pmin(-log10(pmax(discovery_common_fdr, 1e-20)), 20),
    heldout_signed_evidence = ifelse(
      heldout_common_direction == "Up", 1, -1
    ) * pmin(-log10(pmax(heldout_common_fdr, 1e-20)), 20),
    minimum_evidence = pmin(
      abs(discovery_signed_evidence), abs(heldout_signed_evidence)
    )
  )
label_replicated <- replicated %>%
  filter(collection != "GO:BP 2026.1.Hs") %>%
  arrange(desc(minimum_evidence), collection, gene_set) %>%
  slice_head(n = 7L) %>%
  mutate(label = clean_pathway(gene_set))

p3b <- ggplot(
  replicated,
  aes(discovery_signed_evidence, heldout_signed_evidence,
      colour = collection_label, shape = discovery_common_direction)
) +
  geom_abline(
    slope = 1, intercept = 0, colour = figure_colours[["line"]],
    linewidth = 0.35
  ) +
  geom_hline(yintercept = 0, colour = figure_colours[["line"]], linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = figure_colours[["line"]], linewidth = 0.3) +
  geom_point(size = 1.8, alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = label_replicated, aes(label = label),
    family = figure_font, size = 1.55, colour = figure_colours[["ink"]],
    box.padding = 0.25, point.padding = 0.15,
    min.segment.length = 0, segment.size = 0.25,
    segment.colour = figure_colours[["line"]], max.overlaps = Inf,
    seed = 20260830, show.legend = FALSE
  ) +
  scale_colour_manual(values = collection_palette) +
  scale_shape_manual(values = c(Down = 16, Up = 17)) +
  coord_equal() +
  labs(
    x = "Discovery signed −log10 FDR",
    y = "Held-out signed −log10 FDR"
  ) +
  theme_cb() +
  guides(
    shape = "none",
    colour = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.text = element_text(size = 5.0),
    legend.key.width = grid::unit(4.2, "mm")
  )

# c. Directionally consistent but non-replicated Hallmark WNT result.
context_order <- c(
  "discovery_common", "heldout_common",
  "discovery_ABS", "heldout_ABS",
  "discovery_GOB", "heldout_GOB",
  "discovery_TAC", "heldout_TAC"
)
context_labels <- c(
  discovery_common = "Disc.\ncommon", heldout_common = "Held-out\ncommon",
  discovery_ABS = "Disc.\nABS", heldout_ABS = "Held-out\nABS",
  discovery_GOB = "Disc.\nGOB", heldout_GOB = "Held-out\nGOB",
  discovery_TAC = "Disc.\nTAC", heldout_TAC = "Held-out\nTAC"
)
wnt <- pathway_long %>%
  filter(gene_set == "HALLMARK_WNT_BETA_CATENIN_SIGNALING") %>%
  mutate(
    context = factor(context, levels = context_order),
    context_label = factor(
      unname(context_labels[as.character(context)]),
      levels = unname(context_labels[context_order])
    ),
    evidence = -log10(pmax(fdr, 1e-15))
  )

p3c <- ggplot(wnt, aes(context_label, evidence)) +
  geom_hline(
    yintercept = -log10(0.10), linetype = 3, linewidth = 0.4,
    colour = figure_colours[["normal"]]
  ) +
  geom_hline(
    yintercept = -log10(0.05), linetype = 2, linewidth = 0.4,
    colour = figure_colours[["adenoma"]]
  ) +
  geom_segment(
    aes(xend = context_label, y = 0, yend = evidence),
    linewidth = 0.55, colour = figure_colours[["line"]]
  ) +
  geom_point(
    aes(fill = direction), shape = 21, size = 2.3,
    colour = "white", stroke = 0.35
  ) +
  geom_text(
    aes(label = sprintf("%.3f", fdr)),
    nudge_y = 0.08, family = figure_font, size = 1.45,
    colour = figure_colours[["muted"]]
  ) +
  scale_fill_manual(values = c(Up = figure_colours[["adenoma"]])) +
  scale_y_continuous(limits = c(0, 1.65), breaks = c(0, 0.5, 1, 1.5)) +
  labs(x = NULL, y = "Hallmark WNT/β-catenin\n−log10 FDR") +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 4.8),
    legend.position = "none"
  )

# d. Gene-level anchors across discovery and held-out epithelial states.
anchors <- c(
  "ASCL2", "AXIN2", "RNF43", "ZNRF3", "EPHB2", "OLFM4",
  "CA2", "FABP1", "COX6C", "ACAA2"
)
anchor_context <- bind_rows(
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
  filter(gene %in% anchors) %>%
  tidyr::pivot_longer(c(ABS, GOB, TAC), names_to = "state", values_to = "z") %>%
  mutate(
    context = factor(
      paste(partition, state, sep = "\n"),
      levels = c(
        "Discovery\nABS", "Discovery\nGOB", "Discovery\nTAC",
        "Held-out\nABS", "Held-out\nGOB", "Held-out\nTAC"
      )
    ),
    gene = factor(gene, levels = rev(anchors))
  )

p3d <- ggplot(anchor_context, aes(context, gene, fill = pmax(pmin(z, 6), -6))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.1f", z), colour = abs(z) >= 3.5),
    family = figure_font, size = 1.45
  ) +
  scale_fill_gradient2(
    low = figure_colours[["normal"]], mid = "white",
    high = figure_colours[["adenoma"]], midpoint = 0,
    limits = c(-6, 6), breaks = c(-6, 0, 6), name = "Gene z"
  ) +
  scale_colour_manual(
    values = c(`TRUE` = "white", `FALSE` = figure_colours[["ink"]]),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_cb() +
  theme(
    axis.text.x = element_text(size = 5.0),
    axis.text.y = element_text(size = 5.3),
    legend.position = "top", legend.justification = "left",
    panel.border = element_rect(
      fill = NA, colour = figure_colours[["line"]], linewidth = 0.3
    )
  )

figure <- (p3a | p3b) / (p3c | p3d) +
  plot_layout(heights = c(1, 0.88)) +
  plot_annotation(tag_levels = "a")
figure <- tagged(figure)

write_source(stage_counts, source_dir, "figureS3a_pathway_attrition.tsv")
write_source(replicated, source_dir, "figureS3b_replicated_pathways.tsv")
write_source(wnt, source_dir, "figureS3c_wnt_contexts.tsv")
write_source(anchor_context, source_dir, "figureS3d_anchor_gene_replication.tsv")

export_cb_figure(
  figure, out_dir, "figureS3_functional_and_regulatory_structure",
  width_mm = 178, height_mm = 155
)
message("Supplementary Figure 3 exported to ", out_dir)
