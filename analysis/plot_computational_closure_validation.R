#!/usr/bin/env Rscript

# Submission-grade computational-closure figure. No narrative title/subtitle is
# drawn inside the graphic; the complete description belongs in the legend.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  file.path(getwd(), "analysis", "plot_computational_closure_validation.R")

ROOT <- normalizePath(file.path(dirname(script_path), ".."))
DATA_DIR <- file.path(ROOT, "results", "computational_closure_validation")
OUT_DIR <- Sys.getenv(
  "JTM_FIGURE_DIR",
  unset = file.path(ROOT, "figures", "computational_closure_validation")
)
FIGURE_STEM <- Sys.getenv(
  "JTM_FIGURE_STEM",
  unset = "figure_computational_closure_validation"
)
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

read_result <- function(filename) {
  read.delim(
    file.path(DATA_DIR, filename), check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN")
  )
}

write_source <- function(frame, filename) {
  write.table(
    frame, file.path(SOURCE_DIR, filename), sep = "\t", row.names = FALSE,
    quote = FALSE, na = ""
  )
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
  support = "#2F8F83",
  uncertain = "#8C7A6B"
)

theme_submission <- function() {
  theme_classic(base_size = 7, base_family = JTM_FONT) +
    theme(
      axis.title = element_text(size = 7, colour = COL[["ink"]]),
      axis.text = element_text(size = 6.1, colour = COL[["ink"]]),
      axis.line = element_line(linewidth = 0.35, colour = COL[["ink"]]),
      axis.ticks = element_line(linewidth = 0.3, colour = COL[["ink"]]),
      axis.ticks.length = grid::unit(1.1, "mm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 5.8),
      legend.key.height = grid::unit(2.5, "mm"),
      legend.key.width = grid::unit(3.0, "mm"),
      panel.spacing = grid::unit(1.6, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.4, face = "bold"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(1.6, 2.0, 1.6, 2.0, unit = "mm"),
      plot.tag = element_text(size = 9, face = "bold", colour = COL[["ink"]])
    )
}

export_figure <- function(plot, stem, width_mm = 170, height_mm = 155) {
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
    format = names(paths), file = unname(paths), width_mm = width_mm,
    height_mm = height_mm, resolution_dpi = c(NA, NA, 600, 300)
  )
}

unit_effects <- read_result("perturbation_unit_effects.tsv")
effect_summary <- read_result("perturbation_effect_summary.tsv")
matched_tests <- read_result("expression_matched_signature_tests.tsv")
evidence <- read_result("evidence_closure_matrix.tsv")
ulm_summary <- read_result("collectri_ulm_effect_summary.tsv")
signed_tcf <- read_result("collectri_signed_mean_tcf7l2_clone_effects.tsv")
virtual_panel <- read_result("virtual_tf_knockout_panel.tsv")
spatial_units <- read_result("spatial_locked_route_section_effects.tsv")

stopifnot(nrow(filter(evidence, layer == "empirical_perturbation")) == 11)
stopifnot(sum(filter(spatial_units,
                     comparison == "tumor_vs_non_neoplastic_epithelium",
                     feature == "route_score")$difference > 0) == 6)

# a: evidence architecture.
workflow <- data.frame(
  x = 1:5,
  label = c(
    "Locked\n100-gene route",
    "APC loss /\nApc restoration",
    "WNT–TCF /\nASCL2 perturbation",
    "Regulon + graph\nvirtual deletion",
    "Full-route\nspatial test"
  ),
  fill = c(COL[["neutral_pale"]], COL[["route_light"]], "#F2DEC1",
           "#D7E9E5", COL[["wnt_light"]])
)
workflow_arrows <- data.frame(x = 1:4, xend = 2:5, y = 1, yend = 1)
write_source(workflow, "panel_a_evidence_architecture.tsv")

p_a <- ggplot(workflow, aes(x, 1)) +
  geom_segment(
    data = workflow_arrows,
    aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
    linewidth = 0.42, colour = COL[["neutral"]],
    arrow = arrow(length = grid::unit(1.4, "mm"), type = "closed")
  ) +
  geom_label(
    aes(label = label, fill = fill), size = 2.0, lineheight = 0.92,
    linewidth = 0.25, label.padding = grid::unit(1.15, "mm"), colour = COL[["ink"]]
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.55, 5.45), ylim = c(0.72, 1.28), clip = "off") +
  theme_void(base_size = 7, base_family = JTM_FONT) +
  theme(plot.margin = margin(1.5, 5.0, 0.5, 5.0, unit = "mm"),
        plot.tag = element_text(size = 9, face = "bold"),
        plot.tag.position = c(0.01, 0.92))

# b: all empirical route effects, aligned to each frozen expected direction.
comparison_labels <- c(
  "GSE125472\nAPC-KO (+WNT)" = "GSE125472|APC_vs_WT_with_Wnt",
  "GSE125472\nAPC-KO (−WNT)" = "GSE125472|APC_vs_WT_without_Wnt",
  "GSE171910\nWNT off" = "GSE171910|conditional_wnt_silencing",
  "GSE130822\nASCL2-KO" = "GSE130822|ascl2_ko_vs_resting_wt",
  "GSE135328 HT29\nTCF7L2-KO" = "GSE135328_HT29|tcf7l2_ko_vs_wt",
  "GSE135328 HCT116\nTCF7L2-KO" = "GSE135328_HCT116|tcf7l2_ko_vs_wt",
  "GSE114059\nTrametinib" = "GSE114059|trametinib_vs_dmso",
  "GSE114059\nPRI-724 reversal" = "GSE114059|pri724_reversal_of_trametinib",
  "GSE67186\nApc restored" = "GSE67186|apc_restoration_shApc",
  "GSE67186\nApc/Kras restored" = "GSE67186|apc_restoration_shApc_Kras"
)
label_lookup <- setNames(names(comparison_labels), unname(comparison_labels))
ordered_keys <- rev(unname(comparison_labels))

forest_summary <- effect_summary %>%
  filter(feature == "route_score", expected_direction != 0) %>%
  mutate(
    key = paste(dataset, comparison, sep = "|"),
    aligned_mean = mean_difference * expected_direction,
    aligned_low = pmin(bootstrap_mean_ci_low * expected_direction,
                       bootstrap_mean_ci_high * expected_direction),
    aligned_high = pmax(bootstrap_mean_ci_low * expected_direction,
                        bootstrap_mean_ci_high * expected_direction)
  ) %>%
  filter(key %in% unname(comparison_labels)) %>%
  left_join(
    evidence %>% filter(layer == "empirical_perturbation") %>%
      mutate(key = paste(dataset, comparison, sep = "|")) %>%
      select(key, status, specificity_p), by = "key"
  ) %>%
  mutate(
    display = factor(key, levels = ordered_keys, labels = label_lookup[ordered_keys]),
    coverage = ifelse(route_reportable, "≥80% route coverage", "<80% route coverage")
  )

forest_units <- unit_effects %>%
  filter(feature == "route_score", expected_direction != 0) %>%
  mutate(key = paste(dataset, comparison, sep = "|"),
         aligned = difference * expected_direction) %>%
  filter(key %in% unname(comparison_labels)) %>%
  mutate(display = factor(key, levels = ordered_keys, labels = label_lookup[ordered_keys]))
write_source(forest_summary, "panel_b_route_effect_summary.tsv")
write_source(forest_units, "panel_b_route_effects_by_context.tsv")

status_colours <- c(
  supportive_specific = COL[["route"]],
  supportive_direction_only = COL[["wnt"]],
  discordant = COL[["uncertain"]],
  exploratory_low_coverage = COL[["neutral"]]
)
p_b <- ggplot(forest_summary, aes(y = display)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = COL[["neutral"]]) +
  geom_errorbar(
    aes(xmin = aligned_low, xmax = aligned_high), orientation = "y",
    width = 0.20, linewidth = 0.5, colour = "black"
  ) +
  geom_point(
    data = forest_units, aes(x = aligned, y = display), inherit.aes = FALSE,
    position = position_jitter(height = 0.10, width = 0, seed = 20260808),
    size = 1.1, colour = COL[["neutral"]], alpha = 0.75
  ) +
  geom_point(aes(x = aligned_mean, colour = status, shape = coverage),
             size = 2.2, stroke = 0.65, fill = "white") +
  scale_colour_manual(
    values = status_colours,
    labels = c(
      supportive_specific = "Expression-matched support",
      supportive_direction_only = "Direction only",
      discordant = "Discordant",
      exploratory_low_coverage = "Exploratory"
    )
  ) +
  scale_shape_manual(values = c("≥80% route coverage" = 18, "<80% route coverage" = 23)) +
  labs(x = "Direction-aligned route effect (z units)", y = NULL) +
  theme_submission() +
  theme(axis.text.y = element_text(size = 5.5, lineheight = 0.88),
        legend.position = "bottom", legend.box = "vertical") +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 1))

# c: component concordance; positive means movement in the frozen expected route direction.
heat_keys <- unname(comparison_labels[1:8])
route_directions <- effect_summary %>%
  filter(feature == "route_score") %>%
  select(dataset, comparison, route_expected = expected_direction)
component_heat <- effect_summary %>%
  filter(feature %in% c("route_score", "route_up", "route_down", "wnt_stem", "proliferation_control")) %>%
  left_join(route_directions, by = c("dataset", "comparison")) %>%
  mutate(
    key = paste(dataset, comparison, sep = "|"),
    aligned = case_when(
      feature == "route_down" ~ -mean_difference * route_expected,
      TRUE ~ mean_difference * route_expected
    ),
    component = recode(
      feature,
      route_score = "Locked\nroute",
      route_up = "Up\narm",
      route_down = "Differentiation\nloss",
      wnt_stem = "WNT/\nstemness",
      proliferation_control = "Proliferation\nco-movement"
    )
  ) %>%
  filter(key %in% heat_keys) %>%
  mutate(
    context = factor(key, levels = rev(heat_keys), labels = label_lookup[rev(heat_keys)]),
    component = factor(component, levels = c("Locked\nroute", "Up\narm", "Differentiation\nloss",
                                             "WNT/\nstemness", "Proliferation\nco-movement"))
  )
write_source(component_heat, "panel_c_component_concordance.tsv")

p_c <- ggplot(component_heat, aes(component, context, fill = aligned)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = COL[["wnt"]], mid = "white", high = COL[["route"]], midpoint = 0,
                       limits = c(-1.6, 1.6), oob = scales::squish,
                       name = "Aligned\neffect") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7, base_family = JTM_FONT) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 5.2, lineheight = 0.85, angle = 28, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5.4, lineheight = 0.88),
    legend.title = element_text(size = 5.8), legend.text = element_text(size = 5.5),
    legend.position = "bottom", plot.margin = margin(1.6, 2, 1.6, 2, "mm"),
    plot.tag = element_text(size = 9, face = "bold")
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                barwidth = grid::unit(18, "mm"),
                                barheight = grid::unit(2.2, "mm")))

# d: pharmacological calibration retains both the successful WNT/stem reversal
# and the discordant full-route reversal.
pharm <- unit_effects %>%
  filter(dataset == "GSE114059",
         comparison %in% c("trametinib_vs_dmso", "pri724_reversal_of_trametinib"),
         feature %in% c("route_score", "wnt_stem")) %>%
  mutate(
    aligned = difference * expected_direction,
    comparison_label = recode(
      comparison,
      trametinib_vs_dmso = "Trametinib\nvs DMSO",
      pri724_reversal_of_trametinib = "PRI-724 + trametinib\nvs trametinib"
    ),
    endpoint = recode(feature, route_score = "Locked route", wnt_stem = "WNT/stemness"),
    row = factor(
      paste(comparison_label, endpoint, sep = "|"),
      levels = rev(c(
        "Trametinib\nvs DMSO|Locked route", "Trametinib\nvs DMSO|WNT/stemness",
        "PRI-724 + trametinib\nvs trametinib|Locked route",
        "PRI-724 + trametinib\nvs trametinib|WNT/stemness"
      ))
    )
  )
pharm_summary <- pharm %>% group_by(row, comparison_label, endpoint) %>%
  summarise(mean_aligned = mean(aligned), .groups = "drop")
write_source(pharm, "panel_d_pharmacological_calibration.tsv")

p_d <- ggplot(pharm, aes(x = aligned, y = row, colour = endpoint)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = COL[["neutral"]]) +
  geom_point(size = 1.45, alpha = 0.85,
             position = position_jitter(height = 0.08, width = 0, seed = 20260808)) +
  geom_point(data = pharm_summary, aes(x = mean_aligned, y = row, fill = endpoint),
             inherit.aes = FALSE, shape = 18, size = 2.4) +
  scale_colour_manual(values = c("Locked route" = COL[["route"]], "WNT/stemness" = COL[["wnt"]])) +
  scale_fill_manual(values = c("Locked route" = COL[["route"]], "WNT/stemness" = COL[["wnt"]])) +
  scale_y_discrete(labels = function(x) gsub("\\|", "\n", x)) +
  labs(x = "Direction-aligned effect", y = NULL) +
  theme_submission() +
  theme(axis.text.y = element_text(size = 5.4, lineheight = 0.86),
        legend.position = "none")

# e: TCF7L2 knockout consensus across empirical route, two regulon scores and
# topology-only edge deletion. Values are scaled within method only for colour.
tcf_route <- effect_summary %>%
  filter(feature == "route_score", grepl("GSE135328", dataset)) %>%
  transmute(context = ifelse(grepl("HT29", dataset), "HT29 KO", "HCT116 KO"),
            method = "Locked route", value = -mean_difference)
tcf_ulm <- ulm_summary %>%
  filter(tf == "TCF7L2", grepl("GSE135328", dataset)) %>%
  transmute(context = ifelse(grepl("HT29", dataset), "HT29 KO", "HCT116 KO"),
            method = "ULM activity", value = -mean_difference)
tcf_signed <- signed_tcf %>%
  mutate(context = ifelse(cell_line == "HT29", "HT29 KO", "HCT116 KO")) %>%
  group_by(context) %>%
  summarise(value = -mean(difference_vs_WT), .groups = "drop") %>%
  mutate(method = "Signed-mean activity")
tcf_graph <- virtual_panel %>% filter(tf == "TCF7L2") %>%
  transmute(context = "Graph deletion", method = "Edge deletion",
            value = -combined_virtual_ko_route_impact)
tcf_consensus <- bind_rows(tcf_route, tcf_ulm, tcf_signed, tcf_graph) %>%
  group_by(method) %>% mutate(scaled = value / max(abs(value), na.rm = TRUE)) %>%
  ungroup()
tcf_grid <- expand_grid(
  context = c("HT29 KO", "HCT116 KO", "Graph deletion"),
  method = c("Locked route", "ULM activity", "Signed-mean activity", "Edge deletion")
) %>%
  left_join(tcf_consensus, by = c("context", "method")) %>%
  mutate(
    context = factor(context, levels = rev(c("HT29 KO", "HCT116 KO", "Graph deletion"))),
    method = factor(method, levels = c("Locked route", "ULM activity", "Signed-mean activity", "Edge deletion")),
    label = ifelse(is.na(value), "", sprintf("%.2f", value))
  )
write_source(tcf_grid, "panel_e_tcf7l2_regulatory_consensus.tsv")

p_e <- ggplot(tcf_grid, aes(method, context, fill = scaled)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 1.85, colour = "black") +
  scale_fill_gradient2(low = COL[["wnt"]], mid = "white", high = COL[["route"]], midpoint = 0,
                       limits = c(-1, 1), na.value = COL[["neutral_pale"]]) +
  scale_x_discrete(labels = c("Locked route" = "Locked\nroute", "ULM activity" = "ULM\nactivity",
                              "Signed-mean activity" = "Signed\nmean",
                              "Edge deletion" = "Edge\ndeletion")) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7, base_family = JTM_FONT) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 5.1, lineheight = 0.82, angle = 0, hjust = 0.5, vjust = 1),
        axis.text.y = element_text(size = 5.6), legend.position = "none",
        plot.margin = margin(1.6, 2, 1.6, 2, "mm"),
        plot.tag = element_text(size = 9, face = "bold"))

# f: paired section-level full-route recapitulation before and after adjustment.
spatial <- spatial_units %>%
  filter(comparison == "tumor_vs_non_neoplastic_epithelium",
         feature %in% c("route_score", "route_residual_prolif_epithelial")) %>%
  mutate(
    feature_label = recode(
      feature,
      route_score = "Raw route\nP = 0.0313",
      route_residual_prolif_epithelial = "Adjusted route\nP = 0.0313"
    ),
    feature_label = factor(feature_label, levels = c("Raw route\nP = 0.0313", "Adjusted route\nP = 0.0313"))
  )
spatial_median <- spatial %>% group_by(feature_label) %>%
  summarise(median_difference = median(difference), .groups = "drop")
write_source(spatial, "panel_f_spatial_section_effects.tsv")

p_f <- ggplot(spatial, aes(feature_label, difference, group = sample_id)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = COL[["neutral"]]) +
  geom_line(colour = COL[["neutral_light"]], linewidth = 0.45) +
  geom_point(colour = COL[["route"]], size = 1.55, alpha = 0.9) +
  geom_point(data = spatial_median,
             aes(feature_label, median_difference), inherit.aes = FALSE,
             shape = 18, size = 2.7, colour = "black") +
  labs(x = NULL, y = "Tumour − non-neoplastic\nepithelium (section median)") +
  theme_submission() +
  theme(axis.text.x = element_text(size = 5.7, lineheight = 0.9), legend.position = "none")

combined <- p_a / (p_b | p_c) / (p_d | p_e | p_f) +
  plot_layout(heights = c(0.55, 2.75, 2.25), widths = c(1.08, 1.08, 0.84)) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e", "f"))) &
  theme(plot.tag.position = c(0.01, 0.98))

panel_trace <- data.frame(
  figure = "Figure 4",
  panel = letters[1:6],
  claim = c(
    "Frozen validation architecture",
    "Direction-aligned genetic and pharmacological route responses",
    "Component concordance and nonoverlapping proliferation control",
    "PDO pharmacological calibration and full-route boundary",
    "TCF7L2 empirical, regulon and graph-deletion consensus",
    "Paired spatial recapitulation before and after control adjustment"
  ),
  source_table = c(
    "panel_a_evidence_architecture.tsv",
    "panel_b_route_effect_summary.tsv; panel_b_route_effects_by_context.tsv",
    "panel_c_component_concordance.tsv",
    "panel_d_pharmacological_calibration.tsv",
    "panel_e_tcf7l2_regulatory_consensus.tsv",
    "panel_f_spatial_section_effects.tsv"
  ),
  independent_unit = c(
    "evidence layer",
    "donor, model, clone, PDO or oncogenic background",
    "dataset contrast",
    "PDO",
    "clone, cell line or regulatory graph",
    "Visium tissue section"
  ),
  inference = c(
    "design schematic",
    "context-level effects, bootstrap intervals and expression-matched null",
    "direction-aligned descriptive response",
    "context-level descriptive calibration",
    "concordance across related regulatory robustness views",
    "exact paired Wilcoxon signed-rank test"
  ),
  stringsAsFactors = FALSE
)
write_source(panel_trace, "panel_source_trace_computational.tsv")

manifest <- export_figure(combined, FIGURE_STEM)
write_source(manifest, "figure_export_manifest.tsv")

cat("Wrote computational-closure figure to", OUT_DIR, "\n")
