#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

file_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("^--file=", "", file_arg[1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
result_dir <- file.path(root, "results", "cell_state_decomposition_v1")
figure_dir <- Sys.getenv(
  "CELL_STATE_FIGURE_DIR",
  unset = file.path(root, "figures", "research_upgrade_v1")
)
figure_stem <- Sys.getenv(
  "CELL_STATE_FIGURE_STEM",
  unset = "figure_cell_state_decomposition_v1"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

within <- read.delim(file.path(result_dir, "within_cell_state_effects.tsv"), check.names = FALSE)
abundance <- read.delim(file.path(result_dir, "differential_abundance_propeller_style.tsv"), check.names = FALSE)
contributions <- read.delim(file.path(result_dir, "decomposition_state_contributions.tsv"), check.names = FALSE)
decomposition <- read.delim(file.path(result_dir, "decomposition_summary.tsv"), check.names = FALSE)
state_scores <- read.delim(gzfile(file.path(result_dir, "specimen_cell_state_scores.tsv.gz")), check.names = FALSE)
donor_inputs <- read.delim(gzfile(file.path(result_dir, "donor_route_decomposition_inputs.tsv.gz")), check.names = FALSE)

state_levels <- c(
  "Neoplastic", "Absorptive", "Goblet", "Transit-amplifying", "Crypt stem",
  "Other colonic epithelium", "Tuft", "Enteroendocrine", "Abnormal"
)
state_levels_reverse <- rev(state_levels)
route_colours <- c("Normal" = "#3976A8", "Adenoma" = "#D65F43")
component_colours <- c(
  "Total difference" = "#273746",
  "Composition" = "#6A8EAE",
  "Within state" = "#D65F43"
)

base_theme <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    axis.title = element_text(size = 8, colour = "#202124"),
    axis.text = element_text(size = 7, colour = "#30343B"),
    axis.line = element_line(linewidth = 0.35, colour = "#40454B"),
    axis.ticks = element_line(linewidth = 0.3, colour = "#40454B"),
    plot.title = element_text(size = 9, face = "bold", colour = "#202124", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 7, colour = "#5A6168", margin = margin(b = 5)),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold", colour = "#30343B"),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.height = unit(3.5, "mm"),
    legend.key.width = unit(4.5, "mm"),
    plot.margin = margin(5, 7, 5, 5)
  )

# a. Donor-balanced composition in the held-out validation partition.
panel_a_data <- subset(
  abundance,
  dataset == "validation" & transformation == "arcsine_sqrt"
)
panel_a_data$cell_type_label <- factor(panel_a_data$cell_type_label, levels = state_levels_reverse)
panel_a_long <- rbind(
  transform(
    panel_a_data[c("cell_type_label", "mean_proportion_normal")],
    route = "Normal", proportion = mean_proportion_normal
  )[c("cell_type_label", "route", "proportion")],
  transform(
    panel_a_data[c("cell_type_label", "mean_proportion_adenoma")],
    route = "Adenoma", proportion = mean_proportion_adenoma
  )[c("cell_type_label", "route", "proportion")]
)
panel_a_long$route <- factor(panel_a_long$route, levels = c("Normal", "Adenoma"))

panel_a <- ggplot(panel_a_data, aes(y = cell_type_label)) +
  geom_segment(
    aes(x = mean_proportion_normal, xend = mean_proportion_adenoma, yend = cell_type_label),
    linewidth = 0.65, colour = "#C9CED3"
  ) +
  geom_point(
    data = panel_a_long,
    aes(x = proportion, colour = route),
    size = 2.25
  ) +
  scale_colour_manual(values = route_colours) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Epithelial composition changes",
    x = "Mean donor-level proportion", y = NULL
  ) +
  base_theme +
  theme(legend.position = "top", legend.justification = "left")

# b. State-resolved frozen-core effects; sparse estimates remain visible but muted.
panel_b_data <- subset(
  within,
  dataset == "validation" & score == "core_287" &
    min_cells_per_specimen_state == 20 & comparison == "conventional_vs_normal"
)
panel_b_data$cell_type_label <- factor(panel_b_data$cell_type_label, levels = state_levels_reverse)
panel_b_data$support <- ifelse(
  panel_b_data$n_donors_a >= 5 & panel_b_data$n_donors_b >= 5,
  "At least 5 donors per group", "Sparse"
)
panel_b_data$donor_label <- paste0(panel_b_data$n_donors_a, "/", panel_b_data$n_donors_b)
panel_b_data$display_label <- paste0(as.character(panel_b_data$cell_type_label), "  (", panel_b_data$donor_label, ")")
panel_b_data$display_label <- factor(
  panel_b_data$display_label,
  levels = paste0(state_levels_reverse, "  (", panel_b_data$donor_label[match(state_levels_reverse, as.character(panel_b_data$cell_type_label))], ")")
)

panel_b <- ggplot(panel_b_data, aes(y = display_label, x = effect_a_minus_b)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#A7ADB3") +
  geom_errorbar(
    aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high, alpha = support),
    orientation = "y", width = 0, linewidth = 0.65, colour = "#D65F43"
  ) +
  geom_point(aes(alpha = support), size = 2.25, colour = "#D65F43") +
  scale_alpha_manual(values = c("At least 5 donors per group" = 1, "Sparse" = 0.35), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.06))) +
  labs(
    title = "Programme shift within matched states",
    subtitle = "Parentheses show adenoma/normal donors",
    x = "Adenoma − normal 287-gene score", y = NULL
  ) +
  base_theme

# c. Paired donors in the three well-supported epithelial states.
paired_states <- c("ABS", "GOB", "TAC")
paired <- subset(
  state_scores,
  dataset == "validation" & cell_type %in% paired_states & n_cells >= 20 &
    route_group %in% c("normal", "conventional_adenoma")
)
paired <- aggregate(
  paired[["score__core_287__mean"]],
  by = list(
    donor_id = paired$donor_id,
    route_group = paired$route_group,
    cell_type = paired$cell_type,
    cell_type_label = paired$cell_type_label
  ),
  FUN = mean
)
names(paired)[5] <- "score"
pair_counts <- aggregate(route_group ~ donor_id + cell_type, paired, function(x) length(unique(x)))
eligible_pairs <- subset(pair_counts, route_group == 2)[c("donor_id", "cell_type")]
paired <- merge(paired, eligible_pairs, by = c("donor_id", "cell_type"))
paired$route_label <- factor(
  ifelse(paired$route_group == "normal", "Normal", "Adenoma"),
  levels = c("Normal", "Adenoma")
)
paired$cell_type_label <- factor(
  ifelse(paired$cell_type_label == "Transit-amplifying", "Transit-\namplifying", paired$cell_type_label),
  levels = c("Absorptive", "Goblet", "Transit-\namplifying")
)

panel_c <- ggplot(paired, aes(x = route_label, y = score, group = donor_id)) +
  geom_line(linewidth = 0.45, colour = "#AEB4BA", alpha = 0.85) +
  geom_point(aes(colour = route_label), size = 1.8) +
  facet_wrap(~cell_type_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = route_colours, guide = "none") +
  scale_x_discrete(labels = c("Normal" = "N", "Adenoma" = "A")) +
  labs(
    title = "Paired donors retain the shift",
    x = NULL, y = "287-gene score"
  ) +
  base_theme +
  theme(
    panel.spacing.x = unit(4, "mm"),
    strip.text.x = element_text(size = 6.5, lineheight = 0.9),
    strip.clip = "off"
  )

# d. Signed state-specific contributions explain what drives each component.
panel_d_data <- subset(
  contributions,
  dataset == "validation" & score == "core_287" & state_set == "all_states" &
    comparison == "conventional_vs_normal"
)
composition_rows <- transform(
  panel_d_data,
  component = "Composition",
  estimate = composition_contribution,
  low = composition_ci_low,
  high = composition_ci_high
)
within_rows <- transform(
  panel_d_data,
  component = "Within state",
  estimate = within_state_contribution,
  low = within_state_ci_low,
  high = within_state_ci_high
)
panel_d_long <- rbind(composition_rows, within_rows)
panel_d_long$cell_type_label <- factor(panel_d_long$cell_type_label, levels = state_levels_reverse)
panel_d_long$component <- factor(panel_d_long$component, levels = c("Composition", "Within state"))

panel_d <- ggplot(panel_d_long, aes(x = estimate, y = cell_type_label, colour = component)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#A7ADB3") +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 1.9) +
  facet_wrap(~component, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = component_colours[c("Composition", "Within state")], guide = "none") +
  scale_x_continuous(breaks = pretty_breaks(n = 3), labels = label_number(accuracy = 0.01)) +
  labs(
    title = "State-specific component contributions",
    x = "Signed contribution to total difference", y = NULL
  ) +
  base_theme +
  theme(panel.spacing.x = unit(5, "mm"))

# e. Overall decomposition and its prespecified validation checks.
panel_e_source <- subset(
  decomposition,
  score == "core_287" & comparison == "conventional_vs_normal" &
    ((dataset == "validation" & state_set %in% c("all_states", "canonical_states")) |
      (dataset == "discovery" & state_set == "all_states"))
)
panel_e_source$analysis <- with(
  panel_e_source,
  ifelse(
    dataset == "validation" & state_set == "all_states", "Held-out: all states",
    ifelse(dataset == "validation", "Held-out: canonical states", "Discovery: all states")
  )
)
panel_e <- rbind(
  transform(
    panel_e_source,
    component = "Total difference", estimate = total_difference,
    low = total_ci_low, high = total_ci_high
  ),
  transform(
    panel_e_source,
    component = "Composition", estimate = composition_component,
    low = composition_ci_low, high = composition_ci_high
  ),
  transform(
    panel_e_source,
    component = "Within state", estimate = within_state_component,
    low = within_state_ci_low, high = within_state_ci_high
  )
)
panel_e$analysis <- factor(
  panel_e$analysis,
  levels = rev(c("Held-out: all states", "Held-out: canonical states", "Discovery: all states"))
)
panel_e$component <- factor(panel_e$component, levels = names(component_colours))

panel_e_plot <- ggplot(panel_e, aes(y = analysis, x = estimate, colour = component)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#A7ADB3") +
  geom_errorbar(
    aes(xmin = low, xmax = high),
    orientation = "y", position = position_dodge(width = 0.45), width = 0, linewidth = 0.55
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 1.9) +
  scale_colour_manual(values = component_colours) +
  labs(
    title = "Exact decomposition replicates",
    x = "Adenoma − normal 287-gene score", y = NULL
  ) +
  base_theme +
  theme(legend.position = "top", legend.justification = "left")

# f. Compact-score transfer at the donor-route level.
panel_f_data <- subset(
  donor_inputs,
  dataset == "validation" & state_set == "all_states" &
    route_group %in% c("normal", "conventional_adenoma")
)
panel_f_totals <- aggregate(
  contribution ~ donor_id + route_group + score,
  panel_f_data,
  sum
)
panel_f_wide <- reshape(
  panel_f_totals,
  idvar = c("donor_id", "route_group"), timevar = "score", direction = "wide"
)
names(panel_f_wide) <- sub("^contribution\\.", "", names(panel_f_wide))
panel_f_wide$route_label <- factor(
  ifelse(panel_f_wide$route_group == "normal", "Normal", "Adenoma"),
  levels = c("Normal", "Adenoma")
)
rho <- cor(panel_f_wide$core_287, panel_f_wide$signature_12, method = "spearman")

panel_f <- ggplot(panel_f_wide, aes(x = core_287, y = signature_12, colour = route_label)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6, colour = "#8B9298") +
  geom_point(size = 2.0, alpha = 0.9) +
  scale_colour_manual(values = route_colours) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("Spearman ρ = %.3f", rho),
    hjust = -0.05, vjust = 1.25, size = 2.5, family = "Arial", colour = "#30343B"
  ) +
  labs(
    title = "The 12-gene score preserves the state",
    x = "287-gene donor-route score", y = "12-gene donor-route score"
  ) +
  base_theme +
  theme(legend.position = "top", legend.justification = "left")

figure <- ((panel_a | panel_b) + plot_layout(widths = c(0.92, 1.08))) /
  ((panel_c | panel_d) + plot_layout(widths = c(0.95, 1.05))) /
  ((panel_e_plot | panel_f) + plot_layout(widths = c(1.08, 0.92))) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(family = "Arial", size = 11, face = "bold", colour = "#15191D"),
    plot.tag.position = c(0, 1)
  )

pdf_path <- file.path(figure_dir, paste0(figure_stem, ".pdf"))
tiff_path <- file.path(figure_dir, paste0(figure_stem, ".tiff"))
png_path <- file.path(figure_dir, paste0(figure_stem, ".png"))
svg_path <- file.path(figure_dir, paste0(figure_stem, ".svg"))

ggsave(pdf_path, figure, width = 180, height = 226, units = "mm", device = cairo_pdf)
ggsave(tiff_path, figure, width = 180, height = 226, units = "mm", dpi = 600, compression = "lzw")
ggsave(png_path, figure, width = 180, height = 226, units = "mm", dpi = 300, bg = "white")
ggsave(svg_path, figure, width = 180, height = 226, units = "mm", device = svglite::svglite)

cat(sprintf("Saved %s\n", pdf_path))
