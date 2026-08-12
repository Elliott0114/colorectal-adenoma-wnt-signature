#!/usr/bin/env Rscript

# Discovery-only pilot for replacing the fixed 50-up/50-down truncation.
#
# The analysis uses a donor-aware limma model plus the existing donor-block
# bootstrap audit to define an error-controlled stable core.  A forward,
# redundancy-aware score-reconstruction path is then fit within grouped
# leave-one-donor-out cross-validation.  The one-standard-error rule determines
# the number of genes per arm.  Chen validation data are read only after the
# discovery panel has been fixed.

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
})

set.seed(20260710)

root <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(root, "results", "data_adaptive_panel_pilot_v2_6")
fig_dir <- file.path(root, "figures", "data_adaptive_panel_pilot_v2_6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

discovery_meta_path <- file.path(
  root, "results", "route_signature",
  "chen_discovery_specimen_pseudobulk_meta.tsv"
)
discovery_expr_path <- file.path(
  root, "results", "route_signature",
  "chen_discovery_specimen_pseudobulk_expression.tsv.gz"
)
validation_meta_path <- file.path(
  root, "results", "route_signature",
  "chen_validation_specimen_pseudobulk_meta.tsv"
)
validation_expr_path <- file.path(
  root, "results", "route_signature",
  "chen_validation_specimen_pseudobulk_expression.tsv.gz"
)
audit_path <- file.path(
  root, "results", "route_signature_locked",
  "discovery_gene_stability_audit.tsv"
)
reference_100_path <- file.path(
  root, "results", "route_signature_locked",
  "discovery_locked_signature_genes.tsv"
)
compact_10_path <- file.path(
  root, "results", "programme_transparency_v2_5",
  "compact_panel_definition_corrected.tsv"
)

read_tsv_base <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else path
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  read.delim(
    con,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

write_tsv <- function(x, name) {
  write.table(
    x,
    file.path(out_dir, name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

cell_cycle_exclude <- c(
  "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
  "MCM7", "TYMS", "UBE2C", "CENPF", "CDK1", "CCNA2", "CCNB1", "CCNB2",
  "AURKA", "AURKB", "BIRC5", "CDC20", "CDC45", "CDC6", "MKI67IP"
)

is_excluded_gene <- function(gene) {
  upper <- toupper(gene)
  starts_technical <- grepl("^(ENSG|MT-|RPL|RPS|MTRNR)", upper)
  exact_technical <- upper %in% c(
    cell_cycle_exclude,
    "MALAT1", "XIST", "NEAT1", "JUN", "FOS", "HBB", "HBA1", "HBA2"
  )
  starts_technical | exact_technical
}

aggregate_donor_route <- function(meta, expr) {
  keep <- meta$route_group %in% c("normal", "conventional_adenoma")
  meta <- meta[keep, , drop = FALSE]
  expr <- as.matrix(expr[keep, , drop = FALSE])
  storage.mode(expr) <- "double"
  key <- paste(meta$donor_id, meta$route_group, sep = "::")
  levels_key <- sort(unique(key))
  rows <- vector("list", length(levels_key))
  key_meta <- vector("list", length(levels_key))
  for (i in seq_along(levels_key)) {
    idx <- which(key == levels_key[[i]])
    if (length(idx) == 1L) {
      rows[[i]] <- expr[idx, ]
    } else {
      rows[[i]] <- apply(expr[idx, , drop = FALSE], 2, median, na.rm = TRUE)
    }
    key_meta[[i]] <- data.frame(
      donor_id = meta$donor_id[idx[[1]]],
      route_group = meta$route_group[idx[[1]]],
      n_specimens = length(idx),
      stringsAsFactors = FALSE
    )
  }
  aggregated <- do.call(rbind, rows)
  colnames(aggregated) <- colnames(expr)
  rownames(aggregated) <- levels_key
  list(meta = do.call(rbind, key_meta), expr = aggregated)
}

zscore_matrix <- function(x, center = NULL, scale_value = NULL) {
  x <- as.matrix(x)
  if (is.null(center)) center <- colMeans(x, na.rm = TRUE)
  if (is.null(scale_value)) {
    scale_value <- apply(x, 2, sd, na.rm = TRUE)
  }
  scale_value[!is.finite(scale_value) | scale_value == 0] <- NA_real_
  z <- sweep(sweep(x, 2, center, "-"), 2, scale_value, "/")
  list(z = z, center = center, scale = scale_value)
}

safe_cor <- function(x, y, method = "pearson") {
  value <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = method))
  ifelse(is.finite(value), value, NA_real_)
}

greedy_reconstruction_path <- function(x, target, max_k) {
  x <- as.matrix(x)
  max_k <- min(max_k, ncol(x))
  selected <- character(0)
  selected_idx <- integer(0)
  current_sum <- rep(0, nrow(x))
  path_cor <- numeric(max_k)
  for (k in seq_len(max_k)) {
    remaining <- setdiff(seq_len(ncol(x)), selected_idx)
    candidate_scores <- sweep(x[, remaining, drop = FALSE], 1, current_sum, "+") / k
    correlations <- as.numeric(cor(candidate_scores, target, use = "pairwise.complete.obs"))
    correlations[!is.finite(correlations)] <- -Inf
    ordering <- order(-correlations, colnames(x)[remaining])
    chosen <- remaining[ordering[[1]]]
    selected_idx <- c(selected_idx, chosen)
    selected <- c(selected, colnames(x)[chosen])
    current_sum <- current_sum + x[, chosen]
    path_cor[[k]] <- safe_cor(current_sum / k, target)
  }
  data.frame(
    k = seq_len(max_k),
    gene = selected,
    training_correlation = path_cor,
    stringsAsFactors = FALSE
  )
}

bootstrap_oof_curve <- function(oof, donors, arm, max_k, n_boot = 1000L) {
  target <- oof[[paste0("full_", arm)]]
  point <- vapply(
    seq_len(max_k),
    function(k) safe_cor(oof[[paste0(arm, "_k", k)]], target, method = "spearman"),
    numeric(1)
  )
  unique_donors <- sort(unique(donors))
  boot <- matrix(NA_real_, nrow = n_boot, ncol = max_k)
  for (b in seq_len(n_boot)) {
    sampled <- sample(unique_donors, length(unique_donors), replace = TRUE)
    idx <- unlist(lapply(sampled, function(d) which(donors == d)), use.names = FALSE)
    for (k in seq_len(max_k)) {
      boot[b, k] <- safe_cor(
        oof[[paste0(arm, "_k", k)]][idx],
        target[idx],
        method = "spearman"
      )
    }
  }
  se <- apply(boot, 2, sd, na.rm = TRUE)
  lower <- apply(boot, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper <- apply(boot, 2, quantile, probs = 0.975, na.rm = TRUE)
  k_max <- which.max(point)
  threshold <- point[[k_max]] - se[[k_max]]
  selected_k <- min(which(point >= threshold))
  curve <- data.frame(
    arm = arm,
    k = seq_len(max_k),
    oof_spearman = point,
    bootstrap_se = se,
    ci_low = lower,
    ci_high = upper,
    maximum_k = k_max,
    one_se_threshold = threshold,
    selected_k = selected_k,
    stringsAsFactors = FALSE
  )
  list(curve = curve, bootstrap = boot, selected_k = selected_k)
}

score_two_arm <- function(z, up, down) {
  up_present <- intersect(up, colnames(z))
  down_present <- intersect(down, colnames(z))
  if (length(up_present) == 0 || length(down_present) == 0) {
    stop("Both score arms require at least one measured gene")
  }
  rowMeans(z[, up_present, drop = FALSE], na.rm = TRUE) -
    rowMeans(z[, down_present, drop = FALSE], na.rm = TRUE)
}

auc_rank <- function(case, control) {
  case <- case[is.finite(case)]
  control <- control[is.finite(control)]
  ranks <- rank(c(case, control), ties.method = "average")
  u <- sum(ranks[seq_along(case)]) - length(case) * (length(case) + 1) / 2
  u / (length(case) * length(control))
}

panel_metrics <- function(meta, score, full_score, panel_name, gene_count) {
  use <- meta$route_group %in% c("normal", "conventional_adenoma")
  m <- meta[use, , drop = FALSE]
  s <- score[use]
  f <- full_score[use]
  adenoma <- s[m$route_group == "conventional_adenoma"]
  normal <- s[m$route_group == "normal"]
  test <- suppressWarnings(wilcox.test(adenoma, normal, exact = FALSE))
  donor_values <- data.frame(
    donor_id = m$donor_id,
    route_group = m$route_group,
    score = s,
    stringsAsFactors = FALSE
  ) |>
    group_by(donor_id, route_group) |>
    summarise(score = median(score), .groups = "drop") |>
    pivot_wider(names_from = route_group, values_from = score)
  paired <- donor_values |>
    filter(is.finite(normal), is.finite(conventional_adenoma)) |>
    mutate(delta = conventional_adenoma - normal)
  paired_p <- if (nrow(paired) >= 3) {
    suppressWarnings(wilcox.test(paired$delta, mu = 0, exact = FALSE)$p.value)
  } else {
    NA_real_
  }
  data.frame(
    panel = panel_name,
    n_genes = gene_count,
    n_adenoma = length(adenoma),
    n_normal = length(normal),
    median_difference = median(adenoma) - median(normal),
    auc = auc_rank(adenoma, normal),
    p_mann_whitney = test$p.value,
    spearman_vs_full_core = safe_cor(s, f, method = "spearman"),
    n_paired_donors = nrow(paired),
    median_paired_delta = if (nrow(paired)) median(paired$delta) else NA_real_,
    paired_positive_fraction = if (nrow(paired)) mean(paired$delta > 0) else NA_real_,
    p_paired_wilcoxon = paired_p,
    stringsAsFactors = FALSE
  )
}

cluster_bootstrap_auc_difference <- function(meta, score_a, score_b, n_boot = 2000L) {
  use <- meta$route_group %in% c("normal", "conventional_adenoma")
  m <- meta[use, , drop = FALSE]
  a <- score_a[use]
  b <- score_b[use]
  donors <- sort(unique(m$donor_id))
  values <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    sampled <- sample(donors, length(donors), replace = TRUE)
    idx <- unlist(lapply(sampled, function(d) which(m$donor_id == d)), use.names = FALSE)
    group <- m$route_group[idx]
    if (sum(group == "conventional_adenoma") == 0 || sum(group == "normal") == 0) next
    auc_a <- auc_rank(a[idx][group == "conventional_adenoma"], a[idx][group == "normal"])
    auc_b <- auc_rank(b[idx][group == "conventional_adenoma"], b[idx][group == "normal"])
    values[[i]] <- auc_a - auc_b
  }
  values <- values[is.finite(values)]
  data.frame(
    comparison = "adaptive_minus_100_gene_reference",
    n_bootstrap_valid = length(values),
    auc_difference_median = median(values),
    ci_low = unname(quantile(values, 0.025)),
    ci_high = unname(quantile(values, 0.975)),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Discovery-only model and stable evidence core
# -----------------------------------------------------------------------------

message("Loading discovery data")
discovery_meta <- read_tsv_base(discovery_meta_path)
discovery_expr <- read_tsv_base(discovery_expr_path)
stopifnot(nrow(discovery_meta) == nrow(discovery_expr))
discovery_agg <- aggregate_donor_route(discovery_meta, discovery_expr)
rm(discovery_expr)
gc(verbose = FALSE)

audit <- read_tsv_base(audit_path)
as_logical_field <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}
for (field in c("excluded_gene", "ci_excludes_zero", "eligible", "selected")) {
  if (field %in% colnames(audit)) {
    audit[[field]] <- as_logical_field(audit[[field]])
  }
}
expressed_genes <- audit$gene[
  !audit$excluded_gene & audit$mean_expression_discovery > 0.001
]
expressed_genes <- intersect(expressed_genes, colnames(discovery_agg$expr))
y <- t(discovery_agg$expr[, expressed_genes, drop = FALSE])
route <- factor(
  discovery_agg$meta$route_group,
  levels = c("normal", "conventional_adenoma")
)
design <- model.matrix(~ route)
block <- factor(discovery_agg$meta$donor_id)
corfit <- duplicateCorrelation(y, design = design, block = block)
fit <- lmFit(y, design = design, block = block, correlation = corfit$consensus)
fit <- eBayes(fit, robust = TRUE, trend = TRUE)
limma_table <- topTable(
  fit,
  coef = "routeconventional_adenoma",
  number = Inf,
  sort.by = "none"
)
limma_table$gene <- rownames(limma_table)
rownames(limma_table) <- NULL
limma_table <- limma_table |>
  select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, B)

evidence <- audit |>
  inner_join(limma_table, by = "gene") |>
  mutate(
    limma_direction_match = sign(logFC) == sign(discovery_effect_adenoma_minus_normal),
    limma_fdr_pass = adj.P.Val <= 0.05,
    bootstrap_stability_pass = direction_stability >= 0.90 & ci_excludes_zero,
    stable_core = !excluded_gene & mean_expression_discovery > 0.001 &
      limma_fdr_pass & bootstrap_stability_pass & limma_direction_match,
    arm = ifelse(logFC > 0, "up", "down")
  ) |>
  arrange(adj.P.Val, desc(abs(logFC)), gene)

core <- evidence |>
  filter(stable_core)
if (sum(core$arm == "up") < 5 || sum(core$arm == "down") < 5) {
  stop("The stable core contains fewer than five genes in one arm")
}

write_tsv(evidence, "discovery_gene_evidence.tsv")
write_tsv(core, "stable_error_controlled_core.tsv")

message(
  "Stable core: ", nrow(core), " genes (",
  sum(core$arm == "up"), " up; ", sum(core$arm == "down"), " down)"
)

# -----------------------------------------------------------------------------
# Grouped leave-one-donor-out score reconstruction
# -----------------------------------------------------------------------------

core_up <- sort(core$gene[core$arm == "up"])
core_down <- sort(core$gene[core$arm == "down"])
core_genes <- c(core_up, core_down)
max_k <- min(30L, length(core_up), length(core_down))
n <- nrow(discovery_agg$expr)
oof <- data.frame(
  row_id = seq_len(n),
  donor_id = discovery_agg$meta$donor_id,
  route_group = discovery_agg$meta$route_group,
  full_up = NA_real_,
  full_down = NA_real_,
  stringsAsFactors = FALSE
)
for (k in seq_len(max_k)) {
  oof[[paste0("up_k", k)]] <- NA_real_
  oof[[paste0("down_k", k)]] <- NA_real_
}

donors <- sort(unique(discovery_agg$meta$donor_id))
message("Running grouped leave-one-donor-out reconstruction across ", length(donors), " donors")
for (fold in seq_along(donors)) {
  held_donor <- donors[[fold]]
  test_idx <- which(discovery_agg$meta$donor_id == held_donor)
  train_idx <- setdiff(seq_len(n), test_idx)
  train_raw <- discovery_agg$expr[train_idx, core_genes, drop = FALSE]
  test_raw <- discovery_agg$expr[test_idx, core_genes, drop = FALSE]
  train_scaled <- zscore_matrix(train_raw)
  z_train <- train_scaled$z
  z_test <- zscore_matrix(
    test_raw,
    center = train_scaled$center,
    scale_value = train_scaled$scale
  )$z
  full_train_up <- rowMeans(z_train[, core_up, drop = FALSE])
  full_train_down <- rowMeans(z_train[, core_down, drop = FALSE])
  oof$full_up[test_idx] <- rowMeans(z_test[, core_up, drop = FALSE])
  oof$full_down[test_idx] <- rowMeans(z_test[, core_down, drop = FALSE])
  path_up <- greedy_reconstruction_path(
    z_train[, core_up, drop = FALSE], full_train_up, max_k
  )
  path_down <- greedy_reconstruction_path(
    z_train[, core_down, drop = FALSE], full_train_down, max_k
  )
  for (k in seq_len(max_k)) {
    oof[[paste0("up_k", k)]][test_idx] <- rowMeans(
      z_test[, path_up$gene[seq_len(k)], drop = FALSE]
    )
    oof[[paste0("down_k", k)]][test_idx] <- rowMeans(
      z_test[, path_down$gene[seq_len(k)], drop = FALSE]
    )
  }
  if (fold %% 5 == 0 || fold == length(donors)) {
    message("  completed donor fold ", fold, "/", length(donors))
  }
}

up_curve <- bootstrap_oof_curve(oof, oof$donor_id, "up", max_k, n_boot = 1000L)
down_curve <- bootstrap_oof_curve(oof, oof$donor_id, "down", max_k, n_boot = 1000L)
fidelity_curve <- bind_rows(up_curve$curve, down_curve$curve)
selected_k_up <- up_curve$selected_k
selected_k_down <- down_curve$selected_k

full_scaled <- zscore_matrix(discovery_agg$expr[, core_genes, drop = FALSE])$z
final_path_up <- greedy_reconstruction_path(
  full_scaled[, core_up, drop = FALSE],
  rowMeans(full_scaled[, core_up, drop = FALSE]),
  max_k
)
final_path_down <- greedy_reconstruction_path(
  full_scaled[, core_down, drop = FALSE],
  rowMeans(full_scaled[, core_down, drop = FALSE]),
  max_k
)

adaptive_up <- final_path_up$gene[seq_len(selected_k_up)]
adaptive_down <- final_path_down$gene[seq_len(selected_k_down)]
adaptive_panel <- bind_rows(
  data.frame(
    gene = adaptive_up,
    arm = "up",
    selection_order_within_arm = seq_along(adaptive_up),
    selected_k_by_one_se = selected_k_up,
    stringsAsFactors = FALSE
  ),
  data.frame(
    gene = adaptive_down,
    arm = "down",
    selection_order_within_arm = seq_along(adaptive_down),
    selected_k_by_one_se = selected_k_down,
    stringsAsFactors = FALSE
  )
) |>
  left_join(
    evidence |>
      select(
        gene, logFC, adj.P.Val, direction_stability,
        discovery_effect_adenoma_minus_normal,
        bootstrap_ci_low, bootstrap_ci_high
      ),
    by = "gene"
  ) |>
  arrange(desc(arm), selection_order_within_arm)

oof$full_score <- oof$full_up - oof$full_down
oof$adaptive_score <-
  oof[[paste0("up_k", selected_k_up)]] -
  oof[[paste0("down_k", selected_k_down)]]
oof_total_fidelity <- safe_cor(oof$adaptive_score, oof$full_score, method = "spearman")

write_tsv(oof, "discovery_grouped_oof_scores.tsv")
write_tsv(fidelity_curve, "discovery_oof_fidelity_curve.tsv")
write_tsv(final_path_up, "full_discovery_forward_path_up.tsv")
write_tsv(final_path_down, "full_discovery_forward_path_down.tsv")
write_tsv(adaptive_panel, "adaptive_minimum_sufficient_panel.tsv")

# -----------------------------------------------------------------------------
# Read held-out data only after the adaptive panel has been fixed
# -----------------------------------------------------------------------------

message(
  "Adaptive panel fixed before validation: ", length(adaptive_up), " up + ",
  length(adaptive_down), " down"
)
validation_meta <- read_tsv_base(validation_meta_path)
validation_expr <- read_tsv_base(validation_expr_path)
stopifnot(nrow(validation_meta) == nrow(validation_expr))

reference_100 <- read_tsv_base(reference_100_path)
compact_10 <- read_tsv_base(compact_10_path)
reference_up <- reference_100$gene[reference_100$signature_direction == "adenoma_up"]
reference_down <- reference_100$gene[reference_100$signature_direction == "adenoma_down"]
compact_up <- compact_10$gene[compact_10$route_weight > 0]
compact_down <- compact_10$gene[compact_10$route_weight < 0]

needed <- unique(c(
  core_up, core_down, adaptive_up, adaptive_down,
  reference_up, reference_down, compact_up, compact_down
))
missing_validation <- setdiff(needed, colnames(validation_expr))
if (length(missing_validation)) {
  stop("Chen validation is missing selected genes: ", paste(missing_validation, collapse = ", "))
}
validation_z <- zscore_matrix(validation_expr[, needed, drop = FALSE])$z
scores <- data.frame(
  validation_meta,
  full_stable_core = score_two_arm(validation_z, core_up, core_down),
  adaptive_minimum_sufficient = score_two_arm(validation_z, adaptive_up, adaptive_down),
  reference_100_gene = score_two_arm(validation_z, reference_up, reference_down),
  compact_10_gene = score_two_arm(validation_z, compact_up, compact_down),
  stringsAsFactors = FALSE
)

metrics <- bind_rows(
  panel_metrics(
    validation_meta, scores$full_stable_core, scores$full_stable_core,
    "Full stable core", length(core_genes)
  ),
  panel_metrics(
    validation_meta, scores$adaptive_minimum_sufficient, scores$full_stable_core,
    "Adaptive minimum-sufficient", length(adaptive_up) + length(adaptive_down)
  ),
  panel_metrics(
    validation_meta, scores$reference_100_gene, scores$full_stable_core,
    "100-gene reference", 100L
  ),
  panel_metrics(
    validation_meta, scores$compact_10_gene, scores$full_stable_core,
    "10-gene biology-guided", 10L
  )
)

auc_difference <- cluster_bootstrap_auc_difference(
  validation_meta,
  scores$adaptive_minimum_sufficient,
  scores$reference_100_gene,
  n_boot = 2000L
)

write_tsv(scores, "chen_heldout_panel_scores.tsv")
write_tsv(metrics, "chen_heldout_panel_comparison.tsv")
write_tsv(auc_difference, "chen_heldout_auc_difference_bootstrap.tsv")

flow <- data.frame(
  stage_order = 1:5,
  stage = c(
    "Assayed epithelial features",
    "Expressed non-technical genes",
    "limma FDR <= 0.05",
    "FDR + bootstrap stable core",
    "Adaptive minimum-sufficient panel"
  ),
  n_genes = c(
    nrow(audit),
    length(expressed_genes),
    sum(evidence$limma_fdr_pass & !evidence$excluded_gene &
      evidence$mean_expression_discovery > 0.001),
    nrow(core),
    nrow(adaptive_panel)
  ),
  up_genes = c(
    NA_integer_, NA_integer_,
    sum(evidence$limma_fdr_pass & !evidence$excluded_gene &
      evidence$mean_expression_discovery > 0.001 & evidence$logFC > 0),
    sum(core$arm == "up"),
    length(adaptive_up)
  ),
  down_genes = c(
    NA_integer_, NA_integer_,
    sum(evidence$limma_fdr_pass & !evidence$excluded_gene &
      evidence$mean_expression_discovery > 0.001 & evidence$logFC < 0),
    sum(core$arm == "down"),
    length(adaptive_down)
  ),
  stringsAsFactors = FALSE
)
write_tsv(flow, "adaptive_selection_flow.tsv")

summary_table <- data.frame(
  metric = c(
    "duplicate_correlation",
    "stable_core_genes",
    "stable_core_up",
    "stable_core_down",
    "adaptive_up_genes",
    "adaptive_down_genes",
    "adaptive_total_genes",
    "discovery_oof_spearman_vs_full_core",
    "heldout_adaptive_auc",
    "heldout_reference_100_auc",
    "heldout_compact_10_auc",
    "heldout_adaptive_spearman_vs_full_core",
    "heldout_adaptive_paired_positive_fraction",
    "heldout_auc_difference_adaptive_minus_100_ci_low",
    "heldout_auc_difference_adaptive_minus_100_ci_high"
  ),
  value = c(
    corfit$consensus,
    nrow(core),
    sum(core$arm == "up"),
    sum(core$arm == "down"),
    length(adaptive_up),
    length(adaptive_down),
    length(adaptive_up) + length(adaptive_down),
    oof_total_fidelity,
    metrics$auc[metrics$panel == "Adaptive minimum-sufficient"],
    metrics$auc[metrics$panel == "100-gene reference"],
    metrics$auc[metrics$panel == "10-gene biology-guided"],
    metrics$spearman_vs_full_core[metrics$panel == "Adaptive minimum-sufficient"],
    metrics$paired_positive_fraction[metrics$panel == "Adaptive minimum-sufficient"],
    auc_difference$ci_low,
    auc_difference$ci_high
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary_table, "pilot_summary.tsv")

# -----------------------------------------------------------------------------
# R-only exploratory figure
# -----------------------------------------------------------------------------

palette <- c(
  navy = "#315A7D",
  teal = "#3B8C88",
  orange = "#D98B3A",
  red = "#B94A48",
  grey = "#8A9096",
  light = "#DDE3E8",
  dark = "#262A2E"
)

theme_pilot <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = palette[["dark"]]),
      axis.ticks = element_line(linewidth = 0.35, colour = palette[["dark"]]),
      axis.text = element_text(colour = palette[["dark"]]),
      axis.title = element_text(colour = palette[["dark"]]),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 0.5),
      plot.subtitle = element_text(size = base_size - 0.2, colour = palette[["grey"]]),
      panel.grid = element_blank()
    )
}

p_a <- flow |>
  mutate(
    stage = factor(stage, levels = rev(stage)),
    label = format(n_genes, big.mark = ",")
  ) |>
  ggplot(aes(x = n_genes, y = stage)) +
  geom_segment(aes(x = 1, xend = n_genes, yend = stage), linewidth = 2.7,
               colour = palette[["light"]], lineend = "round") +
  geom_point(size = 2.2, colour = palette[["navy"]]) +
  geom_text(aes(label = label), hjust = -0.15, size = 2.3, colour = palette[["dark"]]) +
  scale_x_log10(expand = expansion(mult = c(0.02, 0.23))) +
  labs(
    title = "Discovery-only selection",
    subtitle = "No fixed gene-count truncation",
    x = "Number of genes (log scale)", y = NULL
  ) +
  theme_pilot()

curve_plot <- fidelity_curve |>
  mutate(
    arm_label = recode(arm, up = "Adenoma-up arm", down = "Adenoma-down arm")
  )
p_b <- ggplot(curve_plot, aes(x = k, y = oof_spearman, colour = arm_label)) +
  geom_ribbon(
    aes(ymin = ci_low, ymax = ci_high, fill = arm_label),
    alpha = 0.12, colour = NA
  ) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.1) +
  geom_vline(
    data = distinct(curve_plot, arm_label, selected_k),
    aes(xintercept = selected_k, colour = arm_label),
    linetype = 2, linewidth = 0.45, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Adenoma-up arm" = palette[["orange"]],
    "Adenoma-down arm" = palette[["navy"]]
  )) +
  scale_fill_manual(values = c(
    "Adenoma-up arm" = palette[["orange"]],
    "Adenoma-down arm" = palette[["navy"]]
  )) +
  coord_cartesian(ylim = c(0.45, 1.01)) +
  labs(
    title = "Grouped out-of-fold fidelity",
    subtitle = "Dashed lines: one-SE choices",
    x = "Genes retained per arm", y = "Spearman correlation"
  ) +
  theme_pilot() +
  theme(legend.position = c(0.67, 0.22))

metric_long <- metrics |>
  select(panel, n_genes, auc, spearman_vs_full_core, paired_positive_fraction) |>
  pivot_longer(
    cols = c(auc, spearman_vs_full_core, paired_positive_fraction),
    names_to = "metric", values_to = "value"
  ) |>
  mutate(
    metric = recode(
      metric,
      auc = "Held-out AUC",
      spearman_vs_full_core = "Held-out fidelity",
      paired_positive_fraction = "Positive paired donors"
    ),
    panel = factor(
      panel,
      levels = c(
        "Full stable core", "Adaptive minimum-sufficient",
        "100-gene reference", "10-gene biology-guided"
      )
    )
  )
panel_colours <- c(
  "Full stable core" = palette[["grey"]],
  "Adaptive minimum-sufficient" = palette[["teal"]],
  "100-gene reference" = palette[["navy"]],
  "10-gene biology-guided" = palette[["orange"]]
)
p_c <- ggplot(metric_long, aes(x = value, y = panel, colour = panel)) +
  geom_segment(aes(x = 0.5, xend = value, yend = panel), linewidth = 0.7,
               colour = palette[["light"]]) +
  geom_point(size = 2.2) +
  geom_text(
    aes(label = sprintf("%.3f", value)),
    hjust = -0.25, size = 2.1, colour = palette[["dark"]]
  ) +
  facet_wrap(~ metric, ncol = 1) +
  scale_colour_manual(values = panel_colours) +
  scale_x_continuous(limits = c(0.5, 1.05), breaks = c(0.5, 0.75, 1.0)) +
  labs(
    title = "Untouched Chen held-out comparison",
    x = NULL, y = NULL
  ) +
  theme_pilot() +
  theme(legend.position = "none")

paired_plot <- scores |>
  filter(route_group %in% c("normal", "conventional_adenoma")) |>
  group_by(donor_id, route_group) |>
  summarise(score = median(adaptive_minimum_sufficient), .groups = "drop") |>
  group_by(donor_id) |>
  filter(n_distinct(route_group) == 2) |>
  ungroup() |>
  mutate(
    route_group = factor(
      route_group,
      levels = c("normal", "conventional_adenoma"),
      labels = c("Normal", "Adenoma")
    )
  )
p_d <- ggplot(paired_plot, aes(x = route_group, y = score, group = donor_id)) +
  geom_line(linewidth = 0.45, colour = palette[["grey"]], alpha = 0.75) +
  geom_point(aes(colour = route_group), size = 1.8) +
  scale_colour_manual(values = c("Normal" = palette[["navy"]], "Adenoma" = palette[["red"]])) +
  labs(
    title = "Held-out donor-paired change",
    subtitle = paste0("Adaptive panel; n = ", n_distinct(paired_plot$donor_id), " donors"),
    x = NULL, y = "Adaptive score"
  ) +
  theme_pilot() +
  theme(legend.position = "none")

figure <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(widths = c(1.05, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

base <- file.path(fig_dir, "data_adaptive_panel_pilot")
width_in <- 170 / 25.4
height_in <- 128 / 25.4
svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "sans")
print(figure)
dev.off()
ragg::agg_tiff(
  paste0(base, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(figure)
dev.off()
ragg::agg_png(
  paste0(base, ".png"), width = width_in, height = height_in,
  units = "in", res = 300
)
print(figure)
dev.off()

# -----------------------------------------------------------------------------
# Manifest, QA and concise report
# -----------------------------------------------------------------------------

figure_files <- paste0(base, c(".svg", ".pdf", ".tiff", ".png"))
qa <- data.frame(
  check = c(
    "discovery_validation_separation",
    "both_arms_present",
    "adaptive_panel_smaller_than_100",
    "oof_fidelity_finite",
    "heldout_all_panels_scored",
    "heldout_paired_donors_present",
    "figure_exports_present"
  ),
  pass = c(
    TRUE,
    length(adaptive_up) > 0 && length(adaptive_down) > 0,
    nrow(adaptive_panel) < 100,
    is.finite(oof_total_fidelity),
    nrow(metrics) == 4 && all(is.finite(metrics$auc)),
    all(metrics$n_paired_donors >= 3),
    all(file.exists(figure_files))
  ),
  stringsAsFactors = FALSE
)
write_tsv(qa, "analysis_qa.tsv")

manifest <- list(
  analysis = "data-adaptive minimum-sufficient gene-panel pilot",
  status = "exploratory; manuscript not modified",
  seed = 20260710,
  discovery_only_selection = TRUE,
  discovery_samples = nrow(discovery_agg$expr),
  discovery_donors = length(unique(discovery_agg$meta$donor_id)),
  limma_duplicate_correlation = unname(corfit$consensus),
  eligibility = list(
    limma_bh_fdr = 0.05,
    donor_bootstrap_direction_stability = 0.90,
    donor_bootstrap_ci_excludes_zero = TRUE
  ),
  size_selection = list(
    method = "grouped leave-one-donor-out forward score reconstruction plus one-SE rule",
    maximum_path_length_per_arm = max_k,
    selected_up = selected_k_up,
    selected_down = selected_k_down
  ),
  validation_read_after_panel_fixed = TRUE,
  comparison_panels = c(
    "full stable core", "adaptive minimum-sufficient",
    "100-gene reference", "10-gene biology-guided"
  )
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "analysis_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

adaptive_metric <- metrics[metrics$panel == "Adaptive minimum-sufficient", ]
reference_metric <- metrics[metrics$panel == "100-gene reference", ]
compact_metric <- metrics[metrics$panel == "10-gene biology-guided", ]
preliminary_replace <-
  nrow(adaptive_panel) < 100 &&
  adaptive_metric$spearman_vs_full_core >= 0.90 &&
  adaptive_metric$auc >= reference_metric$auc - 0.02

report <- c(
  "# Data-adaptive gene-panel pilot (v2.6)",
  "",
  "## Scope",
  "",
  "This pilot was run in parallel with the v2.5 submission package. No manuscript, figure or table in v2.5 was overwritten. Gene eligibility, forward reconstruction and panel size were determined from Chen discovery data only. Chen validation data were read after the adaptive panel was fixed.",
  "",
  "## Main results",
  "",
  sprintf("- Stable error-controlled core: %d genes (%d up; %d down).", nrow(core), sum(core$arm == "up"), sum(core$arm == "down")),
  sprintf("- One-SE adaptive panel: %d genes (%d up; %d down).", nrow(adaptive_panel), length(adaptive_up), length(adaptive_down)),
  sprintf("- Discovery grouped out-of-fold fidelity to the full stable core: Spearman rho = %.3f.", oof_total_fidelity),
  sprintf("- Chen held-out AUC: adaptive %.3f; 100-gene reference %.3f; biology-guided 10-gene %.3f.", adaptive_metric$auc, reference_metric$auc, compact_metric$auc),
  sprintf("- Chen held-out adaptive fidelity to the full stable core: Spearman rho = %.3f.", adaptive_metric$spearman_vs_full_core),
  sprintf("- Chen held-out paired positive fraction for the adaptive panel: %.3f (%d donors).", adaptive_metric$paired_positive_fraction, adaptive_metric$n_paired_donors),
  sprintf("- Cluster-bootstrap AUC difference, adaptive minus 100-gene: median %.3f, 95%% CI %.3f to %.3f.", auc_difference$auc_difference_median, auc_difference$ci_low, auc_difference$ci_high),
  "",
  "## Preliminary decision",
  "",
  if (preliminary_replace) {
    "The discovery-only adaptive method passes the prespecified pilot screen for expansion to the five external cohorts. This is not yet sufficient to replace the manuscript panel."
  } else {
    "The discovery-only adaptive method does not pass the prespecified pilot screen. The current manuscript panel should not be replaced on the basis of this result."
  },
  "",
  "## Boundary",
  "",
  "This is a small-sample exploratory feature-selection analysis. The stable core was fixed from the complete discovery dataset before grouped reconstruction, so grouped out-of-fold fidelity is less conservative than a fully nested re-selection pipeline. Replacement would require re-running the complete selector inside each discovery fold and confirming transportability in all external transcriptomic and FFPE resources."
)
writeLines(report, file.path(out_dir, "pilot_report.md"), useBytes = TRUE)

if (!all(qa$pass)) {
  stop(
    "Pilot QA failed: ",
    paste(qa$check[!qa$pass], collapse = ", ")
  )
}

message("Pilot complete")
message("  output: ", out_dir)
message("  figure: ", fig_dir)
message("  adaptive panel: ", length(adaptive_up), " up + ", length(adaptive_down), " down")
message("  held-out AUC: ", sprintf("%.3f", adaptive_metric$auc))
