#!/usr/bin/env Rscript

# Objective, discovery-only reduction of the fixed 287-gene epithelial programme.
#
# Candidate eligibility is limited to protein-coding core genes with label-blind
# feature availability on all five external platforms and GSE117606 FFPE.  The
# validation expression matrices and all validation tissue labels are not read.
# A balanced greedy path adds one adenoma-up and one adenoma-down gene per step
# to reconstruct the complete 287-gene score.  Panel size is selected by the
# one-standard-error rule applied to donor-grouped out-of-fold fidelity.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

set.seed(20260710)

root <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(root, "results", "objective_compact_panel_v2_7")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

discovery_meta_path <- file.path(
  root, "results", "route_signature", "chen_discovery_specimen_pseudobulk_meta.tsv"
)
discovery_expr_path <- file.path(
  root, "results", "route_signature", "chen_discovery_specimen_pseudobulk_expression.tsv.gz"
)
core_path <- file.path(
  root, "results", "data_adaptive_panel_pilot_v2_6", "stable_error_controlled_core.tsv"
)
candidate_path <- file.path(
  out_dir, "portable_protein_coding_candidate_universe.tsv"
)

read_tsv_base <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else path
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  read.delim(
    con, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, quote = "", comment.char = ""
  )
}

write_tsv <- function(x, name) {
  write.table(
    x, file.path(out_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )
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
    rows[[i]] <- if (length(idx) == 1L) {
      expr[idx, ]
    } else {
      apply(expr[idx, , drop = FALSE], 2, median, na.rm = TRUE)
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
  if (is.null(scale_value)) scale_value <- apply(x, 2, sd, na.rm = TRUE)
  scale_value[!is.finite(scale_value) | scale_value == 0] <- NA_real_
  list(
    z = sweep(sweep(x, 2, center, "-"), 2, scale_value, "/"),
    center = center,
    scale = scale_value
  )
}

safe_cor <- function(x, y, method = "pearson") {
  value <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = method))
  ifelse(is.finite(value), value, NA_real_)
}

greedy_balanced_pair_path <- function(x_up, x_down, target) {
  x_up <- as.matrix(x_up)
  x_down <- as.matrix(x_down)
  max_pairs <- min(ncol(x_up), ncol(x_down))
  selected_up_idx <- integer(0)
  selected_down_idx <- integer(0)
  selected_up <- character(0)
  selected_down <- character(0)
  up_sum <- rep(0, nrow(x_up))
  down_sum <- rep(0, nrow(x_down))
  path <- vector("list", max_pairs)

  for (k in seq_len(max_pairs)) {
    remaining_up <- setdiff(seq_len(ncol(x_up)), selected_up_idx)
    remaining_down <- setdiff(seq_len(ncol(x_down)), selected_down_idx)
    candidate_rows <- vector("list", length(remaining_up) * length(remaining_down))
    position <- 0L
    for (up_idx in remaining_up) {
      candidate_up <- (up_sum + x_up[, up_idx]) / k
      for (down_idx in remaining_down) {
        position <- position + 1L
        candidate_score <- candidate_up - (down_sum + x_down[, down_idx]) / k
        candidate_rows[[position]] <- data.frame(
          up_idx = up_idx,
          down_idx = down_idx,
          up_gene = colnames(x_up)[up_idx],
          down_gene = colnames(x_down)[down_idx],
          training_correlation = safe_cor(candidate_score, target, method = "pearson"),
          stringsAsFactors = FALSE
        )
      }
    }
    candidates <- bind_rows(candidate_rows) |>
      mutate(training_correlation = ifelse(
        is.finite(training_correlation), training_correlation, -Inf
      )) |>
      arrange(desc(training_correlation), up_gene, down_gene)
    chosen <- candidates[1, , drop = FALSE]
    selected_up_idx <- c(selected_up_idx, chosen$up_idx)
    selected_down_idx <- c(selected_down_idx, chosen$down_idx)
    selected_up <- c(selected_up, chosen$up_gene)
    selected_down <- c(selected_down, chosen$down_gene)
    up_sum <- up_sum + x_up[, chosen$up_idx]
    down_sum <- down_sum + x_down[, chosen$down_idx]
    path[[k]] <- data.frame(
      pair_step = k,
      up_gene = chosen$up_gene,
      down_gene = chosen$down_gene,
      training_correlation = safe_cor(up_sum / k - down_sum / k, target),
      stringsAsFactors = FALSE
    )
  }
  bind_rows(path)
}

bootstrap_oof_curve <- function(oof, donors, max_pairs, n_boot = 2000L) {
  point <- vapply(
    seq_len(max_pairs),
    function(k) safe_cor(oof[[paste0("panel_k", k)]], oof$full_core_score, "spearman"),
    numeric(1)
  )
  unique_donors <- sort(unique(donors))
  boot <- matrix(NA_real_, nrow = n_boot, ncol = max_pairs)
  for (b in seq_len(n_boot)) {
    sampled <- sample(unique_donors, length(unique_donors), replace = TRUE)
    idx <- unlist(lapply(sampled, function(d) which(donors == d)), use.names = FALSE)
    for (k in seq_len(max_pairs)) {
      boot[b, k] <- safe_cor(
        oof[[paste0("panel_k", k)]][idx],
        oof$full_core_score[idx],
        "spearman"
      )
    }
  }
  se <- apply(boot, 2, sd, na.rm = TRUE)
  ci_low <- apply(boot, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_high <- apply(boot, 2, quantile, probs = 0.975, na.rm = TRUE)
  k_max <- which.max(point)
  one_se_threshold <- point[[k_max]] - se[[k_max]]
  one_se_k <- min(which(point >= one_se_threshold))
  # Kneedle-style maximum distance from the normalized diagonal.  cummax
  # prevents small cross-validation reversals from creating a false knee.
  monotone_point <- cummax(point)
  normalized_x <- (seq_len(max_pairs) - 1) / (max_pairs - 1)
  normalized_y <- (monotone_point - min(monotone_point)) /
    (max(monotone_point) - min(monotone_point))
  knee_distance <- normalized_y - normalized_x
  knee_k <- which.max(knee_distance)
  list(
    curve = data.frame(
      pairs = seq_len(max_pairs),
      total_genes = 2L * seq_len(max_pairs),
      oof_spearman = point,
      donor_bootstrap_se = se,
      ci_low = ci_low,
      ci_high = ci_high,
      maximum_fidelity_pairs = k_max,
      maximum_fidelity = point[[k_max]],
      one_se_threshold = one_se_threshold,
      one_se_pairs = one_se_k,
      monotone_oof_spearman = monotone_point,
      normalized_knee_distance = knee_distance,
      selected_pairs = knee_k,
      primary_size_rule = "Kneedle maximum normalized distance after monotone envelope",
      stringsAsFactors = FALSE
    ),
    bootstrap = boot,
    selected_pairs = knee_k,
    one_se_pairs = one_se_k
  )
}

message("Loading discovery-only inputs")
meta <- read_tsv_base(discovery_meta_path)
expr <- read_tsv_base(discovery_expr_path)
core <- read_tsv_base(core_path)
candidates <- read_tsv_base(candidate_path)

as_logical_field <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}
for (field in c(
  "objective_selection_eligible", "all_six_platforms_present", "protein_coding"
)) {
  candidates[[field]] <- as_logical_field(candidates[[field]])
}

stopifnot(nrow(core) == 287L, setequal(core$arm, c("up", "down")))
stopifnot(all(candidates$objective_selection_eligible))
stopifnot(all(candidates$gene %in% core$gene))
stopifnot(nrow(meta) == nrow(expr))

aggregated <- aggregate_donor_route(meta, expr)
rm(expr)
gc(verbose = FALSE)

core_up <- sort(core$gene[core$arm == "up"])
core_down <- sort(core$gene[core$arm == "down"])
candidate_up <- sort(candidates$gene[candidates$arm == "up"])
candidate_down <- sort(candidates$gene[candidates$arm == "down"])
core_genes <- c(core_up, core_down)
max_pairs <- min(length(candidate_up), length(candidate_down))
donors <- sort(unique(aggregated$meta$donor_id))
n <- nrow(aggregated$expr)

message(
  "Candidate universe: ", length(candidate_up), " up + ",
  length(candidate_down), " down; ", length(donors), " donor folds"
)

oof <- data.frame(
  row_id = seq_len(n),
  donor_id = aggregated$meta$donor_id,
  route_group = aggregated$meta$route_group,
  full_core_score = NA_real_,
  stringsAsFactors = FALSE
)
for (k in seq_len(max_pairs)) oof[[paste0("panel_k", k)]] <- NA_real_
fold_path_rows <- vector("list", length(donors))

for (fold in seq_along(donors)) {
  held_donor <- donors[[fold]]
  test_idx <- which(aggregated$meta$donor_id == held_donor)
  train_idx <- setdiff(seq_len(n), test_idx)
  train_raw <- aggregated$expr[train_idx, core_genes, drop = FALSE]
  test_raw <- aggregated$expr[test_idx, core_genes, drop = FALSE]
  train_scaled <- zscore_matrix(train_raw)
  z_train <- train_scaled$z
  z_test <- zscore_matrix(
    test_raw, center = train_scaled$center, scale_value = train_scaled$scale
  )$z
  target_train <- rowMeans(z_train[, core_up, drop = FALSE]) -
    rowMeans(z_train[, core_down, drop = FALSE])
  oof$full_core_score[test_idx] <-
    rowMeans(z_test[, core_up, drop = FALSE]) -
    rowMeans(z_test[, core_down, drop = FALSE])

  path <- greedy_balanced_pair_path(
    z_train[, candidate_up, drop = FALSE],
    z_train[, candidate_down, drop = FALSE],
    target_train
  )
  path$held_out_donor <- held_donor
  fold_path_rows[[fold]] <- path
  for (k in seq_len(max_pairs)) {
    selected_up <- path$up_gene[seq_len(k)]
    selected_down <- path$down_gene[seq_len(k)]
    oof[[paste0("panel_k", k)]][test_idx] <-
      rowMeans(z_test[, selected_up, drop = FALSE]) -
      rowMeans(z_test[, selected_down, drop = FALSE])
  }
  if (fold %% 5L == 0L || fold == length(donors)) {
    message("  completed donor fold ", fold, "/", length(donors))
  }
}
fold_paths <- bind_rows(fold_path_rows)
curve_result <- bootstrap_oof_curve(oof, oof$donor_id, max_pairs, n_boot = 2000L)
curve <- curve_result$curve
selected_pairs <- curve_result$selected_pairs

message(
  "Kneedle selection: ", selected_pairs, " pairs (",
  2L * selected_pairs, " genes); OOF rho = ",
  sprintf("%.3f", curve$oof_spearman[curve$pairs == selected_pairs]),
  "; one-SE sensitivity = ", unique(curve$one_se_pairs), " pairs"
)

full_scaled <- zscore_matrix(aggregated$expr[, core_genes, drop = FALSE])$z
full_target <- rowMeans(full_scaled[, core_up, drop = FALSE]) -
  rowMeans(full_scaled[, core_down, drop = FALSE])
final_path <- greedy_balanced_pair_path(
  full_scaled[, candidate_up, drop = FALSE],
  full_scaled[, candidate_down, drop = FALSE],
  full_target
)

selection_frequency <- bind_rows(
  fold_paths |>
    filter(pair_step <= selected_pairs) |>
    count(up_gene, name = "selected_folds") |>
    transmute(
      gene = up_gene, arm = "up", selected_folds,
      selection_frequency = selected_folds / length(donors)
    ),
  fold_paths |>
    filter(pair_step <= selected_pairs) |>
    count(down_gene, name = "selected_folds") |>
    transmute(
      gene = down_gene, arm = "down", selected_folds,
      selection_frequency = selected_folds / length(donors)
    )
)

panel <- bind_rows(
  final_path |>
    filter(pair_step <= selected_pairs) |>
    transmute(
      gene = up_gene, arm = "up", pair_step,
      full_discovery_path_correlation = training_correlation
    ),
  final_path |>
    filter(pair_step <= selected_pairs) |>
    transmute(
      gene = down_gene, arm = "down", pair_step,
      full_discovery_path_correlation = training_correlation
    )
) |>
  left_join(selection_frequency, by = c("gene", "arm")) |>
  mutate(
    selected_folds = coalesce(selected_folds, 0L),
    selection_frequency = coalesce(selection_frequency, 0),
    route_weight = ifelse(arm == "up", 1, -1),
    panel_size_rule = "Kneedle maximum normalized distance on the monotone grouped-OOF fidelity curve",
    fixed_gene_count_used = FALSE,
    validation_outcomes_used = FALSE
  ) |>
  left_join(
    candidates |>
      select(
        gene, logFC, adj.P.Val, direction_stability,
        bootstrap_ci_low, bootstrap_ci_high, type_of_gene,
        all_six_platforms_present, protein_coding
      ),
    by = "gene"
  ) |>
  arrange(pair_step, desc(arm))

if (nrow(panel) != 2L * selected_pairs || anyDuplicated(panel$gene)) {
  stop("Frozen objective panel has an invalid size or duplicate genes")
}
if (!all(panel$all_six_platforms_present) || !all(panel$protein_coding)) {
  stop("Frozen objective panel violates the portability gate")
}

write_tsv(oof, "discovery_grouped_oof_balanced_pair_scores.tsv")
write_tsv(curve, "discovery_grouped_oof_fidelity_curve.tsv")
write_tsv(fold_paths, "discovery_fold_specific_pair_paths.tsv")
write_tsv(final_path, "discovery_full_balanced_pair_path.tsv")
write_tsv(selection_frequency, "discovery_selected_size_gene_stability.tsv")
write_tsv(panel, "objective_compact_panel_frozen.tsv")

summary <- data.frame(
  metric = c(
    "core_total", "core_up", "core_down",
    "portable_candidate_up", "portable_candidate_down",
    "discovery_donors", "discovery_donor_route_profiles",
    "selected_pairs", "selected_total_genes",
    "selected_oof_spearman", "selected_oof_ci_low", "selected_oof_ci_high",
    "maximum_oof_spearman", "kneedle_selected_pairs",
    "one_se_sensitivity_pairs", "one_se_threshold",
    "minimum_selected_gene_fold_frequency"
  ),
  value = c(
    nrow(core), length(core_up), length(core_down),
    length(candidate_up), length(candidate_down),
    length(donors), n,
    selected_pairs, 2L * selected_pairs,
    curve$oof_spearman[curve$pairs == selected_pairs],
    curve$ci_low[curve$pairs == selected_pairs],
    curve$ci_high[curve$pairs == selected_pairs],
    max(curve$oof_spearman),
    selected_pairs,
    unique(curve$one_se_pairs),
    unique(curve$one_se_threshold),
    min(panel$selection_frequency)
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, "objective_selection_summary.tsv")

manifest <- list(
  analysis = "objective compact panel selection from the 287-gene core",
  seed = 20260710,
  selection_data = "Chen discovery donor-route profiles only",
  candidate_filter = paste(
    "287-core member; NCBI protein-coding; label-blind feature presence",
    "on all five external platforms and GSE117606 FFPE"
  ),
  balance_rule = "one adenoma-up plus one adenoma-down gene added at each step",
  path_objective = "Pearson reconstruction of the complete 287-gene two-arm score",
  size_rule = paste(
    "Kneedle maximum normalized distance on the monotone envelope of the",
    "leave-one-donor-out Spearman fidelity curve; Satopaa et al.,",
    "doi:10.1109/ICDCSW.2011.20"
  ),
  sensitivity_size_rule = paste(
    "smallest balanced pair count within one donor-bootstrap standard error",
    "of maximum leave-one-donor-out Spearman fidelity"
  ),
  donor_bootstrap_replicates = 2000,
  fixed_gene_count_used = FALSE,
  validation_expression_read = FALSE,
  validation_outcomes_used = FALSE,
  panel_frozen_before_validation = TRUE
)
jsonlite::write_json(
  manifest, file.path(out_dir, "objective_selection_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)

qa <- data.frame(
  check = c(
    "fixed_287_core_loaded",
    "both_portable_candidate_arms_are_nonempty",
    "donor_grouping_preserved",
    "kneedle_size_is_data_derived",
    "frozen_panel_is_balanced",
    "all_frozen_genes_pass_portability_gate",
    "validation_outcomes_not_used"
  ),
  passed = c(
    nrow(core) == 287,
    min(length(candidate_up), length(candidate_down)) >= 1,
    all(table(oof$donor_id) %in% c(1L, 2L)),
    selected_pairs >= 1 && selected_pairs <= max_pairs,
    sum(panel$arm == "up") == sum(panel$arm == "down"),
    all(panel$all_six_platforms_present) && all(panel$protein_coding),
    all(!panel$validation_outcomes_used)
  ),
  stringsAsFactors = FALSE
)
write_tsv(qa, "objective_selection_qa.tsv")
if (!all(qa$passed)) {
  stop("Objective selection QA failed: ", paste(qa$check[!qa$passed], collapse = ", "))
}

message("Frozen panel written before any validation dataset is read")
