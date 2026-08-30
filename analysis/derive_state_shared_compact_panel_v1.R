#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260829)

root <- normalizePath(".", mustWork = TRUE)
pseudobulk_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "discovery_pseudobulk",
  "discovery_state_pseudobulk.rds"
)
common_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "common_effects",
  "cross_state_common_effects.tsv.gz"
)
candidate_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "panel_derivation",
  "portable_state_shared_candidate_universe.tsv"
)
out_dir <- file.path(root, "results", "state_aware_program_v1", "panel_derivation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(file.exists(c(
  pseudobulk_path,
  common_path,
  candidate_path
)))) {
  stop("One or more compact-panel inputs are missing")
}

pseudobulk <- readRDS(pseudobulk_path)
common <- read.delim(common_path, check.names = FALSE)
candidates <- read.delim(candidate_path, check.names = FALSE)
metadata <- as.data.frame(pseudobulk$metadata)
states <- c("ABS", "GOB", "TAC")

strict <- common[common$strict_state_shared, , drop = FALSE]
strict_up <- sort(strict$gene[strict$shared_direction == "up"])
strict_down <- sort(strict$gene[strict$shared_direction == "down"])
candidate_up <- sort(candidates$gene[candidates$arm == "up"])
candidate_down <- sort(candidates$gene[candidates$arm == "down"])
if (!all(c(candidate_up, candidate_down) %in% strict$gene)) {
  stop("Portable candidates are not a subset of the strict shared core")
}
if (length(candidate_up) < 2L || length(candidate_down) < 2L) {
  stop("Both candidate arms require at least two genes")
}

strict_genes <- c(strict_up, strict_down)
expression <- matrix(
  NA_real_,
  nrow = nrow(metadata),
  ncol = length(strict_genes),
  dimnames = list(metadata$pseudobulk_id, strict_genes)
)
for (state in states) {
  column_index <- which(metadata$cell_type == state)
  dge <- DGEList(
    counts = as.matrix(
      pseudobulk$counts[strict_genes, column_index, drop = FALSE]
    )
  )
  dge <- calcNormFactors(dge, method = "TMM")
  state_log_cpm <- t(cpm(dge, log = TRUE, prior.count = 0.5))
  expression[column_index, ] <- state_log_cpm[, strict_genes, drop = FALSE]
}
if (any(!is.finite(expression))) {
  stop("Non-finite TMM log-CPM value in compact-panel input")
}

standardize_by_state <- function(train_index, test_index = integer()) {
  z_train <- matrix(
    NA_real_,
    nrow = length(train_index),
    ncol = ncol(expression),
    dimnames = list(rownames(expression)[train_index], colnames(expression))
  )
  z_test <- matrix(
    NA_real_,
    nrow = length(test_index),
    ncol = ncol(expression),
    dimnames = list(rownames(expression)[test_index], colnames(expression))
  )
  for (state in states) {
    train_local <- which(metadata$cell_type[train_index] == state)
    test_local <- which(metadata$cell_type[test_index] == state)
    centre <- colMeans(expression[train_index[train_local], , drop = FALSE])
    scale_value <- apply(
      expression[train_index[train_local], , drop = FALSE],
      2,
      sd
    )
    if (any(!is.finite(scale_value)) || any(scale_value <= 0)) {
      stop("A strict-core gene has zero training variance within ", state)
    }
    z_train[train_local, ] <- sweep(
      sweep(
        expression[train_index[train_local], , drop = FALSE],
        2,
        centre,
        "-"
      ),
      2,
      scale_value,
      "/"
    )
    if (length(test_local)) {
      z_test[test_local, ] <- sweep(
        sweep(
          expression[test_index[test_local], , drop = FALSE],
          2,
          centre,
          "-"
        ),
        2,
        scale_value,
        "/"
      )
    }
  }
  list(train = z_train, test = z_test)
}

full_score <- function(z) {
  rowMeans(z[, strict_up, drop = FALSE]) -
    rowMeans(z[, strict_down, drop = FALSE])
}

safe_cor <- function(x, y, method = "pearson") {
  value <- suppressWarnings(cor(
    x,
    y,
    method = method,
    use = "pairwise.complete.obs"
  ))
  ifelse(is.finite(value), value, NA_real_)
}

correlate_columns <- function(matrix, target) {
  matrix <- sweep(matrix, 2, colMeans(matrix), "-")
  target <- target - mean(target)
  numerator <- as.numeric(crossprod(target, matrix))
  denominator <- sqrt(sum(target^2) * colSums(matrix^2))
  output <- numerator / denominator
  output[!is.finite(output)] <- -Inf
  output
}

greedy_balanced_pair_path <- function(z, target, max_steps) {
  x_up <- z[, candidate_up, drop = FALSE]
  x_down <- z[, candidate_down, drop = FALSE]
  selected_up <- character()
  selected_down <- character()
  base_score <- rep(0, nrow(z))
  path <- vector("list", max_steps)
  for (step in seq_len(max_steps)) {
    remaining_up <- setdiff(candidate_up, selected_up)
    remaining_down <- setdiff(candidate_down, selected_down)
    pair_grid <- expand.grid(
      up_gene = remaining_up,
      down_gene = remaining_down,
      stringsAsFactors = FALSE
    )
    pair_grid <- pair_grid[
      order(pair_grid$up_gene, pair_grid$down_gene),
      ,
      drop = FALSE
    ]
    candidate_scores <-
      matrix(base_score, nrow = nrow(z), ncol = nrow(pair_grid)) +
      x_up[, pair_grid$up_gene, drop = FALSE] -
      x_down[, pair_grid$down_gene, drop = FALSE]
    correlations <- correlate_columns(candidate_scores, target)
    best <- which.max(correlations)
    chosen_up <- pair_grid$up_gene[best]
    chosen_down <- pair_grid$down_gene[best]
    selected_up <- c(selected_up, chosen_up)
    selected_down <- c(selected_down, chosen_down)
    base_score <- base_score + x_up[, chosen_up] - x_down[, chosen_down]
    path[[step]] <- data.frame(
      pair_step = step,
      up_gene = chosen_up,
      down_gene = chosen_down,
      training_correlation = safe_cor(base_score / step, target),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, path)
}

panel_score_from_path <- function(z, path, pairs) {
  rowMeans(z[, path$up_gene[seq_len(pairs)], drop = FALSE]) -
    rowMeans(z[, path$down_gene[seq_len(pairs)], drop = FALSE])
}

donors <- sort(unique(as.character(metadata$donor_id)))
maximum_pairs <- min(
  length(candidate_up),
  length(candidate_down),
  length(donors) - 2L
)
if (maximum_pairs < 2L) {
  stop("Too few candidate pairs for an objective fidelity curve")
}

message(
  "Running ",
  length(donors),
  " donor-held-out derivations across ",
  maximum_pairs,
  " balanced pair steps"
)
oof <- data.frame(
  pseudobulk_id = metadata$pseudobulk_id,
  donor_id = as.character(metadata$donor_id),
  route = as.character(metadata$route),
  cell_type = as.character(metadata$cell_type),
  full_programme_score = NA_real_,
  stringsAsFactors = FALSE
)
for (pairs in seq_len(maximum_pairs)) {
  oof[[paste0("panel_pairs_", pairs)]] <- NA_real_
}
fold_paths <- vector("list", length(donors))

for (fold in seq_along(donors)) {
  held_donor <- donors[fold]
  test_index <- which(metadata$donor_id == held_donor)
  train_index <- setdiff(seq_len(nrow(metadata)), test_index)
  scaled <- standardize_by_state(train_index, test_index)
  train_target <- full_score(scaled$train)
  test_target <- full_score(scaled$test)
  path <- greedy_balanced_pair_path(
    scaled$train,
    train_target,
    maximum_pairs
  )
  path$held_out_donor <- held_donor
  fold_paths[[fold]] <- path
  oof$full_programme_score[test_index] <- test_target
  for (pairs in seq_len(maximum_pairs)) {
    oof[[paste0("panel_pairs_", pairs)]][test_index] <-
      panel_score_from_path(scaled$test, path, pairs)
  }
  if (fold %% 5L == 0L || fold == length(donors)) {
    message("  completed fold ", fold, "/", length(donors))
  }
}
fold_paths <- do.call(rbind, fold_paths)

n_bootstrap <- 2000L
point_fidelity <- vapply(
  seq_len(maximum_pairs),
  function(pairs) safe_cor(
    oof[[paste0("panel_pairs_", pairs)]],
    oof$full_programme_score,
    method = "spearman"
  ),
  numeric(1)
)
bootstrap_fidelity <- matrix(
  NA_real_,
  nrow = n_bootstrap,
  ncol = maximum_pairs
)
for (bootstrap in seq_len(n_bootstrap)) {
  sampled_donors <- sample(donors, length(donors), replace = TRUE)
  index <- unlist(
    lapply(sampled_donors, function(donor) which(oof$donor_id == donor)),
    use.names = FALSE
  )
  for (pairs in seq_len(maximum_pairs)) {
    bootstrap_fidelity[bootstrap, pairs] <- safe_cor(
      oof[[paste0("panel_pairs_", pairs)]][index],
      oof$full_programme_score[index],
      method = "spearman"
    )
  }
}
fidelity_se <- apply(bootstrap_fidelity, 2, sd, na.rm = TRUE)
fidelity_ci_low <- apply(
  bootstrap_fidelity,
  2,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)
fidelity_ci_high <- apply(
  bootstrap_fidelity,
  2,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)
monotone_fidelity <- cummax(point_fidelity)
normalized_x <- (seq_len(maximum_pairs) - 1) / (maximum_pairs - 1)
fidelity_range <- max(monotone_fidelity) - min(monotone_fidelity)
if (!is.finite(fidelity_range) || fidelity_range <= 0) {
  stop("The out-of-fold fidelity curve has no finite range")
}
normalized_y <-
  (monotone_fidelity - min(monotone_fidelity)) / fidelity_range
knee_distance <- normalized_y - normalized_x
selected_pairs <- which.max(knee_distance)
maximum_index <- which.max(point_fidelity)
one_se_threshold <-
  point_fidelity[maximum_index] - fidelity_se[maximum_index]
one_se_pairs <- min(which(point_fidelity >= one_se_threshold))

fidelity_curve <- data.frame(
  pairs = seq_len(maximum_pairs),
  total_genes = 2L * seq_len(maximum_pairs),
  oof_spearman = point_fidelity,
  donor_bootstrap_se = fidelity_se,
  ci_low = fidelity_ci_low,
  ci_high = fidelity_ci_high,
  monotone_oof_spearman = monotone_fidelity,
  normalized_knee_distance = knee_distance,
  kneedle_selected_pairs = selected_pairs,
  one_se_pairs = one_se_pairs,
  one_se_threshold = one_se_threshold,
  stringsAsFactors = FALSE
)

state_fidelity <- do.call(rbind, lapply(states, function(state) {
  index <- which(oof$cell_type == state)
  data.frame(
    cell_type = state,
    pairs = seq_len(maximum_pairs),
    total_genes = 2L * seq_len(maximum_pairs),
    oof_spearman = vapply(
      seq_len(maximum_pairs),
      function(pairs) safe_cor(
        oof[[paste0("panel_pairs_", pairs)]][index],
        oof$full_programme_score[index],
        method = "spearman"
      ),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}))

message(
  "Kneedle selected ",
  selected_pairs,
  " pairs; OOF rho=",
  sprintf("%.3f", point_fidelity[selected_pairs]),
  "; one-SE sensitivity=",
  one_se_pairs,
  " pairs"
)

all_scaled <- standardize_by_state(seq_len(nrow(metadata)))$train
all_target <- full_score(all_scaled)
final_path <- greedy_balanced_pair_path(
  all_scaled,
  all_target,
  maximum_pairs
)

outer_frequency <- rbind(
  data.frame(
    gene = candidate_up,
    arm = "up",
    selected_folds = vapply(
      candidate_up,
      function(gene) sum(
        fold_paths$pair_step <= selected_pairs & fold_paths$up_gene == gene
      ),
      integer(1)
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    gene = candidate_down,
    arm = "down",
    selected_folds = vapply(
      candidate_down,
      function(gene) sum(
        fold_paths$pair_step <= selected_pairs & fold_paths$down_gene == gene
      ),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
)
outer_frequency$selection_frequency <-
  outer_frequency$selected_folds / length(donors)

message("Running 2,000 whole-donor bootstrap path refits at the selected size")
bootstrap_refits <- parallel::mclapply(
  seq_len(n_bootstrap),
  function(bootstrap) {
    sampled_donors <- sample(donors, length(donors), replace = TRUE)
    index <- unlist(
      lapply(sampled_donors, function(donor) which(metadata$donor_id == donor)),
      use.names = FALSE
    )
    z <- matrix(
      NA_real_,
      nrow = length(index),
      ncol = ncol(expression),
      dimnames = list(NULL, colnames(expression))
    )
    for (state in states) {
      local <- which(metadata$cell_type[index] == state)
      values <- expression[index[local], , drop = FALSE]
      centre <- colMeans(values)
      scale_value <- apply(values, 2, sd)
      if (any(!is.finite(scale_value)) || any(scale_value <= 0)) {
        return(NULL)
      }
      z[local, ] <- sweep(sweep(values, 2, centre, "-"), 2, scale_value, "/")
    }
    path <- greedy_balanced_pair_path(z, full_score(z), selected_pairs)
    data.frame(
      bootstrap = bootstrap,
      pair_step = rep(path$pair_step, 2L),
      gene = c(path$up_gene, path$down_gene),
      arm = rep(c("up", "down"), each = nrow(path)),
      stringsAsFactors = FALSE
    )
  },
  mc.cores = 4L,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE
)
valid_bootstrap <- !vapply(bootstrap_refits, is.null, logical(1))
if (sum(valid_bootstrap) < 0.95 * n_bootstrap) {
  stop("Fewer than 95% of donor bootstrap path refits were valid")
}
bootstrap_selections <- do.call(rbind, bootstrap_refits[valid_bootstrap])
bootstrap_frequency <- aggregate(
  bootstrap ~ gene + arm,
  data = unique(bootstrap_selections[c("bootstrap", "gene", "arm")]),
  FUN = length
)
colnames(bootstrap_frequency)[3L] <- "selected_bootstraps"
bootstrap_frequency$selection_frequency <-
  bootstrap_frequency$selected_bootstraps / sum(valid_bootstrap)

panel <- rbind(
  data.frame(
    gene = final_path$up_gene[seq_len(selected_pairs)],
    arm = "up",
    pair_step = seq_len(selected_pairs),
    route_weight = 1,
    stringsAsFactors = FALSE
  ),
  data.frame(
    gene = final_path$down_gene[seq_len(selected_pairs)],
    arm = "down",
    pair_step = seq_len(selected_pairs),
    route_weight = -1,
    stringsAsFactors = FALSE
  )
)
panel <- merge(
  panel,
  final_path[c("pair_step", "training_correlation")],
  by = "pair_step",
  all.x = TRUE,
  sort = FALSE
)
panel <- merge(
  panel,
  outer_frequency,
  by = c("gene", "arm"),
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_outer")
)
colnames(panel)[colnames(panel) == "selected_folds"] <-
  "selected_donor_heldout_folds"
colnames(panel)[colnames(panel) == "selection_frequency"] <-
  "donor_heldout_selection_frequency"
panel <- merge(
  panel,
  bootstrap_frequency,
  by = c("gene", "arm"),
  all.x = TRUE,
  sort = FALSE
)
colnames(panel)[colnames(panel) == "selection_frequency"] <-
  "donor_bootstrap_selection_frequency"
panel$selected_bootstraps[is.na(panel$selected_bootstraps)] <- 0L
panel$donor_bootstrap_selection_frequency[
  is.na(panel$donor_bootstrap_selection_frequency)
] <- 0
panel <- merge(
  panel,
  candidates[c(
    "gene",
    "common_effect",
    "common_z",
    "common_q_value",
    "max_lfsr",
    "all_six_platforms_present",
    "protein_coding"
  )],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
panel$panel_size_rule <-
  "Kneedle on monotone donor-grouped out-of-fold fidelity curve"
panel$validation_outcomes_used <- FALSE
panel$panel_frozen_before_validation <- TRUE
panel <- panel[order(panel$pair_step, -panel$route_weight, panel$gene), ]

write.table(
  oof,
  gzfile(file.path(out_dir, "discovery_grouped_oof_scores.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  fidelity_curve,
  file.path(out_dir, "discovery_grouped_oof_fidelity_curve.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  state_fidelity,
  file.path(out_dir, "discovery_state_specific_fidelity_curve.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  fold_paths,
  gzfile(file.path(out_dir, "discovery_donor_heldout_pair_paths.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  final_path,
  file.path(out_dir, "discovery_full_balanced_pair_path.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  bootstrap_selections,
  gzfile(file.path(out_dir, "donor_bootstrap_selected_genes_2000.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  merge(
    outer_frequency,
    bootstrap_frequency,
    by = c("gene", "arm"),
    all = TRUE,
    suffixes = c("_donor_heldout", "_donor_bootstrap")
  ),
  file.path(out_dir, "candidate_selection_stability.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  panel,
  file.path(out_dir, "compact_state_shared_panel_frozen.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
summary <- data.frame(
  metric = c(
    "strict_core_total",
    "strict_core_up",
    "strict_core_down",
    "portable_candidate_up",
    "portable_candidate_down",
    "discovery_donors",
    "discovery_profiles",
    "maximum_pairs_evaluated",
    "selected_pairs",
    "selected_total_genes",
    "selected_oof_spearman",
    "selected_oof_ci_low",
    "selected_oof_ci_high",
    "one_se_pairs",
    "valid_gene_selection_bootstraps",
    "minimum_panel_bootstrap_selection_frequency"
  ),
  value = c(
    nrow(strict),
    length(strict_up),
    length(strict_down),
    length(candidate_up),
    length(candidate_down),
    length(donors),
    nrow(metadata),
    maximum_pairs,
    selected_pairs,
    2L * selected_pairs,
    point_fidelity[selected_pairs],
    fidelity_ci_low[selected_pairs],
    fidelity_ci_high[selected_pairs],
    one_se_pairs,
    sum(valid_bootstrap),
    min(panel$donor_bootstrap_selection_frequency)
  ),
  stringsAsFactors = FALSE
)
write.table(
  summary,
  file.path(out_dir, "compact_panel_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = "derive_state_shared_compact_panel_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_hashes = as.list(setNames(
    vapply(
      c(pseudobulk_path, common_path, candidate_path),
      digest::digest,
      character(1),
      algo = "sha256",
      file = TRUE
    ),
    basename(c(pseudobulk_path, common_path, candidate_path))
  )),
  target = "equal-arm strict state-shared programme score",
  normalisation = "state-specific TMM log-CPM",
  cross_validation = "leave complete donor out; state-specific training centring/scaling",
  path = "one frozen-up plus one frozen-down gene per step",
  size_rule = "Kneedle maximum distance on monotone OOF Spearman fidelity curve",
  sensitivity_size_rule = "one-standard-error rule",
  donor_bootstraps = n_bootstrap,
  selected_pairs = selected_pairs,
  selected_genes = nrow(panel),
  validation_expression_read = FALSE,
  validation_outcomes_used = FALSE,
  random_seed = 20260829,
  package_versions = as.list(vapply(
    c("edgeR", "Matrix", "digest", "jsonlite"),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "compact_panel_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Frozen compact readout: ",
  nrow(panel),
  " genes (",
  selected_pairs,
  " balanced pairs)"
)
