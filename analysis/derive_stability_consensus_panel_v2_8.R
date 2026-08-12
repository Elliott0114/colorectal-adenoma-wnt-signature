#!/usr/bin/env Rscript

# Derive a compact consensus panel without fixing its gene count.
#
# The number of genes considered in each donor-held-out fold comes from the
# previously frozen Kneedle point (six balanced pairs).  A gene is retained in
# the consensus only if it appears within that fold-specific path in a strict
# majority (>50%) of the 27 leave-one-donor-out derivations.  Thus neither the
# final gene count nor the up/down arm counts are prescribed.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

root <- normalizePath(getwd(), mustWork = TRUE)
source_dir <- file.path(root, "results", "objective_compact_panel_v2_7")
out_dir <- file.path(root, "results", "stability_consensus_panel_v2_8")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(
    path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE, quote = "", comment.char = ""
  )
}

write_tsv <- function(frame, filename) {
  write.table(
    frame, file.path(out_dir, filename), sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )
}

fold_paths <- read_tsv(file.path(source_dir, "discovery_fold_specific_pair_paths.tsv"))
curve <- read_tsv(file.path(source_dir, "discovery_grouped_oof_fidelity_curve.tsv"))
candidates <- read_tsv(file.path(source_dir, "portable_protein_coding_candidate_universe.tsv"))
as_logical_field <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}
for (field in c("all_six_platforms_present", "protein_coding")) {
  candidates[[field]] <- as_logical_field(candidates[[field]])
}

selected_pairs <- unique(curve$selected_pairs)
if (length(selected_pairs) != 1L || !is.finite(selected_pairs)) {
  stop("The frozen discovery Kneedle size is not uniquely defined")
}
donors <- sort(unique(fold_paths$held_out_donor))
n_folds <- length(donors)
minimum_folds <- floor(n_folds / 2) + 1L

frequency <- bind_rows(
  fold_paths |>
    filter(pair_step <= selected_pairs) |>
    distinct(held_out_donor, gene = up_gene) |>
    count(gene, name = "selected_folds") |>
    mutate(arm = "up"),
  fold_paths |>
    filter(pair_step <= selected_pairs) |>
    distinct(held_out_donor, gene = down_gene) |>
    count(gene, name = "selected_folds") |>
    mutate(arm = "down")
) |>
  mutate(
    total_donor_folds = n_folds,
    selection_frequency = selected_folds / total_donor_folds,
    strict_majority_threshold_folds = minimum_folds,
    retained_by_strict_majority = selected_folds >= minimum_folds
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
  arrange(arm, desc(selected_folds), gene)

panel <- frequency |>
  filter(retained_by_strict_majority) |>
  mutate(
    route_weight = ifelse(arm == "up", 1, -1),
    selection_rule = paste0(
      "Selected within the Kneedle-sized fold path in a strict majority of ",
      n_folds, " donor-held-out derivations (at least ", minimum_folds, ")"
    ),
    fixed_gene_count_used = FALSE,
    validation_outcomes_used = FALSE,
    panel_frozen_before_validation = TRUE
  ) |>
  arrange(desc(selection_frequency), arm, gene)

if (nrow(panel) < 6L || !all(c("up", "down") %in% panel$arm)) {
  stop("Strict-majority consensus produced an unusably small or one-sided panel")
}
if (anyDuplicated(panel$gene)) {
  stop("Consensus panel contains duplicate genes")
}
if (!all(panel$all_six_platforms_present) || !all(panel$protein_coding)) {
  stop("Consensus panel violates the frozen portability gate")
}

write_tsv(frequency, "stability_selection_frequency_all_candidates.tsv")
write_tsv(panel, "stability_consensus_panel_frozen.tsv")

summary <- data.frame(
  metric = c(
    "donor_held_out_derivations", "kneedle_pairs_per_derivation",
    "genes_considered_per_derivation", "strict_majority_minimum_folds",
    "consensus_total_genes", "consensus_up_genes", "consensus_down_genes",
    "minimum_retained_selection_frequency"
  ),
  value = c(
    n_folds, selected_pairs, 2L * selected_pairs, minimum_folds,
    nrow(panel), sum(panel$arm == "up"), sum(panel$arm == "down"),
    min(panel$selection_frequency)
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, "stability_consensus_summary.tsv")

qa <- data.frame(
  check = c(
    "kneedle_size_reused_without_refitting",
    "strict_majority_is_data_independent_threshold",
    "gene_count_not_fixed",
    "both_biological_arms_retained",
    "all_genes_pass_portability_gate",
    "validation_outcomes_not_used"
  ),
  passed = c(
    selected_pairs == 6L,
    minimum_folds == 14L,
    !any(panel$fixed_gene_count_used),
    all(c("up", "down") %in% panel$arm),
    all(panel$all_six_platforms_present) && all(panel$protein_coding),
    all(!panel$validation_outcomes_used)
  ),
  stringsAsFactors = FALSE
)
write_tsv(qa, "stability_consensus_qa.tsv")
if (!all(qa$passed)) {
  stop("Stability-consensus QA failed")
}

manifest <- list(
  analysis = "strict-majority stability consensus from donor-held-out paths",
  source_kneedle_pairs = selected_pairs,
  source_donor_folds = n_folds,
  frequency_threshold = "> 0.5 (strict majority)",
  minimum_selected_folds = minimum_folds,
  fixed_gene_count_used = FALSE,
  validation_expression_read = FALSE,
  validation_outcomes_used = FALSE,
  gene_replacement_after_validation = FALSE
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "stability_consensus_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Frozen stability consensus: ", nrow(panel), " genes (",
  sum(panel$arm == "up"), " up, ", sum(panel$arm == "down"), " down)"
)
