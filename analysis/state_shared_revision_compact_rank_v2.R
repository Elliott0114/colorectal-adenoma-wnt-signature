#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(lme4)
  library(lmerTest)
})

options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260830)

root <- normalizePath(".", mustWork = TRUE)
parent_root <- file.path(root, "results", "state_aware_program_v1")
revision_root <- file.path(root, "results", "state_shared_revision_v2")
out_dir <- file.path(revision_root, "compact_rank")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  discovery = file.path(
    parent_root, "discovery_pseudobulk", "discovery_state_pseudobulk.rds"
  ),
  validation = file.path(
    parent_root, "validation_pseudobulk", "validation_state_pseudobulk.rds"
  ),
  common = file.path(
    parent_root, "common_effects", "cross_state_common_effects.tsv.gz"
  ),
  compact_panel = file.path(
    parent_root, "panel_derivation", "compact_state_shared_panel_frozen.tsv"
  ),
  portable_universe = file.path(
    parent_root,
    "panel_derivation",
    "portable_state_shared_candidate_universe.tsv"
  ),
  donor_disjoint_scores = file.path(
    revision_root,
    "donor_site",
    "donor_disjoint_programme_scores.tsv.gz"
  ),
  contract = file.path(
    root,
    "analysis",
    "contracts",
    "state_shared_revision_validation_v2_2026-08-30.md"
  )
)
if (!all(file.exists(unlist(paths)))) {
  stop("Run donor-disjoint validation before the compact rank analysis")
}

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
input_hashes <- vapply(paths, sha256, character(1))

discovery <- readRDS(paths$discovery)
validation <- readRDS(paths$validation)
common <- read.delim(paths$common, check.names = FALSE)
panel <- read.delim(paths$compact_panel, check.names = FALSE)
portable <- read.delim(paths$portable_universe, check.names = FALSE)
reference_scores <- read.delim(paths$donor_disjoint_scores, check.names = FALSE)

overlap <- intersect(
  unique(as.character(discovery$metadata$donor_id)),
  unique(as.character(validation$metadata$donor_id))
)
metadata_all <- as.data.frame(validation$metadata)
keep <- !metadata_all$donor_id %in% overlap
metadata <- droplevels(metadata_all[keep, , drop = FALSE])
counts <- validation$counts[, keep, drop = FALSE]

strict <- common[common$strict_state_shared %in% TRUE, , drop = FALSE]
up <- sort(strict$gene[strict$shared_direction == "up"])
down <- sort(strict$gene[strict$shared_direction == "down"])
genes <- c(up, down)
if (length(genes) != 1843L || !all(genes %in% rownames(counts))) {
  stop("The frozen programme is incomplete in held-out pseudobulk counts")
}
if (!all(panel$gene %in% genes)) {
  stop("The compact readout is not contained in the frozen programme")
}

log_cpm <- matrix(
  NA_real_,
  nrow = nrow(metadata),
  ncol = length(genes),
  dimnames = list(metadata$pseudobulk_id, genes)
)
for (state in c("ABS", "GOB", "TAC")) {
  index <- which(metadata$cell_type == state)
  dge <- DGEList(counts = as.matrix(counts[genes, index, drop = FALSE]))
  dge <- calcNormFactors(dge, method = "TMM")
  log_cpm[index, ] <- t(cpm(dge, log = TRUE, prior.count = 0.5))
}

within_profile_ranks <- t(apply(log_cpm, 1L, rank, ties.method = "average"))
within_profile_ranks <- within_profile_ranks / ncol(within_profile_ranks)
colnames(within_profile_ranks) <- genes
single_sample_full_rank_score <-
  rowMeans(within_profile_ranks[, up, drop = FALSE]) -
  rowMeans(within_profile_ranks[, down, drop = FALSE])
single_sample_compact_rank_score <-
  rowMeans(within_profile_ranks[, panel$gene[panel$arm == "up"], drop = FALSE]) -
  rowMeans(within_profile_ranks[, panel$gene[panel$arm == "down"], drop = FALSE])

scores <- metadata
scores$single_sample_full_rank_score <- single_sample_full_rank_score
scores$single_sample_compact_rank_score <- single_sample_compact_rank_score
reference_scores <- reference_scores[, c(
  "pseudobulk_id", "full_programme_score", "compact_8_score"
)]
scores <- merge(
  scores,
  reference_scores,
  by = "pseudobulk_id",
  all.x = TRUE,
  sort = FALSE
)
scores <- scores[match(metadata$pseudobulk_id, scores$pseudobulk_id), ]
if (anyNA(scores$full_programme_score)) {
  stop("Rank and reference score tables are misaligned")
}

safe_cor <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
}
cluster_bootstrap_correlation <- function(x, y, donor, iterations = 2000L) {
  donors <- sort(unique(as.character(donor)))
  values <- vapply(seq_len(iterations), function(iteration) {
    sampled <- sample(donors, length(donors), replace = TRUE)
    index <- unlist(lapply(sampled, function(value) which(donor == value)))
    safe_cor(x[index], y[index])
  }, numeric(1))
  c(
    estimate = safe_cor(x, y),
    ci_low = unname(quantile(values, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(values, 0.975, na.rm = TRUE))
  )
}

fidelity_records <- list()
record_index <- 0L
for (scope in c("all", "ABS", "GOB", "TAC")) {
  index <- if (scope == "all") {
    seq_len(nrow(scores))
  } else {
    which(scores$cell_type == scope)
  }
  for (target in c("full_programme_score", "single_sample_full_rank_score")) {
    estimate <- cluster_bootstrap_correlation(
      scores$single_sample_compact_rank_score[index],
      scores[[target]][index],
      scores$donor_id[index]
    )
    record_index <- record_index + 1L
    fidelity_records[[record_index]] <- data.frame(
      scope = scope,
      target = target,
      n_profiles = length(index),
      n_donors = length(unique(scores$donor_id[index])),
      spearman = estimate["estimate"],
      donor_bootstrap_ci_low = estimate["ci_low"],
      donor_bootstrap_ci_high = estimate["ci_high"],
      stringsAsFactors = FALSE
    )
  }
}
fidelity <- do.call(rbind, fidelity_records)

portable_eligible <- tolower(as.character(
  portable$objective_selection_eligible
)) %in% c("true", "t", "1", "yes")
portable_up <- intersect(
  portable$gene[portable$arm == "up" & portable_eligible],
  colnames(within_profile_ranks)
)
portable_down <- intersect(
  portable$gene[portable$arm == "down" & portable_eligible],
  colnames(within_profile_ranks)
)
if (length(portable_up) < 4L || length(portable_down) < 4L) {
  stop("The portable random-panel universe is too small")
}

n_random <- 2000L
random_panel_records <- vector("list", n_random)
random_benchmark <- data.frame(
  random_panel_id = seq_len(n_random),
  spearman_with_full_programme = NA_real_,
  stringsAsFactors = FALSE
)
for (iteration in seq_len(n_random)) {
  selected_up <- sample(portable_up, 4L, replace = FALSE)
  selected_down <- sample(portable_down, 4L, replace = FALSE)
  score <-
    rowMeans(within_profile_ranks[, selected_up, drop = FALSE]) -
    rowMeans(within_profile_ranks[, selected_down, drop = FALSE])
  random_benchmark$spearman_with_full_programme[iteration] <- safe_cor(
    score,
    scores$full_programme_score
  )
  random_panel_records[[iteration]] <- data.frame(
    random_panel_id = iteration,
    arm = c(rep("up", 4L), rep("down", 4L)),
    gene = c(selected_up, selected_down),
    stringsAsFactors = FALSE
  )
}
random_panels <- do.call(rbind, random_panel_records)
observed_fidelity <- fidelity$spearman[
  fidelity$scope == "all" & fidelity$target == "full_programme_score"
]
benchmark_summary <- data.frame(
  observed_compact_spearman = observed_fidelity,
  random_median = median(random_benchmark$spearman_with_full_programme),
  random_q95 = unname(quantile(
    random_benchmark$spearman_with_full_programme,
    0.95
  )),
  random_maximum = max(random_benchmark$spearman_with_full_programme),
  empirical_upper_tail_p = (1 + sum(
    random_benchmark$spearman_with_full_programme >= observed_fidelity
  )) / (n_random + 1),
  stringsAsFactors = FALSE
)

fit_effect <- function(data, scope) {
  local <- if (scope == "all") data else data[data$cell_type == scope, ]
  local$route <- relevel(factor(local$route), ref = "normal")
  local$donor_id <- factor(local$donor_id)
  local$specimen_id <- factor(local$specimen_id)
  local$cell_type <- factor(local$cell_type)
  formula <- if (scope == "all") {
    single_sample_compact_rank_score ~ route + cell_type +
      (1 | donor_id) + (1 | specimen_id)
  } else {
    single_sample_compact_rank_score ~ route + (1 | donor_id)
  }
  fit <- lmerTest::lmer(formula, data = local, REML = FALSE)
  coefficient <- summary(fit)$coefficients["routeconventional_adenoma", ]
  critical <- qt(0.975, coefficient["df"])
  data.frame(
    scope = scope,
    n_profiles = nrow(local),
    n_donors = length(unique(local$donor_id)),
    estimate = coefficient["Estimate"],
    standard_error = coefficient["Std. Error"],
    ci_low = coefficient["Estimate"] - critical * coefficient["Std. Error"],
    ci_high = coefficient["Estimate"] + critical * coefficient["Std. Error"],
    degrees_of_freedom = coefficient["df"],
    p_value = coefficient["Pr(>|t|)"],
    singular_fit = lme4::isSingular(fit, tol = 1e-5),
    stringsAsFactors = FALSE
  )
}
route_effects <- do.call(rbind, lapply(c("all", "ABS", "GOB", "TAC"), function(scope) {
  fit_effect(scores, scope)
}))
route_effects$p_value_BH <- p.adjust(route_effects$p_value, method = "BH")

gates <- data.frame(
  gate = c(
    "compact_rank_fidelity_at_least_0.75",
    "compact_rank_exceeds_random_q95",
    "compact_rank_route_effect_positive_overall",
    "compact_rank_route_effect_positive_all_states"
  ),
  passed = c(
    observed_fidelity >= 0.75,
    observed_fidelity > benchmark_summary$random_q95,
    route_effects$estimate[route_effects$scope == "all"] > 0,
    all(route_effects$estimate[route_effects$scope != "all"] > 0)
  ),
  stringsAsFactors = FALSE
)

write.table(
  scores,
  gzfile(file.path(out_dir, "heldout_single_sample_rank_scores.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  fidelity,
  file.path(out_dir, "heldout_single_sample_rank_fidelity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  route_effects,
  file.path(out_dir, "heldout_single_sample_rank_route_effects.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  random_benchmark,
  file.path(out_dir, "random_eight_gene_benchmark.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  random_panels,
  gzfile(file.path(out_dir, "random_eight_gene_panel_membership.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  benchmark_summary,
  file.path(out_dir, "random_eight_gene_benchmark_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  gates,
  file.path(out_dir, "quality_gates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- list(
  analysis = "state_shared_revision_compact_rank_v2",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = 20260830,
  input_sha256 = as.list(input_hashes),
  score_definition = paste(
    "within-profile percentile ranks across the 1,843-gene programme;",
    "mean rank of four up genes minus mean rank of four down genes"
  ),
  label_dependent_standardisation = FALSE,
  random_panels = n_random,
  random_universe_up = length(portable_up),
  random_universe_down = length(portable_down),
  quality_gates = stats::setNames(as.list(gates$passed), gates$gate),
  package_versions = as.list(vapply(
    c("edgeR", "lme4", "lmerTest"),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "analysis_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Compact single-sample rank analysis completed: observed rho=",
  sprintf("%.3f", observed_fidelity),
  "; random 95th percentile=",
  sprintf("%.3f", benchmark_summary$random_q95)
)
