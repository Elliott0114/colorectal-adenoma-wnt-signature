#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BiocParallel)
  library(edgeR)
  library(limma)
  library(lme4)
  library(lmerTest)
  library(mashr)
  library(Matrix)
  library(rhdf5)
  library(tidyr)
  library(variancePartition)
})

options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260830)

root <- normalizePath(".", mustWork = TRUE)
parent_root <- file.path(root, "results", "state_aware_program_v1")
out_dir <- file.path(root, "results", "state_shared_revision_v2", "donor_site")
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
  prior_validation_effects = file.path(
    parent_root, "validation_models", "state_specific_primary_effects.tsv.gz"
  ),
  compact_panel = file.path(
    parent_root, "panel_derivation", "compact_state_shared_panel_frozen.tsv"
  ),
  validation_h5ad = file.path(
    root,
    "data_sources",
    "Chen_Cell_2021_CELLxGENE",
    "chen_validation_epithelial.h5ad"
  ),
  contract = file.path(
    root,
    "analysis",
    "contracts",
    "state_shared_revision_validation_v2_2026-08-30.md"
  )
)
if (!all(file.exists(unlist(paths)))) {
  stop("At least one input required for donor-disjoint validation is missing")
}

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
input_hashes <- vapply(paths, sha256, character(1))

states <- c("ABS", "GOB", "TAC")
discovery <- readRDS(paths$discovery)
validation <- readRDS(paths$validation)
discovery_metadata <- as.data.frame(discovery$metadata)
validation_metadata_all <- as.data.frame(validation$metadata)
common <- read.delim(paths$common, check.names = FALSE)
prior_validation_effects <- read.delim(
  paths$prior_validation_effects,
  check.names = FALSE
)
panel <- read.delim(paths$compact_panel, check.names = FALSE)

overlap <- intersect(
  unique(as.character(discovery_metadata$donor_id)),
  unique(as.character(validation_metadata_all$donor_id))
)
keep_columns <- !validation_metadata_all$donor_id %in% overlap
validation_metadata <- droplevels(validation_metadata_all[keep_columns, , drop = FALSE])
validation_counts <- validation$counts[, keep_columns, drop = FALSE]
if (!identical(colnames(validation_counts), validation_metadata$pseudobulk_id)) {
  stop("Filtered validation counts and metadata are misaligned")
}
if (length(overlap) == 0L) {
  warning("No overlapping donors were detected; the exclusion changed no profiles")
}

strict <- common[common$strict_state_shared %in% TRUE, , drop = FALSE]
strict_up <- sort(strict$gene[strict$shared_direction == "up"])
strict_down <- sort(strict$gene[strict$shared_direction == "down"])
strict_genes <- c(strict_up, strict_down)
if (length(strict_genes) != 1843L || nrow(panel) != 8L) {
  stop("The frozen programme or compact readout differs from the parent analysis")
}

read_h5ad_categorical <- function(path, field) {
  base <- paste0("/obs/", field)
  codes <- h5read(path, paste0(base, "/codes"))
  categories <- h5read(path, paste0(base, "/categories"))
  output <- rep(NA_character_, length(codes))
  valid <- codes >= 0L
  output[valid] <- as.character(categories[codes[valid] + 1L])
  output
}

specimen_metadata <- data.frame(
  specimen_id = read_h5ad_categorical(paths$validation_h5ad, "HTAN Specimen ID"),
  donor_id_h5ad = read_h5ad_categorical(paths$validation_h5ad, "donor_id"),
  tissue = read_h5ad_categorical(paths$validation_h5ad, "tissue"),
  age_label = read_h5ad_categorical(paths$validation_h5ad, "development_stage"),
  sex = read_h5ad_categorical(paths$validation_h5ad, "sex"),
  polyp_type = read_h5ad_categorical(paths$validation_h5ad, "Polyp_Type"),
  stringsAsFactors = FALSE
)
specimen_metadata <- unique(specimen_metadata)
if (anyDuplicated(specimen_metadata$specimen_id)) {
  consistency <- aggregate(
    cbind(tissue, age_label, sex, donor_id_h5ad) ~ specimen_id,
    specimen_metadata,
    function(x) length(unique(x))
  )
  if (any(consistency[, -1L] != 1L)) {
    stop("A specimen maps to inconsistent anatomical or donor metadata")
  }
  specimen_metadata <- specimen_metadata[
    !duplicated(specimen_metadata$specimen_id),
    ,
    drop = FALSE
  ]
}
specimen_metadata$age_years <- as.numeric(sub("-year-old stage", "", specimen_metadata$age_label))
proximal_sites <- c(
  "ascending colon",
  "transverse colon",
  "hepatic cecum",
  "hepatic flexure of colon"
)
distal_sites <- c("descending colon", "sigmoid colon", "rectum")
specimen_metadata$site_group <- ifelse(
  specimen_metadata$tissue %in% proximal_sites,
  "proximal",
  ifelse(specimen_metadata$tissue %in% distal_sites, "distal", NA_character_)
)
if (anyNA(specimen_metadata$site_group) || anyNA(specimen_metadata$age_years)) {
  stop("Deposited anatomical site or age could not be mapped")
}

calculate_log_cpm <- function(counts, metadata, genes) {
  output <- matrix(
    NA_real_,
    nrow = nrow(metadata),
    ncol = length(genes),
    dimnames = list(metadata$pseudobulk_id, genes)
  )
  for (state in states) {
    index <- which(metadata$cell_type == state)
    dge <- DGEList(counts = as.matrix(counts[genes, index, drop = FALSE]))
    dge <- calcNormFactors(dge, method = "TMM")
    output[index, ] <- t(cpm(dge, log = TRUE, prior.count = 0.5))
  }
  if (any(!is.finite(output))) {
    stop("Non-finite TMM log-CPM values were generated")
  }
  output
}

discovery_expression <- calculate_log_cpm(
  discovery$counts,
  discovery_metadata,
  strict_genes
)
validation_expression <- calculate_log_cpm(
  validation_counts,
  validation_metadata,
  strict_genes
)

centres <- scales <- matrix(
  NA_real_,
  nrow = length(states),
  ncol = length(strict_genes),
  dimnames = list(states, strict_genes)
)
for (state in states) {
  index <- discovery_metadata$cell_type == state
  centres[state, ] <- colMeans(discovery_expression[index, , drop = FALSE])
  scales[state, ] <- apply(discovery_expression[index, , drop = FALSE], 2L, sd)
}
if (any(!is.finite(scales)) || any(scales <= 0)) {
  stop("Invalid discovery scaling parameters")
}

validation_z <- matrix(
  NA_real_,
  nrow = nrow(validation_metadata),
  ncol = length(strict_genes),
  dimnames = list(validation_metadata$pseudobulk_id, strict_genes)
)
for (state in states) {
  index <- validation_metadata$cell_type == state
  validation_z[index, ] <- sweep(
    sweep(validation_expression[index, , drop = FALSE], 2L, centres[state, ], "-"),
    2L,
    scales[state, ],
    "/"
  )
}

equal_arm_score <- function(matrix, up, down) {
  rowMeans(matrix[, up, drop = FALSE]) - rowMeans(matrix[, down, drop = FALSE])
}
score_table <- merge(
  validation_metadata,
  specimen_metadata,
  by = "specimen_id",
  all.x = TRUE,
  sort = FALSE
)
score_table <- score_table[match(validation_metadata$pseudobulk_id, score_table$pseudobulk_id), ]
if (anyNA(score_table$tissue) || any(score_table$donor_id != score_table$donor_id_h5ad)) {
  stop("Pseudobulk and deposited specimen metadata do not agree")
}
score_table$full_programme_score <- equal_arm_score(
  validation_z,
  strict_up,
  strict_down
)
score_table$compact_8_score <- equal_arm_score(
  validation_z,
  panel$gene[panel$arm == "up"],
  panel$gene[panel$arm == "down"]
)

fit_score_model <- function(data, outcome, scope, adjusted) {
  local <- data
  if (scope != "all") {
    local <- local[local$cell_type == scope, , drop = FALSE]
  }
  local$route <- relevel(factor(local$route), ref = "normal")
  local$cell_type <- factor(local$cell_type, levels = states)
  local$site_group <- relevel(factor(local$site_group), ref = "distal")
  local$sex <- factor(local$sex)
  local$donor_id <- factor(local$donor_id)
  local$specimen_id <- factor(local$specimen_id)
  fixed <- if (scope == "all") "route + cell_type" else "route"
  if (adjusted) {
    fixed <- paste(fixed, "+ site_group + age_years + sex")
  }
  random_effects <- if (scope == "all") {
    "+ (1 | donor_id) + (1 | specimen_id)"
  } else {
    "+ (1 | donor_id)"
  }
  formula <- as.formula(paste(outcome, "~", fixed, random_effects))
  fit <- lmerTest::lmer(formula, data = local, REML = FALSE)
  coefficient <- summary(fit)$coefficients["routeconventional_adenoma", ]
  critical <- qt(0.975, coefficient["df"])
  data.frame(
    score = outcome,
    scope = scope,
    model = ifelse(adjusted, "site_age_sex_adjusted", "unadjusted"),
    n_profiles = nrow(local),
    n_specimens = length(unique(local$specimen_id)),
    n_donors = length(unique(local$donor_id)),
    estimate = coefficient["Estimate"],
    standard_error = coefficient["Std. Error"],
    degrees_of_freedom = coefficient["df"],
    ci_low = coefficient["Estimate"] - critical * coefficient["Std. Error"],
    ci_high = coefficient["Estimate"] + critical * coefficient["Std. Error"],
    p_value = coefficient["Pr(>|t|)"],
    singular_fit = lme4::isSingular(fit, tol = 1e-5),
    stringsAsFactors = FALSE
  )
}

score_effects <- do.call(rbind, lapply(
  c("full_programme_score", "compact_8_score"),
  function(outcome) do.call(rbind, lapply(c("all", states), function(scope) {
    rbind(
      fit_score_model(score_table, outcome, scope, FALSE),
      fit_score_model(score_table, outcome, scope, TRUE)
    )
  }))
))
score_effects$p_value_BH <- p.adjust(score_effects$p_value, method = "BH")

same_site_records <- list()
record_index <- 0L
for (outcome in c("full_programme_score", "compact_8_score")) {
  for (state in states) {
    local <- score_table[score_table$cell_type == state, , drop = FALSE]
    donor_site_route <- aggregate(
      local[[outcome]],
      by = list(
        donor_id = local$donor_id,
        tissue = local$tissue,
        route = local$route
      ),
      FUN = mean
    )
    names(donor_site_route)[4L] <- "score"
    paired <- pivot_wider(
      donor_site_route,
      names_from = route,
      values_from = score
    )
    paired <- paired[
      is.finite(paired$normal) & is.finite(paired$conventional_adenoma),
      ,
      drop = FALSE
    ]
    paired$paired_difference <- paired$conventional_adenoma - paired$normal
    donor_difference <- aggregate(paired_difference ~ donor_id, paired, mean)
    difference <- donor_difference$paired_difference
    t_result <- if (length(difference) >= 3L) t.test(difference) else NULL
    wilcoxon_p <- if (length(difference) >= 3L && !all(difference == 0)) {
      suppressWarnings(wilcox.test(difference, exact = TRUE)$p.value)
    } else if (length(difference) >= 3L) {
      1
    } else {
      NA_real_
    }
    sign_p <- if (length(difference)) {
      binom.test(sum(difference > 0), sum(difference != 0), p = 0.5)$p.value
    } else {
      NA_real_
    }
    record_index <- record_index + 1L
    same_site_records[[record_index]] <- data.frame(
      score = outcome,
      cell_type = state,
      n_donors = length(difference),
      n_donor_site_pairs = nrow(paired),
      mean_difference = ifelse(length(difference), mean(difference), NA_real_),
      median_difference = ifelse(length(difference), median(difference), NA_real_),
      ci_low = if (is.null(t_result)) NA_real_ else t_result$conf.int[1L],
      ci_high = if (is.null(t_result)) NA_real_ else t_result$conf.int[2L],
      paired_t_p_value = if (is.null(t_result)) NA_real_ else t_result$p.value,
      wilcoxon_p_value = wilcoxon_p,
      positive_donors = sum(difference > 0),
      nonzero_donors = sum(difference != 0),
      sign_test_p_value = sign_p,
      stringsAsFactors = FALSE
    )
  }
}
same_site_effects <- do.call(rbind, same_site_records)

message("Refitting donor-disjoint state-specific gene effects")
validation_gene_universe <- Reduce(
  intersect,
  lapply(states, function(state) {
    prior_validation_effects$gene[prior_validation_effects$cell_type == state]
  })
)
validation_gene_universe <- rownames(validation_counts)[
  rownames(validation_counts) %in% validation_gene_universe
]
parallel_param <- SnowParam(
  workers = 4L,
  type = "SOCK",
  progressbar = TRUE,
  stop.on.error = FALSE
)

state_effects_list <- list()
for (state in states) {
  column_index <- which(validation_metadata$cell_type == state)
  state_metadata <- droplevels(validation_metadata[column_index, , drop = FALSE])
  state_metadata$route <- relevel(factor(state_metadata$route), ref = "normal")
  state_metadata$donor_id <- factor(state_metadata$donor_id)
  state_counts <- as.matrix(
    validation_counts[validation_gene_universe, column_index, drop = FALSE]
  )
  rownames(state_metadata) <- colnames(state_counts)
  dge <- DGEList(counts = round(state_counts))
  dge <- calcNormFactors(dge, method = "TMM")
  formula <- ~ route + (1 | donor_id)
  voom_object <- voomWithDreamWeights(
    dge,
    formula,
    state_metadata,
    plot = FALSE,
    BPPARAM = parallel_param
  )
  raw_fit <- dream(
    voom_object,
    formula,
    state_metadata,
    ddf = "Satterthwaite",
    BPPARAM = parallel_param
  )
  coefficient_name <- "routeconventional_adenoma"
  coefficient_index <- match(coefficient_name, colnames(raw_fit$coefficients))
  if (is.na(coefficient_index)) {
    stop("The donor-disjoint route coefficient is absent for ", state)
  }
  raw_standard_error <- raw_fit$stdev.unscaled[, coefficient_index] * raw_fit$sigma
  moderated <- variancePartition::eBayes(raw_fit, robust = TRUE)
  table <- variancePartition::topTable(
    moderated,
    coef = coefficient_name,
    number = Inf,
    sort.by = "none"
  )
  table$gene <- rownames(table)
  table$cell_type <- state
  table$raw_se <- as.numeric(raw_standard_error[table$gene])
  state_effects_list[[state]] <- table[, c(
    "gene", "cell_type", "logFC", "raw_se", "P.Value", "adj.P.Val"
  )]
}
state_effects <- do.call(rbind, state_effects_list)

genes <- Reduce(
  intersect,
  lapply(states, function(state) state_effects$gene[state_effects$cell_type == state])
)
extract_matrix <- function(column) {
  matrix <- vapply(states, function(state) {
    local <- state_effects[state_effects$cell_type == state, , drop = FALSE]
    local[[column]][match(genes, local$gene)]
  }, numeric(length(genes)))
  rownames(matrix) <- genes
  colnames(matrix) <- states
  matrix
}
bhat <- extract_matrix("logFC")
shat <- extract_matrix("raw_se")
null_correlation <- estimate_null_correlation_simple(
  mash_set_data(bhat, shat),
  z_thresh = 2,
  est_cor = TRUE
)
ones <- rep(1, length(states))
gls <- t(vapply(seq_along(genes), function(index) {
  covariance <- diag(shat[index, ]) %*% null_correlation %*% diag(shat[index, ])
  precision <- solve(covariance)
  denominator <- as.numeric(t(ones) %*% precision %*% ones)
  effect <- as.numeric(t(ones) %*% precision %*% bhat[index, ] / denominator)
  standard_error <- sqrt(1 / denominator)
  c(common_effect = effect, common_se = standard_error, common_z = effect / standard_error)
}, numeric(3L)))
donor_disjoint_common <- data.frame(
  gene = genes,
  gls,
  common_p_value = 2 * pnorm(-abs(gls[, "common_z"])),
  stringsAsFactors = FALSE
)
for (state in states) {
  donor_disjoint_common[[paste0("logFC_", state)]] <- bhat[, state]
}
discovery_fields <- common[, c(
  "gene", "common_effect", "common_z", "strict_state_shared", "shared_direction"
)]
names(discovery_fields)[2:3] <- c("discovery_common_effect", "discovery_common_z")
donor_disjoint_common <- merge(
  donor_disjoint_common,
  discovery_fields,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
donor_disjoint_common$expected_sign <- ifelse(
  donor_disjoint_common$shared_direction == "up", 1,
  ifelse(donor_disjoint_common$shared_direction == "down", -1, NA_real_)
)
donor_disjoint_common$common_direction_match <-
  sign(donor_disjoint_common$common_effect) == donor_disjoint_common$expected_sign
donor_disjoint_common$all_state_direction_match <- apply(
  cbind(
    sign(donor_disjoint_common$logFC_ABS) == donor_disjoint_common$expected_sign,
    sign(donor_disjoint_common$logFC_GOB) == donor_disjoint_common$expected_sign,
    sign(donor_disjoint_common$logFC_TAC) == donor_disjoint_common$expected_sign
  ),
  1L,
  all
)
strict_testable <- donor_disjoint_common[
  donor_disjoint_common$strict_state_shared %in% TRUE,
  ,
  drop = FALSE
]
replication_summary <- data.frame(
  n_strict_programme_genes = length(strict_genes),
  n_testable_strict_genes = nrow(strict_testable),
  common_direction_match_n = sum(strict_testable$common_direction_match),
  common_direction_match_fraction = mean(strict_testable$common_direction_match),
  all_state_direction_match_n = sum(strict_testable$all_state_direction_match),
  all_state_direction_match_fraction = mean(strict_testable$all_state_direction_match),
  discovery_validation_effect_spearman = cor(
    strict_testable$discovery_common_effect,
    strict_testable$common_effect,
    method = "spearman"
  ),
  discovery_validation_z_spearman = cor(
    strict_testable$discovery_common_z,
    strict_testable$common_z,
    method = "spearman"
  ),
  stringsAsFactors = FALSE
)

same_site_full <- same_site_effects[
  same_site_effects$score == "full_programme_score",
  ,
  drop = FALSE
]
gates <- data.frame(
  gate = c(
    "overall_full_score_positive",
    "all_state_full_scores_positive",
    "strict_common_direction_at_least_0.90",
    "effect_correlation_at_least_0.80",
    "site_adjusted_full_score_positive",
    "same_site_state_effects_positive"
  ),
  passed = c(
    score_effects$estimate[
      score_effects$score == "full_programme_score" &
        score_effects$scope == "all" &
        score_effects$model == "unadjusted"
    ] > 0,
    all(score_effects$estimate[
      score_effects$score == "full_programme_score" &
        score_effects$scope %in% states &
        score_effects$model == "unadjusted"
    ] > 0),
    replication_summary$common_direction_match_fraction >= 0.90,
    replication_summary$discovery_validation_effect_spearman >= 0.80,
    score_effects$estimate[
      score_effects$score == "full_programme_score" &
        score_effects$scope == "all" &
        score_effects$model == "site_age_sex_adjusted"
    ] > 0,
    nrow(same_site_full) == length(states) &&
      all(same_site_full$n_donors >= 5L) &&
      all(same_site_full$mean_difference > 0)
  ),
  stringsAsFactors = FALSE
)

exclusion_audit <- data.frame(
  discovery_donors = length(unique(discovery_metadata$donor_id)),
  validation_donors_before = length(unique(validation_metadata_all$donor_id)),
  overlapping_donors_removed = length(overlap),
  validation_donors_after = length(unique(validation_metadata$donor_id)),
  validation_profiles_before = nrow(validation_metadata_all),
  validation_profiles_after = nrow(validation_metadata),
  validation_specimens_after = length(unique(validation_metadata$specimen_id)),
  stringsAsFactors = FALSE
)

write.table(
  score_table,
  gzfile(file.path(out_dir, "donor_disjoint_programme_scores.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  score_effects,
  file.path(out_dir, "donor_disjoint_score_effects.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  same_site_effects,
  file.path(out_dir, "same_site_paired_score_effects.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  state_effects,
  gzfile(file.path(out_dir, "donor_disjoint_state_gene_effects.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  donor_disjoint_common,
  gzfile(file.path(out_dir, "donor_disjoint_common_gene_effects.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  replication_summary,
  file.path(out_dir, "donor_disjoint_replication_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  exclusion_audit,
  file.path(out_dir, "donor_overlap_exclusion_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  gates,
  file.path(out_dir, "quality_gates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- list(
  analysis = "state_shared_revision_donor_site_v2",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = 20260830,
  input_sha256 = as.list(input_hashes),
  discovery_validation_overlap_count = length(overlap),
  validation_donors_after_exclusion = length(unique(validation_metadata$donor_id)),
  biological_replicate = "donor",
  repeated_units = c("specimen", "broad epithelial state"),
  anatomical_adjustment = c("proximal_or_distal_site", "age_years", "sex"),
  programme_selection_reopened = FALSE,
  quality_gates = stats::setNames(as.list(gates$passed), gates$gate),
  package_versions = as.list(vapply(
    c(
      "BiocParallel", "edgeR", "limma", "lme4", "lmerTest", "mashr",
      "Matrix", "rhdf5", "tidyr", "variancePartition"
    ),
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
  "Donor-disjoint validation completed: ",
  exclusion_audit$validation_donors_after,
  " donors; strict-direction replication=",
  sprintf("%.3f", replication_summary$common_direction_match_fraction),
  "; effect rho=",
  sprintf("%.3f", replication_summary$discovery_validation_effect_spearman)
)
