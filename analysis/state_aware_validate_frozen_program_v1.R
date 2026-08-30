#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(lmerTest)
  library(mashr)
})

options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260829)

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
discovery_pseudobulk_path <- file.path(
  result_root,
  "discovery_pseudobulk",
  "discovery_state_pseudobulk.rds"
)
validation_pseudobulk_path <- file.path(
  result_root,
  "validation_pseudobulk",
  "validation_state_pseudobulk.rds"
)
discovery_common_path <- file.path(
  result_root,
  "common_effects",
  "cross_state_common_effects.tsv.gz"
)
validation_effect_path <- file.path(
  result_root,
  "validation_models",
  "state_specific_primary_effects.tsv.gz"
)
validation_paired_path <- file.path(
  result_root,
  "validation_models",
  "state_specific_paired_sensitivity.tsv.gz"
)
panel_path <- file.path(
  result_root,
  "panel_derivation",
  "compact_state_shared_panel_frozen.tsv"
)
out_dir <- file.path(result_root, "heldout_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  discovery_pseudobulk_path,
  validation_pseudobulk_path,
  discovery_common_path,
  validation_effect_path,
  validation_paired_path,
  panel_path
)
if (!all(file.exists(inputs))) {
  stop("One or more frozen programme-validation inputs are missing")
}

expected_hashes <- c(
  discovery_pseudobulk =
    "2c63422dd5cf8124098974fec62648bb24b546adb407ec41d794c8709cdb4c96",
  validation_pseudobulk =
    "3adbe6ad44ef9baf7fd06d1bf6ec2b077483a3399308ed7ce6456186889b6d19",
  discovery_common =
    "a1ac4b7b67ac279782e04e386971d7463e169cbdbc24d7da4c314e34a4f3e946",
  frozen_panel =
    "c5997e572342a72da8441df312fba4e3461cacfa5e30d0d8590a3f23ae3d96f0"
)
actual_hashes <- c(
  discovery_pseudobulk = digest::digest(
    discovery_pseudobulk_path,
    algo = "sha256",
    file = TRUE
  ),
  validation_pseudobulk = digest::digest(
    validation_pseudobulk_path,
    algo = "sha256",
    file = TRUE
  ),
  discovery_common = digest::digest(
    discovery_common_path,
    algo = "sha256",
    file = TRUE
  ),
  frozen_panel = digest::digest(panel_path, algo = "sha256", file = TRUE)
)
if (!identical(actual_hashes, expected_hashes)) {
  stop("A frozen discovery input or analysis contract changed before validation")
}

states <- c("ABS", "GOB", "TAC")
discovery <- readRDS(discovery_pseudobulk_path)
validation <- readRDS(validation_pseudobulk_path)
discovery_metadata <- as.data.frame(discovery$metadata)
validation_metadata <- as.data.frame(validation$metadata)
discovery_common <- read.delim(discovery_common_path, check.names = FALSE)
validation_effects <- read.delim(validation_effect_path, check.names = FALSE)
validation_paired <- read.delim(validation_paired_path, check.names = FALSE)
panel <- read.delim(panel_path, check.names = FALSE)

if (!setequal(as.character(unique(validation_metadata$cell_type)), states)) {
  stop("Held-out pseudobulk object does not contain exactly ABS, GOB and TAC")
}
if (!setequal(unique(validation_effects$cell_type), states)) {
  stop("Held-out model output has unexpected cell-state labels")
}
if (anyDuplicated(panel$gene) || nrow(panel) != 8L) {
  stop("The frozen compact readout is not the expected eight-gene panel")
}

strict <- discovery_common[discovery_common$strict_state_shared, , drop = FALSE]
strict_up <- sort(strict$gene[strict$shared_direction == "up"])
strict_down <- sort(strict$gene[strict$shared_direction == "down"])
strict_genes <- c(strict_up, strict_down)
if (length(strict_genes) != 1843L) {
  stop("The frozen strict state-shared programme changed size")
}
if (!all(strict_genes %in% rownames(discovery$counts)) ||
    !all(strict_genes %in% rownames(validation$counts))) {
  stop("At least one frozen strict-core gene is absent from a pseudobulk matrix")
}
if (!all(panel$gene %in% strict_genes)) {
  stop("The compact panel is not a subset of the frozen strict programme")
}

calculate_log_cpm <- function(pseudobulk, metadata, genes) {
  output <- matrix(
    NA_real_,
    nrow = nrow(metadata),
    ncol = length(genes),
    dimnames = list(metadata$pseudobulk_id, genes)
  )
  for (state in states) {
    index <- which(metadata$cell_type == state)
    dge <- DGEList(
      counts = as.matrix(pseudobulk$counts[genes, index, drop = FALSE])
    )
    dge <- calcNormFactors(dge, method = "TMM")
    output[index, ] <- t(cpm(dge, log = TRUE, prior.count = 0.5))
  }
  if (any(!is.finite(output))) {
    stop("Non-finite TMM log-CPM value in score input")
  }
  output
}

message("Calculating frozen discovery-standardised programme scores")
discovery_expression <- calculate_log_cpm(
  discovery,
  discovery_metadata,
  strict_genes
)
validation_expression <- calculate_log_cpm(
  validation,
  validation_metadata,
  strict_genes
)

discovery_centres <- matrix(
  NA_real_,
  nrow = length(states),
  ncol = length(strict_genes),
  dimnames = list(states, strict_genes)
)
discovery_scales <- discovery_centres
for (state in states) {
  index <- which(discovery_metadata$cell_type == state)
  discovery_centres[state, ] <- colMeans(
    discovery_expression[index, , drop = FALSE]
  )
  discovery_scales[state, ] <- apply(
    discovery_expression[index, , drop = FALSE],
    2,
    sd
  )
}
if (any(!is.finite(discovery_scales)) || any(discovery_scales <= 0)) {
  stop("A strict-core gene has invalid discovery scaling parameters")
}
scaling_parameters <- do.call(rbind, lapply(states, function(state) {
  data.frame(
    cell_type = state,
    gene = strict_genes,
    discovery_centre = discovery_centres[state, ],
    discovery_scale = discovery_scales[state, ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}))

validation_z <- matrix(
  NA_real_,
  nrow = nrow(validation_metadata),
  ncol = length(strict_genes),
  dimnames = list(validation_metadata$pseudobulk_id, strict_genes)
)
for (state in states) {
  index <- which(validation_metadata$cell_type == state)
  validation_z[index, ] <- sweep(
    sweep(
      validation_expression[index, , drop = FALSE],
      2,
      discovery_centres[state, ],
      "-"
    ),
    2,
    discovery_scales[state, ],
    "/"
  )
}

score_equal_arms <- function(z, up_genes, down_genes) {
  rowMeans(z[, up_genes, drop = FALSE]) -
    rowMeans(z[, down_genes, drop = FALSE])
}

panel_up <- panel$gene[panel$arm == "up"]
panel_down <- panel$gene[panel$arm == "down"]

score_table <- validation_metadata
score_table$full_programme_score <- score_equal_arms(
  validation_z,
  strict_up,
  strict_down
)
score_table$compact_8_score <- score_equal_arms(
  validation_z,
  panel_up,
  panel_down
)

safe_cor <- function(x, y, method = "spearman") {
  value <- suppressWarnings(cor(x, y, method = method, use = "complete.obs"))
  ifelse(is.finite(value), value, NA_real_)
}

clustered_correlation <- function(x, y, donor, n_bootstrap = 2000L) {
  donors <- sort(unique(as.character(donor)))
  bootstrap <- vapply(
    seq_len(n_bootstrap),
    function(iteration) {
      sampled <- sample(donors, length(donors), replace = TRUE)
      index <- unlist(
        lapply(sampled, function(value) which(donor == value)),
        use.names = FALSE
      )
      safe_cor(x[index], y[index])
    },
    numeric(1)
  )
  c(
    spearman = safe_cor(x, y),
    ci_low = unname(quantile(bootstrap, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(bootstrap, 0.975, na.rm = TRUE)),
    bootstrap_se = sd(bootstrap, na.rm = TRUE)
  )
}

message("Evaluating held-out compact-readout fidelity")
fidelity_records <- list()
record_index <- 0L
for (scope in c("all", states)) {
  index <- if (scope == "all") {
    seq_len(nrow(score_table))
  } else {
    which(score_table$cell_type == scope)
  }
  for (readout in "compact_8_score") {
    record_index <- record_index + 1L
    estimate <- clustered_correlation(
      score_table[[readout]][index],
      score_table$full_programme_score[index],
      score_table$donor_id[index]
    )
    fidelity_records[[record_index]] <- data.frame(
      scope = scope,
      readout = readout,
      n_profiles = length(index),
      n_donors = length(unique(score_table$donor_id[index])),
      n_genes = 8L,
      spearman = estimate["spearman"],
      donor_bootstrap_ci_low = estimate["ci_low"],
      donor_bootstrap_ci_high = estimate["ci_high"],
      donor_bootstrap_se = estimate["bootstrap_se"],
      stringsAsFactors = FALSE
    )
  }
}
fidelity <- do.call(rbind, fidelity_records)

fit_score_effect <- function(data, score, scope) {
  local <- data
  local$route <- relevel(factor(local$route), ref = "normal")
  local$donor_id <- factor(local$donor_id)
  local$specimen_id <- factor(local$specimen_id)
  local$cell_type <- factor(local$cell_type, levels = states)
  formula <- if (scope == "all") {
    as.formula(paste(
      score,
      "~ route + cell_type + (1 | donor_id) + (1 | specimen_id)"
    ))
  } else {
    as.formula(paste(score, "~ route + (1 | donor_id)"))
  }
  fit <- lmerTest::lmer(formula, data = local, REML = FALSE)
  coefficient <- summary(fit)$coefficients["routeconventional_adenoma", ]
  critical_value <- qt(0.975, df = coefficient["df"])
  data.frame(
    scope = scope,
    score = score,
    n_profiles = nrow(local),
    n_donors = nlevels(local$donor_id),
    n_adenoma_profiles = sum(local$route == "conventional_adenoma"),
    n_normal_profiles = sum(local$route == "normal"),
    estimate = coefficient["Estimate"],
    standard_error = coefficient["Std. Error"],
    degrees_of_freedom = coefficient["df"],
    ci_low = coefficient["Estimate"] - critical_value * coefficient["Std. Error"],
    ci_high = coefficient["Estimate"] + critical_value * coefficient["Std. Error"],
    t_statistic = coefficient["t value"],
    p_value = coefficient["Pr(>|t|)"],
    singular_fit = lme4::isSingular(fit, tol = 1e-5),
    stringsAsFactors = FALSE
  )
}

score_effects <- list()
record_index <- 0L
for (scope in c("all", states)) {
  local <- if (scope == "all") {
    score_table
  } else {
    score_table[score_table$cell_type == scope, , drop = FALSE]
  }
  for (score in c(
    "full_programme_score",
    "compact_8_score"
  )) {
    record_index <- record_index + 1L
    score_effects[[record_index]] <- fit_score_effect(local, score, scope)
  }
}
score_effects <- do.call(rbind, score_effects)
score_effects$p_value_BH <- p.adjust(score_effects$p_value, method = "BH")

paired_score_records <- list()
record_index <- 0L
for (state in states) {
  for (score in c(
    "full_programme_score",
    "compact_8_score"
  )) {
    local <- score_table[score_table$cell_type == state, , drop = FALSE]
    donor_route <- aggregate(
      local[[score]],
      by = list(donor_id = local$donor_id, route = local$route),
      FUN = mean
    )
    colnames(donor_route)[3L] <- "score"
    paired <- tidyr::pivot_wider(
      donor_route,
      names_from = route,
      values_from = score
    )
    paired <- paired[
      is.finite(paired$normal) & is.finite(paired$conventional_adenoma),
      ,
      drop = FALSE
    ]
    difference <- paired$conventional_adenoma - paired$normal
    t_test <- t.test(difference, mu = 0)
    wilcoxon <- suppressWarnings(wilcox.test(difference, mu = 0, exact = FALSE))
    record_index <- record_index + 1L
    paired_score_records[[record_index]] <- data.frame(
      cell_type = state,
      score = score,
      n_paired_donors = length(difference),
      mean_paired_difference = mean(difference),
      median_paired_difference = median(difference),
      paired_t_p_value = t_test$p.value,
      wilcoxon_p_value = wilcoxon$p.value,
      stringsAsFactors = FALSE
    )
  }
}
paired_score_effects <- do.call(rbind, paired_score_records)
paired_score_effects$paired_t_p_value_BH <- p.adjust(
  paired_score_effects$paired_t_p_value,
  method = "BH"
)
paired_score_effects$wilcoxon_p_value_BH <- p.adjust(
  paired_score_effects$wilcoxon_p_value,
  method = "BH"
)

message("Integrating held-out gene-level effects across states")
validation_genes <- Reduce(
  intersect,
  lapply(states, function(state) {
    validation_effects$gene[validation_effects$cell_type == state]
  })
)
validation_genes <- sort(validation_genes)
extract_effect_matrix <- function(column) {
  output <- vapply(
    states,
    function(state) {
      local <- validation_effects[
        validation_effects$cell_type == state,
        ,
        drop = FALSE
      ]
      local[[column]][match(validation_genes, local$gene)]
    },
    numeric(length(validation_genes))
  )
  rownames(output) <- validation_genes
  colnames(output) <- states
  output
}

bhat <- extract_effect_matrix("logFC")
shat <- extract_effect_matrix("raw_se")
if (any(!is.finite(bhat)) || any(!is.finite(shat)) || any(shat <= 0)) {
  stop("Invalid held-out gene-level effect estimates")
}
null_correlation <- estimate_null_correlation_simple(
  mash_set_data(bhat, shat),
  z_thresh = 2,
  est_cor = TRUE
)
dimnames(null_correlation) <- list(states, states)
if (min(eigen(null_correlation, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) {
  stop("Held-out null-correlation matrix is not positive definite")
}

ones <- rep(1, length(states))
gls <- t(vapply(
  seq_along(validation_genes),
  function(index) {
    covariance <- diag(shat[index, ]) %*%
      null_correlation %*%
      diag(shat[index, ])
    precision <- solve(covariance)
    denominator <- as.numeric(t(ones) %*% precision %*% ones)
    effect <- as.numeric(t(ones) %*% precision %*% bhat[index, ] / denominator)
    standard_error <- sqrt(1 / denominator)
    residual <- bhat[index, ] - effect
    heterogeneity_q <- as.numeric(t(residual) %*% precision %*% residual)
    c(
      common_effect = effect,
      common_se = standard_error,
      common_z = effect / standard_error,
      heterogeneity_q = heterogeneity_q,
      heterogeneity_p = pchisq(
        heterogeneity_q,
        df = length(states) - 1L,
        lower.tail = FALSE
      )
    )
  },
  numeric(5L)
))

validation_common <- data.frame(
  gene = validation_genes,
  gls,
  common_p_value = 2 * pnorm(-abs(gls[, "common_z"])),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
validation_common$common_q_value <- p.adjust(
  validation_common$common_p_value,
  method = "BH"
)
for (state in states) {
  validation_common[[paste0("logFC_", state)]] <- bhat[, state]
  validation_common[[paste0("raw_se_", state)]] <- shat[, state]
}

discovery_columns <- c(
  "gene",
  "common_effect",
  "common_z",
  "common_q_value",
  "strict_state_shared",
  "shared_direction"
)
colnames(discovery_common)[match(
  c("common_effect", "common_z", "common_q_value"),
  colnames(discovery_common)
)] <- c(
  "discovery_common_effect",
  "discovery_common_z",
  "discovery_common_q_value"
)
discovery_columns <- c(
  "gene",
  "discovery_common_effect",
  "discovery_common_z",
  "discovery_common_q_value",
  "strict_state_shared",
  "shared_direction"
)
validation_common <- merge(
  validation_common,
  discovery_common[discovery_columns],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
validation_common$expected_sign <- ifelse(
  validation_common$shared_direction == "up",
  1,
  ifelse(
    validation_common$shared_direction == "down",
    -1,
    sign(validation_common$discovery_common_effect)
  )
)
validation_common$common_direction_match <-
  sign(validation_common$common_effect) == validation_common$expected_sign
validation_common$all_state_directions_match <- vapply(
  seq_len(nrow(validation_common)),
  function(index) {
    expected <- validation_common$expected_sign[index]
    if (!is.finite(expected)) {
      return(NA)
    }
    all(sign(bhat[validation_common$gene[index], ]) == expected)
  },
  logical(1)
)

paired_sign_match <- matrix(
  NA,
  nrow = nrow(validation_common),
  ncol = length(states),
  dimnames = list(validation_common$gene, states)
)
for (state in states) {
  local <- validation_paired[validation_paired$cell_type == state, , drop = FALSE]
  paired_sign_match[, state] <- sign(
    local$logFC[match(validation_common$gene, local$gene)]
  ) == validation_common$expected_sign
}
validation_common$all_paired_state_directions_match <- apply(
  paired_sign_match,
  1,
  function(value) if (anyNA(value)) NA else all(value)
)

testable_discovery <- validation_common[
  is.finite(validation_common$discovery_common_effect),
  ,
  drop = FALSE
]
testable_strict <- validation_common[
  validation_common$strict_state_shared %in% TRUE,
  ,
  drop = FALSE
]

summarise_replication <- function(data, label) {
  expected <- data$expected_sign
  data.frame(
    set = label,
    n_genes = nrow(data),
    common_direction_match_n = sum(data$common_direction_match, na.rm = TRUE),
    common_direction_match_fraction = mean(data$common_direction_match, na.rm = TRUE),
    all_state_direction_match_n = sum(data$all_state_directions_match, na.rm = TRUE),
    all_state_direction_match_fraction = mean(
      data$all_state_directions_match,
      na.rm = TRUE
    ),
    all_paired_state_direction_match_fraction = mean(
      data$all_paired_state_directions_match,
      na.rm = TRUE
    ),
    nominal_common_p_le_0.05_correct_n = sum(
      data$common_direction_match & data$common_p_value <= 0.05,
      na.rm = TRUE
    ),
    heldout_common_q_le_0.05_correct_n = sum(
      data$common_direction_match & data$common_q_value <= 0.05,
      na.rm = TRUE
    ),
    ABS_direction_match_fraction = mean(
      sign(data$logFC_ABS) == expected,
      na.rm = TRUE
    ),
    GOB_direction_match_fraction = mean(
      sign(data$logFC_GOB) == expected,
      na.rm = TRUE
    ),
    TAC_direction_match_fraction = mean(
      sign(data$logFC_TAC) == expected,
      na.rm = TRUE
    ),
    discovery_validation_effect_spearman = safe_cor(
      data$discovery_common_effect,
      data$common_effect
    ),
    discovery_validation_z_spearman = safe_cor(
      data$discovery_common_z,
      data$common_z
    ),
    stringsAsFactors = FALSE
  )
}

replication_summary <- rbind(
  summarise_replication(testable_discovery, "all_discovery_testable"),
  summarise_replication(testable_strict, "strict_state_shared"),
  summarise_replication(
    testable_strict[testable_strict$shared_direction == "up", , drop = FALSE],
    "strict_state_shared_up"
  ),
  summarise_replication(
    testable_strict[testable_strict$shared_direction == "down", , drop = FALSE],
    "strict_state_shared_down"
  )
)

rank_statistic <- validation_common$common_z
names(rank_statistic) <- validation_common$gene
gene_sets <- list(
  frozen_strict_up = intersect(strict_up, names(rank_statistic)),
  frozen_strict_down = intersect(strict_down, names(rank_statistic))
)
enrichment <- cameraPR(
  statistic = rank_statistic,
  index = gene_sets,
  use.ranks = TRUE,
  sort = FALSE
)
enrichment$gene_set <- rownames(enrichment)
expected_enrichment_direction <- c(
  frozen_strict_up = "Up",
  frozen_strict_down = "Down"
)
enrichment$expected_direction <- unname(
  expected_enrichment_direction[enrichment$gene_set]
)
enrichment$direction_matches_expectation <-
  enrichment$Direction == enrichment$expected_direction
if (!"FDR" %in% colnames(enrichment)) {
  enrichment$FDR <- p.adjust(enrichment$PValue, method = "BH")
}
enrichment <- enrichment[, c(
  "gene_set",
  "NGenes",
  "expected_direction",
  "Direction",
  "direction_matches_expectation",
  "PValue",
  "FDR"
)]

validation_panel_fields <- validation_common[, c(
  "gene",
  "common_effect",
  "common_se",
  "common_z",
  "common_p_value",
  "common_q_value",
  "heterogeneity_p",
  "logFC_ABS",
  "logFC_GOB",
  "logFC_TAC",
  "all_paired_state_directions_match"
)]
colnames(validation_panel_fields)[2:7] <- paste0(
  "validation_",
  colnames(validation_panel_fields)[2:7]
)
panel_validation <- merge(
  panel,
  validation_panel_fields,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
panel_validation$expected_sign_panel <- ifelse(
  panel_validation$arm == "up",
  1,
  -1
)
panel_validation$validation_common_direction_match <-
  sign(panel_validation$validation_common_effect) ==
  panel_validation$expected_sign_panel
panel_validation$validation_all_state_directions_match <-
  sign(panel_validation$logFC_ABS) == panel_validation$expected_sign_panel &
  sign(panel_validation$logFC_GOB) == panel_validation$expected_sign_panel &
  sign(panel_validation$logFC_TAC) == panel_validation$expected_sign_panel
panel_validation <- panel_validation[
  order(panel_validation$pair_step, -panel_validation$route_weight),
  ,
  drop = FALSE
]

write.table(
  score_table,
  gzfile(file.path(out_dir, "heldout_programme_scores.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  fidelity,
  file.path(out_dir, "heldout_readout_fidelity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  score_effects,
  file.path(out_dir, "heldout_score_route_effects.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  paired_score_effects,
  file.path(out_dir, "heldout_paired_score_effects.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  validation_common,
  gzfile(file.path(out_dir, "heldout_cross_state_common_effects.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  replication_summary,
  file.path(out_dir, "heldout_gene_level_replication_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  enrichment,
  file.path(out_dir, "heldout_frozen_set_enrichment.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  panel_validation,
  file.path(out_dir, "heldout_compact_panel_gene_validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  scaling_parameters,
  gzfile(file.path(out_dir, "frozen_discovery_scaling_parameters.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(state = states, null_correlation, check.names = FALSE),
  file.path(out_dir, "heldout_estimated_null_correlation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(
    cell_type = rep(states, each = 2L),
    parameter = rep(c("centre", "scale"), times = length(states)),
    gene_count = length(strict_genes),
    minimum = as.numeric(rbind(
      apply(discovery_centres, 1, min),
      apply(discovery_scales, 1, min)
    )),
    maximum = as.numeric(rbind(
      apply(discovery_centres, 1, max),
      apply(discovery_scales, 1, max)
    )),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "discovery_scaling_parameter_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_validate_frozen_program_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  validation_outcomes_used_for_selection = FALSE,
  frozen_inputs = as.list(actual_hashes),
  generated_validation_inputs = list(
    state_effects_path = normalizePath(validation_effect_path),
    state_effects_sha256 = digest::digest(
      validation_effect_path,
      algo = "sha256",
      file = TRUE
    ),
    paired_effects_path = normalizePath(validation_paired_path),
    paired_effects_sha256 = digest::digest(
      validation_paired_path,
      algo = "sha256",
      file = TRUE
    )
  ),
  score = list(
    full_programme_genes = length(strict_genes),
    full_up_genes = length(strict_up),
    full_down_genes = length(strict_down),
    compact_genes = nrow(panel),
    normalisation = "TMM log2 CPM with prior.count 0.5 within state",
    scaling = "discovery state-specific centres and standard deviations",
    scoring = "mean z(up) minus mean z(down); equal arm weights",
    fidelity_interval = "2,000 whole-donor bootstrap replicates"
  ),
  gene_level = list(
    validation_genes = length(validation_genes),
    integration = "correlation-adjusted GLS across ABS, GOB and TAC",
    enrichment = "limma::cameraPR use.ranks=TRUE"
  ),
  random_seed = 20260829,
  package_versions = as.list(vapply(
    c("edgeR", "limma", "lme4", "lmerTest", "mashr", "Matrix"),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "heldout_validation_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Held-out validation completed: ",
  length(validation_genes),
  " genes; compact rho=",
  sprintf(
    "%.3f",
    fidelity$spearman[
      fidelity$scope == "all" & fidelity$readout == "compact_8_score"
    ]
  )
)
