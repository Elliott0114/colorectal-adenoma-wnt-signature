#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(metafor)
})

options(stringsAsFactors = FALSE)
set.seed(20260830)

root <- normalizePath(".", mustWork = TRUE)
parent_dir <- file.path(root, "results", "state_aware_program_v1", "external_validation")
out_dir <- file.path(root, "results", "state_shared_revision_v2", "external_meta")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  cohort_tests = file.path(parent_dir, "external_cohort_tests.tsv"),
  sample_scores = file.path(parent_dir, "external_sample_scores.tsv.gz"),
  contract = file.path(
    root,
    "analysis",
    "contracts",
    "state_shared_revision_validation_v2_2026-08-30.md"
  )
)
if (!all(file.exists(unlist(paths)))) {
  stop("External validation inputs are missing")
}

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
input_hashes <- vapply(paths, sha256, character(1))
cohort_tests <- read.delim(paths$cohort_tests, check.names = FALSE)
sample_scores <- read.delim(paths$sample_scores, check.names = FALSE)

meta_records <- list()
leaveout_records <- list()
record_index <- 0L
leaveout_index <- 0L
for (signature in c("state_shared_1843", "compact_8")) {
  local <- cohort_tests[
    cohort_tests$signature_id == signature &
      cohort_tests$comparison == "adenoma_vs_normal",
    ,
    drop = FALSE
  ]
  local$yi <- local$clustered_standardized_mean_difference
  local$sei <- (
    local$clustered_standardized_ci_high - local$clustered_standardized_ci_low
  ) / (2 * qnorm(0.975))
  if (nrow(local) != 5L || any(!is.finite(local$yi)) || any(local$sei <= 0)) {
    stop("Cohort effect inputs are incomplete for ", signature)
  }
  fit <- rma.uni(
    yi = local$yi,
    sei = local$sei,
    method = "REML",
    test = "knha",
    slab = local$cohort
  )
  prediction <- predict(fit, level = 95)
  record_index <- record_index + 1L
  meta_records[[record_index]] <- data.frame(
    signature_id = signature,
    excluded_cohort = "__NONE__",
    n_cohorts = fit$k,
    pooled_standardized_effect = as.numeric(fit$b),
    standard_error = fit$se,
    ci_low = fit$ci.lb,
    ci_high = fit$ci.ub,
    p_value = fit$pval,
    tau_squared = fit$tau2,
    i_squared = fit$I2,
    q_statistic = fit$QE,
    q_p_value = fit$QEp,
    prediction_interval_low = prediction$pi.lb,
    prediction_interval_high = prediction$pi.ub,
    method = "REML random effects with Hartung-Knapp interval",
    stringsAsFactors = FALSE
  )
  for (excluded in local$cohort) {
    reduced <- local[local$cohort != excluded, , drop = FALSE]
    reduced_fit <- rma.uni(
      yi = reduced$yi,
      sei = reduced$sei,
      method = "REML",
      test = "knha",
      slab = reduced$cohort
    )
    reduced_prediction <- predict(reduced_fit, level = 95)
    leaveout_index <- leaveout_index + 1L
    leaveout_records[[leaveout_index]] <- data.frame(
      signature_id = signature,
      excluded_cohort = excluded,
      n_cohorts = reduced_fit$k,
      pooled_standardized_effect = as.numeric(reduced_fit$b),
      standard_error = reduced_fit$se,
      ci_low = reduced_fit$ci.lb,
      ci_high = reduced_fit$ci.ub,
      p_value = reduced_fit$pval,
      tau_squared = reduced_fit$tau2,
      i_squared = reduced_fit$I2,
      prediction_interval_low = reduced_prediction$pi.lb,
      prediction_interval_high = reduced_prediction$pi.ub,
      stringsAsFactors = FALSE
    )
  }
}
meta_summary <- do.call(rbind, meta_records)
leaveout_summary <- do.call(rbind, leaveout_records)

patient_group_values <- function(data, group_field) {
  aggregate(
    programme_score ~ patient_cluster_id + group,
    transform(data, group = data[[group_field]]),
    mean
  )
}

contrast_groups <- function(data, group_a, group_b, contrast_name) {
  local <- data[data$group %in% c(group_a, group_b), , drop = FALSE]
  a <- local$programme_score[local$group == group_a]
  b <- local$programme_score[local$group == group_b]
  if (length(a) < 3L || length(b) < 3L) {
    return(NULL)
  }
  test <- t.test(a, b, paired = FALSE, var.equal = FALSE)
  mann_whitney <- wilcox.test(a, b, exact = FALSE)
  pooled_sd <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) /
    (length(a) + length(b) - 2))
  hedges_correction <- 1 - 3 / (4 * (length(a) + length(b)) - 9)
  hedges_g <- hedges_correction * (mean(a) - mean(b)) / pooled_sd
  data.frame(
    contrast = contrast_name,
    group_a = group_a,
    group_b = group_b,
    n_a = length(a),
    n_b = length(b),
    mean_a = mean(a),
    mean_b = mean(b),
    mean_difference = mean(a) - mean(b),
    ci_low = test$conf.int[1L],
    ci_high = test$conf.int[2L],
    welch_t_statistic = unname(test$statistic),
    degrees_of_freedom = unname(test$parameter),
    welch_p_value = test$p.value,
    hedges_g = hedges_g,
    mann_whitney_w = unname(mann_whitney$statistic),
    mann_whitney_p_value = mann_whitney$p.value,
    stringsAsFactors = FALSE
  )
}

specificity_records <- list()
specificity_index <- 0L
for (signature in c("state_shared_1843", "compact_8")) {
  gse40362 <- sample_scores[
    sample_scores$signature_id == signature & sample_scores$cohort == "GSE40362",
    ,
    drop = FALSE
  ]
  gse40362 <- patient_group_values(gse40362, "tissue_group")
  for (definition in list(
    c("adenoma", "normal", "adenoma_vs_normal"),
    c("hyperplastic", "normal", "hyperplastic_vs_normal"),
    c("adenoma", "hyperplastic", "adenoma_vs_hyperplastic")
  )) {
    result <- contrast_groups(gse40362, definition[1], definition[2], definition[3])
    if (!is.null(result)) {
      specificity_index <- specificity_index + 1L
      result$cohort <- "GSE40362"
      result$signature_id <- signature
      result$grouping_variable <- "tissue_group"
      specificity_records[[specificity_index]] <- result
    }
  }

  gse41657 <- sample_scores[
    sample_scores$signature_id == signature & sample_scores$cohort == "GSE41657",
    ,
    drop = FALSE
  ]
  gse41657 <- patient_group_values(gse41657, "grade_group")
  for (definition in list(
    c("low_grade", "normal", "low_grade_adenoma_vs_normal"),
    c("high_grade", "normal", "high_grade_adenoma_vs_normal"),
    c("high_grade", "low_grade", "high_vs_low_grade_adenoma"),
    c("crc", "normal", "crc_vs_normal"),
    c("crc", "high_grade", "crc_vs_high_grade_adenoma")
  )) {
    result <- contrast_groups(gse41657, definition[1], definition[2], definition[3])
    if (!is.null(result)) {
      specificity_index <- specificity_index + 1L
      result$cohort <- "GSE41657"
      result$signature_id <- signature
      result$grouping_variable <- "grade_group"
      specificity_records[[specificity_index]] <- result
    }
  }
}
specificity <- do.call(rbind, specificity_records)
specificity$welch_p_value_BH <- ave(
  specificity$welch_p_value,
  specificity$signature_id,
  FUN = function(value) p.adjust(value, method = "BH")
)

full_meta <- meta_summary[meta_summary$signature_id == "state_shared_1843", ]
full_leaveout <- leaveout_summary[
  leaveout_summary$signature_id == "state_shared_1843",
  ,
  drop = FALSE
]
adenoma_hyperplastic <- specificity[
  specificity$signature_id == "state_shared_1843" &
    specificity$contrast == "adenoma_vs_hyperplastic",
  ,
  drop = FALSE
]
adenoma_hyperplastic_supported <- nrow(adenoma_hyperplastic) == 1L &&
  adenoma_hyperplastic$ci_low > 0 &&
  adenoma_hyperplastic$welch_p_value_BH <= 0.05
broad_specificity_supported <- FALSE
gates <- data.frame(
  gate = c(
    "random_effects_full_programme_positive",
    "all_leave_one_cohort_out_full_programme_positive",
    "adenoma_vs_hyperplastic_supported_in_GSE40362",
    "broad_adenoma_specificity_supported"
  ),
  passed = c(
    full_meta$pooled_standardized_effect > 0,
    all(full_leaveout$pooled_standardized_effect > 0),
    adenoma_hyperplastic_supported,
    broad_specificity_supported
  ),
  consequence = c(
    "retain independent-cohort transfer",
    "retain leave-one-cohort-out robustness",
    ifelse(
      adenoma_hyperplastic_supported,
      "report this as a one-cohort histology boundary result",
      "do not claim separation from hyperplastic polyps"
    ),
    "do not describe the programme as broadly adenoma-specific without serrated and inflammatory controls"
  ),
  stringsAsFactors = FALSE
)

write.table(
  meta_summary,
  file.path(out_dir, "random_effects_meta_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  leaveout_summary,
  file.path(out_dir, "random_effects_leave_one_cohort_out.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  specificity,
  file.path(out_dir, "histology_and_grade_boundary_analysis.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  gates,
  file.path(out_dir, "quality_gates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- list(
  analysis = "state_shared_revision_external_meta_v2",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = 20260830,
  input_sha256 = as.list(input_hashes),
  meta_analysis = list(
    effect = "patient-clustered standardized mean difference per cohort",
    tau_estimator = "REML",
    interval = "Hartung-Knapp",
    prediction_interval = TRUE
  ),
  boundary_analyses = c(
    "GSE40362 hyperplastic polyps",
    "GSE41657 low-grade adenoma, high-grade adenoma and CRC"
  ),
  specificity_is_primary = FALSE,
  quality_gates = stats::setNames(as.list(gates$passed), gates$gate),
  package_versions = list(metafor = as.character(packageVersion("metafor"))),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "analysis_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Random-effects synthesis completed: full-programme pooled effect=",
  sprintf("%.3f", full_meta$pooled_standardized_effect),
  "; one-cohort adenoma-vs-hyperplastic gate=",
  adenoma_hyperplastic_supported,
  "; broad specificity gate=FALSE"
)
