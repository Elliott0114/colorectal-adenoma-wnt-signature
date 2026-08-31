#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
  library(metafor)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
source_root <- file.path(result_root, "functional_architecture_v1")
run_root_env <- Sys.getenv("STATE_AWARE_MODULE_RUN_ROOT", unset = "")
out_root <- if (nzchar(run_root_env)) {
  normalizePath(run_root_env, mustWork = FALSE)
} else {
  source_root
}
wgcna_root <- file.path(source_root, "consensus_wgcna")
external_root <- file.path(out_root, "module_external_validation")
selection_path_env <- Sys.getenv("STATE_AWARE_MODULE_SELECTION_PATH", unset = "")
summary_path <- if (nzchar(selection_path_env)) {
  normalizePath(selection_path_env, mustWork = TRUE)
} else {
  file.path(wgcna_root, "module_internal_gate_summary.tsv")
}
route_column <- Sys.getenv(
  "STATE_AWARE_MODULE_ROUTE_COLUMN", unset = "internal_gate_pass"
)
effect_path <- file.path(external_root, "module_external_cohort_effects.tsv")
addendum_path <- file.path(
  root, "analysis", "contracts",
  "state_aware_functional_architecture_downstream_addendum_v1_2026-08-30.md"
)
expected_addendum_hash <- "32343afe117d09007066fbe01f8fbe7cf4a11ee628dbe6a121f0f814968da3bb"

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)
if (sha256(addendum_path) != expected_addendum_hash) {
  stop("The frozen downstream validation addendum changed")
}
if (!file.exists(summary_path)) {
  stop("The frozen module summary is missing")
}

module_summary <- read.delim(summary_path, check.names = FALSE)
if (!route_column %in% colnames(module_summary)) {
  stop("Routing column is missing: ", route_column)
}
as_bool <- function(value) {
  if (is.logical(value)) value else tolower(as.character(value)) == "true"
}
if (!file.exists(effect_path)) {
  if (!any(as_bool(module_summary[[route_column]]))) {
    module_summary$external_gate_pass <- FALSE
    module_summary$routing_status <- "audit_analysis_route_failed"
    write.table(
      module_summary, file.path(out_root, "module_validation.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    quit(save = "no", status = 0)
  }
  stop("External cohort effects are missing for internally eligible modules")
}

effects <- read.delim(effect_path, check.names = FALSE)

fit_meta <- function(frame) {
  frame <- frame[
    is.finite(frame$oriented_adenoma_effect_sd) &
      is.finite(frame$standard_error) & frame$standard_error > 0,
    , drop = FALSE
  ]
  if (nrow(frame) < 2L) {
    return(data.frame(
      k = nrow(frame), pooled_effect = NA_real_, pooled_se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
      tau2_reml = NA_real_, i2_percent = NA_real_, stringsAsFactors = FALSE
    ))
  }
  fit <- rma.uni(
    yi = frame$oriented_adenoma_effect_sd,
    sei = frame$standard_error,
    method = "REML", test = "knha"
  )
  data.frame(
    k = nrow(frame), pooled_effect = unname(fit$b[1L]),
    pooled_se = fit$se, ci_low = fit$ci.lb, ci_high = fit$ci.ub,
    p_value = fit$pval, tau2_reml = fit$tau2, i2_percent = fit$I2,
    stringsAsFactors = FALSE
  )
}

meta_records <- list()
loo_records <- list()
meta_counter <- 0L
loo_counter <- 0L
for (module in unique(effects$module)) {
  local <- effects[effects$module == module, , drop = FALSE]
  pooled <- fit_meta(local)
  loo_values <- numeric(0)
  if (nrow(local) >= 3L) {
    for (cohort in sort(unique(local$cohort))) {
      loo <- fit_meta(local[local$cohort != cohort, , drop = FALSE])
      loo_counter <- loo_counter + 1L
      loo_records[[loo_counter]] <- data.frame(
        module = module, excluded_cohort = cohort, loo,
        stringsAsFactors = FALSE
      )
      loo_values <- c(loo_values, loo$pooled_effect)
    }
  }
  loo_stable <- length(loo_values) == nrow(local) &&
    all(is.finite(loo_values)) && all(loo_values > 0)
  external_pass <- nrow(local) >= 3L && is.finite(pooled$ci_low) &&
    pooled$ci_low > 0 && loo_stable
  meta_counter <- meta_counter + 1L
  meta_records[[meta_counter]] <- data.frame(
    module = module,
    n_coverage_eligible_cohorts = nrow(local),
    pooled,
    all_leave_one_cohort_out_positive = loo_stable,
    external_gate_pass = external_pass,
    stringsAsFactors = FALSE
  )
}

meta_table <- do.call(rbind, meta_records)
loo_table <- if (length(loo_records)) {
  do.call(rbind, loo_records)
} else {
  data.frame()
}
validation <- merge(module_summary, meta_table, by = "module", all.x = TRUE)
validation$external_gate_pass[is.na(validation$external_gate_pass)] <- FALSE
route_pass <- as_bool(validation[[route_column]])
validation$routing_status <- ifelse(
  !route_pass,
  "audit_analysis_route_failed",
  ifelse(
    validation$external_gate_pass,
    "pending_technical_orthogonal_and_perturbation_gates",
    "audit_external_gate_failed"
  )
)

write.table(
  meta_table, file.path(external_root, "module_external_random_effects.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  loo_table, file.path(external_root, "module_external_leave_one_cohort_out.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  validation, file.path(out_root, "module_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_module_meta_analysis_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  model = "metafor REML random effects with Knapp-Hartung inference",
  route_column = route_column,
  external_gate = paste(
    "at least three coverage-eligible cohorts; pooled 95% CI above zero;",
    "all leave-one-cohort-out estimates above zero"
  ),
  modules_tested = unique(effects$module),
  modules_passing_external_gate = meta_table$module[meta_table$external_gate_pass],
  input_sha256 = list(
    module_summary = sha256(summary_path),
    cohort_effects = sha256(effect_path),
    downstream_addendum = sha256(addendum_path)
  ),
  package_versions = list(
    R = as.character(getRversion()),
    metafor = as.character(packageVersion("metafor")),
    digest = as.character(packageVersion("digest")),
    jsonlite = as.character(packageVersion("jsonlite"))
  ),
  session_info = capture.output(sessionInfo())
)
write_json(
  manifest, file.path(external_root, "module_external_meta_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)

message(
  "REML-KH module meta-analysis complete: ", nrow(meta_table), " tested; ",
  sum(meta_table$external_gate_pass), " passed the external gate"
)
