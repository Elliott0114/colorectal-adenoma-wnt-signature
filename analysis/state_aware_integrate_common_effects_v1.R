#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mashr)
})

options(stringsAsFactors = FALSE)
set.seed(20260829)

root <- normalizePath(".", mustWork = TRUE)
model_dir <- file.path(root, "results", "state_aware_program_v1", "discovery_models")
input_path <- file.path(model_dir, "state_specific_primary_effects.tsv.gz")
paired_path <- file.path(model_dir, "state_specific_paired_sensitivity.tsv.gz")
contract_path <- file.path(
  root,
  "analysis",
  "contracts",
  "state_aware_program_rederivation_v1_2026-08-29.md"
)
out_dir <- file.path(root, "results", "state_aware_program_v1", "common_effects")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_contract_sha256 <-
  "0f39a03154408d14bb0bbe0c1e4f55d498e94d403bd6b45d4be5dd8918555f69"
actual_contract_sha256 <- digest::digest(contract_path, algo = "sha256", file = TRUE)
if (!identical(actual_contract_sha256, expected_contract_sha256)) {
  stop("The frozen parent analysis contract changed before cross-state integration")
}
if (!file.exists(input_path) || !file.exists(paired_path)) {
  stop("Discovery state-model outputs are incomplete")
}

effects <- read.delim(input_path, check.names = FALSE)
paired <- read.delim(paired_path, check.names = FALSE)
states <- c("ABS", "GOB", "TAC")
if (!setequal(unique(effects$cell_type), states)) {
  stop("Unexpected state labels in primary model output")
}
if (anyDuplicated(effects[c("gene", "cell_type")])) {
  stop("Duplicated gene-state estimates in primary model output")
}

genes <- Reduce(
  intersect,
  lapply(states, function(state) effects$gene[effects$cell_type == state])
)
genes <- sort(genes)
if (length(genes) < 1000L) {
  stop("Unexpectedly small common gene universe")
}

extract_matrix <- function(column) {
  output <- vapply(
    states,
    function(state) {
      state_table <- effects[effects$cell_type == state, , drop = FALSE]
      state_table[[column]][match(genes, state_table$gene)]
    },
    numeric(length(genes))
  )
  rownames(output) <- genes
  colnames(output) <- states
  output
}

bhat <- extract_matrix("logFC")
shat <- extract_matrix("raw_se")
if (any(!is.finite(bhat)) || any(!is.finite(shat)) || any(shat <= 0)) {
  stop("Invalid effect estimates or standard errors supplied to mashr")
}

message("Estimating residual correlation across ", nrow(bhat), " genes")
data_independent <- mash_set_data(bhat, shat)
fit_one_by_one <- mash_1by1(data_independent)
strong <- get_significant_results(fit_one_by_one, thresh = 0.05)
if (length(strong) < ncol(bhat) + 2L) {
  warning(
    "Fewer than ",
    ncol(bhat) + 2L,
    " strong effects were available; canonical covariance matrices only will be used"
  )
}

null_correlation <- estimate_null_correlation_simple(
  data_independent,
  z_thresh = 2,
  est_cor = TRUE
)
dimnames(null_correlation) <- list(states, states)
eigenvalues <- eigen(null_correlation, symmetric = TRUE, only.values = TRUE)$values
if (min(eigenvalues) <= 1e-8) {
  stop("Estimated null-correlation matrix is not positive definite")
}
data_correlated <- mash_set_data(bhat, shat, V = null_correlation)

covariances <- cov_canonical(data_correlated)
covariance_mode <- "canonical_only"
if (length(strong) >= ncol(bhat) + 2L) {
  pca_covariances <- cov_pca(
    data_correlated,
    npc = ncol(bhat),
    subset = strong
  )
  ed_covariances <- cov_ed(
    data_correlated,
    pca_covariances,
    subset = strong
  )
  covariances <- c(covariances, ed_covariances)
  covariance_mode <- "canonical_plus_data_driven"
}

message("Fitting multivariate adaptive shrinkage model")
mash_fit <- mash(
  data_correlated,
  Ulist = covariances,
  prior = "nullbiased",
  verbose = FALSE,
  seed = 20260829
)
posterior_mean <- get_pm(mash_fit)
lfsr <- get_lfsr(mash_fit)
rownames(posterior_mean) <- genes
colnames(posterior_mean) <- states
rownames(lfsr) <- genes
colnames(lfsr) <- states

ones <- rep(1, length(states))
gls_result <- t(vapply(
  seq_len(nrow(bhat)),
  function(index) {
    covariance <- diag(shat[index, ]) %*%
      null_correlation %*%
      diag(shat[index, ])
    precision <- solve(covariance)
    denominator <- as.numeric(t(ones) %*% precision %*% ones)
    common_effect <- as.numeric(
      t(ones) %*% precision %*% bhat[index, ] / denominator
    )
    common_se <- sqrt(1 / denominator)
    residual <- bhat[index, ] - common_effect
    heterogeneity_q <- as.numeric(t(residual) %*% precision %*% residual)
    c(
      common_effect = common_effect,
      common_se = common_se,
      common_z = common_effect / common_se,
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
rownames(gls_result) <- genes

output <- data.frame(
  gene = genes,
  gls_result,
  common_p_value = 2 * pnorm(-abs(gls_result[, "common_z"])),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
output$common_q_value <- p.adjust(output$common_p_value, method = "BH")

for (state in states) {
  output[[paste0("logFC_", state)]] <- bhat[, state]
  output[[paste0("raw_se_", state)]] <- shat[, state]
  output[[paste0("posterior_mean_", state)]] <- posterior_mean[, state]
  output[[paste0("lfsr_", state)]] <- lfsr[, state]
}

posterior_sign <- sign(posterior_mean)
same_nonzero_sign <- apply(
  posterior_sign,
  1,
  function(value) all(value > 0) || all(value < 0)
)
output$same_posterior_direction <- same_nonzero_sign
output$max_lfsr <- apply(lfsr, 1, max)
output$strict_state_shared <-
  output$common_q_value <= 0.05 &
  output$same_posterior_direction &
  output$max_lfsr <= 0.05
output$relaxed_state_shared <-
  output$common_q_value <= 0.05 &
  output$same_posterior_direction &
  output$max_lfsr <= 0.10
output$shared_direction <- ifelse(
  output$strict_state_shared | output$relaxed_state_shared,
  ifelse(rowMeans(posterior_mean) > 0, "up", "down"),
  "not_shared"
)

paired_direction <- matrix(
  NA_real_,
  nrow = length(genes),
  ncol = length(states),
  dimnames = list(genes, states)
)
for (state in states) {
  state_table <- paired[paired$cell_type == state, , drop = FALSE]
  paired_direction[, state] <- sign(
    state_table$logFC[match(genes, state_table$gene)]
  )
}
output$paired_all_states_match_common <- vapply(
  seq_along(genes),
  function(index) {
    value <- paired_direction[index, ]
    all(is.finite(value)) && all(value == sign(output$common_effect[index]))
  },
  logical(1)
)

output <- output[order(-abs(output$common_z), output$gene), , drop = FALSE]
output$absolute_rank <- seq_len(nrow(output))
output$signed_rank <- rank(-output$common_z, ties.method = "average")

write.table(
  output,
  gzfile(file.path(out_dir, "cross_state_common_effects.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(
    state = states,
    null_correlation,
    check.names = FALSE
  ),
  file.path(out_dir, "estimated_null_correlation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
saveRDS(
  list(
    fit = mash_fit,
    bhat = bhat,
    shat = shat,
    posterior_mean = posterior_mean,
    lfsr = lfsr,
    null_correlation = null_correlation,
    strong_indices = strong
  ),
  file.path(out_dir, "mashr_fit.rds"),
  compress = "xz"
)

summary_table <- data.frame(
  metric = c(
    "genes_tested",
    "one_by_one_lfsr_le_0.05",
    "common_q_le_0.05",
    "same_posterior_direction",
    "strict_state_shared",
    "strict_state_shared_up",
    "strict_state_shared_down",
    "relaxed_state_shared",
    "heterogeneity_q_le_0.05"
  ),
  value = c(
    nrow(output),
    length(strong),
    sum(output$common_q_value <= 0.05),
    sum(output$same_posterior_direction),
    sum(output$strict_state_shared),
    sum(output$strict_state_shared & output$shared_direction == "up"),
    sum(output$strict_state_shared & output$shared_direction == "down"),
    sum(output$relaxed_state_shared),
    sum(output$heterogeneity_p <= 0.05)
  ),
  stringsAsFactors = FALSE
)
write.table(
  summary_table,
  file.path(out_dir, "common_effect_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_integrate_common_effects_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_path = normalizePath(input_path),
  input_sha256 = digest::digest(input_path, algo = "sha256", file = TRUE),
  paired_input_sha256 = digest::digest(paired_path, algo = "sha256", file = TRUE),
  contract_sha256 = actual_contract_sha256,
  genes_tested = nrow(output),
  states = states,
  correlation_method = "mashr::estimate_null_correlation_simple(z_thresh=2)",
  common_effect_method = "correlation-adjusted inverse-covariance GLS",
  mash_covariance_mode = covariance_mode,
  strict_definition = "common BH q<=0.05; concordant posterior sign; LFSR<=0.05 in all states",
  relaxed_definition = "common BH q<=0.05; concordant posterior sign; LFSR<=0.10 in all states",
  random_seed = 20260829,
  package_versions = as.list(vapply(
    c("mashr", "ashr", "digest", "jsonlite"),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  summary = summary_table,
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "common_effect_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Cross-state integration complete: ",
  sum(output$strict_state_shared),
  " strict and ",
  sum(output$relaxed_state_shared),
  " relaxed shared genes"
)
