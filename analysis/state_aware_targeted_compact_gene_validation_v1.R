#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(lmerTest)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
pseudobulk_path <- file.path(
  result_root,
  "validation_pseudobulk",
  "validation_state_pseudobulk.rds"
)
filter_path <- file.path(
  result_root,
  "validation_models",
  "expression_filter_audit.tsv.gz"
)
panel_path <- file.path(
  result_root,
  "panel_derivation",
  "compact_state_shared_panel_frozen.tsv"
)
null_correlation_path <- file.path(
  result_root,
  "heldout_validation",
  "heldout_estimated_null_correlation.tsv"
)
out_dir <- file.path(result_root, "heldout_validation")

expected_hashes <- c(
  pseudobulk =
    "3adbe6ad44ef9baf7fd06d1bf6ec2b077483a3399308ed7ce6456186889b6d19",
  panel =
    "c5997e572342a72da8441df312fba4e3461cacfa5e30d0d8590a3f23ae3d96f0"
)
actual_hashes <- c(
  pseudobulk = digest::digest(pseudobulk_path, algo = "sha256", file = TRUE),
  panel = digest::digest(panel_path, algo = "sha256", file = TRUE)
)
if (!identical(expected_hashes, actual_hashes)) {
  stop("A frozen compact-gene validation input changed")
}

pseudobulk <- readRDS(pseudobulk_path)
metadata <- as.data.frame(pseudobulk$metadata)
filter_audit <- read.delim(filter_path, check.names = FALSE)
panel <- read.delim(panel_path, check.names = FALSE)
null_table <- read.delim(null_correlation_path, check.names = FALSE)
states <- c("ABS", "GOB", "TAC")
panel_genes <- panel$gene[order(panel$pair_step, -panel$route_weight)]
shared_genes <- filter_audit$gene[filter_audit$in_shared_universe]
shared_genes <- rownames(pseudobulk$counts)[
  rownames(pseudobulk$counts) %in% shared_genes
]
if (length(panel_genes) != 8L || length(shared_genes) != 7115L) {
  stop("Unexpected compact-panel or validation-universe size")
}

records <- list()
record_index <- 0L
paired_records <- list()
paired_index <- 0L
for (state in states) {
  index <- which(metadata$cell_type == state)
  local_metadata <- droplevels(metadata[index, , drop = FALSE])
  local_metadata$route <- relevel(factor(local_metadata$route), ref = "normal")
  local_metadata$donor_id <- factor(local_metadata$donor_id)
  dge_reference <- DGEList(
    counts = as.matrix(pseudobulk$counts[shared_genes, index, drop = FALSE])
  )
  dge_reference <- calcNormFactors(dge_reference, method = "TMM")
  dge_panel <- DGEList(
    counts = as.matrix(pseudobulk$counts[panel_genes, index, drop = FALSE])
  )
  dge_panel$samples$lib.size <- dge_reference$samples$lib.size
  dge_panel$samples$norm.factors <- dge_reference$samples$norm.factors
  log_cpm <- t(cpm(dge_panel, log = TRUE, prior.count = 0.5))

  for (gene in panel_genes) {
    local_metadata$expression <- log_cpm[, gene]
    fit <- lmerTest::lmer(
      expression ~ route + (1 | donor_id),
      data = local_metadata,
      REML = FALSE
    )
    coefficient <- summary(fit)$coefficients[
      "routeconventional_adenoma",
      ,
      drop = TRUE
    ]
    count <- as.numeric(pseudobulk$counts[gene, index])
    record_index <- record_index + 1L
    records[[record_index]] <- data.frame(
      gene = gene,
      arm = panel$arm[match(gene, panel$gene)],
      pair_step = panel$pair_step[match(gene, panel$gene)],
      cell_type = state,
      n_profiles = length(index),
      n_donors = nlevels(local_metadata$donor_id),
      normal_detection_fraction = mean(count[local_metadata$route == "normal"] > 0),
      adenoma_detection_fraction = mean(
        count[local_metadata$route == "conventional_adenoma"] > 0
      ),
      normal_median_log_cpm = median(
        local_metadata$expression[local_metadata$route == "normal"]
      ),
      adenoma_median_log_cpm = median(
        local_metadata$expression[
          local_metadata$route == "conventional_adenoma"
        ]
      ),
      log_cpm_route_effect = coefficient["Estimate"],
      standard_error = coefficient["Std. Error"],
      degrees_of_freedom = coefficient["df"],
      t_statistic = coefficient["t value"],
      p_value = coefficient["Pr(>|t|)"],
      singular_fit = lme4::isSingular(fit, tol = 1e-5),
      stringsAsFactors = FALSE
    )

    donor_route <- aggregate(
      local_metadata$expression,
      by = list(
        donor_id = as.character(local_metadata$donor_id),
        route = as.character(local_metadata$route)
      ),
      FUN = mean
    )
    colnames(donor_route)[3L] <- "expression"
    paired <- tidyr::pivot_wider(
      donor_route,
      names_from = route,
      values_from = expression
    )
    paired <- paired[
      is.finite(paired$normal) & is.finite(paired$conventional_adenoma),
      ,
      drop = FALSE
    ]
    difference <- paired$conventional_adenoma - paired$normal
    paired_index <- paired_index + 1L
    paired_records[[paired_index]] <- data.frame(
      gene = gene,
      arm = panel$arm[match(gene, panel$gene)],
      cell_type = state,
      n_paired_donors = length(difference),
      mean_paired_difference = mean(difference),
      median_paired_difference = median(difference),
      paired_direction_matches_arm =
        sign(mean(difference)) == ifelse(
          panel$arm[match(gene, panel$gene)] == "up",
          1,
          -1
        ),
      stringsAsFactors = FALSE
    )
  }
}
targeted <- do.call(rbind, records)
targeted$p_value_BH <- p.adjust(targeted$p_value, method = "BH")
paired <- do.call(rbind, paired_records)

null_correlation <- as.matrix(null_table[, states, drop = FALSE])
rownames(null_correlation) <- null_table$state
null_correlation <- null_correlation[states, states, drop = FALSE]
ones <- rep(1, length(states))
common_records <- lapply(panel_genes, function(gene) {
  local <- targeted[targeted$gene == gene, , drop = FALSE]
  local <- local[match(states, local$cell_type), , drop = FALSE]
  effect <- local$log_cpm_route_effect
  standard_error <- local$standard_error
  covariance <- diag(standard_error) %*%
    null_correlation %*%
    diag(standard_error)
  precision <- solve(covariance)
  denominator <- as.numeric(t(ones) %*% precision %*% ones)
  common_effect <- as.numeric(
    t(ones) %*% precision %*% effect / denominator
  )
  common_se <- sqrt(1 / denominator)
  data.frame(
    gene = gene,
    arm = panel$arm[match(gene, panel$gene)],
    pair_step = panel$pair_step[match(gene, panel$gene)],
    common_effect = common_effect,
    common_se = common_se,
    common_z = common_effect / common_se,
    common_p_value = 2 * pnorm(-abs(common_effect / common_se)),
    common_direction_matches_arm = sign(common_effect) == ifelse(
      panel$arm[match(gene, panel$gene)] == "up",
      1,
      -1
    ),
    all_state_directions_match_arm = all(
      sign(effect) == ifelse(
        panel$arm[match(gene, panel$gene)] == "up",
        1,
        -1
      )
    ),
    stringsAsFactors = FALSE
  )
})
common <- do.call(rbind, common_records)
common$common_q_value <- p.adjust(common$common_p_value, method = "BH")

write.table(
  targeted,
  file.path(out_dir, "heldout_compact_panel_targeted_state_effects.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  paired,
  file.path(out_dir, "heldout_compact_panel_targeted_paired_effects.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  common,
  file.path(out_dir, "heldout_compact_panel_targeted_common_effects.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_targeted_compact_gene_validation_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  frozen_input_hashes = as.list(actual_hashes),
  genes = panel_genes,
  normalisation = "validation shared-universe TMM factors",
  expression = "log2 CPM with prior.count 0.5",
  targeted_model = "logCPM ~ route + (1 | donor_id); lmerTest Satterthwaite",
  use = "secondary targeted audit for frozen genes, including genes outside the genome-wide shared expression filter",
  validation_outcomes_used_for_panel_revision = FALSE,
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "heldout_compact_panel_targeted_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Targeted compact-gene validation completed; common direction matches=",
  sum(common$common_direction_matches_arm),
  "/",
  nrow(common)
)
