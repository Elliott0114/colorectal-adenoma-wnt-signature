#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
input_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "purity_audit",
  "specimen_state_purity_metrics.tsv"
)
out_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "purity_audit",
  "paired_donor_purity_tests.tsv"
)
if (!file.exists(input_path)) {
  stop("Specimen-state purity metrics are missing")
}

data <- read.delim(input_path, check.names = FALSE)
metrics <- c(
  "epithelial_positive_fraction",
  "immune_positive_fraction",
  "stromal_positive_fraction",
  "endothelial_positive_fraction",
  "any_nonepithelial_positive_fraction",
  "INPP5D_positive_fraction",
  "TIMP1_positive_fraction",
  "S100P_positive_fraction",
  "LCN2_positive_fraction"
)
states <- c("ABS", "GOB", "TAC")

rows <- list()
counter <- 1L
for (state in states) {
  state_data <- data[data$cell_type == state, , drop = FALSE]
  for (metric in metrics) {
    donor_route <- aggregate(
      state_data[[metric]],
      by = list(donor_id = state_data$donor_id, route = state_data$route),
      FUN = median,
      na.rm = TRUE
    )
    colnames(donor_route)[3L] <- "value"
    wide <- reshape(
      donor_route,
      idvar = "donor_id",
      timevar = "route",
      direction = "wide"
    )
    required <- c("value.normal", "value.conventional_adenoma")
    if (!all(required %in% colnames(wide))) {
      next
    }
    paired <- wide[complete.cases(wide[required]), , drop = FALSE]
    difference <-
      paired$value.conventional_adenoma - paired$value.normal
    p_value <- if (length(difference) >= 3L && any(difference != 0)) {
      suppressWarnings(wilcox.test(
        paired$value.conventional_adenoma,
        paired$value.normal,
        paired = TRUE,
        exact = FALSE,
        alternative = "two.sided"
      )$p.value)
    } else {
      NA_real_
    }
    rows[[counter]] <- data.frame(
      cell_type = state,
      metric = metric,
      n_paired_donors = nrow(paired),
      median_normal = median(paired$value.normal, na.rm = TRUE),
      median_adenoma = median(
        paired$value.conventional_adenoma,
        na.rm = TRUE
      ),
      median_paired_difference = median(difference, na.rm = TRUE),
      positive_difference_fraction = mean(difference > 0, na.rm = TRUE),
      p_paired_wilcoxon = p_value,
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}
result <- do.call(rbind, rows)
result$q_value <- p.adjust(result$p_paired_wilcoxon, method = "BH")
write.table(
  result,
  out_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message("Paired-donor purity tests written for ", nrow(result), " state-metric pairs")
