#!/usr/bin/env Rscript

# Independent base-R audit of the primary Kitagawa identity from saved inputs.

file_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("^--file=", "", file_arg[1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input_path <- file.path(root, "results", "cell_state_decomposition_v1", "donor_route_decomposition_inputs.tsv.gz")
reported_path <- file.path(root, "results", "cell_state_decomposition_v1", "decomposition_summary.tsv")
output_path <- file.path(root, "results", "cell_state_decomposition_v1", "independent_r_audit.tsv")

inputs <- read.delim(gzfile(input_path), check.names = FALSE, stringsAsFactors = FALSE)
reported <- read.delim(reported_path, check.names = FALSE, stringsAsFactors = FALSE)

audit_one <- function(dataset, score, state_set, comparison, group_a, group_b) {
  # Avoid NSE name collisions in subset() by applying explicit indices.
  x <- inputs[
    inputs$dataset == dataset & inputs$score == score & inputs$state_set == state_set &
      inputs$route_group %in% c(group_a, group_b),
  ]
  group_state <- aggregate(cbind(proportion, contribution) ~ route_group + cell_type, x, mean)
  p <- reshape(group_state[c("route_group", "cell_type", "proportion")], idvar = "cell_type", timevar = "route_group", direction = "wide")
  ctab <- reshape(group_state[c("route_group", "cell_type", "contribution")], idvar = "cell_type", timevar = "route_group", direction = "wide")
  merged <- merge(p, ctab, by = "cell_type", all = TRUE)
  p_a <- merged[[paste0("proportion.", group_a)]]
  p_b <- merged[[paste0("proportion.", group_b)]]
  c_a <- merged[[paste0("contribution.", group_a)]]
  c_b <- merged[[paste0("contribution.", group_b)]]
  mu_a <- c_a / p_a
  mu_b <- c_b / p_b
  total <- sum(c_a) - sum(c_b)
  composition <- sum((p_a - p_b) * (mu_a + mu_b) / 2)
  within_state <- sum((mu_a - mu_b) * (p_a + p_b) / 2)
  data.frame(
    dataset = dataset,
    score = score,
    state_set = state_set,
    comparison = comparison,
    total_r = total,
    composition_r = composition,
    within_state_r = within_state,
    closure_error_r = total - composition - within_state,
    stringsAsFactors = FALSE
  )
}

rows <- lapply(seq_len(nrow(reported)), function(i) {
  audit_one(
    reported$dataset[i], reported$score[i], reported$state_set[i], reported$comparison[i],
    reported$group_a[i], reported$group_b[i]
  )
})
audit <- do.call(rbind, rows)
audit <- merge(
  audit,
  reported[c("dataset", "score", "state_set", "comparison", "total_difference", "composition_component", "within_state_component")],
  by = c("dataset", "score", "state_set", "comparison"),
  all.x = TRUE
)
audit$total_abs_error <- abs(audit$total_r - audit$total_difference)
audit$composition_abs_error <- abs(audit$composition_r - audit$composition_component)
audit$within_state_abs_error <- abs(audit$within_state_r - audit$within_state_component)
audit$pass <- apply(
  audit[c("closure_error_r", "total_abs_error", "composition_abs_error", "within_state_abs_error")],
  1,
  function(values) all(abs(values) <= 1e-10)
)

write.table(audit, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
if (!all(audit$pass)) {
  stop("Independent R decomposition audit failed")
}
cat(sprintf("Independent R audit passed for %d decompositions; maximum absolute error %.3e\n", nrow(audit), max(abs(unlist(audit[c("closure_error_r", "total_abs_error", "composition_abs_error", "within_state_abs_error")])))))
