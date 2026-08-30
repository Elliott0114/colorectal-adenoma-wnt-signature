#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260830)

root <- normalizePath(".", mustWork = TRUE)
fine_dir <- file.path(root, "results", "state_shared_revision_v2", "fine_states")
out_dir <- file.path(root, "results", "state_shared_revision_v2", "fine_state_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  assignments = file.path(fine_dir, "cell_fine_state_assignments.tsv.gz"),
  fine_state_manifest = file.path(fine_dir, "analysis_manifest.json"),
  contract = file.path(
    root,
    "analysis",
    "contracts",
    "state_shared_revision_validation_v2_2026-08-30.md"
  )
)
if (!all(file.exists(unlist(paths)))) {
  stop("Run the programme-blind fine-state assignment before modelling")
}

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
input_hashes <- vapply(paths, sha256, character(1))
cells <- read.delim(paths$assignments, check.names = FALSE)

states <- c("ABS", "GOB", "TAC")
k_values <- c(3L, 4L, 5L)
minimum_cells <- 20L
n_bootstrap <- 5000L

proximal_sites <- c(
  "ascending colon",
  "transverse colon",
  "hepatic cecum",
  "hepatic flexure of colon"
)
distal_sites <- c("descending colon", "sigmoid colon", "rectum")
cells$site_group <- ifelse(
  cells$tissue %in% proximal_sites,
  "proximal",
  ifelse(cells$tissue %in% distal_sites, "distal", NA_character_)
)
if (anyNA(cells$site_group)) {
  stop("An anatomical site could not be mapped to proximal or distal colon")
}

summarise_strata <- function(data, k) {
  fine_column <- paste0("fine_state_k", k)
  output <- aggregate(
    data$programme_rank_score,
    by = list(
      partition = data$partition,
      donor_id = data$donor_id,
      specimen_id = data$specimen_id,
      route = data$route,
      tissue = data$tissue,
      site_group = data$site_group,
      broad_state = data$broad_state,
      fine_state = data[[fine_column]]
    ),
    FUN = function(value) c(n_cells = length(value), mean_score = mean(value))
  )
  values <- if (is.matrix(output$x)) {
    output$x
  } else {
    do.call(rbind, output$x)
  }
  output$x <- NULL
  output$n_cells <- values[, "n_cells"]
  output$mean_score <- values[, "mean_score"]
  output$k <- k
  output
}

strata <- do.call(rbind, lapply(k_values, function(k) summarise_strata(cells, k)))
eligible_strata <- strata[strata$n_cells >= minimum_cells, , drop = FALSE]

fit_fine_state_effect <- function(data, k, scope, weighted = FALSE) {
  local <- data[
    data$partition == "validation" & data$k == k,
    ,
    drop = FALSE
  ]
  if (scope != "all") {
    local <- local[local$broad_state == scope, , drop = FALSE]
  }
  local$route <- relevel(factor(local$route), ref = "normal")
  local$site_group <- relevel(factor(local$site_group), ref = "distal")
  local$donor_id <- factor(local$donor_id)
  local$specimen_id <- factor(local$specimen_id)
  local$fine_state <- factor(local$fine_state)
  local$broad_state <- factor(local$broad_state, levels = states)
  fixed <- if (scope == "all") {
    "route + broad_state:fine_state + site_group"
  } else {
    "route + fine_state + site_group"
  }
  formula <- as.formula(paste(
    "mean_score ~",
    fixed,
    "+ (1 | donor_id) + (1 | specimen_id)"
  ))
  fit <- if (weighted) {
    lmerTest::lmer(
      formula,
      data = local,
      weights = n_cells,
      REML = FALSE
    )
  } else {
    lmerTest::lmer(formula, data = local, REML = FALSE)
  }
  coefficient <- summary(fit)$coefficients["routeconventional_adenoma", ]
  critical <- qt(0.975, coefficient["df"])
  data.frame(
    k = k,
    scope = scope,
    model = ifelse(weighted, "cell_count_weighted", "unweighted"),
    n_strata = nrow(local),
    n_specimens = length(unique(local$specimen_id)),
    n_donors = length(unique(local$donor_id)),
    n_adenoma_strata = sum(local$route == "conventional_adenoma"),
    n_normal_strata = sum(local$route == "normal"),
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

fine_state_effects <- do.call(rbind, lapply(k_values, function(k) {
  do.call(rbind, lapply(c("all", states), function(scope) {
    rbind(
      fit_fine_state_effect(eligible_strata, k, scope, FALSE),
      fit_fine_state_effect(eligible_strata, k, scope, TRUE)
    )
  }))
}))
fine_state_effects$p_value_BH <- p.adjust(
  fine_state_effects$p_value,
  method = "BH"
)

make_donor_strata <- function(data, partition, k, scope, broad_only = FALSE) {
  local <- data[data$partition == partition, , drop = FALSE]
  if (scope != "all") {
    local <- local[local$broad_state == scope, , drop = FALSE]
  }
  fine_column <- paste0("fine_state_k", k)
  local$stratum <- if (broad_only) {
    local$broad_state
  } else if (scope == "all") {
    paste(local$broad_state, local[[fine_column]], sep = "::")
  } else {
    local[[fine_column]]
  }
  summary <- aggregate(
    local$programme_rank_score,
    by = list(
      donor_id = local$donor_id,
      route = local$route,
      stratum = local$stratum
    ),
    FUN = function(value) c(n_cells = length(value), score_sum = sum(value))
  )
  values <- if (is.matrix(summary$x)) {
    summary$x
  } else {
    do.call(rbind, summary$x)
  }
  summary$x <- NULL
  summary$n_cells <- values[, "n_cells"]
  summary$score_sum <- values[, "score_sum"]
  totals <- aggregate(n_cells ~ donor_id + route, summary, sum)
  names(totals)[3L] <- "donor_route_cells"
  summary <- merge(summary, totals, by = c("donor_id", "route"), all.x = TRUE)
  summary$proportion <- summary$n_cells / summary$donor_route_cells
  summary$mean_score <- summary$score_sum / summary$n_cells
  summary
}

decompose_once <- function(donor_strata, donor_weights) {
  all_strata <- sort(unique(donor_strata$stratum))
  group_values <- list()
  for (route_value in c("normal", "conventional_adenoma")) {
    local <- donor_strata[donor_strata$route == route_value, , drop = FALSE]
    local$weight <- donor_weights[match(local$donor_id, names(donor_weights))]
    route_donors <- unique(local$donor_id)
    denominator <- sum(donor_weights[route_donors])
    values <- data.frame(
      stratum = all_strata,
      proportion = 0,
      contribution = 0,
      stringsAsFactors = FALSE
    )
    for (index in seq_along(all_strata)) {
      stratum_value <- all_strata[index]
      current <- local[local$stratum == stratum_value, , drop = FALSE]
      if (nrow(current)) {
        values$proportion[index] <- sum(current$weight * current$proportion) /
          denominator
        values$contribution[index] <- sum(
          current$weight * current$proportion * current$mean_score
        ) / denominator
      }
    }
    values$mean_score <- ifelse(
      values$proportion > 0,
      values$contribution / values$proportion,
      0
    )
    group_values[[route_value]] <- values
  }
  normal <- group_values$normal
  adenoma <- group_values$conventional_adenoma
  total <- sum(adenoma$contribution) - sum(normal$contribution)
  composition <- sum(
    (adenoma$proportion - normal$proportion) *
      (adenoma$mean_score + normal$mean_score) / 2
  )
  within <- sum(
    (adenoma$mean_score - normal$mean_score) *
      (adenoma$proportion + normal$proportion) / 2
  )
  c(
    total = total,
    composition = composition,
    within = within,
    reconstruction_error = total - composition - within
  )
}

run_decomposition <- function(partition, k, scope, broad_only = FALSE) {
  donor_strata <- make_donor_strata(cells, partition, k, scope, broad_only)
  donors <- sort(unique(donor_strata$donor_id))
  observed_weights <- stats::setNames(rep(1, length(donors)), donors)
  observed <- decompose_once(donor_strata, observed_weights)
  bootstrap <- t(vapply(seq_len(n_bootstrap), function(iteration) {
    sampled <- sample(donors, length(donors), replace = TRUE)
    weights <- table(factor(sampled, levels = donors))
    names(weights) <- donors
    decompose_once(donor_strata, as.numeric(weights) |> stats::setNames(donors))
  }, numeric(4L)))
  decomposition_type <- ifelse(broad_only, "broad_state", "fine_state")
  summary <- data.frame(
    partition = partition,
    k = ifelse(broad_only, NA_integer_, k),
    scope = scope,
    decomposition_type = decomposition_type,
    component = names(observed),
    estimate = as.numeric(observed),
    ci_low = apply(bootstrap, 2L, quantile, 0.025, na.rm = TRUE),
    ci_high = apply(bootstrap, 2L, quantile, 0.975, na.rm = TRUE),
    n_donors = length(donors),
    n_normal_donors = length(unique(
      donor_strata$donor_id[donor_strata$route == "normal"]
    )),
    n_adenoma_donors = length(unique(
      donor_strata$donor_id[donor_strata$route == "conventional_adenoma"]
    )),
    stringsAsFactors = FALSE
  )
  bootstrap_frame <- data.frame(
    partition = partition,
    k = ifelse(broad_only, NA_integer_, k),
    scope = scope,
    decomposition_type = decomposition_type,
    iteration = seq_len(n_bootstrap),
    bootstrap,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(summary = summary, bootstrap = bootstrap_frame)
}

decomposition_records <- list()
bootstrap_records <- list()
record_index <- 0L
for (partition in c("discovery", "validation")) {
  for (k in k_values) {
    for (scope in c("all", states)) {
      record_index <- record_index + 1L
      result <- run_decomposition(partition, k, scope, FALSE)
      decomposition_records[[record_index]] <- result$summary
      bootstrap_records[[record_index]] <- result$bootstrap
    }
  }
  record_index <- record_index + 1L
  result <- run_decomposition(partition, 4L, "all", TRUE)
  decomposition_records[[record_index]] <- result$summary
  bootstrap_records[[record_index]] <- result$bootstrap
}
decomposition <- do.call(rbind, decomposition_records)
decomposition_bootstrap <- do.call(rbind, bootstrap_records)

primary_decomposition <- decomposition[
  decomposition$partition == "validation" &
    decomposition$k == 4L &
    decomposition$scope == "all" &
    decomposition$decomposition_type == "fine_state",
  ,
  drop = FALSE
]
component_value <- stats::setNames(
  primary_decomposition$estimate,
  primary_decomposition$component
)
within_fraction <- component_value["within"] / component_value["total"]
gates <- data.frame(
  gate = c(
    "fine_state_adjusted_effect_positive_at_k3_k4_k5",
    "fine_state_adjusted_effect_positive_within_all_broad_states_at_k4",
    "primary_within_component_positive",
    "primary_within_component_larger_than_composition",
    "primary_within_fraction_exceeds_half",
    "decomposition_reconstruction_error_below_1e-10"
  ),
  passed = c(
    all(fine_state_effects$estimate[
      fine_state_effects$scope == "all" &
        fine_state_effects$model == "unweighted"
    ] > 0),
    all(fine_state_effects$estimate[
      fine_state_effects$k == 4L &
        fine_state_effects$scope %in% states &
        fine_state_effects$model == "unweighted"
    ] > 0),
    component_value["within"] > 0,
    component_value["within"] > component_value["composition"],
    within_fraction > 0.5,
    max(abs(decomposition$estimate[
      decomposition$component == "reconstruction_error"
    ])) < 1e-10
  ),
  stringsAsFactors = FALSE
)

fine_state_inventory <- aggregate(
  n_cells ~ partition + k + broad_state + fine_state + route,
  strata,
  sum
)
donor_inventory <- aggregate(
  donor_id ~ partition + k + broad_state + fine_state + route,
  strata,
  function(value) length(unique(value))
)
names(donor_inventory)[6L] <- "n_donors"
fine_state_inventory <- merge(
  fine_state_inventory,
  donor_inventory,
  by = c("partition", "k", "broad_state", "fine_state", "route")
)

write.table(
  strata,
  gzfile(file.path(out_dir, "specimen_fine_state_scores.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  fine_state_effects,
  file.path(out_dir, "fine_state_adjusted_route_effects.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  decomposition,
  file.path(out_dir, "programme_composition_decomposition.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  decomposition_bootstrap,
  gzfile(file.path(out_dir, "programme_composition_decomposition_bootstrap.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  fine_state_inventory,
  file.path(out_dir, "fine_state_inventory.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  gates,
  file.path(out_dir, "quality_gates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- list(
  analysis = "state_shared_revision_fine_state_models_v2",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = 20260830,
  input_sha256 = as.list(input_hashes),
  biological_replicate = "donor",
  cell_role = "measurement unit only",
  pseudobulk_stratum = "specimen by broad state by programme-blind fine state",
  minimum_cells_per_stratum = minimum_cells,
  primary_fine_state_resolution = 4L,
  sensitivity_resolutions = c(3L, 5L),
  decomposition_bootstrap_replicates = n_bootstrap,
  primary_within_fraction = unname(within_fraction),
  quality_gates = stats::setNames(as.list(gates$passed), gates$gate),
  package_versions = as.list(vapply(
    c("lme4", "lmerTest"),
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
  "Fine-state modelling completed: primary total=",
  sprintf("%.4f", component_value["total"]),
  "; composition=",
  sprintf("%.4f", component_value["composition"]),
  "; within=",
  sprintf("%.4f", component_value["within"]),
  "; within fraction=",
  sprintf("%.3f", within_fraction)
)
