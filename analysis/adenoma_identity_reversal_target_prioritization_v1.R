#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(decoupleR)
  library(digest)
  library(edgeR)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)
set.seed(20260830)

root <- normalizePath(".", mustWork = TRUE)
contract_path <- file.path(
  root,
  "analysis",
  "contracts",
  "adenoma_identity_reversal_target_prioritization_v1_2026-08-30.md"
)
common_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "common_effects",
  "cross_state_common_effects.tsv.gz"
)
heldout_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "heldout_validation",
  "heldout_cross_state_common_effects.tsv.gz"
)
pseudobulk_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "discovery_pseudobulk",
  "discovery_state_pseudobulk.rds"
)
collectri_path <- file.path(
  root,
  "data_sources",
  "regulatory_priors_v2",
  "interactions_collectri_9606_2023-09-17.tsv.gz"
)
progeny_path <- file.path(
  root,
  "data_sources",
  "regulatory_priors_v2",
  "annotations_PROGENy_9606_2023-09-17.tsv.gz"
)
perturbation_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "extended_validation_full_programme",
  "perturbation_spatial",
  "perturbation_full_programme_unit_effects.tsv"
)
external_perturbation_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "external_validation",
  "perturbation_sample_scores.tsv"
)
out_dir <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "identity_reversal_target_prioritization_v1"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_hashes <- c(
  contract = "08b9824cdf8605e14394cd0162393673fd383248db98d89cd2bbc4c88d97284a",
  common = "a1ac4b7b67ac279782e04e386971d7463e169cbdbc24d7da4c314e34a4f3e946",
  heldout = "ab9918f1481fcb52923133cafc64c59d16eabe3995bf46f7714da742acf95fc5",
  pseudobulk = "2c63422dd5cf8124098974fec62648bb24b546adb407ec41d794c8709cdb4c96",
  collectri = "193e03a026fb80bcdaef2513fe93d40d83ff4e4309be83a2768292ee81269fe9",
  progeny = "fd9f1b8c25b98c61412baef2f226a652423eb580e37712ffb7b93141cbcc91e7",
  perturbation = "4653f1a5c2f4d28621981994a6d4763643afb760388f5171e247df75f0715eed",
  external_perturbation = "24efafa1d030ca34e03b700e779ccbda13e692f0bd16b46d23e3163159cc7c09"
)
input_paths <- c(
  contract = contract_path,
  common = common_path,
  heldout = heldout_path,
  pseudobulk = pseudobulk_path,
  collectri = collectri_path,
  progeny = progeny_path,
  perturbation = perturbation_path,
  external_perturbation = external_perturbation_path
)
observed_hashes <- vapply(
  input_paths,
  digest,
  character(1),
  algo = "sha256",
  file = TRUE
)
if (!identical(unname(observed_hashes), unname(expected_hashes))) {
  mismatch <- names(expected_hashes)[observed_hashes != expected_hashes]
  stop("Frozen-input hash mismatch: ", paste(mismatch, collapse = ", "))
}

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

write_tsv <- function(x, path, compress = FALSE) {
  target <- file.path(out_dir, path)
  if (compress) {
    connection <- gzfile(target, open = "wt")
    on.exit(close(connection), add = TRUE)
    write.table(x, connection, sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    write.table(x, target, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

message("Reading frozen discovery and donor-disjoint rankings")
common <- read.delim(common_path, check.names = FALSE)
heldout <- read.delim(heldout_path, check.names = FALSE)
common$strict_state_shared <- as_flag(common$strict_state_shared)
programme <- common[common$strict_state_shared, c("gene", "shared_direction")]
programme$programme_sign <- ifelse(programme$shared_direction == "up", 1, -1)
if (nrow(programme) != 1843L || sum(programme$programme_sign == 1) != 884L) {
  stop("Frozen 1,843-gene programme does not match the contract")
}
universe <- common$gene

message("Preparing the complete signed CollecTRI prior")
collectri <- read.delim(
  collectri_path,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)
stim <- suppressWarnings(as.integer(collectri$is_stimulation))
inhib <- suppressWarnings(as.integer(collectri$is_inhibition))
source_symbol <- collectri$source_genesymbol
is_complex <- grepl("^COMPLEX:", collectri$source)
source_symbol[is_complex & grepl("JUN|FOS", source_symbol)] <- "AP1"
source_symbol[is_complex & grepl("REL|NFKB", source_symbol)] <- "NFKB"
mor <- ifelse(stim == 1L & inhib == 0L, 1, ifelse(inhib == 1L & stim == 0L, -1, NA))
network <- data.frame(
  source = source_symbol,
  target = collectri$target_genesymbol,
  mor = mor,
  stringsAsFactors = FALSE
)
network <- network[
  nzchar(network$source) & nzchar(network$target) & is.finite(network$mor),
]
network <- aggregate(mor ~ source + target, data = network, FUN = mean)
network <- network[network$mor != 0, ]

run_ulm_context <- function(statistics, context, prior, minsize = 20L) {
  statistics <- statistics[is.finite(statistics)]
  if (is.null(names(statistics)) || anyDuplicated(names(statistics))) {
    stop("Invalid named statistic for ", context)
  }
  context_network <- prior[prior$target %in% names(statistics), ]
  result <- run_ulm(
    mat = matrix(
      statistics,
      ncol = 1L,
      dimnames = list(names(statistics), context)
    ),
    network = context_network,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = minsize
  )
  result$q_value <- p.adjust(result$p_value, method = "BH")
  target_counts <- table(context_network$source)
  result$n_targets_measured <- as.integer(target_counts[result$source])
  result$context <- context
  result[, c(
    "context",
    "source",
    "score",
    "p_value",
    "q_value",
    "n_targets_measured"
  )]
}

state_z <- function(frame, state) {
  estimate <- frame[[paste0("logFC_", state)]]
  standard_error <- frame[[paste0("raw_se_", state)]]
  value <- estimate / standard_error
  names(value) <- frame$gene
  value
}

context_statistics <- list(
  discovery_common = setNames(common$common_z, common$gene),
  heldout_common = setNames(heldout$common_z, heldout$gene),
  discovery_ABS = state_z(common, "ABS"),
  discovery_GOB = state_z(common, "GOB"),
  discovery_TAC = state_z(common, "TAC"),
  heldout_ABS = state_z(heldout, "ABS"),
  heldout_GOB = state_z(heldout, "GOB"),
  heldout_TAC = state_z(heldout, "TAC")
)
regulator_context <- do.call(
  rbind,
  Map(
    function(statistics, context) run_ulm_context(statistics, context, network),
    context_statistics,
    names(context_statistics)
  )
)
row.names(regulator_context) <- NULL
write_tsv(regulator_context, "collectri_complete_context_activity.tsv")

context_wide <- lapply(split(regulator_context, regulator_context$context), function(x) {
  context <- unique(x$context)
  result <- x[, c("source", "score", "p_value", "q_value", "n_targets_measured")]
  colnames(result)[-1L] <- paste0(colnames(result)[-1L], "__", context)
  result
})
regulators <- Reduce(function(x, y) merge(x, y, by = "source", all = TRUE), context_wide)

disc_score <- regulators$score__discovery_common
held_score <- regulators$score__heldout_common
regulators$activity_direction <- ifelse(disc_score > 0, "activated_in_adenoma", "suppressed_in_adenoma")
regulators$proposed_intervention <- ifelse(disc_score > 0, "inhibit", "restore_or_activate")
regulators$replicated <-
  regulators$q_value__discovery_common <= 0.05 &
  regulators$q_value__heldout_common <= 0.10 &
  sign(disc_score) == sign(held_score)
disc_state_scores <- as.matrix(regulators[, c(
  "score__discovery_ABS",
  "score__discovery_GOB",
  "score__discovery_TAC"
)])
held_state_scores <- as.matrix(regulators[, c(
  "score__heldout_ABS",
  "score__heldout_GOB",
  "score__heldout_TAC"
)])
regulators$discovery_state_consistent <- apply(
  disc_state_scores,
  1L,
  function(x) all(is.finite(x)) && all(sign(x) == sign(x[1L]))
)
regulators$heldout_state_consistent <- apply(
  held_state_scores,
  1L,
  function(x) all(is.finite(x)) && all(sign(x) == sign(x[1L]))
)

message("Projecting signed regulons onto both frozen programme arms")
programme_sign <- setNames(programme$programme_sign, programme$gene)
projection_rows <- vector("list", nrow(regulators))
edge_rows <- vector("list", nrow(regulators))
for (i in seq_len(nrow(regulators))) {
  regulator <- regulators$source[i]
  activity_sign <- sign(regulators$score__discovery_common[i])
  edges <- network[network$source == regulator & network$target %in% universe, ]
  edges$predicted_disease_sign <- activity_sign * sign(edges$mor)
  edges$programme_sign <- unname(programme_sign[edges$target])
  edges$in_programme <- is.finite(edges$programme_sign)
  edges$concordant <- edges$in_programme &
    edges$predicted_disease_sign == edges$programme_sign
  edges$regulator <- regulator
  edge_rows[[i]] <- edges[, c(
    "regulator",
    "target",
    "mor",
    "predicted_disease_sign",
    "programme_sign",
    "in_programme",
    "concordant"
  )]

  predicted_up <- unique(edges$target[edges$predicted_disease_sign > 0])
  predicted_down <- unique(edges$target[edges$predicted_disease_sign < 0])
  programme_up <- programme$gene[programme$programme_sign > 0]
  programme_down <- programme$gene[programme$programme_sign < 0]
  n_concordant_up <- length(intersect(predicted_up, programme_up))
  n_concordant_down <- length(intersect(predicted_down, programme_down))
  n_programme_targets <- sum(edges$in_programme)
  n_concordant <- sum(edges$concordant)

  enrichment_test <- function(predicted, selected) {
    overlap <- length(intersect(predicted, selected))
    matrix_values <- matrix(
      c(
        overlap,
        length(predicted) - overlap,
        length(selected) - overlap,
        length(universe) - length(union(predicted, selected))
      ),
      nrow = 2L
    )
    result <- fisher.test(matrix_values, alternative = "greater")
    c(odds_ratio = unname(result$estimate), p_value = result$p.value)
  }
  up_test <- enrichment_test(predicted_up, programme_up)
  down_test <- enrichment_test(predicted_down, programme_down)
  any_test <- enrichment_test(unique(edges$target), programme$gene)
  projection_rows[[i]] <- data.frame(
    source = regulator,
    n_targets_universe = nrow(edges),
    n_targets_programme = n_programme_targets,
    n_concordant_programme = n_concordant,
    signed_precision = ifelse(n_programme_targets > 0, n_concordant / n_programme_targets, NA),
    n_concordant_up = n_concordant_up,
    n_concordant_down = n_concordant_down,
    balanced_concordant_count = min(n_concordant_up, n_concordant_down),
    up_arm_coverage = n_concordant_up / sum(programme$programme_sign > 0),
    down_arm_coverage = n_concordant_down / sum(programme$programme_sign < 0),
    balanced_arm_coverage = min(
      n_concordant_up / sum(programme$programme_sign > 0),
      n_concordant_down / sum(programme$programme_sign < 0)
    ),
    up_enrichment_odds_ratio = up_test["odds_ratio"],
    up_enrichment_p_value = up_test["p_value"],
    down_enrichment_odds_ratio = down_test["odds_ratio"],
    down_enrichment_p_value = down_test["p_value"],
    programme_target_enrichment_odds_ratio = any_test["odds_ratio"],
    programme_target_enrichment_p_value = any_test["p_value"],
    stringsAsFactors = FALSE
  )
}
projection <- do.call(rbind, projection_rows)
signed_edges <- do.call(rbind, edge_rows)
projection$up_enrichment_q_value <- p.adjust(projection$up_enrichment_p_value, method = "BH")
projection$down_enrichment_q_value <- p.adjust(projection$down_enrichment_p_value, method = "BH")
projection$programme_target_enrichment_q_value <- p.adjust(
  projection$programme_target_enrichment_p_value,
  method = "BH"
)
regulators <- merge(regulators, projection, by = "source", all.x = TRUE)

message("Calculating expression feasibility in the discovery pseudobulk profiles")
pseudobulk <- readRDS(pseudobulk_path)
dge <- DGEList(counts = pseudobulk$counts)
dge <- calcNormFactors(dge, method = "TMM")
cpm_values <- cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)
route <- as.character(pseudobulk$metadata$route)
expression_rows <- lapply(regulators$source, function(regulator) {
  if (!regulator %in% rownames(cpm_values)) {
    return(data.frame(
      source = regulator,
      mean_cpm_normal = NA_real_,
      mean_cpm_adenoma = NA_real_,
      fraction_adenoma_profiles_cpm_ge_1 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  values <- cpm_values[regulator, ]
  data.frame(
    source = regulator,
    mean_cpm_normal = mean(values[route == "normal"]),
    mean_cpm_adenoma = mean(values[route == "conventional_adenoma"]),
    fraction_adenoma_profiles_cpm_ge_1 = mean(
      values[route == "conventional_adenoma"] >= 1
    ),
    stringsAsFactors = FALSE
  )
})
expression_feasibility <- do.call(rbind, expression_rows)
regulators <- merge(regulators, expression_feasibility, by = "source", all.x = TRUE)

pareto_frontier <- function(frame, columns) {
  if (!nrow(frame)) {
    return(logical())
  }
  values <- as.matrix(frame[, columns, drop = FALSE])
  keep <- complete.cases(values)
  result <- rep(FALSE, nrow(frame))
  for (i in which(keep)) {
    dominated <- FALSE
    for (j in which(keep)) {
      if (i == j) next
      if (all(values[j, ] >= values[i, ]) && any(values[j, ] > values[i, ])) {
        dominated <- TRUE
        break
      }
    }
    result[i] <- !dominated
  }
  result
}
regulators$absolute_discovery_score <- abs(regulators$score__discovery_common)
regulators$absolute_heldout_score <- abs(regulators$score__heldout_common)
regulators$pareto_replicated <- FALSE
replicated_index <- which(regulators$replicated %in% TRUE)
regulators$pareto_replicated[replicated_index] <- pareto_frontier(
  regulators[replicated_index, ],
  c(
    "absolute_discovery_score",
    "absolute_heldout_score",
    "signed_precision",
    "balanced_concordant_count"
  )
)
regulators$single_gene_entity <- grepl("^[A-Za-z0-9.-]+$", regulators$source) &
  !regulators$source %in% c("AP1", "NFKB")
regulators$priority_tier <- ifelse(
  regulators$replicated & regulators$pareto_replicated &
    regulators$discovery_state_consistent &
    regulators$balanced_concordant_count >= 3L &
    regulators$single_gene_entity,
  "A: replicated non-dominated two-arm node",
  ifelse(
    regulators$replicated & regulators$pareto_replicated,
    "B: replicated non-dominated node with a boundary",
    ifelse(regulators$replicated, "C: replicated dominated node", "not replicated")
  )
)
regulators$direct_empirical_perturbation <- ifelse(
  regulators$source == "ASCL2",
  "ASCL2 knockout resource available",
  ifelse(
    regulators$source == "TCF7L2",
    "TCF7L2-edited resource available",
    "none in current perturbation collection"
  )
)
regulators <- regulators[order(
  factor(regulators$priority_tier, levels = c(
    "A: replicated non-dominated two-arm node",
    "B: replicated non-dominated node with a boundary",
    "C: replicated dominated node",
    "not replicated"
  )),
  -regulators$balanced_concordant_count,
  -regulators$absolute_heldout_score
), ]
write_tsv(regulators, "regulator_replication_and_programme_projection.tsv")
write_tsv(signed_edges, "collectri_signed_programme_edges.tsv.gz", compress = TRUE)

shortlist <- regulators[
  regulators$replicated & regulators$pareto_replicated,
  c(
    "source",
    "activity_direction",
    "proposed_intervention",
    "priority_tier",
    "score__discovery_common",
    "q_value__discovery_common",
    "score__heldout_common",
    "q_value__heldout_common",
    "discovery_state_consistent",
    "heldout_state_consistent",
    "n_targets_programme",
    "n_concordant_up",
    "n_concordant_down",
    "balanced_concordant_count",
    "signed_precision",
    "up_enrichment_q_value",
    "down_enrichment_q_value",
    "mean_cpm_adenoma",
    "fraction_adenoma_profiles_cpm_ge_1",
    "direct_empirical_perturbation"
  )
]
write_tsv(shortlist, "regulator_pareto_shortlist.tsv")

message("Inferring PROGENy pathway activities")
progeny_long <- read.delim(
  progeny_path,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)
pathway_values <- progeny_long[
  progeny_long$label == "pathway",
  c("record_id", "genesymbol", "value")
]
colnames(pathway_values)[3L] <- "pathway"
weight_values <- progeny_long[
  progeny_long$label == "weight",
  c("record_id", "genesymbol", "value")
]
colnames(weight_values)[3L] <- "weight"
pvalue_values <- progeny_long[
  progeny_long$label == "p_value",
  c("record_id", "genesymbol", "value")
]
colnames(pvalue_values)[3L] <- "prior_p_value"
progeny_wide <- Reduce(
  function(x, y) merge(x, y, by = c("record_id", "genesymbol")),
  list(pathway_values, weight_values, pvalue_values)
)
progeny_wide$weight <- as.numeric(progeny_wide$weight)
progeny_wide$prior_p_value <- as.numeric(progeny_wide$prior_p_value)
progeny_top <- do.call(rbind, lapply(split(progeny_wide, progeny_wide$pathway), function(x) {
  x <- x[order(x$prior_p_value, -abs(x$weight)), ]
  head(x, 500L)
}))
progeny_network <- data.frame(
  source = progeny_top$pathway,
  target = progeny_top$genesymbol,
  mor = progeny_top$weight,
  stringsAsFactors = FALSE
)
progeny_network <- aggregate(
  mor ~ source + target,
  data = progeny_network,
  FUN = mean
)
progeny_network <- progeny_network[is.finite(progeny_network$mor), ]
progeny_activity <- rbind(
  run_ulm_context(
    setNames(common$common_z, common$gene),
    "discovery_common",
    progeny_network,
    minsize = 20L
  ),
  run_ulm_context(
    setNames(heldout$common_z, heldout$gene),
    "heldout_common",
    progeny_network,
    minsize = 20L
  )
)
write_tsv(progeny_activity, "progeny_discovery_heldout_activity.tsv")

message("Building the empirical two-component perturbation map")
unit_effects <- read.delim(perturbation_path, check.names = FALSE)
unit_wide <- reshape(
  unit_effects[, c(
    "dataset",
    "species",
    "model_system",
    "comparison",
    "unit_id",
    "feature",
    "difference"
  )],
  idvar = c("dataset", "species", "model_system", "comparison", "unit_id"),
  timevar = "feature",
  direction = "wide"
)
colnames(unit_wide) <- sub("^difference\\.", "", colnames(unit_wide))

external_scores <- read.delim(external_perturbation_path, check.names = FALSE)
external_scores <- external_scores[external_scores$signature_id == "state_shared_1843", ]
features <- c("route_up", "route_down", "wnt_stem", "proliferation_control")

extract_condition <- function(frame, values) {
  keep <- rep(TRUE, nrow(frame))
  for (name in names(values)) {
    keep <- keep & as.character(frame[[name]]) == values[[name]]
  }
  frame[keep, , drop = FALSE]
}

gse125 <- external_scores[external_scores$dataset == "GSE125472", ]
gse125 <- aggregate(
  gse125[, features],
  by = gse125[, c("donor_id", "genotype", "wnt_rspo")],
  FUN = mean
)
gse125_rows <- list()
gse125_comparisons <- list(
  list(
    comparison = "wnt_rspo_withdrawal_in_WT",
    target = c(genotype = "WT", wnt_rspo = "without"),
    reference = c(genotype = "WT", wnt_rspo = "with")
  ),
  list(
    comparison = "wnt_rspo_withdrawal_in_APC_KO",
    target = c(genotype = "APC", wnt_rspo = "without"),
    reference = c(genotype = "APC", wnt_rspo = "with")
  )
)
counter <- 1L
for (comparison in gse125_comparisons) {
  target <- extract_condition(gse125, comparison$target)
  reference <- extract_condition(gse125, comparison$reference)
  merged <- merge(target, reference, by = "donor_id", suffixes = c("_target", "_reference"))
  for (i in seq_len(nrow(merged))) {
    row <- data.frame(
      dataset = "GSE125472",
      species = "human",
      model_system = "isogenic_human_colon_organoid",
      comparison = comparison$comparison,
      unit_id = merged$donor_id[i],
      stringsAsFactors = FALSE
    )
    for (feature in features) {
      row[[feature]] <- merged[[paste0(feature, "_target")]][i] -
        merged[[paste0(feature, "_reference")]][i]
    }
    gse125_rows[[counter]] <- row
    counter <- counter + 1L
  }
}
gse125_wide <- do.call(rbind, gse125_rows)

gse135 <- external_scores[external_scores$dataset == "GSE135328", ]
gse135 <- aggregate(
  gse135[, features],
  by = gse135[, c("cell_line", "clone_id", "genotype")],
  FUN = mean
)
gse135_rows <- list()
counter <- 1L
for (cell_line in unique(gse135$cell_line)) {
  cell_frame <- gse135[gse135$cell_line == cell_line, ]
  reference <- cell_frame[cell_frame$genotype == "WT", , drop = FALSE]
  if (nrow(reference) != 1L) next
  targets <- cell_frame[cell_frame$genotype != "WT", , drop = FALSE]
  for (i in seq_len(nrow(targets))) {
    row <- data.frame(
      dataset = "GSE135328",
      species = "human",
      model_system = "crc_cell_line_clone",
      comparison = ifelse(
        targets$genotype[i] == "KO",
        "tcf7l2_KO_vs_WT",
        "tcf7l2_heterozygous_vs_WT"
      ),
      unit_id = targets$clone_id[i],
      stringsAsFactors = FALSE
    )
    for (feature in features) {
      row[[feature]] <- targets[[feature]][i] - reference[[feature]][1L]
    }
    gse135_rows[[counter]] <- row
    counter <- counter + 1L
  }
}
gse135_wide <- do.call(rbind, gse135_rows)

two_axis_units <- rbind(
  unit_wide[, intersect(colnames(unit_wide), colnames(gse125_wide))],
  gse125_wide,
  gse135_wide
)
for (feature in features) {
  if (!feature %in% colnames(two_axis_units)) {
    two_axis_units[[feature]] <- NA_real_
  }
}
two_axis_units$wnt_stem_suppression <- -two_axis_units$wnt_stem
two_axis_units$adenoma_up_arm_suppression <- -two_axis_units$route_up
two_axis_units$mature_function_restoration <- two_axis_units$route_down
two_axis_units$proliferation_suppression <- -two_axis_units$proliferation_control
two_axis_units$conservative_joint_reversal <- pmin(
  two_axis_units$wnt_stem_suppression,
  two_axis_units$mature_function_restoration
)
two_axis_units$quadrant <- ifelse(
  two_axis_units$wnt_stem_suppression > 0 &
    two_axis_units$mature_function_restoration > 0,
  "both axes favourable",
  ifelse(
    two_axis_units$wnt_stem_suppression > 0,
    "WNT/stem only",
    ifelse(
      two_axis_units$mature_function_restoration > 0,
      "mature-function only",
      "neither axis favourable"
    )
  )
)
two_axis_units$interpretation_role <- ifelse(
  grepl("control", two_axis_units$comparison, ignore.case = TRUE),
  "control",
  ifelse(
    two_axis_units$comparison == "tcf7l2_heterozygous_vs_WT",
    "dosage sensitivity",
    "causal or pathway perturbation"
  )
)
write_tsv(two_axis_units, "perturbation_two_component_unit_effects.tsv")

bootstrap_mean <- function(values, seed) {
  values <- values[is.finite(values)]
  if (!length(values)) return(c(mean = NA, low = NA, high = NA))
  if (length(values) == 1L) return(c(mean = values, low = values, high = values))
  set.seed(seed)
  boot <- replicate(5000L, mean(sample(values, replace = TRUE)))
  c(
    mean = mean(values),
    low = unname(quantile(boot, 0.025)),
    high = unname(quantile(boot, 0.975))
  )
}

group_key <- interaction(
  two_axis_units$dataset,
  two_axis_units$species,
  two_axis_units$model_system,
  two_axis_units$comparison,
  two_axis_units$interpretation_role,
  drop = TRUE,
  lex.order = TRUE
)
summary_rows <- lapply(seq_along(split(two_axis_units, group_key)), function(i) {
  frame <- split(two_axis_units, group_key)[[i]]
  x <- bootstrap_mean(frame$wnt_stem_suppression, 20260830L + i)
  y <- bootstrap_mean(frame$mature_function_restoration, 20261830L + i)
  up <- bootstrap_mean(frame$adenoma_up_arm_suppression, 20262830L + i)
  data.frame(
    dataset = frame$dataset[1L],
    species = frame$species[1L],
    model_system = frame$model_system[1L],
    comparison = frame$comparison[1L],
    interpretation_role = frame$interpretation_role[1L],
    n_units = nrow(frame),
    mean_wnt_stem_suppression = x["mean"],
    wnt_stem_ci_low = x["low"],
    wnt_stem_ci_high = x["high"],
    mean_mature_function_restoration = y["mean"],
    mature_function_ci_low = y["low"],
    mature_function_ci_high = y["high"],
    mean_adenoma_up_arm_suppression = up["mean"],
    adenoma_up_arm_ci_low = up["low"],
    adenoma_up_arm_ci_high = up["high"],
    conservative_joint_reversal = min(x["mean"], y["mean"]),
    quadrant = ifelse(
      x["mean"] > 0 & y["mean"] > 0,
      "both axes favourable",
      ifelse(
        x["mean"] > 0,
        "WNT/stem only",
        ifelse(y["mean"] > 0, "mature-function only", "neither axis favourable")
      )
    ),
    stringsAsFactors = FALSE
  )
})
two_axis_summary <- do.call(rbind, summary_rows)
two_axis_summary <- two_axis_summary[order(
  -two_axis_summary$conservative_joint_reversal,
  two_axis_summary$dataset,
  two_axis_summary$comparison
), ]
write_tsv(two_axis_summary, "perturbation_two_component_summary.tsv")

eligible_separability <- two_axis_summary$interpretation_role ==
  "causal or pathway perturbation"
has_wnt_only <- any(
  two_axis_summary$mean_wnt_stem_suppression[eligible_separability] > 0 &
    two_axis_summary$mean_mature_function_restoration[eligible_separability] <= 0
)
has_both <- any(
  two_axis_summary$mean_wnt_stem_suppression[eligible_separability] > 0 &
    two_axis_summary$mean_mature_function_restoration[eligible_separability] > 0
)
separability_supported <- has_wnt_only && has_both

manifest <- list(
  analysis = "adenoma_identity_reversal_target_prioritization_v1",
  analysis_date = "2026-08-30",
  contract_sha256 = observed_hashes[["contract"]],
  input_hashes = as.list(observed_hashes[names(observed_hashes) != "contract"]),
  programme = list(
    n_genes = nrow(programme),
    n_up = sum(programme$programme_sign > 0),
    n_down = sum(programme$programme_sign < 0),
    gene_reselection = FALSE,
    compact_readout_used_for_selection = FALSE
  ),
  regulator_analysis = list(
    method = "decoupleR ULM with complete signed CollecTRI prior",
    minimum_targets = 20,
    n_regulators_tested_discovery = sum(
      regulator_context$context == "discovery_common"
    ),
    n_replicated = sum(regulators$replicated, na.rm = TRUE),
    n_pareto_replicated = sum(regulators$pareto_replicated, na.rm = TRUE),
    n_tier_a = sum(
      regulators$priority_tier == "A: replicated non-dominated two-arm node",
      na.rm = TRUE
    )
  ),
  perturbation_analysis = list(
    primary_x = "negative change in fixed WNT/stem score",
    primary_y = "change in fixed adenoma-down arm score",
    n_comparison_groups = nrow(two_axis_summary),
    has_wnt_only_perturbation = has_wnt_only,
    has_dual_favourable_perturbation = has_both,
    separability_supported = separability_supported
  ),
  software = list(
    R = R.version.string,
    decoupleR = as.character(packageVersion("decoupleR")),
    edgeR = as.character(packageVersion("edgeR"))
  )
)
writeLines(
  toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
  file.path(out_dir, "analysis_manifest.json"),
  useBytes = TRUE
)

summary_lines <- c(
  paste0("regulators_tested\t", manifest$regulator_analysis$n_regulators_tested_discovery),
  paste0("replicated_regulators\t", manifest$regulator_analysis$n_replicated),
  paste0("pareto_replicated_regulators\t", manifest$regulator_analysis$n_pareto_replicated),
  paste0("tier_a_regulators\t", manifest$regulator_analysis$n_tier_a),
  paste0("perturbation_groups\t", nrow(two_axis_summary)),
  paste0("has_wnt_only_perturbation\t", has_wnt_only),
  paste0("has_dual_favourable_perturbation\t", has_both),
  paste0("separability_supported\t", separability_supported)
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.txt"), useBytes = TRUE)
message("Completed identity-reversal target prioritization: ", out_dir)
