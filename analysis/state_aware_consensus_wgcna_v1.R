#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(doParallel)
  library(edgeR)
  library(igraph)
  library(limma)
  library(Matrix)
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
set.seed(20260830)
allowWGCNAThreads(nThreads = 8)

root <- normalizePath(".", mustWork = TRUE)
parent_root <- file.path(root, "results", "state_aware_program_v1")
out_root <- file.path(parent_root, "functional_architecture_v1")
out_dir <- file.path(out_root, "consensus_wgcna")
audit_dir <- file.path(out_dir, "audit")
tom_dir <- file.path(audit_dir, "individual_tom")
dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  discovery_ranking = file.path(
    parent_root, "common_effects", "cross_state_common_effects.tsv.gz"
  ),
  heldout_ranking = file.path(
    parent_root, "heldout_validation", "heldout_cross_state_common_effects.tsv.gz"
  ),
  discovery_pseudobulk = file.path(
    parent_root, "discovery_pseudobulk", "discovery_state_pseudobulk.rds"
  ),
  heldout_pseudobulk = file.path(
    parent_root, "validation_pseudobulk", "validation_state_pseudobulk.rds"
  ),
  pathway_communities = file.path(
    out_root, "pathway_replication", "replicated_pathway_communities.tsv"
  ),
  pathway_leading_edges = file.path(
    out_root, "pathway_replication", "replicated_pathway_leading_edges.tsv"
  ),
  contract = file.path(
    root, "analysis", "contracts",
    "state_aware_functional_architecture_v1_2026-08-30.md"
  )
)
expected_hashes <- c(
  discovery_ranking = "a1ac4b7b67ac279782e04e386971d7463e169cbdbc24d7da4c314e34a4f3e946",
  heldout_ranking = "ab9918f1481fcb52923133cafc64c59d16eabe3995bf46f7714da742acf95fc5",
  discovery_pseudobulk = "2c63422dd5cf8124098974fec62648bb24b546adb407ec41d794c8709cdb4c96",
  heldout_pseudobulk = "3adbe6ad44ef9baf7fd06d1bf6ec2b077483a3399308ed7ce6456186889b6d19",
  contract = "f19df7fe6ffc421477c280d64327c2ba5ffa951ab300a6e26c3be897bc9c3c9b"
)
if (!all(file.exists(unlist(paths)))) {
  stop("At least one consensus-WGCNA input is missing")
}
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
observed_hashes <- vapply(
  paths[names(expected_hashes)], sha256, character(1)
)
if (!identical(unname(observed_hashes), unname(expected_hashes[names(observed_hashes)]))) {
  stop("A frozen consensus-WGCNA input changed")
}

ranking <- read.delim(paths$discovery_ranking, check.names = FALSE)
heldout_ranking <- read.delim(paths$heldout_ranking, check.names = FALSE)
if (nrow(ranking) != 8221L || sum(ranking$strict_state_shared) != 1843L ||
    sum(ranking$strict_state_shared & ranking$shared_direction == "up") != 884L ||
    sum(ranking$strict_state_shared & ranking$shared_direction == "down") != 959L) {
  stop("The testable universe or strict programme changed")
}
universe <- ranking$gene
strict_genes <- ranking$gene[ranking$strict_state_shared]

discovery_object <- readRDS(paths$discovery_pseudobulk)
heldout_object <- readRDS(paths$heldout_pseudobulk)
discovery_metadata <- as.data.frame(discovery_object$metadata)
heldout_metadata_all <- as.data.frame(heldout_object$metadata)
if (!identical(colnames(discovery_object$counts), discovery_metadata$pseudobulk_id) ||
    !identical(colnames(heldout_object$counts), heldout_metadata_all$pseudobulk_id)) {
  stop("Pseudobulk columns and metadata are not aligned")
}

overlap_donors <- intersect(
  unique(as.character(discovery_metadata$donor_id)),
  unique(as.character(heldout_metadata_all$donor_id))
)
heldout_keep <- !heldout_metadata_all$donor_id %in% overlap_donors
heldout_metadata <- heldout_metadata_all[heldout_keep, , drop = FALSE]
heldout_counts <- heldout_object$counts[, heldout_keep, drop = FALSE]
if (length(overlap_donors) != 1L || length(unique(heldout_metadata$donor_id)) != 23L) {
  stop("The frozen donor-disjoint exclusion did not yield 23 held-out donors")
}

missing_discovery <- setdiff(universe, rownames(discovery_object$counts))
if (length(missing_discovery)) {
  stop("Discovery pseudobulk is missing testable-universe genes")
}
network_genes <- universe

aggregate_state_profiles <- function(counts, metadata, genes, partition) {
  state_results <- list()
  sample_audit <- list()
  audit_index <- 0L
  for (state in c("ABS", "GOB", "TAC")) {
    local_meta <- metadata[metadata$cell_type == state, , drop = FALSE]
    local_index <- match(local_meta$pseudobulk_id, colnames(counts))
    local_counts <- counts[genes, local_index, drop = FALSE]
    group_id <- paste(local_meta$donor_id, local_meta$route, sep = "__")
    groups <- unique(group_id)
    aggregate_counts <- matrix(
      0, nrow = length(genes), ncol = length(groups),
      dimnames = list(genes, groups)
    )
    aggregate_meta <- vector("list", length(groups))
    for (group_index in seq_along(groups)) {
      members <- which(group_id == groups[group_index])
      aggregate_counts[, group_index] <- Matrix::rowSums(
        local_counts[, members, drop = FALSE]
      )
      aggregate_meta[[group_index]] <- data.frame(
        profile_id = groups[group_index],
        donor_id = as.character(local_meta$donor_id[members[1L]]),
        route = as.character(local_meta$route[members[1L]]),
        state = state,
        n_specimens = length(members),
        n_cells = sum(local_meta$n_cells[members]),
        library_size_input = sum(local_meta$library_size[members]),
        stringsAsFactors = FALSE
      )
    }
    aggregate_meta <- do.call(rbind, aggregate_meta)
    aggregate_meta$route <- factor(
      aggregate_meta$route,
      levels = c("normal", "conventional_adenoma")
    )
    if (anyNA(aggregate_meta$route)) {
      stop("Unexpected route label in ", partition, " ", state)
    }
    dge <- DGEList(counts = aggregate_counts)
    dge <- calcNormFactors(dge, method = "TMM")
    design <- model.matrix(~ route, data = aggregate_meta)
    voom_object <- voom(dge, design = design, plot = FALSE)
    fit <- lmFit(voom_object, design)
    fitted <- fit$coefficients %*% t(design)
    residual_expression <- voom_object$E - fitted

    donors <- unique(aggregate_meta$donor_id)
    donor_expression <- matrix(
      NA_real_, nrow = length(donors), ncol = length(genes),
      dimnames = list(donors, genes)
    )
    for (donor_index in seq_along(donors)) {
      profile_columns <- aggregate_meta$donor_id == donors[donor_index]
      donor_expression[donor_index, ] <- rowMeans(
        residual_expression[, profile_columns, drop = FALSE]
      )
    }
    profile_expression <- t(voom_object$E)
    rownames(profile_expression) <- aggregate_meta$profile_id
    state_results[[state]] <- list(
      residual_donor_expression = donor_expression,
      profile_expression = profile_expression,
      profile_metadata = aggregate_meta,
      norm_factors = dge$samples$norm.factors
    )
    for (profile_index in seq_len(nrow(aggregate_meta))) {
      audit_index <- audit_index + 1L
      sample_audit[[audit_index]] <- data.frame(
        partition = partition,
        aggregate_meta[profile_index, ],
        effective_library_size = sum(aggregate_counts[, profile_index]) *
          dge$samples$norm.factors[profile_index],
        stringsAsFactors = FALSE
      )
    }
  }
  list(states = state_results, audit = do.call(rbind, sample_audit))
}

message("Building independent donor residual matrices")
discovery_prepared <- aggregate_state_profiles(
  discovery_object$counts, discovery_metadata, network_genes, "discovery"
)
heldout_genes <- intersect(network_genes, rownames(heldout_counts))
heldout_prepared <- aggregate_state_profiles(
  heldout_counts, heldout_metadata, heldout_genes, "heldout"
)

expected_discovery_donors <- c(ABS = 27L, GOB = 26L, TAC = 25L)
observed_discovery_donors <- vapply(
  discovery_prepared$states,
  function(value) nrow(value$residual_donor_expression), integer(1)
)
if (!identical(observed_discovery_donors, expected_discovery_donors)) {
  stop("Discovery donor counts do not match the frozen contract")
}

good_gene_by_state <- lapply(
  discovery_prepared$states,
  function(value) {
    matrix_value <- value$residual_donor_expression
    apply(matrix_value, 2L, function(column) {
      all(is.finite(column)) && stats::sd(column) > 0
    })
  }
)
good_genes <- Reduce(`&`, good_gene_by_state)
excluded_zero_variance <- names(good_genes)[!good_genes]
network_genes_used <- names(good_genes)[good_genes]
if (length(network_genes_used) < 0.95 * length(network_genes)) {
  stop("More than 5% of the testable universe has no residual variance")
}

multi_expr <- lapply(
  discovery_prepared$states,
  function(value) list(data = as.data.frame(
    value$residual_donor_expression[, network_genes_used, drop = FALSE]
  ))
)

message("Selecting a common signed soft-thresholding power")
powers <- 1:30
soft_records <- list()
for (state in names(multi_expr)) {
  soft <- pickSoftThreshold(
    multi_expr[[state]]$data,
    powerVector = powers,
    networkType = "signed",
    corFnc = bicor,
    corOptions = list(use = "p", maxPOutliers = 0.05),
    moreNetworkConcepts = TRUE,
    verbose = 0
  )$fitIndices
  soft$state <- state
  soft_records[[state]] <- soft
}
soft_threshold <- do.call(rbind, soft_records)
power_gate <- do.call(rbind, lapply(powers, function(power) {
  local <- soft_threshold[soft_threshold$Power == power, , drop = FALSE]
  data.frame(
    power = power,
    n_r2_0.80 = sum(local$SFT.R.sq >= 0.80, na.rm = TRUE),
    minimum_r2 = min(local$SFT.R.sq, na.rm = TRUE),
    minimum_mean_connectivity = min(local$mean.k., na.rm = TRUE),
    passes = sum(local$SFT.R.sq >= 0.80, na.rm = TRUE) >= 2L &&
      min(local$SFT.R.sq, na.rm = TRUE) >= 0.70 &&
      min(local$mean.k., na.rm = TRUE) >= 5,
    stringsAsFactors = FALSE
  )
}))
passing_powers <- power_gate$power[power_gate$passes]
selected_power <- if (length(passing_powers)) min(passing_powers) else 16L
power_selection_method <- if (length(passing_powers)) {
  "first prespecified topology/connectivity gate"
} else {
  "prespecified signed-network fallback for 20-30 donors"
}

run_consensus <- function(
  consensus_quantile, tag, expression_sets = multi_expr
) {
  blockwiseConsensusModules(
    expression_sets,
    maxBlockSize = 9000,
    randomSeed = 20260830,
    individualTOMInfo = NULL,
    corType = "bicor",
    maxPOutliers = 0.05,
    power = selected_power,
    networkType = "signed",
    TOMType = "signed",
    saveIndividualTOMs = FALSE,
    individualTOMFileNames = file.path(
      tom_dir, paste0(tag, "-set%s-block%b.RData")
    ),
    saveConsensusTOMs = FALSE,
    networkCalibration = "single quantile",
    calibrationQuantile = 0.95,
    consensusQuantile = consensus_quantile,
    useMean = FALSE,
    deepSplit = 2,
    minModuleSize = 30,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    nThreads = 8,
    verbose = 3
  )
}

message("Constructing the primary state-consensus network")
primary_network <- run_consensus(0.25, "primary")
if (!length(primary_network$colors) || all(primary_network$colors == 0L)) {
  stop("The primary consensus network detected no non-grey module")
}
message("Running the two frozen consensus-quantile sensitivities")
sensitivity_q0 <- run_consensus(
  0, "q0"
)
sensitivity_q50 <- run_consensus(
  0.5, "q50"
)

primary_numeric <- primary_network$colors
primary_colors <- labels2colors(primary_numeric)
names(primary_numeric) <- names(primary_colors) <- network_genes_used
module_sizes_numeric <- sort(table(primary_numeric[primary_numeric != 0L]), decreasing = TRUE)
numeric_order <- as.integer(names(module_sizes_numeric))
module_label_map <- setNames(
  sprintf("M%02d", seq_along(numeric_order)), as.character(numeric_order)
)
module_labels <- ifelse(
  primary_numeric == 0L, "grey", module_label_map[as.character(primary_numeric)]
)
names(module_labels) <- network_genes_used
color_to_label <- unique(data.frame(
  module_color = primary_colors,
  module = module_labels,
  stringsAsFactors = FALSE
))

module_jaccard <- function(primary_labels, alternative_numeric, alternative_name) {
  alternative_labels <- paste0("A", alternative_numeric)
  records <- list()
  counter <- 0L
  for (module in setdiff(unique(primary_labels), "grey")) {
    genes <- names(primary_labels)[primary_labels == module]
    candidates <- setdiff(unique(alternative_labels), "A0")
    if (!length(candidates)) {
      counter <- counter + 1L
      records[[counter]] <- data.frame(
        module = module, sensitivity = alternative_name,
        best_alternative = NA_character_, maximum_jaccard = 0,
        stringsAsFactors = FALSE
      )
      next
    }
    values <- vapply(candidates, function(candidate) {
      alternative_genes <- names(primary_labels)[alternative_labels == candidate]
      length(intersect(genes, alternative_genes)) /
        length(union(genes, alternative_genes))
    }, numeric(1))
    best <- which.max(values)
    counter <- counter + 1L
    records[[counter]] <- data.frame(
      module = module, sensitivity = alternative_name,
      best_alternative = candidates[best], maximum_jaccard = values[best],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, records)
}
names(sensitivity_q0$colors) <- names(sensitivity_q50$colors) <- network_genes_used
parameter_stability <- rbind(
  module_jaccard(module_labels, sensitivity_q0$colors, "consensus_quantile_0"),
  module_jaccard(module_labels, sensitivity_q50$colors, "consensus_quantile_0.5")
)

sample_connectivity_records <- list()
sample_counter <- 0L
outliers_by_state <- list()
for (state in names(multi_expr)) {
  expression <- as.matrix(multi_expr[[state]]$data)
  sample_cor <- bicor(t(expression), use = "p", maxPOutliers = 0.05)
  sample_adjacency <- ((1 + sample_cor) / 2)^selected_power
  connectivity <- rowSums(sample_adjacency, na.rm = TRUE) - 1
  connectivity_z <- as.numeric(scale(connectivity))
  flagged <- connectivity_z < -3
  outliers_by_state[[state]] <- rownames(expression)[flagged]
  for (row_index in seq_len(nrow(expression))) {
    sample_counter <- sample_counter + 1L
    sample_connectivity_records[[sample_counter]] <- data.frame(
      state = state, donor_id = rownames(expression)[row_index],
      connectivity = connectivity[row_index],
      connectivity_z = connectivity_z[row_index],
      audit_outlier = flagged[row_index],
      stringsAsFactors = FALSE
    )
  }
}
sample_connectivity <- do.call(rbind, sample_connectivity_records)

calculate_kme <- function(data_expression, colors) {
  eigengenes <- orderMEs(moduleEigengenes(
    data_expression, colors = colors, excludeGrey = TRUE
  )$eigengenes)
  kme <- signedKME(
    data_expression, eigengenes,
    corFnc = "bicor", corOptions = "use='p', maxPOutliers=0.05"
  )
  list(eigengenes = eigengenes, kme = kme)
}

message("Calculating state-resolved module membership")
discovery_kme <- list()
heldout_kme <- list()
for (state in c("ABS", "GOB", "TAC")) {
  discovery_matrix <- as.matrix(multi_expr[[state]]$data)
  discovery_kme[[state]] <- calculate_kme(discovery_matrix, primary_colors)
  shared <- intersect(network_genes_used, colnames(
    heldout_prepared$states[[state]]$residual_donor_expression
  ))
  heldout_matrix <- heldout_prepared$states[[state]]$residual_donor_expression[, shared, drop = FALSE]
  heldout_kme[[state]] <- calculate_kme(
    heldout_matrix, primary_colors[shared]
  )
}

extract_gene_kme <- function(kme_list, state, partition) {
  frame <- kme_list[[state]]$kme
  records <- vector("list", nrow(frame))
  for (row_index in seq_len(nrow(frame))) {
    gene <- rownames(frame)[row_index]
    module <- module_labels[gene]
    column <- paste0("kME", primary_colors[gene])
    value <- if (module != "grey" && column %in% colnames(frame)) {
      frame[row_index, column]
    } else {
      NA_real_
    }
    records[[row_index]] <- data.frame(
      gene = gene, partition = partition, state = state,
      module = unname(module), kme = value, absolute_kme = abs(value),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, records)
}
kme_long <- do.call(rbind, c(
  lapply(c("ABS", "GOB", "TAC"), function(state) {
    extract_gene_kme(discovery_kme, state, "discovery")
  }),
  lapply(c("ABS", "GOB", "TAC"), function(state) {
    extract_gene_kme(heldout_kme, state, "heldout")
  })
))

message("Running held-out module preservation with 1,000 permutations per state")
run_preservation_state <- function(state) {
  # Network construction is complete before this function.  Use one WGCNA
  # thread per permutation and six deterministic foreach workers per state;
  # the three states are launched by the outer mclapply below.
  disableWGCNAThreads()
  registerDoParallel(cores = 6L)
  on.exit(stopImplicitCluster(), add = TRUE)
  shared <- intersect(
    network_genes_used,
    colnames(heldout_prepared$states[[state]]$residual_donor_expression)
  )
  reference <- as.data.frame(
    discovery_prepared$states[[state]]$residual_donor_expression[, shared, drop = FALSE]
  )
  test <- as.data.frame(
    heldout_prepared$states[[state]]$residual_donor_expression[, shared, drop = FALSE]
  )
  finite_variance <- apply(reference, 2L, sd) > 0 & apply(test, 2L, sd) > 0
  reference <- reference[, finite_variance, drop = FALSE]
  test <- test[, finite_variance, drop = FALSE]
  colors <- primary_colors[colnames(reference)]
  set.seed(20260830 + match(state, c("ABS", "GOB", "TAC")))
  preservation <- modulePreservation(
    multiData = list(reference = list(data = reference), test = list(data = test)),
    multiColor = list(reference = colors, test = rep("grey", ncol(test))),
    dataIsExpr = TRUE,
    networkType = "signed",
    corFnc = "bicor",
    corOptions = "use='p', maxPOutliers=0.05",
    referenceNetworks = 1,
    testNetworks = 2,
    nPermutations = 1000,
    calculateQvalue = FALSE,
    maxModuleSize = 9000,
    quickCor = 0,
    savePermutedStatistics = FALSE,
    parallelCalculation = TRUE,
    randomSeed = 20260830 + match(state, c("ABS", "GOB", "TAC")),
    verbose = 2
  )
  z_table <- preservation$preservation$Z$ref.reference$inColumnsAlsoPresentIn.test
  observed_table <- preservation$preservation$observed$ref.reference$inColumnsAlsoPresentIn.test
  module_colors <- intersect(rownames(z_table), color_to_label$module_color)
  preservation_records <- list()
  preservation_counter <- 0L
  for (module_color in module_colors) {
    preservation_counter <- preservation_counter + 1L
    preservation_records[[preservation_counter]] <- data.frame(
      state = state,
      module_color = module_color,
      module = color_to_label$module[match(module_color, color_to_label$module_color)],
      module_size_tested = z_table[module_color, "moduleSize"],
      zsummary = z_table[module_color, "Zsummary.pres"],
      zdensity = z_table[module_color, "Zdensity.pres"],
      zconnectivity = z_table[module_color, "Zconnectivity.pres"],
      median_rank = observed_table[module_color, "medianRank.pres"],
      stringsAsFactors = FALSE
    )
  }
  list(
    state = state,
    table = do.call(rbind, preservation_records),
    object = preservation
  )
}
preservation_states <- c("ABS", "GOB", "TAC")
preservation_results <- parallel::mclapply(
  preservation_states,
  run_preservation_state,
  mc.cores = 3L,
  mc.preschedule = FALSE,
  mc.set.seed = FALSE
)
if (any(vapply(preservation_results, inherits, logical(1), what = "try-error"))) {
  stop("At least one state-level preservation worker failed")
}
preservation_table <- do.call(
  rbind, lapply(preservation_results, function(value) value$table)
)
preservation_objects <- setNames(
  lapply(preservation_results, function(value) value$object),
  preservation_states
)

ranking_contexts <- list(
  discovery_common = setNames(ranking$common_z, ranking$gene),
  discovery_ABS = setNames(ranking$logFC_ABS / ranking$raw_se_ABS, ranking$gene),
  discovery_GOB = setNames(ranking$logFC_GOB / ranking$raw_se_GOB, ranking$gene),
  discovery_TAC = setNames(ranking$logFC_TAC / ranking$raw_se_TAC, ranking$gene),
  heldout_common = setNames(heldout_ranking$common_z, heldout_ranking$gene),
  heldout_ABS = setNames(heldout_ranking$logFC_ABS / heldout_ranking$raw_se_ABS, heldout_ranking$gene),
  heldout_GOB = setNames(heldout_ranking$logFC_GOB / heldout_ranking$raw_se_GOB, heldout_ranking$gene),
  heldout_TAC = setNames(heldout_ranking$logFC_TAC / heldout_ranking$raw_se_TAC, heldout_ranking$gene)
)
module_gene_sets <- split(
  names(module_labels)[module_labels != "grey"],
  module_labels[module_labels != "grey"]
)
module_enrichment_records <- list()
module_enrichment_counter <- 0L
for (context in names(ranking_contexts)) {
  statistic <- ranking_contexts[[context]]
  statistic <- statistic[is.finite(statistic)]
  local_sets <- lapply(module_gene_sets, intersect, y = names(statistic))
  local_sets <- local_sets[lengths(local_sets) >= 10L]
  result <- cameraPR(
    statistic,
    lapply(local_sets, match, table = names(statistic)),
    use.ranks = TRUE, inter.gene.cor = 0.01,
    sort = FALSE, directional = TRUE
  )
  result$module <- rownames(result)
  rownames(result) <- NULL
  for (row_index in seq_len(nrow(result))) {
    module_enrichment_counter <- module_enrichment_counter + 1L
    module_enrichment_records[[module_enrichment_counter]] <- data.frame(
      module = result$module[row_index], context = context,
      n_genes = result$NGenes[row_index],
      direction = result$Direction[row_index],
      direction_sign = ifelse(result$Direction[row_index] == "Up", 1L, -1L),
      p_value = result$PValue[row_index], fdr = result$FDR[row_index],
      stringsAsFactors = FALSE
    )
  }
}
module_enrichment <- do.call(rbind, module_enrichment_records)

overlap_test <- function(query_genes, target_genes, universe_genes) {
  query <- intersect(query_genes, universe_genes)
  target <- intersect(target_genes, universe_genes)
  overlap <- length(intersect(query, target))
  matrix_value <- matrix(c(
    overlap,
    length(query) - overlap,
    length(target) - overlap,
    length(universe_genes) - length(union(query, target))
  ), nrow = 2L)
  data.frame(
    overlap = overlap,
    query_size = length(query), target_size = length(target),
    odds_ratio = unname(fisher.test(matrix_value, alternative = "greater")$estimate),
    p_value = fisher.test(matrix_value, alternative = "greater")$p.value,
    stringsAsFactors = FALSE
  )
}

strict_overlap_records <- list()
strict_counter <- 0L
strict_sets <- list(
  strict_any = strict_genes,
  strict_up = ranking$gene[ranking$strict_state_shared & ranking$shared_direction == "up"],
  strict_down = ranking$gene[ranking$strict_state_shared & ranking$shared_direction == "down"]
)
for (module in names(module_gene_sets)) {
  for (target_name in names(strict_sets)) {
    strict_counter <- strict_counter + 1L
    value <- overlap_test(module_gene_sets[[module]], strict_sets[[target_name]], universe)
    strict_overlap_records[[strict_counter]] <- data.frame(
      module = module, target = target_name, value,
      stringsAsFactors = FALSE
    )
  }
}
module_overlap <- do.call(rbind, strict_overlap_records)
module_overlap$fdr <- ave(
  module_overlap$p_value, module_overlap$target,
  FUN = function(value) p.adjust(value, method = "BH")
)

pathway_communities <- read.delim(paths$pathway_communities, check.names = FALSE)
pathway_leading <- read.delim(paths$pathway_leading_edges, check.names = FALSE)
community_members <- merge(
  pathway_communities[, c("collection", "gene_set", "community_id", "community_label")],
  pathway_leading[, c("collection", "gene_set", "leading_edge")],
  by = c("collection", "gene_set"), all.x = TRUE
)
community_gene_sets <- lapply(
  split(community_members, community_members$community_id),
  function(frame) {
    leading_edge <- frame$leading_edge
    leading_edge[is.na(leading_edge)] <- ""
    unique(unlist(strsplit(leading_edge, ";", fixed = TRUE)))
  }
)
community_gene_sets <- lapply(
  community_gene_sets,
  function(value) value[!is.na(value) & nzchar(value)]
)
community_overlap_records <- list()
community_counter <- 0L
for (module in names(module_gene_sets)) {
  for (community in names(community_gene_sets)) {
    community_counter <- community_counter + 1L
    value <- overlap_test(
      module_gene_sets[[module]], community_gene_sets[[community]], universe
    )
    community_overlap_records[[community_counter]] <- data.frame(
      module = module, community_id = as.integer(community), value,
      stringsAsFactors = FALSE
    )
  }
}
module_community_overlap <- do.call(rbind, community_overlap_records)
module_community_overlap$fdr <- ave(
  module_community_overlap$p_value, module_community_overlap$community_id,
  FUN = function(value) p.adjust(value, method = "BH")
)
community_labels <- unique(pathway_communities[, c("community_id", "community_label")])
module_community_overlap <- merge(
  module_community_overlap, community_labels,
  by = "community_id", all.x = TRUE
)

module_summary <- do.call(rbind, lapply(names(module_gene_sets), function(module) {
  preservation_local <- preservation_table[preservation_table$module == module, ]
  enrichment_local <- module_enrichment[module_enrichment$module == module, ]
  heldout_common <- enrichment_local[enrichment_local$context == "heldout_common", ]
  heldout_states <- enrichment_local[enrichment_local$context %in%
    c("heldout_ABS", "heldout_GOB", "heldout_TAC"), ]
  stability_local <- parameter_stability[parameter_stability$module == module, ]
  overlap_local <- module_overlap[module_overlap$module == module, ]
  community_local <- module_community_overlap[module_community_overlap$module == module, ]
  preservation_complete <- nrow(preservation_local) == 3L &&
    all(is.finite(preservation_local$zsummary))
  preservation_pass <- preservation_complete &&
    sum(preservation_local$zsummary >= 2) >= 2L &&
    all(preservation_local$zsummary >= 0)
  enrichment_pass <- nrow(heldout_common) == 1L &&
    is.finite(heldout_common$fdr) && heldout_common$fdr <= 0.10 &&
    nrow(heldout_states) == 3L &&
    all(is.finite(heldout_states$direction_sign)) &&
    all(heldout_states$direction_sign == heldout_common$direction_sign)
  overlap_pass <- any(overlap_local$fdr <= 0.05, na.rm = TRUE) ||
    any(community_local$fdr <= 0.05, na.rm = TRUE)
  parameter_stable <- nrow(stability_local) == 2L &&
    all(is.finite(stability_local$maximum_jaccard)) &&
    all(stability_local$maximum_jaccard >= 0.50)
  data.frame(
    module = module,
    module_color = color_to_label$module_color[
      match(module, color_to_label$module)
    ],
    module_size = length(module_gene_sets[[module]]),
    strict_programme_overlap = sum(module_gene_sets[[module]] %in% strict_genes),
    preservation_states_z2 = sum(preservation_local$zsummary >= 2, na.rm = TRUE),
    minimum_preservation_z = if (preservation_complete) {
      min(preservation_local$zsummary)
    } else {
      NA_real_
    },
    heldout_direction = if (nrow(heldout_common) == 1L) heldout_common$direction else NA_character_,
    heldout_common_fdr = if (nrow(heldout_common) == 1L) heldout_common$fdr else NA_real_,
    preservation_pass = preservation_pass,
    heldout_enrichment_pass = enrichment_pass,
    biological_overlap_pass = overlap_pass,
    parameter_stable = parameter_stable,
    internal_gate_pass = preservation_pass && enrichment_pass &&
      overlap_pass && parameter_stable,
    external_gate_status = "pending",
    routing_status = "pending_external_validation",
    stringsAsFactors = FALSE
  )
}))

consensus_modules <- data.frame(
  gene = network_genes_used,
  module = unname(module_labels[network_genes_used]),
  module_color = unname(primary_colors[network_genes_used]),
  module_numeric = unname(primary_numeric[network_genes_used]),
  strict_state_shared = ranking$strict_state_shared[match(network_genes_used, ranking$gene)],
  programme_direction = ranking$shared_direction[match(network_genes_used, ranking$gene)],
  discovery_common_z = ranking$common_z[match(network_genes_used, ranking$gene)],
  heldout_common_z = heldout_ranking$common_z[match(network_genes_used, heldout_ranking$gene)],
  heldout_common_fdr = heldout_ranking$common_q_value[
    match(network_genes_used, heldout_ranking$gene)
  ],
  stringsAsFactors = FALSE
)
for (state in c("ABS", "GOB", "TAC")) {
  discovery_values <- kme_long[
    kme_long$partition == "discovery" & kme_long$state == state,
    c("gene", "kme", "absolute_kme")
  ]
  heldout_values <- kme_long[
    kme_long$partition == "heldout" & kme_long$state == state,
    c("gene", "kme", "absolute_kme")
  ]
  consensus_modules[[paste0("kme_discovery_", state)]] <- discovery_values$kme[
    match(consensus_modules$gene, discovery_values$gene)
  ]
  consensus_modules[[paste0("abs_kme_discovery_", state)]] <- discovery_values$absolute_kme[
    match(consensus_modules$gene, discovery_values$gene)
  ]
  consensus_modules[[paste0("kme_heldout_", state)]] <- heldout_values$kme[
    match(consensus_modules$gene, heldout_values$gene)
  ]
  consensus_modules[[paste0("abs_kme_heldout_", state)]] <- heldout_values$absolute_kme[
    match(consensus_modules$gene, heldout_values$gene)
  ]
}
consensus_modules$module_size <- module_summary$module_size[
  match(consensus_modules$module, module_summary$module)
]

write.table(
  consensus_modules, file.path(out_root, "consensus_modules.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  module_summary, file.path(out_dir, "module_internal_gate_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  preservation_table, file.path(out_dir, "module_preservation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  module_enrichment, file.path(out_dir, "module_rank_enrichment.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  module_overlap, file.path(out_dir, "module_strict_programme_overlap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  module_community_overlap,
  file.path(out_dir, "module_pathway_community_overlap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  parameter_stability, file.path(audit_dir, "module_parameter_stability.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  sample_connectivity, file.path(audit_dir, "sample_connectivity_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  rbind(discovery_prepared$audit, heldout_prepared$audit),
  file.path(audit_dir, "network_profile_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  soft_threshold, file.path(audit_dir, "soft_threshold_by_state.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  power_gate, file.path(audit_dir, "soft_threshold_selection_gate.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(gene = excluded_zero_variance),
  file.path(audit_dir, "zero_residual_variance_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(
  list(
    network_genes = network_genes_used,
    selected_power = selected_power,
    module_labels = module_labels,
    module_colors = primary_colors,
    discovery = discovery_prepared,
    heldout = heldout_prepared,
    discovery_kme = discovery_kme,
    heldout_kme = heldout_kme
  ),
  file.path(out_dir, "network_projection_object.rds"),
  compress = "xz"
)
saveRDS(
  list(
    primary = primary_network,
    sensitivity_q0_colors = sensitivity_q0$colors,
    sensitivity_q50_colors = sensitivity_q50$colors,
    preservation = preservation_objects
  ),
  file.path(audit_dir, "consensus_network_audit_object.rds"),
  compress = FALSE
)

output_paths <- c(
  file.path(out_root, "consensus_modules.tsv"),
  file.path(out_dir, "module_internal_gate_summary.tsv"),
  file.path(out_dir, "module_preservation.tsv"),
  file.path(out_dir, "module_rank_enrichment.tsv"),
  file.path(out_dir, "module_strict_programme_overlap.tsv"),
  file.path(out_dir, "module_pathway_community_overlap.tsv"),
  file.path(audit_dir, "module_parameter_stability.tsv"),
  file.path(audit_dir, "sample_connectivity_audit.tsv")
)
manifest <- list(
  analysis = "state_aware_consensus_wgcna_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = 20260830,
  input_sha256 = as.list(c(
    observed_hashes,
    pathway_communities = sha256(paths$pathway_communities),
    pathway_leading_edges = sha256(paths$pathway_leading_edges)
  )),
  donor_overlap_removed = overlap_donors,
  discovery_donors_by_state = as.list(observed_discovery_donors),
  heldout_donors_by_state = as.list(vapply(
    heldout_prepared$states,
    function(value) nrow(value$residual_donor_expression), integer(1)
  )),
  input_gene_universe = length(universe),
  network_genes_used = length(network_genes_used),
  zero_variance_genes = length(excluded_zero_variance),
  selected_power = selected_power,
  power_selection_method = power_selection_method,
  parameters = list(
    network_type = "signed", tom_type = "signed", correlation = "bicor",
    max_p_outliers = 0.05, min_module_size = 30, deep_split = 2,
    merge_cut_height = 0.25, consensus_quantile = 0.25,
    sensitivity_quantiles = c(0, 0.5), preservation_permutations = 1000,
    preservation_state_workers = 3, preservation_workers_per_state = 6,
    wgcna_threads_per_permutation = 1
  ),
  module_count = nrow(module_summary),
  internally_eligible_modules = module_summary$module[module_summary$internal_gate_pass],
  audit_outlier_count = sum(sample_connectivity$audit_outlier),
  package_versions = as.list(vapply(
    c(
      "WGCNA", "doParallel", "foreach", "edgeR", "limma", "Matrix",
      "igraph", "digest", "jsonlite"
    ),
    function(package) as.character(packageVersion(package)), character(1)
  )),
  output_sha256 = as.list(setNames(vapply(output_paths, sha256, character(1)), basename(output_paths))),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest, file.path(out_dir, "consensus_wgcna_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)

message(
  "Consensus WGCNA complete: ", nrow(module_summary), " modules; ",
  sum(module_summary$internal_gate_pass),
  " pass preservation, held-out enrichment, overlap, and parameter-stability gates"
)
