#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(decoupleR)
  library(igraph)
  library(limma)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)
set.seed(20260829)

root <- normalizePath(".", mustWork = TRUE)
common_path <- file.path(
  root,
  "results",
  "state_aware_program_v1",
  "common_effects",
  "cross_state_common_effects.tsv.gz"
)
hallmark_path <- file.path(
  root,
  "data_sources",
  "reference_gene_sets_v2_5",
  "h.all.v2026.1.Hs.symbols.gmt"
)
regulator_path <- file.path(
  root,
  "data_sources",
  "regulatory_priors",
  "collectri_tcf_ascl2_wnt_controls_2026-08-08.tsv"
)
contract_path <- file.path(
  root,
  "analysis",
  "contracts",
  "state_aware_program_rederivation_v1_2026-08-29.md"
)
out_dir <- file.path(root, "results", "state_aware_program_v1", "interpretation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_contract_sha256 <-
  "0f39a03154408d14bb0bbe0c1e4f55d498e94d403bd6b45d4be5dd8918555f69"
if (!identical(
  digest::digest(contract_path, algo = "sha256", file = TRUE),
  expected_contract_sha256
)) {
  stop("The frozen parent contract changed before biological interpretation")
}
if (!all(file.exists(c(common_path, hallmark_path, regulator_path)))) {
  stop("One or more interpretation inputs are missing")
}

common <- read.delim(common_path, check.names = FALSE)
statistics <- common$common_z
names(statistics) <- common$gene
if (anyDuplicated(names(statistics)) || any(!is.finite(statistics))) {
  stop("Invalid common-effect ranking")
}
strict_genes <- common$gene[common$strict_state_shared]

read_gmt <- function(path) {
  lines <- readLines(path, encoding = "UTF-8")
  fields <- strsplit(lines, "\t", fixed = TRUE)
  gene_sets <- lapply(fields, function(value) unique(value[-c(1L, 2L)]))
  names(gene_sets) <- vapply(fields, `[[`, character(1), 1L)
  descriptions <- vapply(fields, `[[`, character(1), 2L)
  names(descriptions) <- names(gene_sets)
  list(gene_sets = gene_sets, descriptions = descriptions)
}

prepare_msigdb <- function(collection, subcollection) {
  frame <- msigdbr(
    species = "Homo sapiens",
    collection = collection,
    subcollection = subcollection
  )
  list(
    gene_sets = split(frame$gene_symbol, frame$gs_name),
    descriptions = setNames(frame$gs_description, frame$gs_name)
  )
}

run_camera_collection <- function(resource, collection_label) {
  gene_sets <- lapply(
    resource$gene_sets,
    function(genes) intersect(unique(genes), names(statistics))
  )
  sizes <- lengths(gene_sets)
  gene_sets <- gene_sets[sizes >= 15L & sizes <= 500L]
  if (!length(gene_sets)) {
    stop("No eligible gene sets for ", collection_label)
  }
  index <- lapply(gene_sets, match, table = names(statistics))
  result <- cameraPR(
    statistic = statistics,
    index = index,
    use.ranks = TRUE,
    inter.gene.cor = 0.01,
    sort = FALSE,
    directional = TRUE
  )
  result$gene_set <- rownames(result)
  rownames(result) <- NULL
  result$collection <- collection_label
  result$description <- unname(resource$descriptions[result$gene_set])
  result$strict_core_overlap <- vapply(
    result$gene_set,
    function(gene_set) sum(gene_sets[[gene_set]] %in% strict_genes),
    integer(1)
  )
  result$top_contributing_genes <- vapply(
    seq_len(nrow(result)),
    function(index_row) {
      gene_set <- result$gene_set[index_row]
      genes <- gene_sets[[gene_set]]
      gene_statistics <- statistics[genes]
      if (result$Direction[index_row] == "Up") {
        ordered <- names(sort(gene_statistics, decreasing = TRUE))
      } else {
        ordered <- names(sort(gene_statistics, decreasing = FALSE))
      }
      paste(head(ordered, 12L), collapse = ";")
    },
    character(1)
  )
  result <- result[, c(
    "collection",
    "gene_set",
    "description",
    "NGenes",
    "Direction",
    "PValue",
    "FDR",
    "strict_core_overlap",
    "top_contributing_genes"
  )]
  attr(result, "gene_sets") <- gene_sets
  result
}

message("Running Hallmark competitive enrichment")
hallmark_resource <- read_gmt(hallmark_path)
hallmark_result <- run_camera_collection(hallmark_resource, "Hallmark 2026.1.Hs")

message("Running Reactome competitive enrichment")
reactome_resource <- prepare_msigdb("C2", "CP:REACTOME")
reactome_result <- run_camera_collection(reactome_resource, "Reactome 2026.1.Hs")

message("Running GO Biological Process competitive enrichment")
go_resource <- prepare_msigdb("C5", "GO:BP")
go_result <- run_camera_collection(go_resource, "GO:BP 2026.1.Hs")

pathway_results <- rbind(hallmark_result, reactome_result, go_result)
pathway_results$global_fdr <- p.adjust(pathway_results$PValue, method = "BH")
pathway_results <- pathway_results[
  order(pathway_results$collection, pathway_results$FDR, pathway_results$PValue),
  ,
  drop = FALSE
]

pathway_family <- function(name) {
  name <- toupper(name)
  if (grepl("WNT|BETA_CATENIN|STEM|CRYPT", name)) {
    return("WNT/stem-regulatory")
  }
  if (grepl("DIFFERENTIATION|ABSORP|DIGEST|TRANSPORT|BRUSH_BORDER|ENTEROCYTE", name)) {
    return("mature epithelial function")
  }
  if (grepl("OXIDATIVE_PHOSPHORYLATION|MITOCHON|FATTY_ACID|LIPID|BILE_ACID|METABOL", name)) {
    return("metabolic function")
  }
  if (grepl("E2F|G2M|CELL_CYCLE|MITOTIC|MYC_TARGET", name)) {
    return("cell cycle/proliferation")
  }
  if (grepl("INFLAM|INTERFERON|TNFA|NF_KB|STRESS|HYPOXIA|REACTIVE_OXYGEN", name)) {
    return("stress/inflammatory")
  }
  "other"
}
pathway_results$interpretive_family <- vapply(
  pathway_results$gene_set,
  pathway_family,
  character(1)
)

# Build a redundancy network only for presentation. Inferential P values and
# FDR values above always remain at the original gene-set level.
display_candidates <- do.call(rbind, lapply(
  split(pathway_results, pathway_results$collection),
  function(frame) {
    significant <- frame[frame$FDR <= 0.05, , drop = FALSE]
    if (!nrow(significant)) {
      return(significant)
    }
    significant <- significant[order(significant$FDR, significant$PValue), , drop = FALSE]
    head(significant, 200L)
  }
))

resource_by_collection <- list(
  `Hallmark 2026.1.Hs` = hallmark_resource$gene_sets,
  `Reactome 2026.1.Hs` = reactome_resource$gene_sets,
  `GO:BP 2026.1.Hs` = go_resource$gene_sets
)
network_edges <- list()
if (nrow(display_candidates) > 1L) {
  counter <- 1L
  for (collection_label in unique(display_candidates$collection)) {
    subset <- display_candidates[
      display_candidates$collection == collection_label,
      ,
      drop = FALSE
    ]
    if (nrow(subset) < 2L) {
      next
    }
    sets <- resource_by_collection[[collection_label]]
    combinations <- utils::combn(subset$gene_set, 2L, simplify = FALSE)
    for (pair in combinations) {
      first <- intersect(sets[[pair[1L]]], names(statistics))
      second <- intersect(sets[[pair[2L]]], names(statistics))
      jaccard <- length(intersect(first, second)) / length(union(first, second))
      if (is.finite(jaccard) && jaccard >= 0.25) {
        network_edges[[counter]] <- data.frame(
          collection = collection_label,
          from = pair[1L],
          to = pair[2L],
          jaccard = jaccard,
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }
}
if (length(network_edges)) {
  network_edges <- do.call(rbind, network_edges)
} else {
  network_edges <- data.frame(
    collection = character(),
    from = character(),
    to = character(),
    jaccard = numeric(),
    stringsAsFactors = FALSE
  )
}

message("Inferring signed regulator activities from the locked CollecTRI prior")
regulator_prior <- read.delim(regulator_path, check.names = FALSE)
as_logical <- function(value) {
  tolower(trimws(as.character(value))) %in% c("true", "t", "1", "yes")
}
stimulation <- as_logical(regulator_prior$consensus_stimulation)
inhibition <- as_logical(regulator_prior$consensus_inhibition)
regulator_prior$mor <- ifelse(
  stimulation & !inhibition,
  1,
  ifelse(inhibition & !stimulation, -1, NA_real_)
)
network <- regulator_prior[
  is.finite(regulator_prior$mor),
  c("source_genesymbol", "target_genesymbol", "mor"),
  drop = FALSE
]
colnames(network) <- c("source", "target", "mor")
network <- aggregate(mor ~ source + target, data = network, FUN = mean)
network <- network[network$target %in% names(statistics) & network$mor != 0, ]

activity_matrix <- matrix(
  statistics,
  ncol = 1L,
  dimnames = list(names(statistics), "common_effect_z")
)
regulator_result <- run_ulm(
  mat = activity_matrix,
  network = network,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5L
)
regulator_result$q_value <- p.adjust(regulator_result$p_value, method = "BH")
regulator_result$n_targets_in_universe <- vapply(
  regulator_result$source,
  function(source) sum(network$source == source),
  integer(1)
)
regulator_result <- regulator_result[
  order(regulator_result$q_value, -abs(regulator_result$score)),
  ,
  drop = FALSE
]

write.table(
  pathway_results,
  gzfile(file.path(out_dir, "pathway_competitive_enrichment_all.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  display_candidates,
  file.path(out_dir, "pathway_display_candidates.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  network_edges,
  file.path(out_dir, "pathway_redundancy_network_edges.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  regulator_result,
  file.path(out_dir, "collectri_regulator_activity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  network,
  gzfile(file.path(out_dir, "collectri_regulatory_edges_tested.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

collection_summary <- do.call(rbind, lapply(
  split(pathway_results, pathway_results$collection),
  function(frame) {
    data.frame(
      collection = unique(frame$collection),
      tested_gene_sets = nrow(frame),
      fdr_le_0.05 = sum(frame$FDR <= 0.05),
      global_fdr_le_0.05 = sum(frame$global_fdr <= 0.05),
      up_fdr_le_0.05 = sum(frame$FDR <= 0.05 & frame$Direction == "Up"),
      down_fdr_le_0.05 = sum(frame$FDR <= 0.05 & frame$Direction == "Down"),
      stringsAsFactors = FALSE
    )
  }
))
write.table(
  collection_summary,
  file.path(out_dir, "interpretation_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- list(
  analysis = "state_aware_interpret_common_program_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  contract_sha256 = expected_contract_sha256,
  input_sha256 = digest::digest(common_path, algo = "sha256", file = TRUE),
  hallmark_sha256 = digest::digest(hallmark_path, algo = "sha256", file = TRUE),
  regulator_prior_sha256 = digest::digest(
    regulator_path,
    algo = "sha256",
    file = TRUE
  ),
  pathway_test = "limma::cameraPR on the complete signed common-effect ranking",
  pathway_size_range = c(15L, 500L),
  within_collection_fdr = "Benjamini-Hochberg as returned by cameraPR",
  global_fdr = "Benjamini-Hochberg across all three collections",
  regulator_method = "decoupleR::run_ulm on signed locked CollecTRI interactions",
  redundancy_network = "display-only; top 200 significant sets per collection; Jaccard >=0.25",
  random_seed = 20260829,
  db_version = unique(c(
    reactome_result$collection,
    go_result$collection,
    hallmark_result$collection
  )),
  summary = collection_summary,
  package_versions = as.list(vapply(
    c("decoupleR", "igraph", "limma", "msigdbr", "digest", "jsonlite"),
    function(package) as.character(packageVersion(package)),
    character(1)
  )),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "interpretation_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Biological interpretation complete: ",
  sum(pathway_results$FDR <= 0.05),
  " gene sets and ",
  sum(regulator_result$q_value <= 0.05),
  " regulators pass within-family FDR 0.05"
)
