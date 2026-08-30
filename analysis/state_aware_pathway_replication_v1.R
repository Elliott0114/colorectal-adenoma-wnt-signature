#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fgsea)
  library(igraph)
  library(limma)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)
set.seed(20260830)

root <- normalizePath(".", mustWork = TRUE)
parent_root <- file.path(root, "results", "state_aware_program_v1")
out_root <- file.path(parent_root, "functional_architecture_v1")
out_dir <- file.path(out_root, "pathway_replication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  discovery = file.path(
    parent_root, "common_effects", "cross_state_common_effects.tsv.gz"
  ),
  heldout = file.path(
    parent_root, "heldout_validation", "heldout_cross_state_common_effects.tsv.gz"
  )
)
expected_hashes <- c(
  discovery = "a1ac4b7b67ac279782e04e386971d7463e169cbdbc24d7da4c314e34a4f3e946",
  heldout = "ab9918f1481fcb52923133cafc64c59d16eabe3995bf46f7714da742acf95fc5"
)
if (!all(file.exists(unlist(paths)))) {
  stop("At least one pathway-replication input is missing")
}
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
observed_hashes <- vapply(paths[c("discovery", "heldout")], sha256, character(1))
if (!identical(unname(observed_hashes), unname(expected_hashes[names(observed_hashes)]))) {
  stop("A frozen functional-architecture input changed")
}

discovery <- read.delim(paths$discovery, check.names = FALSE)
heldout <- read.delim(paths$heldout, check.names = FALSE)
if (nrow(discovery) != 8221L || sum(discovery$strict_state_shared) != 1843L) {
  stop("The discovery universe or strict programme is not frozen as expected")
}
if (sum(discovery$strict_state_shared & discovery$shared_direction == "up") != 884L ||
    sum(discovery$strict_state_shared & discovery$shared_direction == "down") != 959L) {
  stop("The directional arms changed")
}

finite_named <- function(value, genes, label) {
  keep <- is.finite(value) & !is.na(genes) & nzchar(genes)
  value <- value[keep]
  names(value) <- genes[keep]
  if (anyDuplicated(names(value))) {
    stop("Duplicated genes in ", label)
  }
  value
}

contexts <- list(
  discovery_common = finite_named(discovery$common_z, discovery$gene, "discovery_common"),
  discovery_ABS = finite_named(discovery$logFC_ABS / discovery$raw_se_ABS, discovery$gene, "discovery_ABS"),
  discovery_GOB = finite_named(discovery$logFC_GOB / discovery$raw_se_GOB, discovery$gene, "discovery_GOB"),
  discovery_TAC = finite_named(discovery$logFC_TAC / discovery$raw_se_TAC, discovery$gene, "discovery_TAC"),
  heldout_common = finite_named(heldout$common_z, heldout$gene, "heldout_common"),
  heldout_ABS = finite_named(heldout$logFC_ABS / heldout$raw_se_ABS, heldout$gene, "heldout_ABS"),
  heldout_GOB = finite_named(heldout$logFC_GOB / heldout$raw_se_GOB, heldout$gene, "heldout_GOB"),
  heldout_TAC = finite_named(heldout$logFC_TAC / heldout$raw_se_TAC, heldout$gene, "heldout_TAC")
)

read_gmt <- function(path) {
  fields <- strsplit(readLines(path, encoding = "UTF-8"), "\t", fixed = TRUE)
  sets <- lapply(fields, function(value) unique(value[-c(1L, 2L)]))
  names(sets) <- vapply(fields, `[[`, character(1), 1L)
  descriptions <- vapply(fields, `[[`, character(1), 2L)
  names(descriptions) <- names(sets)
  list(gene_sets = sets, descriptions = descriptions)
}

prepare_msigdb <- function(collection, subcollection = NULL) {
  frame <- if (is.null(subcollection)) {
    msigdbr(species = "Homo sapiens", collection = collection)
  } else {
    msigdbr(
      species = "Homo sapiens",
      collection = collection,
      subcollection = subcollection
    )
  }
  list(
    gene_sets = lapply(split(frame$gene_symbol, frame$gs_name), unique),
    descriptions = setNames(frame$gs_description, frame$gs_name)
  )
}

resources <- list(
  `Hallmark 2026.1.Hs` = prepare_msigdb("H"),
  `Reactome 2026.1.Hs` = prepare_msigdb("C2", "CP:REACTOME"),
  `GO:BP 2026.1.Hs` = prepare_msigdb("C5", "GO:BP")
)

run_camera <- function(resource, collection_label, statistic, context_label) {
  sets <- lapply(
    resource$gene_sets,
    function(genes) intersect(unique(genes), names(statistic))
  )
  sizes <- lengths(sets)
  sets <- sets[sizes >= 10L & sizes <= 500L]
  index <- lapply(sets, match, table = names(statistic))
  result <- cameraPR(
    statistic = statistic,
    index = index,
    use.ranks = TRUE,
    inter.gene.cor = 0.01,
    sort = FALSE,
    directional = TRUE
  )
  result$gene_set <- rownames(result)
  rownames(result) <- NULL
  data.frame(
    collection = collection_label,
    gene_set = result$gene_set,
    description = unname(resource$descriptions[result$gene_set]),
    context = context_label,
    n_genes = result$NGenes,
    direction = result$Direction,
    direction_sign = ifelse(result$Direction == "Up", 1L, -1L),
    p_value = result$PValue,
    fdr = result$FDR,
    stringsAsFactors = FALSE
  )
}

message("Running competitive enrichment in eight frozen contexts")
records <- list()
counter <- 0L
for (collection_label in names(resources)) {
  for (context_label in names(contexts)) {
    counter <- counter + 1L
    records[[counter]] <- run_camera(
      resources[[collection_label]], collection_label,
      contexts[[context_label]], context_label
    )
  }
}
results_long <- do.call(rbind, records)

context_names <- names(contexts)
state_contexts <- c(
  "discovery_ABS", "discovery_GOB", "discovery_TAC",
  "heldout_ABS", "heldout_GOB", "heldout_TAC"
)
pathway_keys <- unique(results_long[, c("collection", "gene_set", "description")])
replication_rows <- vector("list", nrow(pathway_keys))
for (index_row in seq_len(nrow(pathway_keys))) {
  key <- pathway_keys[index_row, ]
  local <- results_long[
    results_long$collection == key$collection &
      results_long$gene_set == key$gene_set,
    , drop = FALSE
  ]
  lookup <- local[match(context_names, local$context), , drop = FALSE]
  discovery_common <- lookup[which(lookup$context == "discovery_common"), , drop = FALSE]
  heldout_common <- lookup[which(lookup$context == "heldout_common"), , drop = FALSE]
  state_rows <- lookup[match(state_contexts, lookup$context), , drop = FALSE]
  reference_sign <- if (nrow(discovery_common) == 1L) {
    discovery_common$direction_sign
  } else {
    NA_integer_
  }
  n_state_direction_match <- sum(state_rows$direction_sign == reference_sign, na.rm = TRUE)
  opposite_supported <- any(
    state_rows$direction_sign == -reference_sign & state_rows$fdr <= 0.10,
    na.rm = TRUE
  )
  replicated <- nrow(discovery_common) == 1L && nrow(heldout_common) == 1L &&
    is.finite(discovery_common$fdr) && discovery_common$fdr <= 0.05 &&
    is.finite(heldout_common$fdr) && heldout_common$fdr <= 0.10 &&
    heldout_common$direction_sign == reference_sign &&
    n_state_direction_match >= 5L && !opposite_supported
  scalar <- function(frame, column, default = NA) {
    if (nrow(frame) == 1L) frame[[column]][1L] else default
  }
  replication_rows[[index_row]] <- data.frame(
    collection = key$collection,
    gene_set = key$gene_set,
    description = key$description,
    discovery_common_direction = scalar(discovery_common, "direction", NA_character_),
    discovery_common_p = scalar(discovery_common, "p_value", NA_real_),
    discovery_common_fdr = scalar(discovery_common, "fdr", NA_real_),
    heldout_common_direction = scalar(heldout_common, "direction", NA_character_),
    heldout_common_p = scalar(heldout_common, "p_value", NA_real_),
    heldout_common_fdr = scalar(heldout_common, "fdr", NA_real_),
    n_state_direction_match = n_state_direction_match,
    opposite_supported_state = opposite_supported,
    replicated = replicated,
    stringsAsFactors = FALSE
  )
}
replication <- do.call(rbind, replication_rows)

message("Extracting discovery leading edges for replicated pathways")
replicated_keys <- replication[replication$replicated, c("collection", "gene_set")]
leading_records <- list()
leading_counter <- 0L
for (collection_label in unique(replicated_keys$collection)) {
  selected <- replicated_keys$gene_set[replicated_keys$collection == collection_label]
  if (!length(selected)) {
    next
  }
  resource <- resources[[collection_label]]
  pathways <- resource$gene_sets[selected]
  fg <- fgseaMultilevel(
    pathways = pathways,
    stats = sort(contexts$discovery_common, decreasing = TRUE),
    minSize = 10L,
    maxSize = 500L,
    eps = 0,
    scoreType = "std"
  )
  if (!nrow(fg)) {
    next
  }
  for (row_index in seq_len(nrow(fg))) {
    leading_counter <- leading_counter + 1L
    leading_records[[leading_counter]] <- data.frame(
      collection = collection_label,
      gene_set = fg$pathway[row_index],
      fgsea_nes = fg$NES[row_index],
      fgsea_p_value = fg$pval[row_index],
      fgsea_fdr = fg$padj[row_index],
      leading_edge = paste(fg$leadingEdge[[row_index]], collapse = ";"),
      leading_edge_size = length(fg$leadingEdge[[row_index]]),
      stringsAsFactors = FALSE
    )
  }
}
if (length(leading_records)) {
  leading <- do.call(rbind, leading_records)
} else {
  leading <- data.frame(
    collection = character(), gene_set = character(), fgsea_nes = numeric(),
    fgsea_p_value = numeric(), fgsea_fdr = numeric(), leading_edge = character(),
    leading_edge_size = integer(), stringsAsFactors = FALSE
  )
}
replication <- merge(
  replication, leading,
  by = c("collection", "gene_set"), all.x = TRUE, sort = FALSE
)

message("Building the replicated leading-edge similarity graph")
nodes <- replication[replication$replicated, , drop = FALSE]
nodes$node_id <- paste(nodes$collection, nodes$gene_set, sep = "::")
leading_lists <- strsplit(ifelse(is.na(nodes$leading_edge), "", nodes$leading_edge), ";", fixed = TRUE)
names(leading_lists) <- nodes$node_id
leading_lists <- lapply(leading_lists, function(value) value[nzchar(value)])
edge_records <- list()
edge_counter <- 0L
if (nrow(nodes) >= 2L) {
  pairs <- combn(nodes$node_id, 2L, simplify = FALSE)
  for (pair in pairs) {
    union_size <- length(union(leading_lists[[pair[1L]]], leading_lists[[pair[2L]]]))
    if (!union_size) {
      next
    }
    jaccard <- length(intersect(
      leading_lists[[pair[1L]]], leading_lists[[pair[2L]]]
    )) / union_size
    if (is.finite(jaccard) && jaccard >= 0.25) {
      edge_counter <- edge_counter + 1L
      edge_records[[edge_counter]] <- data.frame(
        from = pair[1L], to = pair[2L], jaccard = jaccard,
        stringsAsFactors = FALSE
      )
    }
  }
}
if (length(edge_records)) {
  edges <- do.call(rbind, edge_records)
} else {
  edges <- data.frame(from = character(), to = character(), jaccard = numeric())
}

if (nrow(nodes)) {
  graph <- graph_from_data_frame(
    edges, directed = FALSE,
    vertices = data.frame(name = nodes$node_id, stringsAsFactors = FALSE)
  )
  if (ecount(graph)) {
    set.seed(20260830)
    communities <- cluster_leiden(
      graph, objective_function = "modularity", weights = E(graph)$jaccard
    )
    membership_vector <- membership(communities)
  } else {
    membership_vector <- seq_len(vcount(graph))
    names(membership_vector) <- V(graph)$name
  }
  nodes$community_id <- as.integer(membership_vector[nodes$node_id])
  representatives <- do.call(rbind, lapply(
    split(nodes, nodes$community_id),
    function(frame) frame[which.min(frame$discovery_common_fdr), , drop = FALSE]
  ))
  labels <- setNames(representatives$gene_set, representatives$community_id)
  nodes$community_label <- unname(labels[as.character(nodes$community_id)])
} else {
  nodes$community_id <- integer()
  nodes$community_label <- character()
  representatives <- nodes
}

community_map <- nodes[, c(
  "collection", "gene_set", "node_id", "community_id", "community_label"
), drop = FALSE]
replication <- merge(
  replication, community_map[, c("collection", "gene_set", "community_id", "community_label")],
  by = c("collection", "gene_set"), all.x = TRUE, sort = FALSE
)
results_long <- merge(
  results_long,
  replication[, c(
    "collection", "gene_set", "replicated", "community_id", "community_label"
  )],
  by = c("collection", "gene_set"), all.x = TRUE, sort = FALSE
)
results_long$replicated[is.na(results_long$replicated)] <- FALSE

write.table(
  results_long,
  file.path(out_root, "pathway_replication.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  replication,
  file.path(out_dir, "pathway_replication_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  leading,
  file.path(out_dir, "replicated_pathway_leading_edges.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  edges,
  file.path(out_dir, "replicated_pathway_network_edges.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  community_map,
  file.path(out_dir, "replicated_pathway_communities.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  representatives,
  file.path(out_dir, "replicated_pathway_representatives.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

collection_summary <- do.call(rbind, lapply(
  split(replication, replication$collection),
  function(frame) data.frame(
    collection = unique(frame$collection),
    tested = nrow(frame),
    discovery_fdr_0.05 = sum(frame$discovery_common_fdr <= 0.05, na.rm = TRUE),
    heldout_fdr_0.10 = sum(frame$heldout_common_fdr <= 0.10, na.rm = TRUE),
    replicated = sum(frame$replicated),
    stringsAsFactors = FALSE
  )
))
write.table(
  collection_summary,
  file.path(out_dir, "pathway_replication_collection_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

output_paths <- c(
  file.path(out_root, "pathway_replication.tsv"),
  file.path(out_dir, "pathway_replication_summary.tsv"),
  file.path(out_dir, "replicated_pathway_leading_edges.tsv"),
  file.path(out_dir, "replicated_pathway_network_edges.tsv"),
  file.path(out_dir, "replicated_pathway_communities.tsv"),
  file.path(out_dir, "replicated_pathway_representatives.tsv"),
  file.path(out_dir, "pathway_replication_collection_summary.tsv")
)
manifest <- list(
  analysis = "state_aware_pathway_replication_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = 20260830,
  input_sha256 = as.list(c(observed_hashes, hallmark = sha256(paths$hallmark))),
  contexts = names(contexts),
  gene_set_size = c(min = 10L, max = 500L),
  primary_test = "limma::cameraPR with ranks and inter.gene.cor=0.01",
  leading_edge_method = "fgseaMultilevel; visualisation and network only",
  replication_rule = paste(
    "discovery common FDR<=0.05; heldout common FDR<=0.10;",
    ">=5/6 state directions match; no supported opposite state"
  ),
  network = "leading-edge Jaccard >=0.25; fixed-seed Leiden modularity",
  n_replicated = sum(replication$replicated),
  n_communities = length(unique(na.omit(replication$community_id))),
  package_versions = as.list(vapply(
    c("fgsea", "igraph", "limma", "msigdbr", "digest", "jsonlite"),
    function(package) as.character(packageVersion(package)), character(1)
  )),
  output_sha256 = as.list(setNames(vapply(output_paths, sha256, character(1)), basename(output_paths))),
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "pathway_replication_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)

message(
  "Pathway replication complete: ", sum(replication$replicated),
  " replicated pathways in ",
  length(unique(na.omit(replication$community_id))), " communities"
)
