#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(root, "results", "state_aware_program_v1")
architecture_root <- file.path(result_root, "functional_architecture_v1")
wgcna_root <- file.path(architecture_root, "consensus_wgcna")
audit_root <- file.path(wgcna_root, "audit")
qc_root <- file.path(root, "qc", "functional_architecture_v1")
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  discovery_ranking = file.path(
    result_root, "common_effects", "cross_state_common_effects.tsv.gz"
  ),
  heldout_ranking = file.path(
    result_root, "heldout_validation", "heldout_cross_state_common_effects.tsv.gz"
  ),
  pathway = file.path(architecture_root, "pathway_replication.tsv"),
  pathway_summary = file.path(
    architecture_root, "pathway_replication", "pathway_replication_summary.tsv"
  ),
  modules = file.path(architecture_root, "consensus_modules.tsv"),
  module_summary = file.path(wgcna_root, "module_internal_gate_summary.tsv"),
  preservation = file.path(wgcna_root, "module_preservation.tsv"),
  enrichment = file.path(wgcna_root, "module_rank_enrichment.tsv"),
  stability = file.path(audit_root, "module_parameter_stability.tsv"),
  profiles = file.path(audit_root, "network_profile_audit.tsv"),
  community_overlap = file.path(
    wgcna_root, "module_pathway_community_overlap.tsv"
  ),
  routing_addendum = file.path(
    root, "analysis", "contracts",
    "state_aware_functional_architecture_routing_addendum_v1_2026-08-30.md"
  ),
  pathway_manifest = file.path(
    architecture_root, "pathway_replication", "pathway_replication_manifest.json"
  ),
  wgcna_manifest = file.path(wgcna_root, "consensus_wgcna_manifest.json"),
  replay = file.path(qc_root, "wgcna_primary_replay.tsv")
)

checks <- list()
check_counter <- 0L
add_check <- function(check_id, passed, detail) {
  check_counter <<- check_counter + 1L
  checks[[check_counter]] <<- data.frame(
    check_id = check_id,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

missing <- names(paths)[!file.exists(unlist(paths))]
add_check(
  "required_core_outputs_exist",
  !length(missing),
  if (length(missing)) paste(missing, collapse = ";") else "all present"
)
if (length(missing)) {
  output <- do.call(rbind, checks)
  write.table(
    output, file.path(qc_root, "functional_architecture_v1_audit.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  stop("Functional-architecture audit cannot continue; core outputs are missing")
}
add_check(
  "routing_addendum_frozen",
  digest(paths$routing_addendum, algo = "sha256", file = TRUE) ==
    "45120db0d56cc31e0610e63e1b28f1e09638b5fff78b404f34a0347a7cb62ea2",
  basename(paths$routing_addendum)
)

discovery <- read.delim(paths$discovery_ranking, check.names = FALSE)
heldout <- read.delim(paths$heldout_ranking, check.names = FALSE)
add_check(
  "frozen_testable_universe",
  nrow(discovery) == 8221L && length(unique(discovery$gene)) == 8221L,
  paste("genes", nrow(discovery))
)
add_check(
  "frozen_high_confidence_programme",
  sum(discovery$strict_state_shared) == 1843L &&
    sum(discovery$strict_state_shared & discovery$shared_direction == "up") == 884L &&
    sum(discovery$strict_state_shared & discovery$shared_direction == "down") == 959L,
  paste(
    "total", sum(discovery$strict_state_shared),
    "up", sum(discovery$strict_state_shared & discovery$shared_direction == "up"),
    "down", sum(discovery$strict_state_shared & discovery$shared_direction == "down")
  )
)
add_check(
  "heldout_ranking_not_identical_to_discovery",
  !identical(discovery$common_z, heldout$common_z),
  "discovery and held-out common effects remain separate"
)

pathway <- read.delim(paths$pathway, check.names = FALSE)
pathway_summary <- read.delim(paths$pathway_summary, check.names = FALSE)
expected_contexts <- c(
  "discovery_common", "discovery_ABS", "discovery_GOB", "discovery_TAC",
  "heldout_common", "heldout_ABS", "heldout_GOB", "heldout_TAC"
)
add_check(
  "pathway_eight_frozen_contexts",
  setequal(unique(pathway$context), expected_contexts),
  paste(sort(unique(pathway$context)), collapse = ";")
)
recalculated_replication <- with(
  pathway_summary,
  discovery_common_fdr <= 0.05 & heldout_common_fdr <= 0.10 &
    discovery_common_direction == heldout_common_direction &
    n_state_direction_match >= 5L & !opposite_supported_state
)
recalculated_replication[is.na(recalculated_replication)] <- FALSE
add_check(
  "pathway_replication_rule_reproduced",
  identical(as.logical(pathway_summary$replicated), recalculated_replication),
  paste("replicated", sum(pathway_summary$replicated))
)
replicated <- pathway_summary[pathway_summary$replicated, , drop = FALSE]
add_check(
  "replicated_pathways_have_complete_state_direction",
  all(replicated$n_state_direction_match >= 5L) &&
    !any(replicated$opposite_supported_state),
  paste("eligible pathways", nrow(replicated))
)

modules <- read.delim(paths$modules, check.names = FALSE)
module_summary <- read.delim(paths$module_summary, check.names = FALSE)
preservation <- read.delim(paths$preservation, check.names = FALSE)
enrichment <- read.delim(paths$enrichment, check.names = FALSE)
stability <- read.delim(paths$stability, check.names = FALSE)
profiles <- read.delim(paths$profiles, check.names = FALSE)

add_check(
  "one_module_assignment_per_gene",
  !anyDuplicated(modules$gene) && all(nzchar(modules$module)),
  paste("assigned genes", nrow(modules))
)
add_check(
  "network_uses_unselected_testable_universe",
  nrow(modules) >= 0.95 * 8221L &&
    all(modules$gene %in% discovery$gene),
  paste("network genes", nrow(modules), "of 8221")
)
profile_key <- paste(
  profiles$partition, profiles$state, profiles$donor_id, profiles$route,
  sep = "__"
)
add_check(
  "one_donor_route_profile_per_state",
  !anyDuplicated(profile_key),
  paste("network profiles", nrow(profiles))
)
discovery_donors <- unique(profiles$donor_id[profiles$partition == "discovery"])
heldout_donors <- unique(profiles$donor_id[profiles$partition == "heldout"])
add_check(
  "discovery_heldout_donor_disjoint",
  !length(intersect(discovery_donors, heldout_donors)),
  paste("overlap", paste(intersect(discovery_donors, heldout_donors), collapse = ";"))
)

recomputed_internal <- vapply(module_summary$module, function(module) {
  local_preservation <- preservation[preservation$module == module, , drop = FALSE]
  local_enrichment <- enrichment[enrichment$module == module, , drop = FALSE]
  heldout_common <- local_enrichment[local_enrichment$context == "heldout_common", ]
  heldout_states <- local_enrichment[local_enrichment$context %in%
    c("heldout_ABS", "heldout_GOB", "heldout_TAC"), ]
  local_stability <- stability[stability$module == module, , drop = FALSE]
  preservation_pass <- nrow(local_preservation) == 3L &&
    sum(local_preservation$zsummary >= 2, na.rm = TRUE) >= 2L &&
    all(local_preservation$zsummary >= 0)
  enrichment_pass <- nrow(heldout_common) == 1L &&
    heldout_common$fdr <= 0.10 && nrow(heldout_states) == 3L &&
    all(heldout_states$direction_sign == heldout_common$direction_sign)
  stability_pass <- nrow(local_stability) == 2L &&
    all(local_stability$maximum_jaccard >= 0.50)
  biological_pass <- module_summary$biological_overlap_pass[
    match(module, module_summary$module)
  ]
  isTRUE(preservation_pass && enrichment_pass && stability_pass && biological_pass)
}, logical(1))
add_check(
  "module_internal_gate_reproduced",
  identical(
    unname(as.logical(module_summary$internal_gate_pass)),
    unname(recomputed_internal)
  ),
  paste("internally eligible", sum(recomputed_internal))
)

replay <- read.delim(paths$replay, check.names = FALSE)
add_check(
  "deterministic_wgcna_replay_passed",
  nrow(replay) == 4L && all(as.logical(replay$passed)),
  paste("checks", nrow(replay), "passed", sum(as.logical(replay$passed)))
)

pathway_manifest <- fromJSON(paths$pathway_manifest, simplifyVector = FALSE)
wgcna_manifest <- fromJSON(paths$wgcna_manifest, simplifyVector = FALSE)
add_check(
  "fixed_random_seed_recorded",
  isTRUE(as.integer(pathway_manifest$random_seed) == 20260830L) &&
    isTRUE(as.integer(wgcna_manifest$random_seed) == 20260830L),
  "seed 20260830"
)
add_check(
  "preservation_permutations_recorded",
  isTRUE(as.integer(wgcna_manifest$parameters$preservation_permutations) == 1000L),
  paste("n", wgcna_manifest$parameters$preservation_permutations)
)

validation_path <- file.path(architecture_root, "module_validation.tsv")
protein_path <- file.path(architecture_root, "protein_priorities.tsv")
if (file.exists(validation_path)) {
  validation <- read.delim(validation_path, check.names = FALSE)
  add_check(
    "module_validation_one_row_per_module",
    !anyDuplicated(validation$module) && setequal(validation$module, module_summary$module),
    paste("rows", nrow(validation))
  )
  if ("routing_status" %in% names(validation)) {
    final_columns <- c(
      "internal_gate_pass", "external_gate_pass", "technical_gate_pass",
      "functional_community_increment_pass", "perturbation_gate_pass",
      "n_sentinel_proteins", "n_regulatory_nodes",
      "interpretive_increment_pass", "routing_status"
    )
    add_check(
      "final_routing_columns_complete",
      all(final_columns %in% names(validation)),
      paste(setdiff(final_columns, names(validation)), collapse = ";")
    )
    if (all(final_columns %in% names(validation))) {
      community_overlap <- read.delim(paths$community_overlap, check.names = FALSE)
      replicated_direction <- unique(
        pathway_summary[pathway_summary$replicated,
          c("community_id", "discovery_common_direction")]
      )
      direction_counts <- table(
        replicated_direction$community_id,
        replicated_direction$discovery_common_direction
      )
      add_check(
        "replicated_pathway_communities_direction_coherent",
        all(rowSums(direction_counts > 0) == 1L),
        paste("communities", nrow(direction_counts))
      )
      recomputed_functional <- vapply(validation$module, function(module) {
        local <- community_overlap[
          community_overlap$module == module &
            is.finite(community_overlap$fdr) & community_overlap$fdr <= 0.05,
          , drop = FALSE
        ]
        if (!nrow(local)) return(FALSE)
        local <- merge(
          local, replicated_direction,
          by = "community_id", all.x = TRUE
        )
        expected_direction <- validation$heldout_direction[
          match(module, validation$module)
        ]
        any(local$discovery_common_direction == expected_direction, na.rm = TRUE)
      }, logical(1))
      add_check(
        "functional_community_increment_reproduced",
        identical(
          unname(as.logical(validation$functional_community_increment_pass)),
          unname(recomputed_functional)
        ),
        paste("functional increments", sum(recomputed_functional))
      )
      recomputed_increment <- with(
        validation,
        as.logical(functional_community_increment_pass) |
          as.logical(perturbation_gate_pass) |
          n_sentinel_proteins > 0L | n_regulatory_nodes > 0L
      )
      core <- with(
        validation,
        as.logical(internal_gate_pass) & as.logical(external_gate_pass) &
          as.logical(technical_gate_pass)
      )
      expected_routing <- ifelse(
        core & recomputed_increment, "Main",
        ifelse(core, "Supplement", "Audit")
      )
      add_check(
        "interpretive_increment_reproduced",
        identical(
          as.logical(validation$interpretive_increment_pass),
          recomputed_increment
        ),
        paste("increments", sum(recomputed_increment))
      )
      add_check(
        "final_module_routing_reproduced",
        identical(as.character(validation$routing_status), expected_routing),
        paste(sort(unique(expected_routing)), collapse = ";")
      )
    }
    add_check(
      "module_routing_uses_fixed_categories",
      all(validation$routing_status %in% c("Main", "Supplement", "Audit")),
      paste(sort(unique(validation$routing_status)), collapse = ";")
    )
  }
}
if (file.exists(protein_path) && file.info(protein_path)$size > 0) {
  protein <- tryCatch(
    read.delim(protein_path, check.names = FALSE),
    error = function(error) data.frame()
  )
  add_check(
    "protein_roles_are_separated",
    !nrow(protein) || all(protein$priority_role %in%
      c("measurement_sentinel", "regulatory_node")),
    paste("n priorities", nrow(protein))
  )
}

output <- do.call(rbind, checks)
write.table(
  output, file.path(qc_root, "functional_architecture_v1_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
manifest <- list(
  analysis = "audit_functional_architecture_v1",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  passed = all(output$passed),
  n_checks = nrow(output),
  n_failed = sum(!output$passed),
  failed_checks = output$check_id[!output$passed],
  input_sha256 = as.list(setNames(
    vapply(paths, function(path) digest(path, algo = "sha256", file = TRUE), character(1)),
    names(paths)
  )),
  session_info = capture.output(sessionInfo())
)
write_json(
  manifest, file.path(qc_root, "functional_architecture_v1_audit_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)

if (!all(output$passed)) {
  stop(
    "Functional-architecture audit failed: ",
    paste(output$check_id[!output$passed], collapse = ", ")
  )
}
message("Functional-architecture audit passed: ", nrow(output), " checks")
