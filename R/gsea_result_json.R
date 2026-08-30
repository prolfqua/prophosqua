#' Export Enrichment Results as string_gsea GSEAResult JSON
#'
#' The stringdbpy package (string_gsea) defines one JSON structure for
#' enrichment results: per contrast a shared pool of ranked items and, per
#' category, terms that reference pool entries by id. Writing the PTM
#' enrichments (PTM-SEA, kinase-library GSEA, MEA) in that same structure gives
#' every downstream consumer — ptm3d, clusterProfiler channelling, reporting —
#' one format instead of three.
#'
#' The ranked items of the PTM enrichments are flanking sequence windows, not
#' site identifiers. Pool ids are therefore the canonical upper-case sequence
#' window (the form the differential table uses); the form actually submitted
#' to the enrichment — with the PTMsigDB `-p` suffix, or the kinase-library
#' lower-cased phospho residue — is kept in `input_label`. Each term carries
#' the full mapped membership in `gene_ids` and, as an extension to the STRING
#' structure, the leading edge in `leading_edge_ids`.
#'
#' @name gsea_result_json
#' @keywords internal
NULL

#' Canonicalize Window Ids the Way the Differential Table Writes Them
#'
#' @param ids Character vector of ranked-item ids.
#' @return The ids upper-cased, trimmed, and with a trailing `-p` removed.
#' @keywords internal
normalize_window_id <- function(ids) {
  toupper(sub("-p$", "", trimws(as.character(ids))))
}

#' An Empty Named List, Serialized by jsonlite as `{}` Rather Than `[]`
#' @keywords internal
empty_named_list <- function() {
  structure(list(), names = character(0))
}

#' Build the Gene Pool of One Contrast from Its Ranked List
#'
#' @param ranks Named numeric vector, sorted descending, names are the
#'   submitted ranked-item ids.
#' @return Named list of pool entries keyed by canonical window id.
#' @keywords internal
gene_pool_from_ranks <- function(ranks) {
  ids <- normalize_window_id(names(ranks))
  entries <- purrr::map(seq_along(ranks), function(i) {
    list(
      protein_id = ids[[i]],
      label = ids[[i]],
      input_label = names(ranks)[[i]],
      input_value = unname(ranks[[i]]),
      rank = i
    )
  })
  rlang::set_names(entries, ids)
}

#' One Term Entry of the GSEAResult Structure
#'
#' @param term_id,category,description Identity of the term.
#' @param nes Normalized enrichment score.
#' @param fdr Adjusted p-value.
#' @param method Enrichment method label.
#' @param gene_ids Canonical ids of the full mapped membership.
#' @param leading_edge_ids Canonical ids of the leading edge.
#' @param genes_in_set Size of the full (unmapped) set.
#' @return A term as a named list.
#' @keywords internal
gsea_term_entry <- function(
  term_id,
  category,
  description,
  nes,
  fdr,
  method,
  gene_ids,
  leading_edge_ids,
  genes_in_set
) {
  direction <- if (is.na(nes) || nes == 0) {
    "both ends"
  } else if (nes > 0) {
    "top"
  } else {
    "bottom"
  }
  list(
    term_id = term_id,
    category = category,
    description = description,
    enrichment_score = nes,
    direction = direction,
    fdr = fdr,
    method = method,
    genes_mapped = length(gene_ids),
    genes_in_set = genes_in_set,
    gene_ids = as.list(gene_ids),
    leading_edge_ids = as.list(leading_edge_ids)
  )
}

#' Assemble the Per-Contrast Blocks of a GSEAResult
#'
#' @param contrast Contrast name.
#' @param pool Named list from [gene_pool_from_ranks()].
#' @param terms List of term entries for the one category.
#' @param category Category name.
#' @return The `data` block of one contrast.
#' @keywords internal
gsea_contrast_entry <- function(contrast, pool, terms, category) {
  categories <- list()
  categories[[category]] <- list(
    category = category,
    contrast = contrast,
    terms = terms
  )
  list(
    contrast = contrast,
    gene_pool = if (length(pool)) pool else empty_named_list(),
    categories = categories
  )
}

#' Convert clusterProfiler GSEA Results to the GSEAResult Structure
#'
#' Converts the per-contrast `gseaResult` objects that [compute_ptmsea()] and
#' [compute_kinaselib_gsea()] return into the string_gsea `GSEAResult`
#' structure: the ranked list (`@geneList`) becomes the contrast's gene pool,
#' each result row becomes a term whose `gene_ids` are the full mapped set
#' membership (`@geneSets` intersected with the ranked list) and whose
#' `leading_edge_ids` are the `core_enrichment` sites.
#'
#' @param results Named list of `gseaResult` objects, one per contrast.
#' @param category Category name recorded on every term, e.g. `"PTM-SEA"`.
#' @param method Method label recorded on every term.
#' @return A GSEAResult as a named list with `data` and `rank_lists`, ready for
#'   [write_gsea_result_json()].
#' @export
#' @examples
#' \dontrun{
#' res <- compute_ptmsea_example()
#' json_data <- gsea_result_data(res$results, category = "PTM-SEA")
#' write_gsea_result_json(json_data, "PTMSEA_DPA_results.json")
#' }
gsea_result_data <- function(results, category, method = "fgsea") {
  data <- purrr::imap(results, function(res, contrast) {
    ranks <- methods::slot(res, "geneList")
    pool <- gene_pool_from_ranks(ranks)
    gene_sets <- methods::slot(res, "geneSets")
    result <- methods::slot(res, "result")

    terms <- purrr::map(seq_len(nrow(result)), function(i) {
      row <- result[i, ]
      full_set <- gene_sets[[row$ID]]
      mapped <- normalize_window_id(intersect(full_set, names(ranks)))
      leading <- normalize_window_id(strsplit(row$core_enrichment, "/")[[1]])
      gsea_term_entry(
        term_id = row$ID,
        category = category,
        description = row$Description,
        nes = row$NES,
        fdr = row$p.adjust,
        method = method,
        gene_ids = mapped,
        leading_edge_ids = intersect(leading, mapped),
        genes_in_set = length(full_set)
      )
    })
    gsea_contrast_entry(contrast, pool, terms, category)
  })

  rank_lists <- purrr::imap(results, function(res, contrast) {
    ranks <- methods::slot(res, "geneList")
    list(contrast = contrast, entries = as.list(ranks))
  })

  list(data = data, rank_lists = rank_lists)
}

#' Convert Collected MEA Results to the GSEAResult Structure
#'
#' The motif enrichment ran outside R (kinase-library CLI), so its inputs are
#' reassembled: the ranked lists come from the `.rnk` files the pipeline wrote
#' (see [read_mea_ranks()]), the full set membership from the `term2gene.csv`
#' motif-scan assignment, and the leading edge from the `Leading substrates`
#' column of the collected results.
#'
#' @param mea_clean The `mea_clean` table of [compute_mea()].
#' @param ranks Named list of ranked vectors, one per contrast.
#' @param term2gene Motif-scan assignment with `term` and `gene` columns.
#' @param category Category name recorded on every term.
#' @return A GSEAResult as a named list with `data` and `rank_lists`.
#' @export
#' @examples
#' \dontrun{
#' res <- compute_mea("PTM_DPA/KinaseLib")
#' ranks <- read_mea_ranks("PTM_DPA/KinaseLib")
#' term2gene <- read.csv("PTM_DPA/KinaseLib/term2gene.csv")
#' json_data <- mea_gsea_result_data(res$mea_clean, ranks, term2gene)
#' }
mea_gsea_result_data <- function(mea_clean, ranks, term2gene, category = "MEA") {
  term2gene$gene <- normalize_window_id(term2gene$gene)
  members_by_term <- split(term2gene$gene, term2gene$term)
  leading_col <- intersect(c("Leading.substrates", "Leading substrates"), names(mea_clean))[1]

  data <- purrr::imap(ranks, function(contrast_ranks, contrast) {
    pool <- gene_pool_from_ranks(contrast_ranks)
    pool_ids <- names(pool)
    contrast_terms <- mea_clean[mea_clean$contrast == contrast, ]

    terms <- purrr::map(seq_len(nrow(contrast_terms)), function(i) {
      row <- contrast_terms[i, ]
      full_set <- members_by_term[[row$kinase]]
      mapped <- intersect(full_set, pool_ids)
      leading <- normalize_window_id(strsplit(row[[leading_col]], ";")[[1]])
      gsea_term_entry(
        term_id = row$kinase,
        category = category,
        description = row$kinase,
        nes = row$NES,
        fdr = row$FDR,
        method = "mea",
        gene_ids = mapped,
        leading_edge_ids = intersect(leading, mapped),
        genes_in_set = length(full_set)
      )
    })
    gsea_contrast_entry(contrast, pool, terms, category)
  })

  rank_lists <- purrr::imap(ranks, function(contrast_ranks, contrast) {
    list(contrast = contrast, entries = as.list(contrast_ranks))
  })

  list(data = data, rank_lists = rank_lists)
}

#' Read the Per-Contrast Ranked Lists the Pipeline Wrote for MEA
#'
#' The pipeline writes one `<sheet>_MEA_<contrast>.rnk` per contrast next to
#' the `mea_*.csv` results: tab-separated with a header row (`SequenceWindow`,
#' `statistic.site`), sorted descending. Reading it headerless would turn the
#' header into a rank entry and coerce every statistic to a string.
#'
#' @param kinaselib_dir Directory holding the `.rnk` files.
#' @return Named list of named numeric vectors, one per contrast.
#' @export
#' @examples
#' \dontrun{
#' ranks <- read_mea_ranks("PTM_DPA/KinaseLib")
#' names(ranks)
#' }
read_mea_ranks <- function(kinaselib_dir) {
  rnk_files <- list.files(kinaselib_dir, pattern = "_MEA_.*\\.rnk$", full.names = TRUE)
  if (length(rnk_files) == 0) {
    stop("No *_MEA_*.rnk files found in: ", kinaselib_dir, call. = FALSE)
  }
  contrasts <- sub("\\.rnk$", "", sub("^.*_MEA_", "", basename(rnk_files)))
  purrr::map(
    rlang::set_names(rnk_files, contrasts),
    function(path) {
      tab <- utils::read.delim(path, header = TRUE)
      rlang::set_names(as.numeric(tab[[2]]), as.character(tab[[1]]))
    }
  )
}

#' Write a GSEAResult Structure as JSON
#'
#' @param gsea_result A list with `data` and `rank_lists`, as produced by
#'   [gsea_result_data()] or [mea_gsea_result_data()].
#' @param path Output file.
#' @return `path`, invisibly.
#' @export
#' @examples
#' \dontrun{
#' write_gsea_result_json(json_data, "PTMSEA_DPA_results.json")
#' }
write_gsea_result_json <- function(gsea_result, path) {
  jsonlite::write_json(gsea_result, path, auto_unbox = TRUE, digits = NA)
  invisible(path)
}
