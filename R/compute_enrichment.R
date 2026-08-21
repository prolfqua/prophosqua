#' Compute the Enrichment Analyses of the PTM Pipeline
#'
#' The three enrichment analyses — PTM-SEA against PTMsigDB, GSEA against
#' kinase-library motifs, and motif enrichment analysis of the pre-ranked site
#' lists — used to be computed inside the reports that display them. Permutation
#' tests are the slowest part of the pipeline, so a fix to a caption cost a full
#' recomputation, and their result workbooks were written from inside a
#' conditional chunk and could therefore not be declared as rule outputs.
#'
#' These functions do the computing. Each returns everything its report needs,
#' under the names the report uses, and each **always** produces a result —
#' empty tables and `has_results = FALSE` when nothing was enriched, rather than
#' no result at all. That is what makes the workbooks declarable: an enrichment
#' that found nothing is a fact to record, not a missing file.
#'
#' @name compute_enrichment
#' @keywords internal
NULL

#' Compute PTM-SEA against PTMsigDB
#'
#' Ranks the sites of each contrast by `stat_column`, matches them to PTMsigDB
#' signatures on their trimmed flanking sequence, and runs pre-ranked GSEA.
#'
#' `trim_to` must match the width the signature database was trimmed to, or
#' nothing matches: the site identifier PTMsigDB keys on *is* the flanking
#' sequence.
#'
#' @param xlsx_file Combined `PTM_results.xlsx`.
#' @param sheet Sheet to read.
#' @param stat_column Per-site statistic the sites are ranked on.
#' @param ptmsigdb_file Filtered PTMsigDB, as `.rds` or `.gmt`.
#' @param trim_to Flanking width the database uses.
#' @param min_size,max_size Set-size filter for testable signatures.
#' @param n_perm Permutations of the enrichment test.
#' @return A list with `results` (one clusterProfiler object per contrast),
#'   `all_clean` (all results as one table), `pathways`, the summary tables the
#'   report shows, and `has_results`.
#' @export
#' @examples
#' # Needs a combined workbook and a filtered PTMsigDB.
#' \dontrun{
#' res <- compute_ptmsea(
#'   "PTM_results.xlsx", "DPA", "statistic.site",
#'   "data/ptmsigdb/ptmsigdb_filtered_KINASE-PSP_15mer.rds"
#' )
#' res$results_info
#' }
compute_ptmsea <- function(
  xlsx_file,
  sheet,
  stat_column,
  ptmsigdb_file,
  trim_to = 15,
  min_size = 10,
  max_size = 500,
  n_perm = 1000
) {
  data <- readxl::read_xlsx(xlsx_file, sheet = sheet)
  data <- canonicalize_sequence_window(data)

  pathways <- read_ptmsigdb(ptmsigdb_file)

  data_info <- dplyr::tibble(
    Property = c("Mode", "Sheet", "Stat Column", "Rows", "Columns", "Contrasts"),
    Value = c(
      basename(xlsx_file),
      sheet,
      stat_column,
      nrow(data),
      ncol(data),
      paste(unique(data$contrast), collapse = ", ")
    )
  )

  ptmsigdb_summary <- dplyr::tibble(
    Property = c(
      "Source File",
      "Total Signatures",
      "KINASE signatures",
      "PATH signatures",
      "Unique site IDs"
    ),
    Value = c(
      basename(ptmsigdb_file),
      length(pathways),
      sum(grepl("^KINASE-", names(pathways))),
      sum(grepl("^PATH-", names(pathways))),
      length(unique(gsub(";[ud]$", "", unlist(pathways))))
    )
  )

  overlap <- ptmsigdb_overlap(data, pathways, trim_to)

  prep <- ptmsea_data_prep(
    data = data,
    stat_column = stat_column,
    seq_window_col = "SequenceWindow",
    contrast_col = "contrast",
    trim_to = as.character(trim_to)
  )
  prep_info <- dplyr::tibble(
    Contrast = names(prep$ranks),
    Sites = purrr::map_int(prep$ranks, length)
  )

  results <- run_ptmsea(
    ranks_list = prep$ranks,
    pathways = pathways,
    min_size = min_size,
    max_size = max_size,
    n_perm = n_perm,
    # Relaxed so that the report can show near-misses; the significance
    # thresholds it draws are applied when plotting, not here.
    pvalueCutoff = 0.25
  )

  results_info <- dplyr::tibble(
    Contrast = names(results),
    `Total Pathways` = purrr::map_int(results, ~ nrow(.x@result)),
    `FDR < 0.1` = purrr::map_int(results, ~ sum(.x@result$p.adjust < 0.1, na.rm = TRUE)),
    `FDR < 0.05` = purrr::map_int(results, ~ sum(.x@result$p.adjust < 0.05, na.rm = TRUE))
  )
  has_results <- sum(results_info$`Total Pathways`) > 0

  # extract_gsea_results() already yields the right columns with zero rows when
  # nothing was enriched, so an empty enrichment produces an empty table rather
  # than no table. That is what lets the workbook be a declared rule output.
  all_clean <- extract_gsea_results(results) |>
    dplyr::mutate(
      pathway = .data$ID,
      pathway_short = substr(
        gsub("^(KINASE|PERT|PATH|DISEASE)-PSP_", "", .data$ID),
        1,
        40
      )
    )

  list(
    results = results,
    all_clean = all_clean,
    pathways = pathways,
    data_info = data_info,
    ptmsigdb_summary = ptmsigdb_summary,
    overlap_stats = overlap$stats,
    n_overlap = overlap$n_overlap,
    n_our_sites = overlap$n_our_sites,
    prep_info = prep_info,
    results_info = results_info,
    has_results = has_results,
    analysis_inputs = list(
      xlsx_file = xlsx_file,
      sheet = sheet,
      stat_column = stat_column,
      ptmsigdb_file = ptmsigdb_file,
      trim_to = trim_to,
      min_size = min_size,
      max_size = max_size,
      n_perm = n_perm
    )
  )
}

#' Compute Kinase-Library GSEA
#'
#' Ranks the sites of each contrast by their site-level statistic and runs
#' pre-ranked GSEA against the kinase-substrate sets a motif scan assigned.
#'
#' @param xlsx_file Combined `PTM_results.xlsx`.
#' @param sheet Sheet to read.
#' @param term2gene_file Kinase-to-site assignment written by the motif scan.
#' @param min_size,max_size Set-size filter for testable kinase sets.
#' @param n_perm Permutations of the enrichment test.
#' @param max_kinase_sets Subsample the kinase sets to at most this many, for a
#'   quick standalone run. `NULL`, the default, tests all of them.
#' @return A list with `gsea_results`, `all_results`, the summary tables the
#'   report shows, and `has_results`.
#' @export
#' @examples
#' # Needs a combined workbook and a term2gene assignment.
#' \dontrun{
#' res <- compute_kinaselib_gsea(
#'   "PTM_results.xlsx", "DPA", "PTM_DPA/KinaseLib/term2gene.csv"
#' )
#' res$gsea_info
#' }
compute_kinaselib_gsea <- function(
  xlsx_file,
  sheet,
  term2gene_file,
  min_size = 15,
  max_size = 5000,
  n_perm = 1000,
  max_kinase_sets = NULL
) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required to run the kinase-library GSEA.", call. = FALSE)
  }

  data <- readxl::read_xlsx(xlsx_file, sheet = sheet)
  data <- canonicalize_sequence_window(data)
  stat_col <- "statistic.site"

  term2gene <- utils::read.csv(term2gene_file, stringsAsFactors = FALSE)

  data_info <- dplyr::tibble(
    Property = c("Rows", "Columns", "Contrasts"),
    Value = c(
      nrow(data),
      ncol(data),
      paste(unique(data$contrast), collapse = ", ")
    )
  )
  kl_info <- dplyr::tibble(
    Property = c("Total assignments", "Unique kinases", "Unique sequences"),
    Value = c(
      nrow(term2gene),
      dplyr::n_distinct(term2gene$term),
      dplyr::n_distinct(term2gene$gene)
    )
  )

  # Kinase Library writes the phosphorylated residue in lower case, e.g.
  # "PETITIRsGPPSPLP". Our sequence windows are upper case throughout, so
  # without this the two never match and every set comes out empty.
  term2gene$gene <- toupper(term2gene$gene)
  term2gene_df <- term2gene |> dplyr::select("term", "gene")

  if (!is.null(max_kinase_sets)) {
    set.seed(42)
    keep <- sample(
      unique(term2gene_df$term),
      min(max_kinase_sets, dplyr::n_distinct(term2gene_df$term))
    )
    term2gene_df <- term2gene_df |> dplyr::filter(.data$term %in% keep)
    message("Subsampled to ", length(keep), " kinase sets")
  }

  our_sequences <- unique(toupper(trimws(data$SequenceWindow)))
  assigned_sequences <- unique(term2gene$gene)
  overlap_seqs <- intersect(our_sequences, assigned_sequences)

  assignment_stats <- dplyr::tibble(
    Metric = c(
      "Phosphosites in differential analysis",
      "Sites with kinase assignments",
      "Sites usable for GSEA",
      "Coverage (%)"
    ),
    Value = c(
      length(our_sequences),
      length(assigned_sequences),
      length(overlap_seqs),
      round(100 * length(overlap_seqs) / length(our_sequences), 1)
    )
  )

  kinase_sizes <- term2gene_df |>
    dplyr::group_by(.data$term) |>
    dplyr::summarize(size = dplyr::n(), .groups = "drop")
  kinase_overlap <- term2gene_df |>
    dplyr::filter(.data$gene %in% our_sequences) |>
    dplyr::group_by(.data$term) |>
    dplyr::summarize(overlap = dplyr::n(), .groups = "drop")

  kinase_stats <- dplyr::tibble(
    Metric = c(
      "Number of kinases",
      "Mean substrates/kinase",
      "Median substrates/kinase",
      "Range (min-max)",
      "Mean substrates in our data",
      "Kinases with >= 15 substrates"
    ),
    Value = c(
      nrow(kinase_sizes),
      round(mean(kinase_sizes$size)),
      stats::median(kinase_sizes$size),
      paste(min(kinase_sizes$size), "-", max(kinase_sizes$size)),
      round(mean(kinase_overlap$overlap, na.rm = TRUE)),
      sum(kinase_overlap$overlap >= 15, na.rm = TRUE)
    )
  )

  ranks <- prepare_gsea_ranks(
    data,
    stat_col = stat_col,
    seq_col = "SequenceWindow",
    contrast_col = "contrast"
  )
  ranks_info <- dplyr::tibble(
    Contrast = names(ranks),
    Sites = purrr::map_int(ranks, length)
  )

  gsea_results <- purrr::map(
    rlang::set_names(names(ranks)),
    function(ct) {
      clusterProfiler::GSEA(
        geneList = ranks[[ct]],
        TERM2GENE = term2gene_df,
        minGSSize = min_size,
        maxGSSize = max_size,
        pvalueCutoff = 0.25,
        nPermSimple = n_perm,
        verbose = FALSE
      )
    }
  )

  gsea_info <- dplyr::tibble(
    Contrast = names(gsea_results),
    `Significant Kinases (FDR < 0.25)` = purrr::map_int(
      gsea_results,
      ~ sum(.x@result$p.adjust < 0.25, na.rm = TRUE)
    )
  )
  has_results <- sum(purrr::map_int(gsea_results, ~ nrow(.x@result))) > 0

  all_results <- purrr::map_dfr(
    rlang::set_names(names(gsea_results)),
    function(ct) {
      dplyr::as_tibble(gsea_results[[ct]]@result) |>
        dplyr::mutate(contrast = ct) |>
        dplyr::select(
          "contrast",
          kinase = "ID",
          "NES",
          "pvalue",
          FDR = "p.adjust",
          "setSize"
        )
    }
  )

  list(
    gsea_results = gsea_results,
    all_results = all_results,
    term2gene = term2gene,
    term2gene_df = term2gene_df,
    n_our_sequences = length(our_sequences),
    n_overlap_seqs = length(overlap_seqs),
    data_info = data_info,
    kl_info = kl_info,
    assignment_stats = assignment_stats,
    kinase_stats = kinase_stats,
    ranks_info = ranks_info,
    gsea_info = gsea_info,
    has_results = has_results,
    analysis_inputs = list(
      xlsx_file = xlsx_file,
      sheet = sheet,
      term2gene_file = term2gene_file,
      stat_col = stat_col,
      min_size = min_size,
      max_size = max_size,
      n_perm = n_perm
    )
  )
}

#' Collect the Motif Enrichment Results of One Analysis
#'
#' Reads the per-contrast `mea_*.csv` files the kinase-library tool wrote and
#' brings them onto the column names the enrichment reports share. Unlike the
#' two GSEA analyses this does no testing of its own — the enrichment already
#' happened, one contrast per file.
#'
#' @param kinaselib_dir Directory holding the `mea_*.csv` files.
#' @return A list with `mea_clean`, its per-contrast `summary_df`, and
#'   `has_results`.
#' @export
#' @examples
#' # Needs the mea_*.csv files of a pipeline run.
#' \dontrun{
#' res <- compute_mea("PTM_DPA/KinaseLib")
#' res$summary_df
#' }
compute_mea <- function(kinaselib_dir) {
  mea_files <- list.files(
    kinaselib_dir,
    pattern = "^mea_.*\\.csv$",
    full.names = TRUE
  )
  if (length(mea_files) == 0) {
    stop("No MEA result files found in: ", kinaselib_dir, call. = FALSE)
  }

  mea_results <- purrr::map_dfr(
    rlang::set_names(mea_files, gsub("^mea_|\\.csv$", "", basename(mea_files))),
    ~ utils::read.csv(.x, stringsAsFactors = FALSE),
    .id = "contrast"
  )
  mea_results <- canonicalize_mea_columns(mea_results)

  message("Loaded ", nrow(mea_results), " results from ", length(mea_files), " contrasts")

  mea_clean <- prepare_enrichment_data(mea_results, "FDR", 0.1)

  summary_df <- mea_clean |>
    dplyr::group_by(.data$contrast) |>
    dplyr::summarize(
      total_kinases = dplyr::n(),
      sig_up = sum(.data$FDR < 0.1 & .data$NES > 0, na.rm = TRUE),
      sig_down = sum(.data$FDR < 0.1 & .data$NES < 0, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    mea_clean = mea_clean,
    summary_df = summary_df,
    n_files = length(mea_files),
    has_results = nrow(mea_clean) > 0,
    analysis_inputs = list(kinaselib_dir = kinaselib_dir)
  )
}

#' Read a Filtered PTMsigDB from RDS or GMT
#'
#' @param ptmsigdb_file Path ending in `.rds` or `.gmt`.
#' @return Named list of signatures.
#' @keywords internal
read_ptmsigdb <- function(ptmsigdb_file) {
  if (!file.exists(ptmsigdb_file)) {
    stop("PTMsigDB file not found: ", ptmsigdb_file, call. = FALSE)
  }
  if (grepl("\\.rds$", ptmsigdb_file)) {
    pathways <- readRDS(ptmsigdb_file)
    message("Loaded ", length(pathways), " pathways from RDS")
    return(pathways)
  }
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea is required to read a GMT file.", call. = FALSE)
  }
  pathways <- fgsea::gmtPathways(ptmsigdb_file)
  message("Loaded ", length(pathways), " pathways from GMT")
  pathways
}

#' How Much of the Data PTMsigDB Annotates
#'
#' The two sides are matched on trimmed flanking sequence, which is the site
#' identifier PTMsigDB keys on. A low overlap bounds how many signatures can be
#' tested at all, and is the first thing to check when few pathways pass the
#' size filter.
#'
#' @param data Result table with a `SequenceWindow` column.
#' @param pathways PTMsigDB signatures.
#' @param trim_to Flanking width to trim to before matching.
#' @return A list with the summary `stats` table and the two counts behind it.
#' @keywords internal
ptmsigdb_overlap <- function(data, pathways, trim_to) {
  our_sequences <- unique(toupper(trimws(data$SequenceWindow)))
  our_site_ids <- paste0(
    purrr::map_chr(our_sequences, ~ trim_flanking_seq(.x, trim_to = trim_to)),
    "-p"
  )
  n_our_sites <- dplyr::n_distinct(our_site_ids)

  ptmsigdb_ids <- unique(gsub(";[ud]$", "", unique(unlist(pathways))))
  n_ptmsigdb_sites <- length(ptmsigdb_ids)
  n_overlap <- length(intersect(unique(our_site_ids), ptmsigdb_ids))

  stats <- dplyr::tibble(
    Metric = c(
      "Our data (unique sequences)",
      "PTMsigDB (unique site IDs)",
      "Overlap",
      "% of our sites in PTMsigDB",
      "% of PTMsigDB sites in our data"
    ),
    Value = c(
      n_our_sites,
      n_ptmsigdb_sites,
      n_overlap,
      round(100 * n_overlap / n_our_sites, 2),
      round(100 * n_overlap / n_ptmsigdb_sites, 2)
    )
  )

  list(stats = stats, n_overlap = n_overlap, n_our_sites = n_our_sites)
}

#' Bring MEA Output onto the Shared Enrichment Column Names
#'
#' The kinase-library tool writes `Kinase`, `p-value` and a `Subs fraction`
#' field of the form `85/471`. The enrichment reports expect `kinase`, `pvalue`,
#' and the fraction split into its leading-edge count and set size.
#'
#' @param mea_results Concatenated MEA results.
#' @return The same table, renamed and with the fraction split.
#' @keywords internal
canonicalize_mea_columns <- function(mea_results) {
  if ("Kinase" %in% names(mea_results)) {
    mea_results <- dplyr::rename(mea_results, kinase = "Kinase")
  }
  if ("p.value" %in% names(mea_results)) {
    mea_results <- dplyr::rename(mea_results, pvalue = "p.value")
  } else if ("p-value" %in% names(mea_results)) {
    mea_results <- dplyr::rename(mea_results, pvalue = "p-value")
  }

  subs_col <- intersect(c("Subs.fraction", "Subs fraction"), names(mea_results))
  if (length(subs_col) > 0) {
    mea_results <- mea_results |>
      dplyr::mutate(
        n_leading = as.numeric(sub("/.*", "", .data[[subs_col[1]]])),
        set_size = as.numeric(sub(".*/", "", .data[[subs_col[1]]]))
      )
  }
  mea_results
}

#' Compute the Enrichment Analyses from the Bundled Example Data
#'
#' The three enrichment reports can be rendered without a pipeline run, against
#' the example data the package ships. These helpers assemble that input and
#' then call the same compute function the pipeline calls, so a standalone
#' render exercises the real code path rather than a parallel one.
#'
#' @name compute_enrichment_example
#' @keywords internal
NULL

#' @rdname compute_enrichment_example
compute_ptmsea_example <- function() {
  compute_ptmsea(
    xlsx_file = example_results_workbook(),
    sheet = "DPA",
    stat_column = "statistic.site",
    ptmsigdb_file = unzip_extdata(
      "ptmsigdb_kinase.rds.zip",
      "ptmsigdb_filtered_KINASE_15mer.rds"
    ),
    n_perm = 100
  )
}

#' @rdname compute_enrichment_example
compute_kinaselib_gsea_example <- function() {
  compute_kinaselib_gsea(
    xlsx_file = example_results_workbook(),
    sheet = "DPA",
    term2gene_file = unzip_extdata("term2gene.csv.zip", "term2gene.csv"),
    n_perm = 100,
    max_kinase_sets = 50
  )
}

#' @rdname compute_enrichment_example
compute_mea_example <- function() {
  zip <- system.file("extdata", "mea_results.zip", package = "prophosqua")
  if (!file.exists(zip)) {
    stop("Bundled MEA example data not found.", call. = FALSE)
  }
  dir <- file.path(tempdir(), "mea_example")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(zip, exdir = dir)
  compute_mea(dir)
}

#' Write the Example Result Table as a Workbook
#'
#' The compute functions take a workbook path, because that is what the pipeline
#' has. The example data is an R object, so it is written out once per session
#' to give the standalone path the same input shape.
#'
#' @return Path to the workbook.
#' @keywords internal
example_results_workbook <- function() {
  path <- file.path(tempdir(), "prophosqua_example_results.xlsx")
  if (file.exists(path)) {
    return(path)
  }
  example <- get(utils::data(
    "combined_test_diff_example",
    package = "prophosqua",
    envir = environment()
  ))
  writexl::write_xlsx(list(DPA = example), path)
  path
}

#' Unpack One Bundled Example File
#'
#' @param zip File name under `inst/extdata`.
#' @param member File expected inside it.
#' @return Path to the unpacked file.
#' @keywords internal
unzip_extdata <- function(zip, member) {
  archive <- system.file("extdata", zip, package = "prophosqua")
  if (!file.exists(archive)) {
    stop("Bundled example data not found: ", zip, call. = FALSE)
  }
  dir <- file.path(tempdir(), sub("\\.zip$", "", zip))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(archive, exdir = dir)
  path <- file.path(dir, member)
  if (!file.exists(path)) {
    stop("Bundled archive ", zip, " does not contain ", member, call. = FALSE)
  }
  path
}
