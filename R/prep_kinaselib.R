#' Prepare Kinase-Library Inputs from the Combined PTM Results
#'
#' Writes the two file kinds the kinase-library tool consumes: one list of
#' unique sequence windows for motif scanning, and one ranked file per contrast
#' for motif enrichment analysis.
#'
#' The workbook, sheet and ranking statistic are the same ones PTM-SEA and the
#' KinaseLib GSEA read, so all three enrichment analyses rank the sites
#' identically and their results can be compared.
#'
#' @param xlsx_file The combined `PTM_results.xlsx`.
#' @param output_dir Where to write, normally the analysis `KinaseLib` directory.
#' @param analysis_type `DPA`, `DPU` or `CF`; names the output files.
#' @param sheet Sheet of `xlsx_file` to read.
#' @param stat_column Column to rank the sites on.
#' @return Invisibly, a character vector of the files written.
#' @export
#' @examples
#' # Needs a combined workbook from a pipeline run.
#' \dontrun{
#' prep_kinaselib_inputs(
#'   "PTM_results.xlsx", "PTM_DPA/KinaseLib",
#'   analysis_type = "DPA", sheet = "DPA", stat_column = "statistic.site"
#' )
#' }
prep_kinaselib_inputs <- function(xlsx_file, output_dir, analysis_type, sheet, stat_column) {
  message("=== Preparing kinase-library inputs ===")
  message("Input file: ", xlsx_file)
  message("Output dir: ", output_dir)
  message("Analysis type: ", analysis_type)
  message("Ranking on: ", stat_column)

  message("Reading sheet: ", sheet)
  data <- readxl::read_xlsx(xlsx_file, sheet = sheet)
  message("  Loaded ", nrow(data), " rows")

  data <- canonicalize_sequence_window(data)

  if (!stat_column %in% names(data)) {
    stop(
      "Statistic column '",
      stat_column,
      "' not found. Available: ",
      paste(utils::head(names(data), 10), collapse = ", "),
      call. = FALSE
    )
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  data_clean <- filter_sequence_windows(data)
  message("  After filtering: ", nrow(data_clean), " rows with valid SequenceWindow")

  seqwindows <- data_clean |>
    dplyr::select("SequenceWindow") |>
    dplyr::distinct() |>
    dplyr::arrange(.data$SequenceWindow)

  seqwindows_file <- file.path(
    output_dir,
    paste0(analysis_type, "_seqwindows.tsv")
  )
  utils::write.table(
    seqwindows,
    seqwindows_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  message("Wrote ", nrow(seqwindows), " unique sequences to: ", seqwindows_file)

  contrasts <- unique(data_clean$contrast)
  message("Found ", length(contrasts), " contrasts")

  rnk_files <- character()
  for (ctr in contrasts) {
    ctr_data <- rank_sites_for_mea(data_clean, stat_column, ctr)

    # The contrast name reaches a file name, and a contrast is free to carry
    # characters a file name is not.
    rnk_file <- file.path(
      output_dir,
      paste0(
        analysis_type,
        "_MEA_",
        gsub("[^A-Za-z0-9_-]", "_", ctr),
        ".rnk"
      )
    )
    utils::write.table(
      ctr_data,
      rnk_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    message("  ", ctr, ": ", nrow(ctr_data), " sites -> ", basename(rnk_file))
    rnk_files <- c(rnk_files, rnk_file)
  }

  message("=== Done ===")
  invisible(c(seqwindows_file, rnk_files))
}

#' Give the Flanking Sequence Column its Canonical Name
#'
#' The flanking sequence arrives as `SequenceWindow` from a prolfquapp DEA and
#' as `PTM_FlankingRegion` from some search engines' output. Both mean the same
#' thing.
#'
#' @param data Data frame that should carry a flanking sequence column.
#' @return `data` with the column named `SequenceWindow`.
#' @keywords internal
canonicalize_sequence_window <- function(data) {
  if ("SequenceWindow" %in% names(data)) {
    return(data)
  }
  if ("PTM_FlankingRegion" %in% names(data)) {
    message("  Renamed PTM_FlankingRegion -> SequenceWindow")
    return(dplyr::rename(data, SequenceWindow = "PTM_FlankingRegion"))
  }
  stop("No SequenceWindow or PTM_FlankingRegion column found", call. = FALSE)
}

#' Drop Sequence Windows a Motif Scan Cannot Use
#'
#' A window is usable only if it is a full, uninterrupted stretch of residues.
#' Windows padded with underscores because the site sits near a protein
#' terminus, and windows shorter than seven residues, carry too little context
#' for a motif match and are dropped rather than scanned.
#'
#' @param data Data frame with a `SequenceWindow` column.
#' @return `data` reduced to usable windows, upper-cased.
#' @keywords internal
filter_sequence_windows <- function(data) {
  data |>
    dplyr::filter(
      !is.na(.data$SequenceWindow),
      .data$SequenceWindow != "",
      !grepl("^_", .data$SequenceWindow),
      !grepl("_$", .data$SequenceWindow),
      nchar(.data$SequenceWindow) >= 7
    ) |>
    dplyr::mutate(SequenceWindow = toupper(.data$SequenceWindow))
}

#' Rank the Sites of One Contrast for Motif Enrichment
#'
#' Motif enrichment walks a single ranked list, so a sequence window may appear
#' only once. Where several sites share a window the most extreme statistic is
#' kept — the one that would move the enrichment score most — rather than an
#' average, which would let two opposing sites cancel into an uninformative
#' middle.
#'
#' @param data Data frame of one analysis, already filtered.
#' @param stat_column Column to rank on.
#' @param contrast_name Contrast to select.
#' @return A two-column table, `SequenceWindow` and `statistic.site`, ordered
#'   from the largest statistic down.
#' @keywords internal
rank_sites_for_mea <- function(data, stat_column, contrast_name) {
  data |>
    dplyr::filter(
      .data$contrast == contrast_name,
      !is.na(.data[[stat_column]])
    ) |>
    dplyr::mutate(statistic.site = .data[[stat_column]]) |>
    dplyr::select("SequenceWindow", "statistic.site") |>
    dplyr::group_by(.data$SequenceWindow) |>
    dplyr::slice(which.max(abs(.data$statistic.site))) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data$statistic.site))
}
