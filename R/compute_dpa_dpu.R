#' Compute Differential PTM Abundance and Differential PTM Usage
#'
#' Pairs the site-level results of a phospho DEA run with the protein-level
#' results of a total-proteome run and derives the two integrated views the PTM
#' pipeline reports on:
#'
#' * **DPA**, differential PTM abundance: the site result with its protein
#'   counterpart joined alongside, suffixed `.site` and `.protein`.
#' * **DPU**, differential PTM usage: the protein-normalized site result
#'   computed by [test_diff()], where the effect size is the difference of the
#'   two log2 fold changes and its standard error the root of the summed
#'   squares.
#'
#' Sites whose protein was not quantified keep a DPA row with empty `.protein`
#' columns; they cannot carry a DPU value. `match_rates` reports that split per
#' contrast so a report can state it without recomputing the join.
#'
#' Nothing is re-quantified here. Both DEA runs must already have happened; this
#' function only reads their output.
#'
#' @param phospho_dea_dir Path to the phospho DEA output directory.
#' @param protein_dea_dir Path to the total-proteome DEA output directory.
#' @return A list with
#'   `combined_site_prot`, the DPA table;
#'   `combined_test_diff`, the DPU table;
#'   `match_rates`, sites tested and sites paired with a protein, per contrast.
#' @seealso [compute_cf_dea()] for the alternative that corrects before
#'   modelling rather than after.
#' @export
#' @examples
#' # Needs two prolfquapp DEA output directories; see
#' # tests/testthat/helper-dea_fixture.R for the synthetic pair the tests use.
#' \dontrun{
#' res <- compute_dpa_dpu(
#'   phospho_dea_dir = "DEA_20260814_WUphospho_vsn",
#'   protein_dea_dir = "DEA_20260814_WUtotal_vsn"
#' )
#' res$match_rates
#' }
compute_dpa_dpu <- function(phospho_dea_dir, protein_dea_dir) {
  tot_file <- get_dea_xlsx(protein_dea_dir)
  ptm_file <- get_dea_xlsx(phospho_dea_dir)

  stopifnot(file.exists(tot_file), file.exists(ptm_file))

  required_cols <- c("protein_Id", "protein_length", "contrast")

  tot_res <- load_and_preprocess_data(tot_file, required_cols)
  tot_res <- filter_contaminants(tot_res)
  tot_res <- canonicalize_uniprot_ids(tot_res)

  phospho_res <- load_and_preprocess_data(ptm_file, required_cols)
  phospho_res <- filter_contaminants(phospho_res)

  # The site annotation -- modAA, posInProtein, SequenceWindow -- arrives with
  # the DEA result: every PTM reader attaches it to the analysis rows, so it is
  # already here and joining it again would suffix both copies .x and .y.
  required_site_cols <- c("posInProtein", "modAA", "SequenceWindow")
  missing_site_cols <- setdiff(required_site_cols, colnames(phospho_res))
  if (length(missing_site_cols) > 0) {
    stop(
      "The phospho DEA result carries no site annotation (missing: ",
      paste(missing_site_cols, collapse = ", "),
      "). It was produced by a reader older than the one that attaches it; ",
      "rerun the DEA with the current prolfquappPTMreaders.",
      call. = FALSE
    )
  }

  # The DPA join carries the protein result of the same contrast alongside the
  # site result. description and protein_length are joined on as well so that
  # the pair is only formed within one protein annotation.
  join_column <- c(
    "protein_Id" = "protein_Id",
    "contrast",
    "description",
    "protein_length"
  )

  # A site without an FDR was not tested and can neither be reported nor
  # corrected. Dropping it here keeps DPA and DPU on the same set of sites.
  phospho_res <- phospho_res |> dplyr::filter(!is.na(.data$FDR))

  combined_site_prot <- dplyr::left_join(
    phospho_res,
    tot_res,
    by = join_column,
    suffix = c(".site", ".protein")
  )

  match_rates <- combined_site_prot |>
    dplyr::group_by(.data$contrast) |>
    dplyr::summarize(
      total_sites = dplyr::n(),
      matched_sites = sum(!is.na(.data$diff.protein)),
      match_rate = round(.data$matched_sites / .data$total_sites * 100, 1)
    )

  combined_test_diff <- test_diff(phospho_res, tot_res, join_column = join_column)

  # Downstream reports and combine_ptm_results() key on `site`; a newer DEA
  # writes the column as protein_Id_site.
  if (
    !"site" %in% colnames(combined_test_diff) &&
      "protein_Id_site" %in% colnames(combined_test_diff)
  ) {
    combined_test_diff <- combined_test_diff |>
      dplyr::mutate(site = .data$protein_Id_site)
  }

  list(
    combined_site_prot = combined_site_prot,
    combined_test_diff = combined_test_diff,
    match_rates = match_rates
  )
}

#' Name of the Site Identifier Column
#'
#' A prolfquapp DEA run names the site identifier either `site` or
#' `protein_Id_site`, depending on its version. Every caller that joins site
#' annotation onto a result table has to ask which one it got.
#'
#' @param x Data frame from a DEA result.
#' @return The name of the site column.
#' @keywords internal
site_column <- function(x) {
  candidates <- intersect(c("site", "protein_Id_site"), names(x))
  if (length(candidates) == 0) {
    stop(
      "no site identifier column found: expected `site` or `protein_Id_site`, got ",
      paste(utils::head(names(x), 10), collapse = ", "),
      call. = FALSE
    )
  }
  candidates[[1]]
}
