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
#' contrast.
#'
#' @param phospho_dea_dir Path to the phospho DEA output directory; its
#'   `Results_WU_*/AnnData.h5ad` is read.
#' @param protein_dea_dir Path to the total-proteome DEA output directory; its
#'   `Results_WU_*/AnnData.h5ad` is read.
#' @param remove_contaminants Drop the sites and proteins the DEAs flag as
#'   contaminants (`CON`). Like the DEAs, the default keeps them.
#' @return A list with
#'   `combined_site_prot`, the DPA table;
#'   `combined_test_diff`, the moderated DPU table;
#'   `combined_test_diff_unmoderated`, the unmoderated DPU table;
#'   `n_unmoderated_untestable`, the number of paired rows whose raw degrees of
#'   freedom do not permit an unmoderated t-test;
#'   `match_rates`, sites tested and sites paired with a protein, per contrast.
#' @seealso [compute_cf_dea()] for the alternative that corrects before
#'   modelling rather than after.
#' @export
#' @examples
#' \dontrun{
#' res <- compute_dpa_dpu("DEA_20260814_WUphospho_vsn", "DEA_20260814_WUtotal_vsn")
#' res$match_rates
#' }
compute_dpa_dpu <- function(phospho_dea_dir, protein_dea_dir, remove_contaminants = FALSE) {
  .compute_dpa_dpu_from_pair(.read_dea_pair(phospho_dea_dir, protein_dea_dir, remove_contaminants))
}

.compute_dpa_dpu_from_pair <- function(pair) {
  # A site without an FDR was not tested and can neither be reported nor
  # corrected. Dropping it here keeps DPA and DPU on the same set of sites.
  site <- dplyr::filter(pair$site$differential_results, !is.na(.data$FDR))
  protein <- canonicalize_uniprot_ids(pair$protein$differential_results)
  # description and protein_length are joined on as well, so that a pair is
  # only formed within one protein annotation.
  join_column <- c("protein_Id", "contrast", "description", "protein_length")
  combined_site_prot <- dplyr::left_join(site, protein, by = join_column, suffix = c(".site", ".protein"))
  match_rates <- combined_site_prot |>
    dplyr::group_by(.data$contrast) |>
    dplyr::summarize(
      total_sites = dplyr::n(),
      matched_sites = sum(!is.na(.data$diff.protein)),
      match_rate = round(.data$matched_sites / .data$total_sites * 100, 1)
    )
  unmoderated <- test_diff(site, protein, join_column = join_column, variant = "unmoderated")
  testable <- function(df) is.finite(df) & df > 0
  list(
    combined_site_prot = combined_site_prot,
    combined_test_diff = test_diff(site, protein, join_column = join_column, variant = "moderated"),
    combined_test_diff_unmoderated = unmoderated,
    n_unmoderated_untestable = sum(
      unmoderated$measured_In == "both" &
        !(testable(unmoderated$df.unmoderated.site) & testable(unmoderated$df.unmoderated.protein))
    ),
    match_rates = match_rates
  )
}
