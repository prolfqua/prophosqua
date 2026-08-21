#' Standardize a PTM Result Table
#'
#' DPA, DPU and CF arrive with different column names for the same three
#' quantities, because each is produced by a different route. This brings all
#' three onto the one vocabulary every downstream report reads: `diff.site` for
#' the log2 fold change, `FDR.site` for its adjusted p-value and
#' `statistic.site` for its test statistic.
#'
#' Columns absent from `data` are dropped rather than raising: a DEA run that
#' did not record `SequenceWindow`, say, still produces a usable table, only
#' without the motif reports. That silence is the reason a column can be present
#' in `Result_DPU.xlsx` and missing from every report — check the selection
#' below before suspecting the analysis.
#'
#' @param data Data frame of results for one analysis type.
#' @param analysis_type One of `"dpa"`, `"dpu"` or `"cf"`, case-insensitive.
#' @return `data` reduced to the standard columns, renamed.
#' @export
#' @examples
#' dpu <- data.frame(
#'   protein_Id = "P1", site = "P1~S1", contrast = "a_vs_b",
#'   gene_name.site = "GENE", diff_diff = 1.5, FDR_I = 0.01, tstatistic_I = 4
#' )
#' standardize_ptm_results(dpu, "dpu")
standardize_ptm_results <- function(data, analysis_type) {
  analysis_type <- tolower(analysis_type)

  # DPA and CF export diff.site / FDR.site / statistic.site directly. DPU comes
  # from test_diff(), which names them after the interaction it estimates.
  site_annotation <- c(
    "protein_Id",
    "site",
    "contrast",
    "posInProtein",
    "modAA",
    "SequenceWindow",
    "protein_length"
  )
  configs <- list(
    dpa = list(
      rename = c(gene_name = "gene_name.site"),
      direct_cols = c(
        site_annotation,
        "diff.site",
        "FDR.site",
        "statistic.site",
        "diff.protein",
        "FDR.protein",
        "statistic.protein",
        "estimate_type.site",
        "estimate_type.protein"
      )
    ),
    dpu = list(
      rename = c(
        gene_name = "gene_name.site",
        diff.site = "diff_diff",
        FDR.site = "FDR_I",
        statistic.site = "tstatistic_I"
      ),
      direct_cols = c(
        site_annotation,
        "estimate_type.site",
        "estimate_type.protein"
      )
    ),
    cf = list(
      # CF has a single estimate per site, exported unsuffixed. Its gene_name
      # arrives already named and so sits inside the annotation block, where
      # DPA and DPU carry theirs at the end because it comes from the rename.
      # The order is pinned by a test: downstream readers select by name, but
      # the sheet is also read by people.
      rename = c(estimate_type.site = "estimate_type"),
      direct_cols = c(
        setdiff(site_annotation, "protein_length"),
        "gene_name",
        "protein_length",
        "diff.site",
        "FDR.site",
        "statistic.site"
      )
    )
  )

  if (!analysis_type %in% names(configs)) {
    stop(
      "Unknown analysis_type: ",
      analysis_type,
      ". Must be one of: ",
      paste(names(configs), collapse = ", "),
      call. = FALSE
    )
  }

  config <- configs[[analysis_type]]
  select_spec <- c(
    stats::setNames(config$direct_cols, config$direct_cols),
    config$rename
  )
  existing <- select_spec[select_spec %in% names(data)]

  dplyr::select(data, !!!existing)
}

#' Combine the Three PTM Analyses into One Workbook
#'
#' Reads the DPA, DPU and CF result workbooks, standardizes their columns, and
#' writes them beside the normalized abundances as a single multi-sheet workbook
#' and the same content as an RDS. Every downstream report of the pipeline reads
#' this workbook rather than the three separate ones.
#'
#' Three abundance sheets are produced: protein abundances, site abundances, and
#' site abundances corrected by their protein — the last matching what the
#' CorrectFirst analysis models, so a reader can check a CF result against the
#' values behind it.
#'
#' @param dpa_xlsx,dpu_xlsx,cf_xlsx The three result workbooks.
#' @param protein_parquet,site_parquet Normalized abundances of the two DEA runs.
#' @param output_xlsx,output_rds Where to write the combined result.
#' @return Invisibly, the list that was written.
#' @export
#' @examples
#' # Needs the outputs of a full pipeline run.
#' \dontrun{
#' combine_ptm_results(
#'   "PTM_DPA/Result_DPA.xlsx",
#'   "PTM_DPU/Result_DPU.xlsx",
#'   "PTM_CF_DPU/CorrectFirst_PTM_usage_results.xlsx",
#'   "DEA_total/Results_WU_x/lfqdata_normalized.parquet",
#'   "DEA_phospho/Results_WU_y/lfqdata_normalized.parquet",
#'   "PTM_results.xlsx", "PTM_results.rds"
#' )
#' }
combine_ptm_results <- function(dpa_xlsx, dpu_xlsx, cf_xlsx, protein_parquet, site_parquet, output_xlsx, output_rds) {
  message("Loading DPA results from: ", dpa_xlsx)
  dpa <- standardize_ptm_results(
    readxl::read_xlsx(dpa_xlsx, sheet = "combinedSiteProteinData"),
    "dpa"
  )

  message("Loading DPU results from: ", dpu_xlsx)
  dpu <- standardize_ptm_results(
    readxl::read_xlsx(dpu_xlsx, sheet = "combinedStats"),
    "dpu"
  )

  message("Loading CF results from: ", cf_xlsx)
  cf <- standardize_ptm_results(
    readxl::read_xlsx(cf_xlsx, sheet = "results"),
    "cf"
  )

  message("  DPA: ", nrow(dpa), " rows")
  message("  DPU: ", nrow(dpu), " rows")
  message("  CF:  ", nrow(cf), " rows")

  message("Loading protein abundances from: ", protein_parquet)
  protein_long <- read_normalized_abundances(protein_parquet) |>
    dplyr::filter(!grepl("^rev_", .data$protein_Id)) |>
    canonicalize_uniprot_ids()

  protein_abund <- protein_long |>
    dplyr::select("Name", "protein_Id", "normalized_abundance") |>
    tidyr::pivot_wider(
      names_from = "Name",
      values_from = "normalized_abundance"
    )

  message("Loading site abundances from: ", site_parquet)
  site_long <- read_normalized_abundances(site_parquet)
  site_col <- site_column(site_long)

  site_abund_dpa <- site_long |>
    dplyr::select(
      "Name",
      site = tidyselect::all_of(site_col),
      "protein_Id",
      "normalized_abundance"
    ) |>
    tidyr::pivot_wider(
      names_from = "Name",
      values_from = "normalized_abundance"
    )

  # The CF sheet holds what CorrectFirst models: the site value minus its
  # protein value in the same sample.
  message("Computing corrected abundances for CF...")
  site_abund_cf <- site_long |>
    dplyr::select(
      "Name",
      site = tidyselect::all_of(site_col),
      "protein_Id",
      site_abund = "normalized_abundance"
    ) |>
    dplyr::inner_join(
      protein_long |>
        dplyr::select("Name", "protein_Id", protein_abund = "normalized_abundance"),
      by = c("Name", "protein_Id")
    ) |>
    dplyr::mutate(
      corrected_abundance = .data$site_abund - .data$protein_abund
    ) |>
    dplyr::select("Name", "site", "corrected_abundance") |>
    tidyr::pivot_wider(
      names_from = "Name",
      values_from = "corrected_abundance"
    )

  message("  Protein abundances: ", nrow(protein_abund), " proteins x ", ncol(protein_abund) - 1, " samples")
  message("  Site abundances (DPA): ", nrow(site_abund_dpa), " sites x ", ncol(site_abund_dpa) - 2, " samples")
  message("  Site abundances (CF): ", nrow(site_abund_cf), " sites x ", ncol(site_abund_cf) - 1, " samples")

  result_list <- list(
    DPA = dpa,
    DPU = dpu,
    CF = cf,
    abundances_protein = protein_abund,
    abundances_site_dpa = site_abund_dpa,
    abundances_site_cf = site_abund_cf
  )

  message("Writing Excel to: ", output_xlsx)
  writexl::write_xlsx(result_list, output_xlsx)

  message("Writing RDS to: ", output_rds)
  saveRDS(result_list, output_rds)

  invisible(result_list)
}

#' Read a DEA Parquet with its Sample Column Canonicalized
#'
#' The sample column of a prolfquapp parquet is named by the run's YAML, which
#' sits beside the parquet. Reading the two together is the only way to get a
#' `Name` column that means the same thing across runs.
#'
#' @param parquet Path to `lfqdata_normalized.parquet`.
#' @return The long table with its sample column renamed to `Name`.
#' @keywords internal
read_normalized_abundances <- function(parquet) {
  arrow::read_parquet(parquet) |>
    canonicalize_dea_sample_column(
      file.path(dirname(parquet), "lfqdata.yaml")
    )
}
