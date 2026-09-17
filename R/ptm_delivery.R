.ptm_delivery_tables <- function(statistics) {
  dpa_dpu <- statistics$get_dpa_dpu()
  dpa <- standardize_ptm_results(dpa_dpu$combined_site_prot, "dpa")
  dpu <- standardize_ptm_results(dpa_dpu$combined_test_diff, "dpu")
  cf <- standardize_ptm_results(statistics$get_cf()$results, "cf")
  dpa <- .ptm_delivery_missing(dpa)
  dpu <- .ptm_delivery_missing(dpu)
  cf <- .ptm_delivery_missing(cf)
  pair <- statistics$get_inputs()$get_pair()
  protein_long <- .ptm_abundance_long(pair$protein) |>
    dplyr::filter(!grepl("^rev_", .data$protein_Id)) |>
    canonicalize_uniprot_ids()
  site_long <- .ptm_abundance_long(pair$site)
  names(protein_long)[names(protein_long) == pair$protein$sample_key] <- "Name"
  names(site_long)[names(site_long) == pair$site$sample_key] <- "Name"
  protein_abund <- protein_long |>
    dplyr::select("Name", "protein_Id", "normalized_abundance") |>
    tidyr::pivot_wider(
      names_from = "Name",
      values_from = "normalized_abundance"
    )

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

  list(
    DPA = dpa,
    DPU = dpu,
    CF = cf,
    abundances_protein = protein_abund,
    abundances_site_dpa = site_abund_dpa,
    abundances_site_cf = site_abund_cf
  )
}

.ptm_delivery_missing <- function(data) {
  # Excel's blank string cells were read as NA in the former combined delivery.
  for (column in names(data)[vapply(data, is.character, logical(1))]) {
    data[[column]][!is.na(data[[column]]) & data[[column]] == ""] <- NA_character_
  }
  data
}

.ptm_abundance_long <- function(experiment) {
  data <- experiment$normalized_abundances
  obs <- experiment$obs
  sample_order <- obs[[experiment$sample_key]][order(obs[[experiment$configuration$file_name]])]
  data[order(match(data[[experiment$sample_key]], sample_order)), , drop = FALSE]
}
