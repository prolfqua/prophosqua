.ptm_delivery_tables <- function(statistics) {
  dpa_dpu <- statistics$get_dpa_dpu()
  cf <- statistics$get_cf()
  pair <- statistics$get_inputs()$get_pair()
  protein_abund <- .ptm_abundance_long(pair$protein) |>
    dplyr::filter(!grepl("^rev_", .data$protein_Id)) |>
    canonicalize_uniprot_ids() |>
    dplyr::select(Name = tidyselect::all_of(pair$protein$sample_key), "protein_Id", "normalized_abundance") |>
    tidyr::pivot_wider(names_from = "Name", values_from = "normalized_abundance")
  site_abund_dpa <- .ptm_abundance_long(pair$site) |>
    dplyr::select(Name = tidyselect::all_of(pair$site$sample_key), "site", "protein_Id", "normalized_abundance") |>
    tidyr::pivot_wider(names_from = "Name", values_from = "normalized_abundance")
  # The CF sheet holds the values CorrectFirst modelled.
  samples <- setdiff(names(site_abund_dpa), c("site", "protein_Id"))
  list(
    DPA = standardize_ptm_results(dpa_dpu$combined_site_prot, "dpa"),
    DPU = standardize_ptm_results(dpa_dpu$combined_test_diff, "dpu"),
    CF = standardize_ptm_results(cf$results, "cf"),
    abundances_protein = protein_abund,
    abundances_site_dpa = site_abund_dpa,
    abundances_site_cf = dplyr::select(cf$wide_data, "site", tidyselect::all_of(samples))
  )
}

.ptm_abundance_long <- function(experiment) {
  data <- experiment$normalized_abundances
  obs <- experiment$obs
  sample_order <- obs[[experiment$sample_key]][order(obs[[experiment$configuration$file_name]])]
  data[order(match(data[[experiment$sample_key]], sample_order)), , drop = FALSE]
}

# DPA, DPU and CF name the same three quantities differently, because each is
# produced by a different route. This brings them onto the names every reader
# uses: diff.site, FDR.site and statistic.site. A column the analysis did not
# produce is dropped, so a column missing here is missing from every report.
# The column order is the one people read the sheets in.
standardize_ptm_results <- function(data, analysis_type) {
  site_annotation <- c("protein_Id", "site", "contrast", "posInProtein", "modAA", "SequenceWindow", "protein_length")
  columns <- switch(
    tolower(analysis_type),
    dpa = c(
      site_annotation,
      "diff.site",
      "FDR.site",
      "statistic.site",
      "diff.protein",
      "FDR.protein",
      "statistic.protein",
      "estimate_type.site",
      "estimate_type.protein",
      gene_name = "gene_name.site"
    ),
    dpu = c(
      site_annotation,
      "estimate_type.site",
      "estimate_type.protein",
      gene_name = "gene_name.site",
      diff.site = "diff_diff",
      FDR.site = "FDR_I",
      statistic.site = "tstatistic_I"
    ),
    cf = c(
      setdiff(site_annotation, "protein_length"),
      "gene_name",
      "protein_length",
      "diff.site",
      "FDR.site",
      "statistic.site",
      estimate_type.site = "estimate_type"
    ),
    stop("Unknown analysis_type: ", analysis_type, ". Must be one of: dpa, dpu, cf", call. = FALSE)
  )
  unnamed <- !nzchar(names(columns))
  names(columns)[unnamed] <- columns[unnamed]
  dplyr::select(data, !!!columns[columns %in% names(data)])
}
