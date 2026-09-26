#' Compute PTM Usage the CorrectFirst Way
#'
#' Corrects phosphosite abundances by their protein abundance **before**
#' modelling, then fits one linear model per site on the corrected values and
#' evaluates the contrasts on those models. This is the mirror image of
#' [test_diff()], which models site and protein separately and subtracts the two
#' fold changes afterwards.
#'
#' The correction is a subtraction on the log2 scale, so it is a ratio on the
#' linear scale: a site keeps only the signal its protein does not explain.
#' The sample's median protein abundance is added back, so the corrected value
#' stays on the log2 intensity scale. Sites whose protein was not quantified in
#' the same sample drop out at the join, and a site left without a corrected
#' value in any sample is absent from the result.
#'
#' Two variants are returned under `variants`. `correct_first_protein_imputed`
#' corrects the observed site values with the protein DEA's `imputedData`
#' layer and is fitted like CF: it has a `results` table and its corrected
#' `abundances` (samples x sites). `correct_first_site_protein_imputed` also
#' fills the site values from the site's own lm fit (sites the site DEA
#' refitted at the LOD keep their gaps) and has `abundances` only: its filled
#' values are the site model's predictions, and a fit on them would count them
#' as observations and reuse the site's data. The per-sample median added back
#' is always that of the observed protein values. Every results row counts the
#' imputed protein values its fit used (`n_protein_imputed`).
#'
#' The model formula and the contrasts are the ones the DEAs recorded; the site
#' and protein DEAs must agree on the formula.
#'
#' The model is prolfqua's `lm_impute` facade, as in the prolfquapp DEA behind
#' DPA: a site whose fit fails or has fewer than two residual degrees of
#' freedom, typically one observed in a single group, is refitted with its
#' missing values at the limit of detection, the lower quartile of corrected
#' values seen once in a group, and a variance borrowed from the fully fitted
#' sites. Those rows carry `estimate_type == "lod_imputed"`, and
#' `model_counts` records the split by `estimate_type`.
#'
#' @param phospho_dea_dir Path to the phospho DEA output directory; its
#'   `Results_WU_*/AnnData.h5ad` is read.
#' @param protein_dea_dir Path to the total-proteome DEA output directory; its
#'   `Results_WU_*/AnnData.h5ad` is read.
#' @param remove_contaminants Drop the sites and proteins the DEAs flag as
#'   contaminants (`CON`). Like the DEAs, the default keeps them.
#' @return A list carrying the result table and everything a report needs to
#'   describe it without refitting: `results`, `ptm_data`,
#'   `contrasts`, `wide_data`, `wide_annotation`, `model_counts`, the
#'   measurement and model counts the prose quotes, and `variants`.
#' @seealso [compute_dpa_dpu()] for the correct-last alternative.
#' @export
#' @examples
#' \dontrun{
#' res <- compute_cf_dea("DEA_20260814_WUphospho_vsn", "DEA_20260814_WUtotal_vsn")
#' res$model_counts
#' }
compute_cf_dea <- function(phospho_dea_dir, protein_dea_dir, remove_contaminants = FALSE) {
  .compute_cf_dea_from_pair(.read_dea_pair(phospho_dea_dir, protein_dea_dir, remove_contaminants))
}

.compute_cf_dea_from_pair <- function(pair) {
  modelstr <- .cf_model_string(pair)
  protein_imputed <- .require_imputed_abundances(pair$protein, "total-proteome")
  correct_first <- function(protein_data) {
    .cf_fit(
      .cf_corrected(pair$site$normalized_abundances, protein_data, pair),
      pair$site$site_info,
      pair$site$contrasts,
      modelstr
    )
  }
  result <- correct_first(pair$protein$normalized_abundances)
  protein_imputed_cf <- correct_first(protein_imputed)
  # B's filled site values are the site model's own predictions; a fit would
  # count them as observations and reuse the site's data, so B is a layer only.
  site_protein_imputed <- .cf_corrected(.lm_filled_sites(pair$site), protein_imputed, pair)
  samples <- as.character(pair$site$obs[[pair$site$sample_key]])
  sites <- .cf_sites(
    pair$site$var,
    list(result$wide_data, protein_imputed_cf$wide_data, site_protein_imputed$wide_data)
  )
  result$variants <- list(
    correct_first_protein_imputed = list(
      results = protein_imputed_cf$results,
      abundances = .cf_abundances(protein_imputed_cf$wide_data, samples, sites)
    ),
    correct_first_site_protein_imputed = list(
      abundances = .cf_abundances(site_protein_imputed$wide_data, samples, sites)
    )
  )
  result
}

# The sites CorrectFirst corrects in any variant, in the site experiment's
# order: every site whose protein the total-proteome DEA quantified. They are
# the feature axis of the enriched_CF modality.
.cf_sites <- function(var, wides) {
  corrected <- unlist(lapply(wides, function(wide) as.character(wide$site)))
  var$site[var$site %in% corrected]
}

# CF fits the design the DEAs used; the site and protein DEAs must agree on it.
.cf_model_string <- function(pair) {
  terms <- lapply(
    list(site = pair$site$formula, protein = pair$protein$formula),
    function(formula) {
      if (is.null(formula) || !grepl("~", formula, fixed = TRUE)) {
        stop("A DEA artifact records no model formula; found '", format(formula), "'.", call. = FALSE)
      }
      attr(stats::terms(stats::as.formula(sub("^[^~]*~", "~", formula))), "term.labels")
    }
  )
  if (!setequal(terms$site, terms$protein)) {
    stop(
      "Site and protein DEAs used different models: '~ ",
      paste(terms$site, collapse = " + "),
      "' and '~ ",
      paste(terms$protein, collapse = " + "),
      "'.",
      call. = FALSE
    )
  }
  paste("~", paste(terms$site, collapse = " + "))
}

.require_imputed_abundances <- function(experiment, side) {
  if (is.null(experiment$imputed_abundances)) {
    stop(
      "The ",
      side,
      " AnnData has no imputedData layer; rerun that DEA with ",
      "prolfquapp >= 2.10.5 and the lm_impute model.",
      call. = FALSE
    )
  }
  experiment$imputed_abundances
}

# Site values filled by the site's own lm only: a site whose missing cells
# prolfqua predicted from its observed-data fit (route "fitted") takes those
# values; a site refitted at the LOD keeps its gaps, which the CF model's own
# refit then handles.
.lm_filled_sites <- function(site) {
  keys <- c(site$sample_key, site$feature_keys)
  fitted <- dplyr::filter(site$imputation, .data$route == "fitted")
  filled <- .require_imputed_abundances(site, "site") |>
    dplyr::semi_join(fitted, by = site$feature_keys) |>
    dplyr::select(tidyselect::all_of(keys), fitted_abundance = "normalized_abundance")
  site$normalized_abundances |>
    dplyr::left_join(filled, by = keys) |>
    dplyr::mutate(normalized_abundance = dplyr::coalesce(.data$fitted_abundance, .data$normalized_abundance)) |>
    dplyr::select(-"fitted_abundance")
}

# Corrected abundances, samples x sites, NA for a site without a corrected
# value.
.cf_abundances <- function(wide, samples, sites) {
  abundances <- matrix(NA_real_, nrow = length(samples), ncol = length(sites), dimnames = list(samples, sites))
  corrected <- as.character(wide$site)
  shared <- intersect(sites, corrected)
  abundances[, shared] <- t(as.matrix(wide[match(shared, corrected), samples, drop = FALSE]))
  abundances
}

# Per site, how many of the corrected values rest on a protein value that was
# not observed.
.cf_protein_imputed_counts <- function(corrected, protein_observed, site_sample_col, protein_sample_col) {
  protein_values <- protein_observed |>
    dplyr::select(tidyselect::all_of(protein_sample_col), "protein_Id", protein_observed = "normalized_abundance")
  corrected |>
    dplyr::filter(!is.na(.data$ptm_usage)) |>
    dplyr::left_join(protein_values, by = c(stats::setNames(protein_sample_col, site_sample_col), "protein_Id")) |>
    dplyr::summarize(n_protein_imputed = sum(is.na(.data$protein_observed)), .by = "site")
}

# The one CorrectFirst correction: site minus its protein in the same sample,
# plus the sample's median observed protein abundance.
.cf_corrected <- function(site_data, protein_data, pair) {
  site_sample_col <- pair$site$sample_key
  protein_sample_col <- pair$protein$sample_key
  accession_keyed <- function(data) {
    canonicalize_uniprot_ids(dplyr::filter(data, !grepl("^rev_", .data$protein_Id)))
  }
  tot_d <- accession_keyed(protein_data) |>
    dplyr::select(tidyselect::all_of(protein_sample_col), "protein_Id", "normalized_abundance")
  protein_observed <- accession_keyed(pair$protein$normalized_abundances)
  protein_sample_median <- protein_observed |>
    dplyr::summarize(
      protein_median = stats::median(.data$normalized_abundance, na.rm = TRUE),
      .by = tidyselect::all_of(protein_sample_col)
    )

  ptm_data <- prolfqua::LFQData$new(site_data, pair$site$configuration$clone(deep = TRUE))
  n_site_measurements <- nrow(ptm_data$data_long())
  sample_join <- stats::setNames(protein_sample_col, site_sample_col)
  ptm_data$set_data(dplyr::inner_join(
    ptm_data$data_long(),
    tot_d,
    by = c(sample_join, protein_Id = "protein_Id"),
    suffix = c(".site", ".total")
  ))
  # Adding back the sample's median protein abundance keeps the corrected value
  # on the log2 intensity scale, and corrects each site by how far its protein
  # lies from that sample's median: a sample-wide offset between the separately
  # normalised proteome and site data does not enter the correction.
  ptm_data$set_data(
    ptm_data$data_long() |>
      dplyr::left_join(protein_sample_median, by = sample_join) |>
      dplyr::mutate(
        ptm_usage = .data$normalized_abundance.site - .data$normalized_abundance.total + .data$protein_median
      ) |>
      dplyr::select(-"protein_median")
  )
  n_merged_measurements <- nrow(ptm_data$data_long())
  # A site and its protein can each be measured, yet never in the same sample;
  # such a site has no corrected value, and imputation would invent one.
  ptm_data$set_data(
    dplyr::filter(ptm_data$data_long(), any(!is.na(.data$ptm_usage)), .by = dplyr::all_of(ptm_data$subject_id()))
  )
  ptm_data$get_config()$set_response("ptm_usage")
  c(
    list(ptm_data = ptm_data),
    .cf_wide(ptm_data),
    list(
      n_protein_imputed = .cf_protein_imputed_counts(
        ptm_data$data_long(),
        protein_observed,
        site_sample_col,
        protein_sample_col
      ),
      n_protein_measurements = nrow(tot_d),
      n_site_measurements = n_site_measurements,
      n_merged_measurements = n_merged_measurements
    )
  )
}

# The corrected values as one row per site, samples as columns, and the sample
# annotation beside them.
.cf_wide <- function(ptm_data) {
  wide <- ptm_data$data_wide()
  list(
    wide_data = prolfqua::separate_hierarchy(wide$data, ptm_data$get_config()),
    wide_annotation = wide$annotation
  )
}

.cf_fit <- function(corrected, site_info, contrasts, modelstr) {
  ptm_data <- corrected$ptm_data
  # The default LOD, the median of values seen once in a group, falls inside
  # the corrected-value distribution, because a faint site on a faint protein
  # has an ordinary ratio. The lower quartile puts it about where the default
  # LOD sits in the site intensities.
  lod <- prolfqua::MissingHelpers$new(ptm_data$data_long(), ptm_data$get_config(), prob = 0.25)$get_lod()
  facade <- prolfqua::ContrastsLMImputeFacade$new(ptm_data, modelstr, contrasts, lod = lod, weights = NULL)
  ctr_df <- facade$get_contrasts()
  # DPA, DPU and CF are read interchangeably downstream, so all three name the
  # effect size, its adjusted p-value and its statistic the same way.
  results <- ctr_df |>
    dplyr::left_join(dplyr::select(site_info, -"protein_Id"), by = "site") |>
    dplyr::rename(diff.site = "diff", FDR.site = "FDR", statistic.site = "statistic") |>
    dplyr::left_join(corrected$n_protein_imputed, by = "site")
  # An empty result would leave every downstream report empty but succeeding.
  stopifnot("no site-contrast pair could be estimated" = nrow(results) > 0)
  list(
    results = results,
    ptm_data = ptm_data,
    contrasts = contrasts,
    wide_data = corrected$wide_data,
    wide_annotation = corrected$wide_annotation,
    model_counts = dplyr::summarize(dplyr::group_by(results, .data$estimate_type), Site_contrast_pairs = dplyr::n()),
    n_before = nrow(results),
    n_protein_measurements = corrected$n_protein_measurements,
    n_site_measurements = corrected$n_site_measurements,
    n_merged_measurements = corrected$n_merged_measurements,
    n_models = nrow(facade$model$model_df),
    n_site_contrast = nrow(ctr_df)
  )
}
