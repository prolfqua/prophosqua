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
#' Sites whose protein was not quantified in the same sample drop out at the
#' join and are absent from the result.
#'
#' Site-contrast pairs whose estimate would rest on an imputed value in one of
#' the two groups are dropped from `results`, because a fold change measured
#' against an imputed level is not interpretable. `model_counts` and `n_before`
#' record that split so a report can state what was dropped.
#'
#' @param phospho_dea_dir Path to the phospho DEA output directory.
#' @param protein_dea_dir Path to the total-proteome DEA output directory.
#' @param annot_file Sample annotation defining the groups and the contrasts.
#'   Two layouts are accepted, see the contrast derivation below.
#' @return A list carrying the result table and everything a report needs to
#'   describe it without refitting: `results`, `ptm_data`, `ctr`, `annot`,
#'   `contrasts`, `wide_data`, `wide_annotation`, `model_counts`, and the
#'   measurement and model counts the prose quotes.
#' @seealso [compute_dpa_dpu()] for the correct-last alternative.
#' @export
#' @examples
#' # Needs two prolfquapp DEA output directories and an annotation; see
#' # tests/testthat/helper-dea_fixture.R for the synthetic set the tests use.
#' \dontrun{
#' res <- compute_cf_dea(
#'   phospho_dea_dir = "DEA_20260814_WUphospho_vsn",
#'   protein_dea_dir = "DEA_20260814_WUtotal_vsn",
#'   annot_file = "phospho_dataset.tsv"
#' )
#' res$model_counts
#' }
compute_cf_dea <- function(phospho_dea_dir, protein_dea_dir, annot_file) {
  annot <- readr::read_tsv(annot_file, show_col_types = FALSE)

  # --- the correction baseline: protein abundances -------------------------
  ldata <- arrow::read_parquet(get_dea_parquet(protein_dea_dir))
  ldata <- ldata |>
    dplyr::filter(!grepl("^rev_", .data$protein_Id)) |>
    canonicalize_uniprot_ids()
  protein_sample_col <- get_dea_sample_name_column(protein_dea_dir)

  tot_d <- ldata |>
    dplyr::select(
      tidyselect::all_of(protein_sample_col),
      "protein_Id",
      "normalized_abundance"
    )
  n_protein_measurements <- nrow(tot_d)

  # --- the response: site abundances --------------------------------------
  ldata <- arrow::read_parquet(get_dea_parquet(phospho_dea_dir))
  config <- prolfqua::list_to_AnalysisConfiguration(
    yaml::read_yaml(get_dea_yaml(phospho_dea_dir))
  )
  ptm_data <- prolfqua::LFQData$new(ldata, config)
  ptm_sample_col <- get_dea_sample_name_column(phospho_dea_dir)
  n_site_measurements <- nrow(ptm_data$data_long())

  # --- correct: site minus its protein, in the same sample ----------------
  sample_join <- stats::setNames(protein_sample_col, ptm_sample_col)
  join_columns <- c(sample_join, protein_Id = "protein_Id")

  ptm_data$set_data(dplyr::inner_join(
    ptm_data$data_long(),
    tot_d,
    by = join_columns,
    suffix = c(".site", ".total")
  ))
  ptm_data$set_data(
    ptm_data$data_long() |>
      dplyr::mutate(
        ptm_usage = .data$normalized_abundance.site - .data$normalized_abundance.total
      )
  )
  n_merged_measurements <- nrow(ptm_data$data_long())

  ptm_data$get_config()$set_response("ptm_usage")
  wide <- ptm_data$data_wide()
  wide_data <- prolfqua::separate_hierarchy(wide$data, ptm_data$get_config())

  # Correcting centres the values near zero, which puts half of them below the
  # axis of an intensity plot. The shift is cosmetic and constant across all
  # samples and sites, so every group difference the model estimates is
  # unchanged by it; it is applied before modelling only so that the
  # distribution and PCA a report draws from `ptm_data` show the same values
  # the model saw.
  ptm_data$set_data(
    ptm_data$data_long() |> dplyr::mutate(ptm_usage = .data$ptm_usage + 20)
  )

  # --- model and contrasts ------------------------------------------------
  strategy_lm <- prolfqua::strategy_lm("ptm_usage ~ G_")
  models <- prolfqua::build_model(data = ptm_data, model_strategy = strategy_lm)
  n_models <- nrow(models$model_df)

  contrasts <- derive_contrasts(annot, basename(annot_file))

  ctr <- prolfqua::Contrasts$new(models, contrasts)
  ctr <- prolfqua::ContrastsModerated$new(ctr)
  ctr_df <- ctr$get_contrasts()
  n_site_contrast <- nrow(ctr_df)

  # --- annotate and export shape -----------------------------------------
  site_info <- get_dea_ptm_site_info(phospho_dea_dir) |>
    dplyr::select(-"protein_Id")
  contrast_site_col <- site_column(ctr_df)
  ctr_df <- dplyr::left_join(
    ctr_df,
    site_info,
    by = stats::setNames("site", contrast_site_col)
  )

  # DPA, DPU and CF are read interchangeably downstream, so all three name the
  # effect size, its adjusted p-value and its statistic the same way.
  ctr_df <- ctr_df |>
    dplyr::rename(
      diff.site = "diff",
      FDR.site = "FDR",
      statistic.site = "statistic"
    ) |>
    dplyr::mutate(
      modelName = ifelse(
        .data$modelName %in% c("WaldTest", "WaldTest_moderated"),
        "Linear_Model_Moderated",
        "MissingInOneCondition"
      )
    )

  model_counts <- ctr_df |>
    dplyr::group_by(.data$modelName) |>
    dplyr::summarize(Site_contrast_pairs = dplyr::n())
  n_before <- nrow(ctr_df)

  results <- ctr_df |>
    dplyr::filter(.data$modelName != "MissingInOneCondition")

  # The filter is allowed to drop rows; what must not happen is dropping
  # everything, which would leave every downstream report empty but succeeding.
  stopifnot(
    "no site-contrast pair survived the imputation filter" = nrow(results) > 0
  )

  list(
    results = results,
    ptm_data = ptm_data,
    ctr = ctr,
    annot = annot,
    contrasts = contrasts,
    wide_data = wide_data,
    wide_annotation = wide$annotation,
    model_counts = model_counts,
    n_before = n_before,
    n_protein_measurements = n_protein_measurements,
    n_site_measurements = n_site_measurements,
    n_merged_measurements = n_merged_measurements,
    n_models = n_models,
    n_site_contrast = n_site_contrast
  )
}

#' Derive Contrasts from a Sample Annotation
#'
#' Two annotation layouts are in use and both have to work:
#'
#' * **explicit** — `Contrast` holds the linear-model expression and
#'   `ContrastName` its label. Taken as given.
#' * **dataset** — `Group` names the experimental group and `Control` marks it
#'   as control (`C`) or treatment (`T`). Every treatment-versus-control pair is
#'   generated, named `<treatment>_vs_<control>` and expressed with the `G_`
#'   prefix prolfqua gives group coefficients.
#'
#' Column names of the dataset layout are matched case-insensitively:
#' annotations in the wild spell the second column both `Control` and `CONTROL`.
#'
#' @param annot Data frame read from the annotation file.
#' @param annot_label Name used in error messages, normally the file name.
#' @return A named character vector of contrast expressions.
#' @keywords internal
derive_contrasts <- function(annot, annot_label = "the annotation") {
  if (all(c("Contrast", "ContrastName") %in% colnames(annot))) {
    contrasts <- annot$Contrast
    names(contrasts) <- annot$ContrastName
    contrasts <- contrasts[!is.na(contrasts) & nchar(contrasts) > 0]
  } else if (all(c("group", "control") %in% tolower(colnames(annot)))) {
    group_col <- colnames(annot)[tolower(colnames(annot)) == "group"][[1]]
    control_col <- colnames(annot)[tolower(colnames(annot)) == "control"][[1]]

    levels_df <- annot |>
      dplyr::select(dplyr::all_of(c(group_col, control_col))) |>
      dplyr::distinct()

    control_groups <- levels_df[[group_col]][toupper(levels_df[[control_col]]) == "C"]
    treatment_groups <- levels_df[[group_col]][toupper(levels_df[[control_col]]) == "T"]

    contrasts <- character()
    for (trt in sort(treatment_groups)) {
      for (ctrl in sort(control_groups)) {
        contrasts[paste0(trt, "_vs_", ctrl)] <- paste0("G_", trt, " - G_", ctrl)
      }
    }
  } else {
    stop(
      "Annotation file must have either (Contrast, ContrastName) or ",
      "(Group, Control) columns",
      call. = FALSE
    )
  }

  # An annotation can carry the right columns and still define nothing, for
  # example a template whose Group/Control cells are still empty. Fail here
  # rather than fitting contrasts that do not exist.
  if (length(contrasts) == 0) {
    stop(
      "No contrasts could be derived from ",
      annot_label,
      ": the annotation carries the expected columns but marks no control (C) ",
      "and treatment (T) groups.",
      call. = FALSE
    )
  }

  contrasts
}

#' Compute the CorrectFirst Analysis from the Bundled Example Data
#'
#' As [compute_dpa_dpu_example()], for the CorrectFirst report: it corrects the
#' example site abundances by their protein abundance and fits the model, then
#' drops the two wide tables the pipeline also leaves out of the report's
#' object, so the standalone render reads the same shape a run produces.
#'
#' @return [compute_cf_dea()]'s result without `wide_data` and
#'   `wide_annotation`.
#' @keywords internal
compute_cf_dea_example <- function() {
  dirs <- example_dea_pair()
  res <- compute_cf_dea(dirs$phospho, dirs$protein, dirs$annot_file)
  res[setdiff(names(res), c("wide_data", "wide_annotation"))]
}
