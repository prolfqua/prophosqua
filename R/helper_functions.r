
#' filter contaminants
#' @export
filter_contaminants <- function(data, protein_id_col = "protein_Id") {
  data |>
    dplyr::filter(!grepl("FGCZCont", !!sym(protein_id_col))) |>
    dplyr::filter(!grepl("contam_", !!sym(protein_id_col))) |>
    dplyr::filter(!grepl("^rev_", !!sym(protein_id_col)))
}

#' load_and_preprocess_data
#' @export
#'
load_and_preprocess_data <- function(
    file_path,
    required_cols,
    sheet_name = "diff_exp_analysis") {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Read all columns as text to avoid type guessing warnings
  data <- readxl::read_xlsx(
    path = file_path,
    sheet = sheet_name,
    guess_max = 20000
  )

  # Validate required columns
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  return(data)
}

#' each site goes into own row, if more then on site per multisite
#' @export
explode_multisites <- function(combined_site){
  class(combined_site$startModSite) <- "numeric"
  class(combined_site$endModSite) <- "numeric"
  combined_site <- combined_site |> mutate(PhosSites  = case_when(is.na(PhosSites) ~ paste0("NotLoc",  (startModSite + endModSite)/2) , TRUE ~ PhosSites))
  combined_site_prot_long <- combined_site |> tidyr::separate_longer_delim(PhosSites, delim = ";")
  combined_site_prot_long <- combined_site_prot_long |> tidyr::extract(
    col = PhosSites,
    into = c("modAA", "posInProtein"),
    regex = "([A-Za-z]+)([0-9]+)",
    remove = FALSE
  )
  return(combined_site_prot_long)
}

#' plot protein and PTM expression
#' @export
n_to_c_expression  <- function(
    combined_site_prot_long,
    contrast_name,
    FDR_threshold = 0.05,
    impute_flag = "Imputed_Mean_moderated") {
  if(!contrast_name %in% unique(combined_site_prot_long$contrast)){
    stop(contrast_name, " not in ", paste0(unique(combined_site_prot_long$contrast)))
  }
  combined_site_prot_long <- combined_site_prot_long |>
    filter(contrast == contrast_name)
  combined_site_prot_long <- combined_site_prot_long |>
    group_by(protein_Id) |>
    filter(any(FDR.site < FDR_threshold)) |>
    ungroup()

  required_cols <- c("protein_Id",
                     "contrast",
                     "protein_length",
                     "site", "diff.protein", "diff.site", "FDR.site",
                     "posInProtein", "modAA", "imputation_status"
  )

  plot_data <- combined_site_prot_long |>
    mutate(imputation_status = case_when(modelName.site == impute_flag ~ "imputed", TRUE ~ "observed")) |>
    dplyr::select(all_of(required_cols)) |>
    group_by(.data$protein_Id, .data$contrast, .data$protein_length) |>
    tidyr::nest()

  plot_data$plot <- vector(mode = "list", length = nrow(plot_data))

  pb <- txtProgressBar(min = 0, max = nrow(plot_data), style = 3)
  for (i in seq_len(nrow(plot_data))) {
    setTxtProgressBar(pb, i)
    plot_data$plot[[i]] <- n_to_c_plot(
      plot_data$data[[i]],
      plot_data$protein_Id[[i]],
      plot_data$protein_length[[i]],
      plot_data$contrast[[i]]
    )

  }
  close(pb)

  return(plot_data=plot_data)
}


#' plot PTM occupation, use after adjusting for total proteome
#' @export
n_to_c_usage  <- function(
    data_combined_diff,
    contrast_name,
    FDR_threshold = 0.05,
    impute_flag = "Imputed_Mean_moderated",
    protein_Id = "protein_Id") {
  data_combined_diff <- data_combined_diff |>
    filter(contrast == contrast_name)
  data_combined_diff <- data_combined_diff |>
    dplyr::group_by(!!dplyr::sym(protein_Id)) |>
    filter(any(FDR_I < FDR_threshold)) |>
    dplyr::ungroup()

  required_columns <- c("protein_Id",
                        "contrast",
                        "protein_length",
                        "site",
                        "diff_diff",
                        "FDR_I",
                        "posInProtein",
                        "modAA",
                        "imputation_status")
  data_combined_diff$imputation_status <- ifelse(data_combined_diff$modelName.site == "Imputed_Mean_moderated" | data_combined_diff$modelName.protein == "Imputed_Mean_moderated", "imputed", "observed")

  # Create integrated N-to-C plots
  # Prepare data for plotting
  plot_data <- data_combined_diff |>
    dplyr::select(dplyr::all_of(required_columns)) |>
    dplyr::group_by(.data$protein_Id, .data$contrast, .data$protein_length) |>
    tidyr::nest()

  plot_data$plot <- vector(mode = "list", length = nrow(plot_data))


  # Create plots
  message("Creating plots...")
  pb <- txtProgressBar(min = 0, max = nrow(plot_data), style = 3)
  for (i in seq_len(nrow(plot_data))) {
    setTxtProgressBar(pb, i)
    plot_data$plot[[i]] <- n_to_c_plot_integrated(
      plot_data$data[[i]],
      plot_data$protein_Id[[i]],
      plot_data$protein_length[[i]],
      plot_data$contrast[[i]]
    )
  }
  close(pb)
  return(plot_data)
}


