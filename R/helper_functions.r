

#' get xlsx path
#' @export
dea_xlsx_path <- function(project_dir, WU, date) {
  file.path(
    project_dir,
    paste0("DEA_", date, "_",
           "WU",WU,"_vsn/Results_WU_",WU,"/DE_WU",WU,".xlsx"))
}

#' get xlsx path
#' @export
dea_res_dir <- function(project_dir, WU, date) {
  file.path(
    project_dir,
    paste0("DEA_", date, "_",
           "WU",WU,"_vsn/Results_WU_",WU))
}


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
    fc_threshold = 0,
    impute_flag = "Imputed_Mean_moderated") {
  if(!contrast_name %in% unique(combined_site_prot_long$contrast)){
    stop(contrast_name, " not in ", paste0(unique(combined_site_prot_long$contrast)))
  }
  combined_site_prot_long <- combined_site_prot_long |>
    filter(contrast == contrast_name)
  combined_site_prot_long <- combined_site_prot_long |>
    group_by(protein_Id) |>
    filter(any(FDR.site < FDR_threshold & abs(diff.site) > fc_threshold)) |>
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
    fc_threshold = 0,
    impute_flag = "Imputed_Mean_moderated",
    protein_Id = "protein_Id") {
  data_combined_diff <- data_combined_diff |>
    filter(contrast == contrast_name)
  data_combined_diff <- data_combined_diff |>
    dplyr::group_by(!!dplyr::sym(protein_Id)) |>
    filter(any(FDR_I < FDR_threshold & abs(diff_diff) > fc_threshold)) |>
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

#' Calculate optimal grid layout for multi-panel plots
#' @param n_panels Number of panels to arrange
#' @return List with ncol and nrow for grid layout
#' @keywords internal
#' @examples
#' calculate_optimal_layout(6)  # Returns list(ncol=3, nrow=2)
#' calculate_optimal_layout(8)  # Returns list(ncol=3, nrow=3)
calculate_optimal_layout <- function(n_panels) {
  if (n_panels == 1) {
    return(list(ncol = 1, nrow = 1))
  }

  # Calculate square root and round up for columns
  # Prefer wider layouts (more columns than rows) for better readability
  ncol <- ceiling(sqrt(n_panels))
  nrow <- ceiling(n_panels / ncol)

  return(list(ncol = ncol, nrow = nrow))
}

#' Plot protein and PTM expression across all contrasts in multi-panel layout
#' @param combined_site_prot_long data frame with combined site and protein data
#' @param FDR_threshold FDR threshold for filtering significant sites (default 0.05)
#' @param fc_threshold Fold change threshold for filtering (default 0)
#' @param impute_flag Flag for imputed values (default "Imputed_Mean_moderated")
#' @return data frame with protein_Id, protein_length, n_contrasts, and multi-panel plot
#' @export
#' @examples
#' # data(combined_site_prot_data)
#' # result <- n_to_c_expression_multicontrast(combined_site_prot_data)
#' # print(result$plot[[1]])  # Display first protein's multi-contrast plot
n_to_c_expression_multicontrast <- function(
    combined_site_prot_long,
    FDR_threshold = 0.05,
    fc_threshold = 0,
    impute_flag = "Imputed_Mean_moderated") {

  # Get all unique contrasts
  all_contrasts <- unique(combined_site_prot_long$contrast)
  n_contrasts <- length(all_contrasts)

  message("Found ", n_contrasts, " contrasts: ", paste(all_contrasts, collapse = ", "))

  # Filter proteins that have at least one significant site in ANY contrast
  significant_proteins <- combined_site_prot_long |>
    dplyr::group_by(protein_Id) |>
    dplyr::filter(any(FDR.site < FDR_threshold & abs(diff.site) > fc_threshold)) |>
    dplyr::ungroup()

  # Get unique proteins with their lengths
  proteins_to_plot <- significant_proteins |>
    dplyr::select(protein_Id, protein_length) |>
    dplyr::distinct()

  message("Creating multi-contrast plots for ", nrow(proteins_to_plot), " proteins...")

  # Prepare required columns
  required_cols <- c("protein_Id", "contrast", "protein_length", "site",
                     "diff.protein", "diff.site", "FDR.site",
                     "posInProtein", "modAA", "imputation_status")

  # Add imputation status
  combined_site_prot_long <- combined_site_prot_long |>
    dplyr::mutate(imputation_status = case_when(
      modelName.site == impute_flag ~ "imputed",
      TRUE ~ "observed"
    ))

  # Calculate optimal layout
  layout <- calculate_optimal_layout(n_contrasts)

  # Initialize results
  plot_data <- proteins_to_plot
  plot_data$n_contrasts <- n_contrasts
  plot_data$plot <- vector(mode = "list", length = nrow(plot_data))

  # Create multi-panel plots for each protein
  pb <- txtProgressBar(min = 0, max = nrow(plot_data), style = 3)

  for (i in seq_len(nrow(plot_data))) {
    setTxtProgressBar(pb, i)

    current_protein <- plot_data$protein_Id[[i]]
    current_length <- plot_data$protein_length[[i]]

    # Get data for this protein across all contrasts
    protein_data <- combined_site_prot_long |>
      dplyr::filter(protein_Id == current_protein)

    # Create individual plots for each contrast
    contrast_plots <- list()

    for (j in seq_along(all_contrasts)) {
      contrast <- all_contrasts[j]

      # Get data for this specific contrast
      contrast_data <- protein_data |>
        dplyr::filter(contrast == !!contrast) |>
        dplyr::select(dplyr::all_of(required_cols))

      # Create plot if data exists
      if (nrow(contrast_data) > 0) {
        contrast_plots[[j]] <- n_to_c_plot(
          contrast_data,
          current_protein,
          current_length,
          contrast
        )
      } else {
        # Create empty placeholder plot
        contrast_plots[[j]] <- ggplot2::ggplot() +
          ggplot2::annotate("text",
                           x = 0.5, y = 0.5,
                           label = paste0("No data for\n", contrast),
                           size = 5, color = "gray50") +
          ggplot2::theme_void() +
          ggplot2::labs(title = contrast)
      }
    }

    # Combine plots using patchwork
    combined_plot <- patchwork::wrap_plots(
      contrast_plots,
      ncol = layout$ncol,
      nrow = layout$nrow
    ) +
      patchwork::plot_annotation(
        title = paste0("Protein: ", current_protein, " (length: ", current_length, ")"),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
      )

    plot_data$plot[[i]] <- combined_plot
  }

  close(pb)
  message("Created ", nrow(plot_data), " multi-contrast plots")

  return(plot_data)
}

#' Plot PTM usage across all contrasts in multi-panel layout
#' @param data_combined_diff data frame with combined differential analysis results
#' @param FDR_threshold FDR threshold for filtering significant sites (default 0.05)
#' @param fc_threshold Fold change threshold for filtering (default 0)
#' @param impute_flag Flag for imputed values (default "Imputed_Mean_moderated")
#' @param protein_Id Column name for protein identifier (default "protein_Id")
#' @return data frame with protein_Id, protein_length, n_contrasts, and multi-panel plot
#' @export
#' @examples
#' # data(combined_diff_data)
#' # result <- n_to_c_usage_multicontrast(combined_diff_data)
#' # print(result$plot[[1]])  # Display first protein's multi-contrast plot
n_to_c_usage_multicontrast <- function(
    data_combined_diff,
    FDR_threshold = 0.05,
    fc_threshold = 0,
    impute_flag = "Imputed_Mean_moderated",
    protein_Id = "protein_Id") {

  # Get all unique contrasts
  all_contrasts <- unique(data_combined_diff$contrast)
  n_contrasts <- length(all_contrasts)

  message("Found ", n_contrasts, " contrasts: ", paste(all_contrasts, collapse = ", "))

  # Filter proteins that have at least one significant site in ANY contrast
  significant_proteins <- data_combined_diff |>
    dplyr::group_by(!!dplyr::sym(protein_Id)) |>
    dplyr::filter(any(FDR_I < FDR_threshold & abs(diff_diff) > fc_threshold)) |>
    dplyr::ungroup()

  # Get unique proteins with their lengths
  proteins_to_plot <- significant_proteins |>
    dplyr::select(dplyr::all_of(c(protein_Id, "protein_length"))) |>
    dplyr::distinct()

  message("Creating multi-contrast plots for ", nrow(proteins_to_plot), " proteins...")

  # Prepare required columns
  required_columns <- c("protein_Id", "contrast", "protein_length", "site",
                        "diff_diff", "FDR_I", "posInProtein", "modAA",
                        "imputation_status")

  # Add imputation status
  data_combined_diff$imputation_status <- ifelse(
    data_combined_diff$modelName.site == impute_flag |
      data_combined_diff$modelName.protein == impute_flag,
    "imputed", "observed"
  )

  # Calculate optimal layout
  layout <- calculate_optimal_layout(n_contrasts)

  # Initialize results
  plot_data <- proteins_to_plot
  plot_data$n_contrasts <- n_contrasts
  plot_data$plot <- vector(mode = "list", length = nrow(plot_data))

  # Create multi-panel plots for each protein
  pb <- txtProgressBar(min = 0, max = nrow(plot_data), style = 3)

  for (i in seq_len(nrow(plot_data))) {
    setTxtProgressBar(pb, i)

    current_protein <- plot_data[[protein_Id]][[i]]
    current_length <- plot_data$protein_length[[i]]

    # Get data for this protein across all contrasts
    protein_data <- data_combined_diff |>
      dplyr::filter(!!dplyr::sym(protein_Id) == current_protein)

    # Create individual plots for each contrast
    contrast_plots <- list()

    for (j in seq_along(all_contrasts)) {
      contrast <- all_contrasts[j]

      # Get data for this specific contrast
      contrast_data <- protein_data |>
        dplyr::filter(contrast == !!contrast) |>
        dplyr::select(dplyr::all_of(required_columns))

      # Create plot if data exists
      if (nrow(contrast_data) > 0) {
        contrast_plots[[j]] <- n_to_c_plot_integrated(
          contrast_data,
          current_protein,
          current_length,
          contrast
        )
      } else {
        # Create empty placeholder plot
        contrast_plots[[j]] <- ggplot2::ggplot() +
          ggplot2::annotate("text",
                           x = 0.5, y = 0.5,
                           label = paste0("No data for\n", contrast),
                           size = 5, color = "gray50") +
          ggplot2::theme_void() +
          ggplot2::labs(title = contrast)
      }
    }

    # Combine plots using patchwork
    combined_plot <- patchwork::wrap_plots(
      contrast_plots,
      ncol = layout$ncol,
      nrow = layout$nrow
    ) +
      patchwork::plot_annotation(
        title = paste0("Protein: ", current_protein, " (length: ", current_length, ")"),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
      )

    plot_data$plot[[i]] <- combined_plot
  }

  close(pb)
  message("Created ", nrow(plot_data), " multi-contrast plots")

  return(plot_data)
}


