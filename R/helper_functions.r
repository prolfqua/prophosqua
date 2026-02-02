
#' Get DEA xlsx file path
#'
#' Constructs the path to a DEA results Excel file based on project structure.
#'
#' @param project_dir Path to the project directory
#' @param WU Work unit identifier (e.g., "phospho", "prot")
#' @param date Date string used in folder naming (e.g., "20240101")
#' @return Character string with full path to the Excel file
#' @export
dea_xlsx_path <- function(project_dir, WU, date) {
  file.path(
    project_dir,
    paste0("DEA_", date, "_",
           "WU",WU,"_vsn/Results_WU_",WU,"/DE_WU",WU,".xlsx"))
}

#' Get DEA results directory path
#'
#' Constructs the path to a DEA results directory based on project structure.
#'
#' @param project_dir Path to the project directory
#' @param WU Work unit identifier (e.g., "phospho", "prot")
#' @param date Date string used in folder naming (e.g., "20240101")
#' @return Character string with full path to the results directory
#' @export
dea_res_dir <- function(project_dir, WU, date) {
  file.path(
    project_dir,
    paste0("DEA_", date, "_",
           "WU",WU,"_vsn/Results_WU_",WU))
}


#' Filter contaminant proteins from data
#'
#' Removes rows matching FGCZ contaminants, contam_ prefix, or rev_ prefix.
#'
#' @param data Data frame with protein data
#' @param protein_id_col Name of the protein ID column (default: "protein_Id")
#' @return Filtered data frame
#' @export
filter_contaminants <- function(data, protein_id_col = "protein_Id") {
  data |>
    dplyr::filter(!grepl("FGCZCont", !!rlang::sym(protein_id_col))) |>
    dplyr::filter(!grepl("contam_", !!rlang::sym(protein_id_col))) |>
    dplyr::filter(!grepl("^rev_", !!rlang::sym(protein_id_col)))
}

#' Load and preprocess DEA data from Excel
#'
#' Reads an Excel file and validates that required columns are present.
#'
#' @param file_path Path to the Excel file
#' @param required_cols Character vector of required column names
#' @param sheet_name Name of the sheet to read (default: "diff_exp_analysis")
#' @return Data frame with the loaded data
#' @export
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

#' Explode multisites into individual rows
#'
#' Each phosphosite goes into its own row when multiple sites are present.
#'
#' @param combined_site Data frame with PhosSites column (semicolon-delimited)
#' @return Data frame with one row per site, with modAA and posInProtein columns
#' @importFrom rlang .data
#' @export
explode_multisites <- function(combined_site) {
  class(combined_site$startModSite) <- "numeric"
  class(combined_site$endModSite) <- "numeric"
  combined_site <- combined_site |>
    dplyr::mutate(
      PhosSites = dplyr::case_when(
        is.na(.data$PhosSites) ~ paste0("NotLoc", (.data$startModSite + .data$endModSite) / 2),
        TRUE ~ .data$PhosSites
      )
    )
  combined_site_prot_long <- combined_site |>
    tidyr::separate_longer_delim(.data$PhosSites, delim = ";")
  combined_site_prot_long <- combined_site_prot_long |>
    tidyr::extract(
      col = "PhosSites",
      into = c("modAA", "posInProtein"),
      regex = "([A-Za-z]+)([0-9]+)",
      remove = FALSE
    )
  return(combined_site_prot_long)
}
