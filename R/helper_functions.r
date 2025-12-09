

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
