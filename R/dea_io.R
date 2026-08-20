#' Locate and Normalize prolfquapp DEA Outputs
#'
#' Helpers that resolve the files a prolfquapp DEA run writes, and normalize the
#' two identifiers that differ between runs: the sample column, whose name the
#' run declares in its YAML, and the protein identifier, which may arrive as a
#' bare accession or in FASTA `sp|P12345|PROT_HUMAN` form.
#'
#' A DEA output directory holds a `Results_WU_<workunit>` subdirectory with the
#' results workbook, the normalized abundances, and the analysis configuration.
#'
#' @name dea_io
#' @keywords internal
NULL

#' Get Excel File from DEA Directory
#'
#' Finds the Excel results file within a DEA directory's `Results_WU_*`
#' subfolder, preferring the `DE_`-prefixed workbook that prolfquapp writes.
#' Falls back to searching `dea_dir` itself when no such subfolder exists.
#'
#' @param dea_dir Path to DEA output directory
#' @return Full path to the Excel file
#' @export
#' @examples
#' dea_dir <- file.path(tempdir(), "DEA_example")
#' results <- file.path(dea_dir, "Results_WU_example")
#' dir.create(results, recursive = TRUE, showWarnings = FALSE)
#' file.create(file.path(results, c("DE_example.xlsx", "other.xlsx")))
#'
#' basename(get_dea_xlsx(dea_dir))
get_dea_xlsx <- function(dea_dir) {
  # Find Results_WU_* subdirectory
  results_dirs <- list.dirs(dea_dir, recursive = FALSE, full.names = TRUE)
  results_dir <- results_dirs[grepl("^Results_WU_", basename(results_dirs))]

  if (length(results_dir) == 0) {
    # Fallback: search in main directory
    results_dir <- dea_dir
  } else {
    results_dir <- results_dir[1]
  }

  # Find xlsx file in directory
  xlsx_files <- list.files(results_dir, pattern = "\\.xlsx$", full.names = TRUE)

  if (length(xlsx_files) == 0) {
    stop("No Excel file found in: ", results_dir)
  }

  # Prefer the one with "DE_" prefix (standard prolfquapp output)
  de_files <- xlsx_files[grepl("^DE_", basename(xlsx_files))]
  if (length(de_files) > 0) {
    message("Using: ", basename(de_files[1]))
    return(de_files[1])
  }

  message("Using: ", basename(xlsx_files[1]))
  return(xlsx_files[1])
}

#' Get File from DEA Directory
#'
#' Generic helper to find files within a DEA directory's `Results_WU_*`
#' subfolder.
#'
#' @param dea_dir Path to DEA output directory
#' @param filename Filename to find (e.g., "lfqdata_normalized.parquet")
#' @param description Description for error message (e.g., "parquet file")
#' @return Full path to the file
#' @export
#' @examples
#' dea_dir <- file.path(tempdir(), "DEA_example")
#' results <- file.path(dea_dir, "Results_WU_example")
#' dir.create(results, recursive = TRUE, showWarnings = FALSE)
#' file.create(file.path(results, "lfqdata.yaml"))
#'
#' basename(get_dea_file(dea_dir, "lfqdata.yaml", "yaml file"))
get_dea_file <- function(dea_dir, filename, description = "file") {
  pattern <- file.path(dea_dir, "Results_WU_*", filename)
  matches <- Sys.glob(pattern)
  if (length(matches) == 0) {
    stop("No ", description, " found in: ", dea_dir)
  }
  return(matches[1])
}

#' Get Parquet File from DEA Directory
#'
#' @param dea_dir Path to DEA output directory
#' @return Full path to the parquet file
#' @export
#' @examples
#' dea_dir <- file.path(tempdir(), "DEA_example")
#' results <- file.path(dea_dir, "Results_WU_example")
#' dir.create(results, recursive = TRUE, showWarnings = FALSE)
#' file.create(file.path(results, "lfqdata_normalized.parquet"))
#'
#' basename(get_dea_parquet(dea_dir))
get_dea_parquet <- function(dea_dir) {
  get_dea_file(dea_dir, "lfqdata_normalized.parquet", "parquet file")
}

#' Get YAML Config from DEA Directory
#'
#' @param dea_dir Path to DEA output directory
#' @return Full path to the yaml file
#' @export
#' @examples
#' dea_dir <- file.path(tempdir(), "DEA_example")
#' results <- file.path(dea_dir, "Results_WU_example")
#' dir.create(results, recursive = TRUE, showWarnings = FALSE)
#' file.create(file.path(results, "lfqdata.yaml"))
#'
#' basename(get_dea_yaml(dea_dir))
get_dea_yaml <- function(dea_dir) {
  get_dea_file(dea_dir, "lfqdata.yaml", "yaml file")
}

#' Get the Declared Sample Column from a DEA YAML File
#'
#' prolfquapp records the column holding the sample name in its analysis
#' configuration. Older runs call it `Name`, current ones `sampleName`, so the
#' name is read rather than assumed.
#'
#' @param yaml_file Path to a prolfqua `lfqdata.yaml` file
#' @return The sample column name
#' @export
#' @examples
#' yaml_file <- file.path(tempdir(), "lfqdata.yaml")
#' writeLines("sample_name: sampleName", yaml_file)
#'
#' get_sample_name_column(yaml_file)
get_sample_name_column <- function(yaml_file) {
  config <- yaml::read_yaml(yaml_file)
  sample_col <- config$sample_name

  if (is.null(sample_col) || length(sample_col) != 1 || !nzchar(sample_col)) {
    stop("No valid sample_name declared in: ", yaml_file)
  }

  sample_col
}

#' Get the Declared Sample Column from a DEA Directory
#'
#' @param dea_dir Path to a DEA output directory
#' @return The sample column name
#' @export
#' @examples
#' dea_dir <- file.path(tempdir(), "DEA_example")
#' results <- file.path(dea_dir, "Results_WU_example")
#' dir.create(results, recursive = TRUE, showWarnings = FALSE)
#' writeLines("sample_name: sampleName", file.path(results, "lfqdata.yaml"))
#'
#' get_dea_sample_name_column(dea_dir)
get_dea_sample_name_column <- function(dea_dir) {
  get_sample_name_column(get_dea_yaml(dea_dir))
}

#' Canonicalize the Sample Column in Normalized DEA Data
#'
#' Renames the sample column a DEA run declares to a fixed name, so that data
#' from runs using different schemas can be joined.
#'
#' @param data Normalized DEA data
#' @param yaml_file Path to its `lfqdata.yaml` file
#' @param canonical_name Output column name
#' @return `data` with its declared sample column renamed
#' @export
#' @examples
#' yaml_file <- file.path(tempdir(), "lfqdata.yaml")
#' writeLines("sample_name: sampleName", yaml_file)
#' data <- data.frame(sampleName = c("a", "b"), abundance = c(1, 2))
#'
#' canonicalize_dea_sample_column(data, yaml_file)
canonicalize_dea_sample_column <- function(data, yaml_file, canonical_name = "Name") {
  sample_col <- get_sample_name_column(yaml_file)

  if (!sample_col %in% names(data)) {
    stop("Declared sample column '", sample_col, "' is absent from normalized DEA data")
  }
  if (sample_col == canonical_name) {
    return(data)
  }
  if (canonical_name %in% names(data)) {
    stop("Cannot rename sample column '", sample_col, "' to '", canonical_name, "': both columns are present")
  }

  names(data)[names(data) == sample_col] <- canonical_name
  data
}

#' Canonicalize UniProt Protein Identifiers
#'
#' Converts FASTA-style identifiers such as `sp|P12345|PROT_HUMAN` to their
#' accession while leaving bare accessions unchanged. The mapping must remain
#' one-to-one after upstream decoy filtering, otherwise two distinct input
#' identifiers would be joined as one protein; a violation is an error.
#'
#' @param data Data containing protein identifiers
#' @param id_col Protein identifier column
#' @return `data` with canonicalized protein identifiers
#' @export
#' @examples
#' data <- data.frame(
#'   protein_Id = c("sp|P12345|PROT_HUMAN", "Q67890"),
#'   diff = c(1.2, -0.4)
#' )
#'
#' canonicalize_uniprot_ids(data)
canonicalize_uniprot_ids <- function(data, id_col = "protein_Id") {
  if (!id_col %in% names(data)) {
    stop("Protein identifier column is absent: ", id_col)
  }

  original <- as.character(data[[id_col]])
  canonical <- vapply(
    strsplit(original, "|", fixed = TRUE),
    function(parts) if (length(parts) >= 2) parts[[2]] else parts[[1]],
    character(1)
  )

  id_map <- unique(data.frame(original = original, canonical = canonical))
  ambiguous <- id_map$canonical[duplicated(id_map$canonical)]
  if (length(ambiguous) > 0) {
    stop("Protein identifier canonicalization is not one-to-one for: ", paste(unique(ambiguous), collapse = ", "))
  }

  data[[id_col]] <- canonical
  data
}

#' Load PTM Site Metadata from a DEA Result
#'
#' Uses metadata exported in the normalized-abundance sheet when available.
#' FragPipe single-site DEA currently omits its sequence window there, so the
#' preserved input table is used as the authoritative fallback.
#'
#' @param dea_dir Path to a phosphosite DEA output directory
#' @return One row per site with position, modified residue, and sequence window
#' @export
#' @examplesIf requireNamespace("writexl", quietly = TRUE)
#' dea_dir <- file.path(tempdir(), "DEA_sites")
#' results <- file.path(dea_dir, "Results_WU_sites")
#' dir.create(results, recursive = TRUE, showWarnings = FALSE)
#' writexl::write_xlsx(
#'   list(normalized_abundances = data.frame(
#'     site = c("P12345_S10~PEP", "P12345_T20~PEP"),
#'     posInProtein = c(10L, 20L),
#'     modAA = c("S", "T"),
#'     SequenceWindow = c("AAAAAAASAAAAAAA", "AAAAAAATAAAAAAA"),
#'     protein_Id = "P12345",
#'     gene_name = "GENE",
#'     protein_length = 100L
#'   )),
#'   file.path(results, "DE_sites.xlsx")
#' )
#'
#' get_dea_ptm_site_info(dea_dir)
get_dea_ptm_site_info <- function(dea_dir) {
  full_results <- readxl::read_xlsx(
    get_dea_xlsx(dea_dir),
    sheet = "normalized_abundances"
  )
  seq_col <- grep(
    "FlankingRegion|SequenceWindow",
    names(full_results),
    value = TRUE,
    ignore.case = TRUE
  )
  site_col <- intersect(c("site", "protein_Id_site"), names(full_results))
  pos_col <- intersect(c("posInProtein", "PTM_SiteLocation"), names(full_results))
  aa_col <- intersect(c("modAA", "PTM_SiteAA"), names(full_results))

  if (length(seq_col) > 0 && length(site_col) > 0 && length(pos_col) > 0 && length(aa_col) > 0) {
    site_info <- full_results |>
      dplyr::transmute(
        site = .data[[site_col[[1]]]],
        posInProtein = .data[[pos_col[[1]]]],
        modAA = .data[[aa_col[[1]]]],
        SequenceWindow = .data[[seq_col[[1]]]],
        .data$protein_Id,
        .data$gene_name,
        .data$protein_length
      ) |>
      dplyr::distinct()
  } else {
    input_files <- list.files(
      dea_dir,
      pattern = "^abundance_single-site_None\\.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(input_files) != 1) {
      stop("Expected one preserved FragPipe single-site input in ", dea_dir, "; found ", length(input_files))
    }

    site_info <- readr::read_tsv(
      input_files[[1]],
      col_select = tidyselect::all_of(
        c("Index", "ProteinID", "Peptide", "SequenceWindow", "Start")
      ),
      show_col_types = FALSE
    ) |>
      dplyr::transmute(
        site = paste(.data$Index, .data$Peptide, sep = "~"),
        posInProtein = .data$Start,
        modAA = sub(".*_([A-Za-z])[0-9]+$", "\\1", .data$Index),
        .data$SequenceWindow,
        protein_Id = .data$ProteinID
      ) |>
      dplyr::distinct()

    protein_info <- full_results |>
      dplyr::select("protein_Id", "gene_name", "protein_length") |>
      dplyr::distinct()
    site_info <- dplyr::left_join(site_info, protein_info, by = "protein_Id")
  }

  if (anyDuplicated(site_info$site)) {
    stop("PTM site metadata is not unique by site in: ", dea_dir)
  }
  site_info
}
