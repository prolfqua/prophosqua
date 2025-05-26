#’ @importFrom rlang .data
NULL

.reverse_join_column <- function(join_column){
  reverse_join_column <- vector(mode = "character", length(join_column))
  for (i in seq_along(join_column)) {
    reverse_join_column[i] <- if (names(join_column)[i] != "") {  names(join_column)[i]} else { join_column[i]}
    names(reverse_join_column)[i] <- if (names(join_column)[i] != "") { join_column[i]} else {""}
  }
  return(reverse_join_column)
}






#' extract sequence window from fasta files for plotting.
#' @export
get_sequence_windows <- function(phos_res, fasta_file, rev_pattern = "rev_", window_size = 15) {
  unique_prot_pep_seq <- phos_res |>
    dplyr::filter(.data$AllLocalized == TRUE) |>
    dplyr::select(protein_id = .data$protein_Id, .data$site, phos_sites = .data$PhosSites) |>
    dplyr::distinct()

  unique_prot_pep_seq <- unique_prot_pep_seq |>
    tidyr::separate_longer_delim(.data$phos_sites, delim = ";") |>
    dplyr::mutate(
      pos_in_protein = as.integer(
        stringr::str_remove(.data$phos_sites, "^[A-Z]")
      )
    ) |>
    dplyr::mutate(aa = stringr::str_remove(.data$phos_sites, "\\d+"))

  fasta <- prolfquapp::get_annot_from_fasta(fasta_file, pattern_decoys = rev_pattern, include_seq = TRUE)

  unique_prot_pep_seq <- dplyr::inner_join(unique_prot_pep_seq, fasta, by = c(protein_id = "proteinname"))
  unique_prot_pep_seq <- unique_prot_pep_seq |>
    dplyr::mutate(
      padded_sequence = paste0(strrep("X", window_size), .data$sequence, strrep("X", window_size)),
      pos_in_padded_seq = .data$pos_in_protein + window_size,  # Adjust position due to padding
      pos_start = .data$pos_in_padded_seq - window_size,
      pos_end = .data$pos_in_padded_seq + window_size,
      sequence_window = substr(.data$padded_sequence, start = .data$pos_start, stop = .data$pos_end)
    ) |>
    dplyr::select(-.data$padded_sequence, -.data$pos_in_padded_seq, -.data$pos_start, -.data$pos_end)

  unique_prot_pep_seq$sequence <- NULL
  unique_prot_pep_seq
}
