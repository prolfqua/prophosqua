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
get_sequence_windows <- function(phosRes, fasta_file, rev_pattern = "rev_", window_size = 15) {
  uniqueProtPepSeq <- phosRes |> dplyr::filter(AllLocalized == TRUE) |>
    dplyr::select(protein_Id, site, PhosSites) |> dplyr::distinct()

  uniqueProtPepSeq <- uniqueProtPepSeq |> tidyr::separate_longer_delim(PhosSites, delim = ";")
  uniqueProtPepSeq <- uniqueProtPepSeq |> dplyr::mutate(posInProtein = as.integer(stringr::str_remove(PhosSites, "^[A-Z]")))
  uniqueProtPepSeq <- uniqueProtPepSeq |> dplyr::mutate(AA = stringr::str_remove(PhosSites, "\\d+"))


  fasta <- prolfquapp::get_annot_from_fasta(fasta_file, pattern_decoys = rev_pattern, include_seq = TRUE)

  uniqueProtPepSeq <- dplyr::inner_join(uniqueProtPepSeq, fasta, by = c(protein_Id = "proteinname") )
  uniqueProtPepSeq <- uniqueProtPepSeq |>
    dplyr::mutate(
      padded_sequence = paste0(strrep("X", window_size), sequence, strrep("X", window_size)),
      posInPaddedSeq = posInProtein + window_size,  # Adjust position due to padding
      posStart = posInPaddedSeq - window_size,
      posEnd = posInPaddedSeq + window_size,
      sequence_window = substr(padded_sequence, start = posStart, stop = posEnd)
    ) |>
    dplyr::select(-padded_sequence, -posInPaddedSeq, -posStart, -posEnd)

  uniqueProtPepSeq$sequence <- NULL
  uniqueProtPepSeq
}

