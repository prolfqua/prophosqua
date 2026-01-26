#’ @importFrom rlang .data
NULL


#' extract sequence windows, from sequence at pos in protein
#' @export
#' @param unique_prot_pep_seq data.frame with sequence and pos in protein
#' @param sequence name of column with sequence
#' @param pos_in_protein name of column with position in protein
#' @param flank_size size of the window to extract
#' @examples
#' # Create sample data
#' library(tidyverse)
#' sample_data <- data.frame(
#'   protein_id = c("P12345", "P12345", "Q67890"),
#'   sequence = c("MKFLVLLFNILCLFPVLAADNH", "MKFLVLLFNILCLFPVLAADNH", "AEQKLISEEDLLRKRREQLKHKLEQL"),
#'   pos_in_protein = c(5, 12, 8),
#'   peptide = c("FLV", "ILC", "EED")
#' )
#'
#' # Extract sequence windows with default flank size (7)
#' result <- get_sequence_windows(sample_data)
#' stopifnot(all(result$sequence_window == c( "XXXMKFLVLLFNILC", "VLLFNILCLFPVLAA", "AEQKLISEEDLLRKR")))
#'
#' # Extract sequence windows with custom flank size
#' result_small <- get_sequence_windows(sample_data, flank_size = 3)
#' stopifnot(all(result_small$sequence_window == c( "KFLVLLF",  "NILCLFP",  "LISEEDL")))
#'
#' # Extract sequence windows with different column names
#' sample_data2 <- data.frame(
#'   prot_seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRIY",
#'   position = c(10, 25, 40)
#' )
#'
#' result2 <- get_sequence_windows(sample_data2,
#'                                flank_size = 5,
#'                                sequence = "prot_seq",
#'                                pos_in_protein = "position")
#' print(result2$sequence_window)
get_sequence_windows <- function(unique_prot_pep_seq,
                                 flank_size = 7,
                                 sequence = "sequence",
                                 pos_in_protein = "pos_in_protein"){
  half_window <- flank_size

  unique_prot_pep_seq_2 <- unique_prot_pep_seq |>
    dplyr::mutate(
      padded_sequence = paste0(strrep("X", half_window), !!rlang::sym(sequence), strrep("X", half_window)),
      pos_in_padded_seq = !!rlang::sym(pos_in_protein) + half_window,
      pos_start = !!rlang::sym("pos_in_padded_seq") - half_window,
      pos_end = !!rlang::sym("pos_in_padded_seq") + half_window,  # -1 because we want window_size chars
      sequence_window = substr(!!rlang::sym("padded_sequence"), start = pos_start, stop = pos_end)
    )
  unique_prot_pep_seq <- unique_prot_pep_seq_2 |>
    dplyr::select(-.data$padded_sequence, -.data$pos_in_padded_seq, -.data$pos_start, -.data$pos_end)
  return(unique_prot_pep_seq)
}
