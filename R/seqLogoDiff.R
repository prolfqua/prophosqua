#' Create PWM from amino acid sequences
#'
#' Internal helper function to compute a position weight matrix (PWM) from
#' a vector of amino acid sequences.
#'
#' @param sequences Character vector of amino acid sequences (must be same length)
#' @return A matrix with amino acids as rows and positions as columns
#'
#' @keywords internal
get_pwm <- function(sequences) {
  mat <- do.call(rbind, strsplit(sequences, ""))
  aa <- c(
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
  )
  pwm <- apply(mat, 2, function(x) {
    table(factor(x, levels = aa)) / length(x)
  })
  rownames(pwm) <- aa
  return(pwm)
}


#' Plot Difference Sequence Logo (Upregulated - Downregulated)
#'
#' Creates a difference sequence logo plot that visualizes the difference in
#' amino acid frequencies between upregulated and downregulated phosphorylation
#' sites. This visualization helps identify kinase-specific motifs that are
#' differentially regulated between experimental conditions.
#'
#' @description
#' The function computes position weight matrices (PWMs) for both upregulated
#' and downregulated sites, then calculates their difference (Up - Down).
#' Positive values indicate amino acids enriched in upregulated sites, while
#' negative values indicate enrichment in downregulated sites.
#'
#' This approach is particularly useful for:
#' \itemize{
#'   \item Identifying kinases with opposing activity changes
#'   \item Detecting motif shifts between conditions
#'   \item Highlighting position-specific amino acid preferences
#' }
#'
#' @param sig_sites A data frame containing significant phosphorylation sites
#'   with the following required columns:
#'   \describe{
#'     \item{contrast}{Character. The name of the experimental contrast
#'       (e.g., "KO_vs_WT")}
#'     \item{regulation}{Character. Either "upregulated" or "downregulated"}
#'     \item{SequenceWindow}{Character. The amino acid sequence window
#'       surrounding the phosphorylation site (typically 15 amino acids
#'       centered on the modified residue)}
#'   }
#'
#' @return A ggplot2 object created by \code{ggseqlogo} showing difference
#'   logos for each contrast. Returns \code{NULL} if no contrasts have both
#'   upregulated and downregulated sequences.
#'
#' @details
#' The difference PWM is calculated as:
#' \deqn{PWM_{diff} = PWM_{up} - PWM_{down}}
#'
#' For each position in the sequence window, the function calculates the
#' frequency of each of the 20 standard amino acids. The difference in
#' frequencies is then visualized using the 'custom' method in ggseqlogo,
#' where:
#' \itemize{
#'   \item Letters above the baseline (positive values) are enriched in
#'     upregulated sites
#'   \item Letters below the baseline (negative values) are enriched in
#'     downregulated sites
#'   \item Letter height corresponds to the magnitude of the frequency
#'     difference
#' }
#'
#' @examples
#' # Create example data with phosphorylation site motifs
#' sig_sites <- data.frame(
#'   contrast = rep("KO_vs_WT", 20),
#'   regulation = rep(c("upregulated", "downregulated"), each = 10),
#'   SequenceWindow = c(
#'     # Upregulated: Proline-directed (CDK-like: S/T-P-x-K)
#'     "AAATPASPLTPGKRA", "GGGKSPVSPLTPGKR", "AAATPASPLTPGKRA",
#'     "LLLSQVSPKTPAAAA", "AAATPASPLTPGKRA", "GGGKSPVSPLTPGKR",
#'     "VVKTPASPLTPGKRA", "AAATPASPLTPGKRA", "LLLSQVSPKTPAAAA",
#'     "GGGKSPVSPLTPGKR",
#'
#'     # Downregulated: Basophilic (PKA-like: R-R-x-S)
#'     "LRRRRSLAAAGKPAA", "GGRRRASVAAAGKPA", "LRRRRSLAAAGKPAA",
#'     "GKRRRSSLAAAGKPA", "LRRRRSLAAAGKPAA", "GGRRRASVAAAGKPA",
#'     "DKRRRSSLAAAGKPA", "LRRRRSLAAAGKPAA", "GKRRRSSLAAAGKPA",
#'     "GGRRRASVAAAGKPA"
#'   )
#' )
#'
#' # Generate difference logo
#' diff_plot <- plot_diff_logo(sig_sites)
#' print(diff_plot)
#'
#' @seealso
#' \code{\link[ggseqlogo]{ggseqlogo}} for the underlying plotting function
#'
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom ggplot2 ggtitle
#' @importFrom dplyr filter pull
#' @export
plot_diff_logo <- function(sig_sites) {
  # Iterate over contrasts
  contrasts_list <- unique(sig_sites$contrast)
  diff_list <- list()

  for (cont in contrasts_list) {
    up_seqs <- sig_sites |>
      dplyr::filter(contrast == cont, regulation == "upregulated") |>
      dplyr::pull(SequenceWindow) |>
      toupper()

    down_seqs <- sig_sites |>
      dplyr::filter(contrast == cont, regulation == "downregulated") |>
      dplyr::pull(SequenceWindow) |>
      toupper()

    # Only compute if we have sequences
    if (length(up_seqs) > 0 && length(down_seqs) > 0) {
      pwm_up <- get_pwm(up_seqs)
      pwm_down <- get_pwm(down_seqs)

      # Calculate Difference (Up - Down)
      diff_pwm <- pwm_up - pwm_down

      diff_list[[paste0(cont, " (Up-Down)")]] <- diff_pwm
    }
  }

  if (length(diff_list) > 0) {
    # Plot using ggseqlogo with custom method for difference matrix
    p <- ggseqlogo::ggseqlogo(diff_list, ncol = 2, method = "custom") +
      ggplot2::ggtitle("Difference in Amino Acid Probabilities (Up - Down)")
    return(p)
  } else {
    return(NULL)
  }
}


#' Plot Sequence Logos with Difference Column
#'
#' Creates a 3-column sequence logo plot showing upregulated sites,
#' downregulated sites, and their difference (Up - Down) for each contrast.
#' The difference column uses an independent y-axis scale appropriate for
#' showing positive and negative frequency differences.
#'
#' @param sig_sites A data frame containing significant phosphorylation sites
#'   with required columns: \code{contrast}, \code{regulation}, and
#'   \code{SequenceWindow}.
#'
#' @return A patchwork object combining sequence logos with 3 columns per
#'   contrast row: upregulated, downregulated, and difference logos.
#'   Returns \code{NULL} if no contrasts have both upregulated and
#'   downregulated sequences.
#'
#' @details
#' For each contrast, three sequence logos are displayed:
#' \itemize{
#'   \item \strong{Upregulated}: Amino acid frequencies at sites with increased
#'     phosphorylation (y-axis: 0-1)
#'   \item \strong{Downregulated}: Amino acid frequencies at sites with decreased
#'     phosphorylation (y-axis: 0-1)
#'   \item \strong{Difference}: The difference in frequencies (Up - Down), where
#'     positive values indicate enrichment in upregulated sites (y-axis: symmetric
#'     around 0)
#' }
#'
#' @examples
#' # Create example data with phosphorylation site motifs
#' sig_sites <- data.frame(
#'   contrast = rep("KO_vs_WT", 20),
#'   regulation = rep(c("upregulated", "downregulated"), each = 10),
#'   SequenceWindow = c(
#'     # Upregulated: Proline-directed (CDK-like: S/T-P-x-K)
#'     "AAATPASPLTPGKRA", "GGGKSPVSPLTPGKR", "AAATPASPLTPGKRA",
#'     "LLLSQVSPKTPAAAA", "AAATPASPLTPGKRA", "GGGKSPVSPLTPGKR",
#'     "VVKTPASPLTPGKRA", "AAATPASPLTPGKRA", "LLLSQVSPKTPAAAA",
#'     "GGGKSPVSPLTPGKR",
#'
#'     # Downregulated: Basophilic (PKA-like: R-R-x-S)
#'     "LRRRRSLAAAGKPAA", "GGRRRASVAAAGKPA", "LRRRRSLAAAGKPAA",
#'     "GKRRRSSLAAAGKPA", "LRRRRSLAAAGKPAA", "GGRRRASVAAAGKPA",
#'     "DKRRRSSLAAAGKPA", "LRRRRSLAAAGKPAA", "GKRRRSSLAAAGKPA",
#'     "GGRRRASVAAAGKPA"
#'   )
#' )
#'
#' # Generate 3-column layout: Up, Down, Difference
#' p <- plot_seqlogo_with_diff(sig_sites)
#' print(p)
#'
#' @seealso \code{\link{plot_diff_logo}} for difference logos only
#'
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom dplyr filter pull
#' @importFrom patchwork wrap_plots
#' @export
plot_seqlogo_with_diff <- function(sig_sites) {
  contrasts_list <- unique(sig_sites$contrast)
  up_seq_list <- list()
  down_seq_list <- list()
  diff_list <- list()

  for (cont in contrasts_list) {
    up_seqs <- sig_sites |>
      dplyr::filter(contrast == cont, regulation == "upregulated") |>
      dplyr::pull(SequenceWindow) |>
      toupper()

    down_seqs <- sig_sites |>
      dplyr::filter(contrast == cont, regulation == "downregulated") |>
      dplyr::pull(SequenceWindow) |>
      toupper()

    if (length(up_seqs) > 0 && length(down_seqs) > 0) {
      # Store sequences directly for up/down (ggseqlogo handles scaling correctly)
      up_seq_list[[paste0(cont, "_up")]] <- up_seqs
      down_seq_list[[paste0(cont, "_down")]] <- down_seqs

      # Compute difference PWM for diff plot
      pwm_up <- get_pwm(up_seqs)
      pwm_down <- get_pwm(down_seqs)
      diff_list[[paste0(cont, "_diff")]] <- pwm_up - pwm_down
    }
  }

  if (length(up_seq_list) > 0) {
    # Up/Down plots: use sequences with probability method (correct scaling)
    p_up <- ggseqlogo::ggseqlogo(up_seq_list, ncol = 1, method = "probability", seq_type = "aa") +
      ggplot2::ggtitle("Upregulated") +
      ggplot2::theme(legend.position = "none")
    p_down <- ggseqlogo::ggseqlogo(down_seq_list, ncol = 1, method = "probability", seq_type = "aa") +
      ggplot2::ggtitle("Downregulated") +
      ggplot2::theme(legend.position = "none")

    # Diff plot: use custom matrix but with proper y-axis expansion
    # Keep legend only on this plot
    p_diff <- ggseqlogo::ggseqlogo(diff_list, ncol = 1, method = "custom") +
      ggplot2::ggtitle("Difference (Up-Down)") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.1))

    # Combine with patchwork, collect legend at bottom
    p <- patchwork::wrap_plots(p_up, p_down, p_diff, ncol = 3) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
    return(p)
  } else {
    return(NULL)
  }
}
