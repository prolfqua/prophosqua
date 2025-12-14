#' Make FASTA database summary
#'
#' Generates database statistics including number of sequences and amino acid frequencies.
#' Extracted from prozor package to reduce dependencies.
#'
#' @param resDB A list of sequences from \code{prozor::readPeptideFasta}
#' @param old Logical; use old (slower) method for AA frequency calculation
#' @param as_string Logical; return formatted string (TRUE) or list (FALSE)
#' @return If \code{as_string=TRUE}, an array of strings to pass to \code{cat()}.
#'   If \code{as_string=FALSE}, a list with nrSequences, lengthSummary, and aafreq.
#' @export
#' @examples
#' \dontrun{
#' fasta <- prozor::readPeptideFasta("path/to/file.fasta")
#' cat(make_fasta_summary(fasta))
#' }
make_fasta_summary <- function(resDB, old = FALSE, as_string = TRUE) {
    if (old) {
        bigstr <- paste(resDB, collapse = "")
        vec <- strsplit(bigstr, split = "")[[1]]
        aafreq <- table(vec)
    } else {
        res <- list()
        tmp <- lapply(resDB, function(x) {
            table(strsplit(x, split = "")[[1]])
        })
        for (i in seq_along(tmp)) {
            x <- tmp[[i]]
            for (j in names(x)) {
                if (is.null(res[[j]])) {
                    res[[j]] <- as.numeric(x[j])
                } else {
                    res[[j]] <- res[[j]] + as.numeric(x[j])
                }
            }
        }
        aafreq <- unlist(res)
        aafreq <- aafreq[order(names(aafreq))]
    }
    length_s <- summary(vapply(resDB, seqinr::getLength, numeric(1)))

    if (as_string) {
        aafreq <- paste(utils::capture.output(as.matrix(aafreq)), "\n", sep = "")
        length_s <- paste(utils::capture.output(length_s), "\n", sep = "")
        summaryRes <- c(
            "nr sequences:\n", length(resDB), "\n length summary:\n",
            length_s, "AA frequencies:\n", aafreq
        )
    } else {
        summaryRes <- list()
        summaryRes$nrSequences <- length(resDB)
        summaryRes$lengthSummary <- length_s
        summaryRes$aafreq <- aafreq
    }
    return(summaryRes)
}
