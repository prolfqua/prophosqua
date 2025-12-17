#' @importFrom ggplot2 ggplot geom_segment aes scale_linetype_manual
#' @importFrom ggplot2 annotate scale_color_manual scale_x_continuous
#' @importFrom ggplot2 geom_text labs theme_minimal geom_rect scale_fill_manual guides guide_legend
NULL

#' Example dataset exampleN_C_dat
#'
#' A simple example dataset demonstrating groups and numeric values.
#' @format ## `exampleN_C_dat`
"exampleN_C_dat"

#' Example dataset n_c_integrated_df
#'
#' A simple example dataset demonstrating groups and numeric values.
#' @format ## `n_c_integrated_df`
"n_c_integrated_df"

#' Example dataset exampleN_C_dat_no_prot
#'
#' A simple example dataset demonstrating groups and numeric values.
#' @format ## `exampleN_C_dat_no_prot`
"exampleN_C_dat_no_prot"

#' Example phosphoproteomics differential expression results
#'
#' Combined differential expression results from phosphoproteomics analysis
#' of Atg16l1 knockout vs wild-type macrophages. Contains both site-level (DPA)
#' and protein-corrected (DPU) statistics for 4 contrasts.
#'
#' @format A tibble with 105,824 rows and 56 columns including:
#' \describe{
#'   \item{protein_Id}{UniProt protein identifier}
#'   \item{gene_name.site}{Gene symbol}
#'   \item{modAA}{Modified amino acid (S, T, or Y)}
#'   \item{posInProtein}{Position of modification in protein}
#'   \item{SequenceWindow}{Flanking sequence around modification site}
#'   \item{contrast}{Comparison name (e.g., KO_vs_WT)}
#'   \item{diff.site}{Log2 fold change for site (DPA)}
#'   \item{statistic.site}{t-statistic for site}
#'   \item{p.value.site}{p-value for site}
#'   \item{FDR.site}{FDR-adjusted p-value for site}
#'   \item{diff_diff}{Protein-corrected log2 fold change (DPU)}
#'   \item{tstatistic_I}{t-statistic for DPU}
#'   \item{pValue_I}{p-value for DPU}
#'   \item{FDR_I}{FDR-adjusted p-value for DPU}
#' }
#' @source Maculins et al. (2020) eLife, doi:10.7554/elife.62320
"combined_test_diff_example"
