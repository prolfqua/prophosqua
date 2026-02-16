# Filter significant PTM sites

Filters phosphosite data by FDR and fold-change thresholds and adds a
regulation column indicating direction of change. This prepares data for
downstream analysis like sequence logo visualization with
[`plot_diff_logo`](https://prolfqua.github.io/prophosqua/reference/plot_diff_logo.md)
and
[`plot_seqlogo_with_diff`](https://prolfqua.github.io/prophosqua/reference/plot_seqlogo_with_diff.md).

## Usage

``` r
filter_significant_sites(
  data,
  fdr_col = "FDR.site",
  diff_col = "diff.site",
  fdr_threshold = 0.05,
  fc_threshold = 0.6,
  require_sequence = FALSE
)
```

## Arguments

- data:

  Data frame with PTM results containing FDR and fold-change columns

- fdr_col:

  Character. Name of FDR/adjusted p-value column. Default "FDR.site"

- diff_col:

  Character. Name of log2 fold-change column. Default "diff.site"

- fdr_threshold:

  Numeric. FDR cutoff (sites with FDR \< threshold kept). Default 0.05

- fc_threshold:

  Numeric. Absolute log2 fold-change cutoff (sites with \|FC\| \>
  threshold kept). Default 0.6 (approximately 1.5-fold change)

- require_sequence:

  Logical. If TRUE, filter out rows with invalid SequenceWindow (NA or
  starting/ending with underscore). Default FALSE

## Value

Filtered data frame with added 'regulation' column containing
"upregulated" or "downregulated" based on sign of fold-change

## Details

The function performs two operations:

1.  Filters to significant sites based on FDR and fold-change thresholds

2.  Adds a 'regulation' column based on the sign of the fold-change

The resulting data frame is suitable for use with
[`plot_diff_logo`](https://prolfqua.github.io/prophosqua/reference/plot_diff_logo.md)
and
[`plot_seqlogo_with_diff`](https://prolfqua.github.io/prophosqua/reference/plot_seqlogo_with_diff.md)
which expect a 'regulation' column.

## See also

[`plot_diff_logo`](https://prolfqua.github.io/prophosqua/reference/plot_diff_logo.md),
[`plot_seqlogo_with_diff`](https://prolfqua.github.io/prophosqua/reference/plot_seqlogo_with_diff.md)

## Examples

``` r
# Example with mock data
data <- data.frame(
  contrast = rep("A_vs_B", 6),
  SequenceWindow = c("AAASAAAA", "BBBSBBB", "CCCSCCCC",
                     "DDDSDDDD", "EEESEEEE", "FFFSFFF"),
  FDR.site = c(0.01, 0.03, 0.08, 0.02, 0.15, 0.04),
  diff.site = c(1.2, -0.8, 0.5, -1.5, 0.3, 0.9),
  stringsAsFactors = FALSE
)

# Filter with default thresholds (FDR < 0.05, |FC| > 0.6)
sig_sites <- filter_significant_sites(data)
sig_sites$regulation
#> [1] "upregulated"   "downregulated" "downregulated" "upregulated"  

# Use with DPU column names
data_dpu <- data.frame(
  contrast = "A_vs_B",
  SequenceWindow = c("AAASAAAA", "BBBSBBB"),
  FDR_I = c(0.01, 0.03),
  diff_diff = c(1.2, -0.8)
)
sig_dpu <- filter_significant_sites(data_dpu, fdr_col = "FDR_I", diff_col = "diff_diff")
```
