# Prepare data for PTM-SEA analysis

Converts phosphosite data to a list of named vectors suitable for GSEA.
Each contrast produces one named vector where names are site IDs in
PTMsigDB format (SEQUENCE-p) and values are t-statistics (or other
ranking metric).

## Usage

``` r
ptmsea_data_prep(
  data,
  stat_column,
  seq_window_col = "SequenceWindow",
  contrast_col = "contrast",
  trim_to = c("11", "13", "15")
)
```

## Arguments

- data:

  data.frame containing phosphosite data with statistics per contrast

- stat_column:

  Character. Column name for ranking statistic. Use "statistic.site" for
  DPA or "tstatistic_I" for DPU.

- seq_window_col:

  Character. Column name for sequence windows. Default "SequenceWindow".

- contrast_col:

  Character. Column name for contrasts. Default "contrast".

- trim_to:

  Character. Trim flanking sequences to this width. Options: "11"
  (default), "13", or "15" (no trimming). Trimming increases overlap
  with PTMsigDB.

## Value

List with two elements:

- ranks: Named list of named numeric vectors (one per contrast)

- dropped: Named list of dropped duplicates (one per contrast with
  duplicates)

## Examples

``` r
# Prepare mock data
data <- data.frame(
  contrast = c("A", "A", "B", "B"),
  SequenceWindow = c("AAAAAAASAAAAAAA", "BBBBBBBSBBBBBBB", "AAAAAAASAAAAAAA", "CCCCCCCSCCCCCCC"),
  statistic.site = c(2.5, 1.5, -2.0, -1.0),
  stringsAsFactors = FALSE
)

# Prepare DPA data using t-statistic for ranking
prep_dpa <- ptmsea_data_prep(data, stat_column = "statistic.site")
prep_dpa$ranks # ranked lists per contrast
#> $A
#> AAAAASAAAAA-p BBBBBSBBBBB-p 
#>           2.5           1.5 
#> 
#> $B
#> CCCCCSCCCCC-p AAAAASAAAAA-p 
#>            -1            -2 
#> 
prep_dpa$dropped # duplicates that were dropped
#> list()

# Trim to 13-mer to potentially increase PTMsigDB overlap
prep_dpa_13 <- ptmsea_data_prep(data, stat_column = "statistic.site", trim_to = "13")
```
