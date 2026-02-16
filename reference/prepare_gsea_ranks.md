# Prepare ranked lists for GSEA analysis

Creates named numeric vectors suitable for clusterProfiler::GSEA from
phosphoproteomics data with multiple contrasts. This is a
general-purpose function that works with any GSEA analysis (Kinase
Library, motif enrichment, etc.). For PTMsigDB-specific preparation with
trimming and "-p" suffix, see
[`ptmsea_data_prep`](https://prolfqua.github.io/prophosqua/reference/ptmsea_data_prep.md).

## Usage

``` r
prepare_gsea_ranks(
  data,
  stat_col,
  seq_col = "SequenceWindow",
  contrast_col = "contrast",
  to_uppercase = TRUE,
  add_suffix = NULL
)
```

## Arguments

- data:

  Data frame with phosphosite data containing statistics per contrast

- stat_col:

  Character. Column name for ranking statistic (e.g., "statistic.site"
  for DPA or "tstatistic_I" for DPU)

- seq_col:

  Character. Column name for sequence windows. Default "SequenceWindow"

- contrast_col:

  Character. Column name for contrasts. Default "contrast"

- to_uppercase:

  Logical. Convert sequences to uppercase. Default TRUE

- add_suffix:

  Character or NULL. Suffix to add to sequences (e.g., "-p" for PTMsigDB
  format). Default NULL (no suffix)

## Value

Named list of named numeric vectors, one per contrast. Each vector
contains statistics named by sequence, sorted in descending order.

## Details

For each contrast, the function:

1.  Filters to the specified contrast

2.  Processes sequences (uppercase, trim whitespace)

3.  Removes rows with NA statistics

4.  Keeps first occurrence of duplicate sequences

5.  Sorts by statistic (descending) for GSEA

## See also

[`ptmsea_data_prep`](https://prolfqua.github.io/prophosqua/reference/ptmsea_data_prep.md)
for PTMsigDB-specific preparation with trimming

## Examples

``` r
# Example with mock data
data <- data.frame(
  contrast = rep(c("A_vs_B", "C_vs_D"), each = 5),
  SequenceWindow = c("AAASAAAA", "BBBSBBB", "CCCSCCCC", "DDDSDDDD", "EEESEEEE",
                     "AAASAAAA", "FFFSFFF", "GGGSGGGG", "HHHSHHHH", "IIISIII"),
  statistic.site = c(2.5, 1.8, -0.5, -1.2, 0.3, 1.9, 2.1, -0.8, 0.1, -1.5),
  stringsAsFactors = FALSE
)

# Prepare ranks for Kinase Library GSEA (no suffix needed)
ranks <- prepare_gsea_ranks(data, stat_col = "statistic.site")
names(ranks)  # contrast names
#> [1] "A_vs_B" "C_vs_D"
head(ranks[[1]])  # ranked sequences for first contrast
#> AAASAAAA  BBBSBBB EEESEEEE CCCSCCCC DDDSDDDD 
#>      2.5      1.8      0.3     -0.5     -1.2 

# Prepare ranks for PTMsigDB (with "-p" suffix)
ranks_ptm <- prepare_gsea_ranks(data, stat_col = "statistic.site", add_suffix = "-p")
```
