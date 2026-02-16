# Summarize significant sites by group

Generates summary statistics for significant sites, counting up- and
down-regulated sites by specified grouping variables.

## Usage

``` r
summarize_significant_sites(data, group_cols = c("contrast"))
```

## Arguments

- data:

  Data frame with 'regulation' column (typically output from
  [`filter_significant_sites`](https://prolfqua.github.io/prophosqua/reference/filter_significant_sites.md))

- group_cols:

  Character vector of column names to group by. Default `c("contrast")`

## Value

A tibble with counts pivoted wide, showing upregulated and downregulated
counts per group

## See also

[`filter_significant_sites`](https://prolfqua.github.io/prophosqua/reference/filter_significant_sites.md)

## Examples

``` r
# Example with mock data
data <- data.frame(
  contrast = c(rep("A_vs_B", 4), rep("C_vs_D", 3)),
  regulation = c("upregulated", "upregulated", "downregulated", "upregulated",
                 "downregulated", "downregulated", "upregulated")
)

summarize_significant_sites(data)
#> # A tibble: 2 × 3
#>   contrast downregulated upregulated
#>   <chr>            <int>       <int>
#> 1 A_vs_B               1           3
#> 2 C_vs_D               2           1
```
