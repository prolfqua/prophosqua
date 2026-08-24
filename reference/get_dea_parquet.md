# Get Parquet File from DEA Directory

Get Parquet File from DEA Directory

## Usage

``` r
get_dea_parquet(dea_dir)
```

## Arguments

- dea_dir:

  Path to DEA output directory

## Value

Full path to the parquet file

## Examples

``` r
dea_dir <- file.path(tempdir(), "DEA_example")
results <- file.path(dea_dir, "Results_WU_example")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
file.create(file.path(results, "lfqdata_normalized.parquet"))
#> [1] TRUE

basename(get_dea_parquet(dea_dir))
#> [1] "lfqdata_normalized.parquet"
```
