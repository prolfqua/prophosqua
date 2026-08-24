# Get the Declared Sample Column from a DEA Directory

Get the Declared Sample Column from a DEA Directory

## Usage

``` r
get_dea_sample_name_column(dea_dir)
```

## Arguments

- dea_dir:

  Path to a DEA output directory

## Value

The sample column name

## Examples

``` r
dea_dir <- file.path(tempdir(), "DEA_example")
results <- file.path(dea_dir, "Results_WU_example")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
writeLines("sample_name: sampleName", file.path(results, "lfqdata.yaml"))

get_dea_sample_name_column(dea_dir)
#> [1] "sampleName"
```
