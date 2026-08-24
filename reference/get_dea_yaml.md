# Get YAML Config from DEA Directory

Get YAML Config from DEA Directory

## Usage

``` r
get_dea_yaml(dea_dir)
```

## Arguments

- dea_dir:

  Path to DEA output directory

## Value

Full path to the yaml file

## Examples

``` r
dea_dir <- file.path(tempdir(), "DEA_example")
results <- file.path(dea_dir, "Results_WU_example")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
file.create(file.path(results, "lfqdata.yaml"))
#> [1] TRUE

basename(get_dea_yaml(dea_dir))
#> [1] "lfqdata.yaml"
```
