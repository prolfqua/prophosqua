# Get File from DEA Directory

Generic helper to find files within a DEA directory's \`Results_WU\_\*\`
subfolder.

## Usage

``` r
get_dea_file(dea_dir, filename, description = "file")
```

## Arguments

- dea_dir:

  Path to DEA output directory

- filename:

  Filename to find (e.g., "lfqdata_normalized.parquet")

- description:

  Description for error message (e.g., "parquet file")

## Value

Full path to the file

## Examples

``` r
dea_dir <- file.path(tempdir(), "DEA_example")
results <- file.path(dea_dir, "Results_WU_example")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
file.create(file.path(results, "lfqdata.yaml"))
#> [1] TRUE

basename(get_dea_file(dea_dir, "lfqdata.yaml", "yaml file"))
#> [1] "lfqdata.yaml"
```
