# Get Excel File from DEA Directory

Finds the Excel results file within a DEA directory's \`Results_WU\_\*\`
subfolder, preferring the \`DE\_\`-prefixed workbook that prolfquapp
writes. Falls back to searching \`dea_dir\` itself when no such
subfolder exists.

## Usage

``` r
get_dea_xlsx(dea_dir)
```

## Arguments

- dea_dir:

  Path to DEA output directory

## Value

Full path to the Excel file

## Examples

``` r
dea_dir <- file.path(tempdir(), "DEA_example")
results <- file.path(dea_dir, "Results_WU_example")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
file.create(file.path(results, c("DE_example.xlsx", "other.xlsx")))
#> [1] TRUE TRUE

basename(get_dea_xlsx(dea_dir))
#> Using: DE_example.xlsx
#> [1] "DE_example.xlsx"
```
