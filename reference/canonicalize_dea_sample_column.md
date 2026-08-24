# Canonicalize the Sample Column in Normalized DEA Data

Renames the sample column a DEA run declares to a fixed name, so that
data from runs using different schemas can be joined.

## Usage

``` r
canonicalize_dea_sample_column(data, yaml_file, canonical_name = "Name")
```

## Arguments

- data:

  Normalized DEA data

- yaml_file:

  Path to its \`lfqdata.yaml\` file

- canonical_name:

  Output column name

## Value

\`data\` with its declared sample column renamed

## Examples

``` r
yaml_file <- file.path(tempdir(), "lfqdata.yaml")
writeLines("sample_name: sampleName", yaml_file)
data <- data.frame(sampleName = c("a", "b"), abundance = c(1, 2))

canonicalize_dea_sample_column(data, yaml_file)
#>   Name abundance
#> 1    a         1
#> 2    b         2
```
