# Get the Declared Sample Column from a DEA YAML File

prolfquapp records the column holding the sample name in its analysis
configuration. Older runs call it \`Name\`, current ones \`sampleName\`,
so the name is read rather than assumed.

## Usage

``` r
get_sample_name_column(yaml_file)
```

## Arguments

- yaml_file:

  Path to a prolfqua \`lfqdata.yaml\` file

## Value

The sample column name

## Examples

``` r
yaml_file <- file.path(tempdir(), "lfqdata.yaml")
writeLines("sample_name: sampleName", yaml_file)

get_sample_name_column(yaml_file)
#> [1] "sampleName"
```
