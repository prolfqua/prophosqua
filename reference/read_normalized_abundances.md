# Read a DEA Parquet with its Sample Column Canonicalized

The sample column of a prolfquapp parquet is named by the run's YAML,
which sits beside the parquet. Reading the two together is the only way
to get a \`Name\` column that means the same thing across runs.

## Usage

``` r
read_normalized_abundances(parquet)
```

## Arguments

- parquet:

  Path to \`lfqdata_normalized.parquet\`.

## Value

The long table with its sample column renamed to \`Name\`.
