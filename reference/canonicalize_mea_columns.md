# Bring MEA Output onto the Shared Enrichment Column Names

The kinase-library tool writes \`Kinase\`, \`p-value\` and a \`Subs
fraction\` field of the form \`85/471\`. The enrichment reports expect
\`kinase\`, \`pvalue\`, and the fraction split into its leading-edge
count and set size.

## Usage

``` r
canonicalize_mea_columns(mea_results)
```

## Arguments

- mea_results:

  Concatenated MEA results.

## Value

The same table, renamed and with the fraction split.
