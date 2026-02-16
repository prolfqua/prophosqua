# Explode multisites into individual rows

Each phosphosite goes into its own row when multiple sites are present.

## Usage

``` r
explode_multisites(combined_site)
```

## Arguments

- combined_site:

  Data frame with PhosSites column (semicolon-delimited)

## Value

Data frame with one row per site, with modAA and posInProtein columns
