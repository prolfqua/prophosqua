# Derive the imputation status of each row

Marks a row as imputed when any of the given estimate_type columns holds
impute_flag. Columns that the input does not carry are ignored, which
lets the same plots serve analyses that have no protein level estimate
(CF).

## Usage

``` r
.imputation_status(data, impute_flag, columns)
```

## Arguments

- data:

  data frame with zero or more estimate_type columns

- impute_flag:

  value of estimate_type marking an imputed estimate

- columns:

  estimate_type columns to consider

## Value

character vector of "imputed" / "observed", one element per row
