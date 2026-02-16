# Test if differences of differences are significant (internal)

Test if differences of differences are significant (internal)

## Usage

``` r
.test_diff_diff(
  dataframe_a,
  dataframe_b,
  by,
  diff = c("diff"),
  std_err = c("std.error"),
  df = c("df"),
  suffix_a = ".site",
  suffix_b = ".protein"
)
```

## Arguments

- dataframe_a:

  First data frame (e.g., site-level results)

- dataframe_b:

  Second data frame (e.g., protein-level results)

- by:

  Columns to join by

- diff:

  Column name for difference values

- std_err:

  Column name for standard error values

- df:

  Column name for degrees of freedom

- suffix_a:

  Suffix for columns from dataframe_a

- suffix_b:

  Suffix for columns from dataframe_b

## Value

Data frame with diff_diff test results
