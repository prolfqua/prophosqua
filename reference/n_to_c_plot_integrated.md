# N to C for integrated results

N to C for integrated results

## Usage

``` r
n_to_c_plot_integrated(
  poi_matrix_min,
  protein_name,
  prot_length,
  contrast,
  thr_a = 0.05,
  thr_b = 0.2,
  color_protein = "yellow"
)
```

## Arguments

- POI_matrixMin:

  data.frame

## Examples

``` r
data(n_c_integrated_df)
n_c_integrated_df$imputation_status <- "observed"
n_to_c_plot_integrated(n_c_integrated_df, "A0A1I9LT44", 539, "WTFC")
```
