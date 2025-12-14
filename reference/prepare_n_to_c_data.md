# Prepare data for N-to-C plotting

Prepare data for N-to-C plotting

## Usage

``` r
prepare_n_to_c_data(poi_matrix_min, model_site = "model_site")
```

## Arguments

- poi_matrix_min:

  data.frame with phosphorylation data

## Value

data.frame with prepared data for plotting

## Examples

``` r
# example code
data(exampleN_C_dat)
poi_matrix_min <- prepare_n_to_c_data(exampleN_C_dat)
```
