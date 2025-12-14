# N to C plot using ggplot2

N to C plot using ggplot2

## Usage

``` r
n_to_c_plot(
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

- protein_name:

  name of protein

- contrast:

  name of contrast

- POI_matrixMin:

  data.frame

- protLength:

  protein length

- thrA:

  significance threshold small default 0.05

- thrB:

  significance threshold small large default 0.20

## Examples

``` r
data(exampleN_C_dat)
# Prepare data for plotting
poi_matrix_min <- prepare_n_to_c_data(exampleN_C_dat)

n_to_c_plot(subset(poi_matrix_min, protein_Id == "A0A1I9LPZ1"), "A0A1I9LPZ1", 2160, "H1FC")

n_to_c_plot(subset(poi_matrix_min, protein_Id == "A0A178US29"), "A0A178US29", 806, "H1FC")
```
