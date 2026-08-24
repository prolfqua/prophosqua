# Combine per contrast N to C panels into one figure

Drops the per panel titles so that the protein description is shown once
for the whole figure, and collects the per panel legends into a single
guide area. The contrast of a panel remains readable from its y axis
label.

## Usage

``` r
combine_n_to_c_panels(contrast_plots, protein_data, protein_name, prot_length)
```

## Arguments

- contrast_plots:

  list of per contrast plots, one per contrast

- protein_data:

  data of a single protein across all contrasts, requires the columns
  site and modAA

- protein_name:

  name of protein

- prot_length:

  protein length

## Value

patchwork object with one title, one subtitle and one legend
