# Plot PTM occupation, use after adjusting for total proteome

Plot PTM occupation, use after adjusting for total proteome

## Usage

``` r
n_to_c_usage(
  data_combined_diff,
  contrast_name,
  FDR_threshold = 0.05,
  fc_threshold = 0,
  impute_flag = "Imputed_Mean_moderated",
  protein_Id = "protein_Id"
)
```

## Arguments

- data_combined_diff:

  data frame with combined differential analysis results

- contrast_name:

  name of the contrast to plot

- FDR_threshold:

  FDR threshold for filtering significant sites (default 0.05)

- fc_threshold:

  Fold change threshold for filtering (default 0)

- impute_flag:

  Flag for imputed values (default "Imputed_Mean_moderated")

- protein_Id:

  Column name for protein identifier (default "protein_Id")

## Value

data frame with protein_Id, protein_length, contrast, data, and plot
