# Plot protein and PTM expression

Plot protein and PTM expression

## Usage

``` r
n_to_c_expression(
  combined_site_prot_long,
  contrast_name,
  FDR_threshold = 0.05,
  fc_threshold = 0,
  impute_flag = "Imputed_Mean_moderated"
)
```

## Arguments

- combined_site_prot_long:

  data frame with combined site and protein data

- contrast_name:

  name of the contrast to plot

- FDR_threshold:

  FDR threshold for filtering significant sites (default 0.05)

- fc_threshold:

  Fold change threshold for filtering (default 0)

- impute_flag:

  Flag for imputed values (default "Imputed_Mean_moderated")

## Value

data frame with protein_Id, protein_length, contrast, data, and plot
