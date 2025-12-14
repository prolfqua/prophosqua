# Plot PTM usage across all contrasts in multi-panel layout

Plot PTM usage across all contrasts in multi-panel layout

## Usage

``` r
n_to_c_usage_multicontrast(
  data_combined_diff,
  FDR_threshold = 0.05,
  fc_threshold = 0,
  impute_flag = "Imputed_Mean_moderated",
  protein_Id = "protein_Id",
  max_plots = NULL,
  include_proteins = NULL
)
```

## Arguments

- data_combined_diff:

  data frame with combined differential analysis results

- FDR_threshold:

  FDR threshold for filtering significant sites (default 0.05)

- fc_threshold:

  Fold change threshold for filtering (default 0)

- impute_flag:

  Flag for imputed values (default "Imputed_Mean_moderated")

- protein_Id:

  Column name for protein identifier (default "protein_Id")

- max_plots:

  Maximum number of plots to generate (default NULL = no limit). If
  specified, only the first max_plots proteins will be plotted.

- include_proteins:

  Character vector of protein IDs to always include in the output,
  regardless of max_plots limit (default NULL). These proteins will be
  added even if they exceed the max_plots limit.

## Value

data frame with protein_Id, protein_length, n_contrasts, and multi-panel
plot

## Examples

``` r
# data(combined_diff_data)
# result <- n_to_c_usage_multicontrast(combined_diff_data,
#   max_plots = 50, include_proteins = c("Q64337"))
# print(result$plot[[1]])  # Display first protein's multi-contrast plot
```
