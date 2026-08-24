# Compute PTM Usage the CorrectFirst Way

Corrects phosphosite abundances by their protein abundance
\*\*before\*\* modelling, then fits one linear model per site on the
corrected values and evaluates the contrasts on those models. This is
the mirror image of \[test_diff()\], which models site and protein
separately and subtracts the two fold changes afterwards.

## Usage

``` r
compute_cf_dea(phospho_dea_dir, protein_dea_dir, annot_file)
```

## Arguments

- phospho_dea_dir:

  Path to the phospho DEA output directory.

- protein_dea_dir:

  Path to the total-proteome DEA output directory.

- annot_file:

  Sample annotation defining the groups and the contrasts. Two layouts
  are accepted, see the contrast derivation below.

## Value

A list carrying the result table and everything a report needs to
describe it without refitting: \`results\`, \`ptm_data\`, \`ctr\`,
\`annot\`, \`contrasts\`, \`wide_data\`, \`wide_annotation\`,
\`model_counts\`, and the measurement and model counts the prose quotes.

## Details

The correction is a subtraction on the log2 scale, so it is a ratio on
the linear scale: a site keeps only the signal its protein does not
explain. Sites whose protein was not quantified in the same sample drop
out at the join and are absent from the result.

Site-contrast pairs whose estimate would rest on an imputed value in one
of the two groups are dropped from \`results\`, because a fold change
measured against an imputed level is not interpretable. \`model_counts\`
and \`n_before\` record that split so a report can state what was
dropped.

## See also

\[compute_dpa_dpu()\] for the correct-last alternative.

## Examples

``` r
# Needs two prolfquapp DEA output directories and an annotation; see
# tests/testthat/helper-dea_fixture.R for the synthetic set the tests use.
if (FALSE) { # \dontrun{
res <- compute_cf_dea(
  phospho_dea_dir = "DEA_20260814_WUphospho_vsn",
  protein_dea_dir = "DEA_20260814_WUtotal_vsn",
  annot_file = "phospho_dataset.tsv"
)
res$model_counts
} # }
```
