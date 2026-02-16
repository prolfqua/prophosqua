# Example phosphoproteomics differential expression results

Combined differential expression results from phosphoproteomics analysis
of Atg16l1 knockout vs wild-type macrophages. Contains both site-level
(DPA) and protein-corrected (DPU) statistics for 4 contrasts.

## Usage

``` r
combined_test_diff_example
```

## Format

A tibble with 105,824 rows and 56 columns including:

- protein_Id:

  UniProt protein identifier

- gene_name.site:

  Gene symbol

- modAA:

  Modified amino acid (S, T, or Y)

- posInProtein:

  Position of modification in protein

- SequenceWindow:

  Flanking sequence around modification site

- contrast:

  Comparison name (e.g., KO_vs_WT)

- diff.site:

  Log2 fold change for site (DPA)

- statistic.site:

  t-statistic for site

- p.value.site:

  p-value for site

- FDR.site:

  FDR-adjusted p-value for site

- diff_diff:

  Protein-corrected log2 fold change (DPU)

- tstatistic_I:

  t-statistic for DPU

- pValue_I:

  p-value for DPU

- FDR_I:

  FDR-adjusted p-value for DPU

## Source

Maculins et al. (2020) eLife, doi:10.7554/elife.62320
