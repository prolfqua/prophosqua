# Compute the Enrichment Analyses from the Bundled Example Data

The three enrichment reports can be rendered without a pipeline run,
against the example data the package ships. These helpers assemble that
input and then call the same compute function the pipeline calls, so a
standalone render exercises the real code path rather than a parallel
one.

## Usage

``` r
compute_ptmsea_example()

compute_kinaselib_gsea_example()

compute_mea_example()
```
