# Convert PTMsigDB pathways to TERM2GENE format

Converts pathway list from fgsea::gmtPathways() to data.frame format
required by clusterProfiler::enricher.

## Usage

``` r
ptmsigdb_to_term2gene(pathways)
```

## Arguments

- pathways:

  Named list of pathways (from fgsea::gmtPathways)

## Value

data.frame with columns: term, gene
