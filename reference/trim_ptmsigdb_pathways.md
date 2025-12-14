# Trim PTMsigDB pathways to shorter flanking sequences

Trims all site IDs in PTMsigDB pathways to specified width. Use this
when trimming your data with ptmsea_data_prep(trim_to=...) to ensure
both sides match.

## Usage

``` r
trim_ptmsigdb_pathways(pathways, trim_to = c("11", "13", "15"))
```

## Arguments

- pathways:

  Named list of pathways from fgsea::gmtPathways().

- trim_to:

  Character. Target width: "11" (default), "13", or "15" (no trim).

## Value

Named list of pathways with trimmed site IDs.

## Examples

``` r
# Mock pathways
pathways <- list(
  PathwayA = c("AAAAAAASAAAAAAA-p", "BBBBBBBSBBBBBBB-p"),
  PathwayB = c("CCCCCCCSCCCCCCC-p")
)

# Trim to 11-mer
pathways_11 <- trim_ptmsigdb_pathways(pathways, trim_to = "11")
```
