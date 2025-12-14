# Get default list of kinase motifs

Returns a data.frame of common kinase motifs defined by regular
expressions. Patterns are designed to match 15-mer sequences centered at
the phosphosite.

## Usage

``` r
get_kinase_motifs(extra_motifs = NULL)
```

## Arguments

- extra_motifs:

  Optional data.frame to append to the default list. Must have columns:
  kinase_group, motif_name, pattern, description.

## Value

data.frame with columns: kinase_group, motif_name, pattern, description

## Note

This list is a compiled "starter set" of well-known consensus motifs
from literature (e.g. classical definitions for PKA, PKC, CK2, CDK,
etc.). It does NOT contain the entire proprietary PhosphoSitePlus
database. Users with access to curated motif databases should load them
using \`read_kinase_motifs_csv\` or pass them as \`extra_motifs\`.

## Examples

``` r
motifs <- get_kinase_motifs()
head(motifs)
#>   kinase_group        motif_name          pattern
#> 1   Basophilic PKA_AKT_consensus   R.R..[ST].....
#> 2   Basophilic       PKA_Classic         R.[ST]..
#> 3   Basophilic       PKC_General      ..[ST].[RK]
#> 4   Basophilic          Akt_RRxS       RR.[ST]...
#> 5   Basophilic             CAMK2       R..[ST]...
#> 6   Basophilic Substrate_RK_rich [RK][RK].[ST]...
#>                                          description
#> 1 R-x-R-x-x-S/T (PKA, AKT, etc. - strong basophilic)
#> 2                              R-x-S/T (PKA classic)
#> 3                            S/T-x-R/K (PKC general)
#> 4                         R-R-x-S/T (Akt preference)
#> 5                                 R-x-x-S/T (CaMKII)
#> 6            Basic residues upstream (R/K-R/K-x-S/T)
```
