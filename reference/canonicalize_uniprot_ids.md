# Canonicalize UniProt Protein Identifiers

Converts FASTA-style identifiers such as \`sp\|P12345\|PROT_HUMAN\` to
their accession while leaving bare accessions unchanged. The mapping
must remain one-to-one after upstream decoy filtering, otherwise two
distinct input identifiers would be joined as one protein; a violation
is an error.

## Usage

``` r
canonicalize_uniprot_ids(data, id_col = "protein_Id")
```

## Arguments

- data:

  Data containing protein identifiers

- id_col:

  Protein identifier column

## Value

\`data\` with canonicalized protein identifiers

## Examples

``` r
data <- data.frame(
  protein_Id = c("sp|P12345|PROT_HUMAN", "Q67890"),
  diff = c(1.2, -0.4)
)

canonicalize_uniprot_ids(data)
#>   protein_Id diff
#> 1     P12345  1.2
#> 2     Q67890 -0.4
```
