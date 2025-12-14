# Download PTMsigDB signatures

Downloads and caches PTMsigDB GMT files from the Broad Institute
ssGSEA2.0 repository. PTMsigDB v2.0.0 uses flanking sequences (15 amino
acids centered on the phosphosite), making the analysis
species-invariant within vertebrates.

## Usage

``` r
download_ptmsigdb(
  species = "mouse",
  version = "v2.0.0",
  output_dir = ".",
  force_download = FALSE
)
```

## Arguments

- species:

  Character. Species for signatures ("mouse" or "human"). Default
  "mouse".

- version:

  Character. PTMsigDB version. Default "v2.0.0".

- output_dir:

  Character. Directory for cached files. Default current directory.

- force_download:

  Logical. Force re-download even if cached. Default FALSE.

## Value

Character path to the cached GMT file.

## References

Krug et al. (2019) A Curated Resource for Phosphosite-specific Signature
Analysis. Mol Cell Proteomics. doi:10.1074/mcp.TIR118.000943

## Examples

``` r
if (FALSE) { # \dontrun{
gmt_path <- download_ptmsigdb(species = "mouse", output_dir = "ptmsea_cache")
pathways <- fgsea::gmtPathways(gmt_path)
} # }
```
