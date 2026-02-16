# prophosqua

**Integration of phosphoproteome and total proteome data for
comprehensive PTM analysis**

The `prophosqua` package provides tools for integrating and analyzing
post-translational modification (PTM) data with total proteome
measurements. It enables researchers to distinguish between changes in
protein abundance and changes in modification site usage, providing
deeper insights into cellular signaling and regulation.

## Overview

The integrated PTM analysis is carried out using `prophosqua`, which
implements three complementary statistical approaches:

- **DPA** (Differential PTM Abundance): tests for changes in PTM-site
  intensity between conditions.
- **DPU** (Differential PTM Usage): tests whether the ratio of PTM-site
  to total protein intensity changes, identifying regulation independent
  of protein abundance.
- **CorrectFirst**: applies protein-level correction before testing PTM
  sites, an alternative approach to distinguish PTM-specific regulation
  from protein expression changes.

Building on the [prolfqua](https://github.com/prolfqua/prolfqua) and
[prolfquapp](https://github.com/prolfqua/prolfquapp) packages,
`prophosqua` additionally provides:

- **N-to-C plots** - Visualization of phosphorylation sites along
  protein backbones
- **Sequence logo analysis** - Identification of kinase recognition
  motifs
- **PTM-SEA** - Post-translational modification set enrichment analysis
- **Kinase activity inference** - Kinase library-based analysis from
  phosphoproteomics data
- **Motif enrichment analysis** - MEA visualization
- **Enrichment visualization** - Dot plots, heatmaps, and volcano plots
  for enrichment results

## Installation

### Prerequisites

This package depends on several other packages that should be installed
first:

``` r

# Install prolfqua (core proteomics analysis package)
library(devtools)
devtools::install_github('protviz/prozor', dependencies = TRUE)
devtools::install_github('prolfqua/prolfqua', dependencies = TRUE)
# Install prolfquapp (proteomics analysis workflow package)
devtools::install_github('prolfqua/prolfquapp', dependencies = TRUE, build_vignettes=TRUE)
```

For detailed installation instructions and system requirements, see: -
[prolfqua GitHub repository](https://github.com/prolfqua/prolfqua) -
[prolfquapp GitHub repository](https://github.com/prolfqua/prolfquapp)

### Install prophosqua

``` r

library(devtools)
devtools::install_github('prolfqua/prophosqua', dependencies = TRUE, build_vignettes=TRUE)
```

## Usage

### Basic Workflow

1.  **Run differential expression analysis** using `prolfquapp` for both
    total proteome and PTM data
2.  **Integrate the results** using
    [`prophosqua::test_diff()`](https://prolfqua.github.io/prophosqua/reference/test_diff.md)
3.  **Visualize the data** with N-to-C plots and sequence logos
4.  **Generate comprehensive reports** with statistical summaries

### Example

``` r

library(prophosqua)

# Load your data
tot_res <- load_and_preprocess_data(total_proteome_file, required_cols)
phospho_res <- load_and_preprocess_data(ptm_file, required_cols)

# Integrate analysis
combined_test_diff <- test_diff(phospho_res, tot_res, join_column = join_column)

# Create visualizations
plot_data <- n_to_c_usage(combined_test_diff, contrast_name, FDR_threshold = 0.05)
```

## Vignettes

The package includes vignettes demonstrating the analysis workflow:

- **`MiMBIntegratedPTM.Rmd`** - Integrated analysis of PTM and total
  proteome (DPA, DPU, CorrectFirst)
- **`Analysis_n_to_c.Rmd`** - N-to-C plots for PTM site visualization
- **`Analysis_seqlogo.Rmd`** - Sequence logo analysis
- **`Analysis_PTMSEA.Rmd`** - PTM-SEA analysis
- **`Analysis_KinaseLibrary.Rmd`** - Kinase activity inference from
  phosphoproteomics data
- **`Analysis_MEA.Rmd`** - Motif enrichment analysis visualization
- **`QCReport.qmd`** - FragPipe TMT quality control report

## Citation

If you use this package in your research, please cite:

> Wolski W, Dittmann A, Panse C, Kunz L, Grossmann J. “Integrated
> Analysis of Post-Translational Modifications and Total Proteome:
> Methods for Distinguishing Abundance from Usage Changes.” *Methods in
> Molecular Biology*, 2025 (submitted).

> Grossmann J, Wolski W. “prolfqua/prophosqua: 0.1.0.” Zenodo, 2025.
> DOI:
> [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)

## Related Packages

- [prolfqua](https://github.com/prolfqua/prolfqua) - Core proteomics
  analysis package
- [prolfquapp](https://github.com/prolfqua/prolfquapp) - Proteomics
  analysis workflow package

## Building and Deploying Documentation

### Build pkgdown site locally

``` r

pkgdown::build_site()
```

### Deploy to GitHub Pages

``` bash
# Using ghp-import (via uv)
uvx ghp-import -n -p -f docs
```

Then enable GitHub Pages in repo settings (Settings -\> Pages -\>
Source: `gh-pages` branch).

Site: <https://prolfqua.github.io/prophosqua>

## Contributing

Contributions are welcome! Please visit our [GitHub
repository](https://github.com/prolfqua/prophosqua) for: - Issue
reporting - Feature requests - Code contributions

## License

This package is released under the [MIT
License](https://prolfqua.github.io/prophosqua/LICENSE).
