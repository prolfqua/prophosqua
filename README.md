[![DOI](https://zenodo.org/badge/784111954.svg)](https://doi.org/10.5281/zenodo.15845272)
[![altdoc](https://img.shields.io/badge/docs-altdoc-blue)](https://prolfqua.github.io/prophosqua/)

# prophosqua

**Integration of phosphoproteome and total proteome data for comprehensive PTM analysis**

The `prophosqua` package provides tools for integrating and analyzing post-translational modification (PTM) data with total proteome measurements. It enables researchers to distinguish between changes in protein abundance and changes in modification site usage, providing deeper insights into cellular signaling and regulation.

## Overview

The integrated PTM analysis is carried out using `prophosqua`, which implements three complementary statistical approaches:

- **DPA** (Differential PTM Abundance): tests for changes in PTM-site intensity between conditions.
- **DPU** (Differential PTM Usage): tests whether the ratio of PTM-site to total protein intensity changes, identifying regulation independent of protein abundance.
- **CorrectFirst**: applies protein-level correction before testing PTM sites, an alternative approach to distinguish PTM-specific regulation from protein expression changes.

Building on the [prolfqua](https://github.com/prolfqua/prolfqua) and [prolfquapp](https://github.com/prolfqua/prolfquapp) packages, `prophosqua` additionally provides:

- **N-to-C plots** - Visualization of phosphorylation sites along protein backbones
- **Sequence logo analysis** - Identification of kinase recognition motifs
- **PTM-SEA** - Post-translational modification set enrichment analysis
- **Kinase activity inference** - Kinase library-based analysis from phosphoproteomics data
- **Motif enrichment analysis** - MEA visualization
- **Enrichment visualization** - Dot plots, heatmaps, and volcano plots for enrichment results

## Installation

### Prerequisites

This package depends on several other packages that should be installed first:

```r
# Install prolfqua (core proteomics analysis package)
library(devtools)
devtools::install_github('protviz/prozor', dependencies = TRUE)
devtools::install_github('prolfqua/prolfqua', dependencies = TRUE)
# Install prolfquapp (proteomics analysis workflow package)
devtools::install_github('prolfqua/prolfquapp', dependencies = TRUE, build_vignettes=TRUE)
```

For detailed installation instructions and system requirements, see:
- [prolfqua GitHub repository](https://github.com/prolfqua/prolfqua)
- [prolfquapp GitHub repository](https://github.com/prolfqua/prolfquapp)

### Install prophosqua

```r
library(devtools)
devtools::install_github('prolfqua/prophosqua', dependencies = TRUE, build_vignettes=TRUE)
```

## Usage

### Basic Workflow

1. Run DEA with `prolfquapp` for enriched sites and total protein.
2. Import both `AnnData.h5ad` files, the analysis parameters and the reference data into `PTM_inputs.h5mu`.
3. Compute DPA, DPU and CorrectFirst into `PTM_statistics.h5mu`, then the enrichment of each analysis into files beside it: PTM-SEA, Kinase GSEA and MEA as protsea documents (`.json.gz`), the kinase-library preparations as gzipped JSON.
4. Assemble `PTM_results.h5mu`, which names these files but holds no enrichment, render the reports from it, and export the delivery workbook last. The `ptm-pipeline` workflow coordinates these steps.

Each persisted stage has its own R6 type, and MuData is the persistence boundary:

```r
library(prophosqua)

inputs <- import_ptm_h5mu("phospho/AnnData.h5ad", "total/AnnData.h5ad", "PTM_inputs.h5mu")
statistics <- PTM_statistics$new(inputs)
statistics$write_h5mu("PTM_statistics.h5mu")
restored <- read_ptm_h5mu("PTM_statistics.h5mu", PTM_statistics)
```

`PTM_statistics.h5mu` keeps both DEA experiments as the `enriched` and `total` modalities, with the DPA and DPU results as `varm` frames of `enriched`; CorrectFirst and its imputed variants form the `enriched_CF` modality. The pipeline runs these steps as `ptm.sh import_h5mu`, `ptm_h5mu`, `enrich`, `assemble_h5mu`, `render` and `export_h5mu`.

For a pair of DEA folders outside the pipeline, `compute_dpa_dpu()` and `compute_cf_dea()` return the DPA/DPU and CorrectFirst results in memory.

## Vignettes

The vignettes are the two Quarto reports the PTM pipeline renders from MuData:

- **`ptm_statistics.qmd`** - DPA, DPU and CorrectFirst results from `PTM_statistics.h5mu`
- **`ptm_enrichment.qmd`** - PTM-SEA, Kinase GSEA and MEA for one analysis from `PTM_results.h5mu`

The MiMB manuscript source is kept outside `vignettes/` and is rendered by
`inst/MiMB_build/Snakefile`:

- **`manuscript/_MiMBIntegratedPTM.Rmd`** - Integrated analysis of PTM and total proteome
  (DPA, DPU, CorrectFirst); the MiMB manuscript

The separate quality-control source is a Quarto vignette:

- **`vignettes/_QCReport.qmd`** - FragPipe TMT quality control report

## Citation

If you use this package in your research, please cite:

> Wolski W, Dittmann A, Panse C, Kunz L, Grossmann J.
> "Integrated Analysis of Post-Translational Modifications and Total Proteome: Methods for Distinguishing Abundance from Usage Changes."
> *Methods in Molecular Biology*, 2025 (submitted).

> Grossmann J, Wolski W.
> "prolfqua/prophosqua: 0.1.0." Zenodo, 2025.
> DOI: [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)

## Related Packages

- [prolfqua](https://github.com/prolfqua/prolfqua) - Core proteomics analysis package
- [prolfquapp](https://github.com/prolfqua/prolfquapp) - Proteomics analysis workflow package

## Building and Deploying Documentation
- https://deepwiki.com/wolski/ptm-pipeline - AI attempt for documentation

### Build the altdoc site locally

```r
altdoc::render_docs(freeze = FALSE)
```

### Deploy to GitHub Pages

The `altdoc` GitHub Actions workflow renders the site and deploys `docs/` to
the `gh-pages` branch after every push to `main`.

Site: https://prolfqua.github.io/prophosqua

## Contributing

Contributions are welcome! Please visit our [GitHub repository](https://github.com/prolfqua/prophosqua) for:
- Issue reporting
- Feature requests
- Code contributions

## License

This package is released under the [MIT License](LICENSE).
