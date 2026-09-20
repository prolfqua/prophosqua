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
2. Import both schema 2.0.0 `AnnData.h5ad` files into `PTM_inputs.h5mu`, including the stored design, contrasts, parameters, and reference data.
3. Compute DPA/DPU and CorrectFirst from those paired inputs, then complete the enrichment branches. Every persisted handoff is MuData.
4. Assemble `PTM_results.h5mu`, render reports from it, and export Excel/RDS delivery files last. The `ptm-pipeline` workflow coordinates these steps.

Each completed stage has its own R6 type. A transition returns a new complete object; it does not add optional fields to an earlier object.

```r
library(prophosqua)

inputs <- import_ptm_h5mu("phospho/AnnData.h5ad", "total/AnnData.h5ad", "PTM_inputs.h5mu")
dpa_dpu <- inputs$build(DPA_DPU)
cf <- inputs$build(CF)
statistics <- dpa_dpu$build(PTM_statistics, cf = cf)
statistics$write_h5mu("PTM_statistics.h5mu")
restored <- read_ptm_h5mu("PTM_statistics.h5mu", PTM_statistics)
```

`enriched` and `total` retain both DEA experiments; `cf` contains corrected abundances. DPA statistics belong to `enriched`; DPU and CorrectFirst belong to `cf`, each with its own presence mask. Protein-to-site joins and outer-join result rows retain the existing R behavior.

The former single-site `compute_ptm_results_h5ad()` writer and `ptm.sh ptm_h5ad` command are replaced by `import_h5mu`, `ptm_h5mu`, `enrich_h5mu`, `assemble_h5mu`, and terminal `export_h5mu`. The existing in-memory computational APIs remain available.

## Vignettes

The package includes vignettes demonstrating the analysis workflow:

- **`Analysis_n_to_c.Rmd`** - N-to-C plots for PTM site visualization
- **`Analysis_seqlogo.Rmd`** - Sequence logo analysis
- **`Analysis_PTMSEA.Rmd`** - PTM-SEA analysis
- **`Analysis_KinaseLibrary.Rmd`** - Kinase activity inference from phosphoproteomics data
- **`Analysis_MEA.Rmd`** - Motif enrichment analysis visualization

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
