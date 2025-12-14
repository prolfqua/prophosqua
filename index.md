# prophosqua

**Integration of phosphoproteome and total proteome data for
comprehensive PTM analysis**

The `prophosqua` package provides tools for integrating and analyzing
post-translational modification (PTM) data with total proteome
measurements. It enables researchers to distinguish between changes in
protein abundance and changes in modification site usage, providing
deeper insights into cellular signaling and regulation.

## Overview

`prophosqua` builds upon the `prolfqua` and `prolfquapp` packages to
provide:

- **Differential PTM Feature Expression (DPE)** - Analysis of
  modification site abundance changes
- **Differential PTM Feature Usage (DPU)** - Analysis of modification
  site usage relative to total protein
- **N-to-C plots** - Visualization of phosphorylation sites along
  protein backbones
- **Sequence logo analysis** - Identification of kinase recognition
  motifs
- **Statistical integration** - Combined analysis of phosphoproteome and
  total proteome data

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
devtools::install_github('prolfqua/prolfquapp', dependencies = TRUE)
```

For detailed installation instructions and system requirements, see: -
[prolfqua GitHub repository](https://github.com/prolfqua/prolfqua) -
[prolfquapp GitHub repository](https://github.com/prolfqua/prolfquapp)

### Install prophosqua

``` r

library(devtools)
devtools::install_github('prolfqua/prophosqua', dependencies = TRUE)
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

## Vignettes and Documentation

The package includes comprehensive vignettes that demonstrate the
complete analysis workflow:

### Rendering Vignettes

The following code renders the vignettes and copies them to the example
directory:

``` r

# Render the differential expression analysis vignette

wd <- "inst/PTM_analysis_example"


rmarkdown::render("vignettes/QCReport.qmd",
params=list(wd=wd,
psm="ptm_example-main/qc_example_data/QCmini/psm.tsv",
fasta="ptm_example-main/qc_example_data/fgcz_3702_UP000006548_AraUniprot_1spg_d_20231024.fasta",
workunit="exampleQC",
projectid="ABCD1234"))
file.copy("vignettes/QCReport.html", file.path(wd,"QCReport.html"))

# Render the PTM integration and visualization vignette  
rmarkdown::render("vignettes/Supplementary_Material_v2.Rmd",
params=list(wd = wd, dodea = TRUE), output_format=bookdown::pdf_document2())
file.copy("vignettes/Supplementary_Material_v2.pdf", file.path(wd,"Supplementary_Material_v2.pdf"))

rmarkdown::render("vignettes/Supplementary_Material_v2.Rmd",
params=list(wd = wd, dodea = FALSE), output_format=bookdown::html_document2())
file.copy("vignettes/Supplementary_Material_v2.html", file.path(wd,"Supplementary_Material_v2.html"))
```

This code: 1. **Renders R Markdown files** to HTML using
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
2. **Copies the HTML files** to the `inst/PTM_analysis_example/`
directory for easy access 3. **Creates interactive reports** with code
folding and navigation

### Available Vignettes

- **`Run_DEA_prolfquapp.Rmd`** - Complete workflow for differential
  expression analysis
- **`Visualize_PTM_features.Rmd`** - Integration and visualization of
  PTM data
- **Interactive reports** - Available in the
  `inst/PTM_analysis_example/` directory

## Related Publications

### Core Dependencies

- **prolfqua**: [Proteomics data analysis in
  R](https://www.sciencedirect.com/science/article/pii/S1535947622002857) -
  MSstatsPTM publication
- **prolfquapp**: [Proteomics analysis workflow
  package](https://www.sciencedirect.com/science/article/pii/S1535947623002190) -
  msqrob2PTM publication

### Related Work

- [MSstatsPTM](https://www.sciencedirect.com/science/article/pii/S1535947622002857) -
  Statistical analysis of PTM data
- [msqrob2PTM](https://www.sciencedirect.com/science/article/pii/S1535947623002190) -
  Robust statistical methods for PTM analysis

## GitHub Repositories

- [prolfqua](https://github.com/prolfqua/prolfqua) - Core proteomics
  analysis package
- [prolfquapp](https://github.com/prolfqua/prolfquapp) - Proteomics
  analysis workflow package
- [prophosqua](https://github.com/prolfqua/prophosqua) - PTM integration
  and analysis package

## Poster

[Prophosqua Poster from SIBDays
2024](https://prolfqua.github.io/prophosqua/img/SIBDays24_Poster_Prophosqua_Grossmann.pdf)

## Citation

If you use this package in your research, please cite:

``` bibtex
@article{prophosqua2024,
  title={prophosqua: Integration of phosphoproteome and total proteome data},
  author={Functional Genomics Center Zurich},
  journal={Bioinformatics},
  year={2024}
}
```

## Contributing

Contributions are welcome! Please visit our [GitHub
repository](https://github.com/prolfqua/prophosqua) for: - Issue
reporting - Feature requests - Code contributions - Documentation
improvements

## License

This package is released under the [MIT
License](https://prolfqua.github.io/prophosqua/LICENSE).
