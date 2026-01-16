# MetaboAnalystR LC-MS Polar Metabolomics Pipeline

This project provides a fully reproducible, automated pipeline for LC-MS/MS polar metabolomics data analysis using [MetaboAnalystR](https://github.com/xia-lab/MetaboAnalystR) and the [`targets`](https://docs.ropensci.org/targets/) workflow manager. The pipeline covers raw spectra processing, data cleaning, normalization, statistical analysis, compound annotation, and pathway interpretation, with dynamic Quarto reports at each stage.

## Features
- Modular, declarative pipeline with `targets`
- Full MetaboAnalystR workflow for polar metabolites
- Interactive Quarto reports (PCA, volcano, pathway maps)
- Example data and full reproducibility with `renv`
- GitHub-ready structure

## Quick Start
1. Clone this repo and open in RStudio or VS Code.
2. Run `renv::restore()` to install dependencies.
3. Place your raw LC-MS files in `data/` (mzML/mzXML/netCDF).
4. Run the pipeline: `targets::tar_make()`
5. View reports in `reports/`.

See `packages.R` for required packages and `R/` for modular pipeline code.

## Example Data
Example files are provided in `data/` for demonstration.

## License
MIT
