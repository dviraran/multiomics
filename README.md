# Multi-Omics Analysis Pipelines

A collection of reproducible R pipelines for analyzing omics data, built with the [`{targets}`](https://docs.ropensci.org/targets/) framework.

## Pipelines

| Pipeline | Description | Key Methods |
|----------|-------------|-------------|
| [rnaseq_pipeline](rnaseq_pipeline/) | RNA-seq differential expression analysis | DESeq2, fgsea, GO/KEGG enrichment |
| [proteomics_pipeline](proteomics_pipeline/) | Mass spectrometry proteomics (LFQ/DIA) | limma, VSN normalization, MNAR imputation |
| [metabolomics_pipeline](metabolomics_pipeline/) | LC-MS/GC-MS metabolomics | PQN normalization, drift correction, MetaboAnalyst-style QC |
| [multiomics_pipeline](multiomics_pipeline/) | Multi-omics data integration | MOFA2, DIABLO (mixOmics), SNF |

## Features

- **Reproducible**: Full dependency tracking with `{targets}` - only re-run what changed
- **Configurable**: YAML-based configuration for all parameters
- **Modular**: Functions organized in `R/` directory for easy customization
- **Comprehensive QC**: Built-in quality control checks and visualizations
- **Publication-ready**: Automated HTML reports with R Markdown
- **AI Commentary** (RNA-seq/Proteomics): Optional LLM-generated figure interpretations
- **Auto-install packages**: Prompts to install missing packages on first run

## Quick Start

### Option 1: Run from pipeline directory

```bash
cd <pipeline_name>

# 1. Edit config.yml with your settings
# 2. Place data in data/ directory
# 3. Run the pipeline:
```

```r
library(targets)
tar_make()        # Run pipeline
tar_visnetwork()  # Visualize dependencies
tar_read(result)  # Access any result
```

### Option 2: Run with custom config from anywhere

Use the `pipeline_runner.R` helper to run pipelines with custom config files from any directory:

```r
source("path/to/multiomics/pipeline_runner.R")

# Run with custom config file
run_rnaseq_pipeline("my_project/rnaseq_config.yml")
run_proteomics_pipeline("my_project/proteomics_config.yml")
run_metabolomics_pipeline("my_project/metabolomics_config.yml")
run_multiomics_pipeline("my_project/multiomics_config.yml")

# Or use the generic function
run_omics_pipeline("rnaseq", "path/to/config.yml")

# Clean cache and rerun from scratch
run_rnaseq_pipeline("config.yml", clean = TRUE)
```

### Option 3: Run all pipelines in sequence

```bash
cd examples/scripts
Rscript run_all_pipelines.R --clean
```

## Directory Structure

Each pipeline contains:

```
<pipeline>/
├── _targets.R           # Pipeline definition
├── config.yml           # Configuration
├── R/                   # Modular functions
├── data/                # Input data (+ examples)
├── outputs/             # Results (tables, plots)
└── reports/             # R Markdown templates
```

## Requirements

### Automatic Package Installation

Each pipeline will automatically check for missing packages on first run and prompt you to install them. Just run `tar_make()` and follow the prompts.

### Manual Installation (optional)

If you prefer to install packages manually:

```r
# Core packages
install.packages(c("targets", "tarchetypes", "tidyverse", "yaml"))

# Bioconductor packages (varies by pipeline)
BiocManager::install(c(
  # RNA-seq
  "DESeq2", "fgsea", "org.Hs.eg.db",
  # Proteomics
  "limma", "vsn", "ComplexHeatmap",
  # Multi-omics
  "MOFA2", "mixOmics", "MultiAssayExperiment"
))
```

See individual pipeline READMEs for complete dependency lists.

## Pipeline Details

### RNA-seq Pipeline

Analyze bulk RNA-seq count data from raw counts to pathway enrichment:

- Gene filtering and annotation (Ensembl/Symbol mapping)
- DESeq2 normalization and differential expression
- PCA, sample correlation, MA/volcano plots
- Gene set enrichment (GO, KEGG) via fgsea

### Proteomics Pipeline

Process quantitative proteomics from MaxQuant, FragPipe, DIA-NN, or Spectronaut:

- Contaminant/reverse filtering
- Log2 transformation + median/quantile/VSN normalization
- MNAR-aware imputation (MinProb, QRILC, kNN)
- limma differential abundance
- Pathway enrichment

### Metabolomics Pipeline

Analyze untargeted or targeted metabolomics:

- Blank subtraction and QC sample processing
- Signal drift correction (LOESS, ComBat)
- PQN/median normalization
- Missing value filtering and imputation
- Differential analysis and metabolite set enrichment

### Multi-omics Integration Pipeline

Integrate 2-3 omics layers (RNA + Protein + Metabolites):

- Flexible input: raw or pre-processed data
- Sample and feature harmonization
- **MOFA2**: Unsupervised factor discovery
- **DIABLO**: Supervised multi-block discriminant analysis
- **SNF**: Network-based sample clustering
- Cross-omics concordance (RNA-protein correlation)
- Combined pathway enrichment

## License

MIT

## Citation

If you use these pipelines, please cite the underlying methods:

- **targets**: Landau (2021). The targets R package. JOSS.
- **DESeq2**: Love et al. (2014). Genome Biology.
- **limma**: Ritchie et al. (2015). Nucleic Acids Research.
- **MOFA2**: Argelaguet et al. (2020). Molecular Systems Biology.
- **mixOmics/DIABLO**: Rohart et al. (2017). PLoS Computational Biology.
