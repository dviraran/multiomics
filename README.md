# Multi-Omics Analysis Pipelines

A collection of reproducible R pipelines for analyzing omics data, built with the [`{targets}`](https://docs.ropensci.org/targets/) framework. Includes an interactive **Shiny viewer** for exploring results.

## Pipelines

| Pipeline | Description | Key Methods |
|----------|-------------|-------------|
| [rnaseq_pipeline](src/pipelines/rnaseq_pipeline/) | RNA-seq differential expression analysis | DESeq2, fgsea, GO/KEGG enrichment |
| [proteomics_pipeline](src/pipelines/proteomics_pipeline/) | Mass spectrometry proteomics (LFQ/DIA) | limma, VSN normalization, MNAR imputation |
| [metabolomics_pipeline](src/pipelines/metabolomics_pipeline/) | LC-MS/GC-MS metabolomics | PQN normalization, drift correction, MetaboAnalyst-style QC |
| [multiomics_pipeline](src/pipelines/multiomics_pipeline/) | Multi-omics data integration | MOFA2, DIABLO (mixOmics), SNF |

## Features

- **Reproducible**: Full dependency tracking with `{targets}` - only re-run what changed
- **Configurable**: YAML-based configuration for all parameters
- **Modular**: Functions organized in `R/` directory for easy customization
- **Comprehensive QC**: Built-in quality control checks and visualizations
- **Publication-ready**: Automated HTML reports with R Markdown
- **Interactive Viewer**: Shiny app for exploring results with interactive plots
- **AI Commentary** (RNA-seq/Proteomics): Optional LLM-generated figure interpretations
- **Auto-install packages**: Prompts to install missing packages on first run

## Directory Structure

```
multiomics/
├── src/
│   ├── pipelines/              # Analysis pipelines
│   │   ├── rnaseq_pipeline/
│   │   ├── proteomics_pipeline/
│   │   ├── metabolomics_pipeline/
│   │   └── multiomics_pipeline/
│   ├── shared/                 # Shared R utilities
│   └── scripts/                # Utility scripts
├── viewer/                     # Shiny visualization app
│   ├── shiny_app/              # The Shiny application
│   ├── run_viewer.R            # Launch the viewer
│   └── create_viewer_package.R # Create distributable package
├── docker/                     # Docker configuration
├── examples/                   # Example data and configs
├── vignettes/                  # Tutorials
└── pipeline_runner.R           # User-facing pipeline launcher
```

## Quick Start

### Option 1: Run from pipeline directory

```bash
cd src/pipelines/<pipeline_name>

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

## Viewing Results with Shiny

After running a pipeline, explore your results interactively with the Shiny viewer.

### Quick Launch

```r
source("viewer/run_viewer.R")

# View results from a pipeline output directory
run_viewer("path/to/pipeline/outputs")
```

### Features

The Shiny viewer provides interactive exploration of:

| Tab | Features |
|-----|----------|
| **RNA-seq** | QC metrics, PCA, DE volcano/MA plots, expression heatmaps, GSEA results |
| **Proteomics** | QC, normalization comparison, DE analysis, pathway enrichment |
| **Metabolomics** | QC, drift correction, differential analysis, pathway mapping |
| **Multi-omics** | MOFA factor exploration, DIABLO integration, cross-omics correlations |

### Create a Distributable Package

Share results with collaborators who don't have R expertise:

```r
source("viewer/create_viewer_package.R")

create_viewer_package(
  data_dirs = list(
    rnaseq = "path/to/rnaseq/outputs",
    proteomics = "path/to/proteomics/outputs",
    metabolomics = "path/to/metabolomics/outputs",
    multiomics = "path/to/multiomics/outputs"
  ),
  output_dir = "my_analysis_viewer",
  project_name = "My Multi-Omics Analysis"
)
```

This creates a standalone folder. Collaborators just need to:

```r
# After opening R in the viewer folder:
source("run_viewer.R")
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

# Shiny viewer packages
install.packages(c("shiny", "bslib", "plotly", "DT", "heatmaply"))
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

## Docker Support

Run pipelines in a containerized environment:

```bash
cd docker
docker-compose up -d
```

See [docker/](docker/) for more details.

## License

MIT

## Citation

If you use these pipelines, please cite the underlying methods:

- **targets**: Landau (2021). The targets R package. JOSS.
- **DESeq2**: Love et al. (2014). Genome Biology.
- **limma**: Ritchie et al. (2015). Nucleic Acids Research.
- **MOFA2**: Argelaguet et al. (2020). Molecular Systems Biology.
- **mixOmics/DIABLO**: Rohart et al. (2017). PLoS Computational Biology.
