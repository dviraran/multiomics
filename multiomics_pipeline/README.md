# Multi-Omics Integration Pipeline

A comprehensive R pipeline for integrating transcriptomics, proteomics, and metabolomics data using the `{targets}` framework. Supports both raw and preprocessed inputs, with multiple integration methods including MOFA2 and DIABLO.

## Features

- **Flexible Input Modes**: Accept raw data (counts, intensities) or preprocessed normalized matrices
- **Multi-Omics Support**: Transcriptomics, proteomics, and/or metabolomics in any combination
- **Integration Methods**:
  - **MOFA2** (Multi-Omics Factor Analysis): Unsupervised factor discovery
  - **DIABLO** (mixOmics): Supervised multi-block discriminant analysis
  - **SNF** (Similarity Network Fusion): Network-based sample clustering (optional)
- **Cross-Omics Analysis**: RNA-protein concordance, differential concordance
- **Pathway Enrichment**: Per-omics and combined multi-omics enrichment
- **Reproducible**: Built on `{targets}` for full pipeline reproducibility
- **Comprehensive Reporting**: HTML report with visualizations

## Installation

### Automatic Package Installation

The pipeline will automatically check for missing packages and prompt you to install them when you run `tar_make()`. Just run the pipeline and follow the prompts.

### Manual Installation (optional)

```r
# CRAN packages
install.packages(c(
  "targets", "tarchetypes", "tidyverse", "yaml", "tibble",
  "ggplot2", "patchwork", "rmarkdown", "DT"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "limma", "edgeR",
  "SummarizedExperiment", "MultiAssayExperiment", "S4Vectors",
  "MOFA2", "mixOmics",
  "clusterProfiler", "org.Hs.eg.db", "fgsea",
  "ComplexHeatmap"
))

# Optional
install.packages("SNFtool")  # For SNF integration
```

## Quick Start

1. **Clone/copy the pipeline directory**

2. **Prepare your data** in the `data/` folder:
   - Sample metadata CSV
   - Omics data matrices (raw or preprocessed)

3. **Configure** `config.yml`:
   ```yaml
   global:
     omics_present:
       - "transcriptomics"
       - "proteomics"
     metadata: "data/metadata.csv"
   ```

4. **Run the pipeline**:

   **Option A: Run from pipeline directory**
   ```r
   library(targets)
   tar_make()
   ```

   **Option B: Run with custom config from anywhere**
   ```r
   # From the multiomics root directory
   source("pipeline_runner.R")
   run_multiomics_pipeline("path/to/your/config.yml")

   # Or set environment variable directly
   Sys.setenv(PIPELINE_CONFIG = "/path/to/your/config.yml")
   setwd("multiomics_pipeline")
   targets::tar_make()
   ```

5. **View results** in `outputs/`

## Directory Structure

```
multiomics_pipeline/
├── _targets.R              # Main pipeline definition
├── config.yml              # Configuration file
├── R/                      # R modules
│   ├── 00_utils.R          # Utility functions
│   ├── 01_ingestion.R      # Data loading
│   ├── 02_preprocessing.R  # Per-omics preprocessing
│   ├── 03_mapping.R        # Identifier harmonization
│   ├── 04_harmonize.R      # Sample/feature harmonization, MAE creation
│   ├── 05_feature_selection.R  # Feature selection
│   ├── 06_mofa.R           # MOFA2 integration
│   ├── 07_diablo.R         # DIABLO integration
│   ├── 08_snf.R            # SNF integration
│   ├── 09_concordance.R    # Cross-omics concordance
│   └── 10_enrichment.R     # Pathway enrichment
├── data/                   # Input data
├── outputs/                # Results
│   ├── tables/             # CSV outputs
│   ├── plots/              # Visualizations
│   └── report/             # HTML report
└── reports/                # R Markdown templates
```

## Configuration

### Global Settings

```yaml
global:
  omics_present:
    - "transcriptomics"
    - "proteomics"
    - "metabolomics"  # Uncomment to include
  organism: "Homo sapiens"
  metadata: "data/metadata.csv"
  sample_id_column: "sample_id"
```

### Experimental Design

```yaml
design:
  condition_column: "condition"
  design_formula: "~ condition"
  contrasts:
    - "treated - control"
  reference_level: "control"
  batch_column: null  # Optional batch correction
```

### Input Modes

Each omics can be in "raw" or "preprocessed" mode:

#### Raw Mode (performs full processing)
```yaml
transcriptomics:
  mode: "raw"
  raw:
    counts_matrix: "data/rna_counts_matrix.csv"
    mapping_file: null  # Optional Ensembl -> Symbol mapping
    gmt_file: null      # Optional gene sets
```

#### Preprocessed Mode (skip processing)
```yaml
transcriptomics:
  mode: "preprocessed"
  preprocessed:
    normalized_matrix: "data/rna_normalized.csv"
    de_table: "data/rna_de_results.csv"  # Optional
```

### Integration Methods

```yaml
integration:
  methods:
    - "MOFA2"      # Required
    - "DIABLO"     # Required
    # - "SNF"      # Optional

  mofa2:
    num_factors: 10
    convergence_mode: "fast"
    seed: 42

  diablo:
    ncomp: 3
    design_matrix: "full"
    cv_folds: 5
    cv_repeats: 10
```

### Feature Selection

```yaml
feature_selection:
  transcriptomics:
    strategy: "hybrid"      # variance, significant, hybrid, all
    top_n_variance: 2000
    fdr_threshold: 0.05

  proteomics:
    strategy: "hybrid"
    top_n_variance: 1000
    fdr_threshold: 0.05
```

## Input Data Formats

### Metadata CSV

```csv
sample_id,condition,batch
sample1,control,batch1
sample2,control,batch1
sample3,treated,batch1
sample4,treated,batch2
```

### Transcriptomics (Raw Counts)

```csv
gene_id,sample1,sample2,sample3,sample4
ENSG00000000003,1234,2345,3456,4567
ENSG00000000005,234,345,456,567
```

### Proteomics (Intensity Matrix)

```csv
protein_id,sample1,sample2,sample3,sample4
sp|P12345|PROT1,1000000,1200000,800000,900000
sp|Q12345|PROT2,500000,600000,NA,550000
```

### Metabolomics (Feature Matrix)

```csv
feature_id,sample1,sample2,sample3,sample4
M100.0505_60,50000,55000,45000,48000
M150.0878_90,30000,32000,28000,29000
```

## Pipeline Outputs

### Tables (CSV)

| File | Description |
|------|-------------|
| `sample_alignment.csv` | Sample presence across omics |
| `mae_summary.csv` | MAE structure summary |
| `feature_selection_summary.csv` | Features selected per omics |
| `mofa_factors.csv` | MOFA factor values per sample |
| `mofa_variance_explained.csv` | R² per factor per view |
| `mofa_weights_*.csv` | Feature weights per view |
| `diablo_scores_*.csv` | DIABLO sample scores |
| `diablo_loadings_*.csv` | DIABLO feature loadings |
| `rna_protein_correlations.csv` | Gene-level RNA-protein correlations |
| `combined_enrichment.csv` | Multi-omics pathway enrichment |

### Plots (PNG/PDF)

- MOFA variance heatmap, factor scatter plots, weight plots
- DIABLO sample plots, correlation circles, loadings
- SNF similarity heatmap, MDS plots
- RNA-protein concordance distributions
- Enrichment dotplots

### Report

- `outputs/report/analysis_report.html`: Comprehensive HTML report

## Usage Examples

### Run Full Pipeline

```r
library(targets)
tar_make()
```

### Visualize Pipeline

```r
tar_visnetwork()
```

### Access Results

```r
# Load specific targets
tar_load(mae_data)
tar_load(integration_results)
tar_load(concordance_results)

# Read MOFA factors
mofa_factors <- tar_read(mofa_results)$results$factors

# Get pipeline summary
summary <- tar_read(pipeline_summary)
```

### Selective Re-run

```r
# Invalidate and rerun from preprocessing
tar_invalidate(preprocessed_data)
tar_make()
```

## Integration Methods

### MOFA2

Multi-Omics Factor Analysis discovers latent factors that capture sources of variation across all omics layers. Ideal for:
- Exploratory analysis
- Identifying shared vs omics-specific variation
- Sample stratification

### DIABLO

Data Integration Analysis for Biomarker discovery using Latent variable approaches for Omics. Supervised method ideal for:
- Classification tasks
- Biomarker discovery
- Feature selection
- Cross-omics feature relationships

### SNF (Optional)

Similarity Network Fusion clusters samples based on combined omics similarity. Useful for:
- Subtype discovery
- Patient stratification
- Unsupervised clustering validation

## Cross-Omics Concordance

The pipeline performs:

1. **Expression Concordance**: Correlate RNA and protein levels for matching genes
2. **Differential Concordance**: Compare fold changes between DE genes and DA proteins
3. **Directional Agreement**: Percentage of features changing in same direction

## Non-Model Organisms

For organisms without annotation packages:

1. Provide custom mapping files:
   ```yaml
   transcriptomics:
     raw:
       mapping_file: "data/gene_mapping.csv"  # ensembl_id, gene_symbol
   ```

2. Provide custom GMT files for enrichment:
   ```yaml
   transcriptomics:
     raw:
       gmt_file: "data/custom_pathways.gmt"
   ```

## Troubleshooting

### "No common samples across omics"
- Check that sample IDs match exactly across metadata and data matrices
- Verify `sample_id_column` in config

### MOFA2 fails to converge
- Try reducing `num_factors`
- Use `convergence_mode: "slow"` for better convergence
- Check for features with zero variance

### DIABLO CV error
- Ensure sufficient samples per group (>5 recommended)
- Reduce `ncomp` if too many components requested

### Memory issues
- Reduce `top_n_variance` in feature selection
- Run with fewer omics layers initially

## Citation

If you use this pipeline, please cite the underlying methods:

- **MOFA2**: Argelaguet et al. (2020) Molecular Systems Biology
- **DIABLO/mixOmics**: Rohart et al. (2017) PLoS Computational Biology
- **SNF**: Wang et al. (2014) Nature Methods
- **targets**: Landau (2021) JOSS

## License

MIT License
