# Proteomics Analysis Pipeline

A comprehensive, reproducible mass-spectrometry proteomics analysis pipeline built with `{targets}` in R. Supports label-free (LFQ/DIA) workflows and handles outputs from MaxQuant, FragPipe/MSFragger, Spectronaut, DIA-NN, Skyline, or generic protein/peptide intensity tables.

## Features

- **Modular design** using `{targets}` for reproducibility and caching
- **Flexible input** support for various MS quantification tools
- **Comprehensive QC** with PCA, UMAP, correlation analysis, and outlier detection
- **Robust normalization** (VSN, median, quantile, cyclic loess)
- **Missing value handling** with MNAR-aware imputation (QRILC, MinProb, KNN)
- **Differential abundance** analysis using limma with multiple contrast support
- **Pathway analysis** with ORA and GSEA (supports custom GMT for any organism)
- **Automated HTML report** with interactive visualizations

## Installation

### Prerequisites

R version ≥ 4.1.0

### Required Packages

```r
# Install CRAN packages
install.packages(c(
  "targets",
  "tarchetypes",
  "tidyverse",
  "yaml",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "ggrepel",
  "umap",
  "rmarkdown",
  "DT",
  "plotly",
  "dendextend"
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "limma",
  "vsn",
  "fgsea",
  "impute",
  "imputeLCMD",
  "clusterProfiler",
  "org.Hs.eg.db",      # For human
  "org.Mm.eg.db",      # For mouse (optional)
  "msigdbr",           # For MSigDB gene sets
  "UniProt.ws",        # For ID mapping
  "biomaRt"            # Alternative ID mapping
))
```

### Optional: Set up renv for reproducibility

```r
# Initialize renv in the project directory
renv::init()
renv::snapshot()
```

## Quick Start

### 1. Prepare Your Data

Place your input files in the `data/` directory:

```
data/
├── quant_matrix.csv    # Protein/peptide intensities (required)
├── metadata.csv        # Sample metadata (required)
├── id_mapping.csv      # Custom ID mapping (optional)
└── custom_pathways.gmt # Custom gene sets (optional)
```

**quant_matrix.csv format:**
- First column: feature IDs (protein accessions, UniProt IDs, or gene symbols)
- Remaining columns: sample intensities (one column per sample)
- Values: raw intensities (not counts), can contain zeros/NAs

```csv
protein_id,sample1,sample2,sample3,sample4
sp|P04406|G3P_HUMAN,1234567,2345678,1876543,2109876
sp|P68371|TBB4B_HUMAN,987654,876543,765432,654321
```

**metadata.csv format:**
- Must include sample identifier column matching quant_matrix column names
- Must include condition/group column

```csv
sample_id,condition,batch,replicate
sample1,control,batch1,1
sample2,control,batch1,2
sample3,treated,batch2,1
sample4,treated,batch2,2
```

### 2. Configure the Pipeline

Edit `config.yml` to match your experiment:

```yaml
input:
  quant_matrix: "data/quant_matrix.csv"
  metadata: "data/metadata.csv"
  feature_id_type: "uniprot"  # or "gene_symbol", "protein_accession"
  organism: "Homo sapiens"
  sample_id_column: "sample_id"

design:
  condition_column: "condition"
  design_formula: "~ condition"
  contrasts:
    - "treated - control"
  reference_level: "control"
```

### 3. Run the Pipeline

```r
# Load targets
library(targets)

# Visualize the pipeline (optional)
tar_visnetwork()

# Run the complete pipeline
tar_make()

# Check pipeline status
tar_progress()
```

### 4. View Results

- **HTML Report**: `outputs/report/analysis_report.html`
- **Differential results**: `outputs/tables/differential_*.csv`
- **Pathway results**: `outputs/tables/pathway_*.csv`
- **QC plots**: `outputs/qc/`
- **Analysis plots**: `outputs/plots/`

## Configuration Reference

### Input Settings

| Parameter | Description | Default |
|-----------|-------------|---------|
| `quant_matrix` | Path to intensity matrix | Required |
| `metadata` | Path to sample metadata | Required |
| `feature_id_type` | ID type: `uniprot`, `gene_symbol`, `protein_accession`, `peptide_sequence`, `custom` | `uniprot` |
| `feature_level` | `protein` or `peptide` | `protein` |
| `organism` | Species name (e.g., "Homo sapiens") | Required |
| `sample_id_column` | Column in metadata with sample IDs | `sample_id` |

### Optional Inputs

| Parameter | Description |
|-----------|-------------|
| `mapping_file` | Custom ID → gene symbol mapping (CSV) |
| `contaminants_file` | List of contaminant IDs to remove |
| `gene_set_gmt` | Custom gene sets for pathway analysis (GMT format) |

### Processing Settings

| Parameter | Options | Default |
|-----------|---------|---------|
| `zeros_as_na` | `true`/`false` | `true` |
| `log_transform` | `true`/`false` | `true` |
| `normalization_method` | `vsn`, `median`, `quantile`, `cyclicloess`, `none` | `vsn` |

### Filtering Settings

| Parameter | Description | Default |
|-----------|-------------|---------|
| `global_min_presence` | Min % samples feature must be present in | `0.5` (50%) |
| `group_min_presence` | Min % samples in at least one group | `0.7` (70%) |
| `low_intensity_quantile` | Remove bottom X% by intensity | `0` (disabled) |

### Imputation Settings

| Parameter | Options | Default |
|-----------|---------|---------|
| `method` | `QRILC`, `MinProb`, `MinDet`, `knn`, `random_forest`, `none` | `QRILC` |
| `min_prob_quantile` | Quantile for MinProb/MinDet | `0.01` |

### Differential Analysis

| Parameter | Description | Default |
|-----------|-------------|---------|
| `method` | `limma` or `limma_trend` | `limma` |
| `adj_pvalue_threshold` | FDR significance threshold | `0.05` |
| `log2fc_threshold` | Minimum absolute log2 fold change | `1` |

### Pathway Analysis

| Parameter | Description | Default |
|-----------|-------------|---------|
| `run_pathway_analysis` | Enable/disable | `true` |
| `databases` | Gene set databases (`GO`, `KEGG`, `Reactome`) | `["GO", "KEGG"]` |
| `go_ontologies` | GO ontologies (`BP`, `MF`, `CC`) | `["BP", "MF"]` |
| `min_set_size` | Minimum genes per pathway | `10` |
| `max_set_size` | Maximum genes per pathway | `500` |

## Design Formula and Contrasts

### Simple Two-Group Comparison

```yaml
design:
  design_formula: "~ condition"
  contrasts:
    - "treated - control"
  reference_level: "control"
```

### Multiple Conditions

```yaml
design:
  design_formula: "~ condition"
  contrasts:
    - "drugA - vehicle"
    - "drugB - vehicle"
    - "drugA - drugB"
  reference_level: "vehicle"
```

### With Batch Correction

```yaml
design:
  design_formula: "~ batch + condition"
  batch_column: "batch"
  contrasts:
    - "treated - control"
```

### Paired Samples

```yaml
design:
  design_formula: "~ patient + condition"
  contrasts:
    - "post - pre"
```

## Handling Non-Model Organisms

For organisms without annotation database support:

### Option 1: Provide Custom ID Mapping

Create `data/id_mapping.csv`:

```csv
feature_id,gene_symbol,entrez_id,protein_name
XP_012345,GENEA,12345,Protein A description
XP_012346,GENEB,12346,Protein B description
```

Configure:
```yaml
optional_inputs:
  mapping_file: "data/id_mapping.csv"
```

### Option 2: Provide Custom Gene Sets (GMT)

Download or create GMT file with gene sets matching your feature IDs:

```
PATHWAY_A<TAB>description<TAB>GENEA<TAB>GENEB<TAB>GENEC
PATHWAY_B<TAB>description<TAB>GENED<TAB>GENEE
```

Configure:
```yaml
optional_inputs:
  gene_set_gmt: "data/custom_pathways.gmt"
```

### Option 3: Skip Pathway Analysis

If no mapping is possible:
```yaml
pathway:
  run_pathway_analysis: false
```

The pipeline will still perform QC, normalization, and differential analysis.

## Input Format Examples

### MaxQuant (proteinGroups.txt)

Extract relevant columns:
```r
# Read MaxQuant output
pg <- read.delim("proteinGroups.txt")

# Select intensity columns (LFQ or iBAQ)
intensity_cols <- grep("^LFQ.intensity", colnames(pg), value = TRUE)

# Create quant matrix
quant_matrix <- pg[, c("Protein.IDs", intensity_cols)]
colnames(quant_matrix)[1] <- "protein_id"
colnames(quant_matrix) <- gsub("LFQ.intensity.", "", colnames(quant_matrix))

write.csv(quant_matrix, "data/quant_matrix.csv", row.names = FALSE)
```

### DIA-NN

```r
# Read DIA-NN output
diann <- read.delim("report.pg_matrix.tsv")

# First column is protein IDs, rest are samples
write.csv(diann, "data/quant_matrix.csv", row.names = FALSE)
```

### Spectronaut

Export the "Protein Group Quantity" report in CSV format with samples as columns.

### Generic Format

Any matrix with:
- Rows = proteins/peptides
- Columns = samples
- Values = intensities (not spectral counts)

## Pipeline Outputs

```
outputs/
├── tables/
│   ├── normalized_matrix.csv       # Normalized expression matrix
│   ├── imputed_matrix.csv          # Imputed expression matrix
│   ├── feature_annotations.csv     # ID mapping results
│   ├── differential_*.csv          # DE results per contrast
│   └── pathway_*.csv               # Pathway enrichment results
├── plots/
│   ├── volcano_*.png               # Volcano plots
│   ├── ma_plot_*.png               # MA plots
│   └── pathway_*.png               # Pathway enrichment plots
├── qc/
│   ├── input_summary_*.csv         # Input data summaries
│   ├── filtering_summary.csv       # Filtering log
│   ├── pca_*.png                   # PCA plots
│   ├── correlation_heatmap.png     # Sample correlations
│   ├── outlier_report.csv          # Flagged outliers
│   └── normalization_*.png         # Normalization diagnostics
└── report/
    └── analysis_report.html        # Comprehensive HTML report
```

## Troubleshooting

### "No matching samples between matrix and metadata"

Ensure sample names in your quant_matrix columns exactly match the sample_id column in metadata.

### "Gene symbol mapping failed"

- Check if your organism is supported (human, mouse, rat, zebrafish, yeast, fly)
- Provide a custom `mapping_file` for other organisms
- Ensure UniProt accessions are in standard format (e.g., P04406, not sp|P04406|G3P_HUMAN)

### "Too few significant genes for pathway analysis"

- Relax the `adj_pvalue_threshold` or `log2fc_threshold`
- Check if gene symbol mapping succeeded
- Provide a custom GMT file with gene sets matching your IDs

### Pipeline runs slowly

- For large datasets (>5000 proteins), use `imputation.method: "MinProb"` instead of `random_forest`
- Reduce `pathway.gsea_nperm` from 10000 to 1000 for faster GSEA

### Memory issues

```r
# Increase R memory limit (on Windows)
memory.limit(size = 16000)

# Or run specific targets
tar_make(names = c("ingested_data", "filtered_data"))
```

## Extending the Pipeline

### Adding New Targets

Edit `_targets.R` to add custom analysis steps:

```r
list(
  # ... existing targets ...

  tar_target(
    my_custom_analysis,
    run_my_analysis(imputed_data, config)
  )
)
```

### Adding New Normalization Methods

Edit `R/04_normalization.R` and add to the `apply_normalization()` function:

```r
apply_normalization <- function(mat, method, config) {
  normalized <- switch(method,
    # ... existing methods ...
    "my_method" = normalize_my_method(mat),
    # ...
  )
}
```

## Citation

If you use this pipeline, please cite:

- **limma**: Ritchie et al. (2015) Nucleic Acids Research
- **VSN**: Huber et al. (2002) Bioinformatics
- **fgsea**: Korotkevich et al. (2021) bioRxiv
- **targets**: Landau (2021) JOSS

## License

MIT License

## Contact

For issues and questions, please open a GitHub issue.
