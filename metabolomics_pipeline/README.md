# Metabolomics Analysis Pipeline

A comprehensive, reproducible metabolomics analysis pipeline built with `{targets}` in R. Supports both untargeted LC-MS/GC-MS and targeted metabolomics workflows, with flexible handling of annotations and pathway analysis.

## Features

- **Dual workflow support**: Untargeted (feature-level) and targeted (metabolite-level) data
- **Flexible QC handling**: Automatic detection and use of QC/blank samples
- **Multiple normalization methods**: PQN, TIC, median, quantile, VSN
- **Batch/drift correction**: ComBat, removeBatchEffect, QC-RLSC
- **MNAR-aware imputation**: Half-min, MinProb, QRILC, kNN
- **Comprehensive QC**: PCA, UMAP, correlation analysis, outlier detection
- **Differential analysis**: limma, t-test, Wilcoxon
- **Enrichment analysis**: Pathway ORA, chemical class enrichment, custom GMT sets
- **Automated HTML report** with interactive visualizations

## Installation

### Prerequisites

R version ≥ 4.1.0

### Required Packages

```r
# CRAN packages
install.packages(c(
  "targets",
  "tarchetypes",
  "tidyverse",
  "yaml",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "ggrepel",
  "reshape2",
  "umap",
  "rmarkdown",
  "DT",
  "plotly",
  "dendextend"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "limma",
  "vsn",
  "sva",
  "impute",
  "imputeLCMD"
))

# Optional for PLS-DA
install.packages("mixOmics")
```

## Quick Start

### 1. Prepare Your Data

Place input files in the `data/` directory:

```
data/
├── feature_matrix.csv      # Feature intensities (required)
├── metadata.csv            # Sample metadata (required)
├── feature_metadata.csv    # m/z, RT info (optional, untargeted)
├── annotations.csv         # Metabolite annotations (optional)
├── pathway_mapping.csv     # Pathway annotations (optional)
└── metabolite_sets.gmt     # Custom sets (optional)
```

**feature_matrix.csv format:**
```csv
feature_id,sample1,sample2,sample3,QC1,QC2,Blank1
M100T60,12345,23456,34567,28000,27500,100
M150T90,45678,56789,67890,52000,51000,50
M200T120,78901,89012,90123,82000,81500,200
```

**metadata.csv format:**
```csv
sample_id,sample_type,condition,batch,injection_order
sample1,Sample,control,batch1,1
sample2,Sample,control,batch1,2
sample3,Sample,treated,batch2,3
QC1,QC,QC,batch1,4
QC2,QC,QC,batch2,5
Blank1,Blank,Blank,batch1,6
```

### 2. Configure the Pipeline

Edit `config.yml`:

```yaml
input:
  feature_matrix: "data/feature_matrix.csv"
  metadata: "data/metadata.csv"
  data_type: "untargeted"
  organism: "Homo sapiens"
  sample_id_column: "sample_id"

sample_types:
  sample_type_column: "sample_type"
  sample_value: "Sample"
  qc_value: "QC"
  blank_value: "Blank"

design:
  condition_column: "condition"
  design_formula: "~ condition"
  contrasts:
    - "treated - control"
  reference_level: "control"
```

### 3. Run the Pipeline

```r
library(targets)

# Visualize pipeline
tar_visnetwork()

# Run pipeline
tar_make()

# Check progress
tar_progress()
```

### 4. View Results

- **Report**: `outputs/report/analysis_report.html`
- **Tables**: `outputs/tables/`
- **Plots**: `outputs/plots/`
- **QC outputs**: `outputs/qc/`

## Configuration Reference

### Input Settings

| Parameter | Description | Required |
|-----------|-------------|----------|
| `feature_matrix` | Path to intensity matrix | Yes |
| `metadata` | Path to sample metadata | Yes |
| `data_type` | `untargeted` or `targeted` | Yes |
| `organism` | Species name | Yes |
| `sample_id_column` | Column with sample IDs | Yes |

### Sample Type Settings

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sample_type_column` | Column indicating sample type | `sample_type` |
| `sample_value` | Value for biological samples | `Sample` |
| `qc_value` | Value for QC samples | `QC` |
| `blank_value` | Value for blank samples | `Blank` |

If no sample_type column exists, all samples are treated as biological samples.

### Processing Settings

| Parameter | Options | Default |
|-----------|---------|---------|
| `zeros_as_na` | `true`/`false` | `true` |
| `transform` | `log2`, `log10`, `none` | `log2` |
| `pseudocount` | Numeric | `1` |
| `normalization_method` | `PQN`, `median`, `TIC`, `quantile`, `vsn`, `none` | `PQN` |
| `scaling_method` | `pareto`, `autoscale`, `none` | `pareto` |

### Filtering Settings

| Parameter | Description | Default |
|-----------|-------------|---------|
| `global_min_presence` | Min % samples overall | `0.2` (20%) |
| `group_min_presence` | Min % in any condition | `0.5` (50%) |
| `blank_ratio_threshold` | Max blank/sample ratio | `0.3` |
| `qc_rsd_threshold` | Max RSD in QC samples | `0.3` (30%) |
| `remove_low_variance` | Remove low variance features | `true` |

### Batch Correction

| Parameter | Options | Default |
|-----------|---------|---------|
| `method` | `none`, `combat`, `removeBatchEffect`, `qcrlsc` | `none` |
| `combat_parametric` | Use parametric ComBat | `true` |
| `loess_span` | LOESS span for QC-RLSC | `0.75` |

**Note**: `qcrlsc` requires QC samples and `injection_order_column` in config.

### Imputation Settings

| Parameter | Options | Default |
|-----------|---------|---------|
| `method` | `half_min`, `minprob`, `QRILC`, `knn`, `none` | `half_min` |
| `half_min_fraction` | Fraction of min for half_min | `0.5` |
| `knn_k` | Neighbors for kNN | `10` |

### Differential Analysis

| Parameter | Options | Default |
|-----------|---------|---------|
| `method` | `limma`, `t_test`, `wilcoxon` | `limma` |
| `adj_pvalue_threshold` | FDR cutoff | `0.05` |
| `log2fc_threshold` | Min absolute log2FC | `1` |

## Optional Annotation Files

### feature_metadata.csv (Untargeted)

For untargeted data, provides m/z, RT, and preliminary IDs:

```csv
feature_id,mz,rt,adduct,preliminary_id
M100T60,100.0505,60,M+H,Unknown
M150T90,150.0878,90,M+H,Phenylalanine
M200T120,200.1182,120,M+Na,Unknown
```

### annotations.csv

Maps features to metabolite identifiers and classes:

```csv
feature_id,metabolite_name,hmdb_id,kegg_id,formula,class
M150T90,L-Phenylalanine,HMDB0000159,C00079,C9H11NO2,Amino acids
M180T95,L-Tyrosine,HMDB0000158,C00082,C9H11NO3,Amino acids
```

### pathway_mapping.csv

Maps metabolites to pathways:

```csv
metabolite_id,pathway_id,pathway_name,source
HMDB0000159,hsa00400,Phenylalanine metabolism,KEGG
HMDB0000158,hsa00350,Tyrosine metabolism,KEGG
```

### metabolite_sets.gmt (GMT format)

Custom metabolite sets for enrichment:

```
Amino_acid_metabolism<TAB>Description<TAB>M150T90<TAB>M180T95<TAB>M200T100
Lipid_biosynthesis<TAB>Description<TAB>M300T150<TAB>M350T180
```

## Workflow Without QC/Blanks

The pipeline works gracefully without QC or blank samples:

```yaml
sample_types:
  assume_all_samples: true
```

In this mode:
- All samples are treated as biological samples
- Blank filtering is skipped
- QC-based RSD filtering is skipped
- QC-RLSC drift correction is unavailable

## Workflow Without Annotations

The pipeline runs differential analysis even without annotations:
- Feature IDs are used directly
- Enrichment analysis is skipped with a clear message
- DE results contain feature_id, log2FC, p-values

To enable enrichment for non-annotated data, provide a custom GMT file matching feature IDs.

## Pipeline Outputs

```
outputs/
├── tables/
│   ├── normalized_matrix.csv
│   ├── imputed_matrix.csv
│   ├── batch_corrected_matrix.csv (if enabled)
│   ├── differential_*.csv
│   ├── enrichment_*.csv (if available)
│   └── pca_scores.csv
├── plots/
│   ├── pca_condition.png
│   ├── umap.png
│   ├── volcano_*.png
│   └── enrichment_*.png
├── qc/
│   ├── qc_sample_metrics.csv
│   ├── qc_feature_metrics.csv
│   ├── filtering_summary.csv
│   ├── missingness_summary.csv
│   ├── correlation_heatmap.png
│   └── outlier_report.csv
└── report/
    └── analysis_report.html
```

## Troubleshooting

### "No matching samples"

Ensure sample names in feature_matrix columns exactly match metadata sample_id column.

### "Too few significant features for enrichment"

- Relax `adj_pvalue_threshold` or `log2fc_threshold`
- Check if metabolite IDs match pathway_mapping
- Provide a custom GMT file

### High missingness after filtering

- Lower `global_min_presence` (e.g., 0.1)
- Lower `group_min_presence` (e.g., 0.3)
- Set `blank_ratio_threshold: null` to disable blank filtering

### QC-RLSC not working

Requires:
1. QC samples with `sample_type = "QC"`
2. `injection_order_column` specified in config
3. At least 4 QC samples

### ComBat fails

- Ensure `batch_column` exists in metadata
- Check that each batch has multiple samples
- Try `removeBatchEffect` as alternative

## Extending the Pipeline

### Adding Custom Normalization

Edit `R/04_normalization.R`:

```r
apply_normalization <- function(mat, method, sample_roles, config) {
  normalized <- switch(method,
    # ... existing methods ...
    "my_method" = my_normalization(mat),
    # ...
  )
}
```

### Adding Custom Imputation

Edit `R/06_imputation.R`:

```r
apply_imputation <- function(mat, method, config) {
  imputed <- switch(method,
    # ... existing methods ...
    "my_method" = my_imputation(mat, config),
    # ...
  )
}
```

## Citation

If you use this pipeline, please cite:

- **limma**: Ritchie et al. (2015) Nucleic Acids Research
- **targets**: Landau (2021) JOSS
- **PQN**: Dieterle et al. (2006) Anal Chem

## License

MIT License

## Contact

For issues and questions, please open a GitHub issue.
