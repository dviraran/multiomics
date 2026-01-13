# Multi-Omics Pipeline Examples

This directory contains scripts to run all four multi-omics pipelines on real public data.

## Quick Start (RStudio)

**Important**: Before running any scripts, set your working directory to the `examples/` folder:

```r
# In RStudio: Session > Set Working Directory > Choose Directory
# Or run:
setwd("/path/to/multiomics/examples")
```

### Example 1: TCGA Breast Cancer (RECOMMENDED)

This dataset uses **real gene symbols** that work with pathway analysis (GO, KEGG, Reactome).

```r
# Step 1: Download TCGA breast cancer data from mixOmics
source("scripts/download_breast_tcga_data.R")

# Step 2: Create pipeline configs
source("scripts/setup_breast_tcga_configs.R")

# Step 3: Run all pipelines
source("scripts/run_all_pipelines.R")
```

### Example 2: NCI-60 Human Cancer Cell Lines

```r
# Step 1: Download example data
source("scripts/download_nci60_data.R")

# Step 2: Update all pipeline configs
source("scripts/setup_nci60_configs.R")

# Step 3: Run all pipelines
source("scripts/run_all_pipelines.R")
```

### Example 3: STATegra Mouse B-cell Differentiation

```r
# Step 1: Download mouse data
source("scripts/download_stategra_data.R")

# Step 2: Update all pipeline configs for mouse
source("scripts/setup_stategra_configs.R")

# Step 3: Run all pipelines
source("scripts/run_all_pipelines.R")
```

### Running Individual Pipelines

After setup, you can run pipelines from the root directory using the pipeline runner:

```r
# From multiomics root directory
source("pipeline_runner.R")

# Run with custom config
run_rnaseq_pipeline("examples/data/nci60/rnaseq_config.yml")
run_proteomics_pipeline("examples/data/nci60/proteomics_config.yml")
run_multiomics_pipeline("examples/data/nci60/multiomics_config.yml")
```

Or from individual pipeline directories:

```r
# RNA-seq only
setwd("../rnaseq_pipeline")
targets::tar_make()

# Proteomics only
setwd("../proteomics_pipeline")
targets::tar_make()

# Multi-omics integration
setwd("../multiomics_pipeline")
targets::tar_make()
```

### Reset to Default Configs

```r
# From examples/ directory
source("scripts/reset_configs.R")
```

---

## Dataset 1: TCGA Breast Cancer (RECOMMENDED)

The breast.TCGA dataset from the mixOmics Bioconductor package contains matched multi-omics data from 150 breast cancer patients. **This is the recommended dataset for testing** because it uses real gene symbols that work with pathway analysis.

### Why This Dataset?

- **Real gene symbols** (BRCA1, TP53, ESR1, etc.) - works with GO, KEGG, Reactome
- **Perfectly matched samples** across all omics types
- **Biologically meaningful** comparison (Basal vs Luminal subtypes)
- **Easy download** - single R command from Bioconductor

### Data Contents

| Omics | Features | IDs |
|-------|----------|-----|
| mRNA | 520 genes | Gene symbols |
| miRNA | 184 miRNAs | miRBase IDs |
| Proteomics (RPPA) | 142 proteins | Gene symbols |

### Analysis Design

- **Organism**: Homo sapiens
- **Samples**: 150 breast cancer patients
- **Subtypes**: Basal (45), Her2 (30), LumA (75)
- **Contrast**: `Basal - LumA` (aggressive vs good prognosis)

### Reference

TCGA Network (2012) "Comprehensive molecular portraits of human breast tumours" *Nature* 490:61-70

---

## Dataset 2: NCI-60 Cancer Cell Line Panel (Human)

The NCI-60 is a panel of 60 human cancer cell lines derived from 9 tissue types:
- Breast (5 lines)
- CNS (6 lines)
- Colon (7 lines)
- Leukemia (6 lines)
- Melanoma (9 lines)
- Lung (9 lines)
- Ovarian (7 lines)
- Prostate (2 lines)
- Renal (8 lines)

### Data Sources

| Omics | Source | Reference |
|-------|--------|-----------|
| RNA-seq | [CellMiner](https://discover.nci.nih.gov/cellminer) | Reinhold et al., 2019 |
| Proteomics | [Cell iScience](https://www.sciencedirect.com/science/article/pii/S2589004219304407) | Guo et al., 2019 |
| Metabolomics | NCI DTP | Simulated for demo |

### Analysis Design

- **Organism**: Homo sapiens
- **Condition**: Solid tumors vs. Hematological (Leukemia)
- **Contrast**: `solid - hematological`
- **Samples**: 60 cell lines

---

## Dataset 3: STATegra B-cell Differentiation (Mouse)

The STATegra dataset is a comprehensive multi-omics study of mouse pre-B-cell differentiation using the B3 cell line model.

### Data Sources

| Omics | Source | Accession |
|-------|--------|-----------|
| RNA-seq | ArrayExpress | E-MTAB-2679 |
| Proteomics | PRIDE | PXD001293 |
| Metabolomics | MetaboLights | MTBLS90 |

### Reference

Gomez-Cabrero et al. (2019) "STATegra, a comprehensive multi-omics dataset of B-cell differentiation in mouse" *Scientific Data* 6:256

### Analysis Design

- **Organism**: Mus musculus
- **System**: B3 pre-B-cell line differentiation
- **Condition**: Proliferating (0h) vs. Differentiating (24h)
- **Contrast**: `differentiating - proliferating`
- **Samples**: 6 (3 replicates × 2 conditions)

---

## Scripts

### Download Scripts

| Script | Description |
|--------|-------------|
| `download_breast_tcga_data.R` | Downloads TCGA breast cancer data (RECOMMENDED) |
| `download_nci60_data.R` | Downloads human NCI-60 data |
| `download_stategra_data.R` | Downloads mouse STATegra data |

### Setup Scripts

| Script | Description |
|--------|-------------|
| `setup_breast_tcga_configs.R` | Configures pipelines for TCGA breast cancer |
| `setup_nci60_configs.R` | Configures pipelines for NCI-60 |
| `setup_stategra_configs.R` | Configures pipelines for STATegra |
| `reset_configs.R` | Resets all configs to defaults |

### Runner Script

| Script | Description |
|--------|-------------|
| `run_all_pipelines.R` | Master script to run all pipelines |

#### Run Options (command line)

```bash
# Run all sequentially
Rscript scripts/run_all_pipelines.R

# Run single-omics pipelines in parallel
Rscript scripts/run_all_pipelines.R --parallel

# Skip metabolomics
Rscript scripts/run_all_pipelines.R --skip-metabolomics

# Run only one pipeline
Rscript scripts/run_all_pipelines.R --only=rnaseq
Rscript scripts/run_all_pipelines.R --only=proteomics
Rscript scripts/run_all_pipelines.R --only=multiomics

# Clean and re-run
Rscript scripts/run_all_pipelines.R --clean
```

---

## Output Structure

After running all pipelines, you'll find:

```
rnaseq_pipeline/
  outputs/
    plots/              # QC plots, volcano, MA plots
    tables/             # DE results, enrichment
    report/             # HTML report
    commentary/         # Figure commentary JSON

proteomics_pipeline/
  outputs/
    qc/                 # QC metrics
    plots/              # Analysis figures
    tables/             # DA results
    report/             # HTML report
    commentary/         # Figure commentary

metabolomics_pipeline/
  outputs/
    ...                 # Similar structure

multiomics_pipeline/
  outputs/
    plots/              # MOFA, DIABLO plots
    tables/             # Integration results
    report/             # HTML report
    commentary/         # Figure commentary
```

---

## Customizing for Your Data

### Using Your Own Data

1. Place your data files in `examples/data/mydata/`:
   - RNA-seq: counts matrix (genes × samples)
   - Proteomics: intensity matrix (proteins × samples)
   - Metabolomics: feature matrix (features × samples)
   - Metadata: sample annotations with condition column

2. Copy and modify a setup script (e.g., `setup_nci60_configs.R`)

3. Key config parameters to update:
   ```yaml
   # File paths
   counts_file: "path/to/counts.csv"
   metadata_file: "path/to/metadata.csv"

   # Design
   sample_id_column: "your_sample_column"
   condition_column: "your_condition"
   organism: "Homo sapiens"  # or "Mus musculus"
   contrasts:
     - "treatment - control"
   ```

### Alternative Public Datasets

| Dataset | Description | Download |
|---------|-------------|----------|
| TCGA-BRCA | Breast cancer (TCGA) | [GDC Portal](https://portal.gdc.cancer.gov/) |
| CPTAC | Proteogenomics | [PDC](https://proteomic.datacommons.cancer.gov/) |
| MetaboLights | Metabolomics studies | [EMBL-EBI](https://www.ebi.ac.uk/metabolights/) |

---

## Troubleshooting

### Common Issues

1. **Working directory not set**: Make sure you're in the `examples/` directory
   ```r
   getwd()  # Should end with /examples
   ```

2. **Memory errors**: Reduce dataset size or increase R memory limit
   ```r
   options(future.globals.maxSize = 2000 * 1024^2)  # 2GB
   ```

3. **Missing packages**: The pipelines will automatically prompt you to install missing packages when you run `tar_make()`. Just follow the prompts.

   For manual installation:
   ```r
   # Core packages for all pipelines
   BiocManager::install(c("DESeq2", "limma", "fgsea", "vsn"))

   # Additional packages for multi-omics integration
   BiocManager::install(c(
     "MultiAssayExperiment",  # Multi-omics data structure
     "MOFA2",                 # Multi-Omics Factor Analysis
     "mixOmics",              # DIABLO integration
     "SNFtool",               # Similarity Network Fusion
     "ComplexHeatmap",        # Visualization
     "circlize"               # Circular plots
   ))
   ```

4. **API errors for commentary**: Set environment variables in R
   ```r
   Sys.setenv(ANTHROPIC_API_KEY = "your-key")  # For Claude
   Sys.setenv(OPENAI_API_KEY = "your-key")     # For OpenAI
   ```
   Or use `backend: "none"` in config for data-driven commentary.

### Getting Help

- Check pipeline-specific READMEs in each pipeline directory
- Open an issue on GitHub
- Review the `CLAUDE.md` file for architecture details

---

## Requirements

- R ≥ 4.1.0
- Bioconductor ≥ 3.14
- Python ≥ 3.8 (for AI commentary)
- ~4GB RAM recommended
- ~2GB disk space for outputs
