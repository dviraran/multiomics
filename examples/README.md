# Multi-Omics Pipeline Examples

This directory contains scripts to run all four multi-omics pipelines on real public data.

## Quick Start

```bash
# From the examples/ directory
cd /path/to/multiomics/examples

# Step 1: Download example data (NCI-60 cancer cell line panel)
Rscript scripts/download_nci60_data.R

# Step 2: Update all pipeline configs to use the example data
Rscript scripts/setup_nci60_configs.R

# Step 3: Run all pipelines
Rscript scripts/run_all_pipelines.R
```

## Dataset: NCI-60 Cancer Cell Line Panel

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

- **Condition**: Solid tumors vs. Hematological (Leukemia)
- **Contrast**: `solid - hematological`
- **Organism**: Homo sapiens

## Scripts

### `download_nci60_data.R`
Downloads and processes data into the standard format expected by each pipeline:
- `rna_counts_matrix.csv` - Gene expression counts
- `prot_intensity_matrix.csv` - Protein intensities
- `metab_feature_matrix.csv` - Metabolite abundances
- `metadata.csv` - Sample annotations

### `setup_nci60_configs.R`
Updates all four pipeline `config.yml` files:
- Sets data file paths
- Configures experimental design
- Sets analysis parameters appropriate for cell line data

### `run_all_pipelines.R`
Master script to run all pipelines:

```bash
# Run all sequentially
Rscript scripts/run_all_pipelines.R

# Run single-omics pipelines in parallel
Rscript scripts/run_all_pipelines.R --parallel

# Skip metabolomics (if not needed)
Rscript scripts/run_all_pipelines.R --skip-metabolomics

# Run only one pipeline
Rscript scripts/run_all_pipelines.R --only=rnaseq
Rscript scripts/run_all_pipelines.R --only=proteomics
Rscript scripts/run_all_pipelines.R --only=multiomics

# Clean and re-run everything
Rscript scripts/run_all_pipelines.R --clean
```

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

## Customizing for Your Data

### Using Your Own Data

1. Place your data files in `examples/data/mydata/`:
   - RNA-seq: counts matrix (genes × samples)
   - Proteomics: intensity matrix (proteins × samples)
   - Metabolomics: feature matrix (features × samples)
   - Metadata: sample annotations with condition column

2. Modify `setup_nci60_configs.R` or create a new setup script

3. Key config parameters to update:
   ```yaml
   # File paths
   counts_file: "path/to/counts.csv"
   metadata_file: "path/to/metadata.csv"

   # Design
   sample_id_column: "your_sample_column"
   condition_column: "your_condition"
   contrasts:
     - "treatment - control"
   ```

### Alternative Public Datasets

| Dataset | Description | Download |
|---------|-------------|----------|
| TCGA-BRCA | Breast cancer (TCGA) | [GDC Portal](https://portal.gdc.cancer.gov/) |
| CPTAC | Proteogenomics | [PDC](https://proteomic.datacommons.cancer.gov/) |
| MetaboLights | Metabolomics studies | [EMBL-EBI](https://www.ebi.ac.uk/metabolights/) |

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce dataset size or increase R memory limit
   ```r
   options(future.globals.maxSize = 2000 * 1024^2)  # 2GB
   ```

2. **Missing packages**: Install required Bioconductor packages
   ```r
   BiocManager::install(c("DESeq2", "limma", "fgsea", "MOFA2", "mixOmics"))
   ```

3. **API errors for commentary**: Set environment variables
   ```bash
   export ANTHROPIC_API_KEY="your-key"  # For Claude
   export OPENAI_API_KEY="your-key"     # For OpenAI
   ```
   Or use `backend: "none"` in config for data-driven commentary.

### Getting Help

- Check pipeline-specific READMEs in each pipeline directory
- Open an issue on GitHub
- Review the `CLAUDE.md` file for architecture details

## Requirements

- R ≥ 4.1.0
- Bioconductor ≥ 3.14
- Python ≥ 3.8 (for AI commentary)
- ~4GB RAM recommended
- ~2GB disk space for outputs
