# Bulk RNA-seq Differential Expression & Pathway Analysis Pipeline

A complete, reproducible pipeline for analyzing bulk RNA-seq data using the `{targets}` framework. Supports model organisms (human, mouse, etc.) and non-model organisms with custom annotation files.

## Features

- **End-to-end analysis**: From raw counts to pathway enrichment
- **Modular design**: Each analysis stage is a separate function
- **Reproducible**: Uses `{targets}` for dependency tracking and caching
- **Flexible**: Supports custom designs, multiple contrasts, and non-model organisms
- **Comprehensive outputs**: CSV tables, plots, and an HTML report

## Pipeline Stages

1. **Data Ingestion**: Load and validate counts matrix and metadata
2. **Gene ID Processing**: Strip Ensembl version suffixes, handle duplicates
3. **Annotation**: Map gene IDs to symbols via biomaRt/OrgDb or custom file
4. **Filtering**: Remove low-expression genes using edgeR's filterByExpr
5. **Normalization**: DESeq2 normalization with VST transformation
6. **QC**: Library sizes, PCA, sample correlation, outlier detection
7. **Differential Expression**: DESeq2 with LFC shrinkage
8. **Pathway Analysis**: GSEA (fgsea) and/or ORA with GO/KEGG/custom gene sets
9. **Figure Commentary**: AI-generated or data-driven figure interpretations
10. **Reporting**: HTML report with all results, visualizations, and commentary

## Quick Start

### 1. Install Dependencies

```r
# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "edgeR", "limma", "AnnotationHub", "biomaRt",
  "clusterProfiler", "fgsea", "enrichplot"
))

# Install CRAN packages
install.packages(c(
  "targets", "tarchetypes", "tidyverse", "ggplot2", "pheatmap",
  "ggrepel", "RColorBrewer", "here", "yaml", "janitor", "matrixStats",
  "qs", "DT", "jsonlite", "htmltools"
))

# Install organism annotation packages (as needed)
BiocManager::install("org.Hs.eg.db")  # Human
BiocManager::install("org.Mm.eg.db")  # Mouse
```

### 2. Prepare Input Data

Place your data in the `data/` folder:

**counts_matrix.csv** (genes x samples):
```csv
gene_id,sample1,sample2,sample3,sample4
ENSG00000000003,1234,2345,1456,2567
ENSG00000000005,567,678,589,690
...
```

**metadata.csv**:
```csv
sample_id,condition,batch
sample1,control,batch1
sample2,control,batch1
sample3,treatment,batch2
sample4,treatment,batch2
```

### 3. Configure the Pipeline

Edit `config.yml` to match your data:

```yaml
# Required
counts_file: "data/counts_matrix.csv"
metadata_file: "data/metadata.csv"
organism: "human"

# Sample settings
sample_id_col: "sample_id"
group_col: "condition"

# Design
design_formula: "~ condition"

# Optional: specify contrasts explicitly
contrasts:
  treatment_vs_control:
    factor: "condition"
    numerator: "treatment"
    denominator: "control"
```

### 4. Run the Pipeline

```r
# From R
library(targets)
tar_make()

# Or from terminal
Rscript -e "targets::tar_make()"
```

### 5. View Results

- **HTML Report**: `outputs/analysis_report.html`
- **DE Results**: `outputs/de_results/`
- **Pathway Results**: `outputs/pathway_results/`
- **Plots**: `outputs/plots/`
- **Normalized Counts**: `outputs/normalized_counts.csv`

## Project Structure

```
rnaseq_pipeline/
├── _targets.R              # Pipeline definition
├── config.yml              # Configuration file
├── R/                      # Analysis functions
│   ├── 01_data_ingestion.R
│   ├── 02_annotation.R
│   ├── 03_filtering.R
│   ├── 04_normalization.R
│   ├── 05_qc.R
│   ├── 06_differential_expression.R
│   ├── 07_pathway_analysis.R
│   ├── 08_reporting.R
│   └── 09_commentary.R     # Figure commentary generation
├── scripts/
│   ├── figure_commentary_claude.py  # Claude Vision API script
│   └── figure_commentary_openai.py  # OpenAI Vision API script
├── reports/
│   └── analysis_report.Rmd
├── data/                   # Input data (user-provided)
│   ├── counts_matrix.csv
│   ├── metadata.csv
│   └── (optional) gene_mapping.csv
├── outputs/                # Generated outputs
│   ├── analysis_report.html
│   ├── normalized_counts.csv
│   ├── vst_counts.csv
│   ├── qc_metrics.csv
│   ├── de_results/
│   ├── pathway_results/
│   ├── plots/
│   └── commentary/         # Figure commentary JSON/CSV files
└── _targets/               # Targets cache (auto-generated)
```

## Configuration Reference

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `counts_file` | Path to counts matrix CSV/TSV |
| `metadata_file` | Path to sample metadata CSV/TSV |
| `organism` | Organism name (e.g., "human", "mouse") |

### Design Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sample_id_col` | "sample_id" | Metadata column with sample IDs |
| `group_col` | "condition" | Primary grouping variable |
| `design_formula` | "~ condition" | DESeq2 design formula |

### Contrast Specification

```yaml
contrasts:
  comparison_name:
    factor: "condition"      # Factor in design
    numerator: "treatment"   # Test group
    denominator: "control"   # Reference group
```

### Filtering & Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_count` | 10 | Minimum read count |
| `alpha` | 0.05 | Significance threshold |
| `lfc_threshold` | 0 | Log2FC threshold |
| `pathway_method` | "fgsea" | "fgsea", "ora", or "both" |

## Figure Commentary

The pipeline includes an automated figure commentary system that generates interpretive text for each QC and DE visualization. This helps users understand what each plot shows and identify potential issues.

### Commentary Backends

**Backend A: Claude Vision API**

Uses Claude's vision capabilities to analyze figures and generate context-aware commentary.

```yaml
# config.yml
commentary_enabled: true
commentary_backend: "claude"
commentary_model: "claude-sonnet-4-20250514"
commentary_max_tokens: 1500
```

**Requirements:**
- Set the `ANTHROPIC_API_KEY` environment variable
- Python 3.8+ with the `anthropic` package installed

```bash
# Set API key
export ANTHROPIC_API_KEY="your-api-key-here"

# Install Python dependency
pip install anthropic
```

**Backend B: OpenAI GPT-4 Vision API**

Uses OpenAI's GPT-4 Vision to analyze figures.

```yaml
# config.yml
commentary_enabled: true
commentary_backend: "openai"
commentary_openai_model: "gpt-4o"
commentary_max_tokens: 1500
```

**Requirements:**
- Set the `OPENAI_API_KEY` environment variable
- Python 3.8+ with the `openai` package installed

```bash
# Set API key
export OPENAI_API_KEY="your-api-key-here"

# Install Python dependency
pip install openai
```

**Backend C: Data-Driven Fallback (No LLM, Offline)**

Uses deterministic heuristics based on QC metrics and DE results to generate commentary. No API calls required.

```yaml
# config.yml
commentary_enabled: true
commentary_backend: "none"
```

This backend:
- Analyzes library size distributions and flags low-depth samples
- Evaluates PCA separation between groups
- Identifies correlation-based outliers
- Summarizes DE results (up/down counts, significance rates)

### Enabling/Disabling Commentary

```yaml
# Disable commentary entirely
commentary_enabled: false

# Enable with data-driven fallback (default)
commentary_enabled: true
commentary_backend: "none"

# Enable with Claude Vision
commentary_enabled: true
commentary_backend: "claude"

# Enable with OpenAI GPT-4 Vision
commentary_enabled: true
commentary_backend: "openai"
```

### Commentary Outputs

When enabled, commentary is saved to `outputs/commentary/`:
- `*_commentary.json` - Individual JSON file per figure
- `commentary_all.json` - Combined JSON with all commentary
- `commentary_summary.csv` - CSV summary table
- `figures_metadata.csv` - Metadata table of all figures

### Commentary Structure

Each figure receives commentary with:
1. **What this figure is** - Brief description of the visualization
2. **Key observations** - 3-8 specific patterns observed
3. **Potential issues/checks** - 0-5 quality concerns or verification steps
4. **Recommended next steps** - 0-5 actionable suggestions

### Figures Covered

The commentary system covers all standard pipeline figures:
- Library size barplot
- Detected genes per sample
- Sample correlation heatmap
- PCA plot (PC1 vs PC2)
- PCA scree plot
- Volcano plots (per contrast)
- MA plots (per contrast)
- DE gene heatmaps (per contrast)
- Pathway enrichment plots

### Adding New Figures

To add commentary for new figures:

1. **Register the figure** in `build_figures_table()` (R/09_commentary.R):

```r
figures$my_new_plot <- data.frame(
  figure_id = "my_new_plot",
  filepath = my_plots$new_plot,
  plot_type = "scatter",
  section = "My Section",
  title = "My New Plot",
  description = "Description of what this plot shows",
  x_axis = "X-axis label",
  y_axis = "Y-axis label",
  contrast = NA_character_,
  stringsAsFactors = FALSE
)
```

2. **Add fallback heuristics** (optional, for data-driven backend):

```r
fallback_my_new_plot <- function(data_object) {
  # Analyze data and return commentary list
  list(
    what_is_this = "Description...",
    observations = list("Observation 1", "Observation 2"),
    issues_checks = list("Issue to check"),
    next_steps = list("Suggested action")
  )
}
```

3. **Add static explanation** in `get_static_explanation()`:

```r
my_new_plot = "
  <div class='static-explanation'>
  <details>
  <summary><strong>Why include this plot?</strong></summary>
  <p>Educational text about why this plot is useful...</p>
  </details>
  </div>
"
```

### Caching Behavior

Commentary generation is cached by `{targets}`:
- Commentary regenerates if figures change
- Commentary regenerates if config changes (e.g., backend switch)
- Manually invalidate with: `tar_invalidate(commentary_tbl)`

### Troubleshooting Commentary

**"Commentary unavailable" in report:**
- Check that `commentary_enabled: true` in config.yml
- Verify figures were generated (check `outputs/plots/`)
- Check `outputs/commentary/` for error messages in JSON files

**Claude API errors:**
- Verify `ANTHROPIC_API_KEY` is set correctly
- Check API quota and rate limits
- The pipeline will retry up to `commentary_max_retries` times

**OpenAI API errors:**
- Verify `OPENAI_API_KEY` is set correctly
- Check API quota and rate limits
- Ensure you have access to the vision-capable models (gpt-4o, gpt-4o-mini)

**Fallback commentary is generic:**
- Ensure data objects are passed correctly in `_targets.R`
- Check that the figure_id matches expected patterns

## Non-Model Organisms

For organisms without standard annotation packages:

### Option 1: Custom Mapping File

Create `data/gene_mapping.csv`:
```csv
gene_id,symbol,entrez_id,description
GENE001,ABC1,12345,ABC transporter 1
GENE002,DEF2,12346,Defense protein 2
```

Configure:
```yaml
organism: "my_organism"
mapping_file: "data/gene_mapping.csv"
```

### Option 2: Custom Gene Sets (GMT)

Provide pathway gene sets in GMT format:
```
PathwayA    description    GENE001    GENE002    GENE003
PathwayB    description    GENE004    GENE005
```

Configure:
```yaml
organism: "my_organism"
gmt_file: "data/pathways.gmt"
```

## Targets Commands

```r
# Run the full pipeline
tar_make()

# View pipeline graph
tar_visnetwork()

# Check pipeline status
tar_progress()

# Load a specific target
qc_metrics <- tar_read(qc_metrics)
de_results <- tar_read(de_results_annotated)

# Invalidate and re-run specific target
tar_invalidate(de_results_list)
tar_make()

# Clean all targets
tar_destroy()
```

## Troubleshooting

### Common Issues

1. **Sample names don't match**
   - Ensure `sample_id` column in metadata exactly matches counts column names
   - Check for whitespace or special characters

2. **biomaRt connection fails**
   - Set `annotation_source: "annotationhub"` in config
   - Or provide a custom `mapping_file`

3. **Out of memory**
   - Reduce `pathway_max_size`
   - Filter more aggressively with higher `min_count`

4. **No significant DE genes**
   - Check PCA plot for batch effects
   - Consider adding batch to design formula
   - Verify experimental design

## Citation

If you use this pipeline, please cite the underlying tools:
- DESeq2: Love et al., 2014
- fgsea: Korotkevich et al., 2021
- clusterProfiler: Wu et al., 2021
- targets: Landau, 2021

## License

MIT License
