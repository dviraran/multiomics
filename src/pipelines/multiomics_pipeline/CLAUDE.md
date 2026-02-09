# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R-based multi-omics integration pipeline built on `{targets}` for combining transcriptomics, proteomics, and metabolomics data. Supports raw and preprocessed inputs, multiple integration methods, and automated HTML report generation.

## Common Commands

```r
library(targets)
tar_make()                          # Run full pipeline
tar_make(target_name)               # Run up to a specific target
tar_visnetwork()                    # Visualize pipeline DAG
tar_read(target_name)               # Read a cached target result
tar_load(target_name)               # Load target into environment
tar_invalidate(target_name)         # Mark target stale, then tar_make() to re-run from there
```

Run with a custom config:
```r
Sys.setenv(PIPELINE_CONFIG = "/path/to/config.yml")
tar_make()
```

There is no formal test suite. Debugging scripts (`debug_*.R`, `diagnostic*.R`) exist in the project root for ad-hoc testing of individual modules.

## Architecture

### Pipeline Framework

`_targets.R` is the entry point. It defines all targets (steps) and their dependencies. The `{targets}` framework handles caching, dependency tracking, and selective re-execution. `tar_option_set(error = "continue")` allows partial pipeline completion on failures.

### Data Flow (target execution order)

```
config.yml → config
  → ingested_data (01_ingestion)
  → gene_protein_mapping (03_mapping)
  → preprocessed_data (02_preprocessing)
  → harmonized_data → mae_data (04_harmonize) [MultiAssayExperiment]
  → feature_data (05_feature_selection)
  → foundational_results (05b) + mechanistic_results (05c)
  → mofa_results (06) | diablo_results (07) | snf_results (08) [parallel branches]
  → concordance_results (09) + consensus_results (09b) + stability_results (09c)
  → enrichment_results (10) → multigsea + pathview (13, 15)
  → commentary (11, AI-powered figure captions)
  → report (analysis_report.Rmd → HTML)
```

### Module Layout (R/)

Files are numbered by execution order. Each file contains functions called by targets in `_targets.R`:

- `00_utils.R` — Config loading (`load_config()`), I/O helpers, logging, statistical utilities
- `01_ingestion.R` — Data loading and validation for all omics types
- `02_preprocessing.R` — DESeq2/Limma normalization and differential analysis
- `03_mapping.R` — Gene/protein ID harmonization (Ensembl, UniProt, gene symbols via org.db)
- `04_harmonize.R` — Sample alignment and `MultiAssayExperiment` (MAE) creation
- `05_feature_selection.R` — Strategies: `all`, `variance`, `significant`, `hybrid`
- `05b_foundational_correlations.R` — Cross-omics correlations (runs before integration)
- `05c_mechanistic_inference.R` — Regulatory networks, mediation analysis
- `06_mofa.R` — MOFA2 unsupervised integration
- `07_diablo.R` — DIABLO supervised integration (mixOmics)
- `08_snf.R` — Similarity Network Fusion (optional)
- `09_concordance.R` — RNA-protein concordance analysis
- `09b_integration_consensus.R` — Cross-method comparison (MOFA vs DIABLO vs SNF)
- `09c_stability_analysis.R` — Bootstrap resampling validation
- `10_enrichment.R` — Pathway enrichment (ORA, GSEA, fgsea), Fisher p-value combination
- `11_commentary.R` — AI figure captions (Claude Vision, GPT-4 Vision, or deterministic fallback)
- `13_multigsea_plots.R` — MultiGSEA bubble plots and heatmaps
- `15_kegg_pathview.R` — KEGG pathway diagrams with fold-change overlays

### Configuration

All pipeline behavior is controlled by `config.yml`. Key sections:
- `global` — Which omics are present, organism, metadata path
- `design` — Condition column, contrasts, reference level
- `transcriptomics/proteomics/metabolomics` — Input mode (`raw` vs `preprocessed`), processing parameters
- `integration.methods` — Which methods to run (MOFA2, DIABLO, SNF, spls)
- `enrichment` — Enrichment method, MultiGSEA settings
- `commentary` — AI backend selection (claude/openai/none)

### Key Data Structures

- **MAE** (`MultiAssayExperiment`): Central data object after harmonization. Contains aligned matrices for all omics with shared sample metadata.
- **Preprocessed data**: Named list with per-omics entries, each containing `normalized_matrix`, `de_results`, and processing metadata.
- **Integration results**: Named lists with `results`, `plots`, and `tables` sub-elements.

### Organism Support

Built-in annotation packages for H. sapiens (`org.Hs.eg.db`) and C. elegans (`org.Ce.eg.db`). Non-model organisms are supported via custom mapping files (`gene_protein_mapping.csv`) and GMT files for enrichment.

### Output Structure

```
outputs/
├── tables/          # CSV results (summaries, factors, loadings, enrichment)
├── plots/           # PNG visualizations
├── enrichment/multigsea/pathview/  # KEGG pathway diagrams
├── report/          # Final HTML report
├── stability/       # Bootstrap validation results
├── commentary/      # AI-generated figure descriptions
└── logs/            # Execution logs
```

## Key Dependencies

**Framework**: targets, tarchetypes
**Bioconductor**: DESeq2, limma, MultiAssayExperiment, SummarizedExperiment, MOFA2, mixOmics, clusterProfiler, fgsea, pathview, ComplexHeatmap
**CRAN**: tidyverse, yaml, patchwork, SNFtool (optional)
**Reports**: rmarkdown, knitr

Missing packages are auto-detected and prompted for installation at pipeline runtime.
