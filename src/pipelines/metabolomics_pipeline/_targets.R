# =============================================================================
# Metabolomics Analysis Pipeline - Main targets file
# =============================================================================
# =============================================================================
# Metabolomics Analysis Pipeline - Main targets file
# =============================================================================
# A comprehensive, modular pipeline for LC-MS/GC-MS metabolomics analysis
# Supports both untargeted and targeted workflows
#
# Run with: targets::tar_make()
# Visualize: targets::tar_visnetwork()
# =============================================================================

# =============================================================================
# Package check - prompt to install missing packages
# =============================================================================
check_and_install_packages <- function(packages) {
  # Separate CRAN and Bioconductor packages
  bioc_packages <- c("limma")

  missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    message("\n========================================")
    message("Missing packages detected:")
    message(paste(" -", missing, collapse = "\n"))
    message("========================================\n")

    answer <- readline(prompt = "Would you like to install them? (y/n): ")

    if (tolower(answer) == "y") {
      # Install BiocManager if needed
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
set.seed(1234)
tar_option_set(
  packages = c("MetaboAnalystR", "ggplot2", "plotly", "dplyr", "readr", "tibble", "stringr", "purrr", "tidyr", "yaml", "logger"),
  resources = list(cores = 4),
  format = "rds"
)
      stop("Cannot proceed without required packages: ", paste(missing, collapse = ", "))
    }
  }
}

# Define required packages
required_packages <- c(
  "tidyverse", "readr", "yaml", "tibble",
  "ggplot2", "pheatmap", "RColorBrewer", "ggrepel", "reshape2",
  "limma"
)

# Check packages before loading targets
check_and_install_packages(required_packages)

library(targets)
library(tarchetypes)

# Set options
options(tidyverse.quiet = TRUE)
tar_option_set(
  packages = c(
    # Core
    "tidyverse",
    "readr",
    "yaml",
    "tibble",

    # Visualization
    "ggplot2",
    "pheatmap",
    "RColorBrewer",
    "ggrepel",
    "reshape2",

    # Statistics
    "limma"
  ),
  format = "rds",
  error = "continue"
)

# Source all R functions
source("R/00_utils.R")
source("R/01_data_ingestion.R")
source("R/02_initial_qc.R")
source("R/03_filtering.R")
source("R/04_normalization.R")
source("R/05_batch_correction.R")
source("R/06_imputation.R")
source("R/07_exploratory.R")
source("R/08_differential.R")
source("R/09_enrichment.R")
source("R/10_metaboanalyst_integration.R")
source("R/11_random_forest.R")

# =============================================================================
# Pipeline Definition
# =============================================================================

list(
  # ---------------------------------------------------------------------------
  # Configuration
  # ---------------------------------------------------------------------------
  tar_target(
    config,
    load_config(Sys.getenv("PIPELINE_CONFIG", unset = "config.yml"))
  ),

  # ---------------------------------------------------------------------------
  # Data Ingestion & Validation
  # ---------------------------------------------------------------------------
  tar_target(
    ingested_data,
    ingest_data(config)
  ),

  # ---------------------------------------------------------------------------
  # Initial Quality Control
  # ---------------------------------------------------------------------------
  tar_target(
    initial_qc,
    run_initial_qc(ingested_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Feature Filtering
  # ---------------------------------------------------------------------------
  tar_target(
    filtered_data,
    filter_features(ingested_data, initial_qc, config)
  ),

  # ---------------------------------------------------------------------------
  # Normalization & Transformation
  # ---------------------------------------------------------------------------
  tar_target(
    normalized_data,
    normalize_data(filtered_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Batch/Drift Correction
  # ---------------------------------------------------------------------------
  tar_target(
    batch_corrected_data,
    correct_batch_effects(normalized_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Missing Value Imputation
  # ---------------------------------------------------------------------------
  tar_target(
    imputed_data,
    impute_missing_values(batch_corrected_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Exploratory Multivariate Analysis
  # ---------------------------------------------------------------------------
  tar_target(
    exploratory_results,
    run_exploratory_analysis(imputed_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Differential Abundance Analysis
  # ---------------------------------------------------------------------------
  tar_target(
    de_results,
    run_differential_analysis(imputed_data, exploratory_results, ingested_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Enrichment Analysis
  # ---------------------------------------------------------------------------
  tar_target(
    enrichment_results,
    run_enrichment_analysis(de_results, imputed_data, ingested_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Optional: MetaboAnalyst integration (initialization & core analyses)
  # ---------------------------------------------------------------------------
  tar_target(
    metabo_mset,
    if (isTRUE(config$metaboanalyst$use_metaboanalyst)) {
      init_mset_from_table(config$input$feature_matrix, config)
    } else {
      NULL
    }
  ),

  tar_target(
    metabo_analyses,
    if (!is.null(metabo_mset)) metabo_preprocess(metabo_mset, config) else NULL
  ),

  tar_target(
    metabo_core_results,
    if (!is.null(metabo_mset) && isTRUE(config$metaboanalyst$run_core)) run_metabo_core_analyses(metabo_mset, config) else NULL
  ),

  # ---------------------------------------------------------------------------
  # Random Forest
  # ---------------------------------------------------------------------------
  tar_target(
    rf_results,
    if (isTRUE(config$rf$run_rf)) run_random_forest(imputed_data$matrix, imputed_data$metadata, config) else NULL
  ),

  # ---------------------------------------------------------------------------
  # Export normalized matrix as samples x features TSV
  # ---------------------------------------------------------------------------
  tar_target(
    normalized_samples_tsv,
    if (isTRUE(config$export$samples_by_features_tsv)) export_normalized_samples_tsv(normalized_data, config) else NULL
  ),

  # ---------------------------------------------------------------------------
  # Annotation Summary
  # ---------------------------------------------------------------------------
  tar_target(
    annotation_summary,
    create_annotation_summary(ingested_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Compile All Results
  # ---------------------------------------------------------------------------
  tar_target(
    analysis_results,
    list(
      config = config,
      ingested = ingested_data,
      initial_qc = initial_qc,
      filtered = filtered_data,
      normalized = normalized_data,
      batch_corrected = batch_corrected_data,
      imputed = imputed_data,
      exploratory = exploratory_results,
      differential = de_results,
      enrichment = enrichment_results,
      annotation_summary = annotation_summary
    )
  ),

  # ---------------------------------------------------------------------------
  # Generate Report (optional; disabled with DISABLE_REPORTS=1)
  # ---------------------------------------------------------------------------
  if (Sys.getenv("DISABLE_REPORTS", unset = "0") != "1") {
    tar_render(
      report,
      "reports/analysis_report.Rmd",
      output_dir = "outputs/report",
      params = list(
        results = analysis_results
      )
    )
  } else {
    message("Report generation targets disabled (DISABLE_REPORTS=1)")
  }
)
