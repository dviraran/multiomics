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

      cran_missing <- setdiff(missing, bioc_packages)
      bioc_missing <- intersect(missing, bioc_packages)

      if (length(cran_missing) > 0) {
        message("Installing CRAN packages: ", paste(cran_missing, collapse = ", "))
        install.packages(cran_missing)
      }

      if (length(bioc_missing) > 0) {
        message("Installing Bioconductor packages: ", paste(bioc_missing, collapse = ", "))
        BiocManager::install(bioc_missing, ask = FALSE)
      }

      # Re-check
      still_missing <- missing[!sapply(missing, requireNamespace, quietly = TRUE)]
      if (length(still_missing) > 0) {
        stop("Failed to install: ", paste(still_missing, collapse = ", "))
      }
      message("All packages installed successfully!")
    } else {
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

# =============================================================================
# Pipeline Definition
# =============================================================================

list(
  # ---------------------------------------------------------------------------
  # Configuration
  # ---------------------------------------------------------------------------
  tar_target(
    config,
    load_config("config.yml")
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
  # Generate Report
  # ---------------------------------------------------------------------------
  tar_render(
    report,
    "reports/analysis_report.Rmd",
    output_dir = "outputs/report",
    params = list(
      results = analysis_results
    )
  )
)
