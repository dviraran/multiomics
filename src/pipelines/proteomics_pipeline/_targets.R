# =============================================================================
# Proteomics Analysis Pipeline - Main targets file
# =============================================================================
# A comprehensive, modular pipeline for mass-spec proteomics analysis
# Supports LFQ/DIA workflows from MaxQuant, FragPipe, Spectronaut, DIA-NN, etc.
#
# Run with: targets::tar_make()
# Visualize: targets::tar_visnetwork()
# =============================================================================

# =============================================================================
# Package check - prompt to install missing packages
# =============================================================================
check_and_install_packages <- function(packages) {
  # Separate CRAN and Bioconductor packages
  bioc_packages <- c(
    "limma", "vsn", "fgsea", "clusterProfiler", "org.Hs.eg.db",
    "SummarizedExperiment", "BiocGenerics"
  )

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
  "ggplot2", "pheatmap", "RColorBrewer", "ggrepel",
  "limma", "vsn", "fgsea"
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

    # Statistics
    "limma",
    "vsn",

    # Enrichment
    "fgsea"
  ),
  format = "rds",
  error = "continue"  # Continue pipeline even if some targets fail
)

# Source all R functions
source("R/00_utils.R")
source("R/01_data_ingestion.R")
source("R/02_annotation.R")
source("R/03_filtering.R")
source("R/04_normalization.R")
source("R/05_imputation.R")
source("R/06_qc.R")
source("R/07_differential.R")
source("R/07b_ppi_networks.R")      # NEW: PPI network analysis
source("R/07c_advanced_stats.R")    # NEW: Advanced statistical methods
source("R/08_pathway.R")
source("R/09_commentary.R")

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
  # Feature Annotation
  # ---------------------------------------------------------------------------
  tar_target(
    annotation_data,
    annotate_features(ingested_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Contaminant Removal & Filtering
  # ---------------------------------------------------------------------------
  tar_target(
    filtered_data,
    filter_features(ingested_data, annotation_data, config)
  ),

  # Store pre-processed matrix for comparison
  tar_target(
    raw_matrix,
    ingested_data$matrix
  ),

  # ---------------------------------------------------------------------------
  # Normalization
  # ---------------------------------------------------------------------------
  tar_target(
    normalized_data,
    normalize_data(filtered_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Missing Value Imputation
  # ---------------------------------------------------------------------------
  tar_target(
    imputed_data,
    impute_missing_values(normalized_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Quality Control & Exploratory Analysis
  # ---------------------------------------------------------------------------
  tar_target(
    qc_results,
    run_qc_analysis(imputed_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Differential Abundance Analysis
  # ---------------------------------------------------------------------------
  tar_target(
    de_results,
    run_differential_analysis(imputed_data, qc_results, config)
  ),

  # ---------------------------------------------------------------------------
  # Pathway Analysis
  # ---------------------------------------------------------------------------
  tar_target(
    pathway_results,
    run_pathway_analysis(de_results, imputed_data, config)
  ),

  # ---------------------------------------------------------------------------
  # PPI Network Analysis (NEW)
  # ---------------------------------------------------------------------------
  tar_target(
    ppi_results,
    {
      ppi_config <- config$ppi_analysis %||% list()
      if (ppi_config$run_ppi %||% TRUE) {
        run_ppi_network_analysis(
          da_results = de_results$results[[1]],  # Use first contrast
          normalized_data = imputed_data$matrix,
          config = config
        )
      } else {
        NULL
      }
    }
  ),

  tar_target(
    ppi_network_csv,
    {
      if (!is.null(ppi_results) && !is.null(ppi_results$network$edges)) {
        out_path <- file.path(config$output$output_dir, "ppi_network.csv")
        write.csv(ppi_results$network$edges, out_path, row.names = FALSE)
        out_path
      } else {
        NULL
      }
    },
    format = "file"
  ),

  # ---------------------------------------------------------------------------
  # Advanced Statistical Analysis (NEW)
  # ---------------------------------------------------------------------------
  tar_target(
    advanced_stats,
    {
      stats_config <- config$advanced_stats %||% list()
      if (stats_config$run_advanced_stats %||% TRUE) {
        run_advanced_stats(
          normalized_data = imputed_data$matrix,
          metadata = imputed_data$metadata,
          da_results = de_results$results[[1]],
          config = config
        )
      } else {
        NULL
      }
    }
  ),

  # ---------------------------------------------------------------------------
  # Figure Commentary Generation
  # ---------------------------------------------------------------------------

  # Build table of all figures with metadata
  tar_target(
    figures_tbl,
    build_figures_table(
      qc_results = qc_results,
      de_results = de_results,
      pathway_results = pathway_results,
      ppi_results = ppi_results,
      advanced_stats = advanced_stats,
      config = config
    )
  ),

  # Generate commentary for all figures
  tar_target(
    commentary_tbl,
    generate_all_commentary(
      figures_tbl = figures_tbl,
      config = config,
      qc_results = qc_results,
      de_results = de_results,
      imputed_data = imputed_data,
      ppi_results = ppi_results,
      advanced_stats = advanced_stats,
      output_dir = file.path(config$output$output_dir, "commentary")
    )
  ),

  # Save figures table
  tar_target(
    figures_tbl_csv,
    {
      out_dir <- file.path(config$output$output_dir, "commentary")
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      out_path <- file.path(out_dir, "figures_metadata.csv")
      write_csv(figures_tbl, out_path)
      out_path
    },
    format = "file"
  ),

  # ---------------------------------------------------------------------------
  # Compile All Results
  # ---------------------------------------------------------------------------
  tar_target(
    analysis_results,
    list(
      config = config,
      ingested = ingested_data,
      annotations = annotation_data,
      filtered = filtered_data,
      normalized = normalized_data,
      imputed = imputed_data,
      qc = qc_results,
      differential = de_results,
      pathway = pathway_results,
      ppi = ppi_results,
      advanced_stats = advanced_stats,
      commentary = commentary_tbl
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
