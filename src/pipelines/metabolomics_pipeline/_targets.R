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
source("R/12_metabolite_mapping.R")

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
  # Metabolite Mapping (Optional / Automatic)
  # ---------------------------------------------------------------------------
  tar_target(
    mapped_pathways,
    map_metabolites_to_pathways(ingested_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Enrichment Analysis
  # ---------------------------------------------------------------------------
  tar_target(
    enrichment_results,
    run_enrichment_analysis(de_results, imputed_data, ingested_data, config, mapped_pathways)
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
  # Global Test QEA (MetaboAnalystR) ✅
  # Runs performGlobalTestQEA() when MetaboAnalyst integration is enabled and
  # QEA is requested in the config (config$metaboanalyst$run_qea)
  # ---------------------------------------------------------------------------
  tar_target(
    qea_results,
    if (!is.null(metabo_mset) && isTRUE(config$metaboanalyst$run_qea)) {
      # Prefer normalized matrix from MetaboAnalyst mSet, fallback to imputed matrix
      conc_mat <- extract_normalized_from_mset(metabo_mset, config)
      if (is.null(conc_mat) && !is.null(imputed_data$matrix)) conc_mat <- imputed_data$matrix

      if (is.null(conc_mat)) stop("No concentration matrix available for QEA (neither MetaboAnalyst mSet nor imputed_data$matrix present).")

      conc_df <- as.data.frame(conc_mat)
      conc_df <- tibble::rownames_to_column(conc_df, var = "Compound")

      group_col <- config$design$condition_column
      if (!group_col %in% colnames(ingested_data$metadata)) stop("Condition column '", group_col, "' not found in metadata for QEA.")
      group_vec <- ingested_data$metadata[[group_col]]

      # Save QEA outputs under pipeline output dir (outputs/qea)
      performGlobalTestQEA(
        conc_df,
        group_vec,
        norm_method = config$metaboanalyst$normalization_method %||% "SumNorm",
        output_dir = config$output$output_dir,
        save_outputs = TRUE
      )
    } else {
      NULL
    }
  ),

  # ---------------------------------------------------------------------------
  # QEA integration test (synthetic data)
  # Ensures performGlobalTestQEA runs end-to-end and saves expected outputs
  # ---------------------------------------------------------------------------
  tar_target(
    qea_integration_test,
    if (isTRUE(config$metaboanalyst$run_qea)) {
      set.seed(42)
      n_met <- 100
      n_samp <- 12
      met_names <- paste0("M", sprintf("%03d", 1:n_met))
      groups <- rep(c("Control", "Treatment"), each = n_samp / 2)
      mat <- matrix(rnorm(n_met * n_samp, mean = 5, sd = 2), nrow = n_met)
      # Inject missingness and small signal
      mat[sample(length(mat), size = floor(0.15 * length(mat)))] <- NA
      mat[1:10, groups == "Treatment"] <- mat[1:10, groups == "Treatment"] + 1.5
      df <- data.frame(Compound = met_names, mat, stringsAsFactors = FALSE)
      names(df)[-1] <- paste0("S", 1:n_samp)

      res <- performGlobalTestQEA(df, groups,
        norm_method = "SumNorm",
        output_dir = file.path(config$output$output_dir, "qea", "test"),
        save_outputs = TRUE
      )

      if (is.null(res$pathway_results_df)) stop("QEA integration test failed: no pathway results produced")
      if (!"p.adj" %in% colnames(res$pathway_results_df)) stop("QEA integration test failed: p.adj missing from results")
      if (length(res$saved_paths) == 0) stop("QEA integration test failed: no saved outputs detected")

      TRUE
    } else {
      TRUE
    }
  ),

  # ---------------------------------------------------------------------------
  # Archive QEA outputs (zip) and optionally upload to GitHub release
  # ---------------------------------------------------------------------------
  tar_target(
    qea_archive,
    {
      if (!isTRUE(config$metaboanalyst$run_qea)) {
        return(NULL)
      }

      qea_dir <- file.path(config$output$output_dir, "qea")
      if (!dir.exists(qea_dir)) {
        warning("QEA output directory not found: ", qea_dir)
        return(NULL)
      }

      files <- list.files(qea_dir, recursive = TRUE, full.names = TRUE)
      if (length(files) == 0) {
        warning("No files to archive in: ", qea_dir)
        return(NULL)
      }

      out_zip <- file.path(config$output$output_dir, paste0("qea_outputs_", format(Sys.Date(), "%Y%m%d"), ".zip"))
      # Create zip from within the qea directory so paths are relative
      owd <- getwd()
      on.exit(setwd(owd), add = TRUE)
      setwd(qea_dir)
      rel_files <- list.files(".", recursive = TRUE)

      # Use utils::zip for portability
      utils::zip(out_zip, files = rel_files)

      if (!file.exists(out_zip)) stop("Failed to create qea archive: ", out_zip)

      out_zip
    },
    format = "file"
  ),

  # tar_target(qea_upload, ...) removed by request to disable GitHub uploads
  # If you later want safe upload functionality, re-add a carefully gated target that requires explicit config.enable_upload = TRUE


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
      mapped_pathways = mapped_pathways,
      qea = qea_results,
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
