# =============================================================================
# Multi-Omics Integration Pipeline
# =============================================================================
# A {targets} pipeline for integrating transcriptomics, proteomics, and metabolomics
#
# Usage:
#   targets::tar_make()       # Run the pipeline
#   targets::tar_visnetwork() # Visualize the pipeline
#   targets::tar_read(name)   # Read a target result
#
# =============================================================================

# =============================================================================
# Package check - prompt to install missing packages
# =============================================================================
check_and_install_packages <- function(packages) {
  # Separate CRAN and Bioconductor packages
  bioc_packages <- c(
    "DESeq2", "limma", "edgeR", "SummarizedExperiment", "MultiAssayExperiment",
    "S4Vectors", "MOFA2", "clusterProfiler", "org.Hs.eg.db", "fgsea",
    "ComplexHeatmap"
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
  "tidyverse", "yaml", "tibble",
  "DESeq2", "limma", "edgeR",
  "SummarizedExperiment", "MultiAssayExperiment", "S4Vectors",
  "MOFA2", "mixOmics", "SNFtool",
  "clusterProfiler", "org.Hs.eg.db", "fgsea",
  "ggplot2", "patchwork", "ComplexHeatmap", "circlize"
)

# Check packages before loading targets
check_and_install_packages(required_packages)

library(targets)
library(tarchetypes)

# Source all R functions
source("R/00_utils.R")
source("R/01_ingestion.R")
source("R/02_preprocessing.R")
source("R/03_mapping.R")
source("R/04_harmonize.R")
source("R/05_feature_selection.R")
source("R/06_mofa.R")
source("R/07_diablo.R")
source("R/08_snf.R")
source("R/05b_foundational_correlations.R")
source("R/05c_mechanistic_inference.R")
source("R/09_concordance.R")
source("R/09b_integration_consensus.R")  # NEW: Method comparison
source("R/09c_stability_analysis.R")     # NEW: Bootstrap stability
source("R/10_enrichment.R")
source("R/11_commentary.R")

# Set target options
tar_option_set(
  packages = c(
    "tidyverse", "yaml", "tibble",
    # Bioconductor (optional)
    "DESeq2", "limma", "edgeR",
    "SummarizedExperiment", "MultiAssayExperiment", "S4Vectors",
    # Integration methods (optional)
    "MOFA2", "mixOmics", "SNFtool",
    # Enrichment (optional)
    "clusterProfiler", "org.Hs.eg.db", "fgsea",
    # Visualization
    "ggplot2", "patchwork", "ComplexHeatmap", "circlize"
  ),
  error = "continue"  # Continue pipeline even if some targets fail
)

# =============================================================================
# Define Pipeline
# =============================================================================

list(
  # ---------------------------------------------------------------------------
  # Configuration
  # ---------------------------------------------------------------------------
  tar_target(
    name = config,
    command = load_config(Sys.getenv("PIPELINE_CONFIG", unset = "config.yml"))
  ),

  # ---------------------------------------------------------------------------
  # Data Ingestion
  # ---------------------------------------------------------------------------
  tar_target(
    name = ingested_data,
    command = ingest_all_data(config)
  ),

  # ---------------------------------------------------------------------------
  # Per-Omics Preprocessing
  # ---------------------------------------------------------------------------
  tar_target(
    name = preprocessed_data,
    command = preprocess_all_omics(ingested_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Identifier Harmonization
  # ---------------------------------------------------------------------------
  tar_target(
    name = harmonized_data,
    command = harmonize_identifiers(preprocessed_data, config)
  ),

  # ---------------------------------------------------------------------------
  # MultiAssayExperiment Creation
  # ---------------------------------------------------------------------------
  tar_target(
    name = mae_data,
    command = create_multiassay_experiment(harmonized_data, config)
  ),

  # MAE Summary
  tar_target(
    name = mae_summary,
    command = summarize_mae(mae_data, config)
  ),

  # ---------------------------------------------------------------------------
  # Feature Selection for Integration
  # ---------------------------------------------------------------------------
  tar_target(
    name = feature_data,
    command = select_features_for_integration(mae_data, config)
  ),

  # Feature overlap analysis
  tar_target(
    name = feature_overlap,
    command = get_feature_overlap(mae_data, config)
  ),

  # ---------------------------------------------------------------------------
  # FOUNDATIONAL CROSS-OMICS ANALYSIS (runs BEFORE integration methods)
  # ---------------------------------------------------------------------------
  # These analyses establish basic understanding of cross-omics relationships
  # before running advanced integration methods like MOFA/DIABLO/SNF

  # Cross-omics correlations, sample concordance, pathway overlap
  tar_target(
    name = foundational_results,
    command = {
      fc <- config$foundational %||% list()
      if (fc$run_foundational %||% TRUE) {
        run_foundational_analysis(mae_data, config)
      } else {
        log_message("Foundational analysis disabled in config. Skipping.")
        NULL
      }
    }
  ),

  # Mechanistic inference (regulatory networks, mediation, TF activity)
  tar_target(
    name = mechanistic_results,
    command = {
      mc <- config$mechanistic %||% list()
      if (mc$run_mechanistic %||% TRUE) {
        run_mechanistic_analysis(mae_data, foundational_results, config)
      } else {
        log_message("Mechanistic analysis disabled in config. Skipping.")
        NULL
      }
    }
  ),

  # ---------------------------------------------------------------------------
  # Integration Methods
  # ---------------------------------------------------------------------------

  # MOFA2 Integration
  tar_target(
    name = mofa_results,
    command = {
      if ("MOFA2" %in% config$integration$methods) {
        run_mofa2_integration(feature_data, config)
      } else {
        log_message("MOFA2 not in integration methods, skipping")
        NULL
      }
    }
  ),

  # MOFA factor associations
  tar_target(
    name = mofa_associations,
    command = {
      if (!is.null(mofa_results)) {
        test_mofa_factor_associations(mofa_results$results, mae_data$metadata, config)
      } else {
        NULL
      }
    }
  ),

  # DIABLO Integration
  tar_target(
    name = diablo_results,
    command = {
      if ("DIABLO" %in% config$integration$methods) {
        run_diablo_integration(feature_data, config)
      } else {
        log_message("DIABLO not in integration methods, skipping")
        NULL
      }
    }
  ),

  # SNF Integration (optional)
  tar_target(
    name = snf_results,
    command = {
      if ("SNF" %in% config$integration$methods) {
        run_snf_integration(feature_data, config)
      } else {
        log_message("SNF not in integration methods, skipping")
        NULL
      }
    }
  ),

  # Combined integration results
  tar_target(
    name = integration_results,
    command = list(
      mofa = mofa_results,
      diablo = diablo_results,
      snf = snf_results
    )
  ),

  # ---------------------------------------------------------------------------
  # Cross-Omics Concordance
  # ---------------------------------------------------------------------------
  tar_target(
    name = concordance_results,
    command = run_concordance_analysis(mae_data, integration_results, config)
  ),

  # ---------------------------------------------------------------------------
  # Multi-Omics Enrichment
  # ---------------------------------------------------------------------------
  tar_target(
    name = enrichment_results,
    command = run_multiomics_enrichment(mae_data, integration_results, config)
  ),

  # ---------------------------------------------------------------------------
  # Integration Consensus Analysis (NEW)
  # ---------------------------------------------------------------------------
  tar_target(
    name = consensus_results,
    command = {
      cc <- config$consensus %||% list()
      if (cc$compare_methods %||% TRUE) {
        run_integration_consensus(
          integration_results = integration_results,
          mae_data = mae_data,
          config = config
        )
      } else {
        log_message("Consensus analysis disabled in config. Skipping.")
        NULL
      }
    }
  ),

  # ---------------------------------------------------------------------------
  # Stability Analysis (NEW)
  # ---------------------------------------------------------------------------
  tar_target(
    name = stability_results,
    command = {
      sc <- config$stability %||% list()
      if (sc$run_stability %||% TRUE) {
        run_stability_analysis(
          mae_data = mae_data,
          feature_data = feature_data,
          integration_results = integration_results,
          config = config
        )
      } else {
        log_message("Stability analysis disabled in config. Skipping.")
        NULL
      }
    }
  ),

  # ---------------------------------------------------------------------------
  # Figure Commentary Generation
  # ---------------------------------------------------------------------------

  # Build table of all figures with metadata
  tar_target(
    name = figures_tbl,
    command = build_figures_table(
      mae_data = mae_data,
      foundational_results = foundational_results,
      mechanistic_results = mechanistic_results,
      integration_results = integration_results,
      concordance_results = concordance_results,
      enrichment_results = enrichment_results,
      consensus_results = consensus_results,
      stability_results = stability_results,
      config = config
    )
  ),

  # Generate commentary for all figures
  tar_target(
    name = commentary_tbl,
    command = generate_all_commentary(
      figures_tbl = figures_tbl,
      config = config,
      mae_data = mae_data,
      foundational_results = foundational_results,
      mechanistic_results = mechanistic_results,
      integration_results = integration_results,
      concordance_results = concordance_results,
      consensus_results = consensus_results,
      stability_results = stability_results,
      output_dir = file.path(config$output$output_dir, "commentary")
    )
  ),

  # Save figures table
  tar_target(
    name = figures_tbl_csv,
    command = {
      out_dir <- file.path(config$output$output_dir, "commentary")
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      out_path <- file.path(out_dir, "figures_metadata.csv")
      write_csv(figures_tbl, out_path)
      out_path
    },
    format = "file"
  ),

  # ---------------------------------------------------------------------------
  # Final Summary
  # ---------------------------------------------------------------------------
  tar_target(
    name = pipeline_summary,
    command = {
      log_message("=== Pipeline Complete ===")

      summary <- list(
        config = config,
        n_omics = length(names(mae_data$harmonized_omics)),
        omics_present = names(mae_data$harmonized_omics),
        n_samples = length(mae_data$common_samples),
        integration_methods_run = names(integration_results)[!sapply(integration_results, is.null)],
        mae_summary = mae_summary,
        feature_selection = feature_data$selection_summary,
        foundational_summary = if (!is.null(foundational_results)) foundational_results$summary else NULL,
        mechanistic_summary = if (!is.null(mechanistic_results)) mechanistic_results$summary else NULL,
        consensus_summary = if (!is.null(consensus_results)) consensus_results$summary else NULL,
        stability_summary = if (!is.null(stability_results)) stability_results$summary else NULL,
        commentary = commentary_tbl
      )

      # Save summary
      summary_file <- file.path(config$output$output_dir, "tables", "pipeline_summary.rds")
      saveRDS(summary, summary_file)

      log_message("Summary saved to: ", summary_file)
      log_message("Omics integrated: ", paste(summary$omics_present, collapse = ", "))
      log_message("Samples: ", summary$n_samples)
      log_message("Methods run: ", paste(summary$integration_methods_run, collapse = ", "))

      summary
    }
  ),

  # ---------------------------------------------------------------------------
  # Report Generation
  # ---------------------------------------------------------------------------
  tar_render(
    name = report,
    path = "reports/analysis_report.Rmd",
    output_dir = file.path("outputs", "report"),
    params = list(
      config = config,
      mae_data = mae_data,
      integration_results = integration_results,
      concordance_results = concordance_results,
      enrichment_results = enrichment_results,
      commentary_tbl = commentary_tbl
    )
  )
)
