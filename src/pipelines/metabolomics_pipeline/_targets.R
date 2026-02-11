# =============================================================================
# _targets.R — Metabolomics Pipeline DAG
# =============================================================================
# Run with:  targets::tar_make()
# Visualize: targets::tar_visnetwork()
# =============================================================================

library(targets)
library(tarchetypes)

options(tidyverse.quiet = TRUE)

tar_option_set(
  packages = c(
    "readr", "readxl", "yaml", "tibble", "dplyr",
    "ggplot2", "ggrepel", "pheatmap",
    "limma", "knitr", "rmarkdown", "DT"
  ),
  format = "rds",
  error = "continue"
)

# Source R modules
source("R/00_utils.R")
source("R/01_data_ingestion.R")
source("R/02_normalization.R")
source("R/03_exploratory.R")
source("R/04_differential.R")
source("R/05_random_forest.R")
source("R/06_export.R")
source("R/07_enrichment.R")
source("R/08_plsda.R")

# =============================================================================
# Pipeline DAG
# =============================================================================
#
#   config → ingested_data → normalized_data ─┬─→ exploratory_results
#                                              ├─→ de_results ──────┬─→ enrichment_results
#                                              ├─→ rf_results ──────┼─→ export
#                                              └─→ plsda_results ───┘
#
# =============================================================================

list(
  # -- Configuration ----------------------------------------------------------
  tar_target(
    config,
    load_config(Sys.getenv("PIPELINE_CONFIG", unset = "config.yml"))
  ),

  # -- Data Ingestion ---------------------------------------------------------
  tar_target(
    ingested_data,
    ingest_data(config)
  ),

  # -- Normalization ----------------------------------------------------------
  tar_target(
    normalized_data,
    normalize_data(ingested_data, config)
  ),

  # -- Exploratory Analysis ---------------------------------------------------
  tar_target(
    exploratory_results,
    run_exploratory(normalized_data, config)
  ),

  # -- Differential Expression ------------------------------------------------
  tar_target(
    de_results,
    run_differential(normalized_data, ingested_data, config)
  ),

  # -- Random Forest ----------------------------------------------------------
  tar_target(
    rf_results,
    run_random_forest(normalized_data, config)
  ),

  # -- PLS-DA -----------------------------------------------------------------
  tar_target(
    plsda_results,
    run_plsda(normalized_data, ingested_data, config)
  ),

  # -- Enrichment Analysis ----------------------------------------------------
  tar_target(
    enrichment_results,
    run_enrichment(de_results, normalized_data, ingested_data, config)
  ),

  # -- ssGSEA Pathway Enrichment -----------------------------------------------
  tar_target(
    ssgsea_results,
    run_ssgsea_enrichment(normalized_data, ingested_data, config)
  ),

  # -- Export -----------------------------------------------------------------
  tar_target(
    export,
    export_results(de_results, normalized_data, ingested_data, config)
  ),

  # -- HTML Report -------------------------------------------------------------
  tar_target(
    report,
    {
      # Reference all upstream targets so report runs last
      list(
        config, ingested_data, normalized_data,
        exploratory_results, de_results, rf_results,
        plsda_results, enrichment_results, ssgsea_results, export
      )

      # Collect all generated plot paths to ensure they exist
      plot_paths <- c()

      # Exploratory plots
      if (!is.null(exploratory_results)) {
        if (!is.null(exploratory_results$pca$plots)) {
          plot_paths <- c(plot_paths,
            file.path(config$output$output_dir, "qc", "pca_scores.png"),
            file.path(config$output$output_dir, "qc", "pca_scree.png")
          )
        }
        if (!is.null(exploratory_results$correlation)) {
          plot_paths <- c(plot_paths, exploratory_results$correlation$path)
        }
        if (!is.null(exploratory_results$sample_heatmap)) {
          plot_paths <- c(plot_paths, exploratory_results$sample_heatmap$path)
        }
      }

      # Differential analysis plots
      if (!is.null(de_results)) {
        plot_paths <- c(plot_paths,
          file.path(config$output$output_dir, "plots", "volcano_plot.png"),
          file.path(config$output$output_dir, "plots", "pvalue_histogram.png")
        )
      }

      # Random forest plot
      if (!is.null(rf_results)) {
        plot_paths <- c(plot_paths,
          file.path(config$output$output_dir, "plots", "rf_importance.png")
        )
      }

      # PLS-DA plots
      if (!is.null(plsda_results)) {
        plot_paths <- c(plot_paths,
          file.path(config$output$output_dir, "plots", "plsda_scores.png"),
          file.path(config$output$output_dir, "plots", "plsda_vip_scores.png")
        )
      }

      # Enrichment plot
      if (!is.null(enrichment_results) &&
          !is.null(enrichment_results$plot_path)) {
        plot_paths <- c(plot_paths, enrichment_results$plot_path)
      }

      # Normalization boxplot
      plot_paths <- c(plot_paths,
        file.path(config$output$output_dir, "qc", "normalization_boxplot.png")
      )

      # Verify plots exist (log warnings for missing plots)
      for (p in plot_paths) {
        if (!file.exists(p)) {
          warning("Plot file not found: ", p)
        }
      }

      # Render HTML report
      out_file <- file.path(
        getwd(), config$output$output_dir,
        "analysis_report.html"
      )
      rmarkdown::render(
        "reports/analysis_report.Rmd",
        output_file = out_file,
        params = list(
          config_path = Sys.getenv("PIPELINE_CONFIG", "config.yml"),
          pipeline_root = getwd(),
          config = config,
          ingested_data = ingested_data,
          normalized_data = normalized_data,
          exploratory_results = exploratory_results,
          de_results = de_results,
          rf_results = rf_results,
          plsda_results = plsda_results,
          enrichment_results = enrichment_results,
          ssgsea_results = ssgsea_results,
          export = export,
          plot_paths = plot_paths
        ),
        quiet = FALSE
      )

      log_message("HTML report generated: ", out_file)
      log_message("Report includes ", length(plot_paths), " plots")

      out_file
    }
  )
)
