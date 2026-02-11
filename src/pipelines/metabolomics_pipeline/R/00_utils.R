# =============================================================================
# 00_utils.R — Config loading, save helpers, logging
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

load_config <- function(config_path = "config.yml") {
  if (!file.exists(config_path))
    stop("Configuration file not found: ", config_path)
  config <- yaml::read_yaml(config_path)
  config <- set_config_defaults(config)
  config
}

set_config_defaults <- function(config) {
  # input
  config$input$sample_id_column   <- config$input$sample_id_column   %||% "SampleName"
  config$input$annotation_columns <- config$input$annotation_columns %||%
    c("Molecule", "HMDB", "SMILES", "KEGG", "CAS")
  config$input$feature_sheet      <- config$input$feature_sheet      %||% 1
  config$input$metadata_sheet     <- config$input$metadata_sheet     %||% "SampleSheet"

  # design
  config$design$condition_column <- config$design$condition_column %||% "condition"
  config$design$reference_level  <- config$design$reference_level  %||% NULL

  # normalization
  if (is.null(config$normalization)) config$normalization <- list()
  config$normalization$method      <- config$normalization$method      %||% "metaboanalyst"
  config$normalization$row_norm    <- config$normalization$row_norm    %||% "NULL"
  config$normalization$trans_norm  <- config$normalization$trans_norm  %||% "LogNorm"
  config$normalization$scale_norm  <- config$normalization$scale_norm  %||% "MeanCenter"
  config$normalization$pseudocount <- config$normalization$pseudocount %||% 1

  # qc / exploratory
  if (is.null(config$qc)) config$qc <- list()
  config$qc$run_pca              <- config$qc$run_pca              %||% TRUE
  config$qc$run_correlation      <- config$qc$run_correlation      %||% TRUE
  config$qc$run_sample_heatmap   <- config$qc$run_sample_heatmap   %||% TRUE
  config$qc$top_n_variable       <- config$qc$top_n_variable       %||% 50

  # differential
  if (is.null(config$differential)) config$differential <- list()
  config$differential$method             <- config$differential$method             %||% "limma"
  config$differential$adj_pvalue_threshold <- config$differential$adj_pvalue_threshold %||% 0.05
  config$differential$log2fc_threshold   <- config$differential$log2fc_threshold   %||% 1

  # random forest
  if (is.null(config$rf)) config$rf <- list()
  config$rf$run_rf     <- config$rf$run_rf     %||% TRUE
  config$rf$n_trees    <- config$rf$n_trees    %||% 500
  config$rf$importance <- config$rf$importance  %||% "permutation"
  config$rf$top_n      <- config$rf$top_n      %||% 20
  config$rf$seed       <- config$rf$seed       %||% 1234

  # plsda
  if (is.null(config$plsda)) config$plsda <- list()
  config$plsda$run_plsda <- config$plsda$run_plsda %||% TRUE
  config$plsda$vip_top_n <- config$plsda$vip_top_n %||% 15
  config$plsda$colors    <- config$plsda$colors    %||%
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")

  # enrichment
  if (is.null(config$enrichment)) config$enrichment <- list()
  config$enrichment$run_enrichment <- config$enrichment$run_enrichment %||% FALSE
  config$enrichment$organism       <- config$enrichment$organism       %||% "generic"
  config$enrichment$library        <- config$enrichment$library        %||% "kegg_compound"
  config$enrichment$method         <- config$enrichment$method         %||% "globaltest"

  # output
  if (is.null(config$output)) config$output <- list()
  config$output$output_dir <- config$output$output_dir %||% "outputs"
  config$output$fig_width  <- config$output$fig_width  %||% 10
  config$output$fig_height <- config$output$fig_height %||% 8
  config$output$plot_dpi   <- config$output$plot_dpi   %||% 300

  config
}

create_output_dirs <- function(config) {
  base <- config$output$output_dir
  dirs <- file.path(base, c("tables", "plots", "qc"))
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  invisible(dirs)
}

save_table <- function(df, filename, config, subdir = "tables") {
  path <- file.path(config$output$output_dir, subdir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, path)
  log_message("Saved table: ", path)
  path
}

save_plot <- function(plot, filename, config,
                      width = NULL, height = NULL, subdir = "plots") {
  width  <- width  %||% config$output$fig_width
  height <- height %||% config$output$fig_height
  dpi    <- config$output$plot_dpi
  path   <- file.path(config$output$output_dir, subdir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  log_message("Saved plot: ", path)
  path
}

log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste(..., collapse = ""))
  message(msg)
}
