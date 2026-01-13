# =============================================================================
# Utility Functions for Proteomics Pipeline
# =============================================================================

#' Load and validate configuration file
#'
#' @param config_path Path to config.yml file
#' @return List containing all configuration parameters
load_config <- function(config_path = "config.yml") {
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path,
         "\nPlease create a config.yml file. See README for details.")
  }

  config <- yaml::read_yaml(config_path)

  # Validate required fields
  required_fields <- c("input", "design", "processing", "filtering")
  missing <- setdiff(required_fields, names(config))
  if (length(missing) > 0) {
    stop("Missing required configuration sections: ", paste(missing, collapse = ", "))
  }

  # Set defaults for optional fields
  config <- set_config_defaults(config)

  return(config)
}

#' Set default values for optional configuration parameters
#'
#' @param config Configuration list
#' @return Configuration list with defaults applied
set_config_defaults <- function(config) {
  # Input defaults
  if (is.null(config$input$sample_id_column)) {
    config$input$sample_id_column <- "sample_id"
  }
  if (is.null(config$input$feature_level)) {
    config$input$feature_level <- "protein"
  }

  # Processing defaults
  if (is.null(config$processing$zeros_as_na)) {
    config$processing$zeros_as_na <- TRUE
  }
  if (is.null(config$processing$log_transform)) {
    config$processing$log_transform <- TRUE
  }
  if (is.null(config$processing$normalization_method)) {
    config$processing$normalization_method <- "vsn"
  }

  # Filtering defaults
  if (is.null(config$filtering$global_min_presence)) {
    config$filtering$global_min_presence <- 0.5
  }

  if (is.null(config$filtering$group_min_presence)) {
    config$filtering$group_min_presence <- 0.7
  }

  # Imputation defaults
  if (is.null(config$imputation)) {
    config$imputation <- list(method = "QRILC")
  }
  if (is.null(config$imputation$method)) {
    config$imputation$method <- "QRILC"
  }

  # QC defaults
  if (is.null(config$qc)) {
    config$qc <- list()
  }
  if (is.null(config$qc$n_pca_components)) {
    config$qc$n_pca_components <- 10
  }
  if (is.null(config$qc$outlier_sd_threshold)) {
    config$qc$outlier_sd_threshold <- 3
  }

  # Differential defaults
  if (is.null(config$differential)) {
    config$differential <- list()
  }
  if (is.null(config$differential$method)) {
    config$differential$method <- "limma"
  }
  if (is.null(config$differential$adj_pvalue_threshold)) {
    config$differential$adj_pvalue_threshold <- 0.05
  }
  if (is.null(config$differential$log2fc_threshold)) {
    config$differential$log2fc_threshold <- 1
  }

  # Pathway defaults
  if (is.null(config$pathway)) {
    config$pathway <- list(run_pathway_analysis = TRUE)
  }

  # Output defaults
  if (is.null(config$output)) {
    config$output <- list()
  }
  if (is.null(config$output$output_dir)) {
    config$output$output_dir <- "outputs"
  }

  return(config)
}

#' Create output directories if they don't exist
#'
#' @param config Configuration list
create_output_dirs <- function(config) {
  output_dir <- config$output$output_dir

  dirs <- c(
    output_dir,
    file.path(output_dir, "tables"),
    file.path(output_dir, "plots"),
    file.path(output_dir, "qc"),
    file.path(output_dir, "report")
  )

  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
    }
  }

  invisible(dirs)
}

#' Save data frame to CSV with consistent formatting
#'
#' @param df Data frame to save
#' @param filename Filename (will be placed in output_dir/tables/)
#' @param config Configuration list
#' @param subdir Subdirectory within output_dir (default: "tables")
#' @return Path to saved file
save_table <- function(df, filename, config, subdir = "tables") {
  output_path <- file.path(config$output$output_dir, subdir, filename)
  readr::write_csv(df, output_path)
  message("Saved: ", output_path)
  return(output_path)
}

#' Save plot to file
#'
#' @param plot ggplot object
#' @param filename Filename (will be placed in output_dir/plots/)
#' @param config Configuration list
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param subdir Subdirectory within output_dir (default: "plots")
#' @return Path to saved file
save_plot <- function(plot, filename, config, width = NULL, height = NULL, subdir = "plots") {
  if (is.null(width)) width <- config$output$fig_width %||% 10
  if (is.null(height)) height <- config$output$fig_height %||% 8
  dpi <- config$output$plot_dpi %||% 300

  output_path <- file.path(config$output$output_dir, subdir, filename)

  ggplot2::ggsave(
    output_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )

  message("Saved: ", output_path)
  return(output_path)
}

#' Parse UniProt accession from various formats
#'
#' Handles formats like:
#' - sp|P12345|PROTEIN_NAME
#' - tr|Q12345|PROTEIN_NAME
#' - P12345
#' - P12345-2 (isoforms)
#'
#' @param ids Character vector of protein identifiers
#' @return Character vector of cleaned UniProt accessions
parse_uniprot_accession <- function(ids) {
  # Pattern for UniProt accession (with optional isoform suffix)
  uniprot_pattern <- "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?"

  # Try to extract from sp|XXX|NAME or tr|XXX|NAME format
  extracted <- stringr::str_extract(ids, paste0("(sp|tr)\\|", uniprot_pattern, "\\|"))
  extracted <- stringr::str_replace_all(extracted, "(sp|tr)\\|", "")
  extracted <- stringr::str_replace_all(extracted, "\\|$", "")


  # For entries that didn't match sp/tr format, try direct match
  no_match <- is.na(extracted)
  if (any(no_match)) {
    direct_match <- stringr::str_extract(ids[no_match], uniprot_pattern)
    extracted[no_match] <- direct_match
  }

  return(extracted)
}

#' Calculate missingness statistics
#'
#' @param mat Numeric matrix (features x samples)
#' @return Data frame with missingness statistics per feature
calculate_missingness <- function(mat) {
  n_samples <- ncol(mat)

  data.frame(
    feature_id = rownames(mat),
    n_present = rowSums(!is.na(mat)),
    n_missing = rowSums(is.na(mat)),
    pct_present = rowSums(!is.na(mat)) / n_samples * 100,
    pct_missing = rowSums(is.na(mat)) / n_samples * 100,
    mean_intensity = rowMeans(mat, na.rm = TRUE),
    sd_intensity = apply(mat, 1, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Log message with timestamp
#'
#' @param ... Message components
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(msg)
}

#' Check if required packages are installed
#'
#' @param packages Character vector of package names
#' @return Invisible TRUE if all installed, error otherwise
check_packages <- function(packages) {
  missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nInstall with: install.packages(c('", paste(missing, collapse = "', '"), "'))",
         "\nOr for Bioconductor packages: BiocManager::install(c('", paste(missing, collapse = "', '"), "'))")
  }

  invisible(TRUE)
}

#' Null-coalescing operator
#'
#' @param x Primary value
#' @param y Fallback value
#' @return x if not NULL, otherwise y
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Format numbers for display
#'
#' @param x Numeric vector
#' @param digits Number of significant digits
#' @return Formatted character vector
format_number <- function(x, digits = 3) {
  format(signif(x, digits), scientific = FALSE)
}

#' Create a summary statistics table
#'
#' @param mat Numeric matrix
#' @param metadata Sample metadata data frame
#' @param group_col Column name for grouping
#' @return Summary data frame
summarize_by_group <- function(mat, metadata, group_col) {
  groups <- unique(metadata[[group_col]])

  results <- lapply(groups, function(g) {
    samples <- metadata[[config$input$sample_id_column]][metadata[[group_col]] == g]
    samples <- intersect(samples, colnames(mat))

    if (length(samples) == 0) return(NULL)

    sub_mat <- mat[, samples, drop = FALSE]

    data.frame(
      group = g,
      n_samples = length(samples),
      mean_intensity = mean(sub_mat, na.rm = TRUE),
      median_intensity = median(sub_mat, na.rm = TRUE),
      pct_missing = mean(is.na(sub_mat)) * 100,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}
