## 00_utils.R - helper utilities for the pipeline
# Minimal safe helpers so `_targets.R` can be validated.

check_and_install_packages <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, require, character.only = TRUE))
}

load_config <- function(path = "config.yml") {
  if (!file.exists(path)) stop("Config file not found: ", path)
  yaml::read_yaml(path)
}

safe_write_rds <- function(object, path) {
  tryCatch({
    saveRDS(object, path)
    return(normalizePath(path))
  }, error = function(e) {
    stop("Failed to write RDS to ", path, ": ", e$message)
  })
}
# =============================================================================
# Utility Functions for Metabolomics Pipeline
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
  if (is.null(config$input$data_type)) {
    config$input$data_type <- "untargeted"
  }

  # Sample type defaults
  if (is.null(config$sample_types)) {
    config$sample_types <- list(
      sample_type_column = "sample_type",
      sample_value = "Sample",
      qc_value = "QC",
      blank_value = "Blank",
      assume_all_samples = TRUE
    )
  }

  # Processing defaults
  if (is.null(config$processing$zeros_as_na)) {
    config$processing$zeros_as_na <- TRUE
  }
  if (is.null(config$processing$transform)) {
    config$processing$transform <- "log2"
  }
  if (is.null(config$processing$pseudocount)) {
    config$processing$pseudocount <- 1
  }
  if (is.null(config$processing$normalization_method)) {
    config$processing$normalization_method <- "PQN"
  }
  if (is.null(config$processing$scaling_method)) {
    config$processing$scaling_method <- "pareto"
  }

  # Filtering defaults
  if (is.null(config$filtering$global_min_presence)) {
    config$filtering$global_min_presence <- 0.2
  }
  if (is.null(config$filtering$group_min_presence)) {
    config$filtering$group_min_presence <- 0.5
  }

  # Batch correction defaults
  if (is.null(config$batch_correction)) {
    config$batch_correction <- list(method = "none")
  }

  # Imputation defaults
  if (is.null(config$imputation)) {
    config$imputation <- list(method = "half_min")
  }

  # QC defaults
  if (is.null(config$qc)) {
    config$qc <- list(n_pca_components = 10, outlier_sd_threshold = 3)
  }

  # Differential defaults
  if (is.null(config$differential)) {
    config$differential <- list(
      method = "limma",
      adj_pvalue_threshold = 0.05,
      log2fc_threshold = 1
    )
  }

  # Multivariate defaults
  if (is.null(config$multivariate)) {
    config$multivariate <- list(run_plsda = FALSE)
  }

  # Enrichment defaults
  if (is.null(config$enrichment)) {
    config$enrichment <- list(run_enrichment = TRUE, methods = list("ora"))
  }

  # Output defaults
  if (is.null(config$output)) {
    config$output <- list(output_dir = "outputs")
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
    file.path(output_dir, "report"),
    file.path(output_dir, "logs")
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

  # Ensure target directory exists
  out_dir <- dirname(output_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(df, output_path)
  log_message("Saved: ", output_path)
  return(output_path)
}

#' Save plot to file
#'
#' @param plot ggplot object
#' @param filename Filename
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

  # Ensure target directory exists
  out_dir <- dirname(output_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    output_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )

  log_message("Saved: ", output_path)
  return(output_path)
}

#' Log message with timestamp
#'
#' @param ... Message components
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(msg)
}

#' Null-coalescing operator
#'
#' @param x Primary value
#' @param y Fallback value
#' @return x if not NULL, otherwise y
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Calculate coefficient of variation (CV) / RSD
#'
#' @param x Numeric vector
#' @param na.rm Remove NAs
#' @return CV as proportion (not percentage)
calc_cv <- function(x, na.rm = TRUE) {
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}

#' Calculate RSD (relative standard deviation) as percentage
#'
#' @param x Numeric vector
#' @param na.rm Remove NAs
#' @return RSD as percentage
calc_rsd <- function(x, na.rm = TRUE) {
  calc_cv(x, na.rm = na.rm) * 100
}

#' Sanitize sample names
#'
#' Remove trailing whitespace and problematic characters
#'
#' @param names Character vector of names
#' @return Sanitized names
sanitize_names <- function(names) {
  names <- trimws(names)
  names <- gsub("[^a-zA-Z0-9_.-]", "_", names)
  return(names)
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

#' Check if required packages are installed
#'
#' @param packages Character vector of package names
#' @return Invisible TRUE if all installed, error otherwise
check_packages <- function(packages) {
  missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nInstall with appropriate method (CRAN or BiocManager)")
  }

  invisible(TRUE)
}

#' Pareto scaling
#'
#' @param x Numeric vector
#' @return Pareto-scaled vector
pareto_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sqrt(sd(x, na.rm = TRUE))
}

#' Autoscaling (z-score)
#'
#' @param x Numeric vector
#' @return Autoscaled vector
autoscale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

#' Apply scaling to matrix (row-wise for features)
#'
#' @param mat Numeric matrix (features x samples)
#' @param method Scaling method: "pareto", "autoscale", "none"
#' @return Scaled matrix
apply_scaling <- function(mat, method = "pareto") {
  if (method == "none") {
    return(mat)
  }

  scale_func <- switch(method,
    "pareto" = pareto_scale,
    "autoscale" = autoscale,
    identity
  )

  t(apply(mat, 1, scale_func))
}

#' Format numbers for display
#'
#' @param x Numeric vector
#' @param digits Number of significant digits
#' @return Formatted character vector
format_number <- function(x, digits = 3) {
  format(signif(x, digits), scientific = FALSE)
}

#' Read GMT file for metabolite sets
#'
#' @param gmt_file Path to GMT file
#' @return Named list of sets
read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)

  sets <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      set_name <- parts[1]
      members <- parts[3:length(parts)]
      members <- members[members != ""]
      sets[[set_name]] <- members
    }
  }

  return(sets)
}
