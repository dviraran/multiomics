# =============================================================================
# Pipeline Runner Helper Functions
# =============================================================================
# Source this file to get helper functions for running pipelines with custom
# config files from any directory.
#
# Usage:
#   source("path/to/pipeline_runner.R")
#
#   # Run rnaseq pipeline with custom config
#   run_rnaseq_pipeline("path/to/my_config.yml")
#
#   # Run any pipeline
#   run_omics_pipeline("rnaseq", "path/to/my_config.yml")
#
# =============================================================================

library(targets)

# Find the multiomics root directory
find_multiomics_root <- function() {
  # Try to get path from sourced file location first
  script_path <- tryCatch({
    # Works when file is being sourced
    dirname(sys.frame(1)$ofile)
  }, error = function(e) NULL)

  # Try common locations
  candidates <- c(
    Sys.getenv("MULTIOMICS_ROOT"),
    script_path,  # If sourced, use the script's directory
    ".",
    "..",
    "../.."
  )

  for (candidate in candidates) {
    if (is.null(candidate) || candidate == "") next
    candidate <- normalizePath(candidate, mustWork = FALSE)
    # Check for new structure (src/pipelines/)
    if (dir.exists(file.path(candidate, "src", "pipelines", "rnaseq_pipeline"))) {
      return(candidate)
    }
    # Fallback for old structure
    if (dir.exists(file.path(candidate, "rnaseq_pipeline"))) {
      return(candidate)
    }
  }

  stop("Cannot find multiomics root directory. Set MULTIOMICS_ROOT environment variable.")
}

#' Run an omics pipeline with a custom config file
#'
#' @param pipeline Which pipeline: "rnaseq", "proteomics", "metabolomics", or "multiomics"
#' @param config_file Path to config file (default: config.yml in pipeline directory)
#' @param clean Whether to clean the targets cache first (default: FALSE)
#' @return Invisibly returns the pipeline result
#'
#' @examples
#' # Run with default config
#' run_omics_pipeline("rnaseq")
#'
#' # Run with custom config
#' run_omics_pipeline("rnaseq", "my_analysis/rnaseq_config.yml")
#'
#' # Run with clean cache
#' run_omics_pipeline("rnaseq", clean = TRUE)
run_omics_pipeline <- function(pipeline, config_file = NULL, clean = FALSE) {
  valid_pipelines <- c("rnaseq", "proteomics", "metabolomics", "multiomics")
  if (!pipeline %in% valid_pipelines) {
    stop("Invalid pipeline. Choose from: ", paste(valid_pipelines, collapse = ", "))
  }

  root <- find_multiomics_root()

  # Check for new structure first (src/pipelines/)
  pipeline_dir <- file.path(root, "src", "pipelines", paste0(pipeline, "_pipeline"))
  if (!dir.exists(pipeline_dir)) {
    # Fallback to old structure
    pipeline_dir <- file.path(root, paste0(pipeline, "_pipeline"))
  }

  if (!dir.exists(pipeline_dir)) {
    stop("Pipeline directory not found: ", pipeline_dir)
  }

  # Convert config to absolute path if provided
  if (!is.null(config_file)) {
    if (!file.exists(config_file)) {
      stop("Config file not found: ", config_file)
    }
    config_file <- normalizePath(config_file)
  }

  cat("╔══════════════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║  Running %s pipeline%s\n",
              toupper(pipeline),
              paste0(rep(" ", 52 - nchar(pipeline)), collapse = "")), "║\n")
  cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

  cat("Pipeline directory:", pipeline_dir, "\n")
  if (!is.null(config_file)) {
    cat("Config file:", config_file, "\n")
  }
  cat("\n")

  # Save current state

  original_dir <- getwd()
  original_config <- Sys.getenv("PIPELINE_CONFIG")

  # Set up environment
  setwd(pipeline_dir)
  if (!is.null(config_file)) {
    Sys.setenv(PIPELINE_CONFIG = config_file)
  }

  tryCatch({
    if (clean) {
      cat("Cleaning targets cache...\n")
      tar_destroy(ask = FALSE)
    }

    cat("Starting pipeline...\n\n")
    start_time <- Sys.time()

    tar_make()

    end_time <- Sys.time()
    duration <- difftime(end_time, start_time, units = "mins")

    cat("\n✓ Pipeline completed in", round(as.numeric(duration), 2), "minutes\n")

    invisible(list(status = "success", duration = duration))

  }, error = function(e) {
    cat("\n✗ Pipeline error:", e$message, "\n")
    invisible(list(status = "error", error = e$message))

  }, finally = {
    # Restore state
    Sys.unsetenv("PIPELINE_CONFIG")
    if (original_config != "") {
      Sys.setenv(PIPELINE_CONFIG = original_config)
    }
    setwd(original_dir)
  })
}

# Convenience wrappers for each pipeline
run_rnaseq_pipeline <- function(config_file = NULL, clean = FALSE) {
  run_omics_pipeline("rnaseq", config_file, clean)
}

run_proteomics_pipeline <- function(config_file = NULL, clean = FALSE) {
  run_omics_pipeline("proteomics", config_file, clean)
}

run_metabolomics_pipeline <- function(config_file = NULL, clean = FALSE) {
  run_omics_pipeline("metabolomics", config_file, clean)
}

run_multiomics_pipeline <- function(config_file = NULL, clean = FALSE) {
  run_omics_pipeline("multiomics", config_file, clean)
}

cat("Pipeline runner loaded. Available functions:\n")
cat("  run_omics_pipeline(pipeline, config_file, clean)\n")
cat("  run_rnaseq_pipeline(config_file, clean)\n")
cat("  run_proteomics_pipeline(config_file, clean)\n")
cat("  run_metabolomics_pipeline(config_file, clean)\n")
cat("  run_multiomics_pipeline(config_file, clean)\n")
