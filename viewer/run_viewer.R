# =============================================================================
# Multi-Omics Viewer - Launcher Script
# =============================================================================
# Run this script to start the Multi-Omics Viewer.
#
# Usage:
#   source("run_viewer.R")
#   run_viewer("path/to/data")  # Specify data directory
#   run_viewer()                # Uses ./data by default
#
# Prerequisites:
#   - Run install_viewer.R first (one-time setup)
#   - Data files should be in the specified data folder
# =============================================================================

#' Launch the Multi-Omics Viewer
#'
#' @param data_dir Path to the data directory containing pipeline outputs.
#'   The directory should contain subdirectories like rnaseq/, proteomics/, etc.
#'   Defaults to "data" in the current directory.
#' @param port Port to run the Shiny app on (default: random available port)
#' @param launch_browser Whether to open browser automatically (default: TRUE)
#'
#' @examples
#' run_viewer("path/to/my_analysis/outputs")
#' run_viewer()  # Uses ./data
run_viewer <- function(data_dir = NULL, port = NULL, launch_browser = TRUE) {

    cat("
================================================================================
                     Multi-Omics Viewer
================================================================================
")

    # Determine script location and shiny_app directory
    script_dir <- tryCatch({
        dirname(sys.frame(1)$ofile)
    }, error = function(e) {
        getwd()
    })

    # Find shiny_app directory
    shiny_dir <- file.path(script_dir, "shiny_app")
    if (!dir.exists(shiny_dir)) {
        shiny_dir <- file.path(getwd(), "shiny_app")
    }
    if (!dir.exists(shiny_dir)) {
        shiny_dir <- file.path(getwd(), "viewer", "shiny_app")
    }

    if (!dir.exists(shiny_dir)) {
        stop("Cannot find 'shiny_app' directory. ",
             "Make sure you are running from the viewer directory or multiomics root.")
    }

    # Determine data directory
    if (is.null(data_dir)) {
        # Default locations to check
        candidates <- c(
            file.path(getwd(), "data"),
            file.path(script_dir, "data"),
            file.path(dirname(script_dir), "data")
        )
        for (candidate in candidates) {
            if (dir.exists(candidate)) {
                data_dir <- candidate
                break
            }
        }
        if (is.null(data_dir)) {
            data_dir <- file.path(getwd(), "data")
        }
    }

    # Convert to absolute path
    data_dir <- normalizePath(data_dir, mustWork = FALSE)

    cat("Shiny app location:", shiny_dir, "\n")
    cat("Data directory:", data_dir, "\n\n")

    # Check if data directory exists
    if (!dir.exists(data_dir)) {
        cat("WARNING: Data directory does not exist!\n")
        cat("  Path:", data_dir, "\n")
        cat("  The viewer will start but no data will be available.\n\n")
    } else {
        # Scan for available data
        cat("Scanning for data...\n")
        found_data <- FALSE

        for (dtype in c("rnaseq", "proteomics", "metabolomics", "multiomics")) {
            dtype_dir <- file.path(data_dir, dtype)
            if (dir.exists(dtype_dir)) {
                files <- list.files(dtype_dir, pattern = "\\.(csv|rds)$", recursive = TRUE)
                if (length(files) > 0) {
                    cat("  ", toupper(dtype), ": Found", length(files), "data files\n")
                    found_data <- TRUE
                } else {
                    cat("  ", toupper(dtype), ": Directory exists but no CSV/RDS files found\n")
                }
            }
        }

        if (!found_data) {
            cat("\nWARNING: No data files found!\n")
            cat("Expected structure:\n")
            cat("  ", data_dir, "/\n")
            cat("    rnaseq/       (with .csv files)\n")
            cat("    proteomics/   (with .csv files)\n")
            cat("    metabolomics/ (with .csv files)\n")
            cat("    multiomics/   (with .csv files)\n\n")
        }
        cat("\n")
    }

    # Check if required packages are installed
    required_packages <- c("shiny", "bslib", "plotly", "DT", "dplyr", "ggplot2")
    missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

    if (length(missing_packages) > 0) {
        cat("Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
        cat("Please run: source('install_viewer.R')\n")
        stop("Required packages not installed.")
    }

    # Set the data directory as an environment variable for global.R to use
    Sys.setenv(MULTIOMICS_DATA_DIR = data_dir)

    # Launch the application
    cat("Starting viewer...\n")
    cat("Press Ctrl+C (or Esc in RStudio) to stop the server.\n")
    cat("================================================================================\n\n")

    library(shiny)

    if (!is.null(port)) {
        runApp(shiny_dir, launch.browser = launch_browser, port = port)
    } else {
        runApp(shiny_dir, launch.browser = launch_browser)
    }
}

# If this script is sourced directly (not via run_viewer function),
# provide instructions
if (!interactive()) {
    cat("Multi-Omics Viewer loaded.\n")
    cat("Usage: run_viewer('path/to/data')\n")
} else {
    cat("Multi-Omics Viewer loaded.\n")
    cat("Usage: run_viewer('path/to/data')\n")
    cat("       run_viewer()  # uses ./data\n")
}
