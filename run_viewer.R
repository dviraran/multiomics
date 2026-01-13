# =============================================================================
# Multi-Omics Viewer - Launcher Script
# =============================================================================
# Run this script to start the Multi-Omics Viewer.
#
# Usage:
#   source("run_viewer.R")
#
# Prerequisites:
#   - Run install_viewer.R first (one-time setup)
#   - Data files should be in the data/ folder
# =============================================================================

cat("
================================================================================
                     Multi-Omics Viewer
================================================================================
Starting application...
")

# Check if required packages are installed
required_packages <- c("shiny", "bslib", "plotly", "DT", "dplyr", "ggplot2")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
    cat("\nMissing required packages:", paste(missing_packages, collapse = ", "), "\n")
    cat("Please run: source('install_viewer.R')\n")
    stop("Required packages not installed.")
}

# Check for data directory
data_dir <- file.path(getwd(), "data")
shiny_dir <- file.path(getwd(), "shiny_app")

if (!dir.exists(data_dir)) {
    cat("\nWarning: No 'data' directory found.\n")
    cat("The viewer will start but no data will be loaded.\n")
    cat("Expected data directory:", data_dir, "\n\n")
}

if (!dir.exists(shiny_dir)) {
    cat("\nError: Cannot find 'shiny_app' directory.\n")
    cat("Make sure you are running this script from the viewer root directory.\n")
    stop("shiny_app directory not found.")
}

# Launch the application
cat("Opening viewer in browser...\n")
cat("Press Ctrl+C (or Esc in RStudio) to stop the server.\n")
cat("================================================================================\n\n")

library(shiny)
runApp(shiny_dir, launch.browser = TRUE)
