# =============================================================================
# Global Settings for Multi-Omics Viewer
# =============================================================================
# This file is sourced once when the app starts.
# It loads packages, sets up the theme, and scans for available data.
# LAZY LOADING: Data is NOT loaded here - only paths are discovered.

# Required packages
required_packages <- c(
    "shiny",
    "bslib",
    "plotly",
    "DT",
    "heatmaply",
    "visNetwork",
    "dplyr",
    "tidyr",
    "ggplot2",
    "scales",
    "viridis",
    "RColorBrewer",
    "ggrepel"
)

# Load all packages
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste0("Package '", pkg, "' is required but not installed.\n",
                    "Run source('install_viewer.R') first."))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Source utility functions
source("R/data_loaders.R")
source("R/plot_functions.R")
source("R/utils.R")

# Source modules
module_files <- list.files("R/modules", pattern = "\\.R$",
                           recursive = TRUE, full.names = TRUE)
for (f in module_files) {
    source(f)
}

# =============================================================================
# App Theme (Bootstrap 5 with bslib)
# =============================================================================

app_theme <- bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2c3e50",
    secondary = "#95a5a6",
    success = "#18bc9c",
    info = "#3498db",
    warning = "#f39c12",
    danger = "#e74c3c",
    base_font = font_google("Source Sans Pro"),
    heading_font = font_google("Source Sans Pro"),
    font_scale = 0.95
)

# =============================================================================
# DATA DIRECTORY DETECTION
# =============================================================================

# Priority 1: Environment variable set by run_viewer()
DATA_DIR <- Sys.getenv("MULTIOMICS_DATA_DIR", unset = "")

# Priority 2: Common relative locations (fallback)
if (DATA_DIR == "" || !dir.exists(DATA_DIR)) {
    candidates <- c(
        file.path(dirname(getwd()), "data"),  # ../data (viewer/data from shiny_app)
        file.path(getwd(), "data"),            # ./data
        "data"                                  # relative data
    )

    for (candidate in candidates) {
        if (dir.exists(candidate)) {
            DATA_DIR <- normalizePath(candidate)
            break
        }
    }

    # If still not found, default to ../data
    if (DATA_DIR == "" || !dir.exists(DATA_DIR)) {
        DATA_DIR <- file.path(dirname(getwd()), "data")
    }
}

message("========================================")
message("Multi-Omics Viewer - Data Discovery")
message("========================================")
message("Data directory: ", DATA_DIR)

# =============================================================================
# LAZY LOADING: Scan for available data directories and validate files
# =============================================================================

# Helper function to check for required files in a data type directory
validate_data_dir <- function(dir_path, data_type) {
    if (!dir.exists(dir_path)) {
        return(list(valid = FALSE, files = character(0), message = "Directory not found"))
    }

    # Look for CSV/RDS files
    csv_files <- list.files(dir_path, pattern = "\\.csv$", recursive = TRUE)
    rds_files <- list.files(dir_path, pattern = "\\.rds$", recursive = TRUE)
    all_files <- c(csv_files, rds_files)

    if (length(all_files) == 0) {
        return(list(valid = FALSE, files = character(0),
                    message = "No CSV or RDS files found"))
    }

    # Check for expected files based on data type
    expected_patterns <- switch(data_type,
        "rnaseq" = c("de_results", "pca", "qc", "normalized"),
        "proteomics" = c("differential", "normalized", "qc"),
        "metabolomics" = c("differential", "normalized", "pca"),
        "multiomics" = c("mofa", "diablo", "mae", "correlation")
    )

    found_expected <- sapply(expected_patterns, function(p) {
        any(grepl(p, all_files, ignore.case = TRUE))
    })

    list(
        valid = TRUE,
        files = all_files,
        expected_found = found_expected,
        message = paste(length(all_files), "files found")
    )
}

# Store paths and availability flags (NOT the actual data)
app_data <- list(
    # Paths to data directories
    data_dir = DATA_DIR,
    rnaseq_dir = NULL,
    proteomics_dir = NULL,
    metabolomics_dir = NULL,
    multiomics_dir = NULL,

    # Availability flags (for UI rendering)
    has_rnaseq = FALSE,
    has_proteomics = FALSE,
    has_metabolomics = FALSE,
    has_multiomics = FALSE,

    # File validation results
    rnaseq_files = character(0),
    proteomics_files = character(0),
    metabolomics_files = character(0),
    multiomics_files = character(0),

    # Project name
    project_name = "Multi-Omics Analysis"
)

# Check which data directories exist and have valid data
if (dir.exists(DATA_DIR)) {

    # RNA-seq
    rnaseq_path <- file.path(DATA_DIR, "rnaseq")
    validation <- validate_data_dir(rnaseq_path, "rnaseq")
    if (validation$valid) {
        app_data$rnaseq_dir <- rnaseq_path
        app_data$has_rnaseq <- TRUE
        app_data$rnaseq_files <- validation$files
        message("  [OK] RNA-seq: ", validation$message)
        if (any(!validation$expected_found)) {
            missing <- names(validation$expected_found)[!validation$expected_found]
            message("       Note: Some expected patterns not found: ",
                    paste(missing, collapse = ", "))
        }
    } else {
        message("  [--] RNA-seq: ", validation$message)
    }

    # Proteomics
    prot_path <- file.path(DATA_DIR, "proteomics")
    validation <- validate_data_dir(prot_path, "proteomics")
    if (validation$valid) {
        app_data$proteomics_dir <- prot_path
        app_data$has_proteomics <- TRUE
        app_data$proteomics_files <- validation$files
        message("  [OK] Proteomics: ", validation$message)
    } else {
        message("  [--] Proteomics: ", validation$message)
    }

    # Metabolomics
    metab_path <- file.path(DATA_DIR, "metabolomics")
    validation <- validate_data_dir(metab_path, "metabolomics")
    if (validation$valid) {
        app_data$metabolomics_dir <- metab_path
        app_data$has_metabolomics <- TRUE
        app_data$metabolomics_files <- validation$files
        message("  [OK] Metabolomics: ", validation$message)
    } else {
        message("  [--] Metabolomics: ", validation$message)
    }

    # Multi-omics
    multi_path <- file.path(DATA_DIR, "multiomics")
    validation <- validate_data_dir(multi_path, "multiomics")
    if (validation$valid) {
        app_data$multiomics_dir <- multi_path
        app_data$has_multiomics <- TRUE
        app_data$multiomics_files <- validation$files
        message("  [OK] Multi-omics: ", validation$message)
    } else {
        message("  [--] Multi-omics: ", validation$message)
    }

} else {
    message("  WARNING: Data directory does not exist!")
    message("  Path: ", DATA_DIR)
}

# Check if config file exists with project name
config_file <- file.path(DATA_DIR, "config.yml")
if (file.exists(config_file)) {
    config <- yaml::read_yaml(config_file)
    if (!is.null(config$project_name)) {
        app_data$project_name <- config$project_name
    }
}

# Summary
total_found <- sum(app_data$has_rnaseq, app_data$has_proteomics,
                   app_data$has_metabolomics, app_data$has_multiomics)

message("========================================")
if (total_found == 0) {
    message("WARNING: No valid data found!")
    message("Make sure your data directory contains subdirectories:")
    message("  rnaseq/       - with DE results, PCA data, etc.")
    message("  proteomics/   - with differential results, normalized matrix, etc.")
    message("  metabolomics/ - with normalized matrix, PCA data, etc.")
    message("  multiomics/   - with MOFA, DIABLO results, etc.")
} else {
    message("Found ", total_found, " data type(s). App ready to launch!")
}
message("========================================")

# =============================================================================
# Color Palettes
# =============================================================================

# For categorical variables
categorical_colors <- c(
    "#3498db", "#e74c3c", "#2ecc71", "#9b59b6", "#f39c12",
    "#1abc9c", "#e67e22", "#34495e", "#16a085", "#c0392b"
)

# For heatmaps
heatmap_colors <- colorRampPalette(c("#3498db", "white", "#e74c3c"))(100)

# For significance
significance_colors <- c(
    "up" = "#e74c3c",
    "down" = "#3498db",
    "ns" = "#95a5a6"
)
