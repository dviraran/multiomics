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
# LAZY LOADING: Only scan for available data directories (don't load data)
# =============================================================================

# Data directory is one level up from shiny_app
DATA_DIR <- file.path(dirname(getwd()), "data")
if (!dir.exists(DATA_DIR)) {
    DATA_DIR <- file.path(getwd(), "data")
}

message("Scanning for data in: ", DATA_DIR)

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

    # Project name
    project_name = "Multi-Omics Analysis"
)

# Check which data directories exist
if (dir.exists(file.path(DATA_DIR, "rnaseq"))) {
    app_data$rnaseq_dir <- file.path(DATA_DIR, "rnaseq")
    app_data$has_rnaseq <- TRUE
    message("  Found RNA-seq data directory")
}

if (dir.exists(file.path(DATA_DIR, "proteomics"))) {
    app_data$proteomics_dir <- file.path(DATA_DIR, "proteomics")
    app_data$has_proteomics <- TRUE
    message("  Found Proteomics data directory")
}

if (dir.exists(file.path(DATA_DIR, "metabolomics"))) {
    app_data$metabolomics_dir <- file.path(DATA_DIR, "metabolomics")
    app_data$has_metabolomics <- TRUE
    message("  Found Metabolomics data directory")
}

if (dir.exists(file.path(DATA_DIR, "multiomics"))) {
    app_data$multiomics_dir <- file.path(DATA_DIR, "multiomics")
    app_data$has_multiomics <- TRUE
    message("  Found Multi-omics data directory")
}

# Check if config file exists with project name
config_file <- file.path(DATA_DIR, "config.yml")
if (file.exists(config_file)) {
    config <- yaml::read_yaml(config_file)
    if (!is.null(config$project_name)) {
        app_data$project_name <- config$project_name
    }
}

message("Data scan complete! App ready to launch.")

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
