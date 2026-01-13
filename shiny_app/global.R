# =============================================================================
# Global Settings for Multi-Omics Viewer
# =============================================================================
# This file is sourced once when the app starts.
# It loads packages, sets up the theme, and loads the data.

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
                    "Run source('install.R') first."))
    }
    library(pkg, character.only = TRUE)
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
# Load Data from ../data/ directory
# =============================================================================

# Data directory is one level up from shiny_app
DATA_DIR <- file.path(dirname(getwd()), "data")
if (!dir.exists(DATA_DIR)) {
    DATA_DIR <- file.path(getwd(), "data")
}

# Load all pipeline data
message("Loading data from: ", DATA_DIR)

app_data <- list(
    rnaseq = NULL,
    proteomics = NULL,
    metabolomics = NULL,
    multiomics = NULL,
    project_name = "Multi-Omics Analysis"
)

# Try to load each pipeline's data
if (dir.exists(file.path(DATA_DIR, "rnaseq"))) {
    message("Loading RNA-seq data...")
    app_data$rnaseq <- load_rnaseq_data(file.path(DATA_DIR, "rnaseq"))
}

if (dir.exists(file.path(DATA_DIR, "proteomics"))) {
    message("Loading Proteomics data...")
    app_data$proteomics <- load_proteomics_data(file.path(DATA_DIR, "proteomics"))
}

if (dir.exists(file.path(DATA_DIR, "metabolomics"))) {
    message("Loading Metabolomics data...")
    app_data$metabolomics <- load_metabolomics_data(file.path(DATA_DIR, "metabolomics"))
}

if (dir.exists(file.path(DATA_DIR, "multiomics"))) {
    message("Loading Multi-omics data...")
    app_data$multiomics <- load_multiomics_data(file.path(DATA_DIR, "multiomics"))
}

# Check if config file exists with project name
config_file <- file.path(DATA_DIR, "config.yml")
if (file.exists(config_file)) {
    config <- yaml::read_yaml(config_file)
    if (!is.null(config$project_name)) {
        app_data$project_name <- config$project_name
    }
}

message("Data loading complete!")

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
