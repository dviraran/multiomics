# =============================================================================
# Multi-Omics Viewer - Installation Script
# =============================================================================
# Run this script ONCE to install all required packages.
#
# Usage:
#   source("install_viewer.R")
#
# After installation, run:
#   source("run_viewer.R")
# =============================================================================

cat("
================================================================================
              Multi-Omics Viewer - Package Installation
================================================================================
")

# Function to install packages
install_if_missing <- function(packages, source = "CRAN") {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            cat(paste0("Installing ", pkg, "...\n"))

            if (source == "CRAN") {
                install.packages(pkg, quiet = TRUE)
            } else if (source == "Bioconductor") {
                if (!requireNamespace("BiocManager", quietly = TRUE)) {
                    install.packages("BiocManager", quiet = TRUE)
                }
                BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
            }

            if (requireNamespace(pkg, quietly = TRUE)) {
                cat(paste0("  -> ", pkg, " installed successfully\n"))
            } else {
                cat(paste0("  -> WARNING: ", pkg, " installation may have failed\n"))
            }
        } else {
            cat(paste0("  ", pkg, " - already installed\n"))
        }
    }
}

cat("\n[1/3] Installing core CRAN packages...\n")
cran_packages <- c(
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
    "ggrepel",
    "yaml"
)
install_if_missing(cran_packages, "CRAN")

cat("\n[2/3] Installing Bioconductor packages (optional, for enhanced features)...\n")
# These are optional - the viewer will work without them
bioc_packages <- c()
# Uncomment if you want ComplexHeatmap support:
# bioc_packages <- c("ComplexHeatmap")

if (length(bioc_packages) > 0) {
    install_if_missing(bioc_packages, "Bioconductor")
}

cat("\n[3/3] Verifying installation...\n")

# Test loading core packages
test_packages <- c("shiny", "bslib", "plotly", "DT", "heatmaply", "visNetwork",
                   "dplyr", "ggplot2")

all_ok <- TRUE
for (pkg in test_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(paste0("  ", pkg, " - OK\n"))
    } else {
        cat(paste0("  ", pkg, " - MISSING\n"))
        all_ok <- FALSE
    }
}

cat("\n================================================================================\n")
if (all_ok) {
    cat("
Installation complete!

To launch the Multi-Omics Viewer, run:
    source(\"run_viewer.R\")

The viewer will open in your default web browser.
================================================================================
")
} else {
    cat("
Some packages failed to install. Please check the error messages above
and try installing the missing packages manually.

Common solutions:
1. Make sure you have an internet connection
2. Try running R as administrator (Windows) or with sudo (Mac/Linux)
3. Update R to the latest version
================================================================================
")
}
