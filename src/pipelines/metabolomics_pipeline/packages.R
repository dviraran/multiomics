# packages.R: Install and load all required packages for the pipeline
# Use renv for reproducibility. Pin MetaboAnalystR version for stability.
options(repos = c(CRAN = "https://cloud.r-project.org"))

# List of required CRAN and GitHub packages
cran_packages <- c(
  "targets", "tarchetypes", "ggplot2", "plotly", "dplyr", "readr", "tibble", "quarto", "stringr", "purrr", "tidyr", "yaml", "logger",
  "readxl", "ggrepel", "pheatmap", "knitr", "rmarkdown", "DT",
  # Random forest and feature importance
  "ranger", "vip",
  # Testing
  "testthat"
)
bioc_packages <- c("mixOmics")
github_packages <- c("xia-lab/MetaboAnalystR")

# Install CRAN packages if missing
installed <- rownames(installed.packages())
for (pkg in cran_packages) {
  if (!pkg %in% installed) install.packages(pkg)
}

# Install Bioconductor packages if missing (mixOmics for PLS-DA)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in bioc_packages) {
  if (!pkg %in% installed) BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

# Install MetaboAnalystR from GitHub if missing (optional; controlled by env var INSTALL_METABOANALYST)
install_metabo <- Sys.getenv("INSTALL_METABOANALYST", "1") == "1"
if (install_metabo && !"MetaboAnalystR" %in% installed) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("xia-lab/MetaboAnalystR")
}

# Load all packages (conditionally load MetaboAnalystR)
core_pkgs <- cran_packages
if ("MetaboAnalystR" %in% rownames(installed.packages()) || install_metabo) core_pkgs <- c(core_pkgs, "MetaboAnalystR")
invisible(lapply(core_pkgs, library, character.only = TRUE))
