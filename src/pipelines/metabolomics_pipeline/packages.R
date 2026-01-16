# packages.R: Install and load all required packages for the pipeline
# Use renv for reproducibility. Pin MetaboAnalystR version for stability.

# List of required CRAN and GitHub packages
cran_packages <- c(
  "targets", "tarchetypes", "ggplot2", "plotly", "dplyr", "readr", "tibble", "quarto", "stringr", "purrr", "tidyr", "yaml", "logger",
  # Random forest and feature importance
  "ranger", "vip",
  # Testing
  "testthat"
)
github_packages <- c("xia-lab/MetaboAnalystR")

# Install CRAN packages if missing
installed <- rownames(installed.packages())
for (pkg in cran_packages) {
  if (!pkg %in% installed) install.packages(pkg)
}

# Install MetaboAnalystR from GitHub if missing (optional; controlled by env var INSTALL_METABOANALYST)
install_metabo <- Sys.getenv("INSTALL_METABOANALYST", "0") == "1"
if (install_metabo && !"MetaboAnalystR" %in% installed) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("xia-lab/MetaboAnalystR")
}

# Load all packages (conditionally load MetaboAnalystR)
core_pkgs <- cran_packages
if ("MetaboAnalystR" %in% installed && install_metabo) core_pkgs <- c(core_pkgs, "MetaboAnalystR")
invisible(lapply(core_pkgs, library, character.only = TRUE))
