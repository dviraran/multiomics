# List of required packages for the pipeline
packages <- c(
  "targets", "tarchetypes", "MetaboAnalystR", "ggplot2", "plotly", "dplyr", "readr", "tibble", "quarto", "stringr", "purrr", "tidyr", "renv"
)

# Install missing packages
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) install.packages(pkg)
}

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))
