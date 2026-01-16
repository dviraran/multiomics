#!/usr/bin/env bash
# setup_kallisto_env.sh
# Activate conda env `kallisto_env` and install required R packages non-interactively.

set -euo pipefail

CONDA_BASE="$HOME/miniconda3"
if [ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
  echo "Conda not found at ${CONDA_BASE}. Adjust CONDA_BASE or install miniconda." >&2
  exit 1
fi

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate kallisto_env

echo "Installing CRAN packages in kallisto_env R..."
Rscript -e "options(repos = 'https://cloud.r-project.org'); pkgs <- c('targets','tarchetypes','callr','logger','yaml','tibble','dplyr','ggplot2','plotly','readr','tidyr','purrr','stringr'); missing <- pkgs[!pkgs %in% rownames(installed.packages())]; if(length(missing)) install.packages(missing);"

echo "Installing devtools and MetaboAnalystR from GitHub (if needed)..."
Rscript -e "if(!'devtools' %in% rownames(installed.packages())) install.packages('devtools'); devtools::install_github('xia-lab/MetaboAnalystR', upgrade = 'never')"

echo "Done. Validate pipeline with: Rscript -e 'targets::tar_validate("src/pipelines/metabolomics_pipeline/_targets.R")'"
