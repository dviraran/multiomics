# =============================================================================
# Multi-Omics Pipelines Docker Image
# =============================================================================
#
# This Dockerfile creates a reproducible environment for running all four
# multi-omics pipelines (RNA-seq, proteomics, metabolomics, multi-omics).
#
# Build:
#   docker build -t multiomics-pipelines .
#
# Run interactively:
#   docker run -it -v $(pwd)/data:/data -v $(pwd)/outputs:/outputs multiomics-pipelines
#
# Run a specific pipeline:
#   docker run -v $(pwd):/work multiomics-pipelines Rscript -e "source('pipeline_runner.R'); run_rnaseq_pipeline('/work/config.yml')"

FROM rocker/r-ver:4.3.2

LABEL maintainer="Multi-Omics Pipelines"
LABEL description="Reproducible environment for multi-omics analysis pipelines"
LABEL version="1.0"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libglpk-dev \
    libgmp-dev \
    libhdf5-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Set CRAN mirror
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >> /usr/local/lib/R/etc/Rprofile.site

# Install BiocManager first
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.18', ask = FALSE)"

# Install core CRAN packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'yaml', \
    'targets', \
    'tarchetypes', \
    'future', \
    'future.callr', \
    'R.utils', \
    'rmarkdown', \
    'knitr', \
    'DT', \
    'plotly', \
    'pheatmap', \
    'RColorBrewer', \
    'ggrepel', \
    'scales', \
    'gridExtra', \
    'patchwork', \
    'circlize', \
    'igraph', \
    'optparse' \
))"

# Install Bioconductor packages for RNA-seq
RUN R -e "BiocManager::install(c( \
    'DESeq2', \
    'edgeR', \
    'limma', \
    'fgsea', \
    'clusterProfiler', \
    'org.Hs.eg.db', \
    'org.Mm.eg.db', \
    'AnnotationDbi', \
    'biomaRt', \
    'GSVA', \
    'sva', \
    'vsn' \
), ask = FALSE, update = FALSE)"

# Install Bioconductor packages for multi-omics
RUN R -e "BiocManager::install(c( \
    'MultiAssayExperiment', \
    'MOFA2', \
    'mixOmics', \
    'SNFtool', \
    'ComplexHeatmap' \
), ask = FALSE, update = FALSE)"

# Install additional packages for proteomics/metabolomics
RUN R -e "BiocManager::install(c( \
    'MSnbase', \
    'UniProt.ws', \
    'impute', \
    'pcaMethods' \
), ask = FALSE, update = FALSE)"

# Install optional packages (non-critical, allow failures)
RUN R -e "tryCatch(BiocManager::install('WGCNA', ask = FALSE), error = function(e) message('WGCNA not installed'))"
RUN R -e "tryCatch(install.packages('doParallel'), error = function(e) message('doParallel not installed'))"

# Create app directory
WORKDIR /app

# Copy pipeline code
COPY . /app/

# Create directories for user data
RUN mkdir -p /data /outputs

# Set environment variables
ENV PIPELINE_CONFIG=/data/config.yml
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# Default command: start R
CMD ["R"]
