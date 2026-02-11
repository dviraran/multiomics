# Install MetaboAnalystR and its dependencies
# This script should be run once to set up the environment

message("Installing MetaboAnalystR dependencies...")

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor dependencies
bioc_packages <- c(
  "impute", "pcaMethods", "globaltest", "GlobalAncova", 
  "Rgraphviz", "preprocessCore", "genefilter", "SSPA", 
  "sva", "limma", "KEGGgraph", "siggenes", "BiocParallel",
  "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# Install CRAN dependencies
cran_packages <- c(
  "Rserve", "RColorBrewer", "xtable", "som", "ROCR", 
  "RJSONIO", "gplots", "e1071", "caTools", "igraph", 
  "randomForest", "Cairo", "pls", "pheatmap", "lattice",
  "rmarkdown", "knitr", "data.table", "pROC", "caret",
  "ellipse", "scatterplot3d", "reshape2", "scales"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Install MetaboAnalystR from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

message("Installing MetaboAnalystR from GitHub...")
devtools::install_github("xia-lab/MetaboAnalystR", 
                         dependencies = FALSE,
                         upgrade = "never")

# Verify installation
if (requireNamespace("MetaboAnalystR", quietly = TRUE)) {
  message("\n✓ MetaboAnalystR successfully installed!")
  message("Version: ", packageVersion("MetaboAnalystR"))
} else {
  message("\n✗ MetaboAnalystR installation failed")
}
