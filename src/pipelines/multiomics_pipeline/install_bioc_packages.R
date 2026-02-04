if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c("clusterProfiler", "org.Ce.eg.db", "ggplot2")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

message("All packages installed successfully")
