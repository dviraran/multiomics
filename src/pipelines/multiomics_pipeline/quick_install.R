options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}

# Install only essential packages
BiocManager::install(c("clusterProfiler", "org.Ce.eg.db"),
                    update = FALSE, ask = FALSE, force = FALSE)

message("Done")
