# Helper script to run the metabolomics pipeline with minimal shell quoting issues
if (!requireNamespace("targets", quietly = TRUE)) {
  install.packages("targets", repos = "https://cloud.r-project.org")
}
# Ensure CRAN repo set for reliability
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Run the pipeline
message("Starting targets::tar_make()...")
targets::tar_make()
message("targets::tar_make() finished")
