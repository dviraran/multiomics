source("R/00_utils.R")
config <- load_config("config.yml")
print("Proteomics Config:")
print(config$proteomics$preprocessed)

path <- config$proteomics$preprocessed$da_table
print(paste("Path from config:", path))
print(paste("File exists?", file.exists(path)))
print(paste("Working directory:", getwd()))

# Check the logic in ingestion
de_path <- config$proteomics$preprocessed$de_table
if (is.null(de_path)) de_path <- config$proteomics$preprocessed$da_table
print(paste("Resolved de_path:", de_path))
