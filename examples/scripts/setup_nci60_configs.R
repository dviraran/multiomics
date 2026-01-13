# =============================================================================
# Setup Pipeline Configurations for NCI-60 Data
# =============================================================================
# Updates all pipeline config files to point to the NCI-60 example data
# and sets appropriate analysis parameters.
#
# Usage: Rscript setup_nci60_configs.R
# =============================================================================

library(yaml)

# Base directories
examples_dir <- getwd()  # Should be run from examples/
multiomics_root <- dirname(examples_dir)
data_dir <- file.path(examples_dir, "data", "nci60")

cat("=== Setting Up Pipeline Configurations ===\n\n")
cat("Examples directory:", examples_dir, "\n")
cat("Data directory:", data_dir, "\n\n")

# Verify data exists
if (!file.exists(file.path(data_dir, "rna_counts_matrix.csv"))) {
  stop("Data not found. Run download_nci60_data.R first.")
}

# -----------------------------------------------------------------------------
# Setup RNA-seq Pipeline Config
# -----------------------------------------------------------------------------
cat("Setting up rnaseq_pipeline config...\n")

rnaseq_config_path <- file.path(multiomics_root, "rnaseq_pipeline", "config.yml")
rnaseq_config <- read_yaml(rnaseq_config_path)

# Update paths (relative to rnaseq_pipeline/)
rnaseq_config$counts_file <- file.path("..", "examples", "data", "nci60", "rna_counts_matrix.csv")
rnaseq_config$metadata_file <- file.path("..", "examples", "data", "nci60", "metadata.csv")

# Update analysis settings
rnaseq_config$sample_id_column <- "sample_id"
rnaseq_config$condition_column <- "condition"
rnaseq_config$design_formula <- "~ condition"
rnaseq_config$contrasts <- list("solid - hematological")
rnaseq_config$reference_level <- "hematological"
rnaseq_config$organism <- "Homo sapiens"

# QC parameters for cell lines
rnaseq_config$min_counts <- 10
rnaseq_config$min_samples <- 3

# Commentary settings
if (is.null(rnaseq_config$commentary_enabled)) {
  rnaseq_config$commentary_enabled <- TRUE
  rnaseq_config$commentary_backend <- "none"
}

write_yaml(rnaseq_config, rnaseq_config_path)
cat("  Updated:", rnaseq_config_path, "\n")

# -----------------------------------------------------------------------------
# Setup Proteomics Pipeline Config
# -----------------------------------------------------------------------------
cat("\nSetting up proteomics_pipeline config...\n")

prot_config_path <- file.path(multiomics_root, "proteomics_pipeline", "config.yml")
prot_config <- read_yaml(prot_config_path)

# Update input paths
prot_config$input$quant_matrix <- file.path("..", "examples", "data", "nci60", "prot_intensity_matrix.csv")
prot_config$input$metadata <- file.path("..", "examples", "data", "nci60", "metadata.csv")
prot_config$input$sample_id_column <- "sample_id"
prot_config$input$organism <- "Homo sapiens"
prot_config$input$feature_id_type <- "uniprot"
prot_config$input$feature_level <- "protein"

# Update design
prot_config$design$condition_column <- "condition"
prot_config$design$design_formula <- "~ condition"
prot_config$design$contrasts <- list("solid - hematological")
prot_config$design$reference_level <- "hematological"

# Processing for cell line data
prot_config$processing$zeros_as_na <- TRUE
prot_config$processing$log_transform <- TRUE
prot_config$processing$normalization_method <- "median"

# Filtering
prot_config$filtering$global_min_presence <- 0.5
prot_config$filtering$group_min_presence <- 0.7

# Imputation
prot_config$imputation$method <- "QRILC"

write_yaml(prot_config, prot_config_path)
cat("  Updated:", prot_config_path, "\n")

# -----------------------------------------------------------------------------
# Setup Metabolomics Pipeline Config
# -----------------------------------------------------------------------------
cat("\nSetting up metabolomics_pipeline config...\n")

metab_config_path <- file.path(multiomics_root, "metabolomics_pipeline", "config.yml")
metab_config <- read_yaml(metab_config_path)

# Update input paths
metab_config$input$feature_matrix <- file.path("..", "examples", "data", "nci60", "metab_feature_matrix.csv")
metab_config$input$metadata <- file.path("..", "examples", "data", "nci60", "metadata.csv")
metab_config$input$feature_metadata <- file.path("..", "examples", "data", "nci60", "metab_feature_metadata.csv")
metab_config$input$sample_id_column <- "sample_id"

# Update design
metab_config$design$condition_column <- "condition"
metab_config$design$design_formula <- "~ condition"
metab_config$design$contrasts <- list("solid - hematological")
metab_config$design$reference_level <- "hematological"

# Processing
metab_config$processing$zeros_as_na <- TRUE
metab_config$processing$log_transform <- TRUE
metab_config$processing$normalization_method <- "PQN"

# QC
metab_config$qc$run_pca <- TRUE
metab_config$qc$n_pca_components <- 10

write_yaml(metab_config, metab_config_path)
cat("  Updated:", metab_config_path, "\n")

# -----------------------------------------------------------------------------
# Setup Multi-omics Pipeline Config
# -----------------------------------------------------------------------------
cat("\nSetting up multiomics_pipeline config...\n")

multi_config_path <- file.path(multiomics_root, "multiomics_pipeline", "config.yml")
multi_config <- read_yaml(multi_config_path)

# Global settings
multi_config$global$omics_present <- list("transcriptomics", "proteomics")
multi_config$global$organism <- "Homo sapiens"
multi_config$global$metadata <- file.path("..", "examples", "data", "nci60", "metadata.csv")
multi_config$global$sample_id_column <- "sample_id"

# Design
multi_config$design$condition_column <- "condition"
multi_config$design$design_formula <- "~ condition"
multi_config$design$contrasts <- list("solid - hematological")
multi_config$design$reference_level <- "hematological"

# Input files - point to outputs from single-omics pipelines or raw data
multi_config$transcriptomics$mode <- "raw"
multi_config$transcriptomics$raw$counts_matrix <- file.path("..", "examples", "data", "nci60", "rna_counts_matrix.csv")

multi_config$proteomics$mode <- "raw"
multi_config$proteomics$raw$intensity_matrix <- file.path("..", "examples", "data", "nci60", "prot_intensity_matrix.csv")

# Integration methods
multi_config$integration$methods <- list("MOFA2", "DIABLO")

# Feature selection
multi_config$feature_selection$transcriptomics$strategy <- "hybrid"
multi_config$feature_selection$transcriptomics$top_n_variance <- 2000
multi_config$feature_selection$proteomics$strategy <- "hybrid"
multi_config$feature_selection$proteomics$top_n_variance <- 1000

write_yaml(multi_config, multi_config_path)
cat("  Updated:", multi_config_path, "\n")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=== Configuration Setup Complete ===\n\n")

cat("All pipeline configs have been updated with:\n")
cat("- Data paths pointing to examples/data/nci60/\n")
cat("- Condition: solid vs hematological tumors\n")
cat("- Organism: Homo sapiens\n")
cat("- Commentary: enabled (data-driven mode)\n")

cat("\nNext step: Run all pipelines with:\n")
cat("  Rscript scripts/run_all_pipelines.R\n")
