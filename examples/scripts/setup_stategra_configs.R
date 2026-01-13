# =============================================================================
# Setup Pipeline Configurations for STATegra Mouse Data
# =============================================================================
# Updates all pipeline config files to point to the STATegra mouse B-cell
# differentiation data and sets appropriate analysis parameters.
#
# Usage: Rscript setup_stategra_configs.R
# =============================================================================

library(yaml)

# Base directories
examples_dir <- getwd()  # Should be run from examples/
multiomics_root <- dirname(examples_dir)
data_dir <- file.path(examples_dir, "data", "stategra")

cat("=== Setting Up Pipeline Configurations for STATegra (Mouse) ===\n\n")
cat("Examples directory:", examples_dir, "\n")
cat("Data directory:", data_dir, "\n\n")

# Verify data exists
if (!file.exists(file.path(data_dir, "rna_counts_matrix.csv"))) {
  stop("Data not found. Run download_stategra_data.R first.")
}

# -----------------------------------------------------------------------------
# Setup RNA-seq Pipeline Config
# -----------------------------------------------------------------------------
cat("Setting up rnaseq_pipeline config...\n")

rnaseq_config_path <- file.path(multiomics_root, "rnaseq_pipeline", "config.yml")
rnaseq_config <- read_yaml(rnaseq_config_path)

# Update paths (relative to rnaseq_pipeline/)
rnaseq_config$counts_file <- file.path("..", "examples", "data", "stategra", "rna_counts_matrix.csv")
rnaseq_config$metadata_file <- file.path("..", "examples", "data", "stategra", "metadata.csv")

# Update analysis settings for mouse
rnaseq_config$sample_id_column <- "sample_id"
rnaseq_config$condition_column <- "condition"
rnaseq_config$design_formula <- "~ condition"
rnaseq_config$contrasts <- list("differentiating - proliferating")
rnaseq_config$reference_level <- "proliferating"
rnaseq_config$organism <- "Mus musculus"

# QC parameters for cell line data
rnaseq_config$min_counts <- 10
rnaseq_config$min_samples <- 2  # Smaller sample size

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
prot_config$input$quant_matrix <- file.path("..", "examples", "data", "stategra", "prot_intensity_matrix.csv")
prot_config$input$metadata <- file.path("..", "examples", "data", "stategra", "metadata.csv")
prot_config$input$sample_id_column <- "sample_id"
prot_config$input$organism <- "Mus musculus"
prot_config$input$feature_id_type <- "uniprot"
prot_config$input$feature_level <- "protein"

# Update design
prot_config$design$condition_column <- "condition"
prot_config$design$design_formula <- "~ condition"
prot_config$design$contrasts <- list("differentiating - proliferating")
prot_config$design$reference_level <- "proliferating"

# Processing for cell line data
prot_config$processing$zeros_as_na <- TRUE
prot_config$processing$log_transform <- TRUE
prot_config$processing$normalization_method <- "median"

# Filtering - more relaxed for small sample size
prot_config$filtering$global_min_presence <- 0.5
prot_config$filtering$group_min_presence <- 0.6

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
metab_config$input$feature_matrix <- file.path("..", "examples", "data", "stategra", "metab_feature_matrix.csv")
metab_config$input$metadata <- file.path("..", "examples", "data", "stategra", "metadata.csv")
metab_config$input$feature_metadata <- file.path("..", "examples", "data", "stategra", "metab_feature_metadata.csv")
metab_config$input$sample_id_column <- "sample_id"

# Update design
metab_config$design$condition_column <- "condition"
metab_config$design$design_formula <- "~ condition"
metab_config$design$contrasts <- list("differentiating - proliferating")
metab_config$design$reference_level <- "proliferating"

# Processing
metab_config$processing$zeros_as_na <- TRUE
metab_config$processing$log_transform <- TRUE
metab_config$processing$normalization_method <- "PQN"

# QC
metab_config$qc$run_pca <- TRUE
metab_config$qc$n_pca_components <- 5  # Smaller for 6 samples

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
multi_config$global$organism <- "Mus musculus"
multi_config$global$metadata <- file.path("..", "examples", "data", "stategra", "metadata.csv")
multi_config$global$sample_id_column <- "sample_id"

# Design
multi_config$design$condition_column <- "condition"
multi_config$design$design_formula <- "~ condition"
multi_config$design$contrasts <- list("differentiating - proliferating")
multi_config$design$reference_level <- "proliferating"

# Input files - point to raw data
multi_config$transcriptomics$mode <- "raw"
multi_config$transcriptomics$raw$counts_matrix <- file.path("..", "examples", "data", "stategra", "rna_counts_matrix.csv")

multi_config$proteomics$mode <- "raw"
multi_config$proteomics$raw$intensity_matrix <- file.path("..", "examples", "data", "stategra", "prot_intensity_matrix.csv")

# Integration methods
multi_config$integration$methods <- list("MOFA2", "DIABLO")

# Feature selection - adjusted for smaller dataset
multi_config$feature_selection$transcriptomics$strategy <- "hybrid"
multi_config$feature_selection$transcriptomics$top_n_variance <- 1500
multi_config$feature_selection$proteomics$strategy <- "hybrid"
multi_config$feature_selection$proteomics$top_n_variance <- 800

write_yaml(multi_config, multi_config_path)
cat("  Updated:", multi_config_path, "\n")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=== Configuration Setup Complete ===\n\n")

cat("All pipeline configs have been updated with:\n")
cat("- Data paths pointing to examples/data/stategra/\n")
cat("- Condition: proliferating (0h) vs differentiating (24h)\n")
cat("- Organism: Mus musculus\n")
cat("- Commentary: enabled (data-driven mode)\n")

cat("\nDataset details:\n")
cat("- System: B3 pre-B-cell differentiation\n")
cat("- Samples: 6 (3 replicates x 2 conditions)\n")
cat("- Comparison: differentiating - proliferating\n")

cat("\nNext step: Run all pipelines with:\n")
cat("  Rscript scripts/run_all_pipelines.R\n")
