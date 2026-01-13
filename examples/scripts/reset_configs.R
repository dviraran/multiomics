#!/usr/bin/env Rscript
# =============================================================================
# Reset Pipeline Configurations to Defaults
# =============================================================================
# Resets all pipeline config.yml files to use their default example data paths.
# Useful when switching back from a custom dataset.
#
# Usage: Rscript reset_configs.R
# =============================================================================

library(yaml)

examples_dir <- getwd()
multiomics_root <- dirname(examples_dir)

cat("=== Resetting Pipeline Configurations to Defaults ===\n\n")

# -----------------------------------------------------------------------------
# Reset RNA-seq Pipeline
# -----------------------------------------------------------------------------
cat("Resetting rnaseq_pipeline config...\n")

rnaseq_config_path <- file.path(multiomics_root, "rnaseq_pipeline", "config.yml")
rnaseq_config <- read_yaml(rnaseq_config_path)

rnaseq_config$counts_file <- "data/counts_matrix.csv"
rnaseq_config$metadata_file <- "data/metadata.csv"
rnaseq_config$sample_id_column <- "sample_id"
rnaseq_config$condition_column <- "condition"
rnaseq_config$design_formula <- "~ condition"
rnaseq_config$contrasts <- list("treated - control")
rnaseq_config$reference_level <- "control"

write_yaml(rnaseq_config, rnaseq_config_path)
cat("  Reset:", rnaseq_config_path, "\n")

# -----------------------------------------------------------------------------
# Reset Proteomics Pipeline
# -----------------------------------------------------------------------------
cat("\nResetting proteomics_pipeline config...\n")

prot_config_path <- file.path(multiomics_root, "proteomics_pipeline", "config.yml")
prot_config <- read_yaml(prot_config_path)

prot_config$input$quant_matrix <- "data/example_quant_matrix.csv"
prot_config$input$metadata <- "data/example_metadata.csv"
prot_config$design$condition_column <- "condition"
prot_config$design$contrasts <- list("treated - control")
prot_config$design$reference_level <- "control"

write_yaml(prot_config, prot_config_path)
cat("  Reset:", prot_config_path, "\n")

# -----------------------------------------------------------------------------
# Reset Metabolomics Pipeline
# -----------------------------------------------------------------------------
cat("\nResetting metabolomics_pipeline config...\n")

metab_config_path <- file.path(multiomics_root, "metabolomics_pipeline", "config.yml")
metab_config <- read_yaml(metab_config_path)

metab_config$input$feature_matrix <- "data/example_feature_matrix.csv"
metab_config$input$metadata <- "data/example_metadata.csv"
metab_config$input$feature_metadata <- "data/example_feature_metadata.csv"
metab_config$design$condition_column <- "condition"
metab_config$design$contrasts <- list("treated - control")
metab_config$design$reference_level <- "control"

write_yaml(metab_config, metab_config_path)
cat("  Reset:", metab_config_path, "\n")

# -----------------------------------------------------------------------------
# Reset Multi-omics Pipeline
# -----------------------------------------------------------------------------
cat("\nResetting multiomics_pipeline config...\n")

multi_config_path <- file.path(multiomics_root, "multiomics_pipeline", "config.yml")
multi_config <- read_yaml(multi_config_path)

multi_config$global$metadata <- "data/metadata.csv"
multi_config$transcriptomics$raw$counts_matrix <- "data/rna_counts_matrix.csv"
multi_config$proteomics$raw$intensity_matrix <- "data/prot_intensity_matrix.csv"
multi_config$design$condition_column <- "condition"
multi_config$design$contrasts <- list("treated - control")
multi_config$design$reference_level <- "control"

write_yaml(multi_config, multi_config_path)
cat("  Reset:", multi_config_path, "\n")

cat("\n=== All configs reset to defaults ===\n")
