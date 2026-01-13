#!/usr/bin/env Rscript
# =============================================================================
# Setup Pipeline Configurations for TCGA Breast Cancer Data
# =============================================================================
# Creates config files in examples/breast_tcga/configs/ for the breast.TCGA
# dataset from mixOmics. This dataset has REAL gene symbols that work with
# pathway analysis.
#
# Usage: Rscript setup_breast_tcga_configs.R
#        Then: Rscript run_all_pipelines.R --config-dir=../breast_tcga/configs
# =============================================================================

library(yaml)

# Find the multiomics root directory
find_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) > 0) {
    script_path <- sub("^--file=", "", script_arg[1])
    script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
    candidate <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
    if (file.exists(file.path(candidate, "pipeline_runner.R"))) {
      return(candidate)
    }
  }

  # Try common locations
  candidates <- c("..", "../..", ".", Sys.getenv("MULTIOMICS_ROOT"))
  for (candidate in candidates) {
    if (candidate == "") next
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (file.exists(file.path(candidate, "pipeline_runner.R"))) {
      return(candidate)
    }
  }
  stop("Cannot find multiomics root directory.")
}

multiomics_root <- find_root()
examples_dir <- file.path(multiomics_root, "examples")
data_dir <- file.path(examples_dir, "data", "breast_tcga")
config_dir <- file.path(examples_dir, "breast_tcga", "configs")

cat("=== Creating Pipeline Configurations for TCGA Breast Cancer ===\n\n")
cat("Multiomics root:", multiomics_root, "\n")
cat("Data directory:", data_dir, "\n")
cat("Config output:", config_dir, "\n\n")

# Verify data exists
if (!file.exists(file.path(data_dir, "mrna_counts_matrix.csv"))) {
  stop("Data not found. Run download_breast_tcga_data.R first.")
}

# Create config directory
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# RNA-seq Pipeline Config
# -----------------------------------------------------------------------------
cat("Creating rnaseq_config.yml...\n")

# Note: breast.TCGA mRNA data is already normalized, but we created pseudo-counts
# The gene IDs are REAL gene symbols (like BRCA1, TP53, ESR1, etc.)

rnaseq_config <- list(
  # Input files
  counts_file = file.path(data_dir, "mrna_counts_matrix.csv"),
  metadata_file = file.path(data_dir, "metadata.csv"),

  # Sample and condition settings
  sample_id_column = "sample_id",
  condition_column = "subtype",
  design_formula = "~ subtype",
  # Compare Basal (aggressive) vs LumA (better prognosis)
  contrasts = list("Basal - LumA"),
  reference_level = "LumA",

  # Organism settings - IMPORTANT: gene symbols, not Ensembl IDs
  organism = "Homo sapiens",
  gene_id_type = "symbol",  # Already gene symbols!
  strip_ensembl_version = FALSE,
  annotation_source = "auto",

  # QC parameters
  min_count = 10,
  min_samples = 3,
  alpha = 0.05,
  lfc_threshold = 0,
  outlier_sd_threshold = 3,

  # Pathway analysis - should work with real gene symbols

  pathway_database = list("GO", "KEGG"),
  pathway_method = "fgsea",
  pathway_min_size = 10,
  pathway_max_size = 500,

  # Commentary
  commentary_enabled = TRUE,
  commentary_backend = "none"
)

write_yaml(rnaseq_config, file.path(config_dir, "rnaseq_config.yml"))

# -----------------------------------------------------------------------------
# Proteomics Pipeline Config
# -----------------------------------------------------------------------------
cat("Creating proteomics_config.yml...\n")

# The RPPA protein IDs are based on gene symbols (with phospho modifications)
# We created a mapping file to help with this

proteomics_config <- list(
  input = list(
    quant_matrix = file.path(data_dir, "protein_intensity_matrix.csv"),
    metadata = file.path(data_dir, "metadata.csv"),
    sample_id_column = "sample_id",
    organism = "Homo sapiens",
    # RPPA uses gene symbol-based IDs
    feature_id_type = "gene_symbol",
    feature_level = "protein",
    # Mapping file to clean up phospho-protein names
    mapping_file = file.path(data_dir, "protein_id_mapping.csv")
  ),

  design = list(
    condition_column = "subtype",
    design_formula = "~ subtype",
    contrasts = list("Basal - LumA"),
    reference_level = "LumA"
  ),

  processing = list(
    zeros_as_na = TRUE,
    log_transform = FALSE,  # RPPA data is already normalized
    normalization_method = "none"  # Already normalized
  ),

  filtering = list(
    global_min_presence = 0.5,
    group_min_presence = 0.7,
    low_intensity_quantile = 0
  ),

  imputation = list(
    method = "knn",  # KNN works well for RPPA
    knn_k = 10
  ),

  qc = list(
    n_pca_components = 10,
    outlier_sd_threshold = 3,
    min_sample_correlation = 0.7,
    run_umap = TRUE,
    umap_n_neighbors = 15,
    umap_min_dist = 0.1
  ),

  differential = list(
    method = "limma",
    adj_pvalue_threshold = 0.05,
    log2fc_threshold = 0.5  # Lower threshold for RPPA
  ),

  pathway = list(
    run_pathway_analysis = TRUE,
    databases = list("GO", "KEGG"),
    go_ontologies = list("BP", "MF"),
    ora_pvalue_threshold = 0.05,
    min_set_size = 5,  # Smaller sets due to limited proteins
    max_set_size = 500
  ),

  output = list(
    output_dir = "outputs",
    report_format = "html"
  ),

  commentary = list(
    enabled = TRUE,
    backend = "none"
  )
)

write_yaml(proteomics_config, file.path(config_dir, "proteomics_config.yml"))

# -----------------------------------------------------------------------------
# Multi-omics Pipeline Config
# -----------------------------------------------------------------------------
cat("Creating multiomics_config.yml...\n")

# Integration of mRNA + protein (RPPA)
# Both have gene symbol-based IDs, so RNA-protein concordance should work

multiomics_config <- list(
  global = list(
    omics_present = list("transcriptomics", "proteomics"),
    organism = "Homo sapiens",
    metadata = file.path(data_dir, "metadata.csv"),
    sample_id_column = "sample_id"
  ),

  design = list(
    condition_column = "subtype",
    design_formula = "~ subtype",
    contrasts = list("Basal - LumA"),
    reference_level = "LumA"
  ),

  transcriptomics = list(
    mode = "preprocessed",  # Use normalized data directly
    preprocessed = list(
      normalized_matrix = file.path(data_dir, "mrna_normalized_matrix.csv")
    ),
    id_type = "symbol"
  ),

  proteomics = list(
    mode = "preprocessed",  # RPPA is already normalized
    preprocessed = list(
      normalized_matrix = file.path(data_dir, "protein_intensity_matrix.csv")
    ),
    id_type = "symbol"
  ),

  harmonization = list(
    sample_mode = "intersection",  # All 150 samples have all data
    duplicate_strategy = "keep_max_mean"
  ),

  feature_selection = list(
    transcriptomics = list(
      strategy = "all",  # Use all 200 genes
      top_n_variance = 200
    ),
    proteomics = list(
      strategy = "all",  # Use all 142 proteins
      top_n_variance = 150
    )
  ),

  integration = list(
    methods = list("MOFA2", "DIABLO"),
    mofa2 = list(
      num_factors = 5,  # Fewer factors for smaller dataset
      convergence_mode = "fast",
      seed = 42
    ),
    diablo = list(
      ncomp = 3,
      design_matrix = "full",
      cv_folds = 5,
      cv_repeats = 5,
      min_samples_per_group = 10
    )
  ),

  enrichment = list(
    run_enrichment = TRUE,
    methods = "ora",
    ora_pvalue = 0.05,
    min_set_size = 5,
    max_set_size = 500
  ),

  output = list(
    output_dir = "outputs",
    report_format = "html"
  ),

  commentary = list(
    enabled = TRUE,
    backend = "none"
  )
)

write_yaml(multiomics_config, file.path(config_dir, "multiomics_config.yml"))

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=== Configuration Setup Complete ===\n\n")

cat("Created config files in:", config_dir, "\n")
cat("  - rnaseq_config.yml\n")
cat("  - proteomics_config.yml\n")
cat("  - multiomics_config.yml\n")

cat("\nDataset details:\n")
cat("  - Source: TCGA Breast Cancer (from mixOmics package)\n")
cat("  - Samples: 150 breast cancer patients\n")
cat("  - Subtypes: Basal (45), Her2 (30), LumA (75)\n")
cat("  - Contrast: Basal vs LumA\n")
cat("  - Organism: Homo sapiens\n")

cat("\nKey features of this dataset:\n")
cat("  - REAL gene symbols (BRCA1, TP53, ESR1, etc.)\n")
cat("  - Works with GO, KEGG, Reactome pathway analysis\n")
cat("  - Matched samples across all omics\n")
cat("  - Biologically meaningful comparison (aggressive vs good prognosis)\n")

cat("\nNext step: Run pipelines with:\n")
cat("  cd ", examples_dir, "\n", sep = "")
cat("  Rscript scripts/run_all_pipelines.R --config-dir=breast_tcga/configs\n")
cat("\nOr from multiomics root:\n")
cat("  source('pipeline_runner.R')\n")
cat("  run_rnaseq_pipeline('examples/breast_tcga/configs/rnaseq_config.yml')\n")
