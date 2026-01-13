#!/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript
# =============================================================================
# Setup Pipeline Configurations for STATegra Mouse Data
# =============================================================================
# Creates config files in examples/stategra/configs/ for the STATegra mouse
# B-cell differentiation data. Does NOT modify pipeline default configs.
#
# Usage: Rscript setup_stategra_configs.R
#        Then: Rscript run_all_pipelines.R --config-dir=../stategra/configs
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
data_dir <- file.path(examples_dir, "data", "stategra")
config_dir <- file.path(examples_dir, "stategra", "configs")

cat("=== Creating Pipeline Configurations for STATegra (Mouse) ===\n\n")
cat("Multiomics root:", multiomics_root, "\n")
cat("Data directory:", data_dir, "\n")
cat("Config output:", config_dir, "\n\n")

# Verify data exists
if (!file.exists(file.path(data_dir, "rna_counts_matrix.csv"))) {
  stop("Data not found. Run download_stategra_data.R first.")
}

# Create config directory
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# RNA-seq Pipeline Config
# -----------------------------------------------------------------------------
cat("Creating rnaseq_config.yml...\n")

rnaseq_config <- list(
  # Input files - absolute paths work best
  counts_file = file.path(data_dir, "rna_counts_matrix.csv"),
  metadata_file = file.path(data_dir, "metadata.csv"),

  # Sample and condition settings
  sample_id_column = "sample_id",
  condition_column = "condition",
  design_formula = "~ condition",
  contrasts = list("differentiating - proliferating"),
  reference_level = "proliferating",

  # Organism settings
  organism = "Mus musculus",
  gene_id_type = "ensembl_gene_id",
  strip_ensembl_version = TRUE,
  annotation_source = "auto",

  # QC parameters for small dataset
  min_count = 10,
  min_samples = 2,
  alpha = 0.05,
  lfc_threshold = 0,
  outlier_sd_threshold = 3,

  # Pathway analysis
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

proteomics_config <- list(
  input = list(
    quant_matrix = file.path(data_dir, "prot_intensity_matrix.csv"),
    metadata = file.path(data_dir, "metadata.csv"),
    sample_id_column = "sample_id",
    organism = "Mus musculus",
    feature_id_type = "uniprot",
    feature_level = "protein"
  ),

  design = list(
    condition_column = "condition",
    design_formula = "~ condition",
    contrasts = list("differentiating - proliferating"),
    reference_level = "proliferating"
  ),

  processing = list(
    zeros_as_na = TRUE,
    log_transform = TRUE,
    normalization_method = "median"
  ),

  filtering = list(
    global_min_presence = 0.5,
    group_min_presence = 0.6,  # Relaxed for small sample size
    low_intensity_quantile = 0
  ),

  imputation = list(
    method = "QRILC",
    min_prob_quantile = 0.01,
    tune_sigma = 1
  ),

  qc = list(
    n_pca_components = 5,  # Smaller for 6 samples
    outlier_sd_threshold = 3,
    min_sample_correlation = 0.8,
    run_umap = FALSE  # Skip for small dataset
  ),

  differential = list(
    method = "limma",
    adj_pvalue_threshold = 0.05,
    log2fc_threshold = 1
  ),

  pathway = list(
    run_pathway_analysis = TRUE,
    databases = list("GO", "KEGG"),
    go_ontologies = list("BP", "MF"),
    ora_pvalue_threshold = 0.05,
    min_set_size = 10,
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
# Metabolomics Pipeline Config
# -----------------------------------------------------------------------------
cat("Creating metabolomics_config.yml...\n")

metabolomics_config <- list(
  input = list(
    feature_matrix = file.path(data_dir, "metab_feature_matrix.csv"),
    metadata = file.path(data_dir, "metadata.csv"),
    feature_metadata = file.path(data_dir, "metab_feature_metadata.csv"),
    sample_id_column = "sample_id"
  ),

  design = list(
    condition_column = "condition",
    design_formula = "~ condition",
    contrasts = list("differentiating - proliferating"),
    reference_level = "proliferating"
  ),

  processing = list(
    zeros_as_na = TRUE,
    log_transform = TRUE,
    normalization_method = "PQN"
  ),

  qc = list(
    run_pca = TRUE,
    n_pca_components = 5  # Smaller for 6 samples
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

write_yaml(metabolomics_config, file.path(config_dir, "metabolomics_config.yml"))

# -----------------------------------------------------------------------------
# Multi-omics Pipeline Config
# -----------------------------------------------------------------------------
cat("Creating multiomics_config.yml...\n")

multiomics_config <- list(
  global = list(
    omics_present = list("transcriptomics", "proteomics"),
    organism = "Mus musculus",
    metadata = file.path(data_dir, "metadata.csv"),
    sample_id_column = "sample_id"
  ),

  design = list(
    condition_column = "condition",
    design_formula = "~ condition",
    contrasts = list("differentiating - proliferating"),
    reference_level = "proliferating"
  ),

  transcriptomics = list(
    mode = "raw",
    raw = list(
      counts_matrix = file.path(data_dir, "rna_counts_matrix.csv")
    ),
    processing = list(
      min_counts = 10,
      min_samples = 2,
      normalization = "vst",
      de_method = "DESeq2"
    )
  ),

  proteomics = list(
    mode = "raw",
    raw = list(
      intensity_matrix = file.path(data_dir, "prot_intensity_matrix.csv")
    ),
    processing = list(
      zeros_as_na = TRUE,
      log_transform = TRUE,
      normalization = "median",
      min_presence = 0.5,
      imputation = "half_min"
    )
  ),

  harmonization = list(
    sample_mode = "intersection",
    duplicate_strategy = "keep_max_mean"
  ),

  feature_selection = list(
    transcriptomics = list(
      strategy = "hybrid",
      top_n_variance = 1500,
      fdr_threshold = 0.05
    ),
    proteomics = list(
      strategy = "hybrid",
      top_n_variance = 800,
      fdr_threshold = 0.05
    )
  ),

  integration = list(
    methods = list("MOFA2", "DIABLO"),
    mofa2 = list(
      num_factors = 5,  # Smaller for 6 samples
      convergence_mode = "fast",
      seed = 42
    ),
    diablo = list(
      ncomp = 2,  # Smaller for 6 samples
      design_matrix = "full",
      cv_folds = 3,  # Reduced for small sample
      cv_repeats = 5,
      min_samples_per_group = 2
    )
  ),

  enrichment = list(
    run_enrichment = TRUE,
    methods = "ora",
    ora_pvalue = 0.05,
    min_set_size = 10,
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
cat("  - metabolomics_config.yml\n")
cat("  - multiomics_config.yml\n")

cat("\nDataset details:\n")
cat("  - System: B3 pre-B-cell differentiation\n")
cat("  - Samples: 6 (3 replicates x 2 conditions)\n")
cat("  - Organism: Mus musculus\n")
cat("  - Comparison: differentiating - proliferating\n")

cat("\nNext step: Run all pipelines with:\n")
cat("  cd examples\n")
cat("  Rscript scripts/run_all_pipelines.R --config-dir=stategra/configs\n")
