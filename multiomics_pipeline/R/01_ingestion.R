# =============================================================================
# Data Ingestion and Validation
# =============================================================================

#' Main ingestion function
ingest_all_data <- function(config) {
  log_message("=== Starting Multi-Omics Data Ingestion ===")
  create_output_dirs(config)

  # Load metadata
  metadata <- load_metadata(config)

  # Load each omics based on config
  omics_present <- config$global$omics_present
  omics_data <- list()

  if ("transcriptomics" %in% omics_present) {
    omics_data$transcriptomics <- load_transcriptomics(config)
  }
  if ("proteomics" %in% omics_present) {
    omics_data$proteomics <- load_proteomics(config)
  }
  if ("metabolomics" %in% omics_present) {
    omics_data$metabolomics <- load_metabolomics(config)
  }

  # Create sample alignment table
  alignment <- create_sample_alignment(metadata, omics_data, config)

  # Save alignment table
  save_table(alignment$alignment_table, "sample_alignment.csv", config)

  log_message("=== Data Ingestion Complete ===")

  list(
    metadata = metadata,
    omics_data = omics_data,
    alignment = alignment
  )
}

#' Load metadata
load_metadata <- function(config) {
  log_message("Loading metadata...")
  meta_path <- config$global$metadata

  if (!file.exists(meta_path)) {
    stop("Metadata file not found: ", meta_path)
  }

  metadata <- readr::read_csv(meta_path, show_col_types = FALSE)
  sample_col <- config$global$sample_id_column

  if (!sample_col %in% colnames(metadata)) {
    stop("Sample ID column '", sample_col, "' not found in metadata")
  }

  metadata[[sample_col]] <- sanitize_names(metadata[[sample_col]])

  # Check for duplicates
  if (any(duplicated(metadata[[sample_col]]))) {
    stop("Duplicate sample IDs found in metadata")
  }

  condition_col <- config$design$condition_column
  if (!condition_col %in% colnames(metadata)) {
    warning("Condition column '", condition_col, "' not found. Creating dummy.")
    metadata[[condition_col]] <- "sample"
  }

  metadata[[condition_col]] <- as.factor(metadata[[condition_col]])

  if (!is.null(config$design$reference_level)) {
    if (config$design$reference_level %in% levels(metadata[[condition_col]])) {
      metadata[[condition_col]] <- relevel(metadata[[condition_col]], ref = config$design$reference_level)
    }
  }

  log_message("Loaded metadata: ", nrow(metadata), " samples")
  return(metadata)
}

#' Load transcriptomics data
load_transcriptomics <- function(config) {
  log_message("Loading transcriptomics data...")
  tc <- config$transcriptomics

  if (tc$mode == "preprocessed") {
    return(load_preprocessed_omics(tc$preprocessed, "transcriptomics"))
  }

  # Raw mode
  counts_path <- tc$raw$counts_matrix
  if (!file.exists(counts_path)) {
    stop("RNA counts matrix not found: ", counts_path)
  }

  counts_df <- readr::read_csv(counts_path, show_col_types = FALSE)
  gene_ids <- as.character(counts_df[[1]])
  sample_cols <- colnames(counts_df)[-1]
  sample_cols <- sanitize_names(sample_cols)

  mat <- as.matrix(counts_df[, -1])
  rownames(mat) <- gene_ids
  colnames(mat) <- sample_cols
  storage.mode(mat) <- "numeric"

  # Load mapping file if provided
  mapping <- NULL
  if (!is.null(tc$raw$mapping_file) && file.exists(tc$raw$mapping_file)) {
    mapping <- readr::read_csv(tc$raw$mapping_file, show_col_types = FALSE)
    log_message("Loaded RNA mapping file")
  }

  # Load GMT if provided
  gmt <- NULL
  if (!is.null(tc$raw$gmt_file) && file.exists(tc$raw$gmt_file)) {
    gmt <- read_gmt(tc$raw$gmt_file)
    log_message("Loaded RNA GMT file: ", length(gmt), " gene sets")
  }

  log_message("Loaded RNA: ", nrow(mat), " genes x ", ncol(mat), " samples")

  list(
    matrix = mat,
    mode = "raw",
    mapping = mapping,
    gmt = gmt,
    de_table = NULL
  )
}

#' Load proteomics data
load_proteomics <- function(config) {
  log_message("Loading proteomics data...")
  pc <- config$proteomics

  if (pc$mode == "preprocessed") {
    return(load_preprocessed_omics(pc$preprocessed, "proteomics"))
  }

  # Raw mode
  intensity_path <- pc$raw$intensity_matrix
  if (!file.exists(intensity_path)) {
    stop("Proteomics intensity matrix not found: ", intensity_path)
  }

  int_df <- readr::read_csv(intensity_path, show_col_types = FALSE)
  protein_ids <- as.character(int_df[[1]])
  sample_cols <- sanitize_names(colnames(int_df)[-1])

  mat <- as.matrix(int_df[, -1])
  rownames(mat) <- protein_ids
  colnames(mat) <- sample_cols
  storage.mode(mat) <- "numeric"

  # Handle zeros
  if (pc$processing$zeros_as_na %||% TRUE) {
    n_zeros <- sum(mat == 0, na.rm = TRUE)
    mat[mat == 0] <- NA
    log_message("Converted ", n_zeros, " zeros to NA")
  }

  # Load mapping file
  mapping <- NULL
  if (!is.null(pc$raw$mapping_file) && file.exists(pc$raw$mapping_file)) {
    mapping <- readr::read_csv(pc$raw$mapping_file, show_col_types = FALSE)
  }

  log_message("Loaded Proteomics: ", nrow(mat), " proteins x ", ncol(mat), " samples")

  list(
    matrix = mat,
    mode = "raw",
    mapping = mapping,
    da_table = NULL
  )
}

#' Load metabolomics data
load_metabolomics <- function(config) {
  log_message("Loading metabolomics data...")
  mc <- config$metabolomics

  if (mc$mode == "preprocessed") {
    return(load_preprocessed_omics(mc$preprocessed, "metabolomics"))
  }

  # Raw mode
  feature_path <- mc$raw$feature_matrix
  if (!file.exists(feature_path)) {
    stop("Metabolomics feature matrix not found: ", feature_path)
  }

  feat_df <- readr::read_csv(feature_path, show_col_types = FALSE)
  feature_ids <- as.character(feat_df[[1]])
  sample_cols <- sanitize_names(colnames(feat_df)[-1])

  mat <- as.matrix(feat_df[, -1])
  rownames(mat) <- feature_ids
  colnames(mat) <- sample_cols
  storage.mode(mat) <- "numeric"

  if (mc$processing$zeros_as_na %||% TRUE) {
    mat[mat == 0] <- NA
  }

  # Load optional files
  feature_metadata <- NULL
  if (!is.null(mc$raw$feature_metadata) && file.exists(mc$raw$feature_metadata)) {
    feature_metadata <- readr::read_csv(mc$raw$feature_metadata, show_col_types = FALSE)
  }

  annotation <- NULL
  if (!is.null(mc$raw$annotation_table) && file.exists(mc$raw$annotation_table)) {
    annotation <- readr::read_csv(mc$raw$annotation_table, show_col_types = FALSE)
  }

  pathway_mapping <- NULL
  if (!is.null(mc$raw$pathway_mapping) && file.exists(mc$raw$pathway_mapping)) {
    pathway_mapping <- readr::read_csv(mc$raw$pathway_mapping, show_col_types = FALSE)
  }

  gmt <- NULL
  if (!is.null(mc$raw$gmt_file) && file.exists(mc$raw$gmt_file)) {
    gmt <- read_gmt(mc$raw$gmt_file)
  }

  log_message("Loaded Metabolomics: ", nrow(mat), " features x ", ncol(mat), " samples")

  list(
    matrix = mat,
    mode = "raw",
    feature_metadata = feature_metadata,
    annotation = annotation,
    pathway_mapping = pathway_mapping,
    gmt = gmt,
    da_table = NULL
  )
}

#' Load preprocessed omics data
load_preprocessed_omics <- function(preproc_config, omics_name) {
  log_message("Loading preprocessed ", omics_name, " data...")

  norm_path <- preproc_config$normalized_matrix
  if (is.null(norm_path) || !file.exists(norm_path)) {
    stop("Preprocessed ", omics_name, " normalized matrix not found")
  }

  norm_df <- readr::read_csv(norm_path, show_col_types = FALSE)
  feature_ids <- as.character(norm_df[[1]])
  sample_cols <- sanitize_names(colnames(norm_df)[-1])

  mat <- as.matrix(norm_df[, -1])
  rownames(mat) <- feature_ids
  colnames(mat) <- sample_cols
  storage.mode(mat) <- "numeric"

  # Load DE/DA table if provided
  de_table <- NULL
  de_path <- preproc_config$de_table %||% preproc_config$da_table
  if (!is.null(de_path) && file.exists(de_path)) {
    de_table <- readr::read_csv(de_path, show_col_types = FALSE)
    log_message("Loaded pre-computed DE/DA table")
  }

  log_message("Loaded preprocessed ", omics_name, ": ", nrow(mat), " features x ", ncol(mat), " samples")

  list(
    matrix = mat,
    mode = "preprocessed",
    de_table = de_table,
    da_table = de_table
  )
}

#' Create sample alignment table
create_sample_alignment <- function(metadata, omics_data, config) {
  sample_col <- config$global$sample_id_column
  meta_samples <- metadata[[sample_col]]

  alignment_df <- data.frame(sample_id = meta_samples, in_metadata = TRUE, stringsAsFactors = FALSE)

  all_omics_samples <- list()

  for (omics_name in names(omics_data)) {
    omics_samples <- colnames(omics_data[[omics_name]]$matrix)
    all_omics_samples[[omics_name]] <- omics_samples
    alignment_df[[paste0("in_", omics_name)]] <- meta_samples %in% omics_samples

    # Check for samples in omics but not metadata
    extra_samples <- setdiff(omics_samples, meta_samples)
    if (length(extra_samples) > 0) {
      warning(omics_name, " has ", length(extra_samples), " samples not in metadata: ",
              paste(head(extra_samples, 3), collapse = ", "))
    }
  }

  # Determine common samples
  sample_mode <- config$harmonization$sample_mode

  if (sample_mode == "intersection") {
    common_samples <- meta_samples
    for (omics_name in names(omics_data)) {
      common_samples <- intersect(common_samples, all_omics_samples[[omics_name]])
    }
  } else {
    common_samples <- meta_samples
    for (omics_name in names(omics_data)) {
      common_samples <- union(common_samples, all_omics_samples[[omics_name]])
    }
  }

  alignment_df$included <- alignment_df$sample_id %in% common_samples

  n_included <- sum(alignment_df$included)
  log_message("Sample alignment (", sample_mode, " mode): ", n_included, "/", length(meta_samples), " samples included")

  list(
    alignment_table = alignment_df,
    common_samples = common_samples,
    mode = sample_mode
  )
}
