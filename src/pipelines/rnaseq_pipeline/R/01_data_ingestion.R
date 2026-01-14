# R/01_data_ingestion.R
# Functions for data ingestion and validation

#' Load parameters from config file
#' @param config_file Path to config YAML file
#' @return List of parameters
load_params <- function(config_file) {
  config <- yaml::read_yaml(config_file)

  # Set defaults for optional parameters
  defaults <- list(
    sample_id_col = "sample_id",
    group_col = "condition",
    design_formula = "~ condition",
    gene_id_type = "ensembl_gene_id",
    strip_ensembl_version = TRUE,
    annotation_source = "auto",
    mapping_file = NULL,
    gmt_file = NULL,
    pathway_database = c("GO", "KEGG"),
    pathway_method = "fgsea",
    pathway_min_size = 10,
    pathway_max_size = 500,
    min_count = 10,
    min_samples = NULL,
    alpha = 0.05,
    lfc_threshold = 0,
    outlier_sd_threshold = 3,
    contrasts = list()
  )

  # Merge with user config
  params <- modifyList(defaults, config)

  # Validate required parameters
  required <- c("counts_file", "metadata_file", "organism")
  missing <- setdiff(required, names(params))
  if (length(missing) > 0) {
    stop("Missing required parameters: ", paste(missing, collapse = ", "))
  }

  # Convert paths to absolute using here
  params$counts_file <- here::here(params$counts_file)
  params$metadata_file <- here::here(params$metadata_file)

  if (!is.null(params$mapping_file)) {
    params$mapping_file <- here::here(params$mapping_file)
  }
  if (!is.null(params$gmt_file)) {
    params$gmt_file <- here::here(params$gmt_file)
  }

  params
}


#' Read counts matrix from file
#' @param file_path Path to counts file (CSV or TSV)
#' @return Matrix of counts
read_counts_matrix <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Counts file not found: ", file_path)
  }

  # Detect delimiter
  ext <- tools::file_ext(file_path)
  if (ext == "tsv" || ext == "txt") {
    counts <- read.delim(file_path, row.names = 1, check.names = FALSE)
  } else {
    counts <- read.csv(file_path, row.names = 1, check.names = FALSE)
  }

  # Convert to matrix
  counts <- as.matrix(counts)

  # Validate counts are numeric
  if (!is.numeric(counts)) {
    # Find which columns are non-numeric to help debug
    non_numeric_cols <- which(!sapply(as.data.frame(counts), function(x) {
      all(grepl("^-?[0-9]*(\\.[0-9]+)?$", as.character(x)) | is.na(x))
    }))
    if (length(non_numeric_cols) > 0) {
      col_names <- colnames(counts)[non_numeric_cols[1:min(5, length(non_numeric_cols))]]
      stop("Counts matrix must contain numeric values.\n",
           "Non-numeric columns found: ", paste(col_names, collapse = ", "), "\n",
           "Check that your counts file has gene IDs as row names (first column) ",
           "and sample names as column headers.")
    }
    stop("Counts matrix must contain numeric values.\n",
         "Check that your counts file has gene IDs as row names (first column) ",
         "and sample names as column headers.")
  }

  # Round to integers if needed (with warning)
  if (any(counts != round(counts), na.rm = TRUE)) {
    warning("Counts contain non-integer values. Rounding to integers.")
    counts <- round(counts)
  }

  storage.mode(counts) <- "integer"

  message("Read counts matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")
  counts
}


#' Read metadata from file
#' @param file_path Path to metadata file
#' @param sample_id_col Expected sample ID column name
#' @return Data frame of metadata
read_metadata <- function(file_path, sample_id_col = "sample_id") {
  if (!file.exists(file_path)) {
    stop("Metadata file not found: ", file_path)
  }

  ext <- tools::file_ext(file_path)
  if (ext == "tsv" || ext == "txt") {
    metadata <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    metadata <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  }

  # Check if first column is unnamed (row names) and sample_id_col doesn't exist
  # This handles files where sample IDs are in an unnamed first column
  if (!sample_id_col %in% colnames(metadata)) {
    first_col <- colnames(metadata)[1]
    # Check if first column name is empty, "X", or looks like row names
    if (first_col == "" || first_col == "X" || grepl("^X\\.?[0-9]*$", first_col)) {
      message("Note: First column appears to be row names. Renaming to '", sample_id_col, "'")
      colnames(metadata)[1] <- sample_id_col
    }
  }

  message("Read metadata: ", nrow(metadata), " samples x ", ncol(metadata), " columns")
  metadata
}


#' Validate and harmonize counts and metadata
#' @param counts Counts matrix
#' @param metadata Metadata data frame
#' @param sample_id_col Column name for sample IDs in metadata
#' @return List with validated counts and metadata
validate_inputs <- function(counts, metadata, sample_id_col = "sample_id") {

  # Check sample ID column exists
  if (!sample_id_col %in% colnames(metadata)) {
    stop("Sample ID column '", sample_id_col, "' not found in metadata. ",
         "Available columns: ", paste(colnames(metadata), collapse = ", "))
  }

  # Get sample IDs
  meta_samples <- as.character(metadata[[sample_id_col]])
  count_samples <- colnames(counts)

  # Check for duplicates
  if (any(duplicated(meta_samples))) {
    stop("Duplicate sample IDs in metadata: ",
         paste(meta_samples[duplicated(meta_samples)], collapse = ", "))
  }
  if (any(duplicated(count_samples))) {
    stop("Duplicate sample names in counts matrix: ",
         paste(count_samples[duplicated(count_samples)], collapse = ", "))
  }

  # Find common samples
  common_samples <- intersect(count_samples, meta_samples)

  if (length(common_samples) == 0) {
    stop("No matching samples between counts and metadata.\n",
         "Counts samples: ", paste(head(count_samples), collapse = ", "), "\n",
         "Metadata samples: ", paste(head(meta_samples), collapse = ", "))
  }

  # Report mismatches
  only_in_counts <- setdiff(count_samples, meta_samples)
  only_in_metadata <- setdiff(meta_samples, count_samples)

  if (length(only_in_counts) > 0) {
    warning("Samples in counts but not metadata (will be dropped): ",
            paste(only_in_counts, collapse = ", "))
  }
  if (length(only_in_metadata) > 0) {
    warning("Samples in metadata but not counts (will be dropped): ",
            paste(only_in_metadata, collapse = ", "))
  }

  # Subset and reorder
  counts <- counts[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, meta_samples), , drop = FALSE]
  rownames(metadata) <- common_samples

  # Check for missing values
  if (any(is.na(counts))) {
    na_count <- sum(is.na(counts))
    warning("Found ", na_count, " NA values in counts matrix. Setting to 0.")
    counts[is.na(counts)] <- 0L
  }

  # Check for negative values
  if (any(counts < 0)) {
    stop("Counts matrix contains negative values")
  }

  message("Validated data: ", nrow(counts), " genes x ", ncol(counts), " samples")

  list(
    counts = counts,
    metadata = metadata
  )
}
