# =============================================================================
# Data Ingestion and Validation
# =============================================================================

#' Load and validate quantification matrix
#'
#' @param config Configuration list
#' @return List with matrix and validation results
load_quant_matrix <- function(config) {
  log_message("Loading quantification matrix...")

quant_path <- config$input$quant_matrix

  if (!file.exists(quant_path)) {
    stop("Quantification matrix not found: ", quant_path)
  }

  # Detect file format and load
  if (grepl("\\.tsv$", quant_path, ignore.case = TRUE)) {
    quant_df <- readr::read_tsv(quant_path, show_col_types = FALSE)
  } else {
    quant_df <- readr::read_csv(quant_path, show_col_types = FALSE)
  }

  # First column should be feature IDs
  feature_ids <- quant_df[[1]]
  sample_cols <- colnames(quant_df)[-1]

  # Convert to matrix
  quant_mat <- as.matrix(quant_df[, -1])
  rownames(quant_mat) <- feature_ids
  colnames(quant_mat) <- sample_cols

  # Ensure numeric
  storage.mode(quant_mat) <- "numeric"

  # Handle zeros as NA if configured
  if (config$processing$zeros_as_na) {
    n_zeros <- sum(quant_mat == 0, na.rm = TRUE)
    quant_mat[quant_mat == 0] <- NA
    log_message("Converted ", n_zeros, " zero values to NA")
  }

  # Validation checks
  validation <- validate_quant_matrix(quant_mat)

  log_message("Loaded matrix: ", nrow(quant_mat), " features x ", ncol(quant_mat), " samples")

  list(
    matrix = quant_mat,
    validation = validation
  )
}

#' Validate quantification matrix
#'
#' @param mat Numeric matrix
#' @return Validation results list
validate_quant_matrix <- function(mat) {
  results <- list(
    n_features = nrow(mat),
    n_samples = ncol(mat),
    n_missing = sum(is.na(mat)),
    pct_missing = mean(is.na(mat)) * 100,
    has_negative = any(mat < 0, na.rm = TRUE),
    intensity_range = range(mat, na.rm = TRUE),
    per_sample_missing = colSums(is.na(mat)),
    per_feature_missing = rowSums(is.na(mat)),
    duplicate_features = sum(duplicated(rownames(mat))),
    issues = character()
  )

  # Check for issues
  if (results$has_negative) {
    results$issues <- c(results$issues,
      "Warning: Matrix contains negative values. This is unusual for intensity data.")
  }

  if (results$pct_missing > 50) {
    results$issues <- c(results$issues,
      paste0("Warning: High missingness (", round(results$pct_missing, 1), "%). Consider relaxing filtering thresholds."))
  }

  if (results$duplicate_features > 0) {
    results$issues <- c(results$issues,
      paste0("Warning: ", results$duplicate_features, " duplicate feature IDs found. Will be handled according to config."))
  }

  # Check for samples with excessive missingness
  high_missing_samples <- names(which(results$per_sample_missing > nrow(mat) * 0.8))
  if (length(high_missing_samples) > 0) {
    results$issues <- c(results$issues,
      paste0("Warning: Samples with >80% missing: ", paste(high_missing_samples, collapse = ", ")))
  }

  return(results)
}

#' Load and validate metadata
#'
#' @param config Configuration list
#' @return Validated metadata data frame
load_metadata <- function(config) {
  log_message("Loading metadata...")

  meta_path <- config$input$metadata

  if (!file.exists(meta_path)) {
    stop("Metadata file not found: ", meta_path)
  }

  metadata <- readr::read_csv(meta_path, show_col_types = FALSE)

  # Validate required columns
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  if (!sample_col %in% colnames(metadata)) {
    stop("Sample ID column '", sample_col, "' not found in metadata. ",
         "Available columns: ", paste(colnames(metadata), collapse = ", "))
  }

  if (!condition_col %in% colnames(metadata)) {
    stop("Condition column '", condition_col, "' not found in metadata. ",
         "Available columns: ", paste(colnames(metadata), collapse = ", "))
  }

  # Check for batch column if specified
  batch_col <- config$design$batch_column
  if (!is.null(batch_col) && !batch_col %in% colnames(metadata)) {
    warning("Batch column '", batch_col, "' not found in metadata. Batch effects will not be modeled.")
    config$design$batch_column <- NULL
  }

  # Ensure condition is a factor
  metadata[[condition_col]] <- as.factor(metadata[[condition_col]])

  # Set reference level if specified
  if (!is.null(config$design$reference_level)) {
    if (config$design$reference_level %in% levels(metadata[[condition_col]])) {
      metadata[[condition_col]] <- relevel(metadata[[condition_col]],
                                           ref = config$design$reference_level)
    } else {
      warning("Reference level '", config$design$reference_level,
              "' not found in condition column. Using first level as reference.")
    }
  }

  log_message("Loaded metadata: ", nrow(metadata), " samples, ",
              length(unique(metadata[[condition_col]])), " conditions")

  return(metadata)
}

#' Match and align samples between matrix and metadata
#'
#' @param quant_data Quantification data list (from load_quant_matrix)
#' @param metadata Metadata data frame
#' @param config Configuration list
#' @return List with aligned matrix and metadata
align_samples <- function(quant_data, metadata, config) {
  log_message("Aligning samples between matrix and metadata...")

  mat <- quant_data$matrix
  sample_col <- config$input$sample_id_column

  matrix_samples <- colnames(mat)
  meta_samples <- metadata[[sample_col]]

  # Find intersection
  common_samples <- intersect(matrix_samples, meta_samples)

  if (length(common_samples) == 0) {
    stop("No matching samples between matrix and metadata!\n",
         "Matrix samples: ", paste(head(matrix_samples, 5), collapse = ", "), "...\n",
         "Metadata samples: ", paste(head(meta_samples, 5), collapse = ", "), "...")
  }

  # Report mismatches
  only_matrix <- setdiff(matrix_samples, meta_samples)
  only_meta <- setdiff(meta_samples, matrix_samples)

  if (length(only_matrix) > 0) {
    warning(length(only_matrix), " samples in matrix but not in metadata: ",
            paste(head(only_matrix, 5), collapse = ", "),
            if (length(only_matrix) > 5) "..." else "")
  }

  if (length(only_meta) > 0) {
    warning(length(only_meta), " samples in metadata but not in matrix: ",
            paste(head(only_meta, 5), collapse = ", "),
            if (length(only_meta) > 5) "..." else "")
  }

  # Subset and reorder
  mat <- mat[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, metadata[[sample_col]]), ]

  log_message("Aligned ", length(common_samples), " samples")

  list(
    matrix = mat,
    metadata = metadata,
    sample_matching = list(
      n_common = length(common_samples),
      only_in_matrix = only_matrix,
      only_in_metadata = only_meta
    )
  )
}

#' Handle duplicate features
#'
#' @param mat Numeric matrix with potential duplicate row names
#' @param config Configuration list
#' @return Matrix with unique feature IDs
handle_duplicates <- function(mat, config) {
  dup_ids <- rownames(mat)[duplicated(rownames(mat))]

  if (length(dup_ids) == 0) {
    log_message("No duplicate features found")
    return(mat)
  }

  log_message("Handling ", length(unique(dup_ids)), " duplicate feature IDs...")

  strategy <- config$duplicates$strategy %||% "keep_most_complete"

  unique_ids <- unique(rownames(mat))
  new_mat <- matrix(NA, nrow = length(unique_ids), ncol = ncol(mat))
  rownames(new_mat) <- unique_ids
  colnames(new_mat) <- colnames(mat)

  for (id in unique_ids) {
    idx <- which(rownames(mat) == id)

    if (length(idx) == 1) {
      new_mat[id, ] <- mat[idx, ]
    } else {
      sub_mat <- mat[idx, , drop = FALSE]

      if (strategy == "keep_most_complete") {
        # Keep row with fewest missing values
        n_missing <- rowSums(is.na(sub_mat))
        best_idx <- which.min(n_missing)
        new_mat[id, ] <- sub_mat[best_idx, ]
      } else if (strategy == "aggregate_median") {
        new_mat[id, ] <- apply(sub_mat, 2, median, na.rm = TRUE)
      } else if (strategy == "aggregate_mean") {
        new_mat[id, ] <- colMeans(sub_mat, na.rm = TRUE)
      } else if (strategy == "keep_first") {
        new_mat[id, ] <- sub_mat[1, ]
      }
    }
  }

  log_message("Reduced from ", nrow(mat), " to ", nrow(new_mat), " features using strategy: ", strategy)

  return(new_mat)
}

#' Create input summary for QC
#'
#' @param aligned_data Aligned data list
#' @param config Configuration list
#' @return Summary data frame
create_input_summary <- function(aligned_data, config) {
  mat <- aligned_data$matrix
  metadata <- aligned_data$metadata
  condition_col <- config$design$condition_column

  # Per-condition summary
  conditions <- unique(metadata[[condition_col]])

  condition_summary <- lapply(conditions, function(cond) {
    sample_col <- config$input$sample_id_column
    samples <- metadata[[sample_col]][metadata[[condition_col]] == cond]
    sub_mat <- mat[, samples, drop = FALSE]

    data.frame(
      condition = cond,
      n_samples = length(samples),
      mean_features_detected = mean(colSums(!is.na(sub_mat))),
      median_intensity = median(sub_mat, na.rm = TRUE),
      pct_missing = mean(is.na(sub_mat)) * 100,
      stringsAsFactors = FALSE
    )
  })

  condition_df <- do.call(rbind, condition_summary)

  # Overall summary
  overall <- data.frame(
    metric = c("Total features", "Total samples", "Overall % missing",
               "Intensity range (min)", "Intensity range (max)",
               "Features with 0% missing", "Features with >50% missing"),
    value = c(
      nrow(mat),
      ncol(mat),
      round(mean(is.na(mat)) * 100, 2),
      round(min(mat, na.rm = TRUE), 2),
      round(max(mat, na.rm = TRUE), 2),
      sum(rowSums(is.na(mat)) == 0),
      sum(rowSums(is.na(mat)) > ncol(mat) * 0.5)
    ),
    stringsAsFactors = FALSE
  )

  list(
    overall = overall,
    by_condition = condition_df
  )
}

#' Main data ingestion function
#'
#' @param config Configuration list
#' @return List with all ingested and validated data
ingest_data <- function(config) {
  log_message("=== Starting Data Ingestion ===")

  # Create output directories
  create_output_dirs(config)

  # Load data
  quant_data <- load_quant_matrix(config)
  metadata <- load_metadata(config)

  # Align samples
  aligned <- align_samples(quant_data, metadata, config)

  # Handle duplicates
  aligned$matrix <- handle_duplicates(aligned$matrix, config)

  # Create summary
  summary <- create_input_summary(aligned, config)

  # Save input summary
  save_table(summary$overall, "input_summary_overall.csv", config, "qc")
  save_table(summary$by_condition, "input_summary_by_condition.csv", config, "qc")

  log_message("=== Data Ingestion Complete ===")

  list(
    matrix = aligned$matrix,
    metadata = aligned$metadata,
    validation = quant_data$validation,
    sample_matching = aligned$sample_matching,
    input_summary = summary
  )
}
