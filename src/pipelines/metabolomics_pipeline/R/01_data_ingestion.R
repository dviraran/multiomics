# =============================================================================
# Data Ingestion and Validation
# =============================================================================

#' Load and validate all input data
#'
#' @param config Configuration list
#' @return List with all ingested data
ingest_data <- function(config) {
  log_message("=== Starting Data Ingestion ===")

  # Create output directories
  create_output_dirs(config)

  # Load feature matrix
  feature_data <- load_feature_matrix(config)

  # Load metadata
  metadata <- load_metadata(config)

  # Align samples
  aligned <- align_samples(feature_data, metadata, config)

  # Handle duplicates
  aligned$matrix <- handle_duplicates(aligned$matrix, config)

  # Identify sample roles (QC, Blank, Sample)
  sample_roles <- identify_sample_roles(aligned$metadata, config)

  # Load optional feature metadata
  feature_metadata <- load_feature_metadata(aligned$matrix, config)

  # Load optional annotations
  annotations <- load_annotations(config)

  # Create input summary
  summary <- create_input_summary(aligned, sample_roles, config)

  # Save summaries
  save_table(summary$overall, "input_summary_overall.csv", config, "qc")
  save_table(summary$by_condition, "input_summary_by_condition.csv", config, "qc")
  if (!is.null(sample_roles$summary)) {
    save_table(sample_roles$summary, "sample_types_summary.csv", config, "qc")
  }

  log_message("=== Data Ingestion Complete ===")

  list(
    matrix = aligned$matrix,
    metadata = aligned$metadata,
    sample_roles = sample_roles,
    feature_metadata = feature_metadata,
    annotations = annotations,
    validation = feature_data$validation,
    sample_matching = aligned$sample_matching,
    input_summary = summary
  )
}

#' Load and validate feature matrix
#'
#' @param config Configuration list
#' @return List with matrix and validation results
load_feature_matrix <- function(config) {
  log_message("Loading feature matrix...")

  matrix_path <- config$input$feature_matrix

  if (!file.exists(matrix_path)) {
    stop("Feature matrix not found: ", matrix_path)
  }

  # Detect file format and load
  if (grepl("\\.tsv$", matrix_path, ignore.case = TRUE)) {
    mat_df <- readr::read_tsv(matrix_path, show_col_types = FALSE)
  } else {
    mat_df <- readr::read_csv(matrix_path, show_col_types = FALSE)
  }

  # First column should be feature IDs
  feature_ids <- as.character(mat_df[[1]])
  sample_cols <- colnames(mat_df)[-1]

  # Sanitize sample names
  sample_cols <- sanitize_names(sample_cols)

  # Convert to matrix
  mat <- as.matrix(mat_df[, -1])
  rownames(mat) <- feature_ids
  colnames(mat) <- sample_cols

  # Ensure numeric
  storage.mode(mat) <- "numeric"

  # Handle zeros as NA if configured
  if (config$processing$zeros_as_na) {
    n_zeros <- sum(mat == 0, na.rm = TRUE)
    mat[mat == 0] <- NA
    log_message("Converted ", n_zeros, " zero values to NA")
  }

  # Validation checks
  validation <- validate_feature_matrix(mat)

  log_message("Loaded matrix: ", nrow(mat), " features x ", ncol(mat), " samples")

  list(
    matrix = mat,
    validation = validation
  )
}

#' Validate feature matrix
#'
#' @param mat Numeric matrix
#' @return Validation results list
validate_feature_matrix <- function(mat) {
  results <- list(
    n_features = nrow(mat),
    n_samples = ncol(mat),
    n_missing = sum(is.na(mat)),
    pct_missing = mean(is.na(mat)) * 100,
    has_negative = any(mat < 0, na.rm = TRUE),
    has_infinite = any(!is.finite(mat) & !is.na(mat)),
    intensity_range = range(mat, na.rm = TRUE),
    per_sample_missing = colSums(is.na(mat)),
    per_feature_missing = rowSums(is.na(mat)),
    duplicate_features = sum(duplicated(rownames(mat))),
    issues = character()
  )

  # Check for issues
  if (results$has_negative) {
    results$issues <- c(results$issues,
      "Warning: Matrix contains negative values. Check if data is already log-transformed.")
  }

  if (results$has_infinite) {
    results$issues <- c(results$issues,
      "Warning: Matrix contains infinite values. These will be treated as NA.")
    # Replace Inf with NA
  }

  if (results$pct_missing > 70) {
    results$issues <- c(results$issues,
      paste0("Warning: Very high missingness (", round(results$pct_missing, 1),
             "%). Consider relaxing filtering thresholds."))
  }

  if (results$duplicate_features > 0) {
    results$issues <- c(results$issues,
      paste0("Warning: ", results$duplicate_features, " duplicate feature IDs found."))
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

  # Sanitize sample names in metadata
  metadata[[sample_col]] <- sanitize_names(metadata[[sample_col]])

  # Check for condition column (may not exist if only QC samples)
  if (!condition_col %in% colnames(metadata)) {
    warning("Condition column '", condition_col, "' not found in metadata. ",
            "Creating dummy condition for all samples.")
    metadata[[condition_col]] <- "sample"
  }

  # Ensure condition is a factor
  metadata[[condition_col]] <- as.factor(metadata[[condition_col]])

  # Set reference level if specified
  if (!is.null(config$design$reference_level)) {
    if (config$design$reference_level %in% levels(metadata[[condition_col]])) {
      metadata[[condition_col]] <- relevel(metadata[[condition_col]],
                                           ref = config$design$reference_level)
    }
  }

  log_message("Loaded metadata: ", nrow(metadata), " samples")

  return(metadata)
}

#' Match and align samples between matrix and metadata
#'
#' @param feature_data Feature data list
#' @param metadata Metadata data frame
#' @param config Configuration list
#' @return List with aligned matrix and metadata
align_samples <- function(feature_data, metadata, config) {
  log_message("Aligning samples between matrix and metadata...")

  mat <- feature_data$matrix
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
            paste(head(only_matrix, 5), collapse = ", "))
  }

  if (length(only_meta) > 0) {
    warning(length(only_meta), " samples in metadata but not in matrix: ",
            paste(head(only_meta, 5), collapse = ", "))
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

#' Handle duplicate feature IDs
#'
#' @param mat Numeric matrix
#' @param config Configuration list
#' @return Matrix with unique feature IDs
handle_duplicates <- function(mat, config) {
  dup_ids <- rownames(mat)[duplicated(rownames(mat))]

  if (length(dup_ids) == 0) {
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
        n_missing <- rowSums(is.na(sub_mat))
        best_idx <- which.min(n_missing)
        new_mat[id, ] <- sub_mat[best_idx, ]
      } else if (strategy == "aggregate_median") {
        new_mat[id, ] <- apply(sub_mat, 2, median, na.rm = TRUE)
      } else if (strategy == "aggregate_mean") {
        new_mat[id, ] <- colMeans(sub_mat, na.rm = TRUE)
      } else {
        new_mat[id, ] <- sub_mat[1, ]
      }
    }
  }

  log_message("Reduced from ", nrow(mat), " to ", nrow(new_mat), " features")

  return(new_mat)
}

#' Identify sample roles (QC, Blank, Sample)
#'
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Sample roles list
identify_sample_roles <- function(metadata, config) {
  sample_type_col <- config$sample_types$sample_type_column
  sample_col <- config$input$sample_id_column

  # Check if sample type column exists
  if (!sample_type_col %in% colnames(metadata)) {
    log_message("No sample type column found. Treating all samples as biological samples.")
    return(list(
      has_sample_types = FALSE,
      samples = metadata[[sample_col]],
      qc_samples = character(),
      blank_samples = character(),
      biological_samples = metadata[[sample_col]],
      summary = NULL
    ))
  }

  sample_types <- metadata[[sample_type_col]]
  sample_ids <- metadata[[sample_col]]

  sample_val <- config$sample_types$sample_value
  qc_val <- config$sample_types$qc_value
  blank_val <- config$sample_types$blank_value

  biological <- sample_ids[sample_types == sample_val]
  qc <- sample_ids[sample_types == qc_val]
  blanks <- sample_ids[sample_types == blank_val]

  summary_df <- data.frame(
    sample_type = c("Biological samples", "QC samples", "Blank samples", "Other"),
    count = c(
      length(biological),
      length(qc),
      length(blanks),
      nrow(metadata) - length(biological) - length(qc) - length(blanks)
    ),
    stringsAsFactors = FALSE
  )

  log_message("Sample types: ", length(biological), " biological, ",
              length(qc), " QC, ", length(blanks), " blank")

  list(
    has_sample_types = TRUE,
    samples = sample_ids,
    qc_samples = qc,
    blank_samples = blanks,
    biological_samples = biological,
    summary = summary_df
  )
}

#' Load optional feature metadata
#'
#' @param mat Feature matrix
#' @param config Configuration list
#' @return Feature metadata data frame or NULL
load_feature_metadata <- function(mat, config) {
  if (is.null(config$optional_inputs$feature_metadata)) {
    return(NULL)
  }

  meta_path <- config$optional_inputs$feature_metadata
  if (!file.exists(meta_path)) {
    log_message("Feature metadata file not found: ", meta_path)
    return(NULL)
  }

  log_message("Loading feature metadata...")
  feature_meta <- readr::read_csv(meta_path, show_col_types = FALSE)

  # Check alignment with matrix
  if ("feature_id" %in% colnames(feature_meta)) {
    matched <- sum(rownames(mat) %in% feature_meta$feature_id)
    log_message("Feature metadata matches ", matched, "/", nrow(mat), " features")
  }

  return(feature_meta)
}

#' Load optional annotations
#'
#' @param config Configuration list
#' @return Annotations list or NULL
load_annotations <- function(config) {
  annotations <- list()

  # Load annotation table
  if (!is.null(config$optional_inputs$annotation_table)) {
    if (file.exists(config$optional_inputs$annotation_table)) {
      log_message("Loading annotation table...")
      annotations$annotation_table <- readr::read_csv(
        config$optional_inputs$annotation_table,
        show_col_types = FALSE
      )
    }
  }

  # Load pathway mapping
  if (!is.null(config$optional_inputs$pathway_mapping)) {
    if (file.exists(config$optional_inputs$pathway_mapping)) {
      log_message("Loading pathway mapping...")
      annotations$pathway_mapping <- readr::read_csv(
        config$optional_inputs$pathway_mapping,
        show_col_types = FALSE
      )
    }
  }

  # Load internal standards
  if (!is.null(config$optional_inputs$internal_standards)) {
    if (file.exists(config$optional_inputs$internal_standards)) {
      log_message("Loading internal standards...")
      annotations$internal_standards <- readr::read_csv(
        config$optional_inputs$internal_standards,
        show_col_types = FALSE
      )
    }
  }

  # Load GMT file
  if (!is.null(config$optional_inputs$gmt_file)) {
    if (file.exists(config$optional_inputs$gmt_file)) {
      log_message("Loading GMT file...")
      annotations$metabolite_sets <- read_gmt(config$optional_inputs$gmt_file)
    }
  }

  if (length(annotations) == 0) {
    return(NULL)
  }

  return(annotations)
}

#' Create input summary
#'
#' @param aligned_data Aligned data list
#' @param sample_roles Sample roles list
#' @param config Configuration list
#' @return Summary list
create_input_summary <- function(aligned_data, sample_roles, config) {
  mat <- aligned_data$matrix
  metadata <- aligned_data$metadata
  condition_col <- config$design$condition_column

  # Per-condition summary (biological samples only)
  bio_samples <- sample_roles$biological_samples
  conditions <- unique(metadata[[condition_col]][metadata[[config$input$sample_id_column]] %in% bio_samples])

  condition_summary <- lapply(conditions, function(cond) {
    sample_col <- config$input$sample_id_column
    samples <- metadata[[sample_col]][metadata[[condition_col]] == cond]
    samples <- intersect(samples, bio_samples)
    samples <- intersect(samples, colnames(mat))

    if (length(samples) == 0) return(NULL)

    sub_mat <- mat[, samples, drop = FALSE]

    data.frame(
      condition = as.character(cond),
      n_samples = length(samples),
      mean_features_detected = mean(colSums(!is.na(sub_mat))),
      median_total_signal = median(colSums(sub_mat, na.rm = TRUE)),
      pct_missing = mean(is.na(sub_mat)) * 100,
      stringsAsFactors = FALSE
    )
  })

  condition_df <- do.call(rbind, condition_summary[!sapply(condition_summary, is.null)])

  # Overall summary
  overall <- data.frame(
    metric = c("Total features", "Total samples", "Biological samples",
               "QC samples", "Blank samples", "Overall % missing",
               "Data type"),
    value = c(
      nrow(mat),
      ncol(mat),
      length(sample_roles$biological_samples),
      length(sample_roles$qc_samples),
      length(sample_roles$blank_samples),
      round(mean(is.na(mat)) * 100, 2),
      config$input$data_type
    ),
    stringsAsFactors = FALSE
  )

  list(
    overall = overall,
    by_condition = condition_df
  )
}
