# =============================================================================
# Batch Correction and Drift Correction
# =============================================================================

#' Apply batch/drift correction
#'
#' @param normalized_data Data from normalization step
#' @param config Configuration list
#' @return Batch-corrected data list
correct_batch_effects <- function(normalized_data, config) {
  log_message("=== Starting Batch Correction ===")

  method <- config$batch_correction$method %||% "none"

  if (method == "none") {
    log_message("Batch correction disabled. Skipping.")
    return(list(
      matrix = normalized_data$matrix,
      metadata = normalized_data$metadata,
      sample_roles = normalized_data$sample_roles,
      method = "none",
      corrected = FALSE
    ))
  }

  mat <- normalized_data$matrix
  metadata <- normalized_data$metadata
  sample_roles <- normalized_data$sample_roles

  # Store pre-correction matrix
  mat_pre <- mat

  # Apply correction method
  if (method == "combat") {
    result <- apply_combat(mat, metadata, config)
  } else if (method == "removeBatchEffect") {
    result <- apply_remove_batch_effect(mat, metadata, config)
  } else if (method == "qcrlsc") {
    result <- apply_qcrlsc(mat, metadata, sample_roles, config)
  } else {
    log_message("Unknown batch correction method: ", method, ". Skipping.")
    return(list(
      matrix = mat,
      metadata = metadata,
      sample_roles = sample_roles,
      method = method,
      corrected = FALSE
    ))
  }

  mat_corrected <- result$matrix

  # Create diagnostic plots
  plots <- create_batch_correction_plots(mat_pre, mat_corrected, metadata, config)

  # Save corrected matrix
  corrected_df <- as.data.frame(mat_corrected)
  corrected_df <- tibble::rownames_to_column(corrected_df, "feature_id")
  save_table(corrected_df, "batch_corrected_matrix.csv", config, "tables")

  log_message("=== Batch Correction Complete ===")

  list(
    matrix = mat_corrected,
    matrix_pre_correction = mat_pre,
    metadata = metadata,
    sample_roles = sample_roles,
    method = method,
    corrected = TRUE,
    correction_details = result$details,
    plots = plots
  )
}

#' Apply ComBat batch correction
#'
#' @param mat Normalized matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Corrected result list
apply_combat <- function(mat, metadata, config) {
  batch_col <- config$design$batch_column
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  if (is.null(batch_col) || !batch_col %in% colnames(metadata)) {
    log_message("Batch column not found. Skipping ComBat.")
    return(list(matrix = mat, details = list(applied = FALSE)))
  }

  if (!requireNamespace("sva", quietly = TRUE)) {
    log_message("sva package not available. Skipping ComBat.")
    return(list(matrix = mat, details = list(applied = FALSE)))
  }

  log_message("Applying ComBat batch correction...")

  # Get batch info
  batch <- metadata[[batch_col]][match(colnames(mat), metadata[[sample_col]])]

  # Create model matrix for condition (to preserve)
  condition <- metadata[[condition_col]][match(colnames(mat), metadata[[sample_col]])]
  mod <- model.matrix(~ condition)

  # Handle missing values (ComBat doesn't handle NAs well)
  mat_imputed <- mat
  for (i in seq_len(nrow(mat_imputed))) {
    row <- mat_imputed[i, ]
    if (any(is.na(row))) {
      mat_imputed[i, is.na(row)] <- min(row, na.rm = TRUE) / 2
    }
  }

  # Apply ComBat
  parametric <- config$batch_correction$combat_parametric %||% TRUE

  tryCatch({
    mat_corrected <- sva::ComBat(
      dat = mat_imputed,
      batch = batch,
      mod = mod,
      par.prior = parametric
    )

    # Restore NAs
    mat_corrected[is.na(mat)] <- NA

    log_message("ComBat correction complete")

    return(list(
      matrix = mat_corrected,
      details = list(
        applied = TRUE,
        n_batches = length(unique(batch)),
        parametric = parametric
      )
    ))
  }, error = function(e) {
    log_message("ComBat failed: ", e$message)
    return(list(matrix = mat, details = list(applied = FALSE, error = e$message)))
  })
}

#' Apply limma removeBatchEffect
#'
#' @param mat Normalized matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Corrected result list
apply_remove_batch_effect <- function(mat, metadata, config) {
  batch_col <- config$design$batch_column
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  if (is.null(batch_col) || !batch_col %in% colnames(metadata)) {
    log_message("Batch column not found. Skipping removeBatchEffect.")
    return(list(matrix = mat, details = list(applied = FALSE)))
  }

  if (!requireNamespace("limma", quietly = TRUE)) {
    log_message("limma package not available. Skipping removeBatchEffect.")
    return(list(matrix = mat, details = list(applied = FALSE)))
  }

  log_message("Applying limma::removeBatchEffect...")

  # Get batch and design info
  batch <- metadata[[batch_col]][match(colnames(mat), metadata[[sample_col]])]
  condition <- metadata[[condition_col]][match(colnames(mat), metadata[[sample_col]])]
  design <- model.matrix(~ condition)

  # Handle NAs
  mat_imputed <- mat
  for (i in seq_len(nrow(mat_imputed))) {
    row <- mat_imputed[i, ]
    if (any(is.na(row))) {
      mat_imputed[i, is.na(row)] <- min(row, na.rm = TRUE) / 2
    }
  }

  tryCatch({
    mat_corrected <- limma::removeBatchEffect(mat_imputed, batch = batch, design = design)

    # Restore NAs
    mat_corrected[is.na(mat)] <- NA

    log_message("removeBatchEffect correction complete")

    return(list(
      matrix = mat_corrected,
      details = list(
        applied = TRUE,
        n_batches = length(unique(batch))
      )
    ))
  }, error = function(e) {
    log_message("removeBatchEffect failed: ", e$message)
    return(list(matrix = mat, details = list(applied = FALSE, error = e$message)))
  })
}

#' Apply QC-RLSC (QC-based LOESS signal correction)
#'
#' @param mat Normalized matrix
#' @param metadata Sample metadata
#' @param sample_roles Sample roles
#' @param config Configuration list
#' @return Corrected result list
apply_qcrlsc <- function(mat, metadata, sample_roles, config) {
  inj_col <- config$design$injection_order_column
  sample_col <- config$input$sample_id_column

  # Check requirements
  if (length(sample_roles$qc_samples) == 0) {
    log_message("No QC samples found. Skipping QC-RLSC.")
    return(list(matrix = mat, details = list(applied = FALSE, reason = "No QC samples")))
  }

  if (is.null(inj_col) || !inj_col %in% colnames(metadata)) {
    log_message("Injection order not found. Skipping QC-RLSC.")
    return(list(matrix = mat, details = list(applied = FALSE, reason = "No injection order")))
  }

  log_message("Applying QC-RLSC drift correction...")

  # Get injection order
  injection_order <- metadata[[inj_col]][match(colnames(mat), metadata[[sample_col]])]

  # Identify QC samples
  qc_idx <- which(colnames(mat) %in% sample_roles$qc_samples)

  if (length(qc_idx) < 4) {
    log_message("Too few QC samples (", length(qc_idx), ") for LOESS. Skipping QC-RLSC.")
    return(list(matrix = mat, details = list(applied = FALSE, reason = "Too few QC samples")))
  }

  span <- config$batch_correction$loess_span %||% 0.75

  # Apply LOESS correction per feature
  mat_corrected <- mat

  for (i in seq_len(nrow(mat))) {
    feature_vals <- mat[i, ]
    qc_vals <- feature_vals[qc_idx]
    qc_order <- injection_order[qc_idx]

    # Skip if too many missing in QC
    if (sum(!is.na(qc_vals)) < 4) next

    # Fit LOESS on QC samples
    tryCatch({
      loess_fit <- loess(qc_vals ~ qc_order, span = span, na.action = na.exclude)

      # Predict for all samples
      predicted <- predict(loess_fit, newdata = injection_order)

      # Calculate reference (median of QC)
      qc_median <- median(qc_vals, na.rm = TRUE)

      # Apply correction
      correction_factor <- qc_median / predicted
      correction_factor[!is.finite(correction_factor)] <- 1

      mat_corrected[i, ] <- feature_vals * correction_factor
    }, error = function(e) {
      # Keep original if LOESS fails
    })
  }

  log_message("QC-RLSC drift correction complete")

  return(list(
    matrix = mat_corrected,
    details = list(
      applied = TRUE,
      n_qc_samples = length(qc_idx),
      loess_span = span
    )
  ))
}

#' Create batch correction diagnostic plots
#'
#' @param mat_pre Matrix before correction
#' @param mat_post Matrix after correction
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return List of plots
create_batch_correction_plots <- function(mat_pre, mat_post, metadata, config) {
  plots <- list()

  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  batch_col <- config$design$batch_column

  # Run PCA before and after
  pca_pre <- run_pca_simple(mat_pre, metadata, config)
  pca_post <- run_pca_simple(mat_post, metadata, config)

  # Combine for plotting
  pca_combined <- rbind(
    data.frame(pca_pre$coordinates, stage = "Before correction"),
    data.frame(pca_post$coordinates, stage = "After correction")
  )

  # PCA colored by condition
  plots$pca_condition <- ggplot2::ggplot(pca_combined,
    ggplot2::aes(x = PC1, y = PC2, color = condition, shape = stage)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::labs(title = "PCA Before/After Batch Correction (by Condition)",
                  x = paste0("PC1 (", round(pca_post$variance_explained[1], 1), "%)"),
                  y = paste0("PC2 (", round(pca_post$variance_explained[2], 1), "%)")) +
    ggplot2::theme_minimal()

  # PCA colored by batch if available
  if (!is.null(batch_col) && batch_col %in% colnames(metadata)) {
    pca_pre$coordinates$batch <- metadata[[batch_col]][match(pca_pre$coordinates$sample_id, metadata[[sample_col]])]
    pca_post$coordinates$batch <- metadata[[batch_col]][match(pca_post$coordinates$sample_id, metadata[[sample_col]])]

    pca_combined_batch <- rbind(
      data.frame(pca_pre$coordinates, stage = "Before correction"),
      data.frame(pca_post$coordinates, stage = "After correction")
    )

    plots$pca_batch <- ggplot2::ggplot(pca_combined_batch,
      ggplot2::aes(x = PC1, y = PC2, color = batch, shape = condition)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::facet_wrap(~ stage) +
      ggplot2::labs(title = "PCA Before/After Batch Correction (by Batch)") +
      ggplot2::theme_minimal()

    save_plot(plots$pca_batch, "batch_correction_pca_batch.png", config, subdir = "qc")
  }

  save_plot(plots$pca_condition, "batch_correction_pca_condition.png", config, subdir = "qc")

  return(plots)
}

#' Simple PCA for batch correction visualization
#'
#' @param mat Matrix
#' @param metadata Metadata
#' @param config Configuration
#' @return PCA results
run_pca_simple <- function(mat, metadata, config) {
  # Remove constant features
  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat_filtered <- mat[vars > 0 & !is.na(vars), ]

  # Simple imputation
  for (i in seq_len(nrow(mat_filtered))) {
    row <- mat_filtered[i, ]
    if (any(is.na(row))) {
      mat_filtered[i, is.na(row)] <- median(row, na.rm = TRUE)
    }
  }

  pca <- prcomp(t(mat_filtered), center = TRUE, scale. = TRUE)

  n_comp <- min(5, ncol(pca$x))
  coords <- as.data.frame(pca$x[, 1:n_comp])
  coords$sample_id <- rownames(coords)

  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  coords$condition <- metadata[[condition_col]][match(coords$sample_id, metadata[[sample_col]])]

  var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

  list(
    coordinates = coords,
    variance_explained = var_explained[1:n_comp]
  )
}
