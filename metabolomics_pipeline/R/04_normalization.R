# =============================================================================
# Normalization and Transformation
# =============================================================================

#' Normalize and transform data
#'
#' @param filtered_data Data from filtering step
#' @param config Configuration list
#' @return Normalized data list
normalize_data <- function(filtered_data, config) {
  log_message("=== Starting Normalization ===")

  mat <- filtered_data$matrix
  metadata <- filtered_data$metadata
  sample_roles <- filtered_data$sample_roles

  # Store pre-normalized for comparison
  mat_raw <- mat

  # Step 1: Log transformation
  transform_method <- config$processing$transform %||% "log2"
  pseudocount <- config$processing$pseudocount %||% 1

  if (transform_method != "none") {
    mat <- apply_transformation(mat, transform_method, pseudocount)
  }

  mat_log <- mat  # Store post-transform, pre-normalization

  # Step 2: Normalization
  norm_method <- config$processing$normalization_method %||% "PQN"
  mat <- apply_normalization(mat, norm_method, sample_roles, config)

  # Create diagnostic plots
  plots <- create_normalization_plots(mat_raw, mat_log, mat, metadata, config)

  # Save normalized matrix
  norm_df <- as.data.frame(mat)
  norm_df <- tibble::rownames_to_column(norm_df, "feature_id")
  save_table(norm_df, "normalized_matrix.csv", config, "tables")

  log_message("=== Normalization Complete ===")

  list(
    matrix = mat,
    matrix_raw = mat_raw,
    matrix_log = mat_log,
    metadata = metadata,
    sample_roles = sample_roles,
    transform_method = transform_method,
    normalization_method = norm_method,
    plots = plots
  )
}

#' Apply log transformation
#'
#' @param mat Intensity matrix
#' @param method Transformation method
#' @param pseudocount Pseudocount to add before log
#' @return Transformed matrix
apply_transformation <- function(mat, method, pseudocount = 1) {
  log_message("Applying ", method, " transformation (pseudocount = ", pseudocount, ")")

  # Check if data looks already log-transformed
  max_val <- max(mat, na.rm = TRUE)
  if (max_val < 30) {
    log_message("Data appears already log-transformed (max = ", round(max_val, 2), "). Skipping transformation.")
    return(mat)
  }

  if (method == "log2") {
    mat <- log2(mat + pseudocount)
  } else if (method == "log10") {
    mat <- log10(mat + pseudocount)
  }

  return(mat)
}

#' Apply normalization method
#'
#' @param mat Log-transformed matrix
#' @param method Normalization method
#' @param sample_roles Sample roles
#' @param config Configuration list
#' @return Normalized matrix
apply_normalization <- function(mat, method, sample_roles, config) {
  log_message("Applying ", method, " normalization")

  normalized <- switch(method,
    "PQN" = normalize_pqn(mat),
    "median" = normalize_median(mat),
    "TIC" = normalize_tic(mat),
    "quantile" = normalize_quantile(mat),
    "vsn" = normalize_vsn(mat),
    "none" = mat,
    {
      warning("Unknown normalization method '", method, "'. Using median.")
      normalize_median(mat)
    }
  )

  return(normalized)
}

#' Probabilistic Quotient Normalization (PQN)
#'
#' Standard normalization for metabolomics data
#'
#' @param mat Log-transformed matrix
#' @return PQN-normalized matrix
normalize_pqn <- function(mat) {
  # Step 1: Calculate reference spectrum (median of all samples)
  reference <- apply(mat, 1, median, na.rm = TRUE)

  # Step 2: Calculate quotients for each sample
  mat_norm <- mat

  for (j in seq_len(ncol(mat))) {
    sample <- mat[, j]
    quotients <- sample / reference
    quotients <- quotients[is.finite(quotients) & quotients > 0]

    if (length(quotients) > 0) {
      scaling_factor <- median(quotients, na.rm = TRUE)
      mat_norm[, j] <- sample / scaling_factor
    }
  }

  log_message("PQN normalization complete")
  return(mat_norm)
}

#' Median normalization
#'
#' @param mat Log-transformed matrix
#' @return Median-normalized matrix
normalize_median <- function(mat) {
  sample_medians <- apply(mat, 2, median, na.rm = TRUE)
  global_median <- median(sample_medians)

  mat_norm <- sweep(mat, 2, sample_medians, "-")
  mat_norm <- mat_norm + global_median

  log_message("Median normalization complete")
  return(mat_norm)
}

#' Total Ion Current (TIC) normalization
#'
#' @param mat Matrix (can be raw or log)
#' @return TIC-normalized matrix
normalize_tic <- function(mat) {
  # If log-transformed, back-transform for TIC
  if (max(mat, na.rm = TRUE) < 50) {
    mat_linear <- 2^mat
    is_log <- TRUE
  } else {
    mat_linear <- mat
    is_log <- FALSE
  }

  # Calculate TIC per sample
  tic <- colSums(mat_linear, na.rm = TRUE)
  median_tic <- median(tic)

  # Scale each sample
  scaling_factors <- median_tic / tic
  mat_norm <- sweep(mat_linear, 2, scaling_factors, "*")

  # Re-log if needed
  if (is_log) {
    mat_norm <- log2(mat_norm)
  }

  log_message("TIC normalization complete")
  return(mat_norm)
}

#' Quantile normalization
#'
#' @param mat Log-transformed matrix
#' @return Quantile-normalized matrix
normalize_quantile <- function(mat) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    warning("limma package not available. Falling back to median normalization.")
    return(normalize_median(mat))
  }

  mat_norm <- limma::normalizeQuantiles(mat)

  log_message("Quantile normalization complete")
  return(mat_norm)
}

#' VSN normalization
#'
#' @param mat Matrix
#' @return VSN-normalized matrix
normalize_vsn <- function(mat) {
  if (!requireNamespace("vsn", quietly = TRUE)) {
    warning("vsn package not available. Falling back to median normalization.")
    return(normalize_median(mat))
  }

  # VSN works on raw intensities
  if (max(mat, na.rm = TRUE) < 50) {
    mat_linear <- 2^mat
  } else {
    mat_linear <- mat
  }

  # Handle NAs
  mat_for_vsn <- mat_linear
  min_val <- min(mat_for_vsn, na.rm = TRUE) / 2
  mat_for_vsn[is.na(mat_for_vsn)] <- min_val

  tryCatch({
    vsn_fit <- vsn::vsn2(mat_for_vsn)
    mat_vsn <- vsn::predict(vsn_fit, mat_for_vsn)

    # Restore NAs
    mat_vsn[is.na(mat_linear)] <- NA

    log_message("VSN normalization complete")
    return(mat_vsn)
  }, error = function(e) {
    warning("VSN failed: ", e$message, ". Falling back to median normalization.")
    return(normalize_median(mat))
  })
}

#' Create normalization diagnostic plots
#'
#' @param mat_raw Raw matrix
#' @param mat_log Log-transformed matrix
#' @param mat_norm Normalized matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return List of plots
create_normalization_plots <- function(mat_raw, mat_log, mat_norm, metadata, config) {
  plots <- list()

  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  # Prepare sample info
  sample_info <- data.frame(
    sample = colnames(mat_norm),
    condition = metadata[[condition_col]][match(colnames(mat_norm), metadata[[sample_col]])]
  )

  # 1. Boxplot comparison
  boxplot_data <- rbind(
    data.frame(
      intensity = as.vector(mat_log),
      sample = rep(colnames(mat_log), each = nrow(mat_log)),
      stage = "Before normalization"
    ),
    data.frame(
      intensity = as.vector(mat_norm),
      sample = rep(colnames(mat_norm), each = nrow(mat_norm)),
      stage = "After normalization"
    )
  )
  boxplot_data <- boxplot_data[!is.na(boxplot_data$intensity), ]
  boxplot_data$condition <- sample_info$condition[match(boxplot_data$sample, sample_info$sample)]

  # Limit samples for readability
  if (length(unique(boxplot_data$sample)) > 30) {
    sample_subset <- sample(unique(boxplot_data$sample), 30)
    boxplot_data <- boxplot_data[boxplot_data$sample %in% sample_subset, ]
  }

  plots$boxplot <- ggplot2::ggplot(boxplot_data,
    ggplot2::aes(x = sample, y = intensity, fill = condition)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::facet_wrap(~ stage, ncol = 1, scales = "free_y") +
    ggplot2::labs(title = "Intensity Distribution Before/After Normalization",
                  x = "Sample", y = "Log2 Intensity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
      legend.position = "bottom"
    )

  # 2. Density plots
  density_data <- rbind(
    data.frame(
      intensity = as.vector(mat_log),
      sample = rep(colnames(mat_log), each = nrow(mat_log)),
      stage = "Before normalization"
    ),
    data.frame(
      intensity = as.vector(mat_norm),
      sample = rep(colnames(mat_norm), each = nrow(mat_norm)),
      stage = "After normalization"
    )
  )
  density_data <- density_data[!is.na(density_data$intensity), ]

  plots$density <- ggplot2::ggplot(density_data,
    ggplot2::aes(x = intensity, color = sample)) +
    ggplot2::geom_density(show.legend = FALSE, alpha = 0.3) +
    ggplot2::facet_wrap(~ stage, ncol = 2) +
    ggplot2::labs(title = "Intensity Density Distribution",
                  x = "Log2 Intensity", y = "Density") +
    ggplot2::theme_minimal()

  # 3. Total signal comparison
  total_before <- colSums(2^mat_log, na.rm = TRUE)
  total_after <- colSums(2^mat_norm, na.rm = TRUE)

  total_df <- data.frame(
    sample = c(names(total_before), names(total_after)),
    total_signal = c(total_before, total_after),
    stage = rep(c("Before", "After"), each = length(total_before))
  )

  plots$total_signal <- ggplot2::ggplot(total_df,
    ggplot2::aes(x = sample, y = log10(total_signal), fill = stage)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = "Total Signal Before/After Normalization",
                  x = "Sample", y = "Log10(Total Signal)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))

  # Save plots
  save_plot(plots$boxplot, "normalization_boxplot.png", config, width = 12, height = 10, subdir = "qc")
  save_plot(plots$density, "normalization_density.png", config, width = 12, height = 6, subdir = "qc")
  save_plot(plots$total_signal, "normalization_total_signal.png", config, subdir = "qc")

  return(plots)
}
