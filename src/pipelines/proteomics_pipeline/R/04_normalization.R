# =============================================================================
# Normalization and Transformation
# =============================================================================

#' Normalize and transform intensity data
#'
#' @param filtered_data Data from filtering step
#' @param config Configuration list
#' @return Normalized data list
normalize_data <- function(filtered_data, config) {
  log_message("=== Starting Normalization ===")

  mat <- filtered_data$matrix
  metadata <- filtered_data$metadata

  # Store pre-normalized data for diagnostics
  mat_pre <- mat

  # Step 1: Log2 transformation (if not already done)
  if (config$processing$log_transform) {
    mat <- log2_transform(mat)
  }

  mat_log <- mat  # Store log-transformed (pre-normalization)

  # Step 2: Normalization
  norm_method <- config$processing$normalization_method
  mat <- apply_normalization(mat, norm_method, config)

  # Create diagnostic plots
  plots <- create_normalization_plots(mat_pre, mat_log, mat, metadata, config)

  # Save normalized matrix
  norm_df <- as.data.frame(mat)
  norm_df <- tibble::rownames_to_column(norm_df, "feature_id")
  save_table(norm_df, "normalized_matrix.csv", config, "tables")

  log_message("=== Normalization Complete ===")

  list(
    matrix = mat,
    matrix_pre_norm = mat_log,
    metadata = metadata,
    annotations = filtered_data$annotations,
    normalization_method = norm_method,
    plots = plots
  )
}

#' Apply log2 transformation
#'
#' @param mat Intensity matrix
#' @return Log2-transformed matrix
log2_transform <- function(mat) {
  # Check if already log-transformed (values typically < 30 for log2)
  max_val <- max(mat, na.rm = TRUE)

  if (max_val > 50) {
    log_message("Applying log2 transformation (max value: ", round(max_val, 2), ")")
    mat <- log2(mat)
  } else {
    log_message("Data appears to be already log-transformed (max value: ", round(max_val, 2), ")")
  }

  return(mat)
}

#' Apply normalization method
#'
#' @param mat Log-transformed matrix
#' @param method Normalization method
#' @param config Configuration list
#' @return Normalized matrix
apply_normalization <- function(mat, method, config) {
  log_message("Applying normalization method: ", method)

  normalized <- switch(method,
    "vsn" = normalize_vsn(mat),
    "median" = normalize_median(mat),
    "quantile" = normalize_quantile(mat),
    "cyclicloess" = normalize_cyclicloess(mat),
    "none" = mat,
    {
      warning("Unknown normalization method '", method, "'. Using median normalization.")
      normalize_median(mat)
    }
  )

  return(normalized)
}

#' VSN normalization
#'
#' Variance Stabilizing Normalization - recommended for proteomics data
#' Addresses mean-variance relationship
#'
#' @param mat Log-transformed matrix
#' @return VSN-normalized matrix
normalize_vsn <- function(mat) {
  if (!requireNamespace("vsn", quietly = TRUE)) {
    warning("vsn package not available. Falling back to median normalization.")
    return(normalize_median(mat))
  }

  # VSN works better on non-log data, so we need to back-transform first
  mat_linear <- 2^mat

  # Replace NA with minimum for VSN fitting, then restore
  mat_for_vsn <- mat_linear
  min_val <- min(mat_for_vsn, na.rm = TRUE) / 2
  mat_for_vsn[is.na(mat_for_vsn)] <- min_val

  # Fit VSN model
  vsn_fit <- vsn::vsn2(mat_for_vsn)
  mat_vsn <- vsn::predict(vsn_fit, mat_for_vsn)

  # Restore NAs
  mat_vsn[is.na(mat_linear)] <- NA

  # VSN output is already log2-like (glog2)
  log_message("VSN normalization complete")

  return(mat_vsn)
}

#' Median normalization
#'
#' Centers each sample to have the same median
#'
#' @param mat Log-transformed matrix
#' @return Median-normalized matrix
normalize_median <- function(mat) {
  # Calculate sample medians
  sample_medians <- apply(mat, 2, median, na.rm = TRUE)
  global_median <- median(sample_medians)

  # Subtract sample median, add global median
  mat_norm <- sweep(mat, 2, sample_medians, "-")
  mat_norm <- mat_norm + global_median

  log_message("Median normalization complete")

  return(mat_norm)
}

#' Quantile normalization
#'
#' Forces all samples to have the same distribution
#'
#' @param mat Log-transformed matrix
#' @return Quantile-normalized matrix
normalize_quantile <- function(mat) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    warning("limma package not available. Falling back to median normalization.")
    return(normalize_median(mat))
  }

  # limma's normalizeQuantiles handles NAs
  mat_norm <- limma::normalizeQuantiles(mat)

  log_message("Quantile normalization complete")

  return(mat_norm)
}

#' Cyclic loess normalization
#'
#' Iterative loess normalization, good for samples with many differences
#'
#' @param mat Log-transformed matrix
#' @return Loess-normalized matrix
normalize_cyclicloess <- function(mat) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    warning("limma package not available. Falling back to median normalization.")
    return(normalize_median(mat))
  }

  # Cyclic loess normalization
  mat_norm <- limma::normalizeCyclicLoess(mat, method = "fast")

  log_message("Cyclic loess normalization complete")

  return(mat_norm)
}

#' Create normalization diagnostic plots
#'
#' @param mat_raw Raw (untransformed) matrix
#' @param mat_log Log-transformed matrix
#' @param mat_norm Normalized matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return List of ggplot objects
create_normalization_plots <- function(mat_raw, mat_log, mat_norm, metadata, config) {
  plots <- list()

  condition_col <- config$design$condition_column
  sample_col <- config$input$sample_id_column

  # Prepare data for plotting
  sample_info <- data.frame(
    sample = colnames(mat_norm),
    condition = metadata[[condition_col]][match(colnames(mat_norm), metadata[[sample_col]])]
  )

  # 1. Boxplot comparison (before/after normalization)
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
  n_samples <- ncol(mat_norm)
  if (n_samples > 30) {
    # Sample subset for visualization
    sample_subset <- sample(unique(boxplot_data$sample), min(30, n_samples))
    boxplot_data <- boxplot_data[boxplot_data$sample %in% sample_subset, ]
  }

  plots$boxplot <- ggplot2::ggplot(boxplot_data,
                                   ggplot2::aes(x = sample, y = intensity, fill = condition)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::facet_wrap(~ stage, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = "Intensity Distribution by Sample",
      x = "Sample",
      y = "Log2 Intensity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6),
      legend.position = "bottom"
    )

  # 2. Density plots (before/after)
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
    ggplot2::geom_density(alpha = 0.3, show.legend = FALSE) +
    ggplot2::facet_wrap(~ stage, ncol = 2) +
    ggplot2::labs(
      title = "Intensity Density Distribution",
      x = "Log2 Intensity",
      y = "Density"
    ) +
    ggplot2::theme_minimal()

  # 3. Mean-variance plot (for normalized data)
  mean_intensity <- rowMeans(mat_norm, na.rm = TRUE)
  var_intensity <- apply(mat_norm, 1, var, na.rm = TRUE)

  mv_data <- data.frame(
    mean = mean_intensity,
    variance = var_intensity
  )
  mv_data <- mv_data[is.finite(mv_data$mean) & is.finite(mv_data$variance), ]

  plots$mean_variance <- ggplot2::ggplot(mv_data, ggplot2::aes(x = mean, y = variance)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.5) +
    ggplot2::geom_smooth(method = "loess", color = "red", se = FALSE) +
    ggplot2::labs(
      title = "Mean-Variance Relationship (Normalized Data)",
      x = "Mean Log2 Intensity",
      y = "Variance"
    ) +
    ggplot2::theme_minimal()

  # Save plots
  save_plot(plots$boxplot, "normalization_boxplot.png", config, width = 12, height = 10, subdir = "qc")
  save_plot(plots$density, "normalization_density.png", config, width = 12, height = 6, subdir = "qc")
  save_plot(plots$mean_variance, "mean_variance_plot.png", config, subdir = "qc")

  return(plots)
}

#' Calculate normalization statistics
#'
#' @param mat_before Matrix before normalization
#' @param mat_after Matrix after normalization
#' @return Statistics data frame
calculate_norm_stats <- function(mat_before, mat_after) {
  # Per-sample statistics
  stats_before <- data.frame(
    sample = colnames(mat_before),
    median_before = apply(mat_before, 2, median, na.rm = TRUE),
    mean_before = colMeans(mat_before, na.rm = TRUE),
    sd_before = apply(mat_before, 2, sd, na.rm = TRUE)
  )

  stats_after <- data.frame(
    sample = colnames(mat_after),
    median_after = apply(mat_after, 2, median, na.rm = TRUE),
    mean_after = colMeans(mat_after, na.rm = TRUE),
    sd_after = apply(mat_after, 2, sd, na.rm = TRUE)
  )

  stats <- merge(stats_before, stats_after, by = "sample")

  return(stats)
}
