# =============================================================================
# Initial Quality Control Metrics
# =============================================================================

#' Run initial QC analysis on raw data
#'
#' @param ingested_data Data from ingestion step
#' @param config Configuration list
#' @return QC metrics list
run_initial_qc <- function(ingested_data, config) {
  log_message("=== Starting Initial QC ===")

  mat <- ingested_data$matrix
  metadata <- ingested_data$metadata
  sample_roles <- ingested_data$sample_roles

  # Calculate per-sample metrics
  sample_metrics <- calculate_sample_metrics(mat, metadata, sample_roles, config)

  # Calculate per-feature metrics
  feature_metrics <- calculate_feature_metrics(mat, sample_roles, config)

  # Calculate intensity distributions
  intensity_stats <- calculate_intensity_stats(mat, config)

  # Run PCA on raw data
  pca_raw <- run_pca_qc(mat, metadata, config, label = "raw")

  # Calculate sample correlations
  correlation_metrics <- calculate_correlation_metrics(mat, config)

  # Injection order analysis (if available)
  drift_analysis <- NULL
  if (!is.null(config$design$injection_order_column)) {
    drift_analysis <- analyze_injection_drift(mat, metadata, config)
  }

  # Internal standards analysis (if available)
  is_analysis <- NULL
  if (!is.null(ingested_data$annotations$internal_standards)) {
    is_analysis <- analyze_internal_standards(mat, ingested_data$annotations$internal_standards, config)
  }

  # Create QC plots
  plots <- create_initial_qc_plots(mat, sample_metrics, feature_metrics, metadata, config)

  # Save metrics
  save_table(sample_metrics, "qc_sample_metrics.csv", config, "qc")
  save_table(feature_metrics, "qc_feature_metrics.csv", config, "qc")

  log_message("=== Initial QC Complete ===")

  list(
    sample_metrics = sample_metrics,
    feature_metrics = feature_metrics,
    intensity_stats = intensity_stats,
    pca = pca_raw,
    correlations = correlation_metrics,
    drift_analysis = drift_analysis,
    internal_standards = is_analysis,
    plots = plots
  )
}

#' Calculate per-sample QC metrics
#'
#' @param mat Feature matrix
#' @param metadata Sample metadata
#' @param sample_roles Sample roles list
#' @param config Configuration list
#' @return Sample metrics data frame
calculate_sample_metrics <- function(mat, metadata, sample_roles, config) {
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  metrics <- data.frame(
    sample_id = colnames(mat),
    stringsAsFactors = FALSE
  )

  # Total ion signal
  metrics$total_signal <- colSums(mat, na.rm = TRUE)

  # Number of detected features
  metrics$n_detected <- colSums(!is.na(mat))
  metrics$pct_detected <- metrics$n_detected / nrow(mat) * 100

  # Missing values
  metrics$n_missing <- colSums(is.na(mat))
  metrics$pct_missing <- metrics$n_missing / nrow(mat) * 100

  # Mean and median intensity
  metrics$mean_intensity <- colMeans(mat, na.rm = TRUE)
  metrics$median_intensity <- apply(mat, 2, median, na.rm = TRUE)

  # Add metadata info
  match_idx <- match(metrics$sample_id, metadata[[sample_col]])
  metrics$condition <- metadata[[condition_col]][match_idx]

  if (sample_roles$has_sample_types) {
    sample_type_col <- config$sample_types$sample_type_column
    metrics$sample_type <- metadata[[sample_type_col]][match_idx]
  } else {
    metrics$sample_type <- "Sample"
  }

  # Injection order if available
  if (!is.null(config$design$injection_order_column)) {
    inj_col <- config$design$injection_order_column
    if (inj_col %in% colnames(metadata)) {
      metrics$injection_order <- metadata[[inj_col]][match_idx]
    }
  }

  return(metrics)
}

#' Calculate per-feature QC metrics
#'
#' @param mat Feature matrix
#' @param sample_roles Sample roles list
#' @param config Configuration list
#' @return Feature metrics data frame
calculate_feature_metrics <- function(mat, sample_roles, config) {
  metrics <- data.frame(
    feature_id = rownames(mat),
    stringsAsFactors = FALSE
  )

  # Overall statistics
  metrics$mean_intensity <- rowMeans(mat, na.rm = TRUE)
  metrics$median_intensity <- apply(mat, 1, median, na.rm = TRUE)
  metrics$sd_intensity <- apply(mat, 1, sd, na.rm = TRUE)
  metrics$cv <- metrics$sd_intensity / metrics$mean_intensity

  # Missingness
  metrics$n_present <- rowSums(!is.na(mat))
  metrics$n_missing <- rowSums(is.na(mat))
  metrics$pct_present <- metrics$n_present / ncol(mat) * 100

  # QC-specific metrics
  if (length(sample_roles$qc_samples) > 0) {
    qc_mat <- mat[, sample_roles$qc_samples, drop = FALSE]
    metrics$qc_mean <- rowMeans(qc_mat, na.rm = TRUE)
    metrics$qc_sd <- apply(qc_mat, 1, sd, na.rm = TRUE)
    metrics$qc_rsd <- (metrics$qc_sd / metrics$qc_mean) * 100
    metrics$qc_n_present <- rowSums(!is.na(qc_mat))
  }

  # Blank-specific metrics
  if (length(sample_roles$blank_samples) > 0) {
    blank_mat <- mat[, sample_roles$blank_samples, drop = FALSE]
    metrics$blank_mean <- rowMeans(blank_mat, na.rm = TRUE)

    # Blank ratio (blank / sample intensity)
    if (length(sample_roles$biological_samples) > 0) {
      bio_mat <- mat[, sample_roles$biological_samples, drop = FALSE]
      bio_mean <- rowMeans(bio_mat, na.rm = TRUE)
      metrics$blank_ratio <- metrics$blank_mean / bio_mean
      metrics$blank_ratio[is.na(metrics$blank_ratio) | is.infinite(metrics$blank_ratio)] <- 0
    }
  }

  return(metrics)
}

#' Calculate intensity distribution statistics
#'
#' @param mat Feature matrix
#' @param config Configuration list
#' @return Intensity statistics list
calculate_intensity_stats <- function(mat, config) {
  values <- as.vector(mat)
  values <- values[!is.na(values)]

  list(
    min = min(values),
    max = max(values),
    mean = mean(values),
    median = median(values),
    q25 = quantile(values, 0.25),
    q75 = quantile(values, 0.75),
    iqr = IQR(values),
    n_zeros = sum(values == 0),
    n_negative = sum(values < 0)
  )
}

#' Run PCA for QC visualization
#'
#' @param mat Feature matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @param label Label for this PCA (e.g., "raw", "normalized")
#' @return PCA results list
run_pca_qc <- function(mat, metadata, config, label = "raw") {
  log_message("Running PCA on ", label, " data...")

  # Remove features with all NA or zero variance
  valid_features <- apply(mat, 1, function(x) {
    sum(!is.na(x)) >= 3 && var(x, na.rm = TRUE) > 0
  })
  mat_filtered <- mat[valid_features, ]

  # Simple imputation for PCA (half-min)
  for (i in seq_len(nrow(mat_filtered))) {
    row <- mat_filtered[i, ]
    if (any(is.na(row))) {
      min_val <- min(row, na.rm = TRUE) / 2
      mat_filtered[i, is.na(row)] <- min_val
    }
  }

  # Log transform for PCA if not already
  if (max(mat_filtered, na.rm = TRUE) > 100) {
    mat_filtered <- log2(mat_filtered + 1)
  }

  n_components <- min(config$qc$n_pca_components %||% 10, ncol(mat_filtered) - 1, nrow(mat_filtered) - 1)

  pca <- prcomp(t(mat_filtered), center = TRUE, scale. = TRUE)

  # Extract coordinates
  pca_coords <- as.data.frame(pca$x[, 1:n_components])
  pca_coords$sample_id <- rownames(pca_coords)

  # Add metadata
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  pca_coords$condition <- metadata[[condition_col]][match(pca_coords$sample_id, metadata[[sample_col]])]

  sample_type_col <- config$sample_types$sample_type_column
  if (sample_type_col %in% colnames(metadata)) {
    pca_coords$sample_type <- metadata[[sample_type_col]][match(pca_coords$sample_id, metadata[[sample_col]])]
  } else {
    pca_coords$sample_type <- "Sample"
  }

  # Variance explained
  var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
  var_explained <- var_explained[1:n_components]
  names(var_explained) <- paste0("PC", 1:n_components)

  list(
    coordinates = pca_coords,
    variance_explained = var_explained,
    loadings = pca$rotation[, 1:n_components],
    label = label
  )
}

#' Calculate sample correlation metrics
#'
#' @param mat Feature matrix
#' @param config Configuration list
#' @return Correlation metrics list
calculate_correlation_metrics <- function(mat, config) {
  # Simple imputation for correlation
  mat_imputed <- mat
  for (i in seq_len(nrow(mat_imputed))) {
    row <- mat_imputed[i, ]
    if (any(is.na(row))) {
      mat_imputed[i, is.na(row)] <- median(row, na.rm = TRUE)
    }
  }

  cor_mat <- cor(mat_imputed, use = "pairwise.complete.obs", method = "spearman")

  # Per-sample mean correlation
  mean_cor <- sapply(seq_len(ncol(cor_mat)), function(i) {
    mean(cor_mat[i, -i], na.rm = TRUE)
  })
  names(mean_cor) <- colnames(mat)

  list(
    correlation_matrix = cor_mat,
    mean_correlation = mean_cor,
    min_correlation = min(cor_mat[lower.tri(cor_mat)], na.rm = TRUE),
    median_correlation = median(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)
  )
}

#' Analyze injection order drift
#'
#' @param mat Feature matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Drift analysis results
analyze_injection_drift <- function(mat, metadata, config) {
  inj_col <- config$design$injection_order_column
  sample_col <- config$input$sample_id_column

  if (!inj_col %in% colnames(metadata)) {
    log_message("Injection order column not found. Skipping drift analysis.")
    return(NULL)
  }

  log_message("Analyzing injection order drift...")

  injection_order <- metadata[[inj_col]][match(colnames(mat), metadata[[sample_col]])]

  # Total signal vs injection order
  total_signal <- colSums(mat, na.rm = TRUE)

  # Correlation between injection order and total signal
  drift_cor <- cor(injection_order, total_signal, use = "complete.obs", method = "spearman")

  # Per-feature drift correlations
  feature_drift <- apply(mat, 1, function(x) {
    if (sum(!is.na(x)) < 5) return(NA)
    cor(injection_order, x, use = "complete.obs", method = "spearman")
  })

  list(
    injection_order = injection_order,
    total_signal = total_signal,
    total_signal_drift_cor = drift_cor,
    feature_drift_correlations = feature_drift,
    n_features_with_drift = sum(abs(feature_drift) > 0.5, na.rm = TRUE)
  )
}

#' Analyze internal standards
#'
#' @param mat Feature matrix
#' @param is_table Internal standards table
#' @param config Configuration list
#' @return Internal standards analysis
analyze_internal_standards <- function(mat, is_table, config) {
  log_message("Analyzing internal standards...")

  is_ids <- is_table$feature_id
  is_in_data <- is_ids[is_ids %in% rownames(mat)]

  if (length(is_in_data) == 0) {
    log_message("No internal standards found in data")
    return(NULL)
  }

  is_mat <- mat[is_in_data, , drop = FALSE]

  # Calculate CV for each IS
  is_cv <- apply(is_mat, 1, calc_cv, na.rm = TRUE) * 100

  # Calculate mean intensity
  is_mean <- rowMeans(is_mat, na.rm = TRUE)

  results <- data.frame(
    feature_id = is_in_data,
    mean_intensity = is_mean,
    cv_percent = is_cv,
    n_detected = rowSums(!is.na(is_mat)),
    stringsAsFactors = FALSE
  )

  # Add expected values if available
  if ("expected_cv" %in% colnames(is_table)) {
    match_idx <- match(results$feature_id, is_table$feature_id)
    results$expected_cv <- is_table$expected_cv[match_idx]
    results$cv_pass <- results$cv_percent <= results$expected_cv
  }

  save_table(results, "internal_standards_summary.csv", config, "qc")

  list(
    summary = results,
    matrix = is_mat
  )
}

#' Create initial QC plots
#'
#' @param mat Feature matrix
#' @param sample_metrics Sample metrics
#' @param feature_metrics Feature metrics
#' @param metadata Metadata
#' @param config Configuration list
#' @return List of plots
create_initial_qc_plots <- function(mat, sample_metrics, feature_metrics, metadata, config) {
  plots <- list()

  # 1. Total signal barplot
  sample_metrics$sample_id <- factor(sample_metrics$sample_id, levels = sample_metrics$sample_id)

  plots$total_signal <- ggplot2::ggplot(sample_metrics,
    ggplot2::aes(x = sample_id, y = total_signal, fill = sample_type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(title = "Total Signal per Sample",
                  x = "Sample", y = "Total Signal", fill = "Sample Type") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))

  # 2. Features detected per sample
  plots$features_detected <- ggplot2::ggplot(sample_metrics,
    ggplot2::aes(x = sample_id, y = n_detected, fill = sample_type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(title = "Features Detected per Sample",
                  x = "Sample", y = "Number of Features", fill = "Sample Type") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))

  # 3. Intensity distribution (density)
  intensity_df <- data.frame(
    intensity = as.vector(mat),
    sample = rep(colnames(mat), each = nrow(mat))
  )
  intensity_df <- intensity_df[!is.na(intensity_df$intensity), ]

  plots$intensity_density <- ggplot2::ggplot(intensity_df,
    ggplot2::aes(x = log10(intensity + 1), color = sample)) +
    ggplot2::geom_density(show.legend = FALSE, alpha = 0.5) +
    ggplot2::labs(title = "Log10 Intensity Distribution",
                  x = "Log10(Intensity + 1)", y = "Density") +
    ggplot2::theme_minimal()

  # 4. Sample missingness
  plots$sample_missing <- ggplot2::ggplot(sample_metrics,
    ggplot2::aes(x = sample_id, y = pct_missing, fill = sample_type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(title = "Missing Values per Sample",
                  x = "Sample", y = "% Missing", fill = "Sample Type") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))

  # 5. Feature missingness distribution
  plots$feature_missing <- ggplot2::ggplot(feature_metrics,
    ggplot2::aes(x = pct_present)) +
    ggplot2::geom_histogram(bins = 50, fill = "#0072B2", color = "white") +
    ggplot2::labs(title = "Feature Detection Rate Distribution",
                  x = "% Samples with Detection", y = "Number of Features") +
    ggplot2::theme_minimal()

  # Save plots
  save_plot(plots$total_signal, "qc_total_signal.png", config, subdir = "qc")
  save_plot(plots$features_detected, "qc_features_detected.png", config, subdir = "qc")
  save_plot(plots$intensity_density, "qc_intensity_density.png", config, subdir = "qc")
  save_plot(plots$sample_missing, "qc_sample_missing.png", config, subdir = "qc")
  save_plot(plots$feature_missing, "qc_feature_missing.png", config, subdir = "qc")

  return(plots)
}
