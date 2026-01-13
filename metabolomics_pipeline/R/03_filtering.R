# =============================================================================
# Feature Filtering
# =============================================================================

#' Apply feature filters
#'
#' @param ingested_data Data from ingestion step
#' @param initial_qc Initial QC results
#' @param config Configuration list
#' @return Filtered data list
filter_features <- function(ingested_data, initial_qc, config) {
  log_message("=== Starting Feature Filtering ===")

  mat <- ingested_data$matrix
  metadata <- ingested_data$metadata
  sample_roles <- ingested_data$sample_roles
  feature_metrics <- initial_qc$feature_metrics

  # Track filtering steps
  filter_log <- data.frame(
    step = character(),
    features_before = integer(),
    features_removed = integer(),
    features_after = integer(),
    stringsAsFactors = FALSE
  )

  initial_n <- nrow(mat)

  # Step 1: Blank filtering (if blanks available)
  if (length(sample_roles$blank_samples) > 0 && !is.null(config$filtering$blank_ratio_threshold)) {
    result <- filter_by_blank_ratio(mat, feature_metrics, config)
    mat <- result$matrix
    filter_log <- rbind(filter_log, data.frame(
      step = "Blank ratio filter",
      features_before = result$n_before,
      features_removed = result$n_removed,
      features_after = nrow(mat),
      stringsAsFactors = FALSE
    ))
  }

  # Step 2: QC RSD filtering (if QC samples available)
  if (length(sample_roles$qc_samples) > 0 && !is.null(config$filtering$qc_rsd_threshold)) {
    # Recalculate feature metrics for remaining features
    feature_metrics_updated <- calculate_feature_metrics_simple(mat, sample_roles)
    result <- filter_by_qc_rsd(mat, feature_metrics_updated, config)
    mat <- result$matrix
    filter_log <- rbind(filter_log, data.frame(
      step = "QC RSD filter",
      features_before = result$n_before,
      features_removed = result$n_removed,
      features_after = nrow(mat),
      stringsAsFactors = FALSE
    ))
  }

  # Step 3: Global missingness filter
  result <- filter_by_global_presence(mat, config)
  mat <- result$matrix
  filter_log <- rbind(filter_log, data.frame(
    step = "Global presence filter",
    features_before = result$n_before,
    features_removed = result$n_removed,
    features_after = nrow(mat),
    stringsAsFactors = FALSE
  ))

  # Step 4: Group-specific missingness filter
  result <- filter_by_group_presence(mat, metadata, sample_roles, config)
  mat <- result$matrix
  filter_log <- rbind(filter_log, data.frame(
    step = "Group presence filter",
    features_before = result$n_before,
    features_removed = result$n_removed,
    features_after = nrow(mat),
    stringsAsFactors = FALSE
  ))

  # Step 5: Low variance filter (optional)
  if (config$filtering$remove_low_variance %||% FALSE) {
    result <- filter_low_variance(mat, config)
    mat <- result$matrix
    filter_log <- rbind(filter_log, data.frame(
      step = "Low variance filter",
      features_before = result$n_before,
      features_removed = result$n_removed,
      features_after = nrow(mat),
      stringsAsFactors = FALSE
    ))
  }

  # Save filter log
  save_table(filter_log, "filtering_summary.csv", config, "qc")

  # Create filtering plots
  plots <- create_filtering_plots(initial_n, filter_log, config)

  log_message("Filtering complete: ", initial_n, " -> ", nrow(mat), " features")
  log_message("=== Feature Filtering Complete ===")

  list(
    matrix = mat,
    metadata = metadata,
    sample_roles = sample_roles,
    filter_log = filter_log,
    plots = plots
  )
}

#' Calculate simple feature metrics
#'
#' @param mat Feature matrix
#' @param sample_roles Sample roles
#' @return Feature metrics data frame
calculate_feature_metrics_simple <- function(mat, sample_roles) {
  metrics <- data.frame(
    feature_id = rownames(mat),
    stringsAsFactors = FALSE
  )

  # QC-specific metrics
  if (length(sample_roles$qc_samples) > 0) {
    qc_mat <- mat[, sample_roles$qc_samples, drop = FALSE]
    metrics$qc_mean <- rowMeans(qc_mat, na.rm = TRUE)
    metrics$qc_sd <- apply(qc_mat, 1, sd, na.rm = TRUE)
    metrics$qc_rsd <- (metrics$qc_sd / metrics$qc_mean) * 100
  }

  return(metrics)
}

#' Filter by blank ratio
#'
#' Remove features where blank/sample ratio exceeds threshold
#'
#' @param mat Feature matrix
#' @param feature_metrics Feature metrics with blank_ratio column
#' @param config Configuration list
#' @return Filtered result list
filter_by_blank_ratio <- function(mat, feature_metrics, config) {
  threshold <- config$filtering$blank_ratio_threshold

  if (!"blank_ratio" %in% colnames(feature_metrics)) {
    return(list(matrix = mat, n_before = nrow(mat), n_removed = 0))
  }

  n_before <- nrow(mat)

  # Match features
  match_idx <- match(rownames(mat), feature_metrics$feature_id)
  blank_ratios <- feature_metrics$blank_ratio[match_idx]

  # Keep features with low blank ratio (or NA which means blank was zero)
  keep <- is.na(blank_ratios) | blank_ratios < threshold

  mat_filtered <- mat[keep, , drop = FALSE]
  n_removed <- n_before - nrow(mat_filtered)

  log_message("Blank ratio filter (<", threshold, "): removed ", n_removed, " features")

  list(
    matrix = mat_filtered,
    n_before = n_before,
    n_removed = n_removed
  )
}

#' Filter by QC RSD
#'
#' Remove features with high RSD in QC samples
#'
#' @param mat Feature matrix
#' @param feature_metrics Feature metrics with qc_rsd column
#' @param config Configuration list
#' @return Filtered result list
filter_by_qc_rsd <- function(mat, feature_metrics, config) {
  threshold <- config$filtering$qc_rsd_threshold * 100  # Convert to percentage

  if (!"qc_rsd" %in% colnames(feature_metrics)) {
    return(list(matrix = mat, n_before = nrow(mat), n_removed = 0))
  }

  n_before <- nrow(mat)

  # Match features
  match_idx <- match(rownames(mat), feature_metrics$feature_id)
  qc_rsds <- feature_metrics$qc_rsd[match_idx]

  # Keep features with low RSD (or NA)
  keep <- is.na(qc_rsds) | qc_rsds <= threshold

  mat_filtered <- mat[keep, , drop = FALSE]
  n_removed <- n_before - nrow(mat_filtered)

  log_message("QC RSD filter (<=", threshold, "%): removed ", n_removed, " features")

  list(
    matrix = mat_filtered,
    n_before = n_before,
    n_removed = n_removed
  )
}

#' Filter by global presence
#'
#' @param mat Feature matrix
#' @param config Configuration list
#' @return Filtered result list
filter_by_global_presence <- function(mat, config) {
  threshold <- config$filtering$global_min_presence

  if (threshold <= 0) {
    return(list(matrix = mat, n_before = nrow(mat), n_removed = 0))
  }

  n_before <- nrow(mat)
  n_samples <- ncol(mat)
  min_present <- ceiling(n_samples * threshold)

  n_present <- rowSums(!is.na(mat))
  keep <- n_present >= min_present

  mat_filtered <- mat[keep, , drop = FALSE]
  n_removed <- n_before - nrow(mat_filtered)

  log_message("Global presence filter (>=", round(threshold * 100), "%): removed ", n_removed, " features")

  list(
    matrix = mat_filtered,
    n_before = n_before,
    n_removed = n_removed
  )
}

#' Filter by group-specific presence
#'
#' @param mat Feature matrix
#' @param metadata Sample metadata
#' @param sample_roles Sample roles
#' @param config Configuration list
#' @return Filtered result list
filter_by_group_presence <- function(mat, metadata, sample_roles, config) {
  threshold <- config$filtering$group_min_presence

  if (threshold <= 0) {
    return(list(matrix = mat, n_before = nrow(mat), n_removed = 0))
  }

  n_before <- nrow(mat)

  condition_col <- config$design$condition_column
  sample_col <- config$input$sample_id_column

  # Only consider biological samples
  bio_samples <- intersect(sample_roles$biological_samples, colnames(mat))

  if (length(bio_samples) == 0) {
    log_message("No biological samples for group filtering. Skipping.")
    return(list(matrix = mat, n_before = n_before, n_removed = 0))
  }

  conditions <- unique(metadata[[condition_col]][metadata[[sample_col]] %in% bio_samples])

  # For each feature, check if it passes threshold in at least one group
  passes_filter <- rep(FALSE, nrow(mat))

  for (cond in conditions) {
    samples <- metadata[[sample_col]][metadata[[condition_col]] == cond]
    samples <- intersect(samples, bio_samples)
    samples <- intersect(samples, colnames(mat))

    if (length(samples) == 0) next

    sub_mat <- mat[, samples, drop = FALSE]
    n_samples <- ncol(sub_mat)
    min_present <- ceiling(n_samples * threshold)

    n_present <- rowSums(!is.na(sub_mat))
    passes_filter <- passes_filter | (n_present >= min_present)
  }

  mat_filtered <- mat[passes_filter, , drop = FALSE]
  n_removed <- n_before - nrow(mat_filtered)

  log_message("Group presence filter (>=", round(threshold * 100), "% in any group): removed ", n_removed, " features")

  list(
    matrix = mat_filtered,
    n_before = n_before,
    n_removed = n_removed
  )
}

#' Filter low variance features
#'
#' @param mat Feature matrix
#' @param config Configuration list
#' @return Filtered result list
filter_low_variance <- function(mat, config) {
  percentile <- config$filtering$variance_percentile %||% 0.05

  if (percentile <= 0) {
    return(list(matrix = mat, n_before = nrow(mat), n_removed = 0))
  }

  n_before <- nrow(mat)

  # Calculate variance per feature
  variances <- apply(mat, 1, var, na.rm = TRUE)

  # Find threshold
  var_threshold <- quantile(variances, probs = percentile, na.rm = TRUE)

  keep <- variances >= var_threshold
  mat_filtered <- mat[keep, , drop = FALSE]
  n_removed <- n_before - nrow(mat_filtered)

  log_message("Low variance filter (bottom ", round(percentile * 100), "%): removed ", n_removed, " features")

  list(
    matrix = mat_filtered,
    n_before = n_before,
    n_removed = n_removed
  )
}

#' Create filtering summary plots
#'
#' @param initial_n Initial number of features
#' @param filter_log Filter log data frame
#' @param config Configuration list
#' @return List of plots
create_filtering_plots <- function(initial_n, filter_log, config) {
  plots <- list()

  if (nrow(filter_log) == 0) {
    return(plots)
  }

  # Waterfall/funnel chart
  filter_log$step <- factor(filter_log$step, levels = filter_log$step)

  plots$filter_funnel <- ggplot2::ggplot(filter_log, ggplot2::aes(x = step, y = features_after)) +
    ggplot2::geom_bar(stat = "identity", fill = "#0072B2") +
    ggplot2::geom_text(ggplot2::aes(label = features_after), vjust = -0.3, size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = paste0("-", features_removed)),
                       vjust = 1.5, color = "red", size = 3) +
    ggplot2::labs(title = "Features Retained After Each Filter",
                  x = "Filtering Step", y = "Number of Features") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  save_plot(plots$filter_funnel, "filtering_funnel.png", config, subdir = "qc")

  return(plots)
}
