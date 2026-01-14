# =============================================================================
# Missing Value Imputation
# =============================================================================

#' Impute missing values
#'
#' @param batch_corrected_data Data from batch correction step
#' @param config Configuration list
#' @return Imputed data list
impute_missing_values <- function(batch_corrected_data, config) {
  log_message("=== Starting Missing Value Imputation ===")

  mat <- batch_corrected_data$matrix
  metadata <- batch_corrected_data$metadata
  sample_roles <- batch_corrected_data$sample_roles

  # Characterize missingness
  missingness <- characterize_missingness(mat, metadata, config)

  # Apply imputation
  method <- config$imputation$method %||% "half_min"
  mat_imputed <- apply_imputation(mat, method, config)

  # Create diagnostic plots
  plots <- create_imputation_plots(mat, mat_imputed, missingness, config)

  # Save imputed matrix
  imputed_df <- as.data.frame(mat_imputed)
  imputed_df <- tibble::rownames_to_column(imputed_df, "feature_id")
  save_table(imputed_df, "imputed_matrix.csv", config, "tables")

  # Save missingness summary
  save_table(missingness$summary, "missingness_summary.csv", config, "qc")
  save_table(missingness$per_feature, "missingness_per_feature.csv", config, "qc")

  log_message("=== Missing Value Imputation Complete ===")

  list(
    matrix = mat_imputed,
    matrix_pre_imputation = mat,
    metadata = metadata,
    sample_roles = sample_roles,
    missingness = missingness,
    imputation_method = method,
    plots = plots
  )
}

#' Characterize missingness pattern
#'
#' @param mat Matrix
#' @param metadata Metadata
#' @param config Configuration
#' @return Missingness characterization
characterize_missingness <- function(mat, metadata, config) {
  log_message("Characterizing missingness patterns...")

  n_features <- nrow(mat)
  n_samples <- ncol(mat)
  n_missing <- sum(is.na(mat))
  pct_missing <- n_missing / (n_features * n_samples) * 100

  # Per-feature missingness
  per_feature <- data.frame(
    feature_id = rownames(mat),
    n_present = rowSums(!is.na(mat)),
    n_missing = rowSums(is.na(mat)),
    pct_missing = rowMeans(is.na(mat)) * 100,
    mean_intensity = rowMeans(mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # Per-sample missingness
  per_sample <- data.frame(
    sample_id = colnames(mat),
    n_present = colSums(!is.na(mat)),
    n_missing = colSums(is.na(mat)),
    pct_missing = colMeans(is.na(mat)) * 100,
    stringsAsFactors = FALSE
  )

  # MNAR analysis
  mnar_analysis <- analyze_mnar_pattern(mat)

  # Summary
  summary_df <- data.frame(
    metric = c(
      "Total values",
      "Missing values",
      "% Missing overall",
      "Features with any missing",
      "Features with >50% missing",
      "Samples with >30% missing",
      "MNAR correlation",
      "Missingness pattern"
    ),
    value = c(
      n_features * n_samples,
      n_missing,
      round(pct_missing, 2),
      sum(per_feature$n_missing > 0),
      sum(per_feature$pct_missing > 50),
      sum(per_sample$pct_missing > 30),
      round(mnar_analysis$correlation, 3),
      mnar_analysis$pattern
    ),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_df,
    per_feature = per_feature,
    per_sample = per_sample,
    mnar_analysis = mnar_analysis
  )
}

#' Analyze MNAR (Missing Not At Random) pattern
#'
#' @param mat Matrix
#' @return MNAR analysis results
analyze_mnar_pattern <- function(mat) {
  # Correlation between mean intensity and missingness
  mean_intensity <- rowMeans(mat, na.rm = TRUE)
  pct_missing <- rowMeans(is.na(mat))

  valid_idx <- is.finite(mean_intensity) & is.finite(pct_missing)
  correlation <- cor(mean_intensity[valid_idx], pct_missing[valid_idx], method = "spearman")

  # Determine pattern
  if (correlation < -0.3) {
    pattern <- "MNAR (left-censored)"
  } else if (correlation > 0.3) {
    pattern <- "Unusual pattern"
  } else {
    pattern <- "MAR/MCAR (random)"
  }

  log_message("Missingness pattern: ", pattern, " (rho = ", round(correlation, 3), ")")

  list(
    correlation = correlation,
    pattern = pattern
  )
}

#' Apply imputation method
#'
#' @param mat Matrix
#' @param method Imputation method
#' @param config Configuration
#' @return Imputed matrix
apply_imputation <- function(mat, method, config) {
  log_message("Applying imputation method: ", method)

  if (method == "none") {
    log_message("Skipping imputation")
    return(mat)
  }

  imputed <- switch(method,
    "half_min" = impute_half_min(mat, config),
    "minprob" = impute_minprob(mat, config),
    "QRILC" = impute_qrilc(mat, config),
    "knn" = impute_knn(mat, config),
    {
      warning("Unknown imputation method '", method, "'. Using half_min.")
      impute_half_min(mat, config)
    }
  )

  # Check for remaining NAs
  remaining_na <- sum(is.na(imputed))
  if (remaining_na > 0) {
    log_message("Warning: ", remaining_na, " missing values remain after imputation")
  } else {
    log_message("Imputation complete. No missing values remain.")
  }

  return(imputed)
}

#' Half-minimum imputation
#'
#' Replace missing values with half the minimum observed value per feature
#'
#' @param mat Matrix
#' @param config Configuration
#' @return Imputed matrix
impute_half_min <- function(mat, config) {
  fraction <- config$imputation$half_min_fraction %||% 0.5

  mat_imputed <- mat

  for (i in seq_len(nrow(mat))) {
    row <- mat[i, ]
    na_idx <- is.na(row)

    if (any(na_idx) && any(!na_idx)) {
      min_val <- min(row, na.rm = TRUE) * fraction
      mat_imputed[i, na_idx] <- min_val
    } else if (all(na_idx)) {
      # Use global minimum if entire row is NA
      global_min <- min(mat, na.rm = TRUE) * fraction
      mat_imputed[i, ] <- global_min
    }
  }

  log_message("Half-min imputation complete (fraction = ", fraction, ")")
  return(mat_imputed)
}

#' MinProb imputation
#'
#' Sample from low-probability distribution
#'
#' @param mat Matrix
#' @param config Configuration
#' @return Imputed matrix
impute_minprob <- function(mat, config) {
  q <- config$imputation$min_prob_quantile %||% 0.01
  tune_sigma <- 0.3

  mat_imputed <- mat

  for (j in seq_len(ncol(mat))) {
    col_data <- mat[, j]
    na_idx <- is.na(col_data)

    if (sum(na_idx) == 0) next

    if (sum(!na_idx) >= 3) {
      col_min <- quantile(col_data, probs = q, na.rm = TRUE)
      col_sd <- sd(col_data, na.rm = TRUE) * tune_sigma
    } else {
      col_min <- quantile(mat, probs = q, na.rm = TRUE)
      col_sd <- sd(mat, na.rm = TRUE) * tune_sigma
    }

    imputed_vals <- rnorm(sum(na_idx), mean = col_min, sd = col_sd)
    mat_imputed[na_idx, j] <- imputed_vals
  }

  log_message("MinProb imputation complete")
  return(mat_imputed)
}

#' QRILC imputation
#'
#' Quantile Regression Imputation of Left-Censored data
#'
#' @param mat Matrix
#' @param config Configuration
#' @return Imputed matrix
impute_qrilc <- function(mat, config) {
  if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
    log_message("imputeLCMD package not available. Using half_min instead.")
    return(impute_half_min(mat, config))
  }

  tryCatch({
    # imputeLCMD expects samples x features
    mat_t <- t(mat)

    result <- imputeLCMD::impute.QRILC(mat_t, tune.sigma = 1)
    mat_imputed <- t(result[[1]])

    log_message("QRILC imputation complete")
    return(mat_imputed)
  }, error = function(e) {
    log_message("QRILC failed: ", e$message, ". Using half_min instead.")
    return(impute_half_min(mat, config))
  })
}

#' kNN imputation
#'
#' @param mat Matrix
#' @param config Configuration
#' @return Imputed matrix
impute_knn <- function(mat, config) {
  k <- config$imputation$knn_k %||% 10

  if (!requireNamespace("impute", quietly = TRUE)) {
    log_message("impute package not available. Using half_min instead.")
    return(impute_half_min(mat, config))
  }

  tryCatch({
    result <- impute::impute.knn(mat, k = k, rowmax = 0.8, colmax = 0.8)
    mat_imputed <- result$data

    log_message("kNN imputation complete (k = ", k, ")")
    return(mat_imputed)
  }, error = function(e) {
    log_message("kNN imputation failed: ", e$message, ". Using half_min instead.")
    return(impute_half_min(mat, config))
  })
}

#' Create imputation diagnostic plots
#'
#' @param mat_before Matrix before imputation
#' @param mat_after Matrix after imputation
#' @param missingness Missingness analysis
#' @param config Configuration
#' @return List of plots
create_imputation_plots <- function(mat_before, mat_after, missingness, config) {
  plots <- list()

  # 1. MNAR diagnostic (missingness vs intensity)
  mnar_df <- data.frame(
    mean_intensity = rowMeans(mat_before, na.rm = TRUE),
    pct_missing = rowMeans(is.na(mat_before)) * 100
  )
  mnar_df <- mnar_df[is.finite(mnar_df$mean_intensity), ]

  plots$mnar_diagnostic <- ggplot2::ggplot(mnar_df,
    ggplot2::aes(x = mean_intensity, y = pct_missing)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.5) +
    ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE) +
    ggplot2::labs(
      title = "Missingness vs Intensity (MNAR Diagnostic)",
      subtitle = paste("Pattern:", missingness$mnar_analysis$pattern),
      x = "Mean Log2 Intensity",
      y = "% Missing"
    ) +
    ggplot2::theme_minimal()

  # 2. Before/after density
  density_before <- data.frame(
    intensity = as.vector(mat_before),
    type = "Observed"
  )
  density_before <- density_before[!is.na(density_before$intensity), ]

  density_after <- data.frame(
    intensity = as.vector(mat_after),
    type = "After Imputation"
  )

  density_df <- rbind(density_before, density_after)

  plots$imputation_density <- ggplot2::ggplot(density_df,
    ggplot2::aes(x = intensity, fill = type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(
      title = "Intensity Distribution Before/After Imputation",
      x = "Log2 Intensity",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = c("Observed" = "#56B4E9", "After Imputation" = "#E69F00"))

  # 3. Feature missingness histogram
  plots$feature_missing <- ggplot2::ggplot(missingness$per_feature,
    ggplot2::aes(x = pct_missing)) +
    ggplot2::geom_histogram(bins = 50, fill = "#0072B2", color = "white") +
    ggplot2::labs(
      title = "Feature Missingness Distribution",
      x = "% Missing per Feature",
      y = "Count"
    ) +
    ggplot2::theme_minimal()

  # Save plots
  save_plot(plots$mnar_diagnostic, "imputation_mnar_diagnostic.png", config, subdir = "qc")
  save_plot(plots$imputation_density, "imputation_density.png", config, subdir = "qc")
  save_plot(plots$feature_missing, "imputation_feature_missing.png", config, subdir = "qc")

  return(plots)
}
