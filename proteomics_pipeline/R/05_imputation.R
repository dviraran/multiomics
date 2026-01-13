# =============================================================================
# Missing Value Characterization and Imputation
# =============================================================================

#' Characterize and impute missing values
#'
#' @param normalized_data Data from normalization step
#' @param config Configuration list
#' @return Imputed data list
impute_missing_values <- function(normalized_data, config) {
  log_message("=== Starting Missing Value Handling ===")

  mat <- normalized_data$matrix
  metadata <- normalized_data$metadata

  # Step 1: Characterize missingness
  missingness_analysis <- characterize_missingness(mat, metadata, config)

  # Step 2: Apply imputation
  imputation_method <- config$imputation$method
  mat_imputed <- apply_imputation(mat, imputation_method, config)

  # Create diagnostic plots
  plots <- create_imputation_plots(mat, mat_imputed, missingness_analysis, config)

  # Save imputed matrix
  imputed_df <- as.data.frame(mat_imputed)
  imputed_df <- tibble::rownames_to_column(imputed_df, "feature_id")
  save_table(imputed_df, "imputed_matrix.csv", config, "tables")

  # Save missingness analysis
  save_table(missingness_analysis$summary, "missingness_summary.csv", config, "qc")
  save_table(missingness_analysis$per_feature, "missingness_per_feature.csv", config, "qc")

  log_message("=== Missing Value Handling Complete ===")

  list(
    matrix = mat_imputed,
    matrix_pre_imputation = mat,
    metadata = metadata,
    annotations = normalized_data$annotations,
    missingness_analysis = missingness_analysis,
    imputation_method = imputation_method,
    plots = plots
  )
}

#' Characterize the missingness pattern
#'
#' Determines if missingness is MCAR (Missing Completely At Random) or
#' MNAR (Missing Not At Random / left-censored)
#'
#' @param mat Normalized matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Missingness characterization list
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
    sample = colnames(mat),
    n_present = colSums(!is.na(mat)),
    n_missing = colSums(is.na(mat)),
    pct_missing = colMeans(is.na(mat)) * 100,
    stringsAsFactors = FALSE
  )

  # MNAR vs MCAR analysis
  # If missing values are more common in low-abundance features, suggests MNAR
  mnar_evidence <- analyze_mnar(mat)

  # Summary
  summary_df <- data.frame(
    metric = c(
      "Total values",
      "Missing values",
      "% Missing overall",
      "Features with any missing",
      "Features with >50% missing",
      "Samples with >20% missing",
      "MNAR evidence (correlation)",
      "Missingness pattern"
    ),
    value = c(
      n_features * n_samples,
      n_missing,
      round(pct_missing, 2),
      sum(per_feature$n_missing > 0),
      sum(per_feature$pct_missing > 50),
      sum(per_sample$pct_missing > 20),
      round(mnar_evidence$correlation, 3),
      mnar_evidence$pattern
    ),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_df,
    per_feature = per_feature,
    per_sample = per_sample,
    mnar_evidence = mnar_evidence
  )
}

#' Analyze evidence for MNAR (left-censored) missingness
#'
#' @param mat Normalized matrix
#' @return List with correlation and pattern assessment
analyze_mnar <- function(mat) {
  # Calculate mean intensity per feature
  mean_intensity <- rowMeans(mat, na.rm = TRUE)

  # Calculate missingness per feature
  pct_missing <- rowMeans(is.na(mat))

  # Correlation between mean intensity and missingness
  # Negative correlation suggests MNAR (low abundance = more missing)
  valid_idx <- is.finite(mean_intensity) & is.finite(pct_missing)
  correlation <- cor(mean_intensity[valid_idx], pct_missing[valid_idx],
                     method = "spearman")

  # Determine pattern
  if (correlation < -0.3) {
    pattern <- "MNAR (left-censored)"
  } else if (correlation > 0.3) {
    pattern <- "Unusual (high abundance = more missing)"
  } else {
    pattern <- "MCAR-like (random)"
  }

  log_message("Missingness pattern: ", pattern, " (rho = ", round(correlation, 3), ")")

  list(
    correlation = correlation,
    pattern = pattern
  )
}

#' Apply selected imputation method
#'
#' @param mat Normalized matrix
#' @param method Imputation method
#' @param config Configuration list
#' @return Imputed matrix
apply_imputation <- function(mat, method, config) {
  log_message("Applying imputation method: ", method)

  if (method == "none") {
    log_message("Skipping imputation (method = none)")
    return(mat)
  }

  imputed <- switch(method,
    "QRILC" = impute_qrilc(mat, config),
    "MinProb" = impute_minprob(mat, config),
    "MinDet" = impute_mindet(mat, config),
    "knn" = impute_knn(mat, config),
    "random_forest" = impute_rf(mat, config),
    {
      warning("Unknown imputation method '", method, "'. Using QRILC.")
      impute_qrilc(mat, config)
    }
  )

  # Verify no NAs remain
  remaining_na <- sum(is.na(imputed))
  if (remaining_na > 0) {
    log_message("Warning: ", remaining_na, " missing values remain after imputation")
  } else {
    log_message("Imputation complete. No missing values remain.")
  }

  return(imputed)
}

#' QRILC imputation (Quantile Regression Imputation of Left-Censored data)
#'
#' Best for MNAR data typical in proteomics
#'
#' @param mat Normalized matrix
#' @param config Configuration list
#' @return Imputed matrix
impute_qrilc <- function(mat, config) {
  if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
    warning("imputeLCMD package not available. Falling back to MinProb.")
    return(impute_minprob(mat, config))
  }

  # QRILC requires complete columns for fitting
  # Transpose for imputeLCMD (expects samples x features)
  mat_t <- t(mat)

  tryCatch({
    # impute.QRILC returns list with imputed matrix
    result <- imputeLCMD::impute.QRILC(mat_t, tune.sigma = 1)
    mat_imputed <- t(result[[1]])

    log_message("QRILC imputation successful")
    return(mat_imputed)
  }, error = function(e) {
    log_message("QRILC failed: ", e$message, ". Falling back to MinProb.")
    return(impute_minprob(mat, config))
  })
}

#' MinProb imputation
#'
#' Replaces missing values with values drawn from a low-probability distribution
#'
#' @param mat Normalized matrix
#' @param config Configuration list
#' @return Imputed matrix
impute_minprob <- function(mat, config) {
  q <- config$imputation$min_prob_quantile %||% 0.01
  tune_sigma <- config$imputation$tune_sigma %||% 1

  mat_imputed <- mat

  for (i in seq_len(ncol(mat))) {
    col_data <- mat[, i]
    na_idx <- is.na(col_data)

    if (sum(na_idx) == 0) next
    if (sum(!na_idx) < 3) {
      # Not enough data to estimate, use global parameters
      global_min <- quantile(mat, probs = q, na.rm = TRUE)
      global_sd <- sd(mat, na.rm = TRUE) * tune_sigma
      mat_imputed[na_idx, i] <- rnorm(sum(na_idx), mean = global_min, sd = global_sd * 0.3)
      next
    }

    # Estimate parameters from observed data
    col_min <- quantile(col_data, probs = q, na.rm = TRUE)
    col_sd <- sd(col_data, na.rm = TRUE) * tune_sigma

    # Draw imputed values
    imputed_vals <- rnorm(sum(na_idx), mean = col_min, sd = col_sd * 0.3)
    mat_imputed[na_idx, i] <- imputed_vals
  }

  log_message("MinProb imputation complete (q = ", q, ")")
  return(mat_imputed)
}

#' MinDet imputation
#'
#' Replaces missing values with a deterministic minimum
#'
#' @param mat Normalized matrix
#' @param config Configuration list
#' @return Imputed matrix
impute_mindet <- function(mat, config) {
  q <- config$imputation$min_prob_quantile %||% 0.01

  mat_imputed <- mat

  for (i in seq_len(ncol(mat))) {
    col_data <- mat[, i]
    na_idx <- is.na(col_data)

    if (sum(na_idx) == 0) next

    # Use column-specific minimum if available, else global
    if (sum(!na_idx) >= 3) {
      min_val <- quantile(col_data, probs = q, na.rm = TRUE)
    } else {
      min_val <- quantile(mat, probs = q, na.rm = TRUE)
    }

    mat_imputed[na_idx, i] <- min_val
  }

  log_message("MinDet imputation complete")
  return(mat_imputed)
}

#' KNN imputation
#'
#' Uses k-nearest neighbors for imputation
#'
#' @param mat Normalized matrix
#' @param config Configuration list
#' @return Imputed matrix
impute_knn <- function(mat, config) {
  if (!requireNamespace("impute", quietly = TRUE)) {
    warning("impute package not available. Falling back to MinProb.")
    return(impute_minprob(mat, config))
  }

  # impute.knn expects features in rows, samples in columns
  tryCatch({
    result <- impute::impute.knn(mat, k = 10, rowmax = 0.8, colmax = 0.8)
    mat_imputed <- result$data

    log_message("KNN imputation successful")
    return(mat_imputed)
  }, error = function(e) {
    log_message("KNN failed: ", e$message, ". Falling back to MinProb.")
    return(impute_minprob(mat, config))
  })
}

#' Random Forest imputation
#'
#' Uses missForest for imputation (can be slow for large datasets)
#'
#' @param mat Normalized matrix
#' @param config Configuration list
#' @return Imputed matrix
impute_rf <- function(mat, config) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    warning("missForest package not available. Falling back to MinProb.")
    return(impute_minprob(mat, config))
  }

  # missForest expects samples in rows
  mat_t <- t(mat)

  # Limit for performance
  if (nrow(mat_t) > 100) {
    log_message("Warning: Random forest imputation may be slow for large datasets")
  }

  tryCatch({
    result <- missForest::missForest(mat_t, maxiter = 10, ntree = 100, verbose = FALSE)
    mat_imputed <- t(result$ximp)

    log_message("Random forest imputation successful (OOB error: ", round(result$OOBerror, 4), ")")
    return(mat_imputed)
  }, error = function(e) {
    log_message("Random forest failed: ", e$message, ". Falling back to MinProb.")
    return(impute_minprob(mat, config))
  })
}

#' Create imputation diagnostic plots
#'
#' @param mat_before Matrix before imputation
#' @param mat_after Matrix after imputation
#' @param missingness_analysis Missingness analysis results
#' @param config Configuration list
#' @return List of ggplot objects
create_imputation_plots <- function(mat_before, mat_after, missingness_analysis, config) {
  plots <- list()

  # 1. Missingness vs intensity (MNAR diagnostic)
  mnar_df <- data.frame(
    mean_intensity = rowMeans(mat_before, na.rm = TRUE),
    pct_missing = rowMeans(is.na(mat_before)) * 100
  )
  mnar_df <- mnar_df[is.finite(mnar_df$mean_intensity), ]

  plots$mnar_diagnostic <- ggplot2::ggplot(mnar_df, ggplot2::aes(x = mean_intensity, y = pct_missing)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.5) +
    ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE) +
    ggplot2::labs(
      title = "Missingness vs Intensity (MNAR Diagnostic)",
      subtitle = paste("Correlation:", round(missingness_analysis$mnar_evidence$correlation, 3)),
      x = "Mean Log2 Intensity",
      y = "% Missing"
    ) +
    ggplot2::theme_minimal()

  # 2. Before/after imputation density
  # Sample for visualization if large
  n_samples <- min(10, ncol(mat_before))
  sample_idx <- sample(seq_len(ncol(mat_before)), n_samples)

  density_before <- data.frame(
    intensity = as.vector(mat_before[, sample_idx]),
    type = "Observed"
  )
  density_after <- data.frame(
    intensity = as.vector(mat_after[, sample_idx]),
    type = "After Imputation"
  )
  density_df <- rbind(
    density_before[!is.na(density_before$intensity), ],
    density_after
  )

  plots$imputation_density <- ggplot2::ggplot(density_df, ggplot2::aes(x = intensity, fill = type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(
      title = "Intensity Distribution Before/After Imputation",
      x = "Log2 Intensity",
      y = "Density",
      fill = "Data"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = c("Observed" = "#56B4E9", "After Imputation" = "#E69F00"))

  # 3. Per-sample missingness barplot
  sample_miss <- missingness_analysis$per_sample
  sample_miss <- sample_miss[order(sample_miss$pct_missing, decreasing = TRUE), ]

  if (nrow(sample_miss) > 50) {
    sample_miss <- sample_miss[1:50, ]  # Top 50 for visibility
  }

  sample_miss$sample <- factor(sample_miss$sample, levels = sample_miss$sample)

  plots$sample_missingness <- ggplot2::ggplot(sample_miss, ggplot2::aes(x = sample, y = pct_missing)) +
    ggplot2::geom_bar(stat = "identity", fill = "#0072B2") +
    ggplot2::labs(
      title = "Missingness by Sample",
      x = "Sample",
      y = "% Missing"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6))

  # Save plots
  save_plot(plots$mnar_diagnostic, "mnar_diagnostic.png", config, subdir = "qc")
  save_plot(plots$imputation_density, "imputation_density.png", config, subdir = "qc")
  save_plot(plots$sample_missingness, "sample_missingness.png", config, subdir = "qc")

  return(plots)
}
