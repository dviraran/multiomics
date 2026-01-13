# =============================================================================
# Exploratory Multivariate Analysis
# =============================================================================

#' Run exploratory multivariate analysis
#'
#' @param imputed_data Data from imputation step
#' @param config Configuration list
#' @return Exploratory analysis results
run_exploratory_analysis <- function(imputed_data, config) {
  log_message("=== Starting Exploratory Analysis ===")

  mat <- imputed_data$matrix
  metadata <- imputed_data$metadata
  sample_roles <- imputed_data$sample_roles

  # Apply scaling for multivariate analysis
  scaling_method <- config$processing$scaling_method %||% "pareto"
  mat_scaled <- apply_scaling(mat, scaling_method)

  # PCA
  pca_results <- run_pca(mat_scaled, metadata, config)

  # UMAP (optional)
  umap_results <- NULL
  if (config$qc$run_umap %||% TRUE) {
    umap_results <- run_umap(mat_scaled, metadata, config)
  }

  # Sample correlations
  correlation_results <- calculate_sample_correlations(mat, metadata, config)

  # Hierarchical clustering
  clustering_results <- run_clustering(mat_scaled, metadata, config)

  # Outlier detection
  outlier_results <- detect_outliers(mat, pca_results, correlation_results, imputed_data, config)

  # PLS-DA (optional, exploratory only)
  plsda_results <- NULL
  if (config$multivariate$run_plsda %||% FALSE) {
    plsda_results <- run_plsda(mat_scaled, metadata, sample_roles, config)
  }

  # Save results
  save_table(pca_results$coordinates, "pca_scores.csv", config, "tables")
  save_table(as.data.frame(pca_results$loadings), "pca_loadings.csv", config, "tables")
  save_table(outlier_results$flagged_samples, "outlier_report.csv", config, "qc")

  log_message("=== Exploratory Analysis Complete ===")

  list(
    pca = pca_results,
    umap = umap_results,
    correlations = correlation_results,
    clustering = clustering_results,
    outliers = outlier_results,
    plsda = plsda_results,
    scaled_matrix = mat_scaled
  )
}

#' Run PCA analysis
#'
#' @param mat Scaled matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return PCA results
run_pca <- function(mat, metadata, config) {
  log_message("Running PCA...")

  n_components <- min(config$qc$n_pca_components %||% 10, ncol(mat) - 1, nrow(mat) - 1)

  # Remove features with zero variance
  var_per_feature <- apply(mat, 1, var, na.rm = TRUE)
  mat_filtered <- mat[var_per_feature > 0 & !is.na(var_per_feature), ]

  pca <- prcomp(t(mat_filtered), center = TRUE, scale. = FALSE)  # Already scaled

  # Extract results
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

  if (!is.null(config$design$batch_column) && config$design$batch_column %in% colnames(metadata)) {
    pca_coords$batch <- metadata[[config$design$batch_column]][match(pca_coords$sample_id, metadata[[sample_col]])]
  }

  # Variance explained
  var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
  var_explained <- var_explained[1:n_components]
  names(var_explained) <- paste0("PC", 1:n_components)

  # Loadings (top features per PC)
  loadings <- pca$rotation[, 1:n_components]

  # Create plots
  plots <- create_pca_plots(pca_coords, var_explained, config)

  list(
    coordinates = pca_coords,
    variance_explained = var_explained,
    loadings = loadings,
    pca_object = pca,
    plots = plots
  )
}

#' Create PCA plots
#'
#' @param pca_coords PCA coordinates
#' @param var_explained Variance explained
#' @param config Configuration
#' @return List of plots
create_pca_plots <- function(pca_coords, var_explained, config) {
  plots <- list()

  # PC1 vs PC2 by condition
  plots$pca_condition <- ggplot2::ggplot(pca_coords,
    ggplot2::aes(x = PC1, y = PC2, color = condition, shape = sample_type)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::labs(
      title = "PCA: PC1 vs PC2",
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      color = "Condition",
      shape = "Sample Type"
    ) +
    ggplot2::theme_minimal()

  # PC1 vs PC3
  if (length(var_explained) >= 3) {
    plots$pca_pc1_pc3 <- ggplot2::ggplot(pca_coords,
      ggplot2::aes(x = PC1, y = PC3, color = condition)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::labs(
        title = "PCA: PC1 vs PC3",
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC3 (", round(var_explained[3], 1), "%)")
      ) +
      ggplot2::theme_minimal()
  }

  # Scree plot
  scree_df <- data.frame(
    PC = factor(names(var_explained), levels = names(var_explained)),
    variance = var_explained,
    cumulative = cumsum(var_explained)
  )

  plots$scree <- ggplot2::ggplot(scree_df, ggplot2::aes(x = PC, y = variance)) +
    ggplot2::geom_bar(stat = "identity", fill = "#0072B2") +
    ggplot2::geom_line(ggplot2::aes(y = cumulative, group = 1), color = "red", size = 1) +
    ggplot2::geom_point(ggplot2::aes(y = cumulative), color = "red", size = 2) +
    ggplot2::labs(title = "PCA Scree Plot", x = "Principal Component", y = "% Variance") +
    ggplot2::theme_minimal()

  # Save plots
  save_plot(plots$pca_condition, "pca_condition.png", config, subdir = "plots")
  save_plot(plots$scree, "pca_scree.png", config, subdir = "qc")

  return(plots)
}

#' Run UMAP analysis
#'
#' @param mat Scaled matrix
#' @param metadata Sample metadata
#' @param config Configuration
#' @return UMAP results
run_umap <- function(mat, metadata, config) {
  if (!requireNamespace("umap", quietly = TRUE)) {
    log_message("umap package not available. Skipping UMAP.")
    return(NULL)
  }

  log_message("Running UMAP...")

  n_neighbors <- config$qc$umap_n_neighbors %||% 15
  min_dist <- config$qc$umap_min_dist %||% 0.1

  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- min(n_neighbors, ncol(mat) - 1)
  umap_config$min_dist <- min_dist

  tryCatch({
    umap_result <- umap::umap(t(mat), config = umap_config)

    umap_coords <- as.data.frame(umap_result$layout)
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    umap_coords$sample_id <- colnames(mat)

    sample_col <- config$input$sample_id_column
    condition_col <- config$design$condition_column
    umap_coords$condition <- metadata[[condition_col]][match(umap_coords$sample_id, metadata[[sample_col]])]

    # Create plot
    p <- ggplot2::ggplot(umap_coords, ggplot2::aes(x = UMAP1, y = UMAP2, color = condition)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::labs(title = "UMAP", color = "Condition") +
      ggplot2::theme_minimal()

    save_plot(p, "umap.png", config, subdir = "plots")

    list(
      coordinates = umap_coords,
      plot = p
    )
  }, error = function(e) {
    log_message("UMAP failed: ", e$message)
    NULL
  })
}

#' Calculate sample correlations
#'
#' @param mat Matrix
#' @param metadata Metadata
#' @param config Configuration
#' @return Correlation results
calculate_sample_correlations <- function(mat, metadata, config) {
  log_message("Calculating sample correlations...")

  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")

  # Mean correlation per sample
  mean_cors <- sapply(seq_len(ncol(cor_mat)), function(i) {
    mean(cor_mat[i, -i], na.rm = TRUE)
  })
  names(mean_cors) <- colnames(mat)

  # Create heatmap
  if (requireNamespace("pheatmap", quietly = TRUE) && ncol(mat) <= 100) {
    sample_col <- config$input$sample_id_column
    condition_col <- config$design$condition_column
    conditions <- metadata[[condition_col]][match(colnames(mat), metadata[[sample_col]])]

    annotation_df <- data.frame(Condition = conditions, row.names = colnames(mat))

    heatmap_file <- file.path(config$output$output_dir, "qc", "correlation_heatmap.png")
    png(heatmap_file, width = 10, height = 10, units = "in", res = 300)
    pheatmap::pheatmap(
      cor_mat,
      annotation_col = annotation_df,
      show_rownames = ncol(mat) <= 50,
      show_colnames = ncol(mat) <= 50,
      main = "Sample Correlation Heatmap"
    )
    dev.off()
    log_message("Saved: ", heatmap_file)
  }

  list(
    correlation_matrix = cor_mat,
    mean_correlation = mean_cors
  )
}

#' Run hierarchical clustering
#'
#' @param mat Scaled matrix
#' @param metadata Metadata
#' @param config Configuration
#' @return Clustering results
run_clustering <- function(mat, metadata, config) {
  log_message("Running hierarchical clustering...")

  dist_mat <- dist(t(mat), method = "euclidean")
  hclust_result <- hclust(dist_mat, method = "ward.D2")

  # Create dendrogram
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  conditions <- metadata[[condition_col]][match(colnames(mat), metadata[[sample_col]])]

  dend_file <- file.path(config$output$output_dir, "qc", "clustering_dendrogram.png")

  if (requireNamespace("dendextend", quietly = TRUE)) {
    dend <- as.dendrogram(hclust_result)

    n_conditions <- length(unique(conditions))
    colors <- RColorBrewer::brewer.pal(max(3, n_conditions), "Set1")[1:n_conditions]
    names(colors) <- unique(conditions)
    label_colors <- colors[as.character(conditions[order.dendrogram(dend)])]

    dend <- dendextend::set(dend, "labels_colors", label_colors)
    dend <- dendextend::set(dend, "labels_cex", 0.6)

    png(dend_file, width = 12, height = 8, units = "in", res = 300)
    plot(dend, main = "Sample Clustering Dendrogram")
    legend("topright", legend = names(colors), fill = colors, border = NA, bty = "n")
    dev.off()
  } else {
    png(dend_file, width = 12, height = 8, units = "in", res = 300)
    plot(hclust_result, main = "Sample Clustering", xlab = "", sub = "")
    dev.off()
  }

  log_message("Saved: ", dend_file)

  list(
    hclust = hclust_result,
    distance_matrix = dist_mat
  )
}

#' Detect potential outliers
#'
#' @param mat Matrix
#' @param pca_results PCA results
#' @param correlation_results Correlation results
#' @param imputed_data Imputed data
#' @param config Configuration
#' @return Outlier detection results
detect_outliers <- function(mat, pca_results, correlation_results, imputed_data, config) {
  log_message("Detecting potential outliers...")

  sd_threshold <- config$qc$outlier_sd_threshold %||% 3
  min_cor <- config$qc$min_sample_correlation %||% 0.8

  samples <- colnames(mat)
  flagged <- data.frame(
    sample_id = samples,
    pca_outlier = FALSE,
    correlation_outlier = FALSE,
    low_signal_outlier = FALSE,
    reason = "",
    stringsAsFactors = FALSE
  )

  # PCA-based detection (Mahalanobis distance)
  pca_coords <- pca_results$coordinates[, c("PC1", "PC2")]
  center <- colMeans(pca_coords)
  cov_mat <- cov(pca_coords)

  if (det(cov_mat) > 1e-10) {
    mahal_dist <- mahalanobis(pca_coords, center, cov_mat)
    mahal_threshold <- qchisq(0.99, df = 2)
    flagged$pca_outlier <- mahal_dist > mahal_threshold
  }

  # Correlation-based detection
  mean_cors <- correlation_results$mean_correlation
  flagged$correlation_outlier <- mean_cors < min_cor

  # Low total signal detection
  total_signal <- colSums(2^mat, na.rm = TRUE)
  signal_threshold <- quantile(total_signal, config$qc$min_total_signal_percentile %||% 0.05)
  flagged$low_signal_outlier <- total_signal < signal_threshold

  # Compile reasons
  for (i in seq_len(nrow(flagged))) {
    reasons <- c()
    if (flagged$pca_outlier[i]) reasons <- c(reasons, "Extreme in PCA")
    if (flagged$correlation_outlier[i]) reasons <- c(reasons, sprintf("Low correlation (%.3f)", mean_cors[i]))
    if (flagged$low_signal_outlier[i]) reasons <- c(reasons, "Low total signal")
    flagged$reason[i] <- paste(reasons, collapse = "; ")
  }

  flagged$is_outlier <- flagged$pca_outlier | flagged$correlation_outlier | flagged$low_signal_outlier

  n_outliers <- sum(flagged$is_outlier)
  log_message("Flagged ", n_outliers, " potential outlier samples (NOT auto-removed)")

  list(
    flagged_samples = flagged,
    n_outliers = n_outliers
  )
}

#' Run PLS-DA (exploratory, with cross-validation)
#'
#' @param mat Scaled matrix
#' @param metadata Metadata
#' @param sample_roles Sample roles
#' @param config Configuration
#' @return PLS-DA results
run_plsda <- function(mat, metadata, sample_roles, config) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    log_message("mixOmics package not available. Skipping PLS-DA.")
    return(NULL)
  }

  log_message("Running PLS-DA (exploratory)...")

  # Only use biological samples
  bio_samples <- intersect(sample_roles$biological_samples, colnames(mat))
  if (length(bio_samples) < 10) {
    log_message("Too few biological samples for PLS-DA")
    return(NULL)
  }

  mat_bio <- mat[, bio_samples]

  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  y <- metadata[[condition_col]][match(bio_samples, metadata[[sample_col]])]

  n_comp <- config$multivariate$plsda_ncomp %||% 3

  tryCatch({
    plsda_result <- mixOmics::plsda(t(mat_bio), y, ncomp = n_comp)

    # Cross-validation
    cv_folds <- config$multivariate$plsda_cv_folds %||% 5
    perf_result <- mixOmics::perf(plsda_result, validation = "Mfold",
                                   folds = cv_folds, nrepeat = 10)

    # Extract scores
    scores <- as.data.frame(plsda_result$variates$X)
    scores$sample_id <- rownames(scores)
    scores$condition <- y

    # Create plot
    p <- ggplot2::ggplot(scores, ggplot2::aes(x = comp1, y = comp2, color = condition)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(
        title = "PLS-DA (Exploratory)",
        subtitle = "Note: For hypothesis generation only, not inference",
        x = "Component 1",
        y = "Component 2"
      ) +
      ggplot2::theme_minimal()

    save_plot(p, "plsda_scores.png", config, subdir = "plots")

    list(
      model = plsda_result,
      performance = perf_result,
      scores = scores,
      plot = p
    )
  }, error = function(e) {
    log_message("PLS-DA failed: ", e$message)
    NULL
  })
}
