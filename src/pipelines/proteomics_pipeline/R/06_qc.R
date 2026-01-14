# =============================================================================
# Quality Control and Exploratory Analysis
# =============================================================================

#' Perform QC and exploratory analysis
#'
#' @param imputed_data Data from imputation step
#' @param config Configuration list
#' @return QC results list
run_qc_analysis <- function(imputed_data, config) {
  log_message("=== Starting QC Analysis ===")

  mat <- imputed_data$matrix
  metadata <- imputed_data$metadata

  # Run PCA
  pca_results <- run_pca(mat, metadata, config)

  # Run UMAP (optional)
  umap_results <- NULL
  if (config$qc$run_umap %||% TRUE) {
    umap_results <- run_umap(mat, metadata, config)
  }

  # Sample correlation analysis
  correlation_results <- analyze_sample_correlation(mat, metadata, config)

  # Clustering analysis
  clustering_results <- run_clustering(mat, metadata, config)

  # Outlier detection
  outlier_results <- detect_outliers(mat, pca_results, correlation_results, config)

  # Batch effect analysis (if batch variable present)
  batch_results <- NULL
  if (!is.null(config$design$batch_column)) {
    batch_results <- analyze_batch_effects(mat, pca_results, metadata, config)
  }

  # Save outlier report
  save_table(outlier_results$flagged_samples, "outlier_report.csv", config, "qc")

  log_message("=== QC Analysis Complete ===")

  list(
    pca = pca_results,
    umap = umap_results,
    correlation = correlation_results,
    clustering = clustering_results,
    outliers = outlier_results,
    batch_effects = batch_results
  )
}

#' Run PCA analysis
#'
#' @param mat Imputed matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return PCA results list
run_pca <- function(mat, metadata, config) {
  log_message("Running PCA...")

  n_components <- min(config$qc$n_pca_components %||% 10, ncol(mat) - 1, nrow(mat) - 1)

  # Transpose for prcomp (samples in rows)
  mat_t <- t(mat)

  # Remove features with zero variance
  var_per_feature <- apply(mat, 1, var, na.rm = TRUE)
  mat_filtered <- mat[var_per_feature > 0, ]

  # Run PCA
  pca <- prcomp(t(mat_filtered), center = TRUE, scale. = TRUE)

  # Extract results
  pca_coords <- as.data.frame(pca$x[, 1:n_components])
  pca_coords$sample <- rownames(pca_coords)

  # Add metadata
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  pca_coords$condition <- metadata[[condition_col]][match(pca_coords$sample, metadata[[sample_col]])]

  if (!is.null(config$design$batch_column)) {
    pca_coords$batch <- metadata[[config$design$batch_column]][match(pca_coords$sample, metadata[[sample_col]])]
  }

  # Variance explained
  var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
  var_explained <- var_explained[1:n_components]
  names(var_explained) <- paste0("PC", 1:n_components)

  # Create PCA plots
  plots <- create_pca_plots(pca_coords, var_explained, config)

  # Save PCA coordinates
  save_table(pca_coords, "pca_coordinates.csv", config, "qc")

  list(
    coordinates = pca_coords,
    variance_explained = var_explained,
    loadings = pca$rotation[, 1:n_components],
    pca_object = pca,
    plots = plots
  )
}

#' Create PCA plots
#'
#' @param pca_coords PCA coordinates data frame
#' @param var_explained Variance explained per PC
#' @param config Configuration list
#' @return List of ggplot objects
create_pca_plots <- function(pca_coords, var_explained, config) {
  plots <- list()

  # PC1 vs PC2 colored by condition
  plots$pca_1_2 <- ggplot2::ggplot(pca_coords, ggplot2::aes(x = PC1, y = PC2, color = condition)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::labs(
      title = "PCA: PC1 vs PC2",
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      color = "Condition"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")

  # PC1 vs PC3
  if (length(var_explained) >= 3) {
    plots$pca_1_3 <- ggplot2::ggplot(pca_coords, ggplot2::aes(x = PC1, y = PC3, color = condition)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::labs(
        title = "PCA: PC1 vs PC3",
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC3 (", round(var_explained[3], 1), "%)"),
        color = "Condition"
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
    ggplot2::labs(
      title = "PCA Scree Plot",
      x = "Principal Component",
      y = "% Variance Explained"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # PCA with batch coloring if available
  if ("batch" %in% colnames(pca_coords)) {
    plots$pca_batch <- ggplot2::ggplot(pca_coords, ggplot2::aes(x = PC1, y = PC2, color = batch, shape = condition)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::labs(
        title = "PCA: Colored by Batch",
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)")
      ) +
      ggplot2::theme_minimal()
  }

  # Save plots
  save_plot(plots$pca_1_2, "pca_pc1_pc2.png", config, subdir = "qc")
  save_plot(plots$scree, "pca_scree.png", config, subdir = "qc")

  if (!is.null(plots$pca_1_3)) {
    save_plot(plots$pca_1_3, "pca_pc1_pc3.png", config, subdir = "qc")
  }
  if (!is.null(plots$pca_batch)) {
    save_plot(plots$pca_batch, "pca_batch.png", config, subdir = "qc")
  }

  return(plots)
}

#' Run UMAP analysis
#'
#' @param mat Imputed matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return UMAP results list
run_umap <- function(mat, metadata, config) {
  if (!requireNamespace("umap", quietly = TRUE)) {
    log_message("umap package not available. Skipping UMAP.")
    return(NULL)
  }

  log_message("Running UMAP...")

  n_neighbors <- config$qc$umap_n_neighbors %||% 15
  min_dist <- config$qc$umap_min_dist %||% 0.1

  # UMAP configuration
  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- min(n_neighbors, ncol(mat) - 1)
  umap_config$min_dist <- min_dist

  tryCatch({
    # Transpose for umap (samples in rows)
    umap_result <- umap::umap(t(mat), config = umap_config)

    umap_coords <- as.data.frame(umap_result$layout)
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    umap_coords$sample <- colnames(mat)

    # Add metadata
    sample_col <- config$input$sample_id_column
    condition_col <- config$design$condition_column
    umap_coords$condition <- metadata[[condition_col]][match(umap_coords$sample, metadata[[sample_col]])]

    # Create plot
    p <- ggplot2::ggplot(umap_coords, ggplot2::aes(x = UMAP1, y = UMAP2, color = condition)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::labs(
        title = "UMAP",
        x = "UMAP 1",
        y = "UMAP 2",
        color = "Condition"
      ) +
      ggplot2::theme_minimal()

    save_plot(p, "umap.png", config, subdir = "qc")
    save_table(umap_coords, "umap_coordinates.csv", config, "qc")

    list(
      coordinates = umap_coords,
      plot = p
    )
  }, error = function(e) {
    log_message("UMAP failed: ", e$message)
    NULL
  })
}

#' Analyze sample correlations
#'
#' @param mat Imputed matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Correlation analysis results
analyze_sample_correlation <- function(mat, metadata, config) {
  log_message("Analyzing sample correlations...")

  # Calculate correlation matrix
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")

  # Summary statistics
  cor_values <- cor_mat[lower.tri(cor_mat)]
  cor_summary <- data.frame(
    metric = c("Min correlation", "Max correlation", "Mean correlation", "Median correlation"),
    value = c(min(cor_values), max(cor_values), mean(cor_values), median(cor_values)),
    stringsAsFactors = FALSE
  )

  # Create heatmap
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  sample_conditions <- metadata[[condition_col]][match(colnames(mat), metadata[[sample_col]])]

  # Create annotation for heatmap
  annotation_df <- data.frame(
    Condition = sample_conditions,
    row.names = colnames(mat)
  )

  # Generate heatmap using pheatmap if available
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    heatmap_file <- file.path(config$output$output_dir, "qc", "correlation_heatmap.png")
    png(heatmap_file, width = 10, height = 10, units = "in", res = 300)
    pheatmap::pheatmap(
      cor_mat,
      annotation_col = annotation_df,
      annotation_row = annotation_df,
      show_rownames = ncol(mat) <= 50,
      show_colnames = ncol(mat) <= 50,
      main = "Sample Correlation Heatmap"
    )
    dev.off()
    log_message("Saved: ", heatmap_file)
  }

  # Per-sample mean correlation (excluding self)
  sample_mean_cor <- sapply(seq_len(ncol(cor_mat)), function(i) {
    mean(cor_mat[i, -i])
  })
  names(sample_mean_cor) <- colnames(mat)

  per_sample_cor <- data.frame(
    sample = colnames(mat),
    mean_correlation = sample_mean_cor,
    condition = sample_conditions,
    stringsAsFactors = FALSE
  )

  save_table(per_sample_cor, "sample_correlations.csv", config, "qc")

  list(
    correlation_matrix = cor_mat,
    summary = cor_summary,
    per_sample = per_sample_cor
  )
}

#' Run hierarchical clustering
#'
#' @param mat Imputed matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Clustering results
run_clustering <- function(mat, metadata, config) {
  log_message("Running hierarchical clustering...")

  # Calculate distance matrix
  dist_mat <- dist(t(mat), method = "euclidean")

  # Hierarchical clustering
  hclust_result <- hclust(dist_mat, method = "ward.D2")

  # Create dendrogram with color coding
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  conditions <- metadata[[condition_col]][match(colnames(mat), metadata[[sample_col]])]

  # Create dendrogram plot
  if (requireNamespace("dendextend", quietly = TRUE)) {
    dend <- as.dendrogram(hclust_result)

    # Get colors for conditions
    n_conditions <- length(unique(conditions))
    condition_colors <- RColorBrewer::brewer.pal(max(3, n_conditions), "Set1")[1:n_conditions]
    names(condition_colors) <- unique(conditions)
    label_colors <- condition_colors[as.character(conditions[order.dendrogram(dend)])]

    dend <- dendextend::set(dend, "labels_colors", label_colors)
    dend <- dendextend::set(dend, "labels_cex", 0.6)

    dend_file <- file.path(config$output$output_dir, "qc", "clustering_dendrogram.png")
    png(dend_file, width = 12, height = 8, units = "in", res = 300)
    plot(dend, main = "Sample Clustering Dendrogram")
    legend("topright", legend = names(condition_colors), fill = condition_colors, border = NA, bty = "n")
    dev.off()
    log_message("Saved: ", dend_file)
  } else {
    # Basic dendrogram
    dend_file <- file.path(config$output$output_dir, "qc", "clustering_dendrogram.png")
    png(dend_file, width = 12, height = 8, units = "in", res = 300)
    plot(hclust_result, main = "Sample Clustering Dendrogram", xlab = "", sub = "")
    dev.off()
    log_message("Saved: ", dend_file)
  }

  list(
    hclust = hclust_result,
    distance_matrix = dist_mat
  )
}

#' Detect potential outlier samples
#'
#' @param mat Imputed matrix
#' @param pca_results PCA results
#' @param correlation_results Correlation results
#' @param config Configuration list
#' @return Outlier detection results
detect_outliers <- function(mat, pca_results, correlation_results, config) {
  log_message("Detecting potential outliers...")

  sd_threshold <- config$qc$outlier_sd_threshold %||% 3
  min_cor_threshold <- config$qc$min_sample_correlation %||% 0.8

  samples <- colnames(mat)
  flagged <- data.frame(
    sample = samples,
    pca_outlier = FALSE,
    correlation_outlier = FALSE,
    reason = "",
    stringsAsFactors = FALSE
  )

  # PCA-based outlier detection (Mahalanobis distance on PC1-PC2)
  pca_coords <- pca_results$coordinates[, c("PC1", "PC2")]
  center <- colMeans(pca_coords)
  cov_mat <- cov(pca_coords)

  if (det(cov_mat) > 1e-10) {
    mahal_dist <- mahalanobis(pca_coords, center, cov_mat)
    mahal_threshold <- qchisq(0.99, df = 2)  # 99% confidence for 2 PCs
    pca_outliers <- mahal_dist > mahal_threshold

    flagged$pca_outlier <- pca_outliers
  }

  # Correlation-based outlier detection
  mean_cors <- correlation_results$per_sample$mean_correlation
  flagged$correlation_outlier <- mean_cors < min_cor_threshold

  # Compile reasons
  for (i in seq_len(nrow(flagged))) {
    reasons <- c()
    if (flagged$pca_outlier[i]) {
      reasons <- c(reasons, "Extreme in PCA space")
    }
    if (flagged$correlation_outlier[i]) {
      reasons <- c(reasons, sprintf("Low correlation (%.3f)", mean_cors[i]))
    }
    flagged$reason[i] <- paste(reasons, collapse = "; ")
  }

  flagged$is_outlier <- flagged$pca_outlier | flagged$correlation_outlier

  n_outliers <- sum(flagged$is_outlier)
  log_message("Flagged ", n_outliers, " potential outlier samples (NOT auto-removed)")

  if (n_outliers > 0) {
    log_message("Outliers: ", paste(flagged$sample[flagged$is_outlier], collapse = ", "))
  }

  list(
    flagged_samples = flagged,
    n_outliers = n_outliers
  )
}

#' Analyze batch effects
#'
#' @param mat Imputed matrix
#' @param pca_results PCA results
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Batch effect analysis results
analyze_batch_effects <- function(mat, pca_results, metadata, config) {
  log_message("Analyzing batch effects...")

  batch_col <- config$design$batch_column
  sample_col <- config$input$sample_id_column

  if (is.null(batch_col) || !batch_col %in% colnames(metadata)) {
    log_message("No batch column specified. Skipping batch analysis.")
    return(NULL)
  }

  batches <- metadata[[batch_col]][match(colnames(mat), metadata[[sample_col]])]

  # Variance explained by batch vs condition
  pca_coords <- pca_results$coordinates

  # ANOVA for PC1 and PC2 vs batch
  pca_coords$batch <- batches

  batch_var_pc1 <- summary(aov(PC1 ~ batch, data = pca_coords))[[1]][1, "Sum Sq"] /
                   sum(summary(aov(PC1 ~ batch, data = pca_coords))[[1]][, "Sum Sq"])
  batch_var_pc2 <- summary(aov(PC2 ~ batch, data = pca_coords))[[1]][1, "Sum Sq"] /
                   sum(summary(aov(PC2 ~ batch, data = pca_coords))[[1]][, "Sum Sq"])

  batch_summary <- data.frame(
    metric = c("Batch variance in PC1", "Batch variance in PC2", "Number of batches"),
    value = c(round(batch_var_pc1 * 100, 2), round(batch_var_pc2 * 100, 2), length(unique(batches))),
    stringsAsFactors = FALSE
  )

  # Create batch effect visualization
  p <- ggplot2::ggplot(pca_coords, ggplot2::aes(x = PC1, y = PC2, color = batch)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::labs(
      title = "PCA Colored by Batch",
      subtitle = sprintf("Batch explains %.1f%% of PC1, %.1f%% of PC2",
                        batch_var_pc1 * 100, batch_var_pc2 * 100),
      x = paste0("PC1 (", round(pca_results$variance_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(pca_results$variance_explained[2], 1), "%)")
    ) +
    ggplot2::theme_minimal()

  save_plot(p, "batch_effect_pca.png", config, subdir = "qc")
  save_table(batch_summary, "batch_effect_summary.csv", config, "qc")

  list(
    summary = batch_summary,
    plot = p
  )
}
