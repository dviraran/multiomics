# R/05_qc.R
# Functions for QC metrics, PCA, and exploratory analysis

#' Compute QC metrics for samples
#' @param counts Raw or filtered counts matrix
#' @param normalized_counts Normalized counts matrix
#' @param metadata Sample metadata
#' @return Data frame of QC metrics per sample
compute_qc_metrics <- function(counts, normalized_counts, metadata) {

  samples <- colnames(counts)

  qc <- data.frame(
    sample = samples,
    total_counts = colSums(counts),
    total_normalized = colSums(normalized_counts),
    detected_genes = colSums(counts > 0),
    detection_rate = colSums(counts > 0) / nrow(counts),
    median_counts = apply(counts, 2, median),
    mean_counts = colMeans(counts),
    max_counts = apply(counts, 2, max),
    stringsAsFactors = FALSE
  )

  # Add coefficient of variation
  qc$cv <- apply(counts, 2, sd) / qc$mean_counts

  # Calculate complexity (genes needed to reach 50% of counts)
  qc$complexity <- sapply(samples, function(s) {
    sorted_counts <- sort(counts[, s], decreasing = TRUE)
    cumsum_counts <- cumsum(sorted_counts)
    threshold <- sum(sorted_counts) * 0.5
    sum(cumsum_counts <= threshold) + 1
  })

  message("Computed QC metrics for ", nrow(qc), " samples")

  qc
}


#' Compute sample-to-sample correlation matrix
#' @param vst_counts VST-transformed counts
#' @param method Correlation method
#' @return Correlation matrix
compute_sample_correlation <- function(vst_counts, method = "pearson") {
  cor_mat <- cor(vst_counts, method = method)
  message("Computed ", method, " correlation matrix: ", nrow(cor_mat), " x ", ncol(cor_mat))
  cor_mat
}


#' Compute PCA on VST counts
#' @param vst_counts VST-transformed counts
#' @param metadata Sample metadata
#' @param sample_id_col Sample ID column name
#' @param n_top_genes Number of most variable genes for PCA
#' @return List with PCA results
compute_pca <- function(vst_counts,
                        metadata,
                        sample_id_col = "sample_id",
                        n_top_genes = 500) {

  # Select top variable genes
  gene_vars <- matrixStats::rowVars(vst_counts)
  top_genes <- order(gene_vars, decreasing = TRUE)[1:min(n_top_genes, length(gene_vars))]
  vst_subset <- vst_counts[top_genes, ]

  # Run PCA
  pca <- prcomp(t(vst_subset), center = TRUE, scale. = TRUE)

  # Extract variance explained
  var_explained <- summary(pca)$importance[2, ] * 100

  # Create results data frame
  pca_data <- as.data.frame(pca$x)
  pca_data$sample <- rownames(pca_data)

  # Merge with metadata
  meta_cols <- setdiff(colnames(metadata), sample_id_col)
  pca_data <- cbind(pca_data, metadata[match(pca_data$sample, rownames(metadata)), meta_cols, drop = FALSE])

  message("PCA computed using top ", length(top_genes), " variable genes")
  message("PC1: ", round(var_explained[1], 1), "%, PC2: ", round(var_explained[2], 1), "%")

  list(
    pca_data = pca_data,
    var_explained = var_explained,
    pca_object = pca,
    n_genes_used = length(top_genes)
  )
}


#' Detect outlier samples based on PCA or correlation
#' @param vst_counts VST-transformed counts
#' @param metadata Sample metadata
#' @param sample_id_col Sample ID column
#' @param sd_threshold SD threshold for outlier detection
#' @return Data frame with outlier flags
detect_outliers <- function(vst_counts,
                            metadata,
                            sample_id_col = "sample_id",
                            sd_threshold = 3) {

  samples <- colnames(vst_counts)

  # Method 1: PCA-based (distance from center in PC1-PC2 space)
  pca <- prcomp(t(vst_counts), center = TRUE, scale. = TRUE)
  pc_scores <- pca$x[, 1:min(5, ncol(pca$x))]

  # Compute Mahalanobis-like distance (simplified: Euclidean in scaled PC space)
  center <- colMeans(pc_scores)
  distances <- sqrt(rowSums((sweep(pc_scores, 2, center))^2))

  # Flag outliers
  dist_mean <- mean(distances)
  dist_sd <- sd(distances)
  pca_outlier <- distances > (dist_mean + sd_threshold * dist_sd)

  # Method 2: Correlation-based (low median correlation with other samples)
  cor_mat <- cor(vst_counts, method = "pearson")
  diag(cor_mat) <- NA
  median_cor <- apply(cor_mat, 1, median, na.rm = TRUE)

  cor_mean <- mean(median_cor)
  cor_sd <- sd(median_cor)
  cor_outlier <- median_cor < (cor_mean - sd_threshold * cor_sd)

  # Combine results
  outliers <- data.frame(
    sample = samples,
    pca_distance = distances,
    pca_outlier = pca_outlier,
    median_correlation = median_cor,
    correlation_outlier = cor_outlier,
    is_outlier = pca_outlier | cor_outlier,
    stringsAsFactors = FALSE
  )

  n_outliers <- sum(outliers$is_outlier)
  if (n_outliers > 0) {
    message("Flagged ", n_outliers, " potential outlier(s): ",
            paste(outliers$sample[outliers$is_outlier], collapse = ", "))
  } else {
    message("No outliers detected (SD threshold = ", sd_threshold, ")")
  }

  outliers
}


#' Generate QC plots
#' @param qc_metrics QC metrics data frame
#' @param sample_correlation Correlation matrix
#' @param pca_result PCA results
#' @param metadata Sample metadata
#' @param group_col Column for color grouping
#' @param output_dir Output directory for plots
#' @return List of plot file paths
generate_qc_plots <- function(qc_metrics,
                               sample_correlation,
                               pca_result,
                               metadata,
                               group_col = "condition",
                               output_dir = "outputs/plots") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  plot_files <- list()

  # 1. Library size barplot
  p1 <- ggplot(qc_metrics, aes(x = reorder(sample, -total_counts), y = total_counts / 1e6)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Library Sizes", x = "Sample", y = "Total Counts (millions)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  ggsave(file.path(output_dir, "library_sizes.png"), p1, width = 10, height = 6)
  plot_files$library_sizes <- file.path(output_dir, "library_sizes.png")

  # 2. Detected genes barplot
  p2 <- ggplot(qc_metrics, aes(x = reorder(sample, -detected_genes), y = detected_genes)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    labs(title = "Detected Genes per Sample", x = "Sample", y = "Number of Genes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  ggsave(file.path(output_dir, "detected_genes.png"), p2, width = 10, height = 6)
  plot_files$detected_genes <- file.path(output_dir, "detected_genes.png")

  # 3. Sample correlation heatmap
  png(file.path(output_dir, "sample_correlation_heatmap.png"), width = 800, height = 800)
  pheatmap::pheatmap(
    sample_correlation,
    main = "Sample-to-Sample Correlation",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    border_color = NA,
    fontsize = 8
  )
  dev.off()
  plot_files$correlation_heatmap <- file.path(output_dir, "sample_correlation_heatmap.png")

  # 4. PCA plot
  pca_data <- pca_result$pca_data
  var_exp <- pca_result$var_explained

  if (group_col %in% colnames(pca_data)) {
    p4 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[group_col]])) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = 15) +
      labs(
        title = "PCA Plot",
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)")
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  } else {
    p4 <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
      geom_point(size = 3, color = "steelblue") +
      ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = 15) +
      labs(
        title = "PCA Plot",
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)")
      ) +
      theme_minimal()
  }

  ggsave(file.path(output_dir, "pca_plot.png"), p4, width = 10, height = 8)
  plot_files$pca <- file.path(output_dir, "pca_plot.png")

  # 5. PCA scree plot
  var_df <- data.frame(
    PC = paste0("PC", 1:min(10, length(var_exp))),
    variance = var_exp[1:min(10, length(var_exp))]
  )
  var_df$PC <- factor(var_df$PC, levels = var_df$PC)

  p5 <- ggplot(var_df, aes(x = PC, y = variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_line(aes(group = 1), color = "red") +
    geom_point(color = "red", size = 2) +
    labs(title = "PCA Variance Explained", x = "Principal Component", y = "% Variance") +
    theme_minimal()

  ggsave(file.path(output_dir, "pca_scree.png"), p5, width = 8, height = 5)
  plot_files$pca_scree <- file.path(output_dir, "pca_scree.png")

  message("Generated ", length(plot_files), " QC plots in ", output_dir)

  plot_files
}
