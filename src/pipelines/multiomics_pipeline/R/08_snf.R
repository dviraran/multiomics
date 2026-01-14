# =============================================================================
# Similarity Network Fusion (SNF) Integration (Optional)
# =============================================================================

#' Run SNF integration
run_snf_integration <- function(feature_data, config) {
  log_message("=== Running SNF Integration ===")

  if (!requireNamespace("SNFtool", quietly = TRUE)) {
    log_message("SNFtool package not available. Skipping SNF integration.")
    return(NULL)
  }

  # Extract parameters
  snf_config <- config$integration$snf
  K <- snf_config$K %||% 20
  alpha <- snf_config$alpha %||% 0.5
  n_clusters <- snf_config$n_clusters %||% 2

  # Get data
  filtered_mae <- feature_data$filtered_mae
  matrices <- extract_matrices_for_integration(filtered_mae)
  metadata <- get_sample_metadata(filtered_mae)

  if (length(matrices) < 2) {
    log_message("SNF requires at least 2 omics layers. Found: ", length(matrices))
    return(NULL)
  }

  log_message("Running SNF with K=", K, ", alpha=", alpha, ", clusters=", n_clusters)

  # Prepare data (samples x features)
  data_list <- list()
  for (omic in names(matrices)) {
    mat <- t(matrices[[omic]])  # Transpose to samples x features
    # Impute NAs with column means
    for (j in seq_len(ncol(mat))) {
      na_idx <- is.na(mat[, j])
      if (any(na_idx)) {
        mat[na_idx, j] <- mean(mat[, j], na.rm = TRUE)
      }
    }
    data_list[[omic]] <- mat
  }

  # Compute similarity matrices for each omics
  log_message("Computing similarity matrices...")
  similarity_matrices <- list()

  for (omic in names(data_list)) {
    # Distance matrix
    dist_mat <- as.matrix(dist(data_list[[omic]], method = "euclidean"))

    # Affinity matrix
    W <- SNFtool::affinityMatrix(dist_mat, K = K, sigma = alpha)
    similarity_matrices[[omic]] <- W

    log_message("  ", omic, ": ", nrow(W), " x ", ncol(W), " similarity matrix")
  }

  # Fuse networks
  log_message("Fusing similarity networks...")
  W_fused <- tryCatch({
    SNFtool::SNF(similarity_matrices, K = K, t = 20)
  }, error = function(e) {
    log_message("Error in SNF: ", e$message)
    return(NULL)
  })

  if (is.null(W_fused)) return(NULL)

  # Spectral clustering on fused network
  log_message("Performing spectral clustering...")
  clusters <- SNFtool::spectralClustering(W_fused, K = n_clusters)

  # Add cluster assignments to results
  sample_ids <- rownames(data_list[[1]])
  cluster_df <- data.frame(
    sample_id = sample_ids,
    snf_cluster = clusters,
    stringsAsFactors = FALSE
  )

  # Add condition if available
  condition_col <- config$design$condition_column
  if (condition_col %in% colnames(metadata)) {
    cluster_df[[condition_col]] <- metadata[sample_ids, condition_col]
  }

  save_table(cluster_df, "snf_clusters.csv", config)

  # Calculate NMI with known labels if available
  nmi_score <- NULL
  if (condition_col %in% colnames(metadata)) {
    true_labels <- as.numeric(factor(metadata[sample_ids, condition_col]))
    nmi_score <- calculate_nmi(clusters, true_labels)
    log_message("SNF cluster-condition NMI: ", round(nmi_score, 3))
  }

  # Save fused similarity matrix
  W_df <- as.data.frame(W_fused)
  colnames(W_df) <- sample_ids
  W_df <- tibble::rownames_to_column(W_df, "sample_id")
  save_table(W_df, "snf_fused_similarity.csv", config)

  # Create visualizations
  snf_results <- list(
    fused_similarity = W_fused,
    individual_similarities = similarity_matrices,
    clusters = clusters,
    cluster_df = cluster_df,
    nmi = nmi_score,
    K = K,
    alpha = alpha
  )

  create_snf_plots(snf_results, metadata, config)

  log_message("SNF integration complete")

  snf_results
}

#' Calculate Normalized Mutual Information
calculate_nmi <- function(pred_labels, true_labels) {
  # Simple NMI calculation
  n <- length(pred_labels)

  # Entropy of predicted labels
  p_pred <- table(pred_labels) / n
  H_pred <- -sum(p_pred * log2(p_pred + 1e-10))

  # Entropy of true labels
  p_true <- table(true_labels) / n
  H_true <- -sum(p_true * log2(p_true + 1e-10))

  # Joint entropy
  joint <- table(pred_labels, true_labels) / n
  H_joint <- -sum(joint * log2(joint + 1e-10))

  # Mutual information
  MI <- H_pred + H_true - H_joint

  # Normalized MI
  NMI <- 2 * MI / (H_pred + H_true + 1e-10)

  NMI
}

#' Create SNF visualizations
create_snf_plots <- function(snf_results, metadata, config) {
  log_message("Creating SNF plots...")

  condition_col <- config$design$condition_column

  # 1. Heatmap of fused similarity matrix
  tryCatch({
    p <- plot_snf_heatmap(snf_results, metadata, condition_col)
    save_plot(p, "snf_fused_heatmap", config, width = 10, height = 10)
  }, error = function(e) log_message("Failed to create SNF heatmap: ", e$message))

  # 2. MDS/UMAP of fused network
  tryCatch({
    p <- plot_snf_mds(snf_results, metadata, condition_col)
    save_plot(p, "snf_mds", config, width = 8, height = 6)
  }, error = function(e) log_message("Failed to create SNF MDS plot: ", e$message))

  # 3. Cluster-condition agreement
  tryCatch({
    if (!is.null(snf_results$cluster_df) && condition_col %in% colnames(snf_results$cluster_df)) {
      p <- plot_snf_cluster_agreement(snf_results$cluster_df, condition_col)
      save_plot(p, "snf_cluster_agreement", config, width = 8, height = 6)
    }
  }, error = function(e) log_message("Failed to create cluster agreement plot: ", e$message))

  log_message("SNF plots complete")
}

#' Plot SNF heatmap
plot_snf_heatmap <- function(snf_results, metadata, condition_col) {
  W <- snf_results$fused_similarity
  clusters <- snf_results$clusters

  # Order by cluster
  ord <- order(clusters)
  W_ordered <- W[ord, ord]

  # Create annotation
  sample_ids <- rownames(W_ordered)
  anno_df <- data.frame(
    cluster = factor(clusters[ord]),
    row.names = sample_ids
  )

  if (condition_col %in% colnames(metadata)) {
    anno_df[[condition_col]] <- metadata[sample_ids, condition_col]
  }

  # Convert to long format for ggplot
  W_df <- as.data.frame(W_ordered)
  W_df$sample1 <- factor(rownames(W_df), levels = rownames(W_df))
  W_long <- tidyr::pivot_longer(W_df, -sample1, names_to = "sample2", values_to = "similarity")
  W_long$sample2 <- factor(W_long$sample2, levels = levels(W_long$sample1))

  ggplot2::ggplot(W_long, ggplot2::aes(x = sample1, y = sample2, fill = similarity)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "white", high = "darkblue", name = "Similarity") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "SNF: Fused Similarity Matrix",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
      axis.text.y = ggplot2::element_text(size = 6)
    )
}

#' Plot SNF MDS
plot_snf_mds <- function(snf_results, metadata, condition_col) {
  W <- snf_results$fused_similarity
  clusters <- snf_results$clusters

  # Convert similarity to distance
  D <- 1 - W / max(W)
  diag(D) <- 0

  # MDS
  mds <- cmdscale(D, k = 2)

  df <- data.frame(
    MDS1 = mds[, 1],
    MDS2 = mds[, 2],
    cluster = factor(clusters),
    sample_id = rownames(W),
    stringsAsFactors = FALSE
  )

  if (condition_col %in% colnames(metadata)) {
    df[[condition_col]] <- metadata[df$sample_id, condition_col]
  }

  # Plot with both cluster shape and condition color
  if (condition_col %in% colnames(df)) {
    ggplot2::ggplot(df, ggplot2::aes(x = MDS1, y = MDS2,
                                      color = .data[[condition_col]],
                                      shape = cluster)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "SNF: MDS of Fused Network",
        x = "MDS 1",
        y = "MDS 2",
        color = condition_col,
        shape = "SNF Cluster"
      )
  } else {
    ggplot2::ggplot(df, ggplot2::aes(x = MDS1, y = MDS2, color = cluster)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "SNF: MDS of Fused Network",
        x = "MDS 1",
        y = "MDS 2",
        color = "SNF Cluster"
      )
  }
}

#' Plot cluster-condition agreement
plot_snf_cluster_agreement <- function(cluster_df, condition_col) {
  # Contingency table
  tbl <- table(cluster_df$snf_cluster, cluster_df[[condition_col]])

  # Convert to data frame
  df <- as.data.frame(tbl)
  colnames(df) <- c("snf_cluster", "condition", "count")

  ggplot2::ggplot(df, ggplot2::aes(x = snf_cluster, y = condition, fill = count)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = count), size = 5) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "SNF Cluster vs Condition Agreement",
      x = "SNF Cluster",
      y = condition_col,
      fill = "Count"
    )
}

#' Estimate optimal number of clusters for SNF
estimate_snf_clusters <- function(W_fused, max_k = 10) {
  # Use eigenvalues of similarity matrix
  eigvals <- eigen(W_fused, only.values = TRUE)$values
  eigvals <- sort(eigvals, decreasing = TRUE)

  # Look for eigengap
  diffs <- diff(eigvals[1:min(max_k, length(eigvals))])
  optimal_k <- which.max(abs(diffs)) + 1

  list(
    optimal_k = optimal_k,
    eigenvalues = eigvals[1:max_k]
  )
}
