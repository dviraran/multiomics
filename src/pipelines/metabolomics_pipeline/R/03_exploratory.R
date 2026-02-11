# =============================================================================
# 03_exploratory.R — PCA, correlation heatmap, sample heatmap
# =============================================================================

run_exploratory <- function(normalized_data, config) {
  log_message("=== Starting Exploratory Analysis ===")

  mat      <- normalized_data$matrix
  metadata <- normalized_data$metadata

  sample_col    <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  cond_vec <- as.character(
    metadata[[condition_col]][match(colnames(mat), metadata[[sample_col]])]
  )
  names(cond_vec) <- colnames(mat)

  results <- list()

  # PCA
  if (isTRUE(config$qc$run_pca)) {
    results$pca <- run_pca(mat, cond_vec, config)
  }

  # Correlation heatmap
  if (isTRUE(config$qc$run_correlation)) {
    results$correlation <- plot_correlation_heatmap(mat, cond_vec, config)
  }

  # Sample heatmap (top variable features)
  if (isTRUE(config$qc$run_sample_heatmap)) {
    top_n <- config$qc$top_n_variable %||% 50
    results$sample_heatmap <- plot_sample_heatmap(mat, cond_vec, top_n, config)
  }

  log_message("=== Exploratory Analysis Complete ===")
  results
}

# ---------- PCA --------------------------------------------------------------

run_pca <- function(mat, cond_vec, config) {
  log_message("Running PCA")

  # Remove features with zero variance or all NA
  keep <- apply(mat, 1, function(x) {
    vals <- x[!is.na(x)]
    length(vals) >= 2 && var(vals) > 0
  })
  mat_pca <- mat[keep, ]

  # Impute remaining NA with row means for PCA
  for (i in seq_len(nrow(mat_pca))) {
    nas <- is.na(mat_pca[i, ])
    if (any(nas)) mat_pca[i, nas] <- mean(mat_pca[i, !nas])
  }

  pca <- prcomp(t(mat_pca), scale. = FALSE)
  pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  pca_df <- data.frame(
    PC1       = pca$x[, 1],
    PC2       = pca$x[, 2],
    sample    = rownames(pca$x),
    condition = cond_vec[rownames(pca$x)],
    stringsAsFactors = FALSE
  )

  # Scores plot
  p_scores <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, y = PC2, color = condition, label = sample)
  ) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 15) +
    ggplot2::labs(
      title = "PCA Scores Plot",
      x = paste0("PC1 (", pct_var[1], "%)"),
      y = paste0("PC2 (", pct_var[2], "%)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(color = "")

  save_plot(p_scores, "pca_scores.png", config, subdir = "qc")

  # Scree plot
  scree_df <- data.frame(
    PC = seq_along(pct_var),
    variance = pct_var
  )
  scree_df <- scree_df[1:min(10, nrow(scree_df)), ]

  p_scree <- ggplot2::ggplot(
    scree_df, ggplot2::aes(x = PC, y = variance)
  ) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(title = "PCA Scree Plot", x = "PC", y = "% Variance") +
    ggplot2::scale_x_continuous(breaks = scree_df$PC) +
    ggplot2::theme_minimal()

  save_plot(p_scree, "pca_scree.png", config, width = 7, height = 5, subdir = "qc")

  log_message("PCA complete")
  list(pca = pca, pct_var = pct_var, plots = list(scores = p_scores, scree = p_scree))
}

# ---------- correlation heatmap ----------------------------------------------

plot_correlation_heatmap <- function(mat, cond_vec, config) {
  log_message("Generating correlation heatmap")

  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")

  annot_df <- data.frame(
    condition = cond_vec[colnames(cor_mat)],
    row.names = colnames(cor_mat)
  )

  # Save to file via pheatmap's filename argument
  out_path <- file.path(config$output$output_dir, "qc",
                         "correlation_heatmap.png")
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  pheatmap::pheatmap(
    cor_mat,
    cluster_rows = TRUE, cluster_cols = TRUE,
    annotation_col = annot_df,
    main = "Sample Correlation Heatmap (Pearson)",
    fontsize_row = 7, fontsize_col = 7,
    filename = out_path,
    width = config$output$fig_width %||% 10,
    height = config$output$fig_height %||% 8
  )

  log_message("Saved correlation heatmap: ", out_path)
  list(cor_matrix = cor_mat, path = out_path)
}

# ---------- sample heatmap ---------------------------------------------------

plot_sample_heatmap <- function(mat, cond_vec, top_n, config) {
  log_message("Generating sample heatmap (top ", top_n, " variable features)")

  # Select top-N most variable features
  row_vars <- apply(mat, 1, var, na.rm = TRUE)
  top_features <- names(sort(row_vars, decreasing = TRUE))[
    seq_len(min(top_n, length(row_vars)))
  ]
  sub_mat <- mat[top_features, , drop = FALSE]

  # Row-scale
  sub_scaled <- t(scale(t(sub_mat)))

  annot_df <- data.frame(
    condition = cond_vec[colnames(sub_scaled)],
    row.names = colnames(sub_scaled)
  )

  out_path <- file.path(config$output$output_dir, "qc",
                         "sample_heatmap.png")
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  pheatmap::pheatmap(
    sub_scaled,
    cluster_rows = TRUE, cluster_cols = TRUE,
    annotation_col = annot_df,
    show_rownames = FALSE,
    scale = "none",
    main = paste0("Top ", length(top_features), " Variable Features"),
    fontsize_col = 7,
    filename = out_path,
    width = config$output$fig_width %||% 10,
    height = config$output$fig_height %||% 8
  )

  log_message("Saved sample heatmap: ", out_path)
  list(path = out_path)
}
