# =============================================================================
# 04_differential.R — limma / t_test / wilcoxon DE analysis
# =============================================================================

run_differential <- function(normalized_data, ingested_data, config) {
  log_message("=== Starting Differential Analysis ===")

  mat      <- normalized_data$matrix
  metadata <- normalized_data$metadata

  sample_col    <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  method        <- config$differential$method %||% "limma"

  # Ensure metadata rows are aligned with matrix columns
  metadata <- metadata[match(colnames(mat), metadata[[sample_col]]), ]
  condition <- metadata[[condition_col]]

  # Get contrast string (e.g. "2024 - 2013")
  contrast_str <- config$design$contrasts
  if (is.list(contrast_str)) contrast_str <- contrast_str[[1]]

  de_results <- switch(method,
    "limma"    = run_limma(mat, condition, contrast_str, config),
    "t_test"   = run_t_test(mat, condition, contrast_str, config),
    "wilcoxon" = run_wilcoxon(mat, condition, contrast_str, config),
    stop("Unknown DE method: ", method)
  )

  # Annotate with feature metadata
  de_results <- annotate_de_results(de_results, ingested_data$feature_metadata)

  # Thresholds
  fc_thresh <- config$differential$log2fc_threshold   %||% 1
  p_thresh  <- config$differential$adj_pvalue_threshold %||% 0.05

  de_results$significant <- abs(de_results$log2FC) >= fc_thresh &
                            de_results$adj.P.Val < p_thresh
  de_results$direction   <- ifelse(de_results$log2FC > 0, "up", "down")
  de_results$direction[!de_results$significant] <- "ns"

  n_sig <- sum(de_results$significant, na.rm = TRUE)
  log_message("Found ", n_sig, " significant features (",
              sum(de_results$direction == "up", na.rm = TRUE), " up, ",
              sum(de_results$direction == "down", na.rm = TRUE), " down)")

  # Plots
  plots <- de_plots(de_results, fc_thresh, p_thresh, config)

  # Save table
  save_table(de_results, "de_results.csv", config, "tables")

  log_message("=== Differential Analysis Complete (method: ", method, ") ===")

  list(table = de_results, method = method, plots = plots)
}

# ---------- limma ------------------------------------------------------------

run_limma <- function(mat, condition, contrast_str, config) {
  log_message("Running limma")

  condition <- factor(condition)
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- levels(condition)

  # Handle NAs: replace with row means for limma
  mat_imp <- mat
  for (i in seq_len(nrow(mat_imp))) {
    nas <- is.na(mat_imp[i, ])
    if (any(nas)) mat_imp[i, nas] <- mean(mat_imp[i, !nas])
  }

  fit <- limma::lmFit(mat_imp, design)

  contrast_matrix <- limma::makeContrasts(contrasts = contrast_str,
                                           levels = design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")

  data.frame(
    feature_id    = rownames(tt),
    log2FC        = tt$logFC,
    avg_intensity = tt$AveExpr,
    statistic     = tt$t,
    P.Value       = tt$P.Value,
    adj.P.Val     = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

# ---------- t-test -----------------------------------------------------------

run_t_test <- function(mat, condition, contrast_str, config) {
  log_message("Running Welch t-tests")

  # Parse contrast to get group names
  groups <- parse_contrast(contrast_str)
  grpB <- groups$numerator
  grpA <- groups$denominator

  idx_A <- which(condition == grpA)
  idx_B <- which(condition == grpB)

  results <- data.frame(
    feature_id    = rownames(mat),
    log2FC        = NA_real_,
    avg_intensity = NA_real_,
    statistic     = NA_real_,
    P.Value       = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(mat))) {
    vals_A <- mat[i, idx_A]
    vals_B <- mat[i, idx_B]
    mean_A <- mean(vals_A, na.rm = TRUE)
    mean_B <- mean(vals_B, na.rm = TRUE)

    results$log2FC[i]        <- mean_B - mean_A  # data already log-transformed
    results$avg_intensity[i] <- mean(c(vals_A, vals_B), na.rm = TRUE)

    tt <- tryCatch(
      t.test(vals_B, vals_A, var.equal = FALSE),
      error = function(e) NULL
    )
    if (!is.null(tt)) {
      results$statistic[i] <- tt$statistic
      results$P.Value[i]   <- tt$p.value
    }
  }

  results$adj.P.Val <- p.adjust(results$P.Value, method = "BH")
  results
}

# ---------- wilcoxon ---------------------------------------------------------

run_wilcoxon <- function(mat, condition, contrast_str, config) {
  log_message("Running Wilcoxon rank-sum tests")

  groups <- parse_contrast(contrast_str)
  grpB <- groups$numerator
  grpA <- groups$denominator

  idx_A <- which(condition == grpA)
  idx_B <- which(condition == grpB)

  results <- data.frame(
    feature_id    = rownames(mat),
    log2FC        = NA_real_,
    avg_intensity = NA_real_,
    statistic     = NA_real_,
    P.Value       = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(mat))) {
    vals_A <- mat[i, idx_A]
    vals_B <- mat[i, idx_B]

    results$log2FC[i]        <- mean(vals_B, na.rm = TRUE) - mean(vals_A, na.rm = TRUE)
    results$avg_intensity[i] <- mean(c(vals_A, vals_B), na.rm = TRUE)

    wt <- tryCatch(
      wilcox.test(vals_B, vals_A, exact = FALSE),
      error = function(e) NULL
    )
    if (!is.null(wt)) {
      results$statistic[i] <- wt$statistic
      results$P.Value[i]   <- wt$p.value
    }
  }

  results$adj.P.Val <- p.adjust(results$P.Value, method = "BH")
  results
}

# ---------- helpers ----------------------------------------------------------

parse_contrast <- function(contrast_str) {
  parts <- trimws(strsplit(contrast_str, "-")[[1]])
  if (length(parts) != 2)
    stop("Cannot parse contrast string: '", contrast_str,
         "'. Expected format: 'groupB - groupA'")
  list(numerator = parts[1], denominator = parts[2])
}

annotate_de_results <- function(de_results, feature_metadata) {
  if (is.null(feature_metadata) || !"feature_id" %in% colnames(feature_metadata))
    return(de_results)

  merge(de_results, feature_metadata, by = "feature_id", all.x = TRUE)
}

# ---------- plots ------------------------------------------------------------

de_plots <- function(de_results, fc_thresh, p_thresh, config) {
  plots <- list()

  # Volcano plot
  vdf <- de_results
  vdf$neg_log10p <- -log10(vdf$adj.P.Val)
  vdf$color <- "ns"
  vdf$color[vdf$adj.P.Val < p_thresh & vdf$log2FC >=  fc_thresh] <- "up"
  vdf$color[vdf$adj.P.Val < p_thresh & vdf$log2FC <= -fc_thresh] <- "down"

  plots$volcano <- ggplot2::ggplot(
    vdf, ggplot2::aes(x = log2FC, y = neg_log10p, color = color)
  ) +
    ggplot2::geom_point(size = 1.5, alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = c(up = "firebrick", down = "steelblue", ns = "grey60")
    ) +
    ggplot2::geom_hline(yintercept = -log10(p_thresh), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
    ggplot2::labs(title = "Volcano Plot", x = "log2 Fold Change",
                  y = "-log10 adjusted p-value", color = "") +
    ggplot2::theme_minimal()

  save_plot(plots$volcano, "volcano_plot.png", config, subdir = "plots")

  # P-value histogram
  plots$pval_hist <- ggplot2::ggplot(
    de_results, ggplot2::aes(x = P.Value)
  ) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggplot2::labs(title = "P-value Distribution", x = "Raw p-value",
                  y = "Count") +
    ggplot2::theme_minimal()

  save_plot(plots$pval_hist, "pvalue_histogram.png", config, subdir = "plots")

  plots
}
