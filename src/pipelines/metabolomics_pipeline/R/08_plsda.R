# =============================================================================
# 08_plsda.R — PLS-DA analysis with VIP scores (mixOmics)
# =============================================================================

run_plsda <- function(normalized_data, ingested_data, config) {
  if (!isTRUE(config$plsda$run_plsda)) {
    log_message("PLS-DA disabled in config — skipping")
    return(NULL)
  }

  log_message("=== Starting PLS-DA Analysis ===")

  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    log_message("mixOmics package not available — skipping PLS-DA")
    log_message("Install with: BiocManager::install('mixOmics')")
    return(NULL)
  }

  mat      <- normalized_data$matrix
  metadata <- normalized_data$metadata

  sample_col    <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  # Align metadata with matrix columns
  metadata <- metadata[match(colnames(mat), metadata[[sample_col]]), ]
  condition <- factor(metadata[[condition_col]])

  # Prepare X matrix: samples (rows) x features (cols)
  X <- t(mat)

  # Impute NAs with column medians
  for (j in seq_len(ncol(X))) {
    nas <- is.na(X[, j])
    if (any(nas)) X[nas, j] <- median(X[!nas, j], na.rm = TRUE)
  }

  # Remove zero-variance features
  col_vars <- apply(X, 2, var, na.rm = TRUE)
  keep <- col_vars > 0
  if (sum(keep) < ncol(X)) {
    log_message("Removed ", sum(!keep), " zero-variance features for PLS-DA")
    X <- X[, keep, drop = FALSE]
  }

  # Fit PLS-DA (2 components for scores plot)
  ncomp <- min(2, ncol(X), nrow(X) - 1)
  log_message("Fitting PLS-DA with ", ncomp, " components on ",
              nrow(X), " samples x ", ncol(X), " features")

  plsda_model <- mixOmics::plsda(X, condition, ncomp = ncomp)

  # Explained variance per component (mixOmics uses prop_expl_var)
  expl_var <- plsda_model$prop_expl_var$X

  # VIP scores
  vip_scores <- mixOmics::vip(plsda_model)

  # PLS-DA config subsection
  plsda_cfg <- config$plsda %||% list()

  # Build feature name lookup from ingested_data
  feat_names <- NULL
  if (!is.null(ingested_data$feature_metadata) &&
      "Molecule" %in% colnames(ingested_data$feature_metadata)) {
    fm <- ingested_data$feature_metadata
    feat_names <- setNames(fm$Molecule, fm$feature_id)
  }

  # Generate plots
  scores_plot <- plot_plsda_scores(plsda_model, condition, expl_var,
                                   plsda_cfg, config)
  vip_plot    <- plot_vip_scores(vip_scores, X, condition, feat_names,
                                 plsda_cfg, config)

  # Save VIP table
  vip_df <- data.frame(
    feature_id = rownames(vip_scores),
    VIP_comp1  = vip_scores[, 1],
    stringsAsFactors = FALSE
  )
  if (ncol(vip_scores) >= 2) vip_df$VIP_comp2 <- vip_scores[, 2]

  # Add molecule names if available
  if (!is.null(feat_names)) {
    vip_df$Molecule <- feat_names[vip_df$feature_id]
  }

  vip_df <- vip_df[order(vip_df$VIP_comp1, decreasing = TRUE), ]
  save_table(vip_df, "plsda_vip_scores.csv", config, "tables")

  log_message("=== PLS-DA Analysis Complete ===")

  list(
    model              = plsda_model,
    vip_scores         = vip_scores,
    explained_variance = expl_var,
    plots              = list(scores = scores_plot, vip = vip_plot)
  )
}

# ---------- PLS-DA Scores Plot ------------------------------------------------

plot_plsda_scores <- function(model, condition, expl_var, plsda_cfg, config) {
  log_message("Generating PLS-DA scores plot")

  scores  <- model$variates$X
  pct_var <- round(expl_var * 100, 1)

  df <- data.frame(
    Comp1     = scores[, 1],
    Comp2     = scores[, 2],
    condition = condition,
    stringsAsFactors = FALSE
  )

  # Configurable colors
  default_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                      "#FF7F00", "#A65628")
  colors <- plsda_cfg$colors %||% default_colors
  n_groups <- nlevels(condition)
  colors <- colors[seq_len(n_groups)]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Comp1, y = Comp2,
                                         color = condition)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.85) +
    ggplot2::scale_color_manual(values = colors)

  # Add 95% confidence ellipses (type="norm" works with 3+ samples)
  grp_counts <- table(condition)
  if (all(grp_counts >= 3)) {
    p <- p + ggplot2::stat_ellipse(type = "norm", level = 0.95,
                                    size = 0.8, show.legend = FALSE)
  } else {
    log_message("Fewer than 3 samples in a group — skipping confidence ellipses")
  }

  p <- p +
    ggplot2::labs(
      title = "PLS-DA Scores Plot",
      x = paste0("Component 1 (", pct_var[1], "%)"),
      y = paste0("Component 2 (", pct_var[2], "%)"),
      color = ""
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold",
                                           size = 14),
      axis.title   = ggplot2::element_text(size = 12),
      legend.text  = ggplot2::element_text(size = 11),
      legend.position = "bottom"
    )

  save_plot(p, "plsda_scores.png", config, subdir = "plots")
  p
}

# ---------- VIP Scores Plot ---------------------------------------------------

plot_vip_scores <- function(vip_scores, X, condition, feat_names,
                            plsda_cfg, config) {
  log_message("Generating VIP scores plot")

  top_n <- plsda_cfg$vip_top_n %||% 15

  # Use component 1 VIP scores
  vip_vec    <- vip_scores[, 1]
  vip_sorted <- sort(vip_vec, decreasing = TRUE)
  top_feats  <- names(vip_sorted)[seq_len(min(top_n, length(vip_sorted)))]

  # Determine which group has the highest mean abundance per feature
  groups <- levels(condition)
  high_group <- vapply(top_feats, function(feat) {
    means <- tapply(X[, feat], condition, mean, na.rm = TRUE)
    groups[which.max(means)]
  }, character(1))

  # Display names: use Molecule names if available
  display_names <- if (!is.null(feat_names)) {
    ifelse(is.na(feat_names[top_feats]) | feat_names[top_feats] == "",
           top_feats, feat_names[top_feats])
  } else {
    top_feats
  }

  df <- data.frame(
    feature      = factor(display_names, levels = rev(display_names)),
    VIP          = vip_sorted[top_feats],
    high_in      = high_group,
    stringsAsFactors = FALSE
  )

  # Colors matching scores plot
  default_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                      "#FF7F00", "#A65628")
  colors <- plsda_cfg$colors %||% default_colors
  n_groups <- length(groups)
  group_colors <- setNames(colors[seq_len(n_groups)], groups)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = VIP, y = feature,
                                         color = high_in)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = VIP, y = feature, yend = feature),
      color = "grey70", size = 0.5
    ) +
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_color_manual(
      values = group_colors,
      labels = paste("High in", names(group_colors))
    ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                         color = "grey40", size = 0.4) +
    ggplot2::labs(
      title = paste0("PLS-DA VIP Scores \u2014 Top ", nrow(df), " Features"),
      x     = "VIP Score",
      y     = NULL,
      color = ""
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold",
                                           size = 14),
      axis.title   = ggplot2::element_text(size = 12),
      axis.text.y  = ggplot2::element_text(size = 9),
      legend.text  = ggplot2::element_text(size = 11),
      legend.position = "bottom"
    )

  save_plot(p, "plsda_vip_scores.png", config, subdir = "plots")
  p
}
