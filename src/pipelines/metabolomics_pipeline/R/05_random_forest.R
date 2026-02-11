# =============================================================================
# 05_random_forest.R — Feature importance via ranger or randomForest
# =============================================================================

run_random_forest <- function(normalized_data, config) {
  if (!isTRUE(config$rf$run_rf)) {
    log_message("Random forest disabled in config — skipping")
    return(NULL)
  }

  log_message("=== Running Random Forest ===")

  mat      <- normalized_data$matrix
  metadata <- normalized_data$metadata

  sample_col    <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  metadata <- metadata[match(colnames(mat), metadata[[sample_col]]), ]
  condition <- factor(metadata[[condition_col]])

  # Prepare data: samples (rows) x features (cols)
  mat_rf <- t(mat)

  # Impute NAs with column medians
  for (j in seq_len(ncol(mat_rf))) {
    nas <- is.na(mat_rf[, j])
    if (any(nas)) mat_rf[nas, j] <- median(mat_rf[!nas, j], na.rm = TRUE)
  }

  n_trees    <- config$rf$n_trees    %||% 500
  importance <- config$rf$importance  %||% "permutation"
  seed       <- config$rf$seed       %||% 1234
  top_n      <- config$rf$top_n      %||% 20

  use_ranger <- requireNamespace("ranger", quietly = TRUE)

  if (use_ranger) {
    log_message("Using ranger")
    set.seed(seed)
    rf_fit <- ranger::ranger(
      x = mat_rf, y = condition,
      num.trees = n_trees, importance = importance, seed = seed
    )
    imp <- ranger::importance(rf_fit)
  } else if (requireNamespace("randomForest", quietly = TRUE)) {
    log_message("ranger unavailable — using randomForest")
    set.seed(seed)
    # randomForest needs clean column names; keep a mapping
    orig_names <- colnames(mat_rf)
    safe_names <- paste0("F", seq_len(ncol(mat_rf)))
    colnames(mat_rf) <- safe_names
    name_map <- setNames(orig_names, safe_names)

    rf_fit <- randomForest::randomForest(
      x = mat_rf, y = condition, ntree = n_trees,
      importance = TRUE
    )
    imp <- randomForest::importance(rf_fit)[, "MeanDecreaseAccuracy"]
    names(imp) <- name_map[names(imp)]
  } else {
    log_message("No RF package available — skipping")
    return(NULL)
  }

  # Importance table
  imp_df <- data.frame(
    feature_id = names(imp),
    importance = unname(imp),
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[order(imp_df$importance, decreasing = TRUE), ]

  save_table(imp_df, "rf_importance.csv", config, "tables")

  # Importance bar plot (top N)
  top_df <- head(imp_df, top_n)
  top_df$feature_id <- factor(top_df$feature_id,
                               levels = rev(top_df$feature_id))

  p <- ggplot2::ggplot(
    top_df, ggplot2::aes(x = feature_id, y = importance)
  ) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Random Forest — Top ", nrow(top_df), " Features"),
      x = NULL, y = "Importance"
    ) +
    ggplot2::theme_minimal()

  save_plot(p, "rf_importance.png", config, subdir = "plots")

  log_message("=== Random Forest Complete ===")

  list(importance = imp_df, model = rf_fit, plot = p)
}
