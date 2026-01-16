# =============================================================================
# Random Forest wrapper using ranger for classification and feature importance
# =============================================================================

#' Run Random Forest classification and report variable importance
#'
#' @param mat Numeric matrix (features x samples)
#' @param metadata Sample metadata data.frame
#' @param config Configuration list
#' @return List containing rf object, importance table, and plots
run_random_forest <- function(mat, metadata, config) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    warning("ranger package not available. Skipping random forest analysis.")
    return(NULL)
  }

  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  # Use biological samples only
  sample_roles <- config$sample_types$sample_type_column
  samples <- metadata[[sample_col]]

  # Align matrix columns to metadata samples
  common <- intersect(colnames(mat), samples)
  if (length(common) < 3) {
    warning("Not enough overlapping samples for RF")
    return(NULL)
  }

  mat_sub <- mat[, common, drop = FALSE]
  y <- factor(metadata[[condition_col]][match(colnames(mat_sub), metadata[[sample_col]])])

  # Prepare data.frame for ranger: features as columns
  df <- as.data.frame(t(mat_sub))
  df$.y <- y

  n_trees <- config$rf$n_trees %||% 500
  importance_type <- ifelse(config$rf$importance %||% "permutation" == "impurity", "impurity", "permutation")

  rf_fit <- tryCatch({
    ranger::ranger(.y ~ ., data = df, num.trees = n_trees, importance = importance_type, seed = config$rf$seed %||% 1234)
  }, error = function(e) {
    warning("Random forest failed: ", e$message); return(NULL)
  })

  if (is.null(rf_fit)) return(NULL)

  imp <- tryCatch({
    impt <- as.data.frame(ranger::importance(rf_fit))
    colnames(impt) <- c("importance")
    impt$feature_id <- rownames(impt)
    impt <- dplyr::arrange(impt, dplyr::desc(importance))
    impt
  }, error = function(e) {
    warning("Failed to extract importance: ", e$message); NULL
  })

  # Plot top features
  plots <- list()
  if (!is.null(imp) && nrow(imp) > 0) {
    top_n <- min(nrow(imp), config$rf$top_n %||% 20)
    p <- ggplot2::ggplot(imp[1:top_n, ], ggplot2::aes(x = reorder(feature_id, importance), y = importance)) +
      ggplot2::geom_col(fill = "#0072B2") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Random Forest Feature Importance", x = "Feature", y = "Importance") +
      ggplot2::theme_minimal()

    plots$importance <- p
    save_plot(p, "rf_feature_importance.png", config, subdir = "plots")
  }

  # Save importance table
  if (!is.null(imp)) {
    save_table(imp, "rf_feature_importance.csv", config, subdir = "tables")
  }

  list(rf = rf_fit, importance = imp, plots = plots)
}
