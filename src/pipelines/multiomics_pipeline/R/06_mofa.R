# =============================================================================
# MOFA2 Integration
# =============================================================================

#' Run MOFA2 integration
run_mofa2_integration <- function(feature_data, config) {
  log_message("=== Running MOFA2 Integration ===")

  if (!requireNamespace("MOFA2", quietly = TRUE)) {
    log_message("MOFA2 package not available. Skipping MOFA2 integration.")
    return(NULL)
  }

  # Extract parameters
  mofa_config <- config$integration$mofa2
  num_factors <- mofa_config$num_factors %||% 10
  convergence_mode <- mofa_config$convergence_mode %||% "fast"
  seed <- mofa_config$seed %||% 42

  # Get matrices
  filtered_mae <- feature_data$filtered_mae
  matrices <- extract_matrices_for_integration(filtered_mae)
  metadata <- get_sample_metadata(filtered_mae)

  # Check we have at least 2 omics
  if (length(matrices) < 2) {
    log_message("MOFA2 requires at least 2 omics layers. Found: ", length(matrices))
    return(NULL)
  }

  # Prepare data for MOFA2
  # MOFA2 expects features x samples matrices
  log_message("Preparing data for MOFA2...")

  # Create MOFA object from matrices
  mofa_data <- list()
  for (omic in names(matrices)) {
    mat <- matrices[[omic]]
    # Center features (MOFA2 expects centered data for Gaussian likelihood)
    mat_centered <- t(scale(t(mat), center = TRUE, scale = FALSE))
    # Handle remaining NAs
    mat_centered[is.na(mat_centered)] <- 0
    mofa_data[[omic]] <- mat_centered
  }

  # Create MOFA object
  mofa_obj <- tryCatch({
    MOFA2::create_mofa(mofa_data)
  }, error = function(e) {
    log_message("Error creating MOFA object: ", e$message)
    return(NULL)
  })

  if (is.null(mofa_obj)) return(NULL)

  # Prepare options using new API (MOFA2 >= 1.10)
  data_opts <- MOFA2::get_default_data_options(mofa_obj)

  model_opts <- MOFA2::get_default_model_options(mofa_obj)
  model_opts$num_factors <- min(num_factors, ncol(matrices[[1]]) - 1)

  train_opts <- MOFA2::get_default_training_options(mofa_obj)
  train_opts$convergence_mode <- convergence_mode
  train_opts$seed <- seed
  train_opts$verbose <- FALSE

  # Use prepare_mofa to set all options at once (new API)
  mofa_obj <- tryCatch({
    MOFA2::prepare_mofa(
      mofa_obj,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )
  }, error = function(e) {
    log_message("Error preparing MOFA object: ", e$message)
    return(NULL)
  })

  if (is.null(mofa_obj)) return(NULL)

  # Run MOFA
  log_message("Training MOFA2 model with ", model_opts$num_factors, " factors...")
  mofa_trained <- tryCatch({
    MOFA2::run_mofa(mofa_obj, use_basilisk = FALSE)
  }, error = function(e) {
    log_message("Error training MOFA model: ", e$message)
    # Try with basilisk
    tryCatch({
      MOFA2::run_mofa(mofa_obj, use_basilisk = TRUE)
    }, error = function(e2) {
      log_message("MOFA training failed: ", e2$message)
      return(NULL)
    })
  })

  if (is.null(mofa_trained)) return(NULL)

  log_message("MOFA2 training complete")

  # Extract results
  mofa_results <- extract_mofa_results(mofa_trained, metadata, config)

  # Create visualizations
  create_mofa_plots(mofa_trained, mofa_results, metadata, config)

  list(
    model = mofa_trained,
    results = mofa_results,
    config = mofa_config
  )
}

#' Extract results from trained MOFA model
extract_mofa_results <- function(mofa_model, metadata, config) {
  log_message("Extracting MOFA2 results...")

  # Factor values (samples x factors)
  factors <- MOFA2::get_factors(mofa_model)[[1]]

  # Weights per view (features x factors per view)
  weights <- MOFA2::get_weights(mofa_model)

  # Variance explained
  var_explained <- MOFA2::get_variance_explained(mofa_model)

  # Per factor per view
  r2_per_factor <- var_explained$r2_per_factor[[1]]

  # Total per view
  r2_total <- var_explained$r2_total[[1]]

  # Save factor values
  factor_df <- as.data.frame(factors)
  factor_df$sample_id <- rownames(factors)
  factor_df <- factor_df[, c("sample_id", setdiff(colnames(factor_df), "sample_id"))]

  # Add metadata
  condition_col <- config$design$condition_column
  if (condition_col %in% colnames(metadata)) {
    factor_df[[condition_col]] <- metadata[factor_df$sample_id, condition_col]
  }
  save_table(factor_df, "mofa_factors.csv", config)

  # Save weights per view
  for (view in names(weights)) {
    w <- weights[[view]]
    w_df <- as.data.frame(w)
    w_df$feature_id <- rownames(w)
    w_df <- w_df[, c("feature_id", setdiff(colnames(w_df), "feature_id"))]
    save_table(w_df, paste0("mofa_weights_", view, ".csv"), config)
  }

  # Save variance explained
  var_df <- as.data.frame(r2_per_factor)
  var_df$view <- rownames(var_df)
  var_df <- var_df[, c("view", setdiff(colnames(var_df), "view"))]
  save_table(var_df, "mofa_variance_explained.csv", config)

  # Identify top features per factor per view
  top_features <- extract_top_mofa_features(weights, n_top = 50)

  list(
    factors = factors,
    weights = weights,
    variance_explained = var_explained,
    r2_per_factor = r2_per_factor,
    r2_total = r2_total,
    top_features = top_features
  )
}

#' Extract top features per factor from MOFA weights
extract_top_mofa_features <- function(weights, n_top = 50) {
  all_top <- list()

  for (view in names(weights)) {
    w <- weights[[view]]
    n_factors <- ncol(w)

    view_top <- list()
    for (k in seq_len(n_factors)) {
      factor_name <- colnames(w)[k]
      loadings <- w[, k]

      # Get top positive and negative
      ord <- order(abs(loadings), decreasing = TRUE)
      top_idx <- head(ord, n_top)

      view_top[[factor_name]] <- data.frame(
        feature_id = rownames(w)[top_idx],
        weight = loadings[top_idx],
        abs_weight = abs(loadings[top_idx]),
        factor = factor_name,
        view = view,
        stringsAsFactors = FALSE
      )
    }

    all_top[[view]] <- do.call(rbind, view_top)
  }

  do.call(rbind, all_top)
}

#' Create MOFA visualization plots
create_mofa_plots <- function(mofa_model, mofa_results, metadata, config) {
  log_message("Creating MOFA2 plots...")

  condition_col <- config$design$condition_column

  # 1. Variance explained heatmap
  tryCatch({
    p <- plot_mofa_variance_heatmap(mofa_results$r2_per_factor)
    save_plot(p, "mofa_variance_heatmap", config, width = 8, height = 6)
  }, error = function(e) log_message("Failed to create variance heatmap: ", e$message))

  # 2. Factor scatter plots
  tryCatch({
    factors <- mofa_results$factors
    if (ncol(factors) >= 2) {
      p <- plot_mofa_factors(factors, metadata, condition_col, c(1, 2))
      save_plot(p, "mofa_factors_1_2", config, width = 8, height = 6)
    }
    if (ncol(factors) >= 3) {
      p <- plot_mofa_factors(factors, metadata, condition_col, c(1, 3))
      save_plot(p, "mofa_factors_1_3", config, width = 8, height = 6)
    }
  }, error = function(e) log_message("Failed to create factor plots: ", e$message))

  # 3. Top weights per factor
  tryCatch({
    for (view in names(mofa_results$weights)) {
      p <- plot_mofa_top_weights(mofa_results$weights[[view]], view, n_top = 20)
      save_plot(p, paste0("mofa_top_weights_", view), config, width = 10, height = 8)
    }
  }, error = function(e) log_message("Failed to create weight plots: ", e$message))

  log_message("MOFA2 plots complete")
}

#' Plot variance explained heatmap
plot_mofa_variance_heatmap <- function(r2_matrix) {
  # Convert to long format
  df <- as.data.frame(r2_matrix)
  df$view <- rownames(df)
  df_long <- tidyr::pivot_longer(df, -view, names_to = "factor", values_to = "r2")

  ggplot2::ggplot(df_long, ggplot2::aes(x = factor, y = view, fill = r2)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", r2)), size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue", name = "R² (%)") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "MOFA2: Variance Explained per Factor",
      x = "Factor",
      y = "View (Omics)"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

#' Plot MOFA factors
plot_mofa_factors <- function(factors, metadata, condition_col, factor_idx = c(1, 2)) {
  f1 <- factor_idx[1]
  f2 <- factor_idx[2]

  df <- data.frame(
    Factor1 = factors[, f1],
    Factor2 = factors[, f2],
    sample_id = rownames(factors),
    stringsAsFactors = FALSE
  )

  if (!is.null(condition_col) && condition_col %in% colnames(metadata)) {
    df$condition <- metadata[df$sample_id, condition_col]
  } else {
    df$condition <- "All"
  }

  ggplot2::ggplot(df, ggplot2::aes(x = Factor1, y = Factor2, color = condition)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0("MOFA2: Factor ", f1, " vs Factor ", f2),
      x = paste0("Factor ", f1),
      y = paste0("Factor ", f2),
      color = condition_col
    )
}

#' Plot top weights for a view
plot_mofa_top_weights <- function(weights, view_name, n_top = 20, n_factors = 3) {
  n_factors <- min(n_factors, ncol(weights))
  plots <- list()

  for (k in seq_len(n_factors)) {
    w <- weights[, k]
    ord <- order(abs(w), decreasing = TRUE)
    top_idx <- head(ord, n_top)

    df <- data.frame(
      feature = rownames(weights)[top_idx],
      weight = w[top_idx],
      stringsAsFactors = FALSE
    )
    df$feature <- factor(df$feature, levels = df$feature[order(abs(df$weight))])
    df$direction <- ifelse(df$weight > 0, "positive", "negative")

    plots[[k]] <- ggplot2::ggplot(df, ggplot2::aes(x = weight, y = feature, fill = direction)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c("positive" = "firebrick", "negative" = "steelblue")) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste0("Factor ", k),
        x = "Weight",
        y = NULL
      ) +
      ggplot2::theme(legend.position = "none")
  }

  # Combine plots
  if (length(plots) > 1 && requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = n_factors) +
      patchwork::plot_annotation(title = paste0("MOFA2 Top Weights: ", view_name))
    return(combined)
  } else {
    return(plots[[1]])
  }
}

#' Test factor associations with phenotypes
test_mofa_factor_associations <- function(mofa_results, metadata, config) {
  log_message("Testing MOFA factor associations with phenotypes...")

  factors <- mofa_results$factors
  condition_col <- config$design$condition_column

  results <- list()

  # Test condition association
  if (condition_col %in% colnames(metadata)) {
    condition <- metadata[rownames(factors), condition_col]

    for (k in seq_len(ncol(factors))) {
      factor_name <- colnames(factors)[k]
      factor_values <- factors[, k]

      # t-test or ANOVA depending on number of groups
      n_groups <- length(unique(condition))

      if (n_groups == 2) {
        test_result <- t.test(factor_values ~ condition)
        results[[factor_name]] <- data.frame(
          factor = factor_name,
          test = "t-test",
          statistic = test_result$statistic,
          pvalue = test_result$p.value,
          variable = condition_col,
          stringsAsFactors = FALSE
        )
      } else if (n_groups > 2) {
        test_result <- summary(aov(factor_values ~ condition))[[1]]
        results[[factor_name]] <- data.frame(
          factor = factor_name,
          test = "ANOVA",
          statistic = test_result$`F value`[1],
          pvalue = test_result$`Pr(>F)`[1],
          variable = condition_col,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) > 0) {
    result_df <- do.call(rbind, results)
    result_df$padj <- p.adjust(result_df$pvalue, method = "BH")
    save_table(result_df, "mofa_factor_associations.csv", config)
    return(result_df)
  }

  return(NULL)
}
