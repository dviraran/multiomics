# =============================================================================
# Univariate Differential Analysis
# =============================================================================

#' Run differential abundance analysis
#'
#' @param imputed_data Data from imputation step
#' @param exploratory_results Exploratory analysis results
#' @param ingested_data Original ingested data (for annotations)
#' @param config Configuration list
#' @return Differential analysis results
run_differential_analysis <- function(imputed_data, exploratory_results, ingested_data, config) {
  log_message("=== Starting Differential Abundance Analysis ===")

  mat <- imputed_data$matrix
  metadata <- imputed_data$metadata
  sample_roles <- imputed_data$sample_roles
  annotations <- ingested_data$annotations

  # Use only biological samples
  bio_samples <- intersect(sample_roles$biological_samples, colnames(mat))

  if (length(bio_samples) == 0) {
    log_message("No biological samples found. Skipping differential analysis.")
    return(NULL)
  }

  mat_bio <- mat[, bio_samples, drop = FALSE]
  sample_col <- config$input$sample_id_column
  metadata_bio <- metadata[metadata[[sample_col]] %in% bio_samples, ]

  # Build design matrix
  design <- build_design_matrix(metadata_bio, config)

  # Build contrasts
  contrasts <- build_contrast_matrix(design, metadata_bio, config)

  # Run statistical analysis
  de_method <- config$differential$method %||% "limma"

  if (de_method == "limma") {
    de_results <- run_limma(mat_bio, design, contrasts, config)
  } else if (de_method == "t_test") {
    de_results <- run_t_tests(mat_bio, metadata_bio, config)
  } else if (de_method == "wilcoxon") {
    de_results <- run_wilcoxon(mat_bio, metadata_bio, config)
  } else {
    log_message("Unknown method '", de_method, "'. Using limma.")
    de_results <- run_limma(mat_bio, design, contrasts, config)
  }

  # Add annotations
  de_results <- annotate_de_results(de_results, annotations, ingested_data$feature_metadata)

  # Create plots
  plots <- create_de_plots(de_results, mat_bio, metadata_bio, config)

  # Save results
  for (contrast_name in names(de_results)) {
    clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))
    save_table(
      de_results[[contrast_name]]$table,
      paste0("differential_", clean_name, ".csv"),
      config, "tables"
    )
  }

  # Summary
  de_summary <- summarize_de_results(de_results, config)
  save_table(de_summary, "differential_summary.csv", config, "tables")

  log_message("=== Differential Abundance Analysis Complete ===")

  list(
    results = de_results,
    design = design,
    contrasts = contrasts,
    summary = de_summary,
    plots = plots
  )
}

#' Build design matrix
#'
#' @param metadata Sample metadata
#' @param config Configuration
#' @return Design matrix
build_design_matrix <- function(metadata, config) {
  formula_str <- config$design$design_formula
  formula <- as.formula(formula_str)

  formula_vars <- all.vars(formula)

  for (var in formula_vars) {
    if (!var %in% colnames(metadata)) {
      stop("Variable '", var, "' not found in metadata")
    }
    if (!is.factor(metadata[[var]])) {
      metadata[[var]] <- as.factor(metadata[[var]])
    }
  }

  design <- model.matrix(formula, data = metadata)
  colnames(design) <- gsub("\\(Intercept\\)", "Intercept", colnames(design))

  log_message("Design matrix: ", ncol(design), " coefficients")

  return(design)
}

#' Build contrast matrix
#'
#' @param design Design matrix
#' @param metadata Metadata
#' @param config Configuration
#' @return Contrast matrix
build_contrast_matrix <- function(design, metadata, config) {
  contrasts_spec <- config$design$contrasts
  condition_col <- config$design$condition_column

  if (is.null(contrasts_spec) || length(contrasts_spec) == 0) {
    levels <- unique(metadata[[condition_col]])
    contrasts_spec <- list()
    for (i in 1:(length(levels) - 1)) {
      for (j in (i + 1):length(levels)) {
        contrasts_spec[[length(contrasts_spec) + 1]] <- paste(levels[j], "-", levels[i])
      }
    }
    log_message("No contrasts specified. Using pairwise: ", paste(contrasts_spec, collapse = ", "))
  }

  coef_names <- colnames(design)
  contrast_list <- list()

  for (contrast_str in contrasts_spec) {
    parts <- strsplit(trimws(contrast_str), "\\s*-\\s*")[[1]]
    if (length(parts) != 2) {
      warning("Invalid contrast: ", contrast_str)
      next
    }

    numerator <- trimws(parts[1])
    denominator <- trimws(parts[2])

    num_coef <- grep(paste0(condition_col, numerator, "$"), coef_names, value = TRUE)
    den_coef <- grep(paste0(condition_col, denominator, "$"), coef_names, value = TRUE)

    contrast_vec <- rep(0, length(coef_names))
    names(contrast_vec) <- coef_names

    if (length(num_coef) == 1) contrast_vec[num_coef] <- 1
    if (length(den_coef) == 1) contrast_vec[den_coef] <- -1

    if (numerator == config$design$reference_level && length(den_coef) == 1) {
      contrast_vec[den_coef] <- -1
    }

    contrast_list[[contrast_str]] <- contrast_vec
  }

  contrast_matrix <- do.call(cbind, contrast_list)
  log_message("Created ", ncol(contrast_matrix), " contrasts")

  return(contrast_matrix)
}

#' Run limma differential analysis
#'
#' @param mat Matrix
#' @param design Design matrix
#' @param contrasts Contrast matrix
#' @param config Configuration
#' @return DE results list
run_limma <- function(mat, design, contrasts, config) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma package required for differential analysis")
  }

  log_message("Running limma analysis...")

  fit <- limma::lmFit(mat, design)
  fit_contrasts <- limma::contrasts.fit(fit, contrasts)
  fit_eb <- limma::eBayes(fit_contrasts)

  results <- list()

  for (contrast_name in colnames(contrasts)) {
    log_message("Processing contrast: ", contrast_name)

    tt <- limma::topTable(
      fit_eb,
      coef = contrast_name,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )

    tt$feature_id <- rownames(tt)
    colnames(tt)[colnames(tt) == "logFC"] <- "log2FC"
    colnames(tt)[colnames(tt) == "AveExpr"] <- "avg_intensity"

    adj_pval <- config$differential$adj_pvalue_threshold %||% 0.05
    log2fc <- config$differential$log2fc_threshold %||% 1

    tt$significant <- tt$adj.P.Val < adj_pval & abs(tt$log2FC) > log2fc
    tt$direction <- ifelse(tt$log2FC > 0, "up", "down")

    tt <- tt[, c("feature_id", "log2FC", "avg_intensity", "t", "P.Value", "adj.P.Val", "B",
                 "significant", "direction")]

    results[[contrast_name]] <- list(
      table = tt,
      n_sig = sum(tt$significant),
      n_up = sum(tt$significant & tt$direction == "up"),
      n_down = sum(tt$significant & tt$direction == "down")
    )

    log_message("  Significant: ", results[[contrast_name]]$n_sig)
  }

  return(results)
}

#' Run t-tests
#'
#' @param mat Matrix
#' @param metadata Metadata
#' @param config Configuration
#' @return DE results
run_t_tests <- function(mat, metadata, config) {
  log_message("Running t-tests...")

  condition_col <- config$design$condition_column
  sample_col <- config$input$sample_id_column
  conditions <- unique(metadata[[condition_col]])

  if (length(conditions) != 2) {
    log_message("t-test requires exactly 2 conditions. Using limma instead.")
    design <- model.matrix(~ metadata[[condition_col]])
    contrasts <- matrix(c(0, 1), ncol = 1)
    colnames(contrasts) <- paste(conditions[2], "-", conditions[1])
    return(run_limma(mat, design, contrasts, config))
  }

  samples_1 <- metadata[[sample_col]][metadata[[condition_col]] == conditions[1]]
  samples_2 <- metadata[[sample_col]][metadata[[condition_col]] == conditions[2]]

  results_list <- lapply(seq_len(nrow(mat)), function(i) {
    x <- mat[i, samples_1]
    y <- mat[i, samples_2]

    if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2) {
      return(data.frame(log2FC = NA, P.Value = NA))
    }

    test <- t.test(y, x)
    data.frame(
      log2FC = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
      P.Value = test$p.value
    )
  })

  tt <- do.call(rbind, results_list)
  tt$feature_id <- rownames(mat)
  tt$adj.P.Val <- p.adjust(tt$P.Value, method = "BH")
  tt$avg_intensity <- rowMeans(mat, na.rm = TRUE)

  adj_pval <- config$differential$adj_pvalue_threshold %||% 0.05
  log2fc <- config$differential$log2fc_threshold %||% 1

  tt$significant <- tt$adj.P.Val < adj_pval & abs(tt$log2FC) > log2fc
  tt$direction <- ifelse(tt$log2FC > 0, "up", "down")

  contrast_name <- paste(conditions[2], "-", conditions[1])

  list(results = setNames(list(list(
    table = tt,
    n_sig = sum(tt$significant, na.rm = TRUE),
    n_up = sum(tt$significant & tt$direction == "up", na.rm = TRUE),
    n_down = sum(tt$significant & tt$direction == "down", na.rm = TRUE)
  )), contrast_name))
}

#' Run Wilcoxon tests
#'
#' @param mat Matrix
#' @param metadata Metadata
#' @param config Configuration
#' @return DE results
run_wilcoxon <- function(mat, metadata, config) {
  log_message("Running Wilcoxon tests...")

  condition_col <- config$design$condition_column
  sample_col <- config$input$sample_id_column
  conditions <- unique(metadata[[condition_col]])

  if (length(conditions) != 2) {
    log_message("Wilcoxon test requires exactly 2 conditions. Using limma instead.")
    design <- model.matrix(~ metadata[[condition_col]])
    contrasts <- matrix(c(0, 1), ncol = 1)
    colnames(contrasts) <- paste(conditions[2], "-", conditions[1])
    return(run_limma(mat, design, contrasts, config))
  }

  samples_1 <- metadata[[sample_col]][metadata[[condition_col]] == conditions[1]]
  samples_2 <- metadata[[sample_col]][metadata[[condition_col]] == conditions[2]]

  results_list <- lapply(seq_len(nrow(mat)), function(i) {
    x <- mat[i, samples_1]
    y <- mat[i, samples_2]

    if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2) {
      return(data.frame(log2FC = NA, P.Value = NA))
    }

    test <- wilcox.test(y, x)
    data.frame(
      log2FC = median(y, na.rm = TRUE) - median(x, na.rm = TRUE),
      P.Value = test$p.value
    )
  })

  tt <- do.call(rbind, results_list)
  tt$feature_id <- rownames(mat)
  tt$adj.P.Val <- p.adjust(tt$P.Value, method = "BH")
  tt$avg_intensity <- rowMeans(mat, na.rm = TRUE)

  adj_pval <- config$differential$adj_pvalue_threshold %||% 0.05
  log2fc <- config$differential$log2fc_threshold %||% 1

  tt$significant <- tt$adj.P.Val < adj_pval & abs(tt$log2FC) > log2fc
  tt$direction <- ifelse(tt$log2FC > 0, "up", "down")

  contrast_name <- paste(conditions[2], "-", conditions[1])

  list(results = setNames(list(list(
    table = tt,
    n_sig = sum(tt$significant, na.rm = TRUE),
    n_up = sum(tt$significant & tt$direction == "up", na.rm = TRUE),
    n_down = sum(tt$significant & tt$direction == "down", na.rm = TRUE)
  )), contrast_name))
}

#' Add annotations to DE results
#'
#' @param de_results DE results
#' @param annotations Annotation data
#' @param feature_metadata Feature metadata
#' @return Annotated results
annotate_de_results <- function(de_results, annotations, feature_metadata) {
  if (is.null(annotations) && is.null(feature_metadata)) {
    return(de_results)
  }

  # Combine annotation sources
  annotation_cols <- NULL

  if (!is.null(annotations$annotation_table)) {
    annotation_cols <- annotations$annotation_table
  }

  if (!is.null(feature_metadata)) {
    if (is.null(annotation_cols)) {
      annotation_cols <- feature_metadata
    } else {
      annotation_cols <- merge(annotation_cols, feature_metadata,
                               by = "feature_id", all = TRUE)
    }
  }

  if (is.null(annotation_cols)) {
    return(de_results)
  }

  for (contrast_name in names(de_results)) {
    tt <- de_results[[contrast_name]]$table

    match_idx <- match(tt$feature_id, annotation_cols$feature_id)

    for (col in setdiff(colnames(annotation_cols), "feature_id")) {
      if (col %in% c("metabolite_name", "hmdb_id", "kegg_id", "class", "mz", "rt")) {
        tt[[col]] <- annotation_cols[[col]][match_idx]
      }
    }

    de_results[[contrast_name]]$table <- tt
  }

  return(de_results)
}

#' Create DE plots
#'
#' @param de_results DE results
#' @param mat Matrix
#' @param metadata Metadata
#' @param config Configuration
#' @return List of plots
create_de_plots <- function(de_results, mat, metadata, config) {
  all_plots <- list()

  adj_pval <- config$differential$adj_pvalue_threshold %||% 0.05
  log2fc <- config$differential$log2fc_threshold %||% 1

  for (contrast_name in names(de_results)) {
    clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))
    tt <- de_results[[contrast_name]]$table

    plots <- list()

    # Volcano plot
    volcano_df <- tt
    volcano_df$neg_log10_pval <- -log10(volcano_df$P.Value)
    volcano_df$status <- "Not significant"
    volcano_df$status[volcano_df$significant & volcano_df$direction == "up"] <- "Up"
    volcano_df$status[volcano_df$significant & volcano_df$direction == "down"] <- "Down"

    plots$volcano <- ggplot2::ggplot(volcano_df,
      ggplot2::aes(x = log2FC, y = neg_log10_pval, color = status)) +
      ggplot2::geom_point(alpha = 0.5, size = 1) +
      ggplot2::scale_color_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2",
                                              "Not significant" = "grey60")) +
      ggplot2::geom_vline(xintercept = c(-log2fc, log2fc), linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(adj_pval), linetype = "dashed") +
      ggplot2::labs(title = paste("Volcano:", contrast_name),
                    x = "Log2 Fold Change", y = "-Log10(P-value)") +
      ggplot2::theme_minimal()

    # P-value histogram
    plots$pval_hist <- ggplot2::ggplot(tt, ggplot2::aes(x = P.Value)) +
      ggplot2::geom_histogram(bins = 50, fill = "#0072B2", color = "white") +
      ggplot2::labs(title = paste("P-value Distribution:", contrast_name),
                    x = "P-value", y = "Count") +
      ggplot2::theme_minimal()

    # Top features boxplot
    top_features <- head(tt[tt$significant, ], 10)
    if (nrow(top_features) > 0) {
      top_ids <- top_features$feature_id
      top_mat <- mat[top_ids, , drop = FALSE]

      sample_col <- config$input$sample_id_column
      condition_col <- config$design$condition_column

      box_df <- reshape2::melt(top_mat)
      colnames(box_df) <- c("feature", "sample", "intensity")
      box_df$condition <- metadata[[condition_col]][match(box_df$sample, metadata[[sample_col]])]

      plots$top_boxplot <- ggplot2::ggplot(box_df,
        ggplot2::aes(x = condition, y = intensity, fill = condition)) +
        ggplot2::geom_boxplot() +
        ggplot2::facet_wrap(~ feature, scales = "free_y") +
        ggplot2::labs(title = paste("Top Features:", contrast_name)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")

      save_plot(plots$top_boxplot, paste0("top_features_", clean_name, ".png"),
                config, height = 10, subdir = "plots")
    }

    save_plot(plots$volcano, paste0("volcano_", clean_name, ".png"), config, subdir = "plots")
    save_plot(plots$pval_hist, paste0("pvalue_dist_", clean_name, ".png"), config, subdir = "plots")

    all_plots[[contrast_name]] <- plots
  }

  return(all_plots)
}

#' Summarize DE results
#'
#' @param de_results DE results
#' @param config Configuration
#' @return Summary data frame
summarize_de_results <- function(de_results, config) {
  summary_rows <- lapply(names(de_results), function(contrast_name) {
    res <- de_results[[contrast_name]]
    data.frame(
      contrast = contrast_name,
      total_tested = nrow(res$table),
      significant = res$n_sig,
      up_regulated = res$n_up,
      down_regulated = res$n_down,
      pct_significant = round(res$n_sig / nrow(res$table) * 100, 2),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_rows)
}
