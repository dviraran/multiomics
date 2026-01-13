# =============================================================================
# Differential Abundance Analysis
# =============================================================================

#' Run differential abundance analysis
#'
#' @param imputed_data Data from imputation step
#' @param qc_data QC results
#' @param config Configuration list
#' @return Differential analysis results
run_differential_analysis <- function(imputed_data, qc_data, config) {
  log_message("=== Starting Differential Abundance Analysis ===")

  mat <- imputed_data$matrix
  metadata <- imputed_data$metadata
  annotations <- imputed_data$annotations

  # Build design matrix
  design <- build_design_matrix(metadata, config)

  # Build contrast matrix
  contrasts <- build_contrast_matrix(design, metadata, config)

  # Run limma analysis
  de_results <- run_limma(mat, design, contrasts, config)

  # Add annotations to results
  de_results <- annotate_de_results(de_results, annotations)

  # Create plots for each contrast
  plots <- create_de_plots(de_results, config)

  # Save results
  for (contrast_name in names(de_results)) {
    clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))
    save_table(
      de_results[[contrast_name]]$table,
      paste0("differential_", clean_name, ".csv"),
      config,
      "tables"
    )
  }

  # Create summary
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

#' Build design matrix from metadata
#'
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Design matrix
build_design_matrix <- function(metadata, config) {
  log_message("Building design matrix...")

  formula_str <- config$design$design_formula
  formula <- as.formula(formula_str)

  # Ensure all formula terms are factors
  formula_vars <- all.vars(formula)

  for (var in formula_vars) {
    if (!var %in% colnames(metadata)) {
      stop("Variable '", var, "' in design formula not found in metadata")
    }
    if (!is.factor(metadata[[var]])) {
      metadata[[var]] <- as.factor(metadata[[var]])
    }
  }

  # Create design matrix
  design <- model.matrix(formula, data = metadata)

  # Clean up column names
  colnames(design) <- gsub("\\(Intercept\\)", "Intercept", colnames(design))

  log_message("Design matrix: ", ncol(design), " coefficients for ", nrow(design), " samples")

  return(design)
}

#' Build contrast matrix
#'
#' @param design Design matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return Contrast matrix
build_contrast_matrix <- function(design, metadata, config) {
  log_message("Building contrast matrix...")

  contrasts_spec <- config$design$contrasts
  condition_col <- config$design$condition_column

  if (is.null(contrasts_spec) || length(contrasts_spec) == 0) {
    # Default: all pairwise contrasts
    levels <- unique(metadata[[condition_col]])
    contrasts_spec <- list()
    for (i in 1:(length(levels)-1)) {
      for (j in (i+1):length(levels)) {
        contrasts_spec[[length(contrasts_spec) + 1]] <- paste(levels[j], "-", levels[i])
      }
    }
    log_message("No contrasts specified. Using pairwise: ", paste(contrasts_spec, collapse = ", "))
  }

  # Get design column names (excluding intercept)
  coef_names <- colnames(design)

  # Parse contrasts and create matrix
  contrast_list <- list()

  for (contrast_str in contrasts_spec) {
    # Parse contrast string like "treated - control"
    parts <- strsplit(trimws(contrast_str), "\\s*-\\s*")[[1]]

    if (length(parts) != 2) {
      warning("Invalid contrast format: ", contrast_str, ". Expected 'A - B'")
      next
    }

    numerator <- trimws(parts[1])
    denominator <- trimws(parts[2])

    # Find matching coefficient names
    # Look for condition_numerator pattern
    num_coef <- grep(paste0(condition_col, numerator, "$"), coef_names, value = TRUE)
    den_coef <- grep(paste0(condition_col, denominator, "$"), coef_names, value = TRUE)

    # Create contrast vector
    contrast_vec <- rep(0, length(coef_names))
    names(contrast_vec) <- coef_names

    if (length(num_coef) == 1) {
      contrast_vec[num_coef] <- 1
    } else if (numerator == config$design$reference_level) {
      # Reference level is not explicitly in the model, it's the intercept
      # Do nothing for numerator
    } else {
      warning("Could not find coefficient for ", numerator)
    }

    if (length(den_coef) == 1) {
      contrast_vec[den_coef] <- -1
    } else if (denominator == config$design$reference_level) {
      # Reference level, handled by the numerator coefficient
    } else {
      warning("Could not find coefficient for ", denominator)
    }

    # Simplify for reference level comparisons
    if (numerator == config$design$reference_level && length(den_coef) == 1) {
      # Comparing ref - other = -other coefficient
      contrast_vec[den_coef] <- -1
    } else if (denominator == config$design$reference_level && length(num_coef) == 1) {
      # Comparing other - ref = other coefficient (already set)
    }

    contrast_list[[contrast_str]] <- contrast_vec
  }

  if (length(contrast_list) == 0) {
    stop("No valid contrasts could be constructed")
  }

  # Convert to matrix
  contrast_matrix <- do.call(cbind, contrast_list)
  colnames(contrast_matrix) <- names(contrast_list)

  log_message("Created ", ncol(contrast_matrix), " contrasts")

  return(contrast_matrix)
}

#' Run limma differential analysis
#'
#' @param mat Expression matrix (features x samples)
#' @param design Design matrix
#' @param contrasts Contrast matrix
#' @param config Configuration list
#' @return List of results per contrast
run_limma <- function(mat, design, contrasts, config) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma package is required for differential analysis")
  }

  log_message("Running limma analysis...")

  # Fit linear model
  fit <- limma::lmFit(mat, design)

  # Apply contrasts
  fit_contrasts <- limma::contrasts.fit(fit, contrasts)

  # Empirical Bayes moderation
  fit_eb <- limma::eBayes(fit_contrasts)

  # Extract results for each contrast
  results <- list()

  for (contrast_name in colnames(contrasts)) {
    log_message("Processing contrast: ", contrast_name)

    # Get topTable
    tt <- limma::topTable(
      fit_eb,
      coef = contrast_name,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )

    # Add feature IDs
    tt$feature_id <- rownames(tt)

    # Rename columns for clarity
    colnames(tt)[colnames(tt) == "logFC"] <- "log2FC"
    colnames(tt)[colnames(tt) == "AveExpr"] <- "avg_intensity"

    # Mark significance
    adj_pval_thresh <- config$differential$adj_pvalue_threshold %||% 0.05
    log2fc_thresh <- config$differential$log2fc_threshold %||% 1

    tt$significant <- tt$adj.P.Val < adj_pval_thresh & abs(tt$log2FC) > log2fc_thresh
    tt$direction <- ifelse(tt$log2FC > 0, "up", "down")

    # Reorder columns
    tt <- tt[, c("feature_id", "log2FC", "avg_intensity", "t", "P.Value", "adj.P.Val", "B",
                 "significant", "direction")]

    results[[contrast_name]] <- list(
      table = tt,
      n_sig = sum(tt$significant),
      n_up = sum(tt$significant & tt$direction == "up"),
      n_down = sum(tt$significant & tt$direction == "down")
    )

    log_message("  Significant: ", results[[contrast_name]]$n_sig,
                " (", results[[contrast_name]]$n_up, " up, ",
                results[[contrast_name]]$n_down, " down)")
  }

  return(results)
}

#' Add annotations to DE results
#'
#' @param de_results DE results list
#' @param annotations Annotation data frame
#' @return Annotated DE results
annotate_de_results <- function(de_results, annotations) {
  for (contrast_name in names(de_results)) {
    tt <- de_results[[contrast_name]]$table

    # Match annotations
    match_idx <- match(tt$feature_id, annotations$feature_id)

    # Add annotation columns if available
    if ("gene_symbol" %in% colnames(annotations)) {
      tt$gene_symbol <- annotations$gene_symbol[match_idx]
    }
    if ("protein_name" %in% colnames(annotations)) {
      tt$protein_name <- annotations$protein_name[match_idx]
    }
    if ("uniprot_accession" %in% colnames(annotations)) {
      tt$uniprot_accession <- annotations$uniprot_accession[match_idx]
    }

    de_results[[contrast_name]]$table <- tt
  }

  return(de_results)
}

#' Create DE visualization plots
#'
#' @param de_results DE results list
#' @param config Configuration list
#' @return List of plots
create_de_plots <- function(de_results, config) {
  all_plots <- list()

  adj_pval_thresh <- config$differential$adj_pvalue_threshold %||% 0.05
  log2fc_thresh <- config$differential$log2fc_threshold %||% 1

  for (contrast_name in names(de_results)) {
    clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))
    tt <- de_results[[contrast_name]]$table

    plots <- list()

    # 1. Volcano plot
    volcano_df <- tt
    volcano_df$neg_log10_pval <- -log10(volcano_df$P.Value)
    volcano_df$status <- "Not significant"
    volcano_df$status[volcano_df$significant & volcano_df$direction == "up"] <- "Up-regulated"
    volcano_df$status[volcano_df$significant & volcano_df$direction == "down"] <- "Down-regulated"
    volcano_df$status <- factor(volcano_df$status,
                                levels = c("Up-regulated", "Down-regulated", "Not significant"))

    # Top genes for labeling
    top_genes <- head(tt[tt$significant, ], 10)

    plots$volcano <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = log2FC, y = neg_log10_pval, color = status)) +
      ggplot2::geom_point(alpha = 0.5, size = 1) +
      ggplot2::scale_color_manual(values = c("Up-regulated" = "#D55E00",
                                              "Down-regulated" = "#0072B2",
                                              "Not significant" = "grey60")) +
      ggplot2::geom_vline(xintercept = c(-log2fc_thresh, log2fc_thresh), linetype = "dashed", color = "grey40") +
      ggplot2::geom_hline(yintercept = -log10(adj_pval_thresh), linetype = "dashed", color = "grey40") +
      ggplot2::labs(
        title = paste("Volcano Plot:", contrast_name),
        x = "Log2 Fold Change",
        y = "-Log10(P-value)",
        color = "Status"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")

    # Add labels for top genes if gene_symbol available
    if ("gene_symbol" %in% colnames(top_genes) && sum(!is.na(top_genes$gene_symbol)) > 0) {
      label_df <- volcano_df[volcano_df$feature_id %in% top_genes$feature_id, ]
      label_df$label <- label_df$gene_symbol
      label_df$label[is.na(label_df$label)] <- label_df$feature_id[is.na(label_df$label)]

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        plots$volcano <- plots$volcano +
          ggrepel::geom_text_repel(
            data = label_df,
            ggplot2::aes(label = label),
            size = 3,
            max.overlaps = 20
          )
      }
    }

    # 2. MA plot
    plots$ma <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = avg_intensity, y = log2FC, color = status)) +
      ggplot2::geom_point(alpha = 0.5, size = 1) +
      ggplot2::scale_color_manual(values = c("Up-regulated" = "#D55E00",
                                              "Down-regulated" = "#0072B2",
                                              "Not significant" = "grey60")) +
      ggplot2::geom_hline(yintercept = 0, color = "black") +
      ggplot2::geom_hline(yintercept = c(-log2fc_thresh, log2fc_thresh), linetype = "dashed", color = "grey40") +
      ggplot2::labs(
        title = paste("MA Plot:", contrast_name),
        x = "Average Log2 Intensity",
        y = "Log2 Fold Change",
        color = "Status"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")

    # 3. P-value distribution
    plots$pval_hist <- ggplot2::ggplot(tt, ggplot2::aes(x = P.Value)) +
      ggplot2::geom_histogram(bins = 50, fill = "#0072B2", color = "white") +
      ggplot2::labs(
        title = paste("P-value Distribution:", contrast_name),
        x = "P-value",
        y = "Count"
      ) +
      ggplot2::theme_minimal()

    # Save plots
    save_plot(plots$volcano, paste0("volcano_", clean_name, ".png"), config, subdir = "plots")
    save_plot(plots$ma, paste0("ma_plot_", clean_name, ".png"), config, subdir = "plots")
    save_plot(plots$pval_hist, paste0("pvalue_dist_", clean_name, ".png"), config, subdir = "plots")

    all_plots[[contrast_name]] <- plots
  }

  return(all_plots)
}

#' Summarize DE results across all contrasts
#'
#' @param de_results DE results list
#' @param config Configuration list
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
