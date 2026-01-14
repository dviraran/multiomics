# =============================================================================
# Batch Effect Detection and Correction
# =============================================================================
# This script detects and corrects batch effects in RNA-seq data.
#
# Detection methods:
#   - PCA visualization by batch variables
#   - PVCA (Principal Variance Component Analysis)
#   - Silhouette analysis comparing batch vs condition clustering
#
# Correction methods:
#   - ComBat-Seq (count-preserving, recommended for RNA-seq)
#   - SVA (Surrogate Variable Analysis)
#   - RUVSeq (Remove Unwanted Variation)
# =============================================================================

#' Main function: Detect and optionally correct batch effects
#' @param vst_counts VST-transformed counts for detection
#' @param counts_filtered Filtered raw counts (for correction)
#' @param metadata Sample metadata
#' @param config Pipeline configuration
#' @return List with batch detection results and optionally corrected counts
run_batch_analysis <- function(vst_counts, counts_filtered, metadata, config) {
  log_message("=== Running Batch Effect Analysis ===")

  bc <- get_batch_config(config)

  if (!bc$detect_batch) {
    log_message("Batch detection disabled in config. Skipping.")
    return(NULL)
  }

  results <- list()

  # 1. Identify potential batch variables
  batch_cols <- identify_batch_columns(metadata, bc, config)

  if (length(batch_cols) == 0) {
    log_message("No batch variables identified. Skipping batch analysis.")
    return(NULL)
  }

  log_message("Analyzing batch variables: ", paste(batch_cols, collapse = ", "))

  # 2. Detect batch effects
  results$detection <- detect_batch_effects(
    vst_counts, metadata, batch_cols, config
  )

  # 3. Quantify batch effects (PVCA)
  results$pvca <- run_pvca_analysis(
    vst_counts, metadata, batch_cols, config
  )

  # 4. Correct batch effects if requested and significant
  if (bc$correct_batch && results$detection$batch_detected) {
    results$correction <- correct_batch_effects(
      counts_filtered, metadata, batch_cols, bc, config
    )
  }

  # 5. Generate visualizations
  plot_batch_analysis(results, metadata, batch_cols, config)

  # 6. Summary
  results$summary <- summarize_batch_results(results, config)

  log_message("=== Batch Analysis Complete ===")
  return(results)
}

#' Get batch correction configuration with defaults
get_batch_config <- function(config) {
  bc <- config$batch_correction %||% list()

  list(
    detect_batch = bc$detect_batch %||% TRUE,
    batch_columns = bc$batch_columns,  # NULL = auto-detect
    correct_batch = bc$correct_batch %||% FALSE,
    correct_method = bc$correct_method %||% "combat_seq",
    n_surrogate_variables = bc$n_surrogate_variables %||% "auto",
    batch_threshold = bc$batch_threshold %||% 0.1  # variance explained threshold
  )
}

#' Identify potential batch columns in metadata
identify_batch_columns <- function(metadata, bc, config) {
  # If explicitly specified, use those
  if (!is.null(bc$batch_columns)) {
    valid_cols <- intersect(bc$batch_columns, colnames(metadata))
    return(valid_cols)
  }

  # Auto-detect potential batch variables
  condition_col <- config$group_col %||% config$condition_column %||% "condition"
  sample_col <- config$sample_id_col %||% "sample_id"

  # Exclude condition and sample ID columns
  exclude_cols <- c(condition_col, sample_col, "sample")

  # Look for common batch-related column names
  batch_patterns <- c("batch", "run", "date", "plate", "lane", "flowcell",
                      "library", "sequencing", "cohort", "site", "center")

  potential_batch <- c()

  for (col in colnames(metadata)) {
    if (col %in% exclude_cols) next

    # Check if column name matches batch patterns
    col_lower <- tolower(col)
    if (any(sapply(batch_patterns, function(p) grepl(p, col_lower)))) {
      potential_batch <- c(potential_batch, col)
      next
    }

    # Check if categorical with reasonable number of levels
    values <- metadata[[col]]
    if (is.factor(values) || is.character(values)) {
      n_levels <- length(unique(values))
      n_samples <- nrow(metadata)

      # Batch variables typically have 2-20 levels, not 1 per sample
      if (n_levels >= 2 && n_levels <= min(20, n_samples / 2)) {
        # Avoid variables that are highly correlated with condition
        if (!col %in% potential_batch) {
          potential_batch <- c(potential_batch, col)
        }
      }
    }
  }

  log_message("Auto-detected potential batch variables: ",
              paste(potential_batch, collapse = ", "))

  return(potential_batch)
}

# =============================================================================
# Batch Effect Detection
# =============================================================================

#' Detect batch effects using multiple approaches
#' @param vst_counts VST-transformed counts
#' @param metadata Sample metadata
#' @param batch_cols Columns to test for batch effects
#' @param config Pipeline configuration
detect_batch_effects <- function(vst_counts, metadata, batch_cols, config) {
  log_message("Detecting batch effects...")

  results <- list()

  # Get condition column
  condition_col <- config$group_col %||% config$condition_column %||% "condition"

  # Run PCA
  pca_result <- prcomp(t(vst_counts), scale. = TRUE, center = TRUE)
  pca_scores <- pca_result$x[, 1:min(10, ncol(pca_result$x))]
  var_explained <- summary(pca_result)$importance[2, 1:min(10, ncol(pca_result$x))]

  results$pca <- list(
    scores = pca_scores,
    var_explained = var_explained
  )

  # For each batch variable, test association with PCs
  batch_stats <- list()

  for (batch_col in batch_cols) {
    if (!batch_col %in% colnames(metadata)) next

    batch_var <- metadata[[batch_col]]
    if (is.null(batch_var)) next

    # Ensure samples are aligned
    common_samples <- intersect(rownames(pca_scores), rownames(metadata))
    batch_var <- metadata[common_samples, batch_col]

    # Test association with top PCs using ANOVA
    pc_associations <- sapply(1:min(5, ncol(pca_scores)), function(pc) {
      pc_values <- pca_scores[common_samples, pc]
      if (length(unique(batch_var)) < 2) return(NA)

      tryCatch({
        aov_result <- aov(pc_values ~ as.factor(batch_var))
        summary(aov_result)[[1]][1, "Pr(>F)"]
      }, error = function(e) NA)
    })

    # R-squared for each PC
    pc_r2 <- sapply(1:min(5, ncol(pca_scores)), function(pc) {
      pc_values <- pca_scores[common_samples, pc]
      if (length(unique(batch_var)) < 2) return(NA)

      tryCatch({
        aov_result <- aov(pc_values ~ as.factor(batch_var))
        ss <- summary(aov_result)[[1]]
        ss[1, "Sum Sq"] / sum(ss[, "Sum Sq"])
      }, error = function(e) NA)
    })

    # Silhouette analysis
    sil_score <- compute_silhouette_score(pca_scores[common_samples, 1:min(5, ncol(pca_scores))],
                                           batch_var)

    # Weighted batch effect (R2 weighted by variance explained)
    weighted_r2 <- sum(pc_r2 * var_explained[1:length(pc_r2)], na.rm = TRUE)

    batch_stats[[batch_col]] <- list(
      pc_pvalues = pc_associations,
      pc_r2 = pc_r2,
      silhouette = sil_score,
      weighted_r2 = weighted_r2,
      n_levels = length(unique(batch_var)),
      significant_pcs = sum(pc_associations < 0.05, na.rm = TRUE)
    )

    log_message("  ", batch_col, ": weighted R2=", round(weighted_r2, 3),
                ", silhouette=", round(sil_score, 3),
                ", significant PCs=", batch_stats[[batch_col]]$significant_pcs)
  }

  # Compare with condition
  if (condition_col %in% colnames(metadata)) {
    cond_var <- metadata[common_samples, condition_col]

    cond_r2 <- sapply(1:min(5, ncol(pca_scores)), function(pc) {
      pc_values <- pca_scores[common_samples, pc]
      tryCatch({
        aov_result <- aov(pc_values ~ as.factor(cond_var))
        ss <- summary(aov_result)[[1]]
        ss[1, "Sum Sq"] / sum(ss[, "Sum Sq"])
      }, error = function(e) NA)
    })

    cond_weighted_r2 <- sum(cond_r2 * var_explained[1:length(cond_r2)], na.rm = TRUE)
    cond_silhouette <- compute_silhouette_score(
      pca_scores[common_samples, 1:min(5, ncol(pca_scores))], cond_var
    )

    results$condition_stats <- list(
      pc_r2 = cond_r2,
      weighted_r2 = cond_weighted_r2,
      silhouette = cond_silhouette
    )

    log_message("  Condition: weighted R2=", round(cond_weighted_r2, 3),
                ", silhouette=", round(cond_silhouette, 3))
  }

  results$batch_stats <- batch_stats

  # Determine if significant batch effects exist
  bc <- get_batch_config(config)
  max_batch_r2 <- max(sapply(batch_stats, function(x) x$weighted_r2), na.rm = TRUE)
  results$batch_detected <- max_batch_r2 > bc$batch_threshold

  if (results$batch_detected) {
    worst_batch <- names(which.max(sapply(batch_stats, function(x) x$weighted_r2)))
    log_message("Batch effects detected. Strongest batch variable: ", worst_batch,
                " (R2=", round(max_batch_r2, 3), ")")
  } else {
    log_message("No significant batch effects detected (threshold: ", bc$batch_threshold, ")")
  }

  # Save detection results
  detection_df <- do.call(rbind, lapply(names(batch_stats), function(bn) {
    bs <- batch_stats[[bn]]
    data.frame(
      batch_variable = bn,
      n_levels = bs$n_levels,
      weighted_r2 = round(bs$weighted_r2, 4),
      silhouette_score = round(bs$silhouette, 4),
      significant_pcs = bs$significant_pcs,
      stringsAsFactors = FALSE
    )
  }))

  output_dir <- config$output_dir %||% "outputs"
  readr::write_csv(detection_df, file.path(output_dir, "batch_effect_detection.csv"))

  results
}

#' Compute silhouette score for clustering
compute_silhouette_score <- function(data_matrix, labels) {
  if (length(unique(labels)) < 2) return(NA)

  tryCatch({
    if (!requireNamespace("cluster", quietly = TRUE)) {
      return(NA)
    }

    dist_mat <- dist(data_matrix)
    sil <- cluster::silhouette(as.numeric(as.factor(labels)), dist_mat)
    mean(sil[, "sil_width"])
  }, error = function(e) NA)
}

# =============================================================================
# PVCA Analysis
# =============================================================================

#' Run Principal Variance Component Analysis
#' @param vst_counts VST-transformed counts
#' @param metadata Sample metadata
#' @param batch_cols Batch columns to analyze
#' @param config Pipeline configuration
run_pvca_analysis <- function(vst_counts, metadata, batch_cols, config) {
  log_message("Running PVCA analysis...")

  # Check if pvca package is available
  if (!requireNamespace("pvca", quietly = TRUE)) {
    log_message("  pvca package not available. Using simplified variance decomposition.")
    return(run_simplified_pvca(vst_counts, metadata, batch_cols, config))
  }

  condition_col <- config$group_col %||% config$condition_column %||% "condition"

  # Prepare metadata for pvca
  all_factors <- c(batch_cols, condition_col)
  all_factors <- intersect(all_factors, colnames(metadata))

  if (length(all_factors) < 2) {
    log_message("  Need at least 2 factors for PVCA. Skipping.")
    return(NULL)
  }

  # Create ExpressionSet
  common_samples <- intersect(colnames(vst_counts), rownames(metadata))
  expr_mat <- vst_counts[, common_samples]
  pdata <- metadata[common_samples, all_factors, drop = FALSE]

  # Convert to factors
  for (col in all_factors) {
    pdata[[col]] <- as.factor(pdata[[col]])
  }

  # Run PVCA
  pvca_result <- tryCatch({
    eset <- Biobase::ExpressionSet(
      assayData = as.matrix(expr_mat),
      phenoData = Biobase::AnnotatedDataFrame(pdata)
    )

    pvca::pvcaBatchAssess(
      abatch = eset,
      batch.factors = all_factors,
      threshold = 0.6
    )
  }, error = function(e) {
    log_message("  PVCA failed: ", e$message)
    NULL
  })

  if (is.null(pvca_result)) {
    return(run_simplified_pvca(vst_counts, metadata, batch_cols, config))
  }

  # Extract results
  var_explained <- pvca_result$dat
  factor_names <- pvca_result$label

  pvca_df <- data.frame(
    factor = factor_names,
    variance_explained = var_explained,
    stringsAsFactors = FALSE
  )

  # Save results
  output_dir <- config$output_dir %||% "outputs"
  readr::write_csv(pvca_df, file.path(output_dir, "pvca_results.csv"))

  log_message("  PVCA variance decomposition:")
  for (i in seq_along(factor_names)) {
    log_message("    ", factor_names[i], ": ", round(var_explained[i] * 100, 1), "%")
  }

  list(
    variance_explained = var_explained,
    factors = factor_names,
    pvca_df = pvca_df
  )
}

#' Simplified PVCA using ANOVA
run_simplified_pvca <- function(vst_counts, metadata, batch_cols, config) {
  log_message("  Running simplified variance decomposition...")

  condition_col <- config$group_col %||% config$condition_column %||% "condition"

  # Run PCA
  pca <- prcomp(t(vst_counts), scale. = TRUE, center = TRUE)
  n_pcs <- min(10, ncol(pca$x))
  pca_scores <- pca$x[, 1:n_pcs]
  var_explained <- summary(pca)$importance[2, 1:n_pcs]

  # Common samples
  common_samples <- intersect(rownames(pca_scores), rownames(metadata))
  pca_sub <- pca_scores[common_samples, ]
  meta_sub <- metadata[common_samples, ]

  # Factors to analyze
  all_factors <- unique(c(batch_cols, condition_col))
  all_factors <- intersect(all_factors, colnames(meta_sub))

  # For each factor, compute variance explained in each PC
  factor_var <- list()

  for (factor_name in all_factors) {
    factor_values <- meta_sub[[factor_name]]
    if (length(unique(factor_values)) < 2) next

    pc_r2 <- sapply(1:n_pcs, function(pc) {
      tryCatch({
        fit <- lm(pca_sub[, pc] ~ as.factor(factor_values))
        summary(fit)$r.squared
      }, error = function(e) 0)
    })

    # Weight by PC variance explained
    weighted_r2 <- sum(pc_r2 * var_explained)
    factor_var[[factor_name]] <- weighted_r2
  }

  # Residual variance
  total_explained <- sum(unlist(factor_var))
  residual <- max(0, 1 - total_explained)
  factor_var$residual <- residual

  # Create results
  pvca_df <- data.frame(
    factor = names(factor_var),
    variance_explained = unlist(factor_var),
    stringsAsFactors = FALSE
  )
  pvca_df <- pvca_df[order(pvca_df$variance_explained, decreasing = TRUE), ]

  # Save results
  output_dir <- config$output_dir %||% "outputs"
  readr::write_csv(pvca_df, file.path(output_dir, "pvca_results.csv"))

  log_message("  Simplified variance decomposition:")
  for (i in seq_len(nrow(pvca_df))) {
    log_message("    ", pvca_df$factor[i], ": ",
                round(pvca_df$variance_explained[i] * 100, 1), "%")
  }

  list(
    variance_explained = pvca_df$variance_explained,
    factors = pvca_df$factor,
    pvca_df = pvca_df
  )
}

# =============================================================================
# Batch Effect Correction
# =============================================================================

#' Correct batch effects
#' @param counts_filtered Filtered raw counts
#' @param metadata Sample metadata
#' @param batch_cols Batch columns
#' @param bc Batch config
#' @param config Full config
correct_batch_effects <- function(counts_filtered, metadata, batch_cols, bc, config) {
  log_message("Correcting batch effects using: ", bc$correct_method)

  # Determine which batch variable to correct (strongest effect)
  # For simplicity, use first batch column or strongest from detection

  batch_col <- batch_cols[1]
  batch_var <- metadata[[batch_col]]

  # Ensure alignment
  common_samples <- intersect(colnames(counts_filtered), rownames(metadata))
  counts_sub <- counts_filtered[, common_samples]
  batch_var <- metadata[common_samples, batch_col]

  condition_col <- config$group_col %||% config$condition_column %||% "condition"
  condition_var <- metadata[common_samples, condition_col]

  results <- list(
    method = bc$correct_method,
    batch_column = batch_col
  )

  if (bc$correct_method == "combat_seq") {
    results$corrected <- correct_combat_seq(counts_sub, batch_var, condition_var)
  } else if (bc$correct_method == "sva") {
    results$corrected <- correct_sva(counts_sub, batch_var, condition_var, bc)
  } else if (bc$correct_method == "ruv") {
    results$corrected <- correct_ruv(counts_sub, batch_var, condition_var, bc)
  } else {
    log_message("Unknown correction method: ", bc$correct_method)
    return(NULL)
  }

  if (!is.null(results$corrected)) {
    # Save corrected counts
    output_dir <- config$output_dir %||% "outputs"
    corrected_df <- as.data.frame(results$corrected)
    corrected_df <- tibble::rownames_to_column(corrected_df, "gene_id")
    readr::write_csv(corrected_df, file.path(output_dir, "batch_corrected_counts.csv"))
    log_message("Saved batch-corrected counts")
  }

  results
}

#' Correct using ComBat-Seq
correct_combat_seq <- function(counts, batch_var, condition_var) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    log_message("  sva package not available for ComBat-Seq")
    return(NULL)
  }

  log_message("  Running ComBat-Seq...")

  corrected <- tryCatch({
    sva::ComBat_seq(
      counts = as.matrix(counts),
      batch = as.factor(batch_var),
      group = as.factor(condition_var)
    )
  }, error = function(e) {
    log_message("  ComBat-Seq failed: ", e$message)
    NULL
  })

  if (!is.null(corrected)) {
    log_message("  ComBat-Seq correction complete")
  }

  corrected
}

#' Correct using SVA
correct_sva <- function(counts, batch_var, condition_var, bc) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    log_message("  sva package not available")
    return(NULL)
  }

  log_message("  Running SVA...")

  # Normalize counts for SVA
  counts_norm <- log2(counts + 1)

  # Create design matrices
  mod <- model.matrix(~as.factor(condition_var))
  mod0 <- model.matrix(~1, data = data.frame(condition = condition_var))

  # Estimate surrogate variables
  n_sv <- bc$n_surrogate_variables
  if (n_sv == "auto") {
    n_sv <- tryCatch({
      sva::num.sv(counts_norm, mod, method = "be")
    }, error = function(e) 2)
  }

  log_message("  Estimating ", n_sv, " surrogate variables")

  svobj <- tryCatch({
    sva::sva(as.matrix(counts_norm), mod, mod0, n.sv = n_sv)
  }, error = function(e) {
    log_message("  SVA failed: ", e$message)
    NULL
  })

  if (is.null(svobj)) return(NULL)

  # Remove SV effects
  corrected <- tryCatch({
    sva::removeBatchEffect(counts_norm, covariates = svobj$sv)
  }, error = function(e) {
    log_message("  Batch effect removal failed: ", e$message)
    NULL
  })

  if (!is.null(corrected)) {
    log_message("  SVA correction complete (", n_sv, " SVs removed)")
  }

  # Note: Returns log-transformed corrected values
  corrected
}

#' Correct using RUVSeq
correct_ruv <- function(counts, batch_var, condition_var, bc) {
  if (!requireNamespace("RUVSeq", quietly = TRUE)) {
    log_message("  RUVSeq package not available")
    return(NULL)
  }

  log_message("  Running RUVSeq...")

  # Create SeqExpressionSet
  eset <- tryCatch({
    RUVSeq::newSeqExpressionSet(
      as.matrix(counts),
      phenoData = data.frame(condition = condition_var, row.names = colnames(counts))
    )
  }, error = function(e) {
    log_message("  Failed to create SeqExpressionSet: ", e$message)
    NULL
  })

  if (is.null(eset)) return(NULL)

  # Use empirical control genes (least variable)
  gene_vars <- apply(counts, 1, var)
  control_genes <- names(sort(gene_vars))[1:min(1000, length(gene_vars))]

  # Run RUVg
  n_k <- min(bc$n_surrogate_variables %||% 2, 5)

  corrected <- tryCatch({
    ruv_result <- RUVSeq::RUVg(eset, control_genes, k = n_k)
    RUVSeq::normCounts(ruv_result)
  }, error = function(e) {
    log_message("  RUVSeq failed: ", e$message)
    NULL
  })

  if (!is.null(corrected)) {
    log_message("  RUVSeq correction complete (k=", n_k, ")")
  }

  corrected
}

# =============================================================================
# Visualizations
# =============================================================================

#' Generate batch effect visualizations
plot_batch_analysis <- function(results, metadata, batch_cols, config) {
  log_message("Generating batch effect plots...")

  output_dir <- config$output_dir %||% "outputs"
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  condition_col <- config$group_col %||% config$condition_column %||% "condition"

  # Plot 1: PCA by batch variables
  if (!is.null(results$detection$pca)) {
    pca_scores <- results$detection$pca$scores
    var_exp <- results$detection$pca$var_explained

    common_samples <- intersect(rownames(pca_scores), rownames(metadata))
    pca_df <- as.data.frame(pca_scores[common_samples, 1:min(5, ncol(pca_scores))])
    pca_df$sample_id <- common_samples

    # Add metadata
    for (col in c(batch_cols, condition_col)) {
      if (col %in% colnames(metadata)) {
        pca_df[[col]] <- metadata[common_samples, col]
      }
    }

    # PCA plot for each batch variable
    for (batch_col in batch_cols) {
      if (!batch_col %in% colnames(pca_df)) next

      p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2,
                                                 color = .data[[batch_col]])) +
        ggplot2::geom_point(size = 3, alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste0("PCA colored by ", batch_col),
          subtitle = paste0("PC1: ", round(var_exp[1] * 100, 1), "%, PC2: ",
                           round(var_exp[2] * 100, 1), "%"),
          x = paste0("PC1 (", round(var_exp[1] * 100, 1), "%)"),
          y = paste0("PC2 (", round(var_exp[2] * 100, 1), "%)"),
          color = batch_col
        )

      ggplot2::ggsave(file.path(plots_dir, paste0("pca_by_", batch_col, ".png")),
                      plot = p, width = 10, height = 8, dpi = 300)
    }

    # PCA by condition
    if (condition_col %in% colnames(pca_df)) {
      p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2,
                                                 color = .data[[condition_col]])) +
        ggplot2::geom_point(size = 3, alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "PCA colored by Condition",
          subtitle = paste0("PC1: ", round(var_exp[1] * 100, 1), "%, PC2: ",
                           round(var_exp[2] * 100, 1), "%"),
          x = paste0("PC1 (", round(var_exp[1] * 100, 1), "%)"),
          y = paste0("PC2 (", round(var_exp[2] * 100, 1), "%)"),
          color = "Condition"
        )

      ggplot2::ggsave(file.path(plots_dir, "pca_by_condition.png"),
                      plot = p, width = 10, height = 8, dpi = 300)
    }
  }

  # Plot 2: PVCA bar plot
  if (!is.null(results$pvca$pvca_df)) {
    pvca_df <- results$pvca$pvca_df

    p <- ggplot2::ggplot(pvca_df, ggplot2::aes(x = reorder(factor, variance_explained),
                                                y = variance_explained * 100)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Principal Variance Component Analysis",
        subtitle = "Variance explained by each factor",
        x = "Factor",
        y = "Variance Explained (%)"
      ) +
      ggplot2::geom_hline(yintercept = 10, linetype = "dashed", color = "red", alpha = 0.7)

    ggplot2::ggsave(file.path(plots_dir, "pvca_variance_decomposition.png"),
                    plot = p, width = 8, height = 6, dpi = 300)
  }

  # Plot 3: Batch effect comparison
  if (!is.null(results$detection$batch_stats)) {
    batch_stats <- results$detection$batch_stats

    comparison_df <- do.call(rbind, lapply(names(batch_stats), function(bn) {
      data.frame(
        variable = bn,
        type = "Batch",
        r_squared = batch_stats[[bn]]$weighted_r2,
        stringsAsFactors = FALSE
      )
    }))

    if (!is.null(results$detection$condition_stats)) {
      comparison_df <- rbind(comparison_df, data.frame(
        variable = "Condition",
        type = "Condition",
        r_squared = results$detection$condition_stats$weighted_r2,
        stringsAsFactors = FALSE
      ))
    }

    p <- ggplot2::ggplot(comparison_df, ggplot2::aes(x = reorder(variable, r_squared),
                                                      y = r_squared * 100,
                                                      fill = type)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Batch vs Condition Effect Comparison",
        subtitle = "Weighted R-squared across top PCs",
        x = "Variable",
        y = "Variance Explained (%)",
        fill = "Type"
      ) +
      ggplot2::scale_fill_manual(values = c("Batch" = "#d73027", "Condition" = "#4575b4"))

    ggplot2::ggsave(file.path(plots_dir, "batch_vs_condition_comparison.png"),
                    plot = p, width = 8, height = 6, dpi = 300)
  }

  log_message("Batch effect plots saved")
}

#' Summarize batch analysis results
summarize_batch_results <- function(results, config) {
  summary_list <- list()

  # Detection summary
  if (!is.null(results$detection)) {
    summary_list$batch_detected <- results$detection$batch_detected

    if (!is.null(results$detection$batch_stats)) {
      max_r2 <- max(sapply(results$detection$batch_stats, function(x) x$weighted_r2))
      worst_batch <- names(which.max(sapply(results$detection$batch_stats,
                                             function(x) x$weighted_r2)))
      summary_list$strongest_batch <- worst_batch
      summary_list$strongest_batch_r2 <- round(max_r2, 4)
    }

    if (!is.null(results$detection$condition_stats)) {
      summary_list$condition_r2 <- round(results$detection$condition_stats$weighted_r2, 4)
    }
  }

  # PVCA summary
  if (!is.null(results$pvca)) {
    summary_list$pvca_factors <- paste(results$pvca$factors, collapse = ", ")
  }

  # Correction summary
  if (!is.null(results$correction)) {
    summary_list$correction_method <- results$correction$method
    summary_list$corrected <- !is.null(results$correction$corrected)
  }

  # Save summary
  summary_df <- data.frame(
    metric = names(summary_list),
    value = as.character(unlist(summary_list)),
    stringsAsFactors = FALSE
  )

  output_dir <- config$output_dir %||% "outputs"
  readr::write_csv(summary_df, file.path(output_dir, "batch_analysis_summary.csv"))

  summary_list
}

# =============================================================================
# Helper function for logging
# =============================================================================

if (!exists("log_message")) {
  log_message <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
    message(msg)
  }
}

if (!exists("%||%")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}
