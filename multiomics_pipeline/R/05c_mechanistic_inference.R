# =============================================================================
# Mechanistic Inference and Regulatory Network Analysis
# =============================================================================
# This script provides deep mechanistic insights from multi-omics data:
#   1. RNA-Protein regulatory analysis (translation efficiency, protein stability)
#   2. Transcription factor activity inference
#   3. Gene regulatory network inference
#   4. Causal mediation analysis (RNA -> Protein -> Phenotype)
# =============================================================================

#' Main function: Run mechanistic analysis
#' @param mae_data MultiAssayExperiment data
#' @param foundational_results Results from foundational analysis
#' @param config Pipeline configuration
#' @return List containing mechanistic analysis results
run_mechanistic_analysis <- function(mae_data, foundational_results, config) {
  log_message("=== Running Mechanistic Inference Analysis ===")

  # Get mechanistic config with defaults
  mc <- get_mechanistic_config(config)

  if (!mc$run_mechanistic) {
    log_message("Mechanistic analysis disabled in config. Skipping.")
    return(NULL)
  }

  results <- list()
  harmonized <- mae_data$harmonized_omics

  # 1. RNA-Protein regulatory analysis
  if ("transcriptomics" %in% names(harmonized) &&
      "proteomics" %in% names(harmonized)) {
    results$rna_protein <- analyze_rna_protein_regulation(
      harmonized, mae_data$gene_mapping, mc, config
    )
  }

  # 2. Transcription factor activity inference
  if ("transcriptomics" %in% names(harmonized) && mc$tf_activity) {
    results$tf_activity <- infer_tf_activity(
      harmonized$transcriptomics, mc, config
    )
  }

  # 3. Gene regulatory network inference
  if (mc$infer_regulatory_network) {
    results$regulatory_network <- infer_regulatory_network(
      harmonized, mae_data$gene_mapping, mc, config
    )
  }

  # 4. Causal mediation analysis
  if (mc$run_mediation && !is.null(mae_data$metadata)) {
    results$mediation <- run_mediation_analysis(
      harmonized, mae_data$metadata, mae_data$gene_mapping, mc, config
    )
  }

  # Generate summary
  results$summary <- summarize_mechanistic_results(results, config)

  log_message("=== Mechanistic Analysis Complete ===")
  return(results)
}

#' Get mechanistic analysis configuration with defaults
get_mechanistic_config <- function(config) {
  mc <- config$mechanistic %||% list()

  list(
    run_mechanistic = mc$run_mechanistic %||% TRUE,
    tf_activity = mc$tf_activity %||% TRUE,
    tf_database = mc$tf_database %||% "dorothea",
    infer_regulatory_network = mc$infer_regulatory_network %||% TRUE,
    network_method = mc$network_method %||% "genie3",
    run_mediation = mc$run_mediation %||% TRUE,
    n_top_regulators = mc$n_top_regulators %||% 50,
    fdr_threshold = mc$fdr_threshold %||% 0.05,
    min_targets = mc$min_targets %||% 10
  )
}

# =============================================================================
# 1. RNA-Protein Regulatory Analysis
# =============================================================================

#' Analyze RNA-Protein regulatory relationships
#' @param harmonized List of harmonized omics data
#' @param gene_mapping Gene mapping table
#' @param mc Mechanistic config
#' @param config Full config
analyze_rna_protein_regulation <- function(harmonized, gene_mapping, mc, config) {
  log_message("Analyzing RNA-Protein regulatory relationships...")

  rna_mat <- harmonized$transcriptomics$normalized_matrix
  prot_mat <- harmonized$proteomics$normalized_matrix

  # Get common samples
  common_samples <- intersect(colnames(rna_mat), colnames(prot_mat))
  if (length(common_samples) < 5) {
    log_message("  Insufficient common samples. Skipping.")
    return(NULL)
  }

  rna_mat <- rna_mat[, common_samples, drop = FALSE]
  prot_mat <- prot_mat[, common_samples, drop = FALSE]

  # Map to common gene symbols
  mapping_result <- map_rna_protein_features(rna_mat, prot_mat, gene_mapping)
  if (is.null(mapping_result)) {
    log_message("  Could not map RNA and protein features. Skipping.")
    return(NULL)
  }

  rna_matched <- mapping_result$rna
  prot_matched <- mapping_result$protein
  gene_symbols <- mapping_result$gene_symbols

  log_message("  Analyzing ", length(gene_symbols), " matched genes")

  results <- list()

  # 1. Per-gene RNA-Protein correlation analysis
  results$correlations <- compute_rna_protein_correlations(
    rna_matched, prot_matched, gene_symbols
  )

  # 2. Translation efficiency inference
  results$translation_efficiency <- compute_translation_efficiency(
    rna_matched, prot_matched, gene_symbols
  )

  # 3. Post-transcriptional regulation detection
  results$post_transcriptional <- detect_post_transcriptional_regulation(
    rna_matched, prot_matched, gene_symbols, mc
  )

  # 4. Protein stability inference
  results$protein_stability <- infer_protein_stability(
    rna_matched, prot_matched, gene_symbols
  )

  # Save results
  save_rna_protein_results(results, config)

  # Visualizations
  plot_rna_protein_regulation(results, config)

  return(results)
}

#' Map RNA and protein features to common identifiers
map_rna_protein_features <- function(rna_mat, prot_mat, gene_mapping) {
  if (is.null(gene_mapping)) {
    # Try direct matching
    common <- intersect(rownames(rna_mat), rownames(prot_mat))
    if (length(common) < 50) return(NULL)

    return(list(
      rna = rna_mat[common, , drop = FALSE],
      protein = prot_mat[common, , drop = FALSE],
      gene_symbols = common
    ))
  }

  # Use gene mapping
  rna_map <- gene_mapping[gene_mapping$omics == "transcriptomics", ]
  prot_map <- gene_mapping[gene_mapping$omics == "proteomics", ]

  common_symbols <- intersect(rna_map$gene_symbol, prot_map$gene_symbol)
  common_symbols <- common_symbols[!is.na(common_symbols) & common_symbols != ""]

  if (length(common_symbols) < 50) return(NULL)

  # Get feature IDs
  rna_features <- rna_map$feature_id[match(common_symbols, rna_map$gene_symbol)]
  prot_features <- prot_map$feature_id[match(common_symbols, prot_map$gene_symbol)]

  # Filter to present features
  valid <- rna_features %in% rownames(rna_mat) & prot_features %in% rownames(prot_mat)
  common_symbols <- common_symbols[valid]
  rna_features <- rna_features[valid]
  prot_features <- prot_features[valid]

  if (length(common_symbols) < 50) return(NULL)

  list(
    rna = rna_mat[rna_features, , drop = FALSE],
    protein = prot_mat[prot_features, , drop = FALSE],
    gene_symbols = common_symbols
  )
}

#' Compute RNA-Protein correlations per gene
compute_rna_protein_correlations <- function(rna_mat, prot_mat, gene_symbols) {
  n_genes <- length(gene_symbols)

  correlations <- numeric(n_genes)
  pvalues <- numeric(n_genes)
  names(correlations) <- names(pvalues) <- gene_symbols

  for (i in seq_len(n_genes)) {
    rna_expr <- rna_mat[i, ]
    prot_expr <- prot_mat[i, ]

    if (sd(rna_expr, na.rm = TRUE) == 0 || sd(prot_expr, na.rm = TRUE) == 0) {
      correlations[i] <- NA
      pvalues[i] <- NA
      next
    }

    ct <- tryCatch({
      cor.test(rna_expr, prot_expr, method = "spearman", use = "pairwise.complete.obs")
    }, error = function(e) NULL)

    if (!is.null(ct)) {
      correlations[i] <- ct$estimate
      pvalues[i] <- ct$p.value
    } else {
      correlations[i] <- NA
      pvalues[i] <- NA
    }
  }

  padj <- p.adjust(pvalues, method = "BH")

  cor_df <- data.frame(
    gene_symbol = gene_symbols,
    rna_protein_cor = correlations,
    pvalue = pvalues,
    padj = padj,
    stringsAsFactors = FALSE
  )

  # Summary
  valid_cors <- correlations[!is.na(correlations)]
  summary_stats <- list(
    n_genes = n_genes,
    n_valid = length(valid_cors),
    mean_cor = mean(valid_cors),
    median_cor = median(valid_cors),
    n_positive = sum(valid_cors > 0),
    n_strong_positive = sum(valid_cors > 0.5),
    n_negative = sum(valid_cors < 0),
    pct_positive = 100 * mean(valid_cors > 0)
  )

  log_message("  RNA-Protein correlations: mean=", round(summary_stats$mean_cor, 3),
              ", ", summary_stats$n_strong_positive, " with r>0.5")

  list(
    correlations = correlations,
    cor_df = cor_df,
    summary = summary_stats
  )
}

#' Compute translation efficiency scores
compute_translation_efficiency <- function(rna_mat, prot_mat, gene_symbols) {
  # Translation efficiency = log2(protein) - log2(RNA)
  # Higher values indicate more efficient translation

  n_genes <- length(gene_symbols)
  n_samples <- ncol(rna_mat)

  # Compute mean expression
  rna_mean <- rowMeans(rna_mat, na.rm = TRUE)
  prot_mean <- rowMeans(prot_mat, na.rm = TRUE)

  # Add small value to avoid log(0)
  rna_mean <- rna_mean + 1
  prot_mean <- prot_mean + 1

  # Translation efficiency (mean across samples)
  te_scores <- log2(prot_mean) - log2(rna_mean)
  names(te_scores) <- gene_symbols

  # Per-sample TE
  te_matrix <- log2(prot_mat + 1) - log2(rna_mat + 1)

  # TE variability across samples
  te_cv <- apply(te_matrix, 1, function(x) {
    sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))
  })

  te_df <- data.frame(
    gene_symbol = gene_symbols,
    translation_efficiency = te_scores,
    te_cv = te_cv,
    rna_mean = log2(rna_mean),
    protein_mean = log2(prot_mean),
    stringsAsFactors = FALSE
  )

  # Classify genes by TE
  te_df$te_class <- cut(te_df$translation_efficiency,
                        breaks = c(-Inf, -1, 1, Inf),
                        labels = c("low_TE", "moderate_TE", "high_TE"))

  # Summary
  summary_stats <- list(
    mean_te = mean(te_scores, na.rm = TRUE),
    sd_te = sd(te_scores, na.rm = TRUE),
    n_high_te = sum(te_df$te_class == "high_TE", na.rm = TRUE),
    n_low_te = sum(te_df$te_class == "low_TE", na.rm = TRUE)
  )

  log_message("  Translation efficiency: mean=", round(summary_stats$mean_te, 3),
              ", ", summary_stats$n_high_te, " high TE, ",
              summary_stats$n_low_te, " low TE")

  list(
    te_scores = te_scores,
    te_matrix = te_matrix,
    te_df = te_df,
    summary = summary_stats
  )
}

#' Detect post-transcriptional regulation
detect_post_transcriptional_regulation <- function(rna_mat, prot_mat, gene_symbols, mc) {
  # Genes with low RNA-protein correlation but high protein variability
  # may be subject to post-transcriptional regulation

  n_genes <- length(gene_symbols)

  # RNA variability
  rna_var <- apply(rna_mat, 1, var, na.rm = TRUE)
  rna_cv <- apply(rna_mat, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))

  # Protein variability
  prot_var <- apply(prot_mat, 1, var, na.rm = TRUE)
  prot_cv <- apply(prot_mat, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))

  # RNA-protein correlation
  cors <- sapply(seq_len(n_genes), function(i) {
    cor(rna_mat[i, ], prot_mat[i, ], use = "pairwise.complete.obs")
  })

  # Residual protein variance not explained by RNA
  residual_var <- sapply(seq_len(n_genes), function(i) {
    fit <- tryCatch({
      lm(prot_mat[i, ] ~ rna_mat[i, ])
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      var(residuals(fit))
    } else {
      NA
    }
  })

  # Fraction of protein variance explained by RNA
  r_squared <- cors^2
  unexplained_var <- 1 - r_squared

  # Post-transcriptional regulation score
  # High score = low correlation + high protein variability
  ptr_score <- unexplained_var * prot_cv / (prot_cv + 0.1)

  ptr_df <- data.frame(
    gene_symbol = gene_symbols,
    rna_cv = rna_cv,
    protein_cv = prot_cv,
    rna_protein_cor = cors,
    r_squared = r_squared,
    unexplained_variance = unexplained_var,
    residual_variance = residual_var,
    ptr_score = ptr_score,
    stringsAsFactors = FALSE
  )

  # Identify likely PTR-regulated genes
  # High PTR score AND low correlation
  ptr_df$likely_ptr <- ptr_df$rna_protein_cor < 0.3 &
    ptr_df$protein_cv > median(ptr_df$protein_cv, na.rm = TRUE)

  n_ptr <- sum(ptr_df$likely_ptr, na.rm = TRUE)
  log_message("  Post-transcriptional regulation: ", n_ptr,
              " genes with likely PTR")

  list(
    ptr_df = ptr_df,
    n_ptr = n_ptr,
    ptr_genes = gene_symbols[ptr_df$likely_ptr & !is.na(ptr_df$likely_ptr)]
  )
}

#' Infer protein stability
infer_protein_stability <- function(rna_mat, prot_mat, gene_symbols) {
  # Protein stability inference based on deviation from expected
  # protein levels given RNA levels

  n_genes <- length(gene_symbols)

  # Fit global RNA-protein relationship
  rna_mean <- rowMeans(rna_mat, na.rm = TRUE)
  prot_mean <- rowMeans(prot_mat, na.rm = TRUE)

  global_fit <- lm(log2(prot_mean + 1) ~ log2(rna_mean + 1))

  # Expected protein from global fit
  expected_prot <- predict(global_fit, newdata = data.frame(`log2(rna_mean + 1)` = log2(rna_mean + 1)))

  # Deviation from expected (stability score)
  # Positive = more stable (higher protein than expected from RNA)
  stability_score <- log2(prot_mean + 1) - expected_prot

  stability_df <- data.frame(
    gene_symbol = gene_symbols,
    rna_mean = rna_mean,
    protein_mean = prot_mean,
    expected_protein = 2^expected_prot - 1,
    stability_score = stability_score,
    stringsAsFactors = FALSE
  )

  # Classify
  stability_df$stability_class <- cut(
    stability_df$stability_score,
    breaks = c(-Inf, -0.5, 0.5, Inf),
    labels = c("unstable", "normal", "stable")
  )

  # Summary
  summary_stats <- list(
    n_stable = sum(stability_df$stability_class == "stable", na.rm = TRUE),
    n_unstable = sum(stability_df$stability_class == "unstable", na.rm = TRUE),
    mean_stability = mean(stability_score, na.rm = TRUE),
    global_r2 = summary(global_fit)$r.squared
  )

  log_message("  Protein stability: ", summary_stats$n_stable, " stable, ",
              summary_stats$n_unstable, " unstable (global R2=",
              round(summary_stats$global_r2, 3), ")")

  list(
    stability_df = stability_df,
    global_fit = global_fit,
    summary = summary_stats
  )
}

#' Save RNA-protein regulation results
save_rna_protein_results <- function(results, config) {
  if (!is.null(results$correlations)) {
    save_table(results$correlations$cor_df, "rna_protein_correlations.csv", config)
  }

  if (!is.null(results$translation_efficiency)) {
    save_table(results$translation_efficiency$te_df, "translation_efficiency.csv", config)
  }

  if (!is.null(results$post_transcriptional)) {
    save_table(results$post_transcriptional$ptr_df, "post_transcriptional_regulation.csv", config)
  }

  if (!is.null(results$protein_stability)) {
    save_table(results$protein_stability$stability_df, "protein_stability.csv", config)
  }
}

#' Plot RNA-protein regulation
plot_rna_protein_regulation <- function(results, config) {
  # Plot 1: RNA-Protein correlation distribution
  if (!is.null(results$correlations)) {
    cor_df <- results$correlations$cor_df

    p1 <- ggplot2::ggplot(cor_df, ggplot2::aes(x = rna_protein_cor)) +
      ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = median(cor_df$rna_protein_cor, na.rm = TRUE),
                          linetype = "dashed", color = "darkgreen") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "RNA-Protein Correlation Distribution",
        subtitle = paste0("Median = ", round(median(cor_df$rna_protein_cor, na.rm = TRUE), 3)),
        x = "Spearman Correlation",
        y = "Count"
      )

    save_plot(p1, "rna_protein_correlation_dist.png", config, width = 8, height = 5)
  }

  # Plot 2: Translation efficiency
  if (!is.null(results$translation_efficiency)) {
    te_df <- results$translation_efficiency$te_df

    p2 <- ggplot2::ggplot(te_df, ggplot2::aes(x = rna_mean, y = protein_mean,
                                               color = te_class)) +
      ggplot2::geom_point(alpha = 0.5, size = 1.5) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "RNA vs Protein Expression",
        subtitle = "Colored by translation efficiency",
        x = "log2(RNA)",
        y = "log2(Protein)",
        color = "Translation\nEfficiency"
      ) +
      ggplot2::scale_color_manual(values = c("low_TE" = "#d73027",
                                              "moderate_TE" = "#fee090",
                                              "high_TE" = "#4575b4"))

    save_plot(p2, "translation_efficiency_scatter.png", config, width = 8, height = 6)

    # TE distribution
    p2b <- ggplot2::ggplot(te_df, ggplot2::aes(x = translation_efficiency)) +
      ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Translation Efficiency Distribution",
        x = "Translation Efficiency (log2 Protein - log2 RNA)",
        y = "Count"
      )

    save_plot(p2b, "translation_efficiency_dist.png", config, width = 8, height = 5)
  }

  # Plot 3: Protein stability
  if (!is.null(results$protein_stability)) {
    stab_df <- results$protein_stability$stability_df

    p3 <- ggplot2::ggplot(stab_df, ggplot2::aes(x = stability_score, fill = stability_class)) +
      ggplot2::geom_histogram(bins = 50, color = "white") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Protein Stability Score Distribution",
        subtitle = "Deviation from expected protein level given RNA",
        x = "Stability Score",
        y = "Count",
        fill = "Stability"
      ) +
      ggplot2::scale_fill_manual(values = c("unstable" = "#d73027",
                                             "normal" = "#fee090",
                                             "stable" = "#4575b4"))

    save_plot(p3, "protein_stability_dist.png", config, width = 8, height = 5)
  }
}

# =============================================================================
# 2. Transcription Factor Activity Inference
# =============================================================================

#' Infer transcription factor activity
#' @param transcriptomics Transcriptomics data
#' @param mc Mechanistic config
#' @param config Full config
infer_tf_activity <- function(transcriptomics, mc, config) {
  log_message("Inferring transcription factor activity...")

  # Check if decoupleR or dorothea is available
  has_decoupler <- requireNamespace("decoupleR", quietly = TRUE)
  has_dorothea <- requireNamespace("dorothea", quietly = TRUE)

  if (!has_decoupler && !has_dorothea) {
    log_message("  Neither decoupleR nor dorothea available. Skipping TF activity.")
    return(NULL)
  }

  expr_mat <- transcriptomics$normalized_matrix

  # Get TF-target regulons
  regulons <- get_tf_regulons(mc$tf_database)
  if (is.null(regulons)) {
    log_message("  Could not load TF regulons. Skipping.")
    return(NULL)
  }

  log_message("  Using ", length(unique(regulons$tf)), " TFs with ",
              nrow(regulons), " TF-target pairs")

  # Infer TF activity
  if (has_decoupler) {
    tf_activity <- infer_tf_activity_decoupler(expr_mat, regulons)
  } else {
    tf_activity <- infer_tf_activity_dorothea(expr_mat, regulons)
  }

  if (is.null(tf_activity)) return(NULL)

  # Test TF activity associations with conditions
  if (!is.null(transcriptomics$metadata)) {
    tf_activity$associations <- test_tf_associations(
      tf_activity$activity_matrix,
      transcriptomics$metadata,
      config
    )
  }

  # Save results
  save_tf_activity_results(tf_activity, config)

  # Visualizations
  plot_tf_activity(tf_activity, config)

  return(tf_activity)
}

#' Get TF-target regulons
get_tf_regulons <- function(database) {
  if (database == "dorothea") {
    if (!requireNamespace("dorothea", quietly = TRUE)) return(NULL)

    # Get high-confidence regulons (A and B)
    regulons <- tryCatch({
      data("dorothea_hs", package = "dorothea", envir = environment())
      dorothea_hs <- get("dorothea_hs", envir = environment())
      dorothea_hs[dorothea_hs$confidence %in% c("A", "B"), ]
    }, error = function(e) NULL)

    if (!is.null(regulons)) {
      # Standardize column names
      regulons <- data.frame(
        tf = regulons$tf,
        target = regulons$target,
        mor = regulons$mor,  # mode of regulation
        confidence = regulons$confidence,
        stringsAsFactors = FALSE
      )
    }

    return(regulons)
  }

  # Other databases could be added here
  return(NULL)
}

#' Infer TF activity using decoupleR
infer_tf_activity_decoupler <- function(expr_mat, regulons) {
  log_message("  Using decoupleR for TF activity inference...")

  # Run univariate linear model (ULM) or VIPER
  activity <- tryCatch({
    decoupleR::run_ulm(
      mat = expr_mat,
      net = regulons,
      .source = "tf",
      .target = "target",
      .mor = "mor",
      minsize = 5
    )
  }, error = function(e) {
    log_message("  decoupleR::run_ulm failed: ", e$message)
    NULL
  })

  if (is.null(activity)) return(NULL)

  # Pivot to wide format
  activity_mat <- tryCatch({
    tidyr::pivot_wider(activity,
                       id_cols = source,
                       names_from = condition,
                       values_from = score)
  }, error = function(e) NULL)

  if (is.null(activity_mat)) return(NULL)

  # Convert to matrix
  tf_names <- activity_mat$source
  activity_mat <- as.matrix(activity_mat[, -1])
  rownames(activity_mat) <- tf_names

  list(
    activity_matrix = activity_mat,
    activity_long = activity,
    method = "decoupleR_ulm"
  )
}

#' Infer TF activity using dorothea directly
infer_tf_activity_dorothea <- function(expr_mat, regulons) {
  log_message("  Using dorothea for TF activity inference...")

  # Simple VIPER-like approach
  # For each TF, compute enrichment score of its targets

  tfs <- unique(regulons$tf)
  n_tfs <- length(tfs)
  n_samples <- ncol(expr_mat)

  activity_mat <- matrix(NA, n_tfs, n_samples)
  rownames(activity_mat) <- tfs
  colnames(activity_mat) <- colnames(expr_mat)

  for (i in seq_along(tfs)) {
    tf <- tfs[i]
    targets <- regulons[regulons$tf == tf, ]

    # Get target expression
    target_genes <- intersect(targets$target, rownames(expr_mat))
    if (length(target_genes) < 5) next

    target_expr <- expr_mat[target_genes, , drop = FALSE]
    target_mor <- targets$mor[match(target_genes, targets$target)]

    # Weight by mode of regulation
    weighted_expr <- sweep(target_expr, 1, target_mor, `*`)

    # TF activity = mean of weighted target expression
    activity_mat[i, ] <- colMeans(weighted_expr, na.rm = TRUE)
  }

  # Scale activity scores
  activity_mat <- t(scale(t(activity_mat)))

  # Remove TFs with no valid activity
  valid_tfs <- !apply(is.na(activity_mat), 1, all)
  activity_mat <- activity_mat[valid_tfs, , drop = FALSE]

  list(
    activity_matrix = activity_mat,
    method = "dorothea_weighted_mean"
  )
}

#' Test TF activity associations with conditions
test_tf_associations <- function(activity_mat, metadata, config) {
  condition_col <- config$design$condition_column %||% "condition"

  if (!condition_col %in% colnames(metadata)) return(NULL)

  # Get conditions for samples
  samples <- intersect(colnames(activity_mat), rownames(metadata))
  if (length(samples) < 5) return(NULL)

  activity_mat <- activity_mat[, samples, drop = FALSE]
  conditions <- metadata[samples, condition_col]

  # Test each TF
  tfs <- rownames(activity_mat)

  results <- lapply(tfs, function(tf) {
    activity <- activity_mat[tf, ]

    # ANOVA or t-test
    if (length(unique(conditions)) == 2) {
      test <- tryCatch({
        t.test(activity ~ conditions)
      }, error = function(e) NULL)

      if (!is.null(test)) {
        return(data.frame(
          tf = tf,
          test = "t.test",
          statistic = test$statistic,
          pvalue = test$p.value,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      test <- tryCatch({
        summary(aov(activity ~ conditions))[[1]]
      }, error = function(e) NULL)

      if (!is.null(test)) {
        return(data.frame(
          tf = tf,
          test = "anova",
          statistic = test[1, "F value"],
          pvalue = test[1, "Pr(>F)"],
          stringsAsFactors = FALSE
        ))
      }
    }

    NULL
  })

  results_df <- do.call(rbind, results)
  if (is.null(results_df)) return(NULL)

  results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
  results_df <- results_df[order(results_df$pvalue), ]

  n_sig <- sum(results_df$padj < 0.05, na.rm = TRUE)
  log_message("  TF associations: ", n_sig, " significantly associated with condition")

  results_df
}

#' Save TF activity results
save_tf_activity_results <- function(tf_activity, config) {
  # Save activity matrix
  activity_df <- as.data.frame(tf_activity$activity_matrix)
  activity_df$tf <- rownames(activity_df)
  save_table(activity_df, "tf_activity_scores.csv", config)

  # Save associations if available
  if (!is.null(tf_activity$associations)) {
    save_table(tf_activity$associations, "tf_condition_associations.csv", config)
  }
}

#' Plot TF activity
plot_tf_activity <- function(tf_activity, config) {
  activity_mat <- tf_activity$activity_matrix

  # Heatmap of top variable TFs
  tf_var <- apply(activity_mat, 1, var, na.rm = TRUE)
  top_tfs <- names(sort(tf_var, decreasing = TRUE))[1:min(50, length(tf_var))]

  mat_plot <- activity_mat[top_tfs, , drop = FALSE]

  # Use pheatmap if available
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    png(file.path(config$output$output_dir, "plots", "tf_activity_heatmap.png"),
        width = 12, height = 10, units = "in", res = 300)

    tryCatch({
      pheatmap::pheatmap(
        mat_plot,
        scale = "row",
        clustering_method = "ward.D2",
        show_colnames = ncol(mat_plot) < 50,
        main = "Transcription Factor Activity",
        color = colorRampPalette(c("blue", "white", "red"))(100)
      )
    }, error = function(e) {
      log_message("  TF heatmap failed: ", e$message)
    })

    dev.off()
  }

  # Bar plot of most significant TFs
  if (!is.null(tf_activity$associations)) {
    assoc <- tf_activity$associations
    assoc_sig <- assoc[assoc$padj < 0.1, ]

    if (nrow(assoc_sig) > 0) {
      assoc_sig <- head(assoc_sig[order(assoc_sig$pvalue), ], 20)

      p <- ggplot2::ggplot(assoc_sig, ggplot2::aes(x = reorder(tf, -log10(pvalue)),
                                                    y = -log10(pvalue))) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        ggplot2::coord_flip() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "TF Activity Association with Condition",
          x = "Transcription Factor",
          y = "-log10(p-value)"
        )

      save_plot(p, "tf_association_barplot.png", config, width = 8, height = 6)
    }
  }
}

# =============================================================================
# 3. Gene Regulatory Network Inference
# =============================================================================

#' Infer gene regulatory network
#' @param harmonized List of harmonized omics data
#' @param gene_mapping Gene mapping table
#' @param mc Mechanistic config
#' @param config Full config
infer_regulatory_network <- function(harmonized, gene_mapping, mc, config) {
  log_message("Inferring gene regulatory network...")

  # Use transcriptomics for network inference
  if (!"transcriptomics" %in% names(harmonized)) {
    log_message("  Transcriptomics data required for network inference. Skipping.")
    return(NULL)
  }

  expr_mat <- harmonized$transcriptomics$normalized_matrix

  # Select top variable genes
  gene_vars <- apply(expr_mat, 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(500, nrow(expr_mat))]
  expr_sub <- expr_mat[top_genes, ]

  log_message("  Using ", length(top_genes), " genes for network inference")

  # Run network inference
  if (mc$network_method == "genie3") {
    network <- infer_network_genie3(expr_sub, mc)
  } else if (mc$network_method == "aracne") {
    network <- infer_network_aracne(expr_sub, mc)
  } else {
    network <- infer_network_correlation(expr_sub, mc)
  }

  if (is.null(network)) return(NULL)

  # Identify hub regulators
  network$hubs <- identify_hub_regulators(network$edges, mc)

  # Save results
  save_table(network$edges, "regulatory_network_edges.csv", config)
  save_table(network$hubs, "hub_regulators.csv", config)

  # Visualizations
  plot_regulatory_network(network, config)

  return(network)
}

#' Infer network using GENIE3
infer_network_genie3 <- function(expr_mat, mc) {
  if (!requireNamespace("GENIE3", quietly = TRUE)) {
    log_message("  GENIE3 not available. Using correlation-based network.")
    return(infer_network_correlation(expr_mat, mc))
  }

  log_message("  Running GENIE3...")

  # Run GENIE3
  weight_mat <- tryCatch({
    GENIE3::GENIE3(
      exprMatrix = expr_mat,
      nCores = 1,
      verbose = FALSE
    )
  }, error = function(e) {
    log_message("  GENIE3 failed: ", e$message)
    NULL
  })

  if (is.null(weight_mat)) {
    return(infer_network_correlation(expr_mat, mc))
  }

  # Convert to edge list
  link_list <- GENIE3::getLinkList(weight_mat)

  # Filter to top edges
  n_edges <- min(5000, nrow(link_list))
  edges <- link_list[1:n_edges, ]
  colnames(edges) <- c("regulator", "target", "weight")

  log_message("  GENIE3 network: ", nrow(edges), " edges")

  list(
    edges = edges,
    method = "GENIE3"
  )
}

#' Infer network using ARACNe
infer_network_aracne <- function(expr_mat, mc) {
  if (!requireNamespace("minet", quietly = TRUE)) {
    log_message("  minet not available. Using correlation-based network.")
    return(infer_network_correlation(expr_mat, mc))
  }

  log_message("  Running ARACNe (via minet)...")

  # Compute mutual information
  mi_mat <- tryCatch({
    minet::build.mim(t(expr_mat), estimator = "spearman")
  }, error = function(e) {
    log_message("  MI estimation failed: ", e$message)
    NULL
  })

  if (is.null(mi_mat)) {
    return(infer_network_correlation(expr_mat, mc))
  }

  # ARACNe network
  aracne_mat <- tryCatch({
    minet::aracne(mi_mat)
  }, error = function(e) NULL)

  if (is.null(aracne_mat)) {
    return(infer_network_correlation(expr_mat, mc))
  }

  # Convert to edge list
  edges <- mat_to_edge_list(aracne_mat, "MI")

  log_message("  ARACNe network: ", nrow(edges), " edges")

  list(
    edges = edges,
    method = "ARACNe"
  )
}

#' Infer network using correlation
infer_network_correlation <- function(expr_mat, mc) {
  log_message("  Using correlation-based network...")

  # Compute correlation matrix
  cor_mat <- cor(t(expr_mat), use = "pairwise.complete.obs", method = "spearman")

  # Convert to edge list (filter by threshold)
  edges <- mat_to_edge_list(cor_mat, "correlation", threshold = 0.5)

  log_message("  Correlation network: ", nrow(edges), " edges")

  list(
    edges = edges,
    method = "correlation"
  )
}

#' Convert matrix to edge list
mat_to_edge_list <- function(mat, weight_name, threshold = 0) {
  genes <- rownames(mat)
  n <- length(genes)

  edges <- data.frame()

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      w <- mat[i, j]
      if (!is.na(w) && abs(w) > threshold) {
        edges <- rbind(edges, data.frame(
          regulator = genes[i],
          target = genes[j],
          weight = w,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  colnames(edges)[3] <- weight_name
  edges <- edges[order(abs(edges[[3]]), decreasing = TRUE), ]

  edges
}

#' Identify hub regulators
identify_hub_regulators <- function(edges, mc) {
  # Count outgoing edges per regulator
  reg_counts <- table(edges$regulator)
  reg_counts <- sort(reg_counts, decreasing = TRUE)

  hub_df <- data.frame(
    gene = names(reg_counts),
    n_targets = as.numeric(reg_counts),
    stringsAsFactors = FALSE
  )

  # Mark top regulators
  hub_df$is_hub <- hub_df$n_targets >= quantile(hub_df$n_targets, 0.9)

  log_message("  Identified ", sum(hub_df$is_hub), " hub regulators")

  hub_df
}

#' Plot regulatory network
plot_regulatory_network <- function(network, config) {
  edges <- network$edges

  # Use igraph for visualization
  if (!requireNamespace("igraph", quietly = TRUE)) return(NULL)

  # Create graph from top edges
  top_edges <- head(edges, 500)

  g <- igraph::graph_from_data_frame(
    top_edges[, c("regulator", "target")],
    directed = TRUE
  )

  # Identify hubs
  if (!is.null(network$hubs)) {
    hub_genes <- network$hubs$gene[network$hubs$is_hub]
    igraph::V(g)$is_hub <- igraph::V(g)$name %in% hub_genes
  }

  # Save as PNG
  png(file.path(config$output$output_dir, "plots", "regulatory_network.png"),
      width = 12, height = 12, units = "in", res = 300)

  tryCatch({
    # Color hubs differently
    v_colors <- ifelse(igraph::V(g)$is_hub, "#d73027", "#4575b4")

    igraph::plot.igraph(
      g,
      vertex.size = ifelse(igraph::V(g)$is_hub, 8, 4),
      vertex.color = v_colors,
      vertex.label = ifelse(igraph::V(g)$is_hub, igraph::V(g)$name, NA),
      vertex.label.cex = 0.7,
      edge.arrow.size = 0.3,
      edge.color = "gray70",
      layout = igraph::layout_with_fr(g),
      main = paste0("Regulatory Network (", network$method, ")")
    )
  }, error = function(e) {
    log_message("  Network plot failed: ", e$message)
  })

  dev.off()
}

# =============================================================================
# 4. Causal Mediation Analysis
# =============================================================================

#' Run mediation analysis: RNA -> Protein -> Phenotype
#' @param harmonized List of harmonized omics data
#' @param metadata Sample metadata
#' @param gene_mapping Gene mapping table
#' @param mc Mechanistic config
#' @param config Full config
run_mediation_analysis <- function(harmonized, metadata, gene_mapping, mc, config) {
  log_message("Running mediation analysis (RNA -> Protein -> Phenotype)...")

  # Check requirements
  if (!"transcriptomics" %in% names(harmonized) ||
      !"proteomics" %in% names(harmonized)) {
    log_message("  Both transcriptomics and proteomics required. Skipping.")
    return(NULL)
  }

  if (!requireNamespace("mediation", quietly = TRUE)) {
    log_message("  mediation package not available. Skipping.")
    return(NULL)
  }

  # Get outcome variable
  outcome_col <- config$design$outcome_column %||% config$design$condition_column
  if (is.null(outcome_col) || !outcome_col %in% colnames(metadata)) {
    log_message("  No outcome variable specified. Skipping mediation.")
    return(NULL)
  }

  # Map features
  rna_mat <- harmonized$transcriptomics$normalized_matrix
  prot_mat <- harmonized$proteomics$normalized_matrix

  mapping_result <- map_rna_protein_features(rna_mat, prot_mat, gene_mapping)
  if (is.null(mapping_result)) return(NULL)

  # Get common samples with outcome
  common_samples <- Reduce(intersect, list(
    colnames(rna_mat),
    colnames(prot_mat),
    rownames(metadata)
  ))

  if (length(common_samples) < 20) {
    log_message("  Insufficient samples for mediation. Skipping.")
    return(NULL)
  }

  rna_matched <- mapping_result$rna[, common_samples, drop = FALSE]
  prot_matched <- mapping_result$protein[, common_samples, drop = FALSE]
  gene_symbols <- mapping_result$gene_symbols
  outcome <- metadata[common_samples, outcome_col]

  # Convert outcome to numeric if factor
  if (is.factor(outcome) || is.character(outcome)) {
    outcome <- as.numeric(as.factor(outcome)) - 1
  }

  log_message("  Testing ", length(gene_symbols), " genes for mediation")

  # Run mediation for each gene (limit to significant genes)
  # First, test for association with outcome
  rna_assoc <- apply(rna_matched, 1, function(x) {
    tryCatch({
      cor.test(x, outcome, method = "spearman")$p.value
    }, error = function(e) 1)
  })

  # Focus on genes associated with outcome
  sig_genes <- gene_symbols[rna_assoc < 0.1]
  sig_genes <- head(sig_genes, 100)  # Limit for computation time

  if (length(sig_genes) < 5) {
    log_message("  Too few genes associated with outcome. Skipping.")
    return(NULL)
  }

  log_message("  Running mediation for ", length(sig_genes), " candidate genes")

  # Run mediation
  mediation_results <- lapply(seq_along(sig_genes), function(i) {
    gene <- sig_genes[i]
    gene_idx <- which(gene_symbols == gene)

    rna_expr <- rna_matched[gene_idx, ]
    prot_expr <- prot_matched[gene_idx, ]

    run_single_mediation(rna_expr, prot_expr, outcome, gene)
  })

  # Combine results
  mediation_df <- do.call(rbind, mediation_results)
  mediation_df <- mediation_df[!is.na(mediation_df$prop_mediated), ]

  if (nrow(mediation_df) == 0) {
    log_message("  No successful mediation analyses.")
    return(NULL)
  }

  # Adjust p-values
  mediation_df$acme_padj <- p.adjust(mediation_df$acme_pval, method = "BH")

  # Sort by proportion mediated
  mediation_df <- mediation_df[order(mediation_df$prop_mediated, decreasing = TRUE), ]

  # Summary
  n_mediated <- sum(mediation_df$prop_mediated > 0.1 &
                      mediation_df$acme_padj < 0.1, na.rm = TRUE)
  log_message("  Found ", n_mediated, " genes with significant mediation")

  # Save results
  save_table(mediation_df, "mediation_analysis_results.csv", config)

  # Visualizations
  plot_mediation_results(mediation_df, config)

  list(
    mediation_df = mediation_df,
    n_mediated = n_mediated
  )
}

#' Run single gene mediation analysis
run_single_mediation <- function(rna_expr, prot_expr, outcome, gene_name) {
  tryCatch({
    # Create data frame
    df <- data.frame(
      rna = as.numeric(rna_expr),
      protein = as.numeric(prot_expr),
      outcome = as.numeric(outcome)
    )
    df <- df[complete.cases(df), ]

    if (nrow(df) < 20) {
      return(data.frame(gene = gene_name, prop_mediated = NA, acme_pval = NA,
                        ade_pval = NA, stringsAsFactors = FALSE))
    }

    # Mediator model: Protein ~ RNA
    med_model <- lm(protein ~ rna, data = df)

    # Outcome model: Outcome ~ RNA + Protein
    out_model <- lm(outcome ~ rna + protein, data = df)

    # Run mediation
    med_result <- mediation::mediate(
      med_model, out_model,
      treat = "rna",
      mediator = "protein",
      boot = FALSE,  # Use analytical estimates for speed
      sims = 100
    )

    data.frame(
      gene = gene_name,
      acme = med_result$d0,  # Average Causal Mediation Effect
      acme_pval = med_result$d0.p,
      ade = med_result$z0,   # Average Direct Effect
      ade_pval = med_result$z0.p,
      total_effect = med_result$tau.coef,
      prop_mediated = med_result$n0,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(gene = gene_name, prop_mediated = NA, acme_pval = NA,
               stringsAsFactors = FALSE)
  })
}

#' Plot mediation results
plot_mediation_results <- function(mediation_df, config) {
  # Filter to significant
  med_sig <- mediation_df[mediation_df$acme_padj < 0.1 & !is.na(mediation_df$acme_padj), ]

  if (nrow(med_sig) == 0) {
    med_sig <- head(mediation_df[!is.na(mediation_df$prop_mediated), ], 20)
  }

  if (nrow(med_sig) == 0) return(NULL)

  # Bar plot of proportion mediated
  med_sig <- head(med_sig[order(med_sig$prop_mediated, decreasing = TRUE), ], 20)

  p1 <- ggplot2::ggplot(med_sig, ggplot2::aes(x = reorder(gene, prop_mediated),
                                               y = prop_mediated)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Mediation Analysis: RNA -> Protein -> Phenotype",
      subtitle = "Proportion of RNA effect mediated by protein",
      x = "Gene",
      y = "Proportion Mediated"
    ) +
    ggplot2::ylim(0, 1)

  save_plot(p1, "mediation_proportion_barplot.png", config, width = 8, height = 6)

  # Effect decomposition
  if (all(c("acme", "ade") %in% colnames(med_sig))) {
    med_long <- tidyr::pivot_longer(med_sig,
                                     cols = c(acme, ade),
                                     names_to = "effect_type",
                                     values_to = "effect_size")

    p2 <- ggplot2::ggplot(med_long, ggplot2::aes(x = reorder(gene, effect_size),
                                                   y = effect_size,
                                                   fill = effect_type)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Mediation Effect Decomposition",
        subtitle = "ACME = mediated effect, ADE = direct effect",
        x = "Gene",
        y = "Effect Size",
        fill = "Effect Type"
      ) +
      ggplot2::scale_fill_manual(values = c("acme" = "#4575b4", "ade" = "#d73027"),
                                  labels = c("acme" = "Mediated (ACME)",
                                            "ade" = "Direct (ADE)"))

    save_plot(p2, "mediation_effect_decomposition.png", config, width = 8, height = 6)
  }
}

# =============================================================================
# Summary Function
# =============================================================================

#' Summarize mechanistic analysis results
summarize_mechanistic_results <- function(results, config) {
  log_message("Generating mechanistic analysis summary...")

  summary_list <- list()

  # RNA-Protein regulation
  if (!is.null(results$rna_protein)) {
    rp <- results$rna_protein
    summary_list$rna_protein <- list(
      mean_correlation = rp$correlations$summary$mean_cor,
      n_high_te = rp$translation_efficiency$summary$n_high_te,
      n_low_te = rp$translation_efficiency$summary$n_low_te,
      n_ptr_genes = rp$post_transcriptional$n_ptr
    )
  }

  # TF activity
  if (!is.null(results$tf_activity)) {
    tf <- results$tf_activity
    n_sig_tfs <- ifelse(!is.null(tf$associations),
                        sum(tf$associations$padj < 0.05, na.rm = TRUE), NA)
    summary_list$tf_activity <- list(
      n_tfs = nrow(tf$activity_matrix),
      n_significant = n_sig_tfs
    )
  }

  # Regulatory network
  if (!is.null(results$regulatory_network)) {
    rn <- results$regulatory_network
    summary_list$regulatory_network <- list(
      method = rn$method,
      n_edges = nrow(rn$edges),
      n_hubs = sum(rn$hubs$is_hub, na.rm = TRUE)
    )
  }

  # Mediation
  if (!is.null(results$mediation)) {
    med <- results$mediation
    summary_list$mediation <- list(
      n_genes_tested = nrow(med$mediation_df),
      n_mediated = med$n_mediated
    )
  }

  # Save summary
  flat_summary <- do.call(rbind, lapply(names(summary_list), function(cat) {
    s <- summary_list[[cat]]
    data.frame(
      category = cat,
      metric = names(s),
      value = as.character(unlist(s)),
      stringsAsFactors = FALSE
    )
  }))

  save_table(flat_summary, "mechanistic_analysis_summary.csv", config)

  log_message("=== Mechanistic Analysis Summary ===")
  for (cat in names(summary_list)) {
    log_message("  ", cat, ": ", paste(names(summary_list[[cat]]), "=",
                                        unlist(summary_list[[cat]]),
                                        collapse = ", "))
  }

  summary_list
}
