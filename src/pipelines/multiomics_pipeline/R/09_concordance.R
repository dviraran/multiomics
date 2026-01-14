# =============================================================================
# Cross-Omics Concordance Analysis
# =============================================================================

#' Run cross-omics concordance analysis
run_concordance_analysis <- function(mae_data, integration_results, config) {
  log_message("=== Running Cross-Omics Concordance Analysis ===")

  harmonized <- mae_data$harmonized_omics
  gene_mapping <- mae_data$gene_mapping

  results <- list()

  # 1. RNA-Protein concordance (if both present)
  if ("transcriptomics" %in% names(harmonized) &&
      "proteomics" %in% names(harmonized)) {
    results$rna_protein <- analyze_rna_protein_concordance(
      harmonized$transcriptomics,
      harmonized$proteomics,
      gene_mapping,
      config
    )
  }

  # 2. Differential expression/abundance concordance
  results$de_concordance <- analyze_de_concordance(harmonized, config)

  # 3. Integration-based concordance (from MOFA/DIABLO)
  if (!is.null(integration_results)) {
    results$integration_concordance <- analyze_integration_concordance(
      integration_results, harmonized, config
    )
  }

  # Summary report
  summarize_concordance(results, config)

  results
}

#' Analyze RNA-Protein concordance
analyze_rna_protein_concordance <- function(rna_data, prot_data, gene_mapping, config) {
  log_message("Analyzing RNA-Protein concordance...")

  rna_mat <- rna_data$normalized_matrix
  prot_mat <- prot_data$normalized_matrix

  # Get gene symbols for matching
  if (is.null(gene_mapping)) {
    log_message("No gene mapping available for RNA-protein comparison")
    return(NULL)
  }

  rna_genes <- gene_mapping[gene_mapping$omics == "transcriptomics", ]
  prot_genes <- gene_mapping[gene_mapping$omics == "proteomics", ]

  # Find common genes
  common_genes <- intersect(rna_genes$gene_symbol, prot_genes$gene_symbol)
  common_genes <- common_genes[!is.na(common_genes)]

  if (length(common_genes) < 10) {
    log_message("Only ", length(common_genes), " genes in common. Skipping concordance.")
    return(NULL)
  }

  log_message("Found ", length(common_genes), " genes in common between RNA and protein")

  # Map to feature IDs
  rna_features <- rna_genes$feature_id[match(common_genes, rna_genes$gene_symbol)]
  prot_features <- prot_genes$feature_id[match(common_genes, prot_genes$gene_symbol)]

  # Get common samples
  common_samples <- intersect(colnames(rna_mat), colnames(prot_mat))

  # Calculate correlations per gene
  correlations <- numeric(length(common_genes))
  names(correlations) <- common_genes

  for (i in seq_along(common_genes)) {
    rna_expr <- rna_mat[rna_features[i], common_samples]
    prot_expr <- prot_mat[prot_features[i], common_samples]

    if (sd(rna_expr, na.rm = TRUE) > 0 && sd(prot_expr, na.rm = TRUE) > 0) {
      correlations[i] <- cor(rna_expr, prot_expr, use = "pairwise.complete.obs")
    } else {
      correlations[i] <- NA
    }
  }

  correlations <- correlations[!is.na(correlations)]

  # Summary statistics
  summary_stats <- list(
    n_genes = length(correlations),
    mean_cor = mean(correlations),
    median_cor = median(correlations),
    sd_cor = sd(correlations),
    n_positive = sum(correlations > 0),
    n_significant = sum(abs(correlations) > 0.5),
    pct_positive = 100 * sum(correlations > 0) / length(correlations)
  )

  log_message("RNA-Protein correlation: mean=", round(summary_stats$mean_cor, 3),
             ", median=", round(summary_stats$median_cor, 3),
             ", ", summary_stats$pct_positive, "% positive")

  # Save correlation table
  cor_df <- data.frame(
    gene_symbol = names(correlations),
    rna_protein_correlation = correlations,
    stringsAsFactors = FALSE
  )
  cor_df <- cor_df[order(abs(cor_df$rna_protein_correlation), decreasing = TRUE), ]
  save_table(cor_df, "rna_protein_correlations.csv", config)

  # Create visualization
  p <- plot_rna_protein_concordance(correlations)
  save_plot(p, "rna_protein_concordance", config, width = 10, height = 6)

  list(
    correlations = correlations,
    summary = summary_stats,
    common_genes = common_genes
  )
}

#' Plot RNA-Protein concordance
plot_rna_protein_concordance <- function(correlations) {
  df <- data.frame(correlation = correlations)

  # Histogram
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = correlation)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = median(correlations), color = "darkgreen",
                        linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "RNA-Protein Correlation Distribution",
      subtitle = paste0("n=", length(correlations), " genes, median=",
                       round(median(correlations), 3)),
      x = "Pearson Correlation",
      y = "Count"
    )

  # Density plot
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = correlation)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.5) +
    ggplot2::geom_rug(alpha = 0.1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Correlation Density",
      x = "Pearson Correlation",
      y = "Density"
    )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p1 + p2
  } else {
    p1
  }
}

#' Analyze differential expression/abundance concordance
analyze_de_concordance <- function(harmonized_omics, config) {
  log_message("Analyzing DE/DA concordance across omics...")

  de_tables <- list()

  # Collect DE/DA tables
  if ("transcriptomics" %in% names(harmonized_omics)) {
    de <- harmonized_omics$transcriptomics$de_table
    if (!is.null(de)) {
      de$omics <- "transcriptomics"
      de$id_col <- de$gene_id
      if ("gene_symbol" %in% colnames(de)) de$id_col <- de$gene_symbol
      de_tables$transcriptomics <- de
    }
  }

  if ("proteomics" %in% names(harmonized_omics)) {
    da <- harmonized_omics$proteomics$da_table
    if (!is.null(da)) {
      da$omics <- "proteomics"
      da$id_col <- da$feature_id
      if ("gene_symbol" %in% colnames(da)) da$id_col <- da$gene_symbol
      de_tables$proteomics <- da
    }
  }

  if ("metabolomics" %in% names(harmonized_omics)) {
    da <- harmonized_omics$metabolomics$da_table
    if (!is.null(da)) {
      da$omics <- "metabolomics"
      da$id_col <- da$feature_id
      if ("compound_name" %in% colnames(da)) da$id_col <- da$compound_name
      de_tables$metabolomics <- da
    }
  }

  if (length(de_tables) < 2) {
    log_message("Need at least 2 omics with DE/DA results for concordance")
    return(NULL)
  }

  # Compare fold changes for matching features (RNA vs Protein)
  concordance_results <- list()

  if ("transcriptomics" %in% names(de_tables) && "proteomics" %in% names(de_tables)) {
    rna_de <- de_tables$transcriptomics
    prot_da <- de_tables$proteomics

    # Match by gene symbol
    common <- intersect(rna_de$id_col, prot_da$id_col)
    common <- common[!is.na(common)]

    if (length(common) >= 10) {
      rna_sub <- rna_de[match(common, rna_de$id_col), ]
      prot_sub <- prot_da[match(common, prot_da$id_col), ]

      fc_cor <- cor(rna_sub$log2FC, prot_sub$log2FC, use = "pairwise.complete.obs")

      # Directional concordance
      same_direction <- sign(rna_sub$log2FC) == sign(prot_sub$log2FC)
      pct_concordant <- 100 * mean(same_direction, na.rm = TRUE)

      concordance_results$rna_protein <- list(
        n_common = length(common),
        fc_correlation = fc_cor,
        pct_directional_concordance = pct_concordant
      )

      log_message("RNA-Protein DE concordance: ",
                 round(fc_cor, 3), " (FC correlation), ",
                 round(pct_concordant, 1), "% directional agreement")

      # Save concordance table
      conc_df <- data.frame(
        gene_symbol = common,
        rna_log2FC = rna_sub$log2FC,
        protein_log2FC = prot_sub$log2FC,
        concordant = same_direction,
        stringsAsFactors = FALSE
      )
      save_table(conc_df, "rna_protein_de_concordance.csv", config)

      # Plot
      p <- plot_de_concordance(rna_sub$log2FC, prot_sub$log2FC, "RNA", "Protein")
      save_plot(p, "rna_protein_de_scatter", config, width = 8, height = 8)
    }
  }

  # Summary of significant features per omics
  sig_summary <- lapply(de_tables, function(dt) {
    padj_col <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(dt))[1]
    if (is.na(padj_col)) return(NULL)

    n_sig <- sum(dt[[padj_col]] < 0.05, na.rm = TRUE)
    n_up <- sum(dt[[padj_col]] < 0.05 & dt$log2FC > 0, na.rm = TRUE)
    n_down <- sum(dt[[padj_col]] < 0.05 & dt$log2FC < 0, na.rm = TRUE)

    data.frame(
      omics = dt$omics[1],
      n_significant = n_sig,
      n_up = n_up,
      n_down = n_down,
      stringsAsFactors = FALSE
    )
  })

  sig_summary_df <- do.call(rbind, sig_summary)
  if (!is.null(sig_summary_df)) {
    save_table(sig_summary_df, "de_da_summary.csv", config)
    log_message("DE/DA summary saved")
  }

  list(
    pairwise = concordance_results,
    sig_summary = sig_summary_df
  )
}

#' Plot DE concordance scatter
plot_de_concordance <- function(fc1, fc2, label1, label2) {
  df <- data.frame(fc1 = fc1, fc2 = fc2)
  df <- df[complete.cases(df), ]

  # Calculate correlation
  r <- cor(df$fc1, df$fc2)

  ggplot2::ggplot(df, ggplot2::aes(x = fc1, y = fc2)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0(label1, " vs ", label2, " Fold Change Concordance"),
      subtitle = paste0("r = ", round(r, 3), ", n = ", nrow(df), " features"),
      x = paste0(label1, " log2(FC)"),
      y = paste0(label2, " log2(FC)")
    ) +
    ggplot2::coord_fixed()
}

#' Analyze concordance from integration results
analyze_integration_concordance <- function(integration_results, harmonized, config) {
  log_message("Analyzing integration-based concordance...")

  results <- list()

  # MOFA-based concordance
  if (!is.null(integration_results$mofa)) {
    mofa_res <- integration_results$mofa$results

    # Cross-omics feature correlations via factor loadings
    weights <- mofa_res$weights
    views <- names(weights)

    if (length(views) >= 2) {
      # For each pair of views, correlate their factor loadings (for shared variance)
      view_cors <- list()

      for (i in 1:(length(views) - 1)) {
        for (j in (i + 1):length(views)) {
          v1 <- views[i]
          v2 <- views[j]

          # Can only correlate if we can map to same genes
          # For now, correlate total variance explained per factor
          r2_1 <- mofa_res$r2_per_factor[v1, ]
          r2_2 <- mofa_res$r2_per_factor[v2, ]

          view_cor <- cor(r2_1, r2_2)
          view_cors[[paste(v1, v2, sep = "_")]] <- view_cor

          log_message("MOFA variance concordance (", v1, " vs ", v2, "): ",
                     round(view_cor, 3))
        }
      }

      results$mofa_variance_concordance <- view_cors
    }
  }

  # DIABLO-based concordance
  if (!is.null(integration_results$diablo)) {
    diablo_res <- integration_results$diablo$results

    # Overlap of selected features between blocks
    selected <- diablo_res$selected_vars
    if (length(selected) >= 2) {
      # This is more meaningful if we can map to genes
      results$diablo_n_selected <- sapply(selected, length)
    }
  }

  results
}

#' Summarize concordance results
summarize_concordance <- function(concordance_results, config) {
  log_message("Generating concordance summary...")

  summary_lines <- c(
    "=== Cross-Omics Concordance Summary ==="
  )

  # RNA-Protein expression concordance
  if (!is.null(concordance_results$rna_protein)) {
    rp <- concordance_results$rna_protein
    summary_lines <- c(summary_lines,
      "",
      "RNA-Protein Expression Concordance:",
      paste0("  - Genes compared: ", rp$summary$n_genes),
      paste0("  - Mean correlation: ", round(rp$summary$mean_cor, 3)),
      paste0("  - Median correlation: ", round(rp$summary$median_cor, 3)),
      paste0("  - Positive correlations: ", rp$summary$pct_positive, "%")
    )
  }

  # DE concordance
  if (!is.null(concordance_results$de_concordance$pairwise$rna_protein)) {
    de_rp <- concordance_results$de_concordance$pairwise$rna_protein
    summary_lines <- c(summary_lines,
      "",
      "RNA-Protein Differential Concordance:",
      paste0("  - Features compared: ", de_rp$n_common),
      paste0("  - Fold change correlation: ", round(de_rp$fc_correlation, 3)),
      paste0("  - Directional concordance: ", round(de_rp$pct_directional_concordance, 1), "%")
    )
  }

  # Write summary
  summary_text <- paste(summary_lines, collapse = "\n")
  summary_file <- file.path(config$output$output_dir, "tables", "concordance_summary.txt")
  writeLines(summary_text, summary_file)

  log_message("Concordance summary saved to ", summary_file)
}
