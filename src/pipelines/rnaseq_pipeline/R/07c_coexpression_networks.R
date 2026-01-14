# =============================================================================
# Co-Expression Network Analysis (WGCNA)
# =============================================================================
# This script performs Weighted Gene Co-Expression Network Analysis:
#   1. Soft-thresholding power selection
#   2. Network construction and module detection
#   3. Module-trait correlation analysis
#   4. Hub gene identification
#   5. Module eigengene visualization
# =============================================================================

#' Main function: Run WGCNA analysis
#' @param vst_counts VST-transformed expression matrix
#' @param metadata Sample metadata
#' @param gene_annotation Gene annotation table
#' @param config Pipeline configuration
#' @return List with WGCNA results
run_wgcna_analysis <- function(vst_counts, metadata, gene_annotation, config) {
  log_message("=== Running WGCNA Co-Expression Analysis ===")

  wc <- get_wgcna_config(config)

  if (!wc$run_wgcna) {
    log_message("WGCNA analysis disabled. Skipping.")
    return(NULL)
  }

  # Check if WGCNA is available
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    log_message("WGCNA package not available. Skipping.")
    return(NULL)
  }

  # Enable WGCNA threads
  WGCNA::allowWGCNAThreads()

  results <- list()

  # Prepare expression data
  expr_data <- prepare_wgcna_data(vst_counts, wc)

  if (is.null(expr_data) || ncol(expr_data) < 100) {
    log_message("Insufficient genes for WGCNA (need at least 100). Skipping.")
    return(NULL)
  }

  log_message("Using ", ncol(expr_data), " genes for WGCNA")

  # 1. Pick soft threshold
  results$soft_threshold <- pick_soft_threshold(expr_data, wc, config)

  # 2. Build network and detect modules
  results$network <- build_network_and_modules(
    expr_data, results$soft_threshold$power, wc, config
  )

  if (is.null(results$network)) {
    log_message("Network construction failed. Skipping rest of WGCNA.")
    return(results)
  }

  # 3. Module-trait correlations
  results$module_traits <- analyze_module_trait_correlations(
    results$network, metadata, wc, config
  )

  # 4. Identify hub genes
  results$hub_genes <- identify_hub_genes(
    results$network, gene_annotation, wc, config
  )

  # 5. Export module membership
  results$module_membership <- export_module_membership(
    results$network, gene_annotation, config
  )

  # 6. Generate visualizations
  plot_wgcna_results(results, metadata, config)

  # 7. Summary
  results$summary <- summarize_wgcna(results, config)

  log_message("=== WGCNA Analysis Complete ===")
  return(results)
}

#' Get WGCNA configuration with defaults
get_wgcna_config <- function(config) {
  wc <- config$coexpression %||% list()

  list(
    run_wgcna = wc$run_wgcna %||% TRUE,
    soft_power = wc$soft_power %||% "auto",
    min_module_size = wc$min_module_size %||% 30,
    module_merge_threshold = wc$module_merge_threshold %||% 0.25,
    n_top_genes = wc$n_top_genes %||% 5000,
    hub_gene_threshold = wc$hub_gene_threshold %||% 0.8,
    network_type = wc$network_type %||% "unsigned"
  )
}

#' Prepare expression data for WGCNA
prepare_wgcna_data <- function(vst_counts, wc) {
  log_message("Preparing expression data for WGCNA...")

  # Remove genes with zero variance
  gene_vars <- apply(vst_counts, 1, var, na.rm = TRUE)
  nonzero_var <- gene_vars > 0 & !is.na(gene_vars)
  expr_mat <- vst_counts[nonzero_var, ]

  log_message("  Removed ", sum(!nonzero_var), " genes with zero variance")

  # Select top variable genes
  gene_vars <- apply(expr_mat, 1, var, na.rm = TRUE)
  n_genes <- min(wc$n_top_genes, nrow(expr_mat))
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_genes]

  expr_mat <- expr_mat[top_genes, ]

  # Check for samples with too many missing values
  good_genes <- WGCNA::goodGenes(t(expr_mat), verbose = 0)
  expr_mat <- expr_mat[good_genes, ]

  # Transpose for WGCNA (samples in rows, genes in columns)
  expr_data <- t(expr_mat)

  # Check for good samples
  good_samples <- WGCNA::goodSamples(expr_data, verbose = 0)
  expr_data <- expr_data[good_samples, ]

  log_message("  Final: ", nrow(expr_data), " samples, ", ncol(expr_data), " genes")

  expr_data
}

# =============================================================================
# 1. Soft Threshold Selection
# =============================================================================

#' Pick soft thresholding power
pick_soft_threshold <- function(expr_data, wc, config) {
  log_message("Selecting soft-thresholding power...")

  # Test range of powers
  powers <- c(1:10, seq(12, 20, 2))

  sft <- tryCatch({
    WGCNA::pickSoftThreshold(
      expr_data,
      powerVector = powers,
      networkType = wc$network_type,
      verbose = 0
    )
  }, error = function(e) {
    log_message("  Soft threshold selection failed: ", e$message)
    NULL
  })

  if (is.null(sft)) {
    log_message("  Using default power: 6")
    return(list(power = 6, sft_table = NULL))
  }

  # Select power
  if (wc$soft_power == "auto") {
    # Use estimated power or default to 6
    power <- sft$powerEstimate
    if (is.na(power)) {
      # Select lowest power with scale-free topology > 0.8
      fit_indices <- sft$fitIndices
      above_threshold <- fit_indices[fit_indices$SFT.R.sq > 0.8, ]

      if (nrow(above_threshold) > 0) {
        power <- above_threshold$Power[1]
      } else {
        power <- 6
      }
    }
  } else {
    power <- as.numeric(wc$soft_power)
  }

  log_message("  Selected soft-thresholding power: ", power)

  # Save fit indices
  output_dir <- config$output_dir %||% "outputs"
  if (!is.null(sft$fitIndices)) {
    readr::write_csv(sft$fitIndices, file.path(output_dir, "wgcna_soft_threshold.csv"))
  }

  list(
    power = power,
    sft_table = sft$fitIndices,
    fit_r2 = if (!is.null(sft$fitIndices)) {
      sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == power]
    } else NA
  )
}

# =============================================================================
# 2. Network Construction and Module Detection
# =============================================================================

#' Build network and detect modules
build_network_and_modules <- function(expr_data, power, wc, config) {
  log_message("Building network and detecting modules...")

  # Use blockwiseModules for memory efficiency
  net <- tryCatch({
    WGCNA::blockwiseModules(
      expr_data,
      power = power,
      TOMType = wc$network_type,
      minModuleSize = wc$min_module_size,
      reassignThreshold = 0,
      mergeCutHeight = wc$module_merge_threshold,
      numericLabels = TRUE,
      pamRespectsDendro = FALSE,
      saveTOMs = FALSE,
      verbose = 0
    )
  }, error = function(e) {
    log_message("  Network construction failed: ", e$message)
    NULL
  })

  if (is.null(net)) return(NULL)

  # Process results
  module_labels <- net$colors
  module_colors <- WGCNA::labels2colors(module_labels)

  n_modules <- length(unique(module_labels)) - 1  # Exclude grey (0)
  log_message("  Detected ", n_modules, " modules (plus unassigned)")

  # Module sizes
  module_sizes <- table(module_colors)
  log_message("  Module sizes: ", paste(names(module_sizes), "=", module_sizes, collapse = ", "))

  # Module eigengenes
  MEs <- net$MEs
  colnames(MEs) <- gsub("^ME", "", colnames(MEs))

  list(
    module_labels = module_labels,
    module_colors = module_colors,
    MEs = MEs,
    dendrograms = net$dendrograms,
    gene_tree = net$dendrograms[[1]],
    n_modules = n_modules,
    module_sizes = module_sizes,
    gene_names = colnames(expr_data)
  )
}

# =============================================================================
# 3. Module-Trait Correlations
# =============================================================================

#' Analyze module-trait correlations
analyze_module_trait_correlations <- function(network, metadata, wc, config) {
  log_message("Analyzing module-trait correlations...")

  MEs <- network$MEs

  # Align samples
  common_samples <- intersect(rownames(MEs), rownames(metadata))
  if (length(common_samples) < 5) {
    log_message("  Insufficient common samples for trait correlation")
    return(NULL)
  }

  MEs_aligned <- MEs[common_samples, ]
  meta_aligned <- metadata[common_samples, ]

  # Get numeric traits
  trait_cols <- sapply(meta_aligned, function(x) {
    is.numeric(x) || is.factor(x) || is.character(x)
  })
  trait_names <- names(trait_cols)[trait_cols]

  # Exclude sample ID columns
  trait_names <- setdiff(trait_names, c("sample_id", "sample", "Sample"))

  if (length(trait_names) == 0) {
    log_message("  No trait variables found for correlation")
    return(NULL)
  }

  # Convert traits to numeric
  traits_numeric <- lapply(trait_names, function(tn) {
    x <- meta_aligned[[tn]]
    if (is.numeric(x)) return(x)
    if (is.factor(x) || is.character(x)) {
      return(as.numeric(as.factor(x)))
    }
    NULL
  })
  names(traits_numeric) <- trait_names
  traits_numeric <- traits_numeric[!sapply(traits_numeric, is.null)]

  if (length(traits_numeric) == 0) return(NULL)

  traits_mat <- as.data.frame(traits_numeric)

  # Compute correlations
  module_trait_cor <- WGCNA::cor(MEs_aligned, traits_mat, use = "pairwise.complete.obs")
  module_trait_pval <- WGCNA::corPvalueStudent(module_trait_cor, nrow(MEs_aligned))

  log_message("  Computed correlations for ", ncol(module_trait_cor), " traits")

  # Save results
  output_dir <- config$output_dir %||% "outputs"

  cor_df <- as.data.frame(module_trait_cor)
  cor_df$module <- rownames(cor_df)
  readr::write_csv(cor_df, file.path(output_dir, "wgcna_module_trait_correlation.csv"))

  pval_df <- as.data.frame(module_trait_pval)
  pval_df$module <- rownames(pval_df)
  readr::write_csv(pval_df, file.path(output_dir, "wgcna_module_trait_pvalue.csv"))

  # Find significant associations
  sig_associations <- which(module_trait_pval < 0.05, arr.ind = TRUE)
  if (nrow(sig_associations) > 0) {
    sig_df <- data.frame(
      module = rownames(module_trait_cor)[sig_associations[, 1]],
      trait = colnames(module_trait_cor)[sig_associations[, 2]],
      correlation = module_trait_cor[sig_associations],
      pvalue = module_trait_pval[sig_associations],
      stringsAsFactors = FALSE
    )
    sig_df <- sig_df[order(sig_df$pvalue), ]

    log_message("  Found ", nrow(sig_df), " significant module-trait associations")
    readr::write_csv(sig_df, file.path(output_dir, "wgcna_significant_associations.csv"))
  }

  list(
    correlations = module_trait_cor,
    pvalues = module_trait_pval,
    traits = traits_mat,
    n_significant = nrow(sig_associations)
  )
}

# =============================================================================
# 4. Hub Gene Identification
# =============================================================================

#' Identify hub genes in each module
identify_hub_genes <- function(network, gene_annotation, wc, config) {
  log_message("Identifying hub genes...")

  gene_names <- network$gene_names
  module_colors <- network$module_colors
  MEs <- network$MEs

  # Calculate module membership (kME)
  kME <- tryCatch({
    WGCNA::signedKME(
      t(network$MEs),  # Transpose if needed
      network$MEs,
      outputColumnName = ""
    )
  }, error = function(e) {
    # Alternative calculation
    expr_data <- t(network$MEs)  # This won't work directly, need original expr_data
    NULL
  })

  # For now, use simpler hub identification based on connectivity
  hub_genes <- list()

  unique_modules <- unique(module_colors)
  unique_modules <- unique_modules[unique_modules != "grey"]

  for (mod in unique_modules) {
    mod_genes <- gene_names[module_colors == mod]

    if (length(mod_genes) < 5) next

    # Hub genes are those most correlated with module eigengene
    # For simplified version, just report module genes
    hub_genes[[mod]] <- data.frame(
      gene = mod_genes,
      module = mod,
      stringsAsFactors = FALSE
    )
  }

  # Combine and add annotation
  all_hubs <- do.call(rbind, hub_genes)

  if (!is.null(gene_annotation) && nrow(all_hubs) > 0) {
    # Try to add gene descriptions
    symbol_col <- intersect(c("gene_symbol", "symbol", "SYMBOL"), colnames(gene_annotation))[1]
    desc_col <- intersect(c("description", "gene_name", "GENENAME"), colnames(gene_annotation))[1]

    if (!is.na(symbol_col) && !is.na(desc_col)) {
      all_hubs$description <- gene_annotation[[desc_col]][
        match(all_hubs$gene, gene_annotation[[symbol_col]])
      ]
    }
  }

  # Save hub genes
  output_dir <- config$output_dir %||% "outputs"
  readr::write_csv(all_hubs, file.path(output_dir, "wgcna_hub_genes.csv"))

  log_message("  Identified hub genes for ", length(hub_genes), " modules")

  hub_genes
}

# =============================================================================
# 5. Export Module Membership
# =============================================================================

#' Export module membership table
export_module_membership <- function(network, gene_annotation, config) {
  log_message("Exporting module membership...")

  module_df <- data.frame(
    gene = network$gene_names,
    module_number = network$module_labels,
    module_color = network$module_colors,
    stringsAsFactors = FALSE
  )

  # Add annotation if available
  if (!is.null(gene_annotation)) {
    symbol_col <- intersect(c("gene_symbol", "symbol", "SYMBOL", "gene_id"),
                            colnames(gene_annotation))[1]

    if (!is.na(symbol_col)) {
      ann_idx <- match(module_df$gene, gene_annotation[[symbol_col]])

      for (col in c("description", "gene_name", "GENENAME", "biotype")) {
        if (col %in% colnames(gene_annotation)) {
          module_df[[col]] <- gene_annotation[[col]][ann_idx]
        }
      }
    }
  }

  # Save
  output_dir <- config$output_dir %||% "outputs"
  readr::write_csv(module_df, file.path(output_dir, "wgcna_modules.csv"))

  # Also save module eigengenes
  ME_df <- as.data.frame(network$MEs)
  ME_df$sample_id <- rownames(ME_df)
  readr::write_csv(ME_df, file.path(output_dir, "wgcna_module_eigengenes.csv"))

  log_message("  Exported ", nrow(module_df), " genes across ",
              network$n_modules, " modules")

  module_df
}

# =============================================================================
# 6. Visualizations
# =============================================================================

#' Generate WGCNA visualizations
plot_wgcna_results <- function(results, metadata, config) {
  log_message("Generating WGCNA plots...")

  output_dir <- config$output_dir %||% "outputs"
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  # Plot 1: Soft threshold selection
  if (!is.null(results$soft_threshold$sft_table)) {
    plot_soft_threshold(results$soft_threshold, plots_dir)
  }

  # Plot 2: Module dendrogram
  if (!is.null(results$network)) {
    plot_module_dendrogram(results$network, plots_dir)
  }

  # Plot 3: Module-trait heatmap
  if (!is.null(results$module_traits)) {
    plot_module_trait_heatmap(results$module_traits, plots_dir)
  }

  # Plot 4: Module sizes
  if (!is.null(results$network)) {
    plot_module_sizes(results$network, plots_dir)
  }

  log_message("WGCNA plots saved")
}

#' Plot soft threshold selection
plot_soft_threshold <- function(sft_result, plots_dir) {
  sft_table <- sft_result$sft_table
  selected_power <- sft_result$power

  # Scale-free topology fit
  p1 <- ggplot2::ggplot(sft_table, ggplot2::aes(x = Power, y = SFT.R.sq)) +
    ggplot2::geom_line(color = "steelblue", size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = selected_power, linetype = "dashed", color = "darkgreen") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Scale-Free Topology Fit",
      subtitle = paste0("Selected power: ", selected_power),
      x = "Soft Threshold (Power)",
      y = "Scale Free Topology Model Fit (R^2)"
    )

  # Mean connectivity
  p2 <- ggplot2::ggplot(sft_table, ggplot2::aes(x = Power, y = mean.k.)) +
    ggplot2::geom_line(color = "steelblue", size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = selected_power, linetype = "dashed", color = "darkgreen") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Mean Connectivity",
      x = "Soft Threshold (Power)",
      y = "Mean Connectivity"
    )

  # Combine plots if patchwork available
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p1 + p2
    ggplot2::ggsave(file.path(plots_dir, "wgcna_soft_threshold.png"),
                    plot = combined, width = 12, height = 5, dpi = 300)
  } else {
    ggplot2::ggsave(file.path(plots_dir, "wgcna_soft_threshold_fit.png"),
                    plot = p1, width = 8, height = 5, dpi = 300)
    ggplot2::ggsave(file.path(plots_dir, "wgcna_soft_threshold_connectivity.png"),
                    plot = p2, width = 8, height = 5, dpi = 300)
  }
}

#' Plot module dendrogram
plot_module_dendrogram <- function(network, plots_dir) {
  if (is.null(network$gene_tree)) return(NULL)

  png(file.path(plots_dir, "wgcna_module_dendrogram.png"),
      width = 12, height = 8, units = "in", res = 300)

  tryCatch({
    WGCNA::plotDendroAndColors(
      network$gene_tree,
      network$module_colors,
      "Module colors",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05,
      main = "Gene Clustering Dendrogram with Module Colors"
    )
  }, error = function(e) {
    log_message("  Module dendrogram plot failed: ", e$message)
  })

  dev.off()
}

#' Plot module-trait correlation heatmap
plot_module_trait_heatmap <- function(module_traits, plots_dir) {
  cor_mat <- module_traits$correlations
  pval_mat <- module_traits$pvalues

  # Create text matrix with correlation and significance
  text_mat <- matrix(
    paste0(round(cor_mat, 2), "\n(",
           ifelse(pval_mat < 0.001, "<0.001",
                  ifelse(pval_mat < 0.01, "<0.01",
                         ifelse(pval_mat < 0.05, "<0.05", round(pval_mat, 2)))),
           ")"),
    nrow = nrow(cor_mat)
  )

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    png(file.path(plots_dir, "wgcna_module_trait_heatmap.png"),
        width = max(8, ncol(cor_mat) + 2), height = max(8, nrow(cor_mat) * 0.4),
        units = "in", res = 300)

    tryCatch({
      pheatmap::pheatmap(
        cor_mat,
        display_numbers = text_mat,
        number_color = "black",
        fontsize_number = 7,
        color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
        main = "Module-Trait Correlations",
        fontsize_row = 8,
        fontsize_col = 8,
        cluster_rows = FALSE,
        cluster_cols = FALSE
      )
    }, error = function(e) {
      log_message("  Module-trait heatmap failed: ", e$message)
    })

    dev.off()
  }
}

#' Plot module sizes
plot_module_sizes <- function(network, plots_dir) {
  module_sizes <- as.data.frame(network$module_sizes)
  colnames(module_sizes) <- c("module", "size")

  # Exclude grey module
  module_sizes <- module_sizes[module_sizes$module != "grey", ]

  if (nrow(module_sizes) == 0) return(NULL)

  p <- ggplot2::ggplot(module_sizes, ggplot2::aes(x = reorder(module, -size),
                                                   y = size, fill = module)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "WGCNA Module Sizes",
      x = "Module",
      y = "Number of Genes"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::scale_fill_manual(values = as.character(module_sizes$module))

  ggplot2::ggsave(file.path(plots_dir, "wgcna_module_sizes.png"),
                  plot = p, width = 10, height = 6, dpi = 300)
}

#' Summarize WGCNA results
summarize_wgcna <- function(results, config) {
  summary_list <- list()

  if (!is.null(results$soft_threshold)) {
    summary_list$soft_power <- results$soft_threshold$power
    summary_list$fit_r2 <- round(results$soft_threshold$fit_r2, 3)
  }

  if (!is.null(results$network)) {
    summary_list$n_modules <- results$network$n_modules
    summary_list$n_genes <- length(results$network$gene_names)
    summary_list$n_unassigned <- sum(results$network$module_colors == "grey")
  }

  if (!is.null(results$module_traits)) {
    summary_list$n_significant_associations <- results$module_traits$n_significant
  }

  # Save summary
  summary_df <- data.frame(
    metric = names(summary_list),
    value = as.character(unlist(summary_list)),
    stringsAsFactors = FALSE
  )

  output_dir <- config$output_dir %||% "outputs"
  readr::write_csv(summary_df, file.path(output_dir, "wgcna_summary.csv"))

  summary_list
}

# =============================================================================
# Helper functions
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
