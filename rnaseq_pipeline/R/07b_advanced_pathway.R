# =============================================================================
# Advanced Pathway Analysis
# =============================================================================
# This script provides enhanced pathway analysis beyond basic FGSEA/ORA:
#   1. GSVA - Gene Set Variation Analysis (sample-level pathway scores)
#   2. MSigDB integration (Hallmark, C2-curated, C5-GO, C7-immunologic)
#   3. Reactome pathway enrichment
#   4. Pathway activity clustering and heatmaps
# =============================================================================

#' Main function: Run advanced pathway analysis
#' @param normalized_counts Normalized expression matrix
#' @param de_results_annotated DE results with gene symbols
#' @param gene_annotation Gene annotation table
#' @param metadata Sample metadata
#' @param config Pipeline configuration
#' @return List with GSVA results and extended enrichment
run_advanced_pathway_analysis <- function(normalized_counts, de_results_annotated,
                                          gene_annotation, metadata, config) {
  log_message("=== Running Advanced Pathway Analysis ===")

  # Wrap entire function in tryCatch to prevent pipeline failure
 tryCatch({
    ap <- get_advanced_pathway_config(config)

    if (!ap$run_advanced_pathway) {
      log_message("Advanced pathway analysis disabled. Skipping.")
      return(NULL)
    }

    results <- list()

    # 1. Run GSVA
    if (ap$run_gsva) {
      results$gsva <- tryCatch(
        run_gsva_analysis(normalized_counts, metadata, ap, config),
        error = function(e) {
          log_message("GSVA analysis failed: ", e$message)
          NULL
        }
      )
    }

    # 2. Extended database enrichment (only if we have DE results)
    if (!is.null(de_results_annotated) && length(de_results_annotated) > 0) {
      results$extended_enrichment <- tryCatch(
        run_extended_enrichment(de_results_annotated, gene_annotation, ap, config),
        error = function(e) {
          log_message("Extended enrichment failed: ", e$message)
          list()
        }
      )
    } else {
      log_message("No DE results available for extended enrichment. Skipping.")
      results$extended_enrichment <- list()
    }

    # 3. Pathway activity comparison across contrasts
    if (length(de_results_annotated) > 1) {
      results$pathway_comparison <- tryCatch(
        compare_pathway_across_contrasts(de_results_annotated, ap, config),
        error = function(e) NULL
      )
    }

    # 4. Generate visualizations
    tryCatch(
      plot_advanced_pathway_results(results, metadata, config),
      error = function(e) {
        log_message("Advanced pathway plots failed: ", e$message)
      }
    )

    # 5. Summary
    results$summary <- tryCatch(
      summarize_advanced_pathway(results, config),
      error = function(e) list()
    )

    log_message("=== Advanced Pathway Analysis Complete ===")
    return(results)

  }, error = function(e) {
    log_message("Advanced pathway analysis failed: ", e$message)
    return(NULL)
  })
}

#' Get advanced pathway configuration with defaults
get_advanced_pathway_config <- function(config) {
  ap <- config$advanced_pathway %||% list()

  list(
    run_advanced_pathway = ap$run_advanced_pathway %||% TRUE,
    run_gsva = ap$run_gsva %||% TRUE,
    gsva_method = ap$gsva_method %||% "gsva",  # "gsva", "ssgsea", "zscore", "plage"
    msigdb_categories = ap$msigdb_categories %||% c("H", "C2", "C5"),
    organism = config$organism %||% "human",
    min_gene_set_size = ap$min_gene_set_size %||% 10,
    max_gene_set_size = ap$max_gene_set_size %||% 500,
    fdr_threshold = ap$fdr_threshold %||% 0.05,
    n_top_pathways = ap$n_top_pathways %||% 50
  )
}

# =============================================================================
# 1. GSVA Analysis
# =============================================================================

#' Run GSVA analysis
#' @param normalized_counts Normalized expression matrix
#' @param metadata Sample metadata
#' @param ap Advanced pathway config
#' @param config Full config
run_gsva_analysis <- function(normalized_counts, metadata, ap, config) {
  log_message("Running GSVA analysis...")

  # Check if GSVA is available
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    log_message("  GSVA package not available. Skipping.")
    return(NULL)
  }

  # Load gene sets
  gene_sets <- load_msigdb_gene_sets(ap$msigdb_categories, ap$organism, ap)

  if (is.null(gene_sets) || length(gene_sets) == 0) {
    log_message("  No gene sets loaded. Skipping GSVA.")
    return(NULL)
  }

  log_message("  Loaded ", length(gene_sets), " gene sets from MSigDB")

  # Ensure we have gene symbols as rownames
  expr_mat <- as.matrix(normalized_counts)

  # Convert to gene symbols if needed (MSigDB uses gene symbols)
  expr_mat <- convert_to_gene_symbols(expr_mat, config)

  # Check overlap with gene sets
  all_gs_genes <- unique(unlist(gene_sets))
  overlap <- length(intersect(rownames(expr_mat), all_gs_genes))
  log_message("  Gene overlap with MSigDB: ", overlap, " / ", length(all_gs_genes), " genes")

  if (overlap < 100) {
    log_message("  WARNING: Low gene overlap. Check if expression matrix uses gene symbols.")
    log_message("  Sample rownames: ", paste(head(rownames(expr_mat), 5), collapse = ", "))
    log_message("  Sample gene set genes: ", paste(head(all_gs_genes, 5), collapse = ", "))
  }

  # Run GSVA (using new API for Bioconductor 3.18+)
  gsva_scores <- tryCatch({
    log_message("  Running GSVA (method: ", ap$gsva_method, ")...")

    # Create parameter object using new GSVA API
    param <- GSVA::gsvaParam(
      exprData = expr_mat,
      geneSets = gene_sets,
      minSize = ap$min_gene_set_size,
      maxSize = ap$max_gene_set_size
    )

    GSVA::gsva(param, verbose = FALSE)
  }, error = function(e) {
    log_message("  GSVA failed: ", e$message)
    NULL
  })

  if (is.null(gsva_scores)) return(NULL)

  log_message("  GSVA computed ", nrow(gsva_scores), " pathway scores")

  # Test pathway associations with condition
  pathway_stats <- NULL
  condition_col <- config$group_col %||% "condition"

  if (condition_col %in% colnames(metadata)) {
    pathway_stats <- test_gsva_associations(
      gsva_scores, metadata, condition_col, config
    )
  }

  # Save GSVA scores
  output_dir <- config$output_dir %||% "outputs"
  gsva_df <- as.data.frame(gsva_scores)
  gsva_df <- tibble::rownames_to_column(gsva_df, "pathway")
  readr::write_csv(gsva_df, file.path(output_dir, "gsva_pathway_scores.csv"))

  if (!is.null(pathway_stats)) {
    readr::write_csv(pathway_stats, file.path(output_dir, "gsva_pathway_statistics.csv"))
  }

  list(
    scores = gsva_scores,
    statistics = pathway_stats,
    n_pathways = nrow(gsva_scores),
    method = ap$gsva_method
  )
}

#' Load MSigDB gene sets
load_msigdb_gene_sets <- function(categories, organism, ap) {
  # Check if msigdbr is available
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    log_message("  msigdbr package not available. Using built-in gene sets if available.")
    return(NULL)
  }

  # Map organism name
  species <- switch(tolower(organism),
    "human" = "Homo sapiens",
    "mouse" = "Mus musculus",
    "rat" = "Rattus norvegicus",
    "Homo sapiens"  # default
  )

  all_sets <- list()

  for (cat in categories) {
    cat_sets <- tryCatch({
      if (cat == "H") {
        # Hallmark
        msigdbr::msigdbr(species = species, category = "H")
      } else if (cat == "C2") {
        # Curated gene sets (canonical pathways)
        msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP")
      } else if (cat == "C5") {
        # GO gene sets (BP only for manageable size)
        msigdbr::msigdbr(species = species, category = "C5", subcategory = "GO:BP")
      } else if (cat == "C7") {
        # Immunologic signatures
        msigdbr::msigdbr(species = species, category = "C7", subcategory = "IMMUNESIGDB")
      } else {
        msigdbr::msigdbr(species = species, category = cat)
      }
    }, error = function(e) {
      log_message("  Failed to load MSigDB category ", cat, ": ", e$message)
      NULL
    })

    if (!is.null(cat_sets) && nrow(cat_sets) > 0) {
      # Convert to list format
      sets <- split(cat_sets$gene_symbol, cat_sets$gs_name)
      all_sets <- c(all_sets, sets)
      log_message("  Loaded ", length(sets), " gene sets from category ", cat)
    }
  }

  # Filter by size
  all_sets <- all_sets[sapply(all_sets, length) >= ap$min_gene_set_size &
                        sapply(all_sets, length) <= ap$max_gene_set_size]

  return(all_sets)
}

#' Test GSVA pathway associations with condition
test_gsva_associations <- function(gsva_scores, metadata, condition_col, config) {
  log_message("  Testing pathway associations with condition...")

  common_samples <- intersect(colnames(gsva_scores), rownames(metadata))
  if (length(common_samples) < 5) return(NULL)

  gsva_sub <- gsva_scores[, common_samples, drop = FALSE]
  conditions <- as.factor(metadata[common_samples, condition_col])

  # Test each pathway
  pathway_names <- rownames(gsva_sub)

  results <- lapply(pathway_names, function(pn) {
    scores <- gsva_sub[pn, ]

    # Test based on number of groups
    if (length(levels(conditions)) == 2) {
      # t-test
      test_result <- tryCatch({
        t.test(scores ~ conditions)
      }, error = function(e) NULL)

      if (!is.null(test_result)) {
        return(data.frame(
          pathway = pn,
          test = "t.test",
          statistic = test_result$statistic,
          pvalue = test_result$p.value,
          mean_diff = diff(tapply(scores, conditions, mean)),
          stringsAsFactors = FALSE
        ))
      }
    } else {
      # ANOVA
      test_result <- tryCatch({
        aov_summary <- summary(aov(scores ~ conditions))[[1]]
        list(
          statistic = aov_summary[1, "F value"],
          pvalue = aov_summary[1, "Pr(>F)"]
        )
      }, error = function(e) NULL)

      if (!is.null(test_result)) {
        return(data.frame(
          pathway = pn,
          test = "anova",
          statistic = test_result$statistic,
          pvalue = test_result$pvalue,
          mean_diff = NA,
          stringsAsFactors = FALSE
        ))
      }
    }

    NULL
  })

  results_df <- do.call(rbind, results)
  if (is.null(results_df) || nrow(results_df) == 0) return(NULL)

  results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
  results_df <- results_df[order(results_df$pvalue), ]

  n_sig <- sum(results_df$padj < 0.05, na.rm = TRUE)
  log_message("  Found ", n_sig, " significantly associated pathways")

  results_df
}

# =============================================================================
# 2. Extended Database Enrichment
# =============================================================================

#' Run enrichment on extended databases
run_extended_enrichment <- function(de_results_annotated, gene_annotation, ap, config) {
  log_message("Running extended database enrichment...")

  results <- list()

  # Check if we have any DE results to process
  if (is.null(de_results_annotated) || length(de_results_annotated) == 0) {
    log_message("  No DE results available. Skipping extended enrichment.")
    return(results)
  }

  # For each contrast
  for (contrast_name in names(de_results_annotated)) {
    de_result <- de_results_annotated[[contrast_name]]

    # Skip if de_result is NULL or empty
    if (is.null(de_result) || nrow(de_result) == 0) {
      log_message("  Contrast ", contrast_name, " has no results. Skipping.")
      next
    }

    log_message("  Processing contrast: ", contrast_name)

    # Get significant genes
    sig_genes <- get_significant_genes(de_result, ap)

    if (length(sig_genes$all) < 10) {
      log_message("  Too few significant genes. Skipping.")
      next
    }

    contrast_results <- list()

    # Reactome enrichment
    if (requireNamespace("ReactomePA", quietly = TRUE)) {
      contrast_results$reactome <- run_reactome_enrichment(sig_genes, ap, config)
    }

    # WikiPathways (if clusterProfiler supports)
    contrast_results$wiki <- run_wiki_enrichment(sig_genes, ap, config)

    # MSigDB category-specific enrichment
    contrast_results$msigdb <- run_msigdb_enrichment(sig_genes, ap, config)

    results[[contrast_name]] <- contrast_results
  }

  # Save results
  save_extended_enrichment(results, config)

  results
}

#' Get significant genes from DE results
get_significant_genes <- function(de_result, ap) {
  # Find padj column
  padj_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(de_result))[1]
  lfc_col <- intersect(c("log2FoldChange", "logFC", "log2FC"), colnames(de_result))[1]
  gene_col <- intersect(c("gene_symbol", "symbol", "SYMBOL", "gene_name"), colnames(de_result))[1]

  if (is.na(padj_col) || is.na(gene_col)) {
    return(list(all = character(0), up = character(0), down = character(0)))
  }

  # Handle NA values safely - use FALSE for any NA comparisons
  padj_values <- de_result[[padj_col]]
  sig_idx <- !is.na(padj_values) & padj_values < ap$fdr_threshold

  all_sig <- de_result[[gene_col]][sig_idx]
  all_sig <- all_sig[!is.na(all_sig) & all_sig != ""]

  # Up and down regulated
  if (!is.na(lfc_col)) {
    lfc_values <- de_result[[lfc_col]]
    # Handle NA in lfc_values as well
    up_idx <- sig_idx & !is.na(lfc_values) & lfc_values > 0
    down_idx <- sig_idx & !is.na(lfc_values) & lfc_values < 0

    up_genes <- de_result[[gene_col]][up_idx]
    down_genes <- de_result[[gene_col]][down_idx]
  } else {
    up_genes <- character(0)
    down_genes <- character(0)
  }

  list(
    all = all_sig,
    up = up_genes[!is.na(up_genes)],
    down = down_genes[!is.na(down_genes)]
  )
}

#' Run Reactome enrichment
run_reactome_enrichment <- function(sig_genes, ap, config) {
  log_message("    Running Reactome enrichment...")

  # Convert to Entrez IDs
  entrez_ids <- tryCatch({
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = sig_genes$all,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      )
    } else {
      NULL
    }
  }, error = function(e) NULL)

  if (is.null(entrez_ids)) return(NULL)

  entrez_ids <- entrez_ids[!is.na(entrez_ids)]

  if (length(entrez_ids) < 10) return(NULL)

  # Run Reactome enrichment
  reactome_result <- tryCatch({
    ReactomePA::enrichPathway(
      gene = entrez_ids,
      organism = "human",
      pvalueCutoff = ap$fdr_threshold,
      pAdjustMethod = "BH",
      readable = TRUE
    )
  }, error = function(e) {
    log_message("    Reactome enrichment failed: ", e$message)
    NULL
  })

  if (is.null(reactome_result)) return(NULL)

  result_df <- as.data.frame(reactome_result)
  log_message("    Found ", nrow(result_df), " enriched Reactome pathways")

  result_df
}

#' Run WikiPathways enrichment
run_wiki_enrichment <- function(sig_genes, ap, config) {
  # WikiPathways via clusterProfiler if available
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) return(NULL)

  # Note: WikiPathways support may require additional setup
  # For now, return NULL - can be implemented with enricher() and WikiPathways GMT
  NULL
}

#' Run MSigDB category-specific enrichment
run_msigdb_enrichment <- function(sig_genes, ap, config) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) return(NULL)
  if (!requireNamespace("msigdbr", quietly = TRUE)) return(NULL)

  results <- list()

  # Hallmark enrichment
  hallmark_sets <- tryCatch({
    msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  }, error = function(e) NULL)

  if (!is.null(hallmark_sets) && nrow(hallmark_sets) > 0) {
    h_result <- tryCatch({
      clusterProfiler::enricher(
        gene = sig_genes$all,
        TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],
        pvalueCutoff = ap$fdr_threshold
      )
    }, error = function(e) NULL)

    if (!is.null(h_result) && nrow(as.data.frame(h_result)) > 0) {
      results$hallmark <- as.data.frame(h_result)
      log_message("    Found ", nrow(results$hallmark), " enriched Hallmark sets")
    }
  }

  results
}

#' Save extended enrichment results
save_extended_enrichment <- function(results, config) {
  output_dir <- config$output_dir %||% "outputs"
  pathway_dir <- file.path(output_dir, "pathway_results")
  if (!dir.exists(pathway_dir)) dir.create(pathway_dir, recursive = TRUE)

  for (contrast_name in names(results)) {
    contrast_results <- results[[contrast_name]]

    for (db_name in names(contrast_results)) {
      if (!is.null(contrast_results[[db_name]]) &&
          nrow(contrast_results[[db_name]]) > 0) {
        filename <- paste0(contrast_name, "_", db_name, "_enrichment.csv")
        readr::write_csv(contrast_results[[db_name]],
                         file.path(pathway_dir, filename))
      }
    }
  }
}

# =============================================================================
# 3. Pathway Comparison Across Contrasts
# =============================================================================

#' Compare pathways across multiple contrasts
compare_pathway_across_contrasts <- function(de_results_annotated, ap, config) {
  log_message("Comparing pathways across contrasts...")

  # This requires enrichment results for each contrast
  # Simple implementation: compare significant pathways

  NULL  # Placeholder - full implementation would track pathway overlap
}

# =============================================================================
# 4. Visualizations
# =============================================================================

#' Generate advanced pathway visualizations
plot_advanced_pathway_results <- function(results, metadata, config) {
  log_message("Generating advanced pathway plots...")

  output_dir <- config$output_dir %||% "outputs"
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  # Plot 1: GSVA heatmap
  if (!is.null(results$gsva) && !is.null(results$gsva$scores)) {
    plot_gsva_heatmap(results$gsva, metadata, config)
  }

  # Plot 2: GSVA top pathways
  if (!is.null(results$gsva$statistics)) {
    plot_gsva_top_pathways(results$gsva, config)
  }

  log_message("Advanced pathway plots saved")
}

#' Plot GSVA heatmap
plot_gsva_heatmap <- function(gsva_result, metadata, config) {
  gsva_scores <- gsva_result$scores

  # Select top variable pathways
  pathway_vars <- apply(gsva_scores, 1, var, na.rm = TRUE)
  top_pathways <- names(sort(pathway_vars, decreasing = TRUE))[1:min(50, length(pathway_vars))]

  mat_plot <- gsva_scores[top_pathways, , drop = FALSE]

  # Use pheatmap if available
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    output_dir <- config$output_dir %||% "outputs"

    # Prepare annotation
    condition_col <- config$group_col %||% "condition"
    ann_col <- NULL

    if (condition_col %in% colnames(metadata)) {
      common_samples <- intersect(colnames(mat_plot), rownames(metadata))
      ann_col <- data.frame(
        Condition = metadata[common_samples, condition_col],
        row.names = common_samples
      )
      mat_plot <- mat_plot[, common_samples, drop = FALSE]
    }

    png(file.path(output_dir, "plots", "gsva_pathway_heatmap.png"),
        width = 12, height = 10, units = "in", res = 300)

    tryCatch({
      pheatmap::pheatmap(
        mat_plot,
        scale = "row",
        clustering_method = "ward.D2",
        show_colnames = ncol(mat_plot) < 50,
        annotation_col = ann_col,
        main = paste0("GSVA Pathway Activity (", gsva_result$method, ")"),
        color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
        fontsize_row = 6
      )
    }, error = function(e) {
      log_message("  GSVA heatmap failed: ", e$message)
    })

    dev.off()
  }
}

#' Plot top significant GSVA pathways
plot_gsva_top_pathways <- function(gsva_result, config) {
  stats_df <- gsva_result$statistics

  # Filter to significant and get top
  sig_pathways <- stats_df[stats_df$padj < 0.1 & !is.na(stats_df$padj), ]

  if (nrow(sig_pathways) == 0) {
    sig_pathways <- head(stats_df[order(stats_df$pvalue), ], 20)
  } else {
    sig_pathways <- head(sig_pathways[order(sig_pathways$pvalue), ], 20)
  }

  if (nrow(sig_pathways) == 0) return(NULL)

  # Shorten pathway names for display
  sig_pathways$pathway_short <- substr(sig_pathways$pathway, 1, 50)

  p <- ggplot2::ggplot(sig_pathways, ggplot2::aes(x = reorder(pathway_short, -log10(pvalue)),
                                                   y = -log10(pvalue))) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Top GSVA Pathway Associations",
      subtitle = "Association with condition",
      x = "Pathway",
      y = "-log10(p-value)"
    ) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))

  output_dir <- config$output_dir %||% "outputs"
  ggplot2::ggsave(file.path(output_dir, "plots", "gsva_top_pathways.png"),
                  plot = p, width = 10, height = 8, dpi = 300)
}

#' Summarize advanced pathway results
summarize_advanced_pathway <- function(results, config) {
  summary_list <- list()

  if (!is.null(results$gsva)) {
    summary_list$gsva_pathways <- results$gsva$n_pathways
    summary_list$gsva_method <- results$gsva$method

    if (!is.null(results$gsva$statistics)) {
      summary_list$gsva_significant <- sum(results$gsva$statistics$padj < 0.05, na.rm = TRUE)
    }
  }

  if (!is.null(results$extended_enrichment)) {
    # Count enriched pathways across databases
    n_reactome <- sum(sapply(results$extended_enrichment, function(x) {
      if (!is.null(x$reactome)) nrow(x$reactome) else 0
    }))
    summary_list$reactome_enriched <- n_reactome
  }

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

#' Convert expression matrix rownames to gene symbols
#' @param expr_mat Expression matrix with any gene identifiers
#' @param config Pipeline configuration
#' @return Expression matrix with gene symbols as rownames
convert_to_gene_symbols <- function(expr_mat, config) {
  current_ids <- rownames(expr_mat)

  # Check if already gene symbols (heuristic: most don't start with ENSG)
  ensembl_pattern <- sum(grepl("^ENSG[0-9]+", current_ids))
  entrez_pattern <- sum(grepl("^[0-9]+$", current_ids))

  if (ensembl_pattern < 0.5 * length(current_ids) && entrez_pattern < 0.5 * length(current_ids)) {
    # Likely already gene symbols
    log_message("  Gene identifiers appear to be symbols already")
    return(expr_mat)
  }

  log_message("  Converting gene identifiers to symbols...")

  # Try org.Hs.eg.db for human
  organism <- config$organism %||% "human"

  if (organism == "human" && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    tryCatch({
      if (ensembl_pattern > entrez_pattern) {
        # Ensembl IDs - strip version numbers first
        clean_ids <- gsub("\\.[0-9]+$", "", current_ids)
        mapping <- AnnotationDbi::select(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = clean_ids,
          keytype = "ENSEMBL",
          columns = "SYMBOL"
        )
      } else {
        # Entrez IDs
        mapping <- AnnotationDbi::select(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = current_ids,
          keytype = "ENTREZID",
          columns = "SYMBOL"
        )
      }

      # Create lookup
      if (ensembl_pattern > entrez_pattern) {
        id_to_symbol <- setNames(mapping$SYMBOL, mapping$ENSEMBL)
        clean_ids <- gsub("\\.[0-9]+$", "", current_ids)
        new_symbols <- id_to_symbol[clean_ids]
      } else {
        id_to_symbol <- setNames(mapping$SYMBOL, mapping$ENTREZID)
        new_symbols <- id_to_symbol[current_ids]
      }

      # Replace NA with original ID
      new_symbols[is.na(new_symbols)] <- current_ids[is.na(new_symbols)]

      # Handle duplicates by keeping first occurrence
      dup_mask <- duplicated(new_symbols)
      if (any(dup_mask)) {
        log_message("  Removing ", sum(dup_mask), " duplicate gene symbols")
        expr_mat <- expr_mat[!dup_mask, , drop = FALSE]
        new_symbols <- new_symbols[!dup_mask]
      }

      rownames(expr_mat) <- new_symbols
      log_message("  Converted ", sum(!is.na(id_to_symbol[clean_ids])), " / ", length(current_ids), " IDs to symbols")

    }, error = function(e) {
      log_message("  Gene symbol conversion failed: ", e$message)
    })
  } else if (organism == "mouse" && requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    tryCatch({
      if (ensembl_pattern > entrez_pattern) {
        clean_ids <- gsub("\\.[0-9]+$", "", current_ids)
        mapping <- AnnotationDbi::select(
          org.Mm.eg.db::org.Mm.eg.db,
          keys = clean_ids,
          keytype = "ENSEMBL",
          columns = "SYMBOL"
        )
        id_to_symbol <- setNames(mapping$SYMBOL, mapping$ENSEMBL)
        new_symbols <- id_to_symbol[clean_ids]
      } else {
        mapping <- AnnotationDbi::select(
          org.Mm.eg.db::org.Mm.eg.db,
          keys = current_ids,
          keytype = "ENTREZID",
          columns = "SYMBOL"
        )
        id_to_symbol <- setNames(mapping$SYMBOL, mapping$ENTREZID)
        new_symbols <- id_to_symbol[current_ids]
      }

      new_symbols[is.na(new_symbols)] <- current_ids[is.na(new_symbols)]

      dup_mask <- duplicated(new_symbols)
      if (any(dup_mask)) {
        log_message("  Removing ", sum(dup_mask), " duplicate gene symbols")
        expr_mat <- expr_mat[!dup_mask, , drop = FALSE]
        new_symbols <- new_symbols[!dup_mask]
      }

      rownames(expr_mat) <- new_symbols

    }, error = function(e) {
      log_message("  Gene symbol conversion failed: ", e$message)
    })
  }

  return(expr_mat)
}
