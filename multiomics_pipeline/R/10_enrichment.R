# =============================================================================
# Multi-Omics Pathway Enrichment
# =============================================================================

#' Run multi-omics pathway enrichment
run_multiomics_enrichment <- function(mae_data, integration_results, config) {
  log_message("=== Running Multi-Omics Pathway Enrichment ===")

  enrich_config <- config$enrichment
  if (!(enrich_config$run_enrichment %||% TRUE)) {
    log_message("Enrichment analysis disabled in config")
    return(NULL)
  }

  harmonized <- mae_data$harmonized_omics
  methods <- enrich_config$methods %||% c("ora")

  results <- list()

  # 1. Per-omics enrichment
  per_omics_results <- list()

  for (omic in names(harmonized)) {
    omic_data <- harmonized[[omic]]
    de_table <- omic_data$de_table %||% omic_data$da_table
    gmt <- omic_data$gmt

    if (is.null(de_table)) {
      log_message("No DE/DA table for ", omic, ", skipping enrichment")
      next
    }

    if (omic == "transcriptomics") {
      # Use gene symbols or gene IDs
      gene_col <- if ("gene_symbol" %in% colnames(de_table)) "gene_symbol" else "gene_id"
      per_omics_results[[omic]] <- run_gene_enrichment(
        de_table, gene_col, gmt, omic, config
      )
    } else if (omic == "metabolomics") {
      # Metabolite enrichment if pathway mapping available
      pathway_mapping <- omic_data$pathway_mapping
      per_omics_results[[omic]] <- run_metabolite_enrichment(
        de_table, pathway_mapping, gmt, config
      )
    }
    # Proteomics can use gene-based enrichment via gene symbols
    else if (omic == "proteomics") {
      feat_anno <- omic_data$feature_annotation
      if (!is.null(feat_anno) && "gene_symbol" %in% colnames(feat_anno)) {
        de_table$gene_symbol <- feat_anno$gene_symbol[match(de_table$feature_id, feat_anno$feature_id)]
        per_omics_results[[omic]] <- run_gene_enrichment(
          de_table, "gene_symbol", NULL, omic, config
        )
      }
    }
  }

  results$per_omics <- per_omics_results

  # 2. Combined multi-omics enrichment
  if (length(per_omics_results) >= 2) {
    results$combined <- combine_enrichment_results(
      per_omics_results,
      enrich_config$combine_method %||% "fisher",
      config
    )
  }

  # 3. Integration-driven enrichment (from MOFA factors)
  if (!is.null(integration_results$mofa)) {
    results$mofa_enrichment <- run_mofa_factor_enrichment(
      integration_results$mofa, harmonized, config
    )
  }

  log_message("=== Enrichment Analysis Complete ===")

  results
}

#' Run gene-based enrichment (ORA)
run_gene_enrichment <- function(de_table, gene_col, gmt, omic_name, config) {
  log_message("Running gene enrichment for ", omic_name, "...")

  enrich_config <- config$enrichment
  pval_thresh <- enrich_config$ora_pvalue %||% 0.05
  min_size <- enrich_config$min_set_size %||% 10
  max_size <- enrich_config$max_set_size %||% 500

  # Get significant genes
  padj_col <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(de_table))[1]
  if (is.na(padj_col)) {
    log_message("No adjusted p-value column found in DE table")
    return(NULL)
  }

  sig_genes <- de_table[[gene_col]][de_table[[padj_col]] < 0.05]
  sig_genes <- unique(sig_genes[!is.na(sig_genes)])
  all_genes <- unique(de_table[[gene_col]][!is.na(de_table[[gene_col]])])

  log_message("  Significant genes: ", length(sig_genes), " / ", length(all_genes))

  if (length(sig_genes) < 5) {
    log_message("  Too few significant genes for enrichment")
    return(NULL)
  }

  # Load gene sets
  gene_sets <- NULL

  # Try custom GMT first
  if (!is.null(gmt)) {
    gene_sets <- gmt
    log_message("  Using custom GMT: ", length(gene_sets), " gene sets")
  }
  # Try clusterProfiler for GO/KEGG
  else if (requireNamespace("clusterProfiler", quietly = TRUE) &&
           requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    log_message("  Running GO enrichment via clusterProfiler...")

    # Try to convert to Entrez IDs
    tryCatch({
      ego <- clusterProfiler::enrichGO(
        gene = sig_genes,
        universe = all_genes,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = pval_thresh,
        minGSSize = min_size,
        maxGSSize = max_size
      )

      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        result <- as.data.frame(ego)
        result$omics <- omic_name
        save_table(result, paste0(omic_name, "_GO_enrichment.csv"), config)

        # Plot
        if (nrow(result) > 0) {
          p <- plot_enrichment_dotplot(result, paste0(omic_name, " GO Enrichment"))
          save_plot(p, paste0(omic_name, "_GO_enrichment"), config, width = 10, height = 8)
        }

        return(list(
          method = "GO_clusterProfiler",
          results = result,
          sig_genes = sig_genes
        ))
      }
    }, error = function(e) {
      log_message("  clusterProfiler enrichment failed: ", e$message)
    })
  }

  # Fallback: simple ORA with custom gene sets
  if (!is.null(gene_sets)) {
    ora_result <- run_simple_ora(sig_genes, all_genes, gene_sets, pval_thresh)
    if (!is.null(ora_result) && nrow(ora_result) > 0) {
      ora_result$omics <- omic_name
      save_table(ora_result, paste0(omic_name, "_ORA_enrichment.csv"), config)

      p <- plot_enrichment_dotplot(ora_result, paste0(omic_name, " ORA Enrichment"))
      save_plot(p, paste0(omic_name, "_ORA_enrichment"), config, width = 10, height = 8)

      return(list(
        method = "ORA",
        results = ora_result,
        sig_genes = sig_genes
      ))
    }
  }

  log_message("  No enrichment results generated")
  return(NULL)
}

#' Simple ORA implementation
run_simple_ora <- function(sig_genes, background, gene_sets, pval_thresh = 0.05) {
  n_bg <- length(background)
  n_sig <- length(sig_genes)

  results <- lapply(names(gene_sets), function(gs_name) {
    gs_genes <- gene_sets[[gs_name]]
    gs_in_bg <- intersect(gs_genes, background)
    n_gs <- length(gs_in_bg)

    if (n_gs < 5) return(NULL)

    overlap <- intersect(sig_genes, gs_in_bg)
    n_overlap <- length(overlap)

    if (n_overlap == 0) return(NULL)

    # Fisher's exact test
    mat <- matrix(c(
      n_overlap,
      n_sig - n_overlap,
      n_gs - n_overlap,
      n_bg - n_gs - n_sig + n_overlap
    ), nrow = 2)

    pval <- fisher.test(mat, alternative = "greater")$p.value

    data.frame(
      term = gs_name,
      overlap = n_overlap,
      term_size = n_gs,
      query_size = n_sig,
      background_size = n_bg,
      pvalue = pval,
      genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results)

  if (is.null(results) || nrow(results) == 0) return(NULL)

  results$padj <- p.adjust(results$pvalue, method = "BH")
  results$fold_enrichment <- (results$overlap / results$query_size) /
                              (results$term_size / results$background_size)

  results <- results[results$padj < pval_thresh, ]
  results <- results[order(results$padj), ]

  results
}

#' Run metabolite enrichment
run_metabolite_enrichment <- function(da_table, pathway_mapping, gmt, config) {
  log_message("Running metabolite enrichment...")

  if (is.null(pathway_mapping) && is.null(gmt)) {
    log_message("  No pathway mapping or GMT for metabolites")
    return(NULL)
  }

  # Get significant metabolites
  padj_col <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(da_table))[1]
  if (is.na(padj_col)) return(NULL)

  sig_features <- da_table$feature_id[da_table[[padj_col]] < 0.05]
  sig_features <- unique(sig_features[!is.na(sig_features)])
  all_features <- unique(da_table$feature_id)

  log_message("  Significant metabolites: ", length(sig_features), " / ", length(all_features))

  if (length(sig_features) < 3) {
    log_message("  Too few significant metabolites")
    return(NULL)
  }

  # Use GMT if available
  if (!is.null(gmt)) {
    ora_result <- run_simple_ora(sig_features, all_features, gmt, 0.05)
    if (!is.null(ora_result) && nrow(ora_result) > 0) {
      ora_result$omics <- "metabolomics"
      save_table(ora_result, "metabolomics_ORA_enrichment.csv", config)

      p <- plot_enrichment_dotplot(ora_result, "Metabolomics ORA Enrichment")
      save_plot(p, "metabolomics_ORA_enrichment", config, width = 10, height = 8)

      return(list(method = "ORA", results = ora_result, sig_features = sig_features))
    }
  }

  # Use pathway mapping if available
  if (!is.null(pathway_mapping)) {
    # Convert to gene set format
    pathway_col <- grep("pathway", colnames(pathway_mapping), ignore.case = TRUE, value = TRUE)[1]
    id_col <- colnames(pathway_mapping)[1]

    if (!is.na(pathway_col)) {
      pathway_list <- split(pathway_mapping[[id_col]], pathway_mapping[[pathway_col]])
      ora_result <- run_simple_ora(sig_features, all_features, pathway_list, 0.05)

      if (!is.null(ora_result) && nrow(ora_result) > 0) {
        ora_result$omics <- "metabolomics"
        save_table(ora_result, "metabolomics_pathway_enrichment.csv", config)
        return(list(method = "pathway", results = ora_result, sig_features = sig_features))
      }
    }
  }

  return(NULL)
}

#' Combine enrichment results across omics
combine_enrichment_results <- function(per_omics_results, method = "fisher", config) {
  log_message("Combining enrichment results across omics (method: ", method, ")...")

  # Collect all enrichment results
  all_results <- list()
  for (omic in names(per_omics_results)) {
    res <- per_omics_results[[omic]]
    if (!is.null(res) && !is.null(res$results)) {
      all_results[[omic]] <- res$results
    }
  }

  if (length(all_results) < 2) {
    log_message("Need at least 2 omics with enrichment results")
    return(NULL)
  }

  # Find common terms
  all_terms <- Reduce(intersect, lapply(all_results, function(x) x$term))

  if (length(all_terms) < 1) {
    log_message("No common enriched terms across omics")
    return(NULL)
  }

  log_message("Found ", length(all_terms), " terms enriched in multiple omics")

  # Combine p-values
  combined_df <- data.frame(term = all_terms, stringsAsFactors = FALSE)

  for (omic in names(all_results)) {
    res <- all_results[[omic]]
    idx <- match(all_terms, res$term)
    combined_df[[paste0("pvalue_", omic)]] <- res$pvalue[idx]
    combined_df[[paste0("padj_", omic)]] <- res$padj[idx]
  }

  # Combine p-values
  pval_cols <- grep("^pvalue_", colnames(combined_df), value = TRUE)
  pvals_matrix <- as.matrix(combined_df[, pval_cols])

  combined_df$combined_pvalue <- apply(pvals_matrix, 1, function(pvals) {
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) < 2) return(NA)

    if (method == "fisher") {
      fisher_combine_pvalues(pvals)
    } else if (method == "stouffer") {
      stouffer_combine_pvalues(pvals)
    } else {
      min(pvals)  # Minimum p-value
    }
  })

  combined_df$combined_padj <- p.adjust(combined_df$combined_pvalue, method = "BH")
  combined_df$n_omics <- rowSums(!is.na(pvals_matrix))

  combined_df <- combined_df[order(combined_df$combined_padj), ]

  save_table(combined_df, "combined_enrichment.csv", config)

  # Plot
  if (sum(combined_df$combined_padj < 0.05, na.rm = TRUE) > 0) {
    sig_combined <- combined_df[combined_df$combined_padj < 0.05, ]
    sig_combined <- head(sig_combined, 20)

    p <- ggplot2::ggplot(sig_combined, ggplot2::aes(
      x = -log10(combined_pvalue),
      y = reorder(term, -log10(combined_pvalue)),
      size = n_omics
    )) +
      ggplot2::geom_point(color = "steelblue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Combined Multi-Omics Enrichment",
        x = "-log10(Combined P-value)",
        y = NULL,
        size = "# Omics"
      )

    save_plot(p, "combined_enrichment", config, width = 10, height = 8)
  }

  list(
    combined = combined_df,
    method = method
  )
}

#' Run enrichment on MOFA factors
run_mofa_factor_enrichment <- function(mofa_results, harmonized, config) {
  log_message("Running enrichment on MOFA factor loadings...")

  weights <- mofa_results$results$weights
  top_features <- mofa_results$results$top_features

  results <- list()

  # For transcriptomics view, run enrichment on top weighted genes per factor
  if ("transcriptomics" %in% names(weights)) {
    rna_weights <- weights$transcriptomics
    rna_anno <- harmonized$transcriptomics$feature_annotation

    for (k in seq_len(ncol(rna_weights))) {
      factor_name <- colnames(rna_weights)[k]

      # Top positive and negative genes
      w <- rna_weights[, k]
      top_pos <- names(sort(w, decreasing = TRUE))[1:100]
      top_neg <- names(sort(w, decreasing = FALSE))[1:100]

      # Map to gene symbols if possible
      if (!is.null(rna_anno) && "gene_symbol" %in% colnames(rna_anno)) {
        idx_pos <- match(top_pos, rna_anno$feature_id)
        top_pos_genes <- rna_anno$gene_symbol[idx_pos]
        top_pos_genes <- top_pos_genes[!is.na(top_pos_genes)]

        idx_neg <- match(top_neg, rna_anno$feature_id)
        top_neg_genes <- rna_anno$gene_symbol[idx_neg]
        top_neg_genes <- top_neg_genes[!is.na(top_neg_genes)]

        # Run enrichment if clusterProfiler available
        if (requireNamespace("clusterProfiler", quietly = TRUE) &&
            requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
            length(top_pos_genes) > 10) {

          tryCatch({
            ego <- clusterProfiler::enrichGO(
              gene = top_pos_genes,
              OrgDb = org.Hs.eg.db::org.Hs.eg.db,
              keyType = "SYMBOL",
              ont = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff = 0.05
            )

            if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
              result <- as.data.frame(ego)
              result$factor <- factor_name
              result$direction <- "positive"
              save_table(result, paste0("mofa_", factor_name, "_pos_enrichment.csv"), config)
              results[[paste0(factor_name, "_pos")]] <- result
            }
          }, error = function(e) NULL)
        }
      }
    }
  }

  if (length(results) > 0) {
    log_message("MOFA factor enrichment complete: ", length(results), " enrichment sets")
  }

  results
}

#' Plot enrichment dotplot
plot_enrichment_dotplot <- function(enrich_df, title, n_terms = 20) {
  if (nrow(enrich_df) == 0) {
    return(ggplot2::ggplot() + ggplot2::ggtitle("No enriched terms"))
  }

  # Select top terms
  plot_df <- head(enrich_df, n_terms)

  # Determine which columns exist
  size_col <- if ("overlap" %in% colnames(plot_df)) "overlap" else "Count"
  pval_col <- if ("padj" %in% colnames(plot_df)) "padj" else "p.adjust"

  if (!size_col %in% colnames(plot_df)) plot_df$overlap <- 1
  if (!pval_col %in% colnames(plot_df)) plot_df$padj <- plot_df$pvalue

  plot_df$term <- factor(plot_df$term, levels = rev(plot_df$term))

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = -log10(padj),
    y = term,
    size = overlap,
    color = -log10(padj)
  )) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = "-log10(Adjusted P-value)",
      y = NULL,
      size = "Gene Count",
      color = "-log10(padj)"
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8)
    )
}
