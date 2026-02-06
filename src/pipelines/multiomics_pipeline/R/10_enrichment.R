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
    # Proteomics can use gene-based enrichment via gene symbols or Wormbase IDs
    else if (omic == "proteomics") {
      # Identify gene identifier column
      gene_cols <- c("gene_symbol", "Wormbase_id", "gene_id")
      target_col <- NULL

      # 1. Check if already in DE table
      for (col in gene_cols) {
        if (col %in% colnames(de_table)) {
          target_col <- col
          break
        }
      }

      # 2. If not, try to map from feature annotation
      if (is.null(target_col)) {
        feat_anno <- omic_data$feature_annotation
        if (!is.null(feat_anno)) {
          for (col in gene_cols) {
            if (col %in% colnames(feat_anno)) {
              de_table[[col]] <- feat_anno[[col]][match(de_table$feature_id, feat_anno$feature_id)]
              target_col <- col
              break
            }
          }
        }
      }

      if (!is.null(target_col)) {
        log_message("  Using ", target_col, " for Proteomics enrichment")
        per_omics_results[[omic]] <- run_gene_enrichment(
          de_table, target_col, NULL, omic, config
        )
      } else {
        log_message("  No gene identifier (Symbol/WormbaseID) found for Proteomics enrichment")
      }
    }
  }

  results$per_omics <- per_omics_results

  # 1.5 Plot RNA vs Proteomics enrichment scatter plot
  if ("transcriptomics" %in% names(per_omics_results) &&
    "proteomics" %in% names(per_omics_results)) {
    plot_rna_protein_enrichment_scatter(
      per_omics_results$transcriptomics,
      per_omics_results$proteomics,
      config
    )
  }

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

#' Run gene-based enrichment (ORA) with multiple collections support
run_gene_enrichment <- function(de_table, gene_col, gmt, omic_name, config) {
  log_message("Running gene enrichment for ", omic_name, "...")

  enrich_config <- config$enrichment
  pval_thresh <- enrich_config$ora_pvalue %||% 0.05
  min_size <- enrich_config$min_set_size %||% 10
  max_size <- enrich_config$max_set_size %||% 500
  top_n_terms <- enrich_config$top_n_terms %||% c(10, 20, 50)

  # Get significant genes
  padj_col <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(de_table))[1]
  if (is.na(padj_col)) {
    log_message("  Warning: No adjusted p-value column found. Trying unadjusted p-values...")
    pval_col <- intersect(c("pvalue", "P.Value", "p.value"), colnames(de_table))[1]
    if (!is.na(pval_col)) {
      log_message("  Applying BH correction to unadjusted p-values")
      de_table$padj_computed <- p.adjust(de_table[[pval_col]], method = "BH")
      padj_col <- "padj_computed"
    } else {
      log_message("  ERROR: No p-value columns found in DE table. Skipping enrichment.")
      return(NULL)
    }
  }

  # Use configurable FDR threshold (defined at line 126)
  fdr_thresh <- enrich_config$fdr_threshold %||% 0.05
  sig_genes <- de_table[[gene_col]][de_table[[padj_col]] < fdr_thresh]
  sig_genes <- unique(sig_genes[!is.na(sig_genes)])
  all_genes <- unique(de_table[[gene_col]][!is.na(de_table[[gene_col]])])

  log_message("  Significant genes: ", length(sig_genes), " / ", length(all_genes))

  if (length(sig_genes) < 5) {
    log_message("  Too few significant genes for enrichment")
    return(NULL)
  }

  # Check if collections are defined in config
  collections <- enrich_config$collections
  if (is.null(collections) || length(collections) == 0) {
    # Fall back to legacy behavior (GO + KEGG)
    return(run_gene_enrichment_legacy(de_table, gene_col, gmt, omic_name, config))
  }

  # Run enrichment for each collection
  all_results <- list()

  for (coll in collections) {
    coll_name <- coll$name
    log_message("  Running enrichment for collection: ", coll_name)

    coll_result <- run_enrichment_for_collection(
      sig_genes, all_genes, coll, omic_name, config,
      pval_thresh, min_size, max_size, top_n_terms
    )

    if (!is.null(coll_result)) {
      all_results[[coll_name]] <- coll_result
    }
  }

  return(all_results)
}

#' Run enrichment for a specific collection
run_enrichment_for_collection <- function(sig_genes, all_genes, collection, omic_name,
                                          config, pval_thresh, min_size, max_size, top_n_terms) {
  coll_name <- collection$name
  coll_type <- collection$type

  enrich_result <- NULL

  # Determine organism and OrgDb
  organism_name <- config$global$organism %||% "human"
  org_db <- NULL
  kegg_code <- NULL

  if (organism_name == "c_elegans") {
    if (requireNamespace("org.Ce.eg.db", quietly = TRUE)) {
      org_db <- org.Ce.eg.db::org.Ce.eg.db
      kegg_code <- "cel"
      species_code <- "cel"
    }
  } else {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      org_db <- org.Hs.eg.db::org.Hs.eg.db
      kegg_code <- "hsa"
      species_code <- "hsa"
    }
  }

  # Determine keyType
  key_type <- "SYMBOL"
  if (organism_name == "c_elegans") {
    if (all(grepl("^WBGene", sig_genes))) {
      key_type <- "WORMBASE"
    }
  } else if (all(grepl("^ENS", sig_genes))) {
    key_type <- "ENSEMBL"
  }

  # Map to ENTREZID
  gene_map <- NULL
  universe_map <- NULL

  if (!is.null(org_db) && requireNamespace("clusterProfiler", quietly = TRUE)) {
    tryCatch({
      gene_map <- clusterProfiler::bitr(sig_genes, fromType = key_type, toType = "ENTREZID", OrgDb = org_db)
      universe_map <- clusterProfiler::bitr(all_genes, fromType = key_type, toType = "ENTREZID", OrgDb = org_db)
    }, error = function(e) {
      log_message("    ID mapping failed: ", e$message)
    })
  }

  if (is.null(gene_map) || nrow(gene_map) == 0) {
    log_message("    No genes mapped for ", coll_name)
    return(NULL)
  }

  # Run enrichment based on collection type
  if (coll_type == "GO") {
    ont <- collection$ont %||% "BP"
    tryCatch({
      ego <- clusterProfiler::enrichGO(
        gene = gene_map$ENTREZID,
        universe = universe_map$ENTREZID,
        OrgDb = org_db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = pval_thresh,
        minGSSize = min_size,
        maxGSSize = max_size
      )

      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        enrich_result <- as.data.frame(ego)
        enrich_result$type <- "GO"
        enrich_result$collection <- coll_name
      }
    }, error = function(e) {
      log_message("    GO enrichment failed: ", e$message)
    })

  } else if (coll_type == "KEGG") {
    tryCatch({
      ekegg <- clusterProfiler::enrichKEGG(
        gene = gene_map$ENTREZID,
        universe = universe_map$ENTREZID,
        organism = kegg_code,
        pvalueCutoff = pval_thresh,
        minGSSize = min_size,
        maxGSSize = max_size
      )

      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        enrich_result <- as.data.frame(ekegg)
        enrich_result$type <- "KEGG"
        enrich_result$collection <- coll_name
      }
    }, error = function(e) {
      log_message("    KEGG enrichment failed: ", e$message)
    })

  } else if (coll_type == "msigdbr") {
    # MSigDB collections (H, C2, etc.)
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      log_message("    msigdbr package not available")
      return(NULL)
    }

    category <- collection$category
    subcategory <- collection$subcategory %||% NULL

    # Determine species for msigdbr
    msigdbr_species <- if (organism_name == "c_elegans") {
      "Caenorhabditis elegans"
    } else {
      "Homo sapiens"
    }

    tryCatch({
      # Get gene sets from msigdbr
      if (is.null(subcategory)) {
        m_df <- msigdbr::msigdbr(species = msigdbr_species, category = category)
      } else {
        m_df <- msigdbr::msigdbr(species = msigdbr_species, category = category, subcategory = subcategory)
      }

      # Convert to term2gene format
      term2gene <- m_df[, c("gs_name", "entrez_gene")]

      # Run enrichment
      emsig <- clusterProfiler::enricher(
        gene = gene_map$ENTREZID,
        universe = universe_map$ENTREZID,
        TERM2GENE = term2gene,
        pvalueCutoff = pval_thresh,
        minGSSize = min_size,
        maxGSSize = max_size
      )

      if (!is.null(emsig) && nrow(as.data.frame(emsig)) > 0) {
        enrich_result <- as.data.frame(emsig)
        enrich_result$type <- "MSigDB"
        enrich_result$collection <- coll_name
      }
    }, error = function(e) {
      log_message("    MSigDB enrichment failed: ", e$message)
    })

  } else if (coll_type == "gmt") {
    # Custom GMT file
    gmt_path <- collection$path
    if (!is.null(gmt_path) && file.exists(gmt_path)) {
      gmt_sets <- read_gmt_file(gmt_path)
      ora_result <- run_simple_ora(sig_genes, all_genes, gmt_sets, pval_thresh)
      if (!is.null(ora_result) && nrow(ora_result) > 0) {
        enrich_result <- ora_result
        enrich_result$type <- "GMT"
        enrich_result$collection <- coll_name
      }
    } else {
      log_message("    GMT file not found: ", gmt_path)
    }
  }

  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    log_message("    No enrichment results for ", coll_name)
    return(NULL)
  }

  # Save results table
  save_table(enrich_result, paste0(omic_name, "_", coll_name, "_enrichment.csv"), config)

  # Generate plots for multiple top_n values
  for (n in top_n_terms) {
    plot_df <- head(enrich_result, n)
    if (nrow(plot_df) > 0) {
      p <- plot_enrichment_dotplot(plot_df, paste0(omic_name, " - ", coll_name, " (Top ", n, ")"), n_terms = n)
      save_plot(p, paste0(omic_name, "_", coll_name, "_dotplot_top", n, ".png"), config, width = 10, height = 8)
    }
  }

  return(list(
    method = coll_type,
    collection = coll_name,
    results = enrich_result,
    sig_genes = sig_genes
  ))
}

#' Legacy enrichment function (fallback when no collections defined)
run_gene_enrichment_legacy <- function(de_table, gene_col, gmt, omic_name, config) {
  enrich_config <- config$enrichment
  pval_thresh <- enrich_config$ora_pvalue %||% 0.05
  min_size <- enrich_config$min_set_size %||% 10
  max_size <- enrich_config$max_set_size %||% 500

  # Use configurable FDR threshold
  fdr_thresh <- enrich_config$fdr_threshold %||% 0.05
  padj_col <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(de_table))[1]

  sig_genes <- de_table[[gene_col]][de_table[[padj_col]] < fdr_thresh]
  sig_genes <- unique(sig_genes[!is.na(sig_genes)])
  all_genes <- unique(de_table[[gene_col]][!is.na(de_table[[gene_col]])])

  # Load gene sets from GMT
  gene_sets <- NULL

  # Try custom GMT first
  if (!is.null(gmt)) {
    gene_sets <- gmt
    log_message("  Using custom GMT: ", length(gene_sets), " gene sets")
  }
  # Try clusterProfiler for GO/KEGG (legacy path)
  else if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    log_message("  Running enrichment via clusterProfiler...")

    enrich_results_list <- list()

    # Determine organism and OrgDb
    organism_name <- config$global$organism %||% "human"
    org_db <- NULL
    kegg_code <- NULL
    naming <- "human"

    if (organism_name == "c_elegans") {
      warning("Organism C. elegans detected: enabling specific database support.")
      if (requireNamespace("org.Ce.eg.db", quietly = TRUE)) {
        org_db <- org.Ce.eg.db::org.Ce.eg.db
        kegg_code <- "cel"
        naming <- "c_elegans"
      }
    } else {
      # Default to Human
      if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        org_db <- org.Hs.eg.db::org.Hs.eg.db
        kegg_code <- "hsa"
      }
    }

    # Determine keyType
    key_type <- "SYMBOL" # Default
    if (naming == "c_elegans") {
      # Check if IDs look like WormBase IDs (WBGene...)
      if (all(grepl("^WBGene", sig_genes))) {
        key_type <- "WORMBASE"
        log_message("    Detected WormBase IDs, using keyType = 'WORMBASE'")
      } else if (gene_col == "Wormbase_id") {
        # Explicit column name
        key_type <- "WORMBASE"
        log_message("    Column is Wormbase_id, using keyType = 'WORMBASE'")
      } else if (all(grepl("^ENS", sig_genes))) {
        key_type <- "ENSEMBL"
      }
    } else {
      if (all(grepl("^ENS", sig_genes))) {
        key_type <- "ENSEMBL"
      }
    }

    if (is.null(org_db)) {
      log_message("    OrgDb for ", organism_name, " not available. Skipping clusterProfiler.")
    } else {
      # explicit mapping to ENTREZID to robustify
      tryCatch(
        {
          gene_map <- clusterProfiler::bitr(sig_genes, fromType = key_type, toType = "ENTREZID", OrgDb = org_db)
          universe_map <- clusterProfiler::bitr(all_genes, fromType = key_type, toType = "ENTREZID", OrgDb = org_db)
        },
        error = function(e) {
          log_message("    ID mapping failed: ", e$message)
          gene_map <- NULL
          universe_map <- NULL
        }
      )

      if (!is.null(gene_map) && nrow(gene_map) > 0) {
        log_message("    Mapped ", nrow(gene_map), " genes to ENTREZID (from ", key_type, ")")

        # 1. GO Enrichment (BP)
        tryCatch(
          {
            log_message("    Running GO Enrichment (BP) using ENTREZID...")
            ego <- clusterProfiler::enrichGO(
              gene = gene_map$ENTREZID,
              universe = universe_map$ENTREZID,
              OrgDb = org_db,
              keyType = "ENTREZID",
              ont = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff = pval_thresh,
              minGSSize = min_size,
              maxGSSize = max_size
            )

            if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
              res_go <- as.data.frame(ego)
              res_go$type <- "GO"
              enrich_results_list[["GO"]] <- res_go
            }
          },
          error = function(e) {
            log_message("    GO enrichment failed: ", e$message)
          }
        )

        # 2. KEGG Enrichment
        # Check config for KEGG (default to TRUE if not specified, to enable pathview)
        if (enrich_config$use_kegg %||% TRUE) {
          log_message("    Running KEGG Enrichment...")
          tryCatch(
            {
              if (!is.null(gene_map) && nrow(gene_map) > 0) {
                ekegg <- clusterProfiler::enrichKEGG(
                  gene = gene_map$ENTREZID,
                  universe = universe_map$ENTREZID,
                  organism = kegg_code,
                  pvalueCutoff = pval_thresh,
                  minGSSize = min_size,
                  maxGSSize = max_size
                )

                if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
                  res_kegg <- as.data.frame(ekegg)
                  res_kegg$type <- "KEGG"
                  enrich_results_list[["KEGG"]] <- res_kegg
                  log_message("    KEGG enrichment found ", nrow(res_kegg), " pathways")
                }
              }
            },
            error = function(e) {
              log_message("    KEGG enrichment failed: ", e$message)
            }
          )
        }
      } else {
        log_message("    Could not map significant genes to ENTREZID. Skipping enrichment.")
      }
    }

    # Combine results
    if (length(enrich_results_list) > 0) {
      result <- do.call(rbind, enrich_results_list)
      # Fill missing columns if any differences

      result$omics <- omic_name
      save_table(result, paste0(omic_name, "_enrichment_combined.csv"), config)

      # Plot top terms (prioritize by pvalue)
      if (nrow(result) > 0) {
        p <- plot_enrichment_dotplot(result, paste0(omic_name, " Enrichment"))
        save_plot(p, paste0(omic_name, "_enrichment_dotplot.png"), config, width = 10, height = 8)
      }

      return(list(
        method = "clusterProfiler_Combined",
        results = result,
        sig_genes = sig_genes
      ))
    }
  }

  # Fallback: simple ORA with custom gene sets
  if (!is.null(gene_sets)) {
    ora_result <- run_simple_ora(sig_genes, all_genes, gene_sets, pval_thresh)
    if (!is.null(ora_result) && nrow(ora_result) > 0) {
      ora_result$omics <- omic_name
      save_table(ora_result, paste0(omic_name, "_ORA_enrichment.csv"), config)

      p <- plot_enrichment_dotplot(ora_result, paste0(omic_name, " ORA Enrichment"))
      save_plot(p, paste0(omic_name, "_ORA_enrichment.png"), config, width = 10, height = 8)

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

    if (n_gs < 5) {
      return(NULL)
    }

    overlap <- intersect(sig_genes, gs_in_bg)
    n_overlap <- length(overlap)

    if (n_overlap == 0) {
      return(NULL)
    }

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

  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }

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
  if (is.na(padj_col)) {
    return(NULL)
  }

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
      save_plot(p, "metabolomics_ORA_enrichment.png", config, width = 10, height = 8)

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
#'
#' This function integrates enrichment results from multiple omics layers by performing
#' a meta-analysis of p-values for common pathways/terms. It identifies pathways that are
#' consistently significant across different biological layers (e.g., RNA and Protein),
#' providing a more robust view of system-level perturbations.
#'
#' The method works as follows:
#' 1. Standardizes identifying terms (e.g., using IDs or Descriptions) across all input results.
#' 2. Identifies the intersection of terms involved in at least 2 omics layers.
#' 3. Combines the p-values for these common terms using established meta-analysis methods:
#'    - Fisher's method: Sum of log p-values (Chi-squared distribution). Good for detecting if *any* omics is significant.
#'    - Stouffer's method: Z-score combination. Good for consensus significance.
#'    - Min-P: Takes the minimum p-value (conservative).
#' 4. Adjusts the combined p-values for multiple hypothesis testing (Benjamini-Hochberg).
#' 5. Visualizes the top combined significantly enriched terms.
combine_enrichment_results <- function(per_omics_results, method = "fisher", config) {
  log_message("Combining enrichment results across omics (method: ", method, ")...")

  # Collect all enrichment results
  all_results <- list()
  for (omic in names(per_omics_results)) {
    res <- per_omics_results[[omic]]
    if (!is.null(res) && !is.null(res$results)) {
      df <- res$results
      # Standardize columns for merging
      # We want to merge on IDs (term_id) but keep Descriptions (term_label)
      if ("ID" %in% colnames(df)) {
        df$term_id <- df$ID
      } else if ("term" %in% colnames(df)) {
        # Fallback if ID not present but term is
        df$term_id <- df$term
      } else {
        df$term_id <- rownames(df)
      }

      if ("Description" %in% colnames(df)) {
        df$term_label <- df$Description
      } else {
        df$term_label <- df$term_id
      }

      # Ensure term column is term_id for backward compatibility if needed,
      # but we will use term_id for intersection
      df$term <- df$term_id

      all_results[[omic]] <- df
    }
  }

  if (length(all_results) < 2) {
    log_message("Need at least 2 omics with enrichment results")
    return(NULL)
  }

  # Find common terms
  # Find common terms based on term_id
  all_terms <- Reduce(intersect, lapply(all_results, function(x) x$term_id))

  if (length(all_terms) < 1) {
    log_message("No common enriched terms across omics")
    return(NULL)
  }

  log_message("Found ", length(all_terms), " terms enriched in multiple omics")

  # Combine p-values
  # Combine p-values
  combined_df <- data.frame(term = all_terms, term_id = all_terms, stringsAsFactors = FALSE)

  # Try to assign term_label from the first omics that has it
  # We iterate to find the first non-NA label for each term
  labels <- rep(NA, length(all_terms))
  for (omic in names(all_results)) {
    res <- all_results[[omic]]
    idx <- match(all_terms, res$term_id)
    # Update NA labels
    na_mask <- is.na(labels)
    if (any(na_mask)) {
      labels[na_mask] <- res$term_label[idx][na_mask]
    }
  }
  combined_df$term_label <- labels
  # Fill remaining NAs with ID
  combined_df$term_label[is.na(combined_df$term_label)] <- combined_df$term_id[is.na(combined_df$term_label)]

  for (omic in names(all_results)) {
    res <- all_results[[omic]]
    idx <- match(all_terms, res$term_id)
    combined_df[[paste0("pvalue_", omic)]] <- res$pvalue[idx]
    combined_df[[paste0("padj_", omic)]] <- res$padj[idx]
  }

  # Combine p-values
  pval_cols <- grep("^pvalue_", colnames(combined_df), value = TRUE)
  pvals_matrix <- as.matrix(combined_df[, pval_cols])

  combined_df$combined_pvalue <- apply(pvals_matrix, 1, function(pvals) {
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) < 2) {
      return(NA)
    }

    if (method == "fisher") {
      fisher_combine_pvalues(pvals)
    } else if (method == "stouffer") {
      stouffer_combine_pvalues(pvals)
    } else {
      min(pvals) # Minimum p-value
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

    save_plot(p, "combined_enrichment.png", config, width = 10, height = 8)
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
          tryCatch(
            {
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
            },
            error = function(e) NULL
          )
        }
      }
    }
  }

  if (length(results) > 0) {
    log_message("MOFA factor enrichment complete: ", length(results), " enrichment sets")
  }

  results
}

#' Plot RNA vs Protein enrichment scatter plot
#'
#' Visualizes the concordance of pathway enrichment results between Transcriptomics (RNA) and
#' Proteomics layers. This helps identify pathways that are consistently regulated at both
#' expression levels versus those with discordant regulation (e.g., buffering or post-transcriptional control).
#'
#' Methodology:
#' 1. Standardizes results from both omics layers (handling different column names from ORA/clusterProfiler).
#' 2. Robust Merging: Uses stable identifiers (`term_id`, e.g., GO IDs) to match pathways between layers,
#'    avoiding mismatches due to slight differences in description text.
#' 3. Labeling: Preserves human-readable descriptions (`term_label`) for the final plot labels.
#' 4. Imputation: Missing values (non-significant in one layer) are imputed to FDR=1 for plotting purposes.
#' 5. Visualization:
#'    - X-axis: -log10(FDR) for RNA
#'    - Y-axis: -log10(FDR) for Protein
#'    - Point Size: Average number of genes in the pathway
#'    - Color: Average fold enrichment
#'    - Labels: Top 20 most significant terms
#'
#' @param rna_results List containing RNA enrichment results
#' @param prot_results List containing Protein enrichment results
#' @param config Pipeline configuration object
plot_rna_protein_enrichment_scatter <- function(rna_results, prot_results, config) {
  log_message("Creating RNA vs Protein enrichment scatter plot...")

  if (is.null(rna_results) || is.null(prot_results)) {
    log_message("  Missing enrichment results for comparison")
    return(NULL)
  }

  rna_df <- rna_results$results
  prot_df <- prot_results$results

  if (is.null(rna_df) || is.null(prot_df) ||
    nrow(rna_df) == 0 || nrow(prot_df) == 0) {
    log_message("  No enrichment results to compare")
    return(NULL)
  }

  tryCatch(
    {
      # Standardize column names
      # clusterProfiler uses: ID, Description, pvalue, p.adjust, Count, GeneRatio
      # Our ORA uses: term, pvalue, padj, overlap

      # Standardize RNA results
      # Standardize RNA results
      # Create term_id for merging and term_label for plotting
      if ("ID" %in% colnames(rna_df)) {
        rna_df$term_id <- rna_df$ID
      } else {
        rna_df$term_id <- rownames(rna_df)
      }

      if ("Description" %in% colnames(rna_df)) {
        rna_df$term_label <- rna_df$Description
      } else {
        rna_df$term_label <- rna_df$term_id
      }

      if ("p.adjust" %in% colnames(rna_df)) {
        rna_df$FDR <- rna_df$p.adjust
      } else if ("padj" %in% colnames(rna_df)) {
        rna_df$FDR <- rna_df$padj
      }

      if ("Count" %in% colnames(rna_df)) {
        rna_df$count <- rna_df$Count
      } else if ("overlap" %in% colnames(rna_df)) {
        rna_df$count <- rna_df$overlap
      }

      # Calculate fold enrichment if not present
      if (!"fold_enrichment" %in% colnames(rna_df)) {
        if ("GeneRatio" %in% colnames(rna_df) && "BgRatio" %in% colnames(rna_df)) {
          gene_ratio <- sapply(strsplit(as.character(rna_df$GeneRatio), "/"), function(x) {
            as.numeric(x[1]) / as.numeric(x[2])
          })
          bg_ratio <- sapply(strsplit(as.character(rna_df$BgRatio), "/"), function(x) {
            as.numeric(x[1]) / as.numeric(x[2])
          })
          rna_df$fold_enrichment <- gene_ratio / bg_ratio
        } else {
          rna_df$fold_enrichment <- 1
        }
      }

      # Standardize Protein results
      # Standardize Protein results
      if ("ID" %in% colnames(prot_df)) {
        prot_df$term_id <- prot_df$ID
      } else {
        prot_df$term_id <- rownames(prot_df)
      }

      if ("Description" %in% colnames(prot_df)) {
        prot_df$term_label <- prot_df$Description
      } else {
        prot_df$term_label <- prot_df$term_id
      }

      if ("p.adjust" %in% colnames(prot_df)) {
        prot_df$FDR <- prot_df$p.adjust
      } else if ("padj" %in% colnames(prot_df)) {
        prot_df$FDR <- prot_df$padj
      }

      if ("Count" %in% colnames(prot_df)) {
        prot_df$count <- prot_df$Count
      } else if ("overlap" %in% colnames(prot_df)) {
        prot_df$count <- prot_df$overlap
      }

      # Calculate fold enrichment if not present
      if (!"fold_enrichment" %in% colnames(prot_df)) {
        if ("GeneRatio" %in% colnames(prot_df) && "BgRatio" %in% colnames(prot_df)) {
          gene_ratio <- sapply(strsplit(as.character(prot_df$GeneRatio), "/"), function(x) {
            as.numeric(x[1]) / as.numeric(x[2])
          })
          bg_ratio <- sapply(strsplit(as.character(prot_df$BgRatio), "/"), function(x) {
            as.numeric(x[1]) / as.numeric(x[2])
          })
          prot_df$fold_enrichment <- gene_ratio / bg_ratio
        } else {
          prot_df$fold_enrichment <- 1
        }
      }

      # Merge on term_id
      rna_subset <- rna_df[, c("term_id", "term_label", "FDR", "count", "fold_enrichment")]
      colnames(rna_subset) <- c("term_id", "term_label_rna", "FDR_RNA", "Count_RNA", "Fold_RNA")

      prot_subset <- prot_df[, c("term_id", "term_label", "FDR", "count", "fold_enrichment")]
      colnames(prot_subset) <- c("term_id", "term_label_prot", "FDR_Prot", "Count_Prot", "Fold_Prot")

      merged_df <- merge(rna_subset, prot_subset, by = "term_id", all = TRUE)

      if (nrow(merged_df) == 0) {
        log_message("  No terms found for comparison")
        return(NULL)
      }

      # Impute missing values
      merged_df$FDR_RNA[is.na(merged_df$FDR_RNA)] <- 1
      merged_df$FDR_Prot[is.na(merged_df$FDR_Prot)] <- 1
      merged_df$Count_RNA[is.na(merged_df$Count_RNA)] <- 0
      merged_df$Count_Prot[is.na(merged_df$Count_Prot)] <- 0
      merged_df$Fold_RNA[is.na(merged_df$Fold_RNA)] <- 0
      merged_df$Fold_Prot[is.na(merged_df$Fold_Prot)] <- 0

      # Calculate -log10(FDR)
      merged_df$logFDR_RNA <- -log10(merged_df$FDR_RNA)
      merged_df$logFDR_Prot <- -log10(merged_df$FDR_Prot)

      # Resolve TermName from labels
      # Prefer RNA label, then Protein label, then ID
      merged_df$TermName <- ifelse(!is.na(merged_df$term_label_rna), merged_df$term_label_rna,
        ifelse(!is.na(merged_df$term_label_prot), merged_df$term_label_prot, merged_df$term_id)
      )

      # Clean term names (truncate only)
      clean_term <- function(term_str) {
        # Truncate long names
        if (nchar(term_str) > 50) {
          term_str <- paste0(substr(term_str, 1, 47), "...")
        }
        return(term_str)
      }

      merged_df$TermName <- sapply(merged_df$TermName, clean_term)

      # Filter to only include terms significant in at least one omics
      sig_threshold <- -log10(0.05)
      merged_df_sig <- merged_df[merged_df$logFDR_RNA >= sig_threshold |
        merged_df$logFDR_Prot >= sig_threshold, ]

      if (nrow(merged_df_sig) == 0) {
        log_message("  No significant terms for scatter plot")
        return(NULL)
      }

      log_message(paste("  Total unique terms:", nrow(merged_df_sig)))

      # Calculate average metrics for plotting
      merged_df_sig$avg_count <- (merged_df_sig$Count_RNA + merged_df_sig$Count_Prot) / 2
      merged_df_sig$avg_fold <- (merged_df_sig$Fold_RNA + merged_df_sig$Fold_Prot) / 2

      # Create scatter plot
      p <- ggplot2::ggplot(merged_df_sig, ggplot2::aes(
        x = logFDR_RNA,
        y = logFDR_Prot,
        label = TermName
      )) +
        ggplot2::geom_point(ggplot2::aes(
          size = avg_count,
          color = avg_fold
        ), alpha = 0.7) +
        ggplot2::geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dashed",
          color = "grey50"
        ) +
        ggplot2::geom_vline(
          xintercept = sig_threshold,
          linetype = "dotted",
          color = "grey70"
        ) +
        ggplot2::geom_hline(
          yintercept = sig_threshold,
          linetype = "dotted",
          color = "grey70"
        )

      # Add text labels for top terms
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        # Label top 20 terms by combined significance
        merged_df_sig$combined_sig <- merged_df_sig$logFDR_RNA + merged_df_sig$logFDR_Prot
        top_terms <- head(merged_df_sig[order(-merged_df_sig$combined_sig), ], 20)

        p <- p + ggrepel::geom_text_repel(
          data = top_terms,
          size = 3,
          max.overlaps = 20,
          box.padding = 0.5,
          force = 2
        )
      }

      p <- p +
        ggplot2::scale_color_viridis_c(name = "Avg Fold\nEnrichment") +
        ggplot2::scale_size_continuous(name = "Avg Gene\nCount") +
        ggplot2::labs(
          title = "RNA-seq vs Proteomics Enrichment Comparison",
          subtitle = paste("Total Unique Terms:", nrow(merged_df_sig)),
          x = "RNA-seq [-log10(FDR)]",
          y = "Proteomics [-log10(FDR)]"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.title = ggplot2::element_text(face = "bold"),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5)
        )

      # Save plot
      save_plot(p, "rna_protein_enrichment_scatter.png", config,
        width = 10, height = 10
      )

      # Save merged data
      save_table(
        merged_df_sig,
        "rna_protein_enrichment_scatter_data.csv",
        config
      )

      log_message("  Scatter plot saved successfully")

      return(p)
    },
    error = function(e) {
      log_message("  Error creating scatter plot: ", e$message)
      return(NULL)
    }
  )
}

#' Read GMT file
read_gmt_file <- function(gmt_path) {
  if (!file.exists(gmt_path)) {
    return(NULL)
  }

  lines <- readLines(gmt_path)
  gene_sets <- list()

  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 3) next

    set_name <- parts[1]
    # parts[2] is typically description, skip it
    genes <- parts[3:length(parts)]
    genes <- genes[nchar(genes) > 0]  # Remove empty strings

    gene_sets[[set_name]] <- genes
  }

  return(gene_sets)
}

#' Plot enrichment dotplot
plot_enrichment_dotplot <- function(enrich_df, title, n_terms = 20) {
  tryCatch(
    {
      if (nrow(enrich_df) == 0) {
        return(ggplot2::ggplot() +
          ggplot2::ggtitle("No enriched terms"))
      }

      # Select top terms
      plot_df <- head(enrich_df, n_terms)

      # Determine which columns exist
      size_col <- if ("overlap" %in% colnames(plot_df)) "overlap" else "Count"
      pval_col <- if ("padj" %in% colnames(plot_df)) "padj" else "p.adjust"

      # Create required columns safely
      if (size_col %in% colnames(plot_df)) {
        plot_df$overlap <- plot_df[[size_col]]
      } else {
        plot_df$overlap <- 1
      }

      if (pval_col %in% colnames(plot_df)) {
        plot_df$padj <- plot_df[[pval_col]]
      } else {
        # Identify SOME pvalue column
        if ("pvalue" %in% colnames(plot_df)) {
          plot_df$padj <- plot_df$pvalue
        } else {
          # Last resort
          plot_df$padj <- 0.05
        }
      }

      # Handle 'term' column (clusterProfiler uses 'Description', generic ORA uses 'term')
      if (!"term" %in% colnames(plot_df)) {
        if ("Description" %in% colnames(plot_df)) {
          plot_df$term <- plot_df$Description
        } else if ("ID" %in% colnames(plot_df)) {
          plot_df$term <- plot_df$ID
        } else {
          plot_df$term <- rownames(plot_df)
        }
      }

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
    },
    error = function(e) {
      log_message("Error in plot_enrichment_dotplot: ", e$message)
      return(ggplot2::ggplot() +
        ggplot2::ggtitle("Plot Generation Failed"))
    }
  )
}

#' Run Pathview for consensus pathways
#'
#' Visualizes the top consensus pathways (enriched in multiple omics) using Pathview.
#' It maps RNA and Protein data to Entrez IDs and Metabolites to KEGG Compound IDs,
#' then overlays their fold-changes on the pathway diagrams.
#'
#' @param enrichment_results List containing the combined enrichment results
#' @param mae_data The MultiAssayExperiment object containing omics data
#' @param config Pipeline configuration object
#' @param top_n Number of top consensus pathways to visualize
run_consensus_pathview <- function(enrichment_results, mae_data, config, top_n = 5) {
  log_message("=== Running Consensus Pathview Visualization ===")

  if (is.null(enrichment_results) || is.null(enrichment_results$combined)) {
    log_message("No combined enrichment results available.")
    return(NULL)
  }

  if (!requireNamespace("pathview", quietly = TRUE)) {
    log_message("Package 'pathview' not installed. Skipping.")
    return(NULL)
  }

  # 1. Identify Top Consensus KEGG Pathways
  combined <- enrichment_results$combined

  # Filter for KEGG-like terms (e.g., hsa01100, cel01100, or just numeric 01100)
  # Look for terms that match KEGG ID pattern: 3-4 letters + 5 digits
  kegg_pattern <- "^[a-z]{3,4}\\d{5}$"
  kegg_indices <- grep(kegg_pattern, combined$term)

  if (length(kegg_indices) == 0) {
    log_message("No KEGG pathway IDs found in combined enrichment terms. Terms might be Descriptions.")
    return(NULL)
  }

  top_pathways <- head(combined$term[kegg_indices], top_n)

  if (length(top_pathways) == 0) {
    log_message("No top pathways found.")
    return(NULL)
  }

  log_message("Top consensus pathways: ", paste(top_pathways, collapse = ", "))

  # 2. Prepare Data (RNA & Proteomics)
  # We construct a matrix where Rows = EntrezIDs, Cols = Layers (RNA, Protein)

  # Helper to extract IDs and LogFC from DE table
  get_omics_data <- function(omic_name) {
    if (!omic_name %in% names(mae_data$harmonized_omics)) {
      return(NULL)
    }

    de_table <- mae_data$harmonized_omics[[omic_name]]$de_table
    if (is.null(de_table)) {
      return(NULL)
    }

    # Identify proper columns
    id_col <- if ("gene_symbol" %in% colnames(de_table)) "gene_symbol" else "gene_id"
    fc_col <- grep("logFC|log2FC", colnames(de_table), value = TRUE)[1]

    if (is.na(fc_col)) {
      return(NULL)
    }

    return(data.frame(ID = de_table[[id_col]], FC = de_table[[fc_col]], stringsAsFactors = FALSE))
  }

  rna_df <- get_omics_data("transcriptomics")
  prot_df <- get_omics_data("proteomics")

  # ID Mapping Setup
  map_to_entrez <- function(ids, config) {
    if (requireNamespace("clusterProfiler", quietly = TRUE)) {
      organism <- config$global$organism %||% "human"
      org_db <- if (organism == "c_elegans") "org.Ce.eg.db" else "org.Hs.eg.db"

      if (requireNamespace(org_db, quietly = TRUE)) {
        db_obj <- get(org_db)
        tryCatch(
          {
            map <- clusterProfiler::bitr(ids, "SYMBOL", "ENTREZID", db_obj)
            return(map)
          },
          error = function(e) {
            return(NULL)
          }
        )
      }
    }
    return(NULL)
  }

  # Process RNA
  rna_data <- NULL
  if (!is.null(rna_df)) {
    map <- map_to_entrez(rna_df$ID, config)
    if (!is.null(map)) {
      merged <- merge(rna_df, map, by.x = "ID", by.y = "SYMBOL")
      # Handle duplicates by mean
      agg <- aggregate(FC ~ ENTREZID, merged, mean)
      rna_data <- setNames(agg$FC, agg$ENTREZID)
    }
  }

  # Process Protein
  prot_data <- NULL
  if (!is.null(prot_df)) {
    map <- map_to_entrez(prot_df$ID, config)
    if (!is.null(map)) {
      merged <- merge(prot_df, map, by.x = "ID", by.y = "SYMBOL")
      agg <- aggregate(FC ~ ENTREZID, merged, mean)
      prot_data <- setNames(agg$FC, agg$ENTREZID)
    }
  }

  # Combine RNA + Protein
  all_entrez <- unique(c(names(rna_data), names(prot_data)))
  gene_data <- NULL

  if (length(all_entrez) > 0) {
    gene_data <- matrix(NA,
      nrow = length(all_entrez), ncol = 2,
      dimnames = list(all_entrez, c("RNA", "Protein"))
    )

    if (!is.null(rna_data)) gene_data[names(rna_data), "RNA"] <- rna_data
    if (!is.null(prot_data)) gene_data[names(prot_data), "Protein"] <- prot_data
  }

  # 3. Prepare Metabolites
  cpd_data <- NULL
  if ("metabolomics" %in% names(mae_data$harmonized_omics)) {
    meta_df <- get_omics_data("metabolomics") # Uses 'gene_id' fallback which might be wrong for metabolites key
    # Re-fetch specific for metabolomics
    de_table <- mae_data$harmonized_omics$metabolomics$de_table
    if (!is.null(de_table)) {
      id_col <- grep("kegg|KEGG", colnames(de_table), value = TRUE)[1]
      fc_col <- grep("logFC|log2FC", colnames(de_table), value = TRUE)[1]

      if (!is.na(id_col) && !is.na(fc_col)) {
        cpd_data <- setNames(de_table[[fc_col]], de_table[[id_col]])
      }
    }
  }

  # 4. Run Pathview
  out_dir <- file.path(config$output$output_dir, "plots", "pathview_consensus")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd))

  species <- gsub("\\d+", "", top_pathways[1])

  results_list <- list()

  for (pid in top_pathways) {
    tryCatch(
      {
        log_message("  Generating plot for ", pid)
        pathview::pathview(
          gene.data = gene_data,
          cpd.data = cpd_data,
          pathway.id = pid,
          species = species,
          gene.idtype = "entrez",
          kegg.native = TRUE,
          multi.state = TRUE, # Important for multi-column gene data
          same.layer = FALSE,
          low = list(gene = "blue", cpd = "blue"),
          mid = list(gene = "gray", cpd = "gray"),
          high = list(gene = "red", cpd = "red"),
          out.suffix = "consensus"
        )
        results_list[[pid]] <- file.path(out_dir, paste0(pid, ".consensus.png"))
      },
      error = function(e) {
        log_message("  Failed for ", pid, ": ", e$message)
      }
    )
  }

  return(results_list)
}
