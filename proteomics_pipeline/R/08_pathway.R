# =============================================================================
# Pathway and Gene Set Enrichment Analysis
# =============================================================================

#' Run pathway enrichment analysis
#'
#' @param de_results Differential expression results
#' @param imputed_data Imputed data with annotations
#' @param config Configuration list
#' @return Pathway analysis results
run_pathway_analysis <- function(de_results, imputed_data, config) {
  log_message("=== Starting Pathway Analysis ===")

  if (!config$pathway$run_pathway_analysis %||% TRUE) {
    log_message("Pathway analysis disabled in config. Skipping.")
    return(NULL)
  }

  annotations <- imputed_data$annotations
  organism <- config$input$organism

  # Check if we have gene symbols for enrichment
  has_gene_symbols <- "gene_symbol" %in% colnames(annotations) &&
                      sum(!is.na(annotations$gene_symbol)) > nrow(annotations) * 0.5

  # Check for custom GMT file
  has_gmt <- !is.null(config$optional_inputs$gene_set_gmt) &&
             file.exists(config$optional_inputs$gene_set_gmt)

  if (!has_gene_symbols && !has_gmt) {
    log_message("Insufficient gene symbol mappings and no custom GMT file provided.")
    log_message("Skipping pathway analysis. Provide a mapping_file or gene_set_gmt in config.")
    return(NULL)
  }

  all_results <- list()

  for (contrast_name in names(de_results$results)) {
    log_message("Running pathway analysis for: ", contrast_name)

    de_table <- de_results$results[[contrast_name]]$table

    # Merge with annotations if not already present
    if (!"gene_symbol" %in% colnames(de_table)) {
      de_table <- merge(de_table, annotations[, c("feature_id", "gene_symbol", "entrez_id")],
                        by = "feature_id", all.x = TRUE)
    }

    contrast_results <- list()

    # ORA (Over-Representation Analysis)
    ora_results <- run_ora(de_table, organism, config)
    if (!is.null(ora_results)) {
      contrast_results$ora <- ora_results
    }

    # GSEA (Gene Set Enrichment Analysis)
    gsea_results <- run_gsea(de_table, organism, config)
    if (!is.null(gsea_results)) {
      contrast_results$gsea <- gsea_results
    }

    # Custom GMT analysis if provided
    if (has_gmt) {
      gmt_results <- run_custom_gmt_analysis(de_table, config)
      if (!is.null(gmt_results)) {
        contrast_results$custom_gmt <- gmt_results
      }
    }

    all_results[[contrast_name]] <- contrast_results

    # Save results
    clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))
    if (!is.null(ora_results)) {
      save_table(ora_results, paste0("pathway_ora_", clean_name, ".csv"), config, "tables")
    }
    if (!is.null(gsea_results)) {
      save_table(gsea_results, paste0("pathway_gsea_", clean_name, ".csv"), config, "tables")
    }
  }

  # Create pathway plots
  plots <- create_pathway_plots(all_results, config)

  log_message("=== Pathway Analysis Complete ===")

  list(
    results = all_results,
    plots = plots
  )
}

#' Run Over-Representation Analysis
#'
#' @param de_table Differential expression table
#' @param organism Organism name
#' @param config Configuration list
#' @return ORA results data frame
run_ora <- function(de_table, organism, config) {
  # Get significant genes
  sig_genes <- de_table$gene_symbol[de_table$significant & !is.na(de_table$gene_symbol)]
  background_genes <- de_table$gene_symbol[!is.na(de_table$gene_symbol)]

  if (length(sig_genes) < 5) {
    log_message("Too few significant genes for ORA (", length(sig_genes), ")")
    return(NULL)
  }

  log_message("Running ORA with ", length(sig_genes), " significant genes")

  # Try clusterProfiler if available and organism is supported
  if (requireNamespace("clusterProfiler", quietly = TRUE) &&
      is_supported_organism(organism)) {
    return(run_ora_clusterprofiler(sig_genes, background_genes, organism, config))
  }

  # Fallback to fgsea with custom GMT if available
  if (!is.null(config$optional_inputs$gene_set_gmt)) {
    return(run_ora_custom(sig_genes, background_genes, config))
  }

  log_message("ORA not available for this organism. Provide a custom GMT file.")
  return(NULL)
}

#' Check if organism is supported by annotation databases
#'
#' @param organism Organism name
#' @return Logical
is_supported_organism <- function(organism) {
  supported <- c("homo sapiens", "human", "mus musculus", "mouse",
                 "rattus norvegicus", "rat", "danio rerio", "zebrafish",
                 "drosophila melanogaster", "saccharomyces cerevisiae", "yeast")
  tolower(organism) %in% supported
}

#' Get OrgDb for organism
#'
#' @param organism Organism name
#' @return OrgDb package name or NULL
get_orgdb <- function(organism) {
  orgdb_map <- list(
    "homo sapiens" = "org.Hs.eg.db",
    "human" = "org.Hs.eg.db",
    "mus musculus" = "org.Mm.eg.db",
    "mouse" = "org.Mm.eg.db",
    "rattus norvegicus" = "org.Rn.eg.db",
    "rat" = "org.Rn.eg.db",
    "danio rerio" = "org.Dr.eg.db",
    "zebrafish" = "org.Dr.eg.db",
    "drosophila melanogaster" = "org.Dm.eg.db",
    "saccharomyces cerevisiae" = "org.Sc.sgd.db",
    "yeast" = "org.Sc.sgd.db"
  )

  orgdb_map[[tolower(organism)]]
}

#' Run ORA using clusterProfiler
#'
#' @param sig_genes Significant gene symbols
#' @param background_genes All gene symbols tested
#' @param organism Organism name
#' @param config Configuration list
#' @return ORA results
run_ora_clusterprofiler <- function(sig_genes, background_genes, organism, config) {
  orgdb_name <- get_orgdb(organism)

  if (is.null(orgdb_name)) {
    log_message("No OrgDb available for ", organism)
    return(NULL)
  }

  # Check if OrgDb is installed
  if (!requireNamespace(orgdb_name, quietly = TRUE)) {
    log_message("OrgDb package ", orgdb_name, " not installed. Skipping clusterProfiler ORA.")
    return(NULL)
  }

  tryCatch({
    orgdb <- get(orgdb_name, envir = asNamespace(orgdb_name))

    databases <- config$pathway$databases %||% c("GO", "KEGG")
    go_ontologies <- config$pathway$go_ontologies %||% c("BP", "MF")
    min_size <- config$pathway$min_set_size %||% 10
    max_size <- config$pathway$max_set_size %||% 500
    pval_thresh <- config$pathway$ora_pvalue_threshold %||% 0.05

    all_results <- list()

    # GO enrichment
    if ("GO" %in% databases) {
      for (ont in go_ontologies) {
        log_message("Running GO ", ont, " enrichment...")

        go_result <- clusterProfiler::enrichGO(
          gene = sig_genes,
          universe = background_genes,
          OrgDb = orgdb,
          keyType = "SYMBOL",
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = pval_thresh,
          qvalueCutoff = 0.2,
          minGSSize = min_size,
          maxGSSize = max_size
        )

        if (!is.null(go_result) && nrow(go_result@result) > 0) {
          result_df <- as.data.frame(go_result)
          result_df$database <- paste0("GO_", ont)
          all_results[[paste0("GO_", ont)]] <- result_df
        }
      }
    }

    # KEGG enrichment
    if ("KEGG" %in% databases) {
      log_message("Running KEGG enrichment...")

      # Get KEGG organism code
      kegg_code <- get_kegg_code(organism)

      if (!is.null(kegg_code)) {
        # Convert symbols to Entrez IDs
        symbol_to_entrez <- clusterProfiler::bitr(
          sig_genes,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = orgdb
        )

        background_entrez <- clusterProfiler::bitr(
          background_genes,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = orgdb
        )

        if (nrow(symbol_to_entrez) > 0) {
          kegg_result <- clusterProfiler::enrichKEGG(
            gene = symbol_to_entrez$ENTREZID,
            universe = background_entrez$ENTREZID,
            organism = kegg_code,
            pAdjustMethod = "BH",
            pvalueCutoff = pval_thresh,
            minGSSize = min_size,
            maxGSSize = max_size
          )

          if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
            result_df <- as.data.frame(kegg_result)
            result_df$database <- "KEGG"
            all_results[["KEGG"]] <- result_df
          }
        }
      }
    }

    # Combine results
    if (length(all_results) > 0) {
      combined <- do.call(rbind, all_results)
      combined <- combined[order(combined$p.adjust), ]
      return(combined)
    }

    return(NULL)
  }, error = function(e) {
    log_message("clusterProfiler ORA failed: ", e$message)
    return(NULL)
  })
}

#' Get KEGG organism code
#'
#' @param organism Organism name
#' @return KEGG code or NULL
get_kegg_code <- function(organism) {
  kegg_map <- list(
    "homo sapiens" = "hsa",
    "human" = "hsa",
    "mus musculus" = "mmu",
    "mouse" = "mmu",
    "rattus norvegicus" = "rno",
    "rat" = "rno",
    "danio rerio" = "dre",
    "zebrafish" = "dre",
    "drosophila melanogaster" = "dme",
    "saccharomyces cerevisiae" = "sce",
    "yeast" = "sce"
  )

  kegg_map[[tolower(organism)]]
}

#' Run ORA with custom gene sets
#'
#' @param sig_genes Significant genes
#' @param background_genes Background genes
#' @param config Configuration list
#' @return ORA results
run_ora_custom <- function(sig_genes, background_genes, config) {
  gmt_file <- config$optional_inputs$gene_set_gmt

  if (!file.exists(gmt_file)) {
    log_message("GMT file not found: ", gmt_file)
    return(NULL)
  }

  # Read GMT file
  pathways <- read_gmt(gmt_file)

  # Run Fisher's exact test for each pathway
  results <- lapply(names(pathways), function(pathway_name) {
    pathway_genes <- pathways[[pathway_name]]

    # Overlap
    sig_in_pathway <- length(intersect(sig_genes, pathway_genes))
    sig_not_in_pathway <- length(setdiff(sig_genes, pathway_genes))
    bg_in_pathway <- length(intersect(background_genes, pathway_genes)) - sig_in_pathway
    bg_not_in_pathway <- length(background_genes) - sig_in_pathway - sig_not_in_pathway - bg_in_pathway

    # Contingency table
    mat <- matrix(c(sig_in_pathway, sig_not_in_pathway, bg_in_pathway, bg_not_in_pathway),
                  nrow = 2)

    if (any(mat < 0)) return(NULL)

    test <- fisher.test(mat, alternative = "greater")

    data.frame(
      pathway = pathway_name,
      n_genes = length(pathway_genes),
      n_overlap = sig_in_pathway,
      p.value = test$p.value,
      odds_ratio = test$estimate,
      genes = paste(intersect(sig_genes, pathway_genes), collapse = "/"),
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results)
  results <- results[!is.null(results), ]

  if (nrow(results) == 0) return(NULL)

  # Adjust p-values
  results$p.adjust <- p.adjust(results$p.value, method = "BH")

  # Sort by adjusted p-value
  results <- results[order(results$p.adjust), ]

  return(results)
}

#' Read GMT file
#'
#' @param gmt_file Path to GMT file
#' @return Named list of gene sets
read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)

  pathways <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      pathway_name <- parts[1]
      genes <- parts[3:length(parts)]
      genes <- genes[genes != ""]
      pathways[[pathway_name]] <- genes
    }
  }

  return(pathways)
}

#' Run Gene Set Enrichment Analysis
#'
#' @param de_table Differential expression table
#' @param organism Organism name
#' @param config Configuration list
#' @return GSEA results
run_gsea <- function(de_table, organism, config) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    log_message("fgsea package not available. Skipping GSEA.")
    return(NULL)
  }

  # Create ranked list
  ranking_method <- config$pathway$gsea_ranking %||% "t_statistic"

  de_ranked <- de_table[!is.na(de_table$gene_symbol), ]

  if (ranking_method == "t_statistic") {
    ranks <- de_ranked$t
  } else if (ranking_method == "log2fc") {
    ranks <- de_ranked$log2FC
  } else if (ranking_method == "signed_pvalue") {
    ranks <- -log10(de_ranked$P.Value) * sign(de_ranked$log2FC)
  } else {
    ranks <- de_ranked$t
  }

  names(ranks) <- de_ranked$gene_symbol

  # Remove duplicates (keep gene with highest absolute rank)
  ranks <- ranks[order(abs(ranks), decreasing = TRUE)]
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)

  if (length(ranks) < 50) {
    log_message("Too few ranked genes for GSEA (", length(ranks), ")")
    return(NULL)
  }

  log_message("Running GSEA with ", length(ranks), " ranked genes")

  # Get gene sets
  pathways <- get_gene_sets(organism, config)

  if (is.null(pathways) || length(pathways) == 0) {
    log_message("No gene sets available for GSEA")
    return(NULL)
  }

  min_size <- config$pathway$min_set_size %||% 10
  max_size <- config$pathway$max_set_size %||% 500
  nperm <- config$pathway$gsea_nperm %||% 10000

  tryCatch({
    gsea_result <- fgsea::fgsea(
      pathways = pathways,
      stats = ranks,
      minSize = min_size,
      maxSize = max_size,
      nperm = nperm
    )

    # Convert to data frame
    result_df <- as.data.frame(gsea_result)

    # Clean up leadingEdge column
    result_df$leadingEdge <- sapply(result_df$leadingEdge, function(x) paste(x, collapse = "/"))

    # Sort by p-value
    result_df <- result_df[order(result_df$pval), ]

    return(result_df)
  }, error = function(e) {
    log_message("GSEA failed: ", e$message)
    return(NULL)
  })
}

#' Get gene sets for GSEA
#'
#' @param organism Organism name
#' @param config Configuration list
#' @return Named list of gene sets
get_gene_sets <- function(organism, config) {
  # First check for custom GMT
  if (!is.null(config$optional_inputs$gene_set_gmt) &&
      file.exists(config$optional_inputs$gene_set_gmt)) {
    return(read_gmt(config$optional_inputs$gene_set_gmt))
  }

  # Try MSigDB for human
  if (tolower(organism) %in% c("homo sapiens", "human")) {
    if (requireNamespace("msigdbr", quietly = TRUE)) {
      msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
      pathways <- split(msig$gene_symbol, msig$gs_name)
      return(pathways)
    }
  }

  log_message("No gene sets found. Consider providing a custom GMT file.")
  return(NULL)
}

#' Run analysis on custom GMT file
#'
#' @param de_table Differential expression table
#' @param config Configuration list
#' @return Custom GMT analysis results
run_custom_gmt_analysis <- function(de_table, config) {
  gmt_file <- config$optional_inputs$gene_set_gmt

  if (!file.exists(gmt_file)) {
    return(NULL)
  }

  log_message("Running custom GMT analysis from: ", gmt_file)

  pathways <- read_gmt(gmt_file)

  if (length(pathways) == 0) {
    log_message("No pathways found in GMT file")
    return(NULL)
  }

  # Run both ORA and GSEA with custom gene sets
  sig_genes <- de_table$gene_symbol[de_table$significant & !is.na(de_table$gene_symbol)]
  bg_genes <- de_table$gene_symbol[!is.na(de_table$gene_symbol)]

  # ORA
  ora_results <- run_ora_custom(sig_genes, bg_genes, config)

  # GSEA with custom pathways
  if (requireNamespace("fgsea", quietly = TRUE)) {
    de_ranked <- de_table[!is.na(de_table$gene_symbol), ]
    ranks <- de_ranked$t
    names(ranks) <- de_ranked$gene_symbol
    ranks <- ranks[!duplicated(names(ranks))]
    ranks <- sort(ranks, decreasing = TRUE)

    gsea_results <- tryCatch({
      fgsea::fgsea(
        pathways = pathways,
        stats = ranks,
        minSize = 10,
        maxSize = 500,
        nperm = 10000
      )
    }, error = function(e) NULL)
  } else {
    gsea_results <- NULL
  }

  list(
    ora = ora_results,
    gsea = as.data.frame(gsea_results)
  )
}

#' Create pathway visualization plots
#'
#' @param pathway_results Pathway analysis results
#' @param config Configuration list
#' @return List of plots
create_pathway_plots <- function(pathway_results, config) {
  all_plots <- list()

  for (contrast_name in names(pathway_results)) {
    clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))
    results <- pathway_results[[contrast_name]]

    plots <- list()

    # ORA dotplot
    if (!is.null(results$ora) && nrow(results$ora) > 0) {
      ora_top <- head(results$ora[results$ora$p.adjust < 0.1, ], 20)

      if (nrow(ora_top) > 0) {
        # Calculate gene ratio
        if ("GeneRatio" %in% colnames(ora_top)) {
          ora_top$gene_ratio <- sapply(ora_top$GeneRatio, function(x) {
            parts <- strsplit(x, "/")[[1]]
            as.numeric(parts[1]) / as.numeric(parts[2])
          })
        } else if ("n_overlap" %in% colnames(ora_top) && "n_genes" %in% colnames(ora_top)) {
          ora_top$gene_ratio <- ora_top$n_overlap / ora_top$n_genes
        } else {
          ora_top$gene_ratio <- ora_top$Count / nrow(results$ora)
        }

        # Shorten pathway names
        ora_top$short_name <- substr(ora_top$Description %||% ora_top$pathway, 1, 50)

        plots$ora_dotplot <- ggplot2::ggplot(ora_top,
          ggplot2::aes(x = gene_ratio, y = reorder(short_name, gene_ratio), size = Count %||% n_overlap, color = p.adjust)) +
          ggplot2::geom_point() +
          ggplot2::scale_color_gradient(low = "red", high = "blue") +
          ggplot2::labs(
            title = paste("ORA Results:", contrast_name),
            x = "Gene Ratio",
            y = "Pathway",
            size = "Count",
            color = "Adj. P-value"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

        save_plot(plots$ora_dotplot, paste0("pathway_ora_", clean_name, ".png"), config,
                  height = 8, width = 10, subdir = "plots")
      }
    }

    # GSEA enrichment plot (NES barplot)
    if (!is.null(results$gsea) && nrow(results$gsea) > 0) {
      gsea_top <- results$gsea[results$gsea$padj < 0.1, ]
      gsea_top <- head(gsea_top[order(abs(gsea_top$NES), decreasing = TRUE), ], 20)

      if (nrow(gsea_top) > 0) {
        gsea_top$short_name <- substr(gsea_top$pathway, 1, 50)
        gsea_top$direction <- ifelse(gsea_top$NES > 0, "Up", "Down")

        plots$gsea_barplot <- ggplot2::ggplot(gsea_top,
          ggplot2::aes(x = NES, y = reorder(short_name, NES), fill = direction)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::scale_fill_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2")) +
          ggplot2::labs(
            title = paste("GSEA Results:", contrast_name),
            x = "Normalized Enrichment Score (NES)",
            y = "Pathway",
            fill = "Direction"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

        save_plot(plots$gsea_barplot, paste0("pathway_gsea_", clean_name, ".png"), config,
                  height = 8, width = 10, subdir = "plots")
      }
    }

    all_plots[[contrast_name]] <- plots
  }

  return(all_plots)
}
