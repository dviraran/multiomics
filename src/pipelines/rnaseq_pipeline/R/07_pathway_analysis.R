# R/07_pathway_analysis.R
# Functions for pathway and gene set analysis

#' Read GMT file for custom gene sets
#' @param gmt_file Path to GMT file
#' @return Named list of gene sets
read_gmt <- function(gmt_file) {
  if (!file.exists(gmt_file)) {
    stop("GMT file not found: ", gmt_file)
  }

  lines <- readLines(gmt_file)
  gene_sets <- list()

  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      set_name <- parts[1]
      # parts[2] is typically description/URL, skip it
      genes <- parts[3:length(parts)]
      genes <- genes[genes != ""]
      gene_sets[[set_name]] <- genes
    }
  }

  message("Loaded ", length(gene_sets), " gene sets from GMT file")
  gene_sets
}


#' Load gene sets from various sources
#' @param organism Organism name
#' @param pathway_database Vector of databases to use
#' @param gmt_file Optional custom GMT file
#' @param annotation Gene annotation data frame
#' @param gene_id_type Type of gene IDs
#' @return List of gene sets by database
load_gene_sets <- function(organism,
                           pathway_database = c("GO", "KEGG"),
                           gmt_file = NULL,
                           annotation = NULL,
                           gene_id_type = "ensembl_gene_id") {

  gene_sets <- list()
  org_info <- get_organism_info(organism)

  # Load custom GMT if provided
  if (!is.null(gmt_file) && file.exists(gmt_file)) {
    gene_sets$custom <- read_gmt(gmt_file)
    message("Loaded custom gene sets from: ", gmt_file)
  }

  # If organism not supported and no GMT, warn and return empty
  if (!org_info$supported && length(gene_sets) == 0) {
    warning("Organism '", organism, "' not supported for standard pathway databases. ",
            "Provide a GMT file for pathway analysis.")
    return(gene_sets)
  }

  # Try to load standard databases for supported organisms
  if (org_info$supported) {

    # GO terms
    if ("GO" %in% pathway_database) {
      tryCatch({
        if (!is.na(org_info$orgdb) && requireNamespace(org_info$orgdb, quietly = TRUE)) {
          # Load the OrgDb package and get the database object
          orgdb <- getExportedValue(org_info$orgdb, org_info$orgdb)

          # Get GO annotations
          go_bp <- AnnotationDbi::select(
            orgdb,
            keys = AnnotationDbi::keys(orgdb, keytype = "ENSEMBL"),
            columns = c("ENSEMBL", "GO"),
            keytype = "ENSEMBL"
          )

          # Filter to BP (Biological Process)
          go_bp <- go_bp[!is.na(go_bp$GO), ]

          # Convert to gene set list
          go_sets <- split(go_bp$ENSEMBL, go_bp$GO)
          go_sets <- go_sets[lengths(go_sets) >= 10 & lengths(go_sets) <= 500]

          gene_sets$GO <- go_sets
          message("Loaded ", length(go_sets), " GO gene sets")
        }
      }, error = function(e) {
        warning("Failed to load GO terms: ", e$message)
      })
    }

    # KEGG pathways
    if ("KEGG" %in% pathway_database && !is.na(org_info$kegg)) {
      tryCatch({
        # Get KEGG pathway-to-gene mapping
        kegg_gene <- clusterProfiler::download_KEGG(org_info$kegg, keggType = "KEGG")

        if (!is.null(kegg_gene) && nrow(kegg_gene$KEGGPATHID2EXTID) > 0) {
          kegg_df <- kegg_gene$KEGGPATHID2EXTID

          # Need to convert Entrez to Ensembl if annotation available
          if (!is.null(annotation) && "entrez_id" %in% colnames(annotation)) {
            entrez_to_ensembl <- annotation[!is.na(annotation$entrez_id),
                                            c("gene_id", "entrez_id")]
            kegg_df <- merge(kegg_df, entrez_to_ensembl,
                             by.x = "to", by.y = "entrez_id")

            kegg_sets <- split(kegg_df$gene_id, kegg_df$from)
          } else {
            # Use Entrez IDs directly
            kegg_sets <- split(kegg_df$to, kegg_df$from)
          }

          kegg_sets <- kegg_sets[lengths(kegg_sets) >= 10 & lengths(kegg_sets) <= 500]
          gene_sets$KEGG <- kegg_sets
          message("Loaded ", length(kegg_sets), " KEGG pathways")
        }
      }, error = function(e) {
        warning("Failed to load KEGG pathways: ", e$message)
      })
    }

    # Reactome
    if ("Reactome" %in% pathway_database) {
      tryCatch({
        if (requireNamespace("reactome.db", quietly = TRUE)) {
          # Implementation for Reactome
          message("Reactome support: checking availability...")
          # This would require reactome.db package
        }
      }, error = function(e) {
        warning("Reactome not available: ", e$message)
      })
    }
  }

  if (length(gene_sets) == 0) {
    warning("No gene sets loaded. Pathway analysis will be skipped.")
  }

  gene_sets
}


#' Run pathway analysis using fgsea or ORA
#' @param de_results List of DE results
#' @param gene_sets Gene sets from load_gene_sets
#' @param annotation Gene annotation
#' @param method Analysis method: "fgsea", "ora", or "both"
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @return List of pathway results by contrast and database
run_pathway_analysis <- function(de_results,
                                  gene_sets,
                                  annotation = NULL,
                                  method = "fgsea",
                                  min_size = 10,
                                  max_size = 500) {

  if (length(gene_sets) == 0) {
    message("No gene sets available. Skipping pathway analysis.")
    return(list())
  }

  pathway_results <- list()

  for (contrast_name in names(de_results)) {
    res <- de_results[[contrast_name]]

    if (is.null(res) || nrow(res) == 0) next

    message("Running pathway analysis for: ", contrast_name)
    contrast_results <- list()

    for (db_name in names(gene_sets)) {
      gs <- gene_sets[[db_name]]

      if (length(gs) == 0) next

      message("  Database: ", db_name, " (", length(gs), " gene sets)")

      tryCatch({
        if (method %in% c("fgsea", "both")) {
          # GSEA-style analysis with fgsea

          # Create ranked gene list (use stat or compute from log2FC and pvalue)
          if ("stat" %in% colnames(res)) {
            ranks <- setNames(res$stat, res$gene_id)
          } else {
            # Alternative: sign(log2FC) * -log10(pvalue)
            ranks <- setNames(
              sign(res$log2FoldChange) * -log10(res$pvalue + 1e-300),
              res$gene_id
            )
          }

          ranks <- ranks[!is.na(ranks)]
          ranks <- sort(ranks, decreasing = TRUE)

          # Run fgsea
          fgsea_res <- fgsea::fgsea(
            pathways = gs,
            stats = ranks,
            minSize = min_size,
            maxSize = max_size,
            nPermSimple = 10000
          )

          fgsea_df <- as.data.frame(fgsea_res)
          fgsea_df$contrast <- contrast_name
          fgsea_df$database <- db_name
          fgsea_df$method <- "fgsea"

          # Convert leading edge genes to comma-separated string
          fgsea_df$leadingEdge <- sapply(fgsea_df$leadingEdge, paste, collapse = ",")

          contrast_results[[paste0(db_name, "_fgsea")]] <- fgsea_df

          n_sig <- sum(fgsea_df$padj < 0.05, na.rm = TRUE)
          message("    fgsea: ", n_sig, " significant pathways (padj < 0.05)")
        }

        if (method %in% c("ora", "both")) {
          # Over-representation analysis

          # Get significant genes
          sig_genes <- res$gene_id[res$significant & res$direction == "up"]
          sig_genes_down <- res$gene_id[res$significant & res$direction == "down"]
          background <- res$gene_id

          # Run ORA for upregulated
          if (length(sig_genes) >= 5) {
            ora_up <- run_ora(sig_genes, gs, background, min_size, max_size)
            ora_up$contrast <- contrast_name
            ora_up$database <- db_name
            ora_up$method <- "ora"
            ora_up$direction <- "up"
            contrast_results[[paste0(db_name, "_ora_up")]] <- ora_up
          }

          # Run ORA for downregulated
          if (length(sig_genes_down) >= 5) {
            ora_down <- run_ora(sig_genes_down, gs, background, min_size, max_size)
            ora_down$contrast <- contrast_name
            ora_down$database <- db_name
            ora_down$method <- "ora"
            ora_down$direction <- "down"
            contrast_results[[paste0(db_name, "_ora_down")]] <- ora_down
          }
        }

      }, error = function(e) {
        warning("Pathway analysis failed for ", db_name, ": ", e$message)
      })
    }

    pathway_results[[contrast_name]] <- contrast_results
  }

  pathway_results
}


#' Run over-representation analysis
#' @param sig_genes Significant gene IDs
#' @param gene_sets List of gene sets
#' @param background Background gene IDs
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @return Data frame of ORA results
run_ora <- function(sig_genes, gene_sets, background, min_size = 10, max_size = 500) {

  # Filter gene sets by size
  gs_sizes <- lengths(gene_sets)
  gs_filtered <- gene_sets[gs_sizes >= min_size & gs_sizes <= max_size]

  if (length(gs_filtered) == 0) {
    return(data.frame())
  }

  # Run Fisher's exact test for each gene set
  results <- lapply(names(gs_filtered), function(gs_name) {
    gs_genes <- gs_filtered[[gs_name]]

    # Intersect with background
    gs_genes <- intersect(gs_genes, background)
    sig_in_gs <- length(intersect(sig_genes, gs_genes))
    sig_not_gs <- length(sig_genes) - sig_in_gs
    gs_not_sig <- length(gs_genes) - sig_in_gs
    neither <- length(background) - length(sig_genes) - gs_not_sig

    # Fisher's exact test
    mat <- matrix(c(sig_in_gs, gs_not_sig, sig_not_gs, neither), nrow = 2)

    if (any(mat < 0)) return(NULL)

    test <- fisher.test(mat, alternative = "greater")

    data.frame(
      pathway = gs_name,
      size = length(gs_genes),
      overlap = sig_in_gs,
      pvalue = test$p.value,
      odds_ratio = test$estimate,
      stringsAsFactors = FALSE
    )
  })

  results_df <- do.call(rbind, results[!sapply(results, is.null)])

  if (nrow(results_df) > 0) {
    results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
    results_df <- results_df[order(results_df$pvalue), ]
  }

  results_df
}


#' Save pathway results to files
#' @param pathway_results Pathway analysis results
#' @param output_dir Output directory
#' @return Vector of output file paths
save_pathway_results <- function(pathway_results, output_dir = "outputs/pathway_results") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_files <- c()

  for (contrast_name in names(pathway_results)) {
    contrast_res <- pathway_results[[contrast_name]]

    for (analysis_name in names(contrast_res)) {
      res_df <- contrast_res[[analysis_name]]

      if (is.null(res_df) || nrow(res_df) == 0) next

      clean_contrast <- gsub("[^a-zA-Z0-9_-]", "_", contrast_name)
      clean_analysis <- gsub("[^a-zA-Z0-9_-]", "_", analysis_name)

      output_file <- file.path(output_dir,
                               paste0("pathway_", clean_contrast, "_", clean_analysis, ".csv"))

      write.csv(res_df, output_file, row.names = FALSE)
      output_files <- c(output_files, output_file)
    }
  }

  message("Saved ", length(output_files), " pathway result file(s) to ", output_dir)

  output_files
}


#' Generate pathway visualization plots
#' @param pathway_results Pathway analysis results
#' @param output_dir Output directory
#' @return List of plot file paths
generate_pathway_plots <- function(pathway_results, output_dir = "outputs/plots") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  plot_files <- list()

  for (contrast_name in names(pathway_results)) {
    contrast_res <- pathway_results[[contrast_name]]
    clean_contrast <- gsub("[^a-zA-Z0-9_-]", "_", contrast_name)

    for (analysis_name in names(contrast_res)) {
      res_df <- contrast_res[[analysis_name]]

      if (is.null(res_df) || nrow(res_df) == 0) next

      clean_analysis <- gsub("[^a-zA-Z0-9_-]", "_", analysis_name)

      # Get top pathways
      if ("padj" %in% colnames(res_df)) {
        top_pathways <- head(res_df[order(res_df$padj), ], 20)
      } else if ("pvalue" %in% colnames(res_df)) {
        top_pathways <- head(res_df[order(res_df$pvalue), ], 20)
      } else {
        next
      }

      if (nrow(top_pathways) < 3) next

      # Truncate long pathway names
      top_pathways$pathway_short <- substr(top_pathways$pathway, 1, 50)

      # Dotplot
      if ("NES" %in% colnames(top_pathways)) {
        # fgsea results
        p <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_short, NES))) +
          geom_point(aes(size = size, color = -log10(padj))) +
          scale_color_gradient(low = "blue", high = "red") +
          labs(
            title = paste("Top Pathways:", contrast_name),
            subtitle = analysis_name,
            x = "Normalized Enrichment Score",
            y = ""
          ) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 8))
      } else {
        # ORA results
        p <- ggplot(top_pathways, aes(x = -log10(pvalue), y = reorder(pathway_short, -log10(pvalue)))) +
          geom_point(aes(size = overlap, color = -log10(padj))) +
          scale_color_gradient(low = "blue", high = "red") +
          labs(
            title = paste("Top Pathways:", contrast_name),
            subtitle = analysis_name,
            x = "-log10(p-value)",
            y = ""
          ) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 8))
      }

      plot_file <- file.path(output_dir, paste0("pathway_", clean_contrast, "_", clean_analysis, ".png"))
      ggsave(plot_file, p, width = 10, height = 8)
      plot_files[[paste0(clean_contrast, "_", clean_analysis)]] <- plot_file
    }
  }

  message("Generated pathway plots in ", output_dir)

  plot_files
}
