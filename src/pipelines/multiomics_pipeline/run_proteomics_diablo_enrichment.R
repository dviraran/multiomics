#!/usr/bin/env Rscript
# Run proteomics enrichment for DIABLO components

library(clusterProfiler)
library(ggplot2)

# Source utility functions
source("R/00_utils.R")

# Load configuration
config <- yaml::read_yaml("config.yml")

# Load proteomics DIABLO loadings
prot_loadings <- read.csv("outputs/tables/diablo_loadings_proteomics.csv", stringsAsFactors = FALSE)

# Load gene-protein mapping
gene_prot_map <- read.csv(config$global$gene_protein_mapping, stringsAsFactors = FALSE)

log_message("Loaded ", nrow(prot_loadings), " protein loadings")
log_message("Loaded ", nrow(gene_prot_map), " gene-protein mappings")

# Set up organism database
organism_name <- config$global$organism %||% "c_elegans"
if (organism_name == "c_elegans") {
  if (!requireNamespace("org.Ce.eg.db", quietly = TRUE)) {
    stop("org.Ce.eg.db package required for C. elegans enrichment")
  }
  org_db <- org.Ce.eg.db::org.Ce.eg.db
  kegg_code <- "cel"
} else {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db package required for human enrichment")
  }
  org_db <- org.Hs.eg.db::org.Hs.eg.db
  kegg_code <- "hsa"
}

# Function to run enrichment for one component
run_component_enrichment <- function(loadings_df, comp_col, gene_map,
                                    org_db, kegg_code, config) {
  comp_name <- comp_col
  log_message("Processing ", comp_name, "...")

  # Get top 100 features by absolute loading
  loadings_df$abs_loading <- abs(loadings_df[[comp_col]])
  loadings_df <- loadings_df[order(loadings_df$abs_loading,
                                   decreasing = TRUE), ]
  top_proteins <- head(loadings_df$feature_id, 100)

  log_message("  Top 100 proteins selected")

  # Map to gene IDs
  if ("protein_id" %in% colnames(gene_map)) {
    mapped_rows <- gene_map[gene_map$protein_id %in% top_proteins, ]
    mapped <- unique(mapped_rows$gene_id)
  } else {
    mapped <- character(0)
  }

  if (length(mapped) == 0) {
    log_message("  No protein IDs mapped. Trying direct Wormbase ID...")
    # Try using protein IDs directly as they might be gene IDs
    mapped <- top_proteins[grepl("^WBGene", top_proteins)]
  }

  log_message("  Mapped to ", length(mapped), " gene IDs")

  if (length(mapped) < 5) {
    log_message("  Too few genes mapped (", length(mapped), "), skipping enrichment")
    return(NULL)
  }

  # Determine key type
  key_type <- "SYMBOL"
  if (all(grepl("^WBGene", mapped[1:min(5, length(mapped))]))) {
    key_type <- "WORMBASE"
  } else if (all(grepl("^ENS", mapped[1:min(5, length(mapped))]))) {
    key_type <- "ENSEMBL"
  }

  log_message("  Using keyType: ", key_type)

  # Convert to Entrez IDs
  tryCatch({
    gene_map_entrez <- clusterProfiler::bitr(
      mapped,
      fromType = key_type,
      toType = "ENTREZID",
      OrgDb = org_db
    )

    if (is.null(gene_map_entrez) || nrow(gene_map_entrez) < 5) {
      log_message("  Too few genes converted to Entrez IDs")
      return(NULL)
    }

    genes <- gene_map_entrez$ENTREZID
    log_message("  Converted to ", length(genes), " Entrez IDs")

    # Run GO Enrichment (Biological Process)
    tryCatch({
      c_go <- clusterProfiler::enrichGO(
        gene = genes,
        OrgDb = org_db,
        ont = "BP",
        pvalueCutoff = 0.05
      )

      if (!is.null(c_go) && nrow(as.data.frame(c_go)) > 0) {
        res_go <- as.data.frame(c_go)
        res_go$block <- "proteomics"
        res_go$component <- comp_name
        res_go$type <- "GO"

        out_name <- paste0("diablo_enrichment_go_proteomics_", comp_name)
        save_table(res_go, paste0(out_name, ".csv"), config)
        log_message("  Saved GO enrichment: ", out_name, ".csv (", nrow(res_go), " terms)")

        # Create plot
        p <- clusterProfiler::dotplot(c_go, showCategory = 15) +
          ggtitle(paste("DIABLO Proteomics", comp_name, "- GO Enrichment"))

        save_plot(p, paste0(out_name, ".png"), config, width = 10, height = 8)
        log_message("  Saved plot: ", out_name, ".png")

        return(res_go)
      } else {
        log_message("  No significant GO terms found")
      }
    }, error = function(e) {
      log_message("  GO enrichment failed: ", e$message)
    })

    # Run KEGG Enrichment
    tryCatch({
      c_kegg <- clusterProfiler::enrichKEGG(
        gene = genes,
        organism = kegg_code,
        pvalueCutoff = 0.05
      )

      if (!is.null(c_kegg) && nrow(as.data.frame(c_kegg)) > 0) {
        res_kegg <- as.data.frame(c_kegg)
        res_kegg$block <- "proteomics"
        res_kegg$component <- comp_name
        res_kegg$type <- "KEGG"

        out_name <- paste0("diablo_enrichment_kegg_proteomics_", comp_name)
        save_table(res_kegg, paste0(out_name, ".csv"), config)
        log_message("  Saved KEGG enrichment: ", out_name, ".csv (", nrow(res_kegg), " pathways)")

        # Create plot
        p <- clusterProfiler::dotplot(c_kegg, showCategory = 15) +
          ggtitle(paste("DIABLO Proteomics", comp_name, "- KEGG Enrichment"))

        save_plot(p, paste0(out_name, ".png"), config, width = 10, height = 8)
        log_message("  Saved plot: ", out_name, ".png")
      } else {
        log_message("  No significant KEGG pathways found")
      }
    }, error = function(e) {
      log_message("  KEGG enrichment failed: ", e$message)
    })

  }, error = function(e) {
    log_message("  Gene ID conversion failed: ", e$message)
  })

  return(NULL)
}

# Run enrichment for each component
log_message("=== Running Proteomics DIABLO Enrichment ===")

for (comp_col in c("comp1", "comp2")) {
  run_component_enrichment(prot_loadings, comp_col, gene_prot_map, org_db, kegg_code, config)
}

log_message("=== Proteomics enrichment complete ===")
