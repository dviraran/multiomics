#!/usr/bin/env Rscript
# Extract top proteomics features from DIABLO components for enrichment

# Load necessary data
config <- yaml::read_yaml("config.yml")
prot_loadings <- read.csv("outputs/tables/diablo_loadings_proteomics.csv",
                          stringsAsFactors = FALSE)
gene_prot_map <- read.csv(config$global$gene_protein_mapping,
                          stringsAsFactors = FALSE)

message("Loaded ", nrow(prot_loadings), " protein loadings")
message("Loaded ", nrow(gene_prot_map), " gene-protein mappings")

# Function to extract and save top features
extract_top_features <- function(loadings_df, comp_col, gene_map, n_top = 100) {
  comp_name <- comp_col
  message("\nProcessing ", comp_name, "...")

  # Get top features by absolute loading
  loadings_df$abs_loading <- abs(loadings_df[[comp_col]])
  loadings_df$signed_loading <- loadings_df[[comp_col]]
  loadings_df <- loadings_df[order(loadings_df$abs_loading,
                                   decreasing = TRUE), ]
  top_features <- head(loadings_df, n_top)

  # Try to map to gene IDs
  # Check which column to use for merging
  if ("uniprot_id" %in% colnames(gene_map)) {
    top_features <- merge(top_features, gene_map,
                         by.x = "feature_id",
                         by.y = "uniprot_id",
                         all.x = TRUE)
  } else if ("protein_id" %in% colnames(gene_map)) {
    top_features <- merge(top_features, gene_map,
                         by.x = "feature_id",
                         by.y = "protein_id",
                         all.x = TRUE)
  }

  # Save results
  out_file <- paste0("outputs/tables/diablo_proteomics_top_features_",
                    comp_name, ".csv")
  write.csv(top_features, out_file, row.names = FALSE)
  message("  Saved ", out_file)

  # Also save just gene IDs for web enrichment
  if ("gene_id" %in% colnames(top_features)) {
    genes <- unique(top_features$gene_id[!is.na(top_features$gene_id)])
    gene_file <- paste0("outputs/tables/diablo_proteomics_genes_",
                       comp_name, ".txt")
    writeLines(genes, gene_file)
    message("  Saved ", length(genes), " genes to ", gene_file)
    message("  You can upload this file to g:Profiler, Enrichr, or MetaScape")
  }

  return(top_features)
}

# Process both components
comp1_features <- extract_top_features(prot_loadings, "comp1", gene_prot_map)
comp2_features <- extract_top_features(prot_loadings, "comp2", gene_prot_map)

message("\n=== Extraction complete ===")
message("Next steps:")
message("1. Upload gene lists to web enrichment tools:")
message("   - g:Profiler: https://biit.cs.ut.ee/gprofiler/")
message("   - Enrichr: https://maayanlab.cloud/Enrichr/")
message("   - MetaScape: https://metascape.org/")
message("2. Or install clusterProfiler and run enrichment locally")
