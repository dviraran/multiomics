# =============================================================================
# Test script for KEGG Pathview with Amir_Sapir data
# =============================================================================
# RNA data: C. elegans gene symbols with log2FC
# Metabolomics data: HMDB IDs with log2FC
# =============================================================================

# Source the main pathview script
source("/home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/15_kegg_pathview.R")

# Install and load required packages
install_load_packages()

# Also need C. elegans annotation package
if (!requireNamespace("org.Ce.eg.db", quietly = TRUE)) {
    message("Installing C. elegans annotation package...")
    BiocManager::install("org.Ce.eg.db", ask = FALSE, update = FALSE)
}
library(org.Ce.eg.db)

# -----------------------------------------------------------------------------
# Define paths
# -----------------------------------------------------------------------------
rna_file <- "/home/ozsol/multiAnalysis/projects/Amir_Sapir/RNA_12.5_vs_0.tab"
metab_file <- "/home/ozsol/multiAnalysis/projects/Amir_Sapir/Metabolomics_12.5_vs_0.tab"
output_dir <- "/home/ozsol/multiAnalysis/projects/Amir_Sapir/pathview_output"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load and prepare gene data (C. elegans gene symbols -> Entrez IDs)
# -----------------------------------------------------------------------------
message("\n=== Loading RNA data ===")
rna_df <- readr::read_tsv(rna_file, show_col_types = FALSE)
message(sprintf("Read %d genes from RNA file", nrow(rna_df)))

# Map C. elegans gene symbols to Entrez IDs
gene_symbols <- rna_df$GeneSymbol
gene_fc_values <- rna_df$log2FC

# Try mapping with org.Ce.eg.db
entrez_ids <- tryCatch(
    {
        AnnotationDbi::mapIds(
            org.Ce.eg.db,
            keys = gene_symbols,
            column = "ENTREZID",
            keytype = "SYMBOL",
            multiVals = "first"
        )
    },
    error = function(e) {
        message("Direct SYMBOL mapping failed, trying GENENAME...")
        # Try with GENENAME as some C. elegans genes use different naming
        tryCatch(
            {
                AnnotationDbi::mapIds(
                    org.Ce.eg.db,
                    keys = gene_symbols,
                    column = "ENTREZID",
                    keytype = "GENENAME",
                    multiVals = "first"
                )
            },
            error = function(e2) {
                message("GENENAME mapping also failed: ", e2$message)
                rep(NA, length(gene_symbols))
            }
        )
    }
)

# Create gene fold-change vector
valid_idx <- !is.na(entrez_ids)
gene_fc <- setNames(gene_fc_values[valid_idx], entrez_ids[valid_idx])
gene_fc <- gene_fc[!duplicated(names(gene_fc))]

message(sprintf(
    "Successfully mapped %d/%d genes to Entrez IDs",
    length(gene_fc), length(gene_symbols)
))

# -----------------------------------------------------------------------------
# Load and prepare metabolite data (HMDB -> KEGG)
# -----------------------------------------------------------------------------
message("\n=== Loading Metabolomics data ===")
metab_df <- readr::read_tsv(metab_file, show_col_types = FALSE)
message(sprintf("Read %d metabolites from metabolomics file", nrow(metab_df)))

# Use the new HMDB mapping function
cpd_fc <- load_metabolite_data(
    file_path = metab_file,
    id_col = "HMDB",
    fc_col = "log2FC",
    id_type = "hmdb"
)

message(sprintf("Successfully mapped %d metabolites to KEGG IDs", length(cpd_fc)))

# -----------------------------------------------------------------------------
# Test with a common metabolic pathway
# For C. elegans, use species code "cel"
# -----------------------------------------------------------------------------
message("\n=== Running pathview ===")

# Change to output directory
setwd(output_dir)

# Test pathways for C. elegans:
# cel00010 - Glycolysis / Gluconeogenesis
# cel00020 - TCA cycle
# cel01100 - Metabolic pathways (overview)

test_pathway <- "cel00010" # Glycolysis

message(sprintf("\nGenerating pathway visualization for: %s", test_pathway))
message(sprintf("Genes with data: %d", length(gene_fc)))
message(sprintf("Metabolites with data: %d", length(cpd_fc)))

# Show some example data
message("\nTop 5 genes by absolute fold-change:")
top_genes <- head(gene_fc[order(abs(gene_fc), decreasing = TRUE)], 5)
print(top_genes)

message("\nTop 5 metabolites by absolute fold-change:")
top_cpds <- head(cpd_fc[order(abs(cpd_fc), decreasing = TRUE)], 5)
print(top_cpds)

# Run pathview
pv_result <- tryCatch(
    {
        pathview::pathview(
            gene.data = gene_fc,
            cpd.data = cpd_fc,
            pathway.id = test_pathway,
            species = "cel",
            gene.idtype = "entrez",
            cpd.idtype = "kegg",
            kegg.native = TRUE,
            limit = list(gene = 2, cpd = 2),
            low = c("#3366CC", "#3366CC"),
            mid = c("#CCCCCC", "#CCCCCC"),
            high = c("#CC3333", "#CC3333"),
            out.suffix = "FC_overlay",
            na.col = "transparent"
        )
    },
    error = function(e) {
        message("pathview error: ", e$message)
        NULL
    }
)

# Report results
message("\n=== Results ===")
output_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
if (length(output_files) > 0) {
    message("Generated PNG files:")
    for (f in output_files) {
        message("  ", f)
    }
} else {
    message("No PNG files generated. Check for errors above.")
}

message("\nDone!")
