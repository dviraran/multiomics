# =============================================================================
# Download TCGA Breast Cancer Multi-Omics Data (from mixOmics)
# =============================================================================
# Downloads the breast.TCGA dataset from the mixOmics Bioconductor package.
# This dataset contains REAL gene symbols that work with pathway analysis.
#
# Data contents:
# - mRNA expression (200 genes, gene symbols)
# - miRNA expression (184 miRNAs)
# - Proteomics/RPPA (142 proteins, gene symbols)
# - 3 breast cancer subtypes: Basal, Her2, LumA
#
# Reference: TCGA Network (2012) Nature
#
# Usage: source("scripts/download_breast_tcga_data.R") (from examples/ directory)
# =============================================================================

# Set CRAN mirror if not set
if (getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Output directory
data_dir <- "data/breast_tcga"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

cat("=== Downloading TCGA Breast Cancer Multi-Omics Data ===\n\n")

# Try to load breast.TCGA data
breast.TCGA <- NULL

# Method 1: Try loading from mixOmics if installed
if (requireNamespace("mixOmics", quietly = TRUE)) {
  cat("Loading breast.TCGA from mixOmics package...\n")
  data("breast.TCGA", package = "mixOmics")
} else {
  # Method 2: Download and extract from Bioconductor
  cat("mixOmics not installed. Downloading data directly from Bioconductor...\n")

  temp_dir <- tempdir()
  pkg_file <- file.path(temp_dir, "mixomics.tar.gz")
  rda_file <- file.path(temp_dir, "mixOmics", "data", "breast.TCGA.rda")

  # Check if already downloaded
  if (!file.exists(rda_file)) {
    cat("  Downloading mixOmics package (18MB)...\n")
    download.file(
      "https://bioconductor.org/packages/3.18/bioc/src/contrib/mixOmics_6.26.0.tar.gz",
      pkg_file,
      mode = "wb",
      quiet = TRUE
    )

    cat("  Extracting data...\n")
    untar(pkg_file, exdir = temp_dir)
  }

  # Load the data
  cat("  Loading breast.TCGA.rda...\n")
  load(rda_file)
}

# Verify data loaded
if (!exists("breast.TCGA") || is.null(breast.TCGA)) {
  stop("Failed to load breast.TCGA data")
}

cat("Data loaded successfully!\n\n")

# -----------------------------------------------------------------------------
# Extract and save mRNA expression data
# -----------------------------------------------------------------------------
cat("Processing mRNA expression data...\n")

mrna_data <- breast.TCGA$data.train$mrna
cat("  Original dimensions:", nrow(mrna_data), "samples x", ncol(mrna_data), "genes\n")

# Transpose to genes x samples format (standard for our pipelines)
mrna_matrix <- t(mrna_data)

# For multi-omics: save as-is (normalized values)
write.csv(mrna_matrix, file.path(data_dir, "mrna_normalized_matrix.csv"))
cat("  Saved normalized mRNA matrix:", nrow(mrna_matrix), "genes x", ncol(mrna_matrix), "samples\n")

# For RNA-seq pipeline: create pseudo-counts (exponentiating and rounding)
# The original data appears to be log2-transformed
mrna_counts <- round(2^(mrna_matrix + 5))  # Shift to positive range and exponentiate
mrna_counts[mrna_counts < 0] <- 0
write.csv(mrna_counts, file.path(data_dir, "mrna_counts_matrix.csv"))
cat("  Saved pseudo-counts matrix for RNA-seq pipeline\n")

# Show some gene names to confirm they're real
cat("  Sample gene symbols:", paste(head(rownames(mrna_matrix), 10), collapse = ", "), "...\n")

# -----------------------------------------------------------------------------
# Extract and save miRNA expression data
# -----------------------------------------------------------------------------
cat("\nProcessing miRNA expression data...\n")

mirna_data <- breast.TCGA$data.train$mirna
mirna_matrix <- t(mirna_data)

write.csv(mirna_matrix, file.path(data_dir, "mirna_normalized_matrix.csv"))
cat("  Saved miRNA matrix:", nrow(mirna_matrix), "miRNAs x", ncol(mirna_matrix), "samples\n")
cat("  Sample miRNA IDs:", paste(head(rownames(mirna_matrix), 5), collapse = ", "), "...\n")

# -----------------------------------------------------------------------------
# Extract and save proteomics (RPPA) data
# -----------------------------------------------------------------------------
cat("\nProcessing proteomics (RPPA) data...\n")

protein_data <- breast.TCGA$data.train$protein
protein_matrix <- t(protein_data)

write.csv(protein_matrix, file.path(data_dir, "protein_intensity_matrix.csv"))
cat("  Saved protein matrix:", nrow(protein_matrix), "proteins x", ncol(protein_matrix), "samples\n")
cat("  Sample protein IDs:", paste(head(rownames(protein_matrix), 10), collapse = ", "), "...\n")

# Create a protein ID mapping file (RPPA name to gene symbol where possible)
protein_names <- rownames(protein_matrix)
protein_mapping <- data.frame(
  protein_id = protein_names,
  # Extract base gene symbol (remove phospho suffixes)
  gene_symbol = gsub("_p[STY][0-9]+.*$", "", protein_names),
  # Flag if it's a phospho-protein
  is_phospho = grepl("_p[STY]", protein_names),
  stringsAsFactors = FALSE
)
# Clean up common patterns
protein_mapping$gene_symbol <- gsub("-", "", protein_mapping$gene_symbol)
protein_mapping$gene_symbol <- gsub("_", "", protein_mapping$gene_symbol)

write.csv(protein_mapping, file.path(data_dir, "protein_id_mapping.csv"), row.names = FALSE)
cat("  Saved protein ID mapping file\n")

# -----------------------------------------------------------------------------
# Create metadata file
# -----------------------------------------------------------------------------
cat("\nCreating metadata file...\n")

subtype <- breast.TCGA$data.train$subtype
sample_ids <- rownames(breast.TCGA$data.train$mrna)

metadata <- data.frame(
  sample_id = sample_ids,
  subtype = as.character(subtype),
  # Create binary condition for simple differential analysis
  condition = ifelse(subtype == "Basal", "basal",
              ifelse(subtype == "Her2", "her2", "luminal")),
  stringsAsFactors = FALSE
)

# Add sample counts per subtype
cat("  Subtype distribution:\n")
print(table(metadata$subtype))

write.csv(metadata, file.path(data_dir, "metadata.csv"), row.names = FALSE)
cat("  Saved metadata for", nrow(metadata), "samples\n")

# -----------------------------------------------------------------------------
# Create gene annotation file (for pathway analysis)
# -----------------------------------------------------------------------------
cat("\nCreating gene annotation file...\n")

gene_symbols <- rownames(mrna_matrix)

# Create basic annotation
gene_annotation <- data.frame(
  gene_symbol = gene_symbols,
  entrez_id = NA_character_,
  stringsAsFactors = FALSE
)

# Try to get Entrez IDs using org.Hs.eg.db if available
if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)) {
  tryCatch({
    entrez_map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                         keys = gene_symbols,
                                         columns = "ENTREZID",
                                         keytype = "SYMBOL")
    # Merge
    for (i in seq_len(nrow(gene_annotation))) {
      match_idx <- which(entrez_map$SYMBOL == gene_annotation$gene_symbol[i])
      if (length(match_idx) > 0) {
        gene_annotation$entrez_id[i] <- entrez_map$ENTREZID[match_idx[1]]
      }
    }
    cat("  Mapped", sum(!is.na(gene_annotation$entrez_id)), "of",
        nrow(gene_annotation), "genes to Entrez IDs\n")
  }, error = function(e) {
    cat("  Note: Could not map to Entrez IDs:", conditionMessage(e), "\n")
  })
} else {
  cat("  Note: org.Hs.eg.db not available, skipping Entrez ID mapping\n")
}

write.csv(gene_annotation, file.path(data_dir, "gene_annotation.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=== Download Complete ===\n")
cat("Data saved to:", normalizePath(data_dir), "\n\n")

cat("Files created:\n")
for (f in list.files(data_dir, full.names = FALSE)) {
  cat("  -", f, "\n")
}

cat("\nDataset summary:\n")
cat("  - 150 samples (breast cancer patients)\n")
cat("  - 3 subtypes: Basal (", sum(metadata$subtype == "Basal"),
    "), Her2 (", sum(metadata$subtype == "Her2"),
    "), LumA (", sum(metadata$subtype == "LumA"), ")\n", sep = "")
cat("  - mRNA:", nrow(mrna_matrix), "genes (real gene symbols)\n")
cat("  - miRNA:", nrow(mirna_matrix), "miRNAs\n")
cat("  - Proteins:", nrow(protein_matrix), "(RPPA, gene symbol-based)\n")

cat("\nThis dataset uses REAL gene symbols that will work with:\n")
cat("  - Pathway analysis (GO, KEGG, Reactome)\n")
cat("  - Gene set enrichment analysis (GSEA)\n")
cat("  - ID mapping to Entrez/Ensembl\n")

cat("\n\nNext steps:\n")
cat("1. Run: source('scripts/setup_breast_tcga_configs.R')\n")
cat("2. Run: source('scripts/run_all_pipelines.R')\n")
