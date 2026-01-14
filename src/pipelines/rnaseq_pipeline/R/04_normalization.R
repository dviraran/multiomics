# R/04_normalization.R
# Functions for DESeq2 normalization and transformation

#' Create DESeq2 object from counts and metadata
#' @param counts Filtered counts matrix
#' @param metadata Sample metadata
#' @param design_formula Design formula for DESeq2
#' @return DESeqDataSet object
create_deseq_object <- function(counts, metadata, design_formula) {

  # Ensure counts are integers
  storage.mode(counts) <- "integer"

  # Ensure metadata row order matches counts column order
  metadata <- metadata[match(colnames(counts), rownames(metadata)), , drop = FALSE]

  # Convert character columns to factors for DESeq2
  for (col in colnames(metadata)) {
    if (is.character(metadata[[col]])) {
      metadata[[col]] <- factor(metadata[[col]])
    }
  }

  # Create DESeqDataSet
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = design_formula
  )

  message("Created DESeq2 object: ", nrow(dds), " genes x ", ncol(dds), " samples")
  message("Design formula: ", as.character(design_formula)[2])

  dds
}


#' Normalize DESeq2 object (estimate size factors and dispersions)
#' @param dds DESeqDataSet object
#' @return Normalized DESeqDataSet
normalize_deseq <- function(dds) {

  message("Running DESeq2 normalization...")

  # Estimate size factors
  dds <- DESeq2::estimateSizeFactors(dds)

  # Report size factors
  sf <- DESeq2::sizeFactors(dds)
  message("Size factors range: ", round(min(sf), 3), " - ", round(max(sf), 3))

  # Estimate dispersions
  dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)

  message("Normalization complete")

  dds
}


#' Get normalized counts from DESeq2 object
#' @param dds Normalized DESeqDataSet
#' @return Matrix of normalized counts
get_normalized_counts <- function(dds) {
  counts <- DESeq2::counts(dds, normalized = TRUE)
  message("Extracted normalized counts: ", nrow(counts), " x ", ncol(counts))
  counts
}


#' Get variance-stabilized transformed counts
#' @param dds Normalized DESeqDataSet
#' @param blind Whether to blind transformation to experimental design
#' @return Matrix of VST-transformed values
get_vst_counts <- function(dds, blind = TRUE) {
  n_genes <- nrow(dds)
  n_samples <- ncol(dds)

  # Use VST for datasets with many samples, rlog for small datasets
  if (n_samples < 30) {
    message("Using rlog transformation (n < 30 samples)")
    vst <- DESeq2::rlog(dds, blind = blind)
  } else {
    message("Using VST transformation (n >= 30 samples)")
    # Adjust nsub for small gene sets (default is 1000)
    nsub <- min(1000, n_genes)
    vst <- DESeq2::vst(dds, blind = blind, nsub = nsub)
  }

  vst_counts <- SummarizedExperiment::assay(vst)
  message("Generated variance-stabilized counts: ", nrow(vst_counts), " x ", ncol(vst_counts))

  vst_counts
}
