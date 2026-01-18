# =============================================================================
# Per-Omics Preprocessing
# =============================================================================

#' Preprocess all omics data
preprocess_all_omics <- function(ingested_data, config) {
  log_message("=== Starting Per-Omics Preprocessing ===")

  metadata <- ingested_data$metadata
  omics_data <- ingested_data$omics_data
  processed <- list()

  if ("transcriptomics" %in% names(omics_data)) {
    processed$transcriptomics <- preprocess_transcriptomics(
      omics_data$transcriptomics, metadata, config
    )
  }

  if ("proteomics" %in% names(omics_data)) {
    processed$proteomics <- preprocess_proteomics(
      omics_data$proteomics, metadata, config
    )
  }

  if ("metabolomics" %in% names(omics_data)) {
    processed$metabolomics <- preprocess_metabolomics(
      omics_data$metabolomics, metadata, config
    )
  }

  log_message("=== Preprocessing Complete ===")

  list(
    processed_omics = processed,
    metadata = metadata,
    alignment = ingested_data$alignment
  )
}

#' Preprocess transcriptomics
preprocess_transcriptomics <- function(rna_data, metadata, config) {
  log_message("Preprocessing transcriptomics...")
  tc <- config$transcriptomics

  if (rna_data$mode == "preprocessed") {
    log_message("Using preprocessed RNA data")
    return(list(
      normalized_matrix = rna_data$matrix,
      de_table = rna_data$de_table,
      mapping = rna_data$mapping,
      gmt = rna_data$gmt
    ))
  }

  # Raw mode processing
  counts <- rna_data$matrix
  sample_col <- config$global$sample_id_column
  condition_col <- config$design$condition_column

  # Align samples
  common_samples <- intersect(colnames(counts), metadata[[sample_col]])
  counts <- counts[, common_samples]
  meta_aligned <- metadata[match(common_samples, metadata[[sample_col]]), ]

  # Filter low expression
  min_counts <- tc$processing$min_counts %||% 10
  min_samples <- tc$processing$min_samples %||% 3
  keep <- rowSums(counts >= min_counts) >= min_samples
  counts_filtered <- counts[keep, ]
  log_message("Filtered RNA: ", nrow(counts), " -> ", nrow(counts_filtered), " genes")

  # Normalization and DE
  norm_method <- tc$processing$normalization %||% "vst"
  de_method <- tc$processing$de_method %||% "DESeq2"

  if (de_method == "DESeq2" && requireNamespace("DESeq2", quietly = TRUE)) {
    result <- run_deseq2(counts_filtered, meta_aligned, config)
    normalized_matrix <- result$normalized
    de_table <- result$de_table
  } else if (requireNamespace("limma", quietly = TRUE) && requireNamespace("edgeR", quietly = TRUE)) {
    result <- run_limma_voom(counts_filtered, meta_aligned, config)
    normalized_matrix <- result$normalized
    de_table <- result$de_table
  } else {
    log_message("No DE package available. Simple log2(CPM+1) normalization")
    lib_sizes <- colSums(counts_filtered)
    cpm <- sweep(counts_filtered, 2, lib_sizes, "/") * 1e6
    normalized_matrix <- log2(cpm + 1)
    de_table <- NULL
  }

  # Save results
  norm_df <- as.data.frame(normalized_matrix)
  norm_df <- tibble::rownames_to_column(norm_df, "gene_id")
  save_table(norm_df, "rna_normalized_matrix.csv", config)

  if (!is.null(de_table)) {
    save_table(de_table, "rna_de_results.csv", config)
  }

  list(
    normalized_matrix = normalized_matrix,
    de_table = de_table,
    mapping = rna_data$mapping,
    gmt = rna_data$gmt
  )
}

#' Run DESeq2 analysis
run_deseq2 <- function(counts, metadata, config) {
  condition_col <- config$design$condition_column
  formula <- as.formula(config$design$design_formula)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = metadata,
    design = formula
  )

  dds <- DESeq2::DESeq(dds)

  if (nrow(dds) < 1000) {
    log_message("Warning: Small dataset (< 1000 features). Using rlog instead of vst.")
    normalized <- DESeq2::rlog(dds, blind = FALSE)
  } else {
    normalized <- DESeq2::vst(dds, blind = FALSE)
  }

  normalized_matrix <- SummarizedExperiment::assay(normalized)

  # Get DE results for each contrast
  contrasts <- config$design$contrasts
  all_results <- list()

  for (contrast_str in contrasts) {
    parts <- strsplit(trimws(contrast_str), "\\s*-\\s*")[[1]]
    if (length(parts) == 2) {
      res <- DESeq2::results(dds, contrast = c(condition_col, parts[1], parts[2]))
      res_df <- as.data.frame(res)
      res_df$gene_id <- rownames(res_df)
      res_df$contrast <- contrast_str
      all_results[[contrast_str]] <- res_df
    }
  }

  de_table <- do.call(rbind, all_results)
  de_table <- de_table[, c("gene_id", "log2FoldChange", "baseMean", "pvalue", "padj", "contrast")]
  colnames(de_table)[2:5] <- c("log2FC", "baseMean", "pvalue", "adj.P.Val")

  list(normalized = normalized_matrix, de_table = de_table)
}

#' Run limma-voom analysis
run_limma_voom <- function(counts, metadata, config) {
  condition_col <- config$design$condition_column
  formula <- as.formula(config$design$design_formula)

  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)

  design <- model.matrix(formula, data = metadata)
  v <- limma::voom(dge, design)

  fit <- limma::lmFit(v, design)

  # Build contrasts
  contrasts_spec <- config$design$contrasts
  all_results <- list()

  for (contrast_str in contrasts_spec) {
    parts <- strsplit(trimws(contrast_str), "\\s*-\\s*")[[1]]
    if (length(parts) == 2) {
      coef_name <- paste0(condition_col, parts[1])
      if (coef_name %in% colnames(design)) {
        fit2 <- limma::eBayes(fit)
        tt <- limma::topTable(fit2, coef = coef_name, number = Inf)
        tt$gene_id <- rownames(tt)
        tt$contrast <- contrast_str
        all_results[[contrast_str]] <- tt
      }
    }
  }

  de_table <- do.call(rbind, all_results)
  de_table <- de_table[, c("gene_id", "logFC", "AveExpr", "P.Value", "adj.P.Val", "contrast")]
  colnames(de_table)[2:4] <- c("log2FC", "baseMean", "pvalue")

  list(normalized = v$E, de_table = de_table)
}

#' Preprocess proteomics
preprocess_proteomics <- function(prot_data, metadata, config) {
  log_message("Preprocessing proteomics...")
  pc <- config$proteomics

  if (prot_data$mode == "preprocessed") {
    log_message("Using preprocessed proteomics data")
    return(list(
      normalized_matrix = prot_data$matrix,
      da_table = prot_data$da_table,
      mapping = prot_data$mapping
    ))
  }

  mat <- prot_data$matrix
  sample_col <- config$global$sample_id_column

  # Align samples
  common_samples <- intersect(colnames(mat), metadata[[sample_col]])
  mat <- mat[, common_samples]
  meta_aligned <- metadata[match(common_samples, metadata[[sample_col]]), ]

  # Filter by presence
  min_presence <- pc$processing$min_presence %||% 0.5
  keep <- rowMeans(!is.na(mat)) >= min_presence
  mat <- mat[keep, ]
  log_message("Filtered proteins: ", sum(!keep), " removed, ", nrow(mat), " remaining")

  # Log transform
  if (pc$processing$log_transform %||% TRUE) {
    max_val <- max(mat, na.rm = TRUE)
    if (max_val > 50) {
      mat <- log2(mat + 1)
    }
  }

  # Normalize
  norm_method <- pc$processing$normalization %||% "median"
  if (norm_method == "median") {
    sample_medians <- apply(mat, 2, median, na.rm = TRUE)
    global_median <- median(sample_medians)
    mat <- sweep(mat, 2, sample_medians, "-") + global_median
  } else if (norm_method == "quantile" && requireNamespace("limma", quietly = TRUE)) {
    mat <- limma::normalizeQuantiles(mat)
  }

  # Imputation
  impute_method <- pc$processing$imputation %||% "half_min"
  if (impute_method == "half_min") {
    for (i in seq_len(nrow(mat))) {
      row <- mat[i, ]
      na_idx <- is.na(row)
      if (any(na_idx) && any(!na_idx)) {
        mat[i, na_idx] <- min(row, na.rm = TRUE) * 0.5
      }
    }
  }

  # DA analysis with limma
  da_table <- NULL
  if (requireNamespace("limma", quietly = TRUE)) {
    da_table <- run_limma_da(mat, meta_aligned, config, "protein")
  }

  # Save
  norm_df <- as.data.frame(mat)
  norm_df <- tibble::rownames_to_column(norm_df, "protein_id")
  save_table(norm_df, "prot_normalized_matrix.csv", config)
  if (!is.null(da_table)) {
    save_table(da_table, "prot_da_results.csv", config)
  }

  list(
    normalized_matrix = mat,
    da_table = da_table,
    mapping = prot_data$mapping
  )
}

#' Preprocess metabolomics
preprocess_metabolomics <- function(metab_data, metadata, config) {
  log_message("Preprocessing metabolomics...")
  mc <- config$metabolomics

  if (metab_data$mode == "preprocessed") {
    log_message("Using preprocessed metabolomics data")
    return(list(
      normalized_matrix = metab_data$matrix,
      da_table = metab_data$da_table,
      feature_metadata = metab_data$feature_metadata,
      annotation = metab_data$annotation,
      pathway_mapping = metab_data$pathway_mapping,
      gmt = metab_data$gmt
    ))
  }

  mat <- metab_data$matrix
  sample_col <- config$global$sample_id_column

  # Align samples
  common_samples <- intersect(colnames(mat), metadata[[sample_col]])
  mat <- mat[, common_samples]
  meta_aligned <- metadata[match(common_samples, metadata[[sample_col]]), ]

  # Filter
  min_presence <- mc$processing$min_presence %||% 0.3
  keep <- rowMeans(!is.na(mat)) >= min_presence
  mat <- mat[keep, ]
  log_message("Filtered metabolites: ", nrow(mat), " remaining")

  # Transform
  transform <- mc$processing$transform %||% "log2"
  if (transform == "log2") {
    max_val <- max(mat, na.rm = TRUE)
    if (max_val > 50) {
      mat <- log2(mat + 1)
    }
  }

  # Normalize (PQN)
  norm_method <- mc$processing$normalization %||% "PQN"
  if (norm_method == "PQN") {
    reference <- apply(mat, 1, median, na.rm = TRUE)
    for (j in seq_len(ncol(mat))) {
      quotients <- mat[, j] / reference
      quotients <- quotients[is.finite(quotients) & quotients > 0]
      if (length(quotients) > 0) {
        mat[, j] <- mat[, j] / median(quotients, na.rm = TRUE)
      }
    }
  } else if (norm_method == "median") {
    sample_medians <- apply(mat, 2, median, na.rm = TRUE)
    mat <- sweep(mat, 2, sample_medians, "-") + median(sample_medians)
  }

  # Imputation
  impute_method <- mc$processing$imputation %||% "half_min"
  if (impute_method == "half_min") {
    for (i in seq_len(nrow(mat))) {
      row <- mat[i, ]
      na_idx <- is.na(row)
      if (any(na_idx) && any(!na_idx)) {
        mat[i, na_idx] <- min(row, na.rm = TRUE) * 0.5
      }
    }
  }

  # DA analysis
  da_table <- NULL
  if (requireNamespace("limma", quietly = TRUE)) {
    da_table <- run_limma_da(mat, meta_aligned, config, "metabolite")
  }

  # Save
  norm_df <- as.data.frame(mat)
  norm_df <- tibble::rownames_to_column(norm_df, "feature_id")
  save_table(norm_df, "metab_normalized_matrix.csv", config)
  if (!is.null(da_table)) {
    save_table(da_table, "metab_da_results.csv", config)
  }

  list(
    normalized_matrix = mat,
    da_table = da_table,
    feature_metadata = metab_data$feature_metadata,
    annotation = metab_data$annotation,
    pathway_mapping = metab_data$pathway_mapping,
    gmt = metab_data$gmt
  )
}

#' Run limma DA
run_limma_da <- function(mat, metadata, config, feature_type = "feature") {
  condition_col <- config$design$condition_column
  formula <- as.formula(config$design$design_formula)

  design <- model.matrix(formula, data = metadata)
  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit)

  contrasts_spec <- config$design$contrasts
  all_results <- list()

  for (contrast_str in contrasts_spec) {
    parts <- strsplit(trimws(contrast_str), "\\s*-\\s*")[[1]]
    if (length(parts) == 2) {
      coef_name <- paste0(condition_col, parts[1])
      if (coef_name %in% colnames(design)) {
        tt <- limma::topTable(fit, coef = coef_name, number = Inf)
        tt$feature_id <- rownames(tt)
        tt$contrast <- contrast_str
        all_results[[contrast_str]] <- tt
      }
    }
  }

  if (length(all_results) == 0) {
    # Fallback: use second coefficient
    if (ncol(design) >= 2) {
      tt <- limma::topTable(fit, coef = 2, number = Inf)
      tt$feature_id <- rownames(tt)
      tt$contrast <- "coef2"
      all_results[["coef2"]] <- tt
    }
  }

  if (length(all_results) == 0) {
    return(NULL)
  }

  da_table <- do.call(rbind, all_results)
  da_table <- da_table[, c("feature_id", "logFC", "AveExpr", "P.Value", "adj.P.Val", "contrast")]
  colnames(da_table)[2:4] <- c("log2FC", "avgExpr", "pvalue")

  return(da_table)
}
