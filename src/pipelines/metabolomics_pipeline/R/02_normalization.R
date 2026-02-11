# =============================================================================
# 02_normalization.R — MetaboAnalystR normalization + manual fallback
# =============================================================================

normalize_data <- function(ingested_data, config) {
  log_message("=== Starting Normalization ===")

  mat      <- ingested_data$matrix
  metadata <- ingested_data$metadata
  method   <- config$normalization$method %||% "metaboanalyst"

  mat_raw <- mat  # keep original for QC plots

  if (method == "metaboanalyst") {
    result <- tryCatch(
      normalize_metaboanalyst(mat, metadata, config),
      error = function(e) {
        log_message("MetaboAnalystR failed (", conditionMessage(e),
                    "). Falling back to manual normalization.")
        normalize_manual(mat, config)
      }
    )
  } else {
    result <- normalize_manual(mat, config)
  }

  mat_norm <- result$matrix

  # QC boxplot before/after
  plots <- normalization_qc_plots(mat_raw, mat_norm, metadata, config)

  # Save normalized matrix
  norm_df <- as.data.frame(mat_norm)
  norm_df <- tibble::rownames_to_column(norm_df, "feature_id")
  save_table(norm_df, "normalized_matrix.csv", config, "tables")

  log_message("=== Normalization Complete (method: ", result$method, ") ===")

  list(
    matrix      = mat_norm,
    matrix_raw  = mat_raw,
    metadata    = metadata,
    method_used = result$method,
    plots       = plots
  )
}

# ---------- MetaboAnalystR path ----------------------------------------------

normalize_metaboanalyst <- function(mat, metadata, config) {
  if (!requireNamespace("MetaboAnalystR", quietly = TRUE))
    stop("MetaboAnalystR is not installed")

  log_message("Running MetaboAnalystR normalization chain")

  sample_col    <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  # Build samples x features data.frame with Sample + Group columns
  # (MetaboAnalystR expects: rows = samples, cols = features, first two cols = Sample, Group)
  df_t <- as.data.frame(t(mat), check.names = FALSE)
  conditions <- metadata[[condition_col]][match(rownames(df_t),
                                                 metadata[[sample_col]])]
  df_t <- cbind(
    data.frame(Sample = rownames(df_t), Group = as.character(conditions),
               stringsAsFactors = FALSE),
    df_t
  )

  # Write to temp file
  tmp_dir <- tempfile("metabo_norm_")
  dir.create(tmp_dir, recursive = TRUE)
  combined_path <- file.path(tmp_dir, "combined_data.txt")
  write.table(df_t, combined_path, sep = "\t", row.names = FALSE, quote = FALSE)

  # MetaboAnalystR changes working directory internally — protect ours

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  mSet <- MetaboAnalystR::InitDataObjects("conc", "stat", FALSE)
  mSet <- MetaboAnalystR::Read.TextData(mSet, combined_path, "rowu", "disc")
  mSet <- MetaboAnalystR::SanityCheckData(mSet)
  mSet <- MetaboAnalystR::ReplaceMin(mSet)
  mSet <- MetaboAnalystR::PreparePrenormData(mSet)

  row_norm   <- config$normalization$row_norm   %||% "NULL"
  trans_norm <- config$normalization$trans_norm  %||% "LogNorm"
  scale_norm <- config$normalization$scale_norm  %||% "MeanCenter"

  mSet <- MetaboAnalystR::Normalization(
    mSet, row_norm, trans_norm, scale_norm,
    "S10T0", ratio = FALSE, ratioNum = 20
  )

  # Extract normalised data (samples x features) and transpose to features x samples
  norm_mat <- as.matrix(mSet$dataSet$norm)
  mat_out  <- t(norm_mat)

  # Restore original feature names — MetaboAnalystR strips special chars like |
  # Match by stripping the same chars from original names
  orig_ids <- rownames(mat)
  mangled  <- gsub("[^A-Za-z0-9]", "", orig_ids)
  mset_ids <- gsub("[^A-Za-z0-9]", "", rownames(mat_out))
  idx      <- match(mset_ids, mangled)
  if (!anyNA(idx)) {
    rownames(mat_out) <- orig_ids[idx]
  } else {
    log_message("Warning: could not fully restore original feature IDs (",
                sum(is.na(idx)), " unmatched)")
    matched <- !is.na(idx)
    rownames(mat_out)[matched] <- orig_ids[idx[matched]]
  }

  # Clean up temp files
  unlink(tmp_dir, recursive = TRUE)

  log_message("MetaboAnalystR normalization complete: ", nrow(mat_out),
              " features x ", ncol(mat_out), " samples")

  list(matrix = mat_out, method = "metaboanalyst")
}

# ---------- manual fallback --------------------------------------------------

normalize_manual <- function(mat, config) {
  log_message("Running manual normalization (log2 + mean-center)")

  pseudocount <- config$normalization$pseudocount %||% 1

  # log2 transform
  mat_log <- log2(mat + pseudocount)

  # Mean-center per feature (row-wise)
  row_means <- rowMeans(mat_log, na.rm = TRUE)
  mat_norm  <- mat_log - row_means

  log_message("Manual normalization complete")
  list(matrix = mat_norm, method = "manual")
}

# ---------- QC plots ---------------------------------------------------------

normalization_qc_plots <- function(mat_raw, mat_norm, metadata, config) {
  plots <- list()

  # Boxplot before/after
  n_samples <- ncol(mat_raw)
  sample_names <- colnames(mat_raw)
  sample_col    <- config$input$sample_id_column
  condition_col <- config$design$condition_column
  cond_vec <- as.character(
    metadata[[condition_col]][match(sample_names, metadata[[sample_col]])]
  )

  before_df <- data.frame(
    intensity = as.vector(log2(mat_raw + 1)),
    sample    = rep(sample_names, each = nrow(mat_raw)),
    stage     = "Before normalization",
    condition = rep(cond_vec, each = nrow(mat_raw)),
    stringsAsFactors = FALSE
  )
  after_df <- data.frame(
    intensity = as.vector(mat_norm),
    sample    = rep(colnames(mat_norm), each = nrow(mat_norm)),
    stage     = "After normalization",
    condition = rep(cond_vec, each = nrow(mat_norm)),
    stringsAsFactors = FALSE
  )
  box_df <- rbind(before_df, after_df)
  box_df <- box_df[!is.na(box_df$intensity), ]

  plots$boxplot <- ggplot2::ggplot(
    box_df,
    ggplot2::aes(x = sample, y = intensity, fill = condition)
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.3) +
    ggplot2::facet_wrap(~stage, ncol = 1, scales = "free_y") +
    ggplot2::labs(title = "Intensity Distribution Before/After Normalization",
                  x = "Sample", y = "Intensity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7),
      legend.position = "bottom"
    )

  save_plot(plots$boxplot, "normalization_boxplot.png", config,
            width = 12, height = 10, subdir = "qc")

  plots
}
