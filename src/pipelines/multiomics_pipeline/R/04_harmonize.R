# =============================================================================
# Sample and Feature Harmonization + MultiAssayExperiment Creation
# =============================================================================

#' Harmonize samples across omics and create MultiAssayExperiment
create_multiassay_experiment <- function(harmonized_data, config) {
  log_message("=== Creating MultiAssayExperiment ===")

  harmonized <- harmonized_data$harmonized_omics
  metadata <- harmonized_data$metadata
  sample_col <- config$global$sample_id_column
  sample_mode <- config$harmonization$sample_mode %||% "intersection"

  # Collect sample IDs from each omics
  omics_samples <- list()
  for (omic in names(harmonized)) {
    mat <- harmonized[[omic]]$normalized_matrix
    omics_samples[[omic]] <- colnames(mat)
  }

  # Determine common/union samples
  if (sample_mode == "intersection") {
    common_samples <- Reduce(intersect, omics_samples)
    log_message("Using intersection mode: ", length(common_samples), " common samples")
  } else {
    common_samples <- Reduce(union, omics_samples)
    log_message("Using union mode: ", length(common_samples), " total samples")
  }

  if (length(common_samples) < 3) {
    stop("Fewer than 3 samples in common across omics. Check sample IDs.")
  }

  # Align metadata
  meta_aligned <- metadata[metadata[[sample_col]] %in% common_samples, ]
  rownames(meta_aligned) <- meta_aligned[[sample_col]]
  meta_aligned <- meta_aligned[common_samples, ]

  # Create experiment list
  exp_list <- list()
  sample_map_list <- list()

  for (omic in names(harmonized))
{
    mat <- harmonized[[omic]]$normalized_matrix

    # Subset to common samples (fill with NA for union mode)
    if (sample_mode == "intersection") {
      mat <- mat[, common_samples, drop = FALSE]
    } else {
      # Union mode: create full matrix with NA for missing
      full_mat <- matrix(
        NA,
        nrow = nrow(mat),
        ncol = length(common_samples),
        dimnames = list(rownames(mat), common_samples)
      )
      overlap <- intersect(colnames(mat), common_samples)
      full_mat[, overlap] <- mat[, overlap]
      mat <- full_mat
    }

    # Create SummarizedExperiment
    if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      # Get feature annotation
      feat_anno <- harmonized[[omic]]$feature_annotation
      if (!is.null(feat_anno)) {
        rownames(feat_anno) <- feat_anno$feature_id
        feat_anno <- feat_anno[rownames(mat), , drop = FALSE]
        row_data <- S4Vectors::DataFrame(feat_anno)
      } else {
        row_data <- S4Vectors::DataFrame(feature_id = rownames(mat))
      }

      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(normalized = mat),
        rowData = row_data
      )
      exp_list[[omic]] <- se
    } else {
      # Fallback: store as matrix
      exp_list[[omic]] <- mat
    }

    # Create sample map for MAE
    sample_map_list[[omic]] <- data.frame(
      primary = colnames(mat),
      colname = colnames(mat),
      stringsAsFactors = FALSE
    )
  }

  # Create MultiAssayExperiment
  mae <- NULL
  if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    # Build sample map
    sample_map <- MultiAssayExperiment::listToMap(sample_map_list)

    # Create col_data
    col_data <- S4Vectors::DataFrame(meta_aligned)
    rownames(col_data) <- common_samples

    mae <- MultiAssayExperiment::MultiAssayExperiment(
      experiments = exp_list,
      colData = col_data,
      sampleMap = sample_map
    )

    log_message("Created MultiAssayExperiment with ", length(exp_list), " omics")
    for (omic in names(exp_list)) {
      if (inherits(exp_list[[omic]], "SummarizedExperiment")) {
        log_message("  - ", omic, ": ", nrow(exp_list[[omic]]), " features x ",
                   ncol(exp_list[[omic]]), " samples")
      } else {
        log_message("  - ", omic, ": ", nrow(exp_list[[omic]]), " features x ",
                   ncol(exp_list[[omic]]), " samples")
      }
    }
  } else {
    log_message("MultiAssayExperiment not available, using list structure")
    mae <- list(
      experiments = exp_list,
      colData = meta_aligned,
      sample_map = sample_map_list
    )
  }

  # Save sample alignment info
  alignment_df <- data.frame(
    sample_id = common_samples,
    stringsAsFactors = FALSE
  )
  for (omic in names(omics_samples)) {
    alignment_df[[paste0("in_", omic)]] <- common_samples %in% omics_samples[[omic]]
  }
  save_table(alignment_df, "sample_alignment.csv", config)

  list(
    mae = mae,
    harmonized_omics = harmonized,
    metadata = meta_aligned,
    common_samples = common_samples,
    gene_mapping = create_gene_mapping(harmonized_data)
  )
}

#' Handle duplicate features after mapping
resolve_duplicates <- function(mat, feature_ids, strategy = "keep_max_mean") {
  if (length(unique(feature_ids)) == length(feature_ids)) {
    rownames(mat) <- feature_ids
    return(mat)
  }

  log_message("Resolving ", length(feature_ids) - length(unique(feature_ids)),
             " duplicate feature IDs using strategy: ", strategy)

  unique_ids <- unique(feature_ids)
  result <- matrix(
    NA,
    nrow = length(unique_ids),
    ncol = ncol(mat),
    dimnames = list(unique_ids, colnames(mat))
  )

  for (uid in unique_ids) {
    idx <- which(feature_ids == uid)

    if (length(idx) == 1) {
      result[uid, ] <- mat[idx, ]
    } else {
      sub_mat <- mat[idx, , drop = FALSE]

      if (strategy == "sum") {
        result[uid, ] <- colSums(sub_mat, na.rm = TRUE)
      } else if (strategy == "median") {
        result[uid, ] <- apply(sub_mat, 2, median, na.rm = TRUE)
      } else if (strategy == "keep_first") {
        result[uid, ] <- sub_mat[1, ]
      } else if (strategy == "keep_max_mean") {
        row_means <- rowMeans(sub_mat, na.rm = TRUE)
        best_idx <- which.max(row_means)
        result[uid, ] <- sub_mat[best_idx, ]
      }
    }
  }

  result
}

#' Extract matrices from MAE for integration methods
extract_matrices_for_integration <- function(mae_data, omics_to_use = NULL) {
  mae <- mae_data$mae

  if (inherits(mae, "MultiAssayExperiment")) {
    exp_names <- names(MultiAssayExperiment::experiments(mae))

    if (!is.null(omics_to_use)) {
      exp_names <- intersect(exp_names, omics_to_use)
    }

    matrices <- list()
    for (exp_name in exp_names) {
      exp <- MultiAssayExperiment::experiments(mae)[[exp_name]]
      if (inherits(exp, "SummarizedExperiment")) {
        matrices[[exp_name]] <- SummarizedExperiment::assay(exp, "normalized")
      } else {
        matrices[[exp_name]] <- exp
      }
    }

    return(matrices)
  } else {
    # List structure fallback
    return(mae$experiments)
  }
}

#' Get sample metadata from MAE
get_sample_metadata <- function(mae_data) {
  mae <- mae_data$mae

  if (inherits(mae, "MultiAssayExperiment")) {
    return(as.data.frame(MultiAssayExperiment::colData(mae)))
  } else {
    return(mae$colData)
  }
}

#' Summarize MAE structure
summarize_mae <- function(mae_data, config) {
  log_message("=== MultiAssayExperiment Summary ===")

  mae <- mae_data$mae
  summary_df <- data.frame(
    omics = character(),
    n_features = integer(),
    n_samples = integer(),
    missing_pct = numeric(),
    stringsAsFactors = FALSE
  )

  matrices <- extract_matrices_for_integration(mae_data)

  for (omic in names(matrices)) {
    mat <- matrices[[omic]]
    missing_pct <- 100 * sum(is.na(mat)) / length(mat)

    summary_df <- rbind(summary_df, data.frame(
      omics = omic,
      n_features = nrow(mat),
      n_samples = ncol(mat),
      missing_pct = round(missing_pct, 2),
      stringsAsFactors = FALSE
    ))

    log_message(omic, ": ", nrow(mat), " features x ", ncol(mat), " samples (",
               round(missing_pct, 1), "% missing)")
  }

  save_table(summary_df, "mae_summary.csv", config)

  summary_df
}
