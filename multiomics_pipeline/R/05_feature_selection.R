# =============================================================================
# Feature Selection for Integration
# =============================================================================

#' Select features for integration across all omics
select_features_for_integration <- function(mae_data, config) {
  log_message("=== Feature Selection for Integration ===")

  matrices <- extract_matrices_for_integration(mae_data)
  harmonized <- mae_data$harmonized_omics
  selected <- list()

  for (omic in names(matrices)) {
    mat <- matrices[[omic]]

    # Get omic-specific config
    if (omic == "transcriptomics") {
      fs_config <- config$feature_selection$transcriptomics
      da_table <- harmonized$transcriptomics$de_table
    } else if (omic == "proteomics") {
      fs_config <- config$feature_selection$proteomics
      da_table <- harmonized$proteomics$da_table
    } else if (omic == "metabolomics") {
      fs_config <- config$feature_selection$metabolomics
      da_table <- harmonized$metabolomics$da_table
    } else {
      fs_config <- list(strategy = "variance", top_n_variance = 1000)
      da_table <- NULL
    }

    strategy <- fs_config$strategy %||% "variance"
    top_n <- fs_config$top_n_variance %||% 1000
    fdr_thresh <- fs_config$fdr_threshold %||% 0.05

    selected[[omic]] <- select_features_single_omics(
      mat = mat,
      da_table = da_table,
      strategy = strategy,
      top_n_variance = top_n,
      fdr_threshold = fdr_thresh,
      omic_name = omic
    )
  }

  # Create filtered MAE
  filtered_mae <- filter_mae_features(mae_data, selected)

  # Save selection summary
  selection_summary <- summarize_feature_selection(selected, config)

  list(
    filtered_mae = filtered_mae,
    selected_features = selected,
    selection_summary = selection_summary,
    original_mae = mae_data
  )
}

#' Select features for a single omics dataset
select_features_single_omics <- function(mat, da_table, strategy, top_n_variance,
                                          fdr_threshold, omic_name) {
  log_message("Selecting features for ", omic_name, " using strategy: ", strategy)

  n_features <- nrow(mat)
  all_features <- rownames(mat)

  # Calculate variance for all features
  feature_var <- apply(mat, 1, var, na.rm = TRUE)
  feature_var[is.na(feature_var)] <- 0

  # Get significant features from DA table
  sig_features <- NULL
  if (!is.null(da_table) && nrow(da_table) > 0) {
    # Handle different column names for adjusted p-value
    padj_col <- intersect(c("adj.P.Val", "padj", "FDR", "q.value"),
                          colnames(da_table))[1]
    id_col <- intersect(c("feature_id", "gene_id", "protein_id"),
                        colnames(da_table))[1]

    if (!is.na(padj_col) && !is.na(id_col)) {
      sig_idx <- which(da_table[[padj_col]] < fdr_threshold)
      sig_features <- unique(da_table[[id_col]][sig_idx])
      sig_features <- intersect(sig_features, all_features)
    }
  }

  # Apply selection strategy
  if (strategy == "all") {
    selected <- all_features
    log_message("  Using all ", length(selected), " features")

  } else if (strategy == "variance") {
    # Top N by variance
    n_select <- min(top_n_variance, n_features)
    selected <- names(sort(feature_var, decreasing = TRUE))[1:n_select]
    log_message("  Selected top ", length(selected), " by variance")

  } else if (strategy == "significant") {
    # Only significant features
    if (is.null(sig_features) || length(sig_features) == 0) {
      log_message("  No significant features found, falling back to variance")
      n_select <- min(top_n_variance, n_features)
      selected <- names(sort(feature_var, decreasing = TRUE))[1:n_select]
    } else {
      selected <- sig_features
      log_message("  Selected ", length(selected), " significant features (FDR < ",
                 fdr_threshold, ")")
    }

  } else if (strategy == "hybrid") {
    # Union of significant + top variance
    n_variance <- min(top_n_variance, n_features)
    top_var_features <- names(sort(feature_var, decreasing = TRUE))[1:n_variance]

    if (!is.null(sig_features) && length(sig_features) > 0) {
      selected <- union(sig_features, top_var_features)
      log_message("  Hybrid: ", length(sig_features), " significant + ",
                 length(top_var_features), " top variance = ",
                 length(selected), " unique features")
    } else {
      selected <- top_var_features
      log_message("  Hybrid (no sig): using top ", length(selected), " by variance")
    }

  } else {
    # Default to variance
    n_select <- min(top_n_variance, n_features)
    selected <- names(sort(feature_var, decreasing = TRUE))[1:n_select]
    log_message("  Unknown strategy, using variance: ", length(selected), " features")
  }

  # Create selection result
  list(
    selected_features = selected,
    n_original = n_features,
    n_selected = length(selected),
    strategy = strategy,
    n_significant = length(sig_features),
    feature_variance = feature_var[selected]
  )
}

#' Filter MAE to selected features
filter_mae_features <- function(mae_data, selected_features) {
  mae <- mae_data$mae

  if (inherits(mae, "MultiAssayExperiment")) {
    # Filter each experiment
    for (omic in names(selected_features)) {
      if (omic %in% names(MultiAssayExperiment::experiments(mae))) {
        sel <- selected_features[[omic]]$selected_features
        mae <- mae[sel, , omic]
      }
    }
  } else {
    # List structure
    for (omic in names(selected_features)) {
      if (omic %in% names(mae$experiments)) {
        sel <- selected_features[[omic]]$selected_features
        mae$experiments[[omic]] <- mae$experiments[[omic]][sel, , drop = FALSE]
      }
    }
  }

  # Update mae_data
  mae_data$mae <- mae
  mae_data
}

#' Summarize feature selection
summarize_feature_selection <- function(selected_features, config) {
  summary_df <- data.frame(
    omics = character(),
    strategy = character(),
    n_original = integer(),
    n_selected = integer(),
    n_significant = integer(),
    pct_selected = numeric(),
    stringsAsFactors = FALSE
  )

  for (omic in names(selected_features)) {
    sel <- selected_features[[omic]]
    summary_df <- rbind(summary_df, data.frame(
      omics = omic,
      strategy = sel$strategy,
      n_original = sel$n_original,
      n_selected = sel$n_selected,
      n_significant = sel$n_significant %||% 0,
      pct_selected = round(100 * sel$n_selected / sel$n_original, 1),
      stringsAsFactors = FALSE
    ))
  }

  save_table(summary_df, "feature_selection_summary.csv", config)

  for (i in seq_len(nrow(summary_df))) {
    log_message(summary_df$omics[i], ": ", summary_df$n_original[i], " -> ",
               summary_df$n_selected[i], " features (",
               summary_df$pct_selected[i], "%)")
  }

  summary_df
}

#' Get feature overlap between omics (for gene-centric comparison)
get_feature_overlap <- function(mae_data, config) {
  gene_mapping <- mae_data$gene_mapping

  if (is.null(gene_mapping)) {
    log_message("No gene mapping available for feature overlap analysis")
    return(NULL)
  }

  # Split by omics
  by_omics <- split(gene_mapping, gene_mapping$omics)

  if (length(by_omics) < 2) {
    return(NULL)
  }

  # Get unique genes per omics
  genes_per_omics <- lapply(by_omics, function(x) unique(x$gene_symbol))

  # Calculate overlaps
  omics_names <- names(genes_per_omics)
  overlap_matrix <- matrix(
    0,
    nrow = length(omics_names),
    ncol = length(omics_names),
    dimnames = list(omics_names, omics_names)
  )

  for (i in seq_along(omics_names)) {
    for (j in seq_along(omics_names)) {
      overlap_matrix[i, j] <- length(intersect(
        genes_per_omics[[i]],
        genes_per_omics[[j]]
      ))
    }
  }

  # Log summary
  log_message("Gene overlap between omics:")
  for (i in 1:(length(omics_names) - 1)) {
    for (j in (i + 1):length(omics_names)) {
      log_message("  ", omics_names[i], " & ", omics_names[j], ": ",
                 overlap_matrix[i, j], " genes")
    }
  }

  # Save
  overlap_df <- as.data.frame(overlap_matrix)
  overlap_df <- tibble::rownames_to_column(overlap_df, "omics")
  save_table(overlap_df, "gene_overlap_matrix.csv", config)

  list(
    overlap_matrix = overlap_matrix,
    genes_per_omics = genes_per_omics
  )
}
