# R/03_filtering.R
# Functions for filtering low-expression genes

#' Filter low-expression genes using edgeR's filterByExpr
#' @param counts Counts matrix
#' @param metadata Sample metadata
#' @param group_col Column in metadata for grouping
#' @param min_count Minimum count threshold
#' @param min_samples Minimum samples with min_count (NULL = automatic)
#' @return List with filtered counts and statistics
filter_low_expression <- function(counts,
                                   metadata,
                                   group_col = "condition",
                                   min_count = 10,
                                   min_samples = NULL) {

  n_genes_before <- nrow(counts)

  # Use edgeR's filterByExpr for robust filtering
  if (requireNamespace("edgeR", quietly = TRUE)) {

    # Create DGEList object
    dge <- edgeR::DGEList(counts = counts)

    # Get group factor
    if (group_col %in% colnames(metadata)) {
      group <- factor(metadata[[group_col]])
    } else {
      warning("Group column '", group_col, "' not found. Using all samples as one group.")
      group <- factor(rep("all", ncol(counts)))
    }

    # Determine min_samples if not specified
    if (is.null(min_samples)) {
      # Default: smallest group size
      min_samples <- min(table(group))
    }

    # Apply filterByExpr
    keep <- edgeR::filterByExpr(
      dge,
      group = group,
      min.count = min_count,
      min.total.count = min_count * 2,
      large.n = 10,
      min.prop = 0.7
    )

  } else {
    # Fallback: simple filtering
    warning("edgeR not available. Using simple filtering approach.")

    if (is.null(min_samples)) {
      min_samples <- floor(ncol(counts) * 0.2)  # At least 20% of samples
    }

    keep <- rowSums(counts >= min_count) >= min_samples
  }

  # Filter counts
  counts_filtered <- counts[keep, , drop = FALSE]
  n_genes_after <- nrow(counts_filtered)
  n_removed <- n_genes_before - n_genes_after

  # Compile statistics
  stats <- data.frame(
    metric = c(
      "genes_before_filtering",
      "genes_after_filtering",
      "genes_removed",
      "percent_retained",
      "min_count_threshold",
      "min_samples_threshold"
    ),
    value = c(
      n_genes_before,
      n_genes_after,
      n_removed,
      round(n_genes_after / n_genes_before * 100, 2),
      min_count,
      min_samples
    ),
    stringsAsFactors = FALSE
  )

  message("Filtering: ", n_genes_before, " -> ", n_genes_after, " genes ",
          "(removed ", n_removed, ", ", round(n_removed/n_genes_before*100, 1), "%)")

  list(
    counts = counts_filtered,
    stats = stats,
    keep = keep
  )
}
