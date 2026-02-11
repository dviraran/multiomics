# =============================================================================
# 06_export.R — Final annotated results table
# =============================================================================

export_results <- function(de_results, normalized_data, ingested_data, config) {
  log_message("=== Exporting Final Results ===")

  de_table <- de_results$table
  mat_norm <- normalized_data$matrix

  # Per-sample normalized values as columns
  norm_df <- as.data.frame(t(mat_norm), check.names = FALSE)
  norm_df <- data.frame(
    feature_id = colnames(norm_df),
    t(as.matrix(norm_df)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  # Actually we want features as rows, samples as columns

  per_sample <- as.data.frame(mat_norm, check.names = FALSE)
  per_sample$feature_id <- rownames(per_sample)

  # Merge DE results with per-sample values
  final <- merge(de_table, per_sample, by = "feature_id", all.x = TRUE)

  # Sort by adjusted p-value
  final <- final[order(final$adj.P.Val, na.last = TRUE), ]

  # Save
  save_table(final, "de_results_annotated.csv", config, "tables")

  # Also save normalised matrix separately
  norm_out <- as.data.frame(mat_norm, check.names = FALSE)
  norm_out <- tibble::rownames_to_column(norm_out, "feature_id")
  save_table(norm_out, "normalized_matrix.csv", config, "tables")

  log_message("=== Export Complete ===")

  final
}
