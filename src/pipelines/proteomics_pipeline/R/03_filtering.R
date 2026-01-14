# =============================================================================
# Contaminant Removal and Feature Filtering
# =============================================================================

#' Remove contaminants and apply feature filters
#'
#' @param ingested_data Data from ingestion step
#' @param annotation_data Annotation data
#' @param config Configuration list
#' @return Filtered data list
filter_features <- function(ingested_data, annotation_data, config) {
  log_message("=== Starting Feature Filtering ===")

  mat <- ingested_data$matrix
  metadata <- ingested_data$metadata
  annotations <- annotation_data$annotations

  # Track filtering steps
  filter_log <- data.frame(
    step = character(),
    features_before = integer(),
    features_removed = integer(),
    features_after = integer(),
    stringsAsFactors = FALSE
  )

  initial_n <- nrow(mat)

  # Step 1: Remove contaminants
  result <- remove_contaminants(mat, annotations, config)
  mat <- result$matrix
  annotations <- annotations[annotations$feature_id %in% rownames(mat), ]

  filter_log <- rbind(filter_log, data.frame(
    step = "Contaminant removal",
    features_before = initial_n,
    features_removed = result$n_removed,
    features_after = nrow(mat),
    stringsAsFactors = FALSE
  ))

  # Step 2: Global missingness filter
  before_n <- nrow(mat)
  result <- filter_by_global_presence(mat, config)
  mat <- result$matrix

  filter_log <- rbind(filter_log, data.frame(
    step = "Global presence filter",
    features_before = before_n,
    features_removed = result$n_removed,
    features_after = nrow(mat),
    stringsAsFactors = FALSE
  ))

  # Step 3: Group-aware missingness filter
  before_n <- nrow(mat)
  result <- filter_by_group_presence(mat, metadata, config)
  mat <- result$matrix

  filter_log <- rbind(filter_log, data.frame(
    step = "Group presence filter",
    features_before = before_n,
    features_removed = result$n_removed,
    features_after = nrow(mat),
    stringsAsFactors = FALSE
  ))

  # Step 4: Low intensity filter (optional)
  if (config$filtering$low_intensity_quantile > 0) {
    before_n <- nrow(mat)
    result <- filter_low_intensity(mat, config)
    mat <- result$matrix

    filter_log <- rbind(filter_log, data.frame(
      step = "Low intensity filter",
      features_before = before_n,
      features_removed = result$n_removed,
      features_after = nrow(mat),
      stringsAsFactors = FALSE
    ))
  }

  # Update annotations to match filtered matrix
  annotations <- annotations[annotations$feature_id %in% rownames(mat), ]

  # Save filter log
  save_table(filter_log, "filtering_summary.csv", config, "qc")

  log_message("Filtering complete: ", initial_n, " -> ", nrow(mat), " features")
  log_message("=== Feature Filtering Complete ===")

  list(
    matrix = mat,
    metadata = metadata,
    annotations = annotations,
    filter_log = filter_log
  )
}

#' Remove contaminant proteins
#'
#' @param mat Intensity matrix
#' @param annotations Annotation data frame
#' @param config Configuration list
#' @return List with filtered matrix and removal count
remove_contaminants <- function(mat, annotations, config) {
  feature_ids <- rownames(mat)
  n_initial <- length(feature_ids)

  # Get contaminant patterns from config
  patterns <- config$filtering$contaminant_patterns %||%
    c("CON__", "REV__", "^DECOY_", "HUMAN_CONTAM")

  # Load custom contaminants file if provided
  custom_contaminants <- NULL
  if (!is.null(config$optional_inputs$contaminants_file) &&
      file.exists(config$optional_inputs$contaminants_file)) {
    custom_contaminants <- readLines(config$optional_inputs$contaminants_file)
    custom_contaminants <- custom_contaminants[nchar(custom_contaminants) > 0]
    log_message("Loaded ", length(custom_contaminants), " custom contaminant IDs")
  }

  # Identify contaminants by pattern matching
  is_contaminant <- rep(FALSE, n_initial)

  for (pattern in patterns) {
    matches <- grepl(pattern, feature_ids, ignore.case = TRUE)
    is_contaminant <- is_contaminant | matches
  }

  # Add common keratin proteins
  keratin_patterns <- c("^KRT[0-9]", "^KRT1[0-9]", "KERATIN", "^K2C[0-9]", "^K1C[0-9]")
  for (pattern in keratin_patterns) {
    # Check in gene symbols if available
    if ("gene_symbol" %in% colnames(annotations)) {
      matches <- grepl(pattern, annotations$gene_symbol, ignore.case = TRUE)
      is_contaminant <- is_contaminant | matches
    }
    # Also check in feature IDs
    matches <- grepl(pattern, feature_ids, ignore.case = TRUE)
    is_contaminant <- is_contaminant | matches
  }

  # Check against custom contaminant list
  if (!is.null(custom_contaminants)) {
    # Match against feature IDs
    matches <- feature_ids %in% custom_contaminants
    is_contaminant <- is_contaminant | matches

    # Match against UniProt accessions
    if ("uniprot_accession" %in% colnames(annotations)) {
      matches <- annotations$uniprot_accession %in% custom_contaminants
      is_contaminant <- is_contaminant | matches
    }
  }

  n_removed <- sum(is_contaminant)
  mat_filtered <- mat[!is_contaminant, , drop = FALSE]

  if (n_removed > 0) {
    removed_ids <- feature_ids[is_contaminant]
    log_message("Removed ", n_removed, " contaminants")

    # Save list of removed contaminants
    contam_df <- data.frame(
      feature_id = removed_ids,
      stringsAsFactors = FALSE
    )
    if ("gene_symbol" %in% colnames(annotations)) {
      contam_df$gene_symbol <- annotations$gene_symbol[match(removed_ids, annotations$feature_id)]
    }
    # Note: We'll save this in the parent function
  }

  list(
    matrix = mat_filtered,
    n_removed = n_removed,
    removed_ids = if (n_removed > 0) removed_ids else character()
  )
}

#' Filter by global presence (present in X% of all samples)
#'
#' @param mat Intensity matrix
#' @param config Configuration list
#' @return List with filtered matrix and removal count
filter_by_global_presence <- function(mat, config) {
  threshold <- config$filtering$global_min_presence

  if (threshold <= 0) {
    return(list(matrix = mat, n_removed = 0))
  }

  n_samples <- ncol(mat)
  min_present <- ceiling(n_samples * threshold)

  n_present <- rowSums(!is.na(mat))
  keep <- n_present >= min_present

  n_removed <- sum(!keep)
  mat_filtered <- mat[keep, , drop = FALSE]

  log_message("Global filter (>=", round(threshold * 100), "% presence): ",
              "removed ", n_removed, " features")

  list(
    matrix = mat_filtered,
    n_removed = n_removed
  )
}

#' Filter by group-specific presence
#'
#' Features must be present in X% of samples in at least one condition group
#'
#' @param mat Intensity matrix
#' @param metadata Sample metadata
#' @param config Configuration list
#' @return List with filtered matrix and removal count
filter_by_group_presence <- function(mat, metadata, config) {
  threshold <- config$filtering$group_min_presence

  if (threshold <= 0) {
    return(list(matrix = mat, n_removed = 0))
  }

  condition_col <- config$design$condition_column
  sample_col <- config$input$sample_id_column
  conditions <- unique(metadata[[condition_col]])

  # For each feature, check if it passes threshold in at least one group
  passes_filter <- rep(FALSE, nrow(mat))

  for (cond in conditions) {
    samples <- metadata[[sample_col]][metadata[[condition_col]] == cond]
    samples <- intersect(samples, colnames(mat))

    if (length(samples) == 0) next

    sub_mat <- mat[, samples, drop = FALSE]
    n_samples <- ncol(sub_mat)
    min_present <- ceiling(n_samples * threshold)

    n_present <- rowSums(!is.na(sub_mat))
    passes_filter <- passes_filter | (n_present >= min_present)
  }

  n_removed <- sum(!passes_filter)
  mat_filtered <- mat[passes_filter, , drop = FALSE]

  log_message("Group filter (>=", round(threshold * 100), "% in any group): ",
              "removed ", n_removed, " features")

  list(
    matrix = mat_filtered,
    n_removed = n_removed
  )
}

#' Filter low intensity features
#'
#' @param mat Intensity matrix
#' @param config Configuration list
#' @return List with filtered matrix and removal count
filter_low_intensity <- function(mat, config) {
  quantile_threshold <- config$filtering$low_intensity_quantile

  if (quantile_threshold <= 0) {
    return(list(matrix = mat, n_removed = 0))
  }

  # Calculate mean intensity per feature (ignoring NAs)
  mean_intensities <- rowMeans(mat, na.rm = TRUE)

  # Find threshold
  intensity_cutoff <- quantile(mean_intensities, probs = quantile_threshold, na.rm = TRUE)

  keep <- mean_intensities >= intensity_cutoff
  n_removed <- sum(!keep)
  mat_filtered <- mat[keep, , drop = FALSE]

  log_message("Low intensity filter (bottom ", round(quantile_threshold * 100), "%): ",
              "removed ", n_removed, " features")

  list(
    matrix = mat_filtered,
    n_removed = n_removed
  )
}

#' Generate filtering diagnostic plots
#'
#' @param mat_before Matrix before filtering
#' @param mat_after Matrix after filtering
#' @param config Configuration list
#' @return List of ggplot objects
create_filtering_plots <- function(mat_before, mat_after, config) {
  plots <- list()

  # Missingness distribution before/after
  miss_before <- data.frame(
    pct_missing = rowMeans(is.na(mat_before)) * 100,
    stage = "Before filtering"
  )
  miss_after <- data.frame(
    pct_missing = rowMeans(is.na(mat_after)) * 100,
    stage = "After filtering"
  )
  miss_df <- rbind(miss_before, miss_after)

  plots$missingness_dist <- ggplot2::ggplot(miss_df, ggplot2::aes(x = pct_missing, fill = stage)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    ggplot2::labs(
      title = "Missingness Distribution",
      x = "% Missing per Feature",
      y = "Count",
      fill = "Stage"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = c("Before filtering" = "#E69F00", "After filtering" = "#56B4E9"))

  # Intensity distribution before/after
  int_before <- data.frame(
    mean_intensity = rowMeans(mat_before, na.rm = TRUE),
    stage = "Before filtering"
  )
  int_after <- data.frame(
    mean_intensity = rowMeans(mat_after, na.rm = TRUE),
    stage = "After filtering"
  )
  int_df <- rbind(int_before, int_after)

  plots$intensity_dist <- ggplot2::ggplot(int_df, ggplot2::aes(x = mean_intensity, fill = stage)) +
    ggplot2::geom_density(alpha = 0.7) +
    ggplot2::labs(
      title = "Intensity Distribution",
      x = "Mean Intensity per Feature",
      y = "Density",
      fill = "Stage"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = c("Before filtering" = "#E69F00", "After filtering" = "#56B4E9"))

  # Save plots
  save_plot(plots$missingness_dist, "filtering_missingness.png", config, subdir = "qc")
  save_plot(plots$intensity_dist, "filtering_intensity.png", config, subdir = "qc")

  return(plots)
}
