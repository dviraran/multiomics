#' Pre-flight Data Validation Checks
#'
#' This module provides functions to validate data files before running pipelines.
#' Catches data format issues, sample mismatches, and potential problems early.
#'
#' @description
#' Features:
#' - Validate matrix format (numeric, no NA columns, etc.)
#' - Check sample ID matching between data and metadata
#' - Detect potential batch effects
#' - Sample size warnings
#' - Organism/ID type validation

# ==============================================================================
# MAIN PREFLIGHT CHECK FUNCTION
# ==============================================================================

#' Run all pre-flight checks for a pipeline
#'
#' @param config Configuration list
#' @param pipeline_type One of "rnaseq", "proteomics", "metabolomics", "multiomics"
#' @param base_path Base path for resolving file paths
#' @param verbose Print progress
#'
#' @return List with check results
#' @export
run_preflight_checks <- function(config, pipeline_type, base_path = ".", verbose = TRUE) {

    result <- list(
        passed = TRUE,
        errors = character(0),
        warnings = character(0),
        info = character(0),
        stats = list()
    )

    if (verbose) message("\n=== Running Pre-flight Checks ===\n")

    # Run checks based on pipeline type
    result <- switch(pipeline_type,
        "rnaseq" = preflight_rnaseq(config, base_path, result, verbose),
        "proteomics" = preflight_proteomics(config, base_path, result, verbose),
        "metabolomics" = preflight_metabolomics(config, base_path, result, verbose),
        "multiomics" = preflight_multiomics(config, base_path, result, verbose),
        result
    )

    result$passed <- length(result$errors) == 0

    # Print summary
    if (verbose) {
        message("\n=== Pre-flight Check Summary ===")

        if (length(result$errors) > 0) {
            message("\nERRORS (must fix):")
            for (e in result$errors) message(sprintf("  ✗ %s", e))
        }

        if (length(result$warnings) > 0) {
            message("\nWARNINGS:")
            for (w in result$warnings) message(sprintf("  ⚠ %s", w))
        }

        if (length(result$info) > 0) {
            message("\nINFO:")
            for (i in result$info) message(sprintf("  ℹ %s", i))
        }

        if (result$passed) {
            message("\n✓ All pre-flight checks passed")
        } else {
            message("\n✗ Pre-flight checks failed - please fix errors before running")
        }
    }

    result
}

# ==============================================================================
# RNASEQ PREFLIGHT
# ==============================================================================

preflight_rnaseq <- function(config, base_path, result, verbose) {

    # Load counts
    counts_file <- resolve_path(config$counts_file, base_path)
    if (!is.null(counts_file) && file.exists(counts_file)) {
        if (verbose) message("Checking counts file...")

        counts <- tryCatch({
            read_data_file(counts_file)
        }, error = function(e) {
            result$errors <<- c(result$errors, sprintf("Failed to read counts file: %s", e$message))
            NULL
        })

        if (!is.null(counts)) {
            result <- check_expression_matrix(counts, "counts", result, verbose)
            result$stats$n_genes <- nrow(counts)
            result$stats$n_samples_counts <- ncol(counts)

            # Check for raw counts (should be integers or close to it)
            numeric_cols <- sapply(counts, is.numeric)
            if (any(numeric_cols)) {
                sample_vals <- unlist(counts[, numeric_cols, drop = FALSE][1:min(100, nrow(counts)), ])
                sample_vals <- sample_vals[!is.na(sample_vals) & sample_vals > 0]
                if (length(sample_vals) > 0) {
                    frac_integer <- mean(sample_vals == floor(sample_vals))
                    if (frac_integer < 0.9) {
                        result$warnings <- c(result$warnings,
                            "Counts appear to be non-integer (may be normalized already)")
                    }
                }
            }
        }
    }

    # Load metadata
    metadata_file <- resolve_path(config$metadata_file, base_path)
    if (!is.null(metadata_file) && file.exists(metadata_file)) {
        if (verbose) message("Checking metadata file...")

        metadata <- tryCatch({
            read_data_file(metadata_file)
        }, error = function(e) {
            result$errors <<- c(result$errors, sprintf("Failed to read metadata file: %s", e$message))
            NULL
        })

        if (!is.null(metadata)) {
            result <- check_metadata(metadata, config, result, verbose)
            result$stats$n_samples_metadata <- nrow(metadata)

            # Check sample matching
            if (!is.null(counts)) {
                result <- check_sample_matching(counts, metadata, config, result, verbose)
            }

            # Check group sizes
            result <- check_group_sizes(metadata, config, result, verbose)
        }
    }

    # Check gene IDs
    if (!is.null(counts)) {
        gene_ids <- get_feature_ids(counts)
        result <- check_gene_ids(gene_ids, config, result, verbose)
    }

    result
}

# ==============================================================================
# PROTEOMICS PREFLIGHT
# ==============================================================================

preflight_proteomics <- function(config, base_path, result, verbose) {

    # Get file path (nested or flat)
    intensity_file <- config$input$quant_matrix %||%
                      config$intensity_file %||%
                      config$input$intensity_file

    intensity_file <- resolve_path(intensity_file, base_path)

    if (!is.null(intensity_file) && file.exists(intensity_file)) {
        if (verbose) message("Checking intensity file...")

        intensity <- tryCatch({
            read_data_file(intensity_file)
        }, error = function(e) {
            result$errors <<- c(result$errors, sprintf("Failed to read intensity file: %s", e$message))
            NULL
        })

        if (!is.null(intensity)) {
            result <- check_expression_matrix(intensity, "intensity", result, verbose)
            result$stats$n_proteins <- nrow(intensity)

            # Check for zeros (common in proteomics)
            numeric_cols <- sapply(intensity, is.numeric)
            if (any(numeric_cols)) {
                n_zeros <- sum(intensity[, numeric_cols] == 0, na.rm = TRUE)
                n_total <- sum(numeric_cols) * nrow(intensity)
                pct_zeros <- 100 * n_zeros / n_total

                if (pct_zeros > 50) {
                    result$warnings <- c(result$warnings,
                        sprintf("%.1f%% zeros in intensity matrix - check if zeros should be NA", pct_zeros))
                }

                result$stats$pct_zeros <- pct_zeros
            }

            # Check for missing values
            pct_na <- 100 * sum(is.na(intensity)) / (nrow(intensity) * ncol(intensity))
            result$stats$pct_missing <- pct_na

            if (pct_na > 30) {
                result$info <- c(result$info,
                    sprintf("%.1f%% missing values - imputation will be applied", pct_na))
            }
        }
    }

    # Load metadata
    metadata_file <- config$input$metadata %||% config$metadata_file
    metadata_file <- resolve_path(metadata_file, base_path)

    if (!is.null(metadata_file) && file.exists(metadata_file)) {
        if (verbose) message("Checking metadata file...")

        metadata <- tryCatch({
            read_data_file(metadata_file)
        }, error = function(e) {
            result$errors <<- c(result$errors, sprintf("Failed to read metadata: %s", e$message))
            NULL
        })

        if (!is.null(metadata)) {
            result <- check_metadata(metadata, config, result, verbose)

            # Check sample matching
            if (!is.null(intensity)) {
                result <- check_sample_matching(intensity, metadata, config, result, verbose)
            }
        }
    }

    result
}

# ==============================================================================
# METABOLOMICS PREFLIGHT
# ==============================================================================

preflight_metabolomics <- function(config, base_path, result, verbose) {

    # Get file path
    feature_file <- config$input$feature_matrix %||% config$feature_file
    feature_file <- resolve_path(feature_file, base_path)

    if (!is.null(feature_file) && file.exists(feature_file)) {
        if (verbose) message("Checking feature matrix...")

        features <- tryCatch({
            read_data_file(feature_file)
        }, error = function(e) {
            result$errors <<- c(result$errors, sprintf("Failed to read feature file: %s", e$message))
            NULL
        })

        if (!is.null(features)) {
            result <- check_expression_matrix(features, "features", result, verbose)
            result$stats$n_features <- nrow(features)

            # Check for negative values (shouldn't exist in metabolomics)
            numeric_cols <- sapply(features, is.numeric)
            if (any(numeric_cols)) {
                n_negative <- sum(features[, numeric_cols] < 0, na.rm = TRUE)
                if (n_negative > 0) {
                    result$warnings <- c(result$warnings,
                        sprintf("%d negative values detected - unusual for metabolomics data", n_negative))
                }
            }
        }
    }

    # Load metadata
    metadata_file <- config$input$metadata %||% config$metadata_file
    metadata_file <- resolve_path(metadata_file, base_path)

    if (!is.null(metadata_file) && file.exists(metadata_file)) {
        if (verbose) message("Checking metadata file...")

        metadata <- tryCatch({
            read_data_file(metadata_file)
        }, error = function(e) {
            result$errors <<- c(result$errors, sprintf("Failed to read metadata: %s", e$message))
            NULL
        })

        if (!is.null(metadata)) {
            result <- check_metadata(metadata, config, result, verbose)

            # Check for sample type column (QC, blank, sample)
            sample_type_col <- config$sample_types$sample_type_column %||% "sample_type"
            if (sample_type_col %in% names(metadata)) {
                sample_types <- unique(metadata[[sample_type_col]])
                result$info <- c(result$info,
                    sprintf("Sample types found: %s", paste(sample_types, collapse = ", ")))

                # Check for QC samples if batch correction enabled
                if (!is.null(config$batch_correction$method) &&
                    config$batch_correction$method == "loess") {
                    if (!any(grepl("qc", sample_types, ignore.case = TRUE))) {
                        result$warnings <- c(result$warnings,
                            "Loess batch correction enabled but no QC samples detected")
                    }
                }
            }
        }
    }

    result
}

# ==============================================================================
# MULTIOMICS PREFLIGHT
# ==============================================================================

preflight_multiomics <- function(config, base_path, result, verbose) {

    omics_present <- config$omics_present

    if (is.null(omics_present) || length(omics_present) < 2) {
        result$errors <- c(result$errors, "Need at least 2 omics types")
        return(result)
    }

    sample_ids_per_omics <- list()

    # Check each omics type
    for (omics in omics_present) {
        if (verbose) message(sprintf("\nChecking %s data...", omics))

        omics_config <- config[[omics]]
        if (is.null(omics_config)) {
            result$errors <- c(result$errors, sprintf("Missing config section for: %s", omics))
            next
        }

        # Get data file
        data_file <- omics_config$counts_file %||%
                     omics_config$intensity_file %||%
                     omics_config$feature_matrix

        data_file <- resolve_path(data_file, base_path)

        if (!is.null(data_file) && file.exists(data_file)) {
            data <- tryCatch({
                read_data_file(data_file)
            }, error = function(e) {
                result$errors <<- c(result$errors,
                    sprintf("Failed to read %s data: %s", omics, e$message))
                NULL
            })

            if (!is.null(data)) {
                result$stats[[paste0(omics, "_features")]] <- nrow(data)
                result$stats[[paste0(omics, "_samples")]] <- ncol(data) - 1  # Assuming first col is ID

                # Store sample IDs
                sample_ids_per_omics[[omics]] <- names(data)[-1]
            }
        }
    }

    # Check sample overlap
    if (length(sample_ids_per_omics) >= 2) {
        if (verbose) message("\nChecking sample overlap across omics...")

        all_samples <- Reduce(union, sample_ids_per_omics)
        common_samples <- Reduce(intersect, sample_ids_per_omics)

        result$stats$total_unique_samples <- length(all_samples)
        result$stats$common_samples <- length(common_samples)

        if (length(common_samples) == 0) {
            result$errors <- c(result$errors,
                "No common samples found across omics types - check sample ID formats")
        } else {
            result$info <- c(result$info,
                sprintf("Found %d common samples across %d omics types",
                        length(common_samples), length(omics_present)))

            # Warn if low overlap
            min_samples <- min(sapply(sample_ids_per_omics, length))
            if (length(common_samples) < min_samples * 0.5) {
                result$warnings <- c(result$warnings,
                    sprintf("Low sample overlap: only %d of %d samples are shared",
                            length(common_samples), min_samples))
            }
        }
    }

    # Check minimum sample size for integration
    if (!is.null(result$stats$common_samples) && result$stats$common_samples < 10) {
        result$warnings <- c(result$warnings,
            sprintf("Only %d common samples - integration results may be unstable",
                    result$stats$common_samples))
    }

    result
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Resolve a file path relative to base path
resolve_path <- function(file_path, base_path) {
    if (is.null(file_path)) return(NULL)

    if (file.exists(file_path)) {
        return(file_path)
    }

    full_path <- file.path(base_path, file_path)
    if (file.exists(full_path)) {
        return(full_path)
    }

    file_path
}

#' Read a data file (CSV or TSV)
read_data_file <- function(file_path) {
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
        read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (grepl("\\.tsv$|\\.txt$", file_path, ignore.case = TRUE)) {
        read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
        # Try CSV first, then TSV
        tryCatch({
            read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
        }, error = function(e) {
            read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
        })
    }
}

#' Check expression/intensity matrix format
check_expression_matrix <- function(data, name, result, verbose) {

    if (verbose) message(sprintf("  Matrix dimensions: %d rows x %d columns", nrow(data), ncol(data)))

    # Check for numeric columns
    numeric_cols <- sapply(data, is.numeric)
    n_numeric <- sum(numeric_cols)

    if (n_numeric == 0) {
        result$errors <- c(result$errors,
            sprintf("%s matrix has no numeric columns", name))
    } else if (n_numeric < ncol(data) - 1) {
        result$info <- c(result$info,
            sprintf("%s matrix: %d numeric columns (samples), %d non-numeric (ID columns)",
                    name, n_numeric, ncol(data) - n_numeric))
    }

    # Check for empty rows/columns
    if (any(numeric_cols)) {
        numeric_data <- data[, numeric_cols, drop = FALSE]

        # Rows with all NA
        all_na_rows <- apply(numeric_data, 1, function(x) all(is.na(x)))
        if (sum(all_na_rows) > 0) {
            result$warnings <- c(result$warnings,
                sprintf("%d rows with all missing values in %s matrix", sum(all_na_rows), name))
        }

        # Columns with all NA
        all_na_cols <- apply(numeric_data, 2, function(x) all(is.na(x)))
        if (sum(all_na_cols) > 0) {
            result$errors <- c(result$errors,
                sprintf("%d columns with all missing values in %s matrix", sum(all_na_cols), name))
        }
    }

    result
}

#' Check metadata format
check_metadata <- function(metadata, config, result, verbose) {

    if (verbose) message(sprintf("  Metadata: %d samples, %d columns", nrow(metadata), ncol(metadata)))

    # Check for sample ID column
    sample_id_col <- config$sample_id_column %||%
                     config$sample_id_col %||%
                     config$input$sample_id_column %||%
                     "sample_id"

    if (!(sample_id_col %in% names(metadata))) {
        # Try to find a likely sample ID column
        potential_id_cols <- c("sample_id", "SampleID", "Sample_ID", "sample", "Sample",
                               "sampleid", "ID", "id", "name", "Name")
        found <- intersect(potential_id_cols, names(metadata))

        if (length(found) > 0) {
            result$info <- c(result$info,
                sprintf("Sample ID column '%s' not found, but found: %s",
                        sample_id_col, paste(found, collapse = ", ")))
        } else {
            result$warnings <- c(result$warnings,
                sprintf("Sample ID column '%s' not found in metadata", sample_id_col))
        }
    }

    # Check for group/condition column
    group_col <- config$group_col %||%
                 config$condition_column %||%
                 config$input$condition_column %||%
                 "condition"

    if (!(group_col %in% names(metadata))) {
        potential_group_cols <- c("condition", "Condition", "group", "Group", "treatment",
                                  "Treatment", "class", "Class", "phenotype")
        found <- intersect(potential_group_cols, names(metadata))

        if (length(found) > 0) {
            result$info <- c(result$info,
                sprintf("Group column '%s' not found, but found: %s",
                        group_col, paste(found, collapse = ", ")))
        } else {
            result$warnings <- c(result$warnings,
                sprintf("Group/condition column '%s' not found in metadata", group_col))
        }
    } else {
        # Report group levels
        groups <- unique(metadata[[group_col]])
        result$info <- c(result$info,
            sprintf("Groups found: %s (%d levels)", paste(groups, collapse = ", "), length(groups)))
    }

    result
}

#' Check sample ID matching between data and metadata
check_sample_matching <- function(data, metadata, config, result, verbose) {

    # Get sample IDs from data (column names, excluding first which is usually feature ID)
    data_samples <- names(data)
    # Exclude obvious non-sample columns
    exclude_patterns <- c("^gene", "^feature", "^protein", "^id$", "^name$",
                          "^symbol$", "^description$", "^entrez")
    is_sample <- !grepl(paste(exclude_patterns, collapse = "|"), data_samples, ignore.case = TRUE)
    data_samples <- data_samples[is_sample]

    # Get sample IDs from metadata
    sample_id_col <- config$sample_id_column %||%
                     config$sample_id_col %||%
                     config$input$sample_id_column

    if (is.null(sample_id_col) || !(sample_id_col %in% names(metadata))) {
        # Try first column
        sample_id_col <- names(metadata)[1]
    }

    metadata_samples <- as.character(metadata[[sample_id_col]])

    # Check overlap
    common <- intersect(data_samples, metadata_samples)
    only_data <- setdiff(data_samples, metadata_samples)
    only_metadata <- setdiff(metadata_samples, data_samples)

    if (verbose) {
        message(sprintf("  Sample matching: %d in both, %d only in data, %d only in metadata",
                        length(common), length(only_data), length(only_metadata)))
    }

    if (length(common) == 0) {
        result$errors <- c(result$errors,
            "No matching sample IDs between data and metadata - check ID formats")

        # Show examples
        if (length(data_samples) > 0 && length(metadata_samples) > 0) {
            result$info <- c(result$info,
                sprintf("Data sample IDs look like: %s",
                        paste(head(data_samples, 3), collapse = ", ")))
            result$info <- c(result$info,
                sprintf("Metadata sample IDs look like: %s",
                        paste(head(metadata_samples, 3), collapse = ", ")))
        }
    } else if (length(only_data) > 0 || length(only_metadata) > 0) {
        if (length(only_data) > 0) {
            result$warnings <- c(result$warnings,
                sprintf("%d samples in data but not in metadata (will be excluded)",
                        length(only_data)))
        }
        if (length(only_metadata) > 0) {
            result$info <- c(result$info,
                sprintf("%d samples in metadata but not in data", length(only_metadata)))
        }
    }

    result$stats$n_matched_samples <- length(common)

    result
}

#' Check group sizes for statistical validity
check_group_sizes <- function(metadata, config, result, verbose) {

    group_col <- config$group_col %||%
                 config$condition_column %||%
                 config$input$condition_column

    if (is.null(group_col) || !(group_col %in% names(metadata))) {
        return(result)
    }

    group_sizes <- table(metadata[[group_col]])

    if (verbose) {
        message("  Group sizes:")
        for (g in names(group_sizes)) {
            message(sprintf("    %s: %d", g, group_sizes[g]))
        }
    }

    # Check for small groups
    small_groups <- names(group_sizes)[group_sizes < 3]
    if (length(small_groups) > 0) {
        result$warnings <- c(result$warnings,
            sprintf("Groups with n < 3: %s (statistical tests may be unreliable)",
                    paste(small_groups, collapse = ", ")))
    }

    # Check for very unbalanced groups
    if (length(group_sizes) >= 2) {
        ratio <- max(group_sizes) / min(group_sizes)
        if (ratio > 5) {
            result$warnings <- c(result$warnings,
                sprintf("Highly unbalanced groups (ratio %.1f:1) - consider this in interpretation",
                        ratio))
        }
    }

    result$stats$group_sizes <- as.list(group_sizes)

    result
}

#' Get feature IDs from a data frame
get_feature_ids <- function(data) {
    # First column is usually feature ID
    first_col <- data[[1]]

    # Check if it looks like IDs (character, mostly unique)
    if (is.character(first_col) || is.factor(first_col)) {
        return(as.character(first_col))
    }

    # Check row names
    if (!is.null(rownames(data)) && !all(rownames(data) == as.character(1:nrow(data)))) {
        return(rownames(data))
    }

    # Fall back to first column
    as.character(first_col)
}

#' Check gene ID format and type
check_gene_ids <- function(gene_ids, config, result, verbose) {

    # Source organism detection if available
    tryCatch({
        shared_dir <- file.path(dirname(dirname(getwd())), "shared", "R")
        if (file.exists(file.path(shared_dir, "organism_detection.R"))) {
            source(file.path(shared_dir, "organism_detection.R"))

            id_type <- detect_id_type(gene_ids)
            if (verbose) message(sprintf("  Detected ID type: %s", id_type))
            result$stats$detected_id_type <- id_type

            # Check for version suffixes
            if (has_ensembl_version(gene_ids)) {
                result$info <- c(result$info,
                    "Ensembl IDs have version suffixes (e.g., .1) - will be stripped during annotation")
            }

            # Detect organism
            org_result <- detect_organism_from_ids(gene_ids)
            if (org_result$organism != "unknown") {
                if (verbose) message(sprintf("  Detected organism: %s (confidence: %s)",
                                              org_result$organism, org_result$confidence))
                result$stats$detected_organism <- org_result$organism

                # Check if matches config
                config_organism <- config$organism
                if (!is.null(config_organism) &&
                    tolower(org_result$organism) != tolower(config_organism) &&
                    org_result$confidence == "high") {
                    result$warnings <- c(result$warnings,
                        sprintf("Detected organism (%s) differs from config (%s)",
                                org_result$organism, config_organism))
                }
            }
        }
    }, error = function(e) {
        # Silently skip organism detection if utilities not available
    })

    # Check for duplicate IDs
    n_dup <- sum(duplicated(gene_ids))
    if (n_dup > 0) {
        result$info <- c(result$info,
            sprintf("%d duplicate feature IDs (will be aggregated)", n_dup))
    }

    result
}

# ==============================================================================
# NULL COALESCING
# ==============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
