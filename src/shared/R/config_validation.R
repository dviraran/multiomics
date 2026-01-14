#' Configuration Validation Utilities
#'
#' This module provides functions to validate pipeline configuration files
#' before running analyses, catching errors early and providing helpful messages.
#'
#' @description
#' Features:
#' - Schema validation for each pipeline type
#' - File existence checks
#' - Parameter range validation
#' - Cross-parameter consistency checks
#' - Helpful error messages with suggestions

# ==============================================================================
# MAIN VALIDATION FUNCTION
# ==============================================================================

#' Validate a pipeline configuration
#'
#' @param config List containing the configuration (from yaml::read_yaml)
#' @param pipeline_type One of "rnaseq", "proteomics", "metabolomics", "multiomics"
#' @param base_path Base path for resolving relative file paths
#' @param verbose Print validation progress
#'
#' @return List with:
#'   - valid: Logical, whether config passed all checks
#'   - errors: Character vector of critical errors (will fail)
#'   - warnings: Character vector of warnings (may cause issues)
#'   - info: Character vector of informational messages
#'
#' @export
validate_config <- function(config, pipeline_type, base_path = ".", verbose = TRUE) {

    result <- list(
        valid = TRUE,
        errors = character(0),
        warnings = character(0),
        info = character(0)
    )

    if (verbose) message(sprintf("\n=== Validating %s pipeline configuration ===\n", pipeline_type))

    # Validate based on pipeline type
    result <- switch(pipeline_type,
        "rnaseq" = validate_rnaseq_config(config, base_path, result, verbose),
        "proteomics" = validate_proteomics_config(config, base_path, result, verbose),
        "metabolomics" = validate_metabolomics_config(config, base_path, result, verbose),
        "multiomics" = validate_multiomics_config(config, base_path, result, verbose),
        {
            result$errors <- c(result$errors, sprintf("Unknown pipeline type: %s", pipeline_type))
            result$valid <- FALSE
            result
        }
    )

    # Set valid flag based on errors
    result$valid <- length(result$errors) == 0

    # Print summary
    if (verbose) {
        message("\n=== Validation Summary ===")

        if (length(result$errors) > 0) {
            message("\nERRORS (must fix):")
            for (e in result$errors) message(sprintf("  ✗ %s", e))
        }

        if (length(result$warnings) > 0) {
            message("\nWARNINGS (may cause issues):")
            for (w in result$warnings) message(sprintf("  ⚠ %s", w))
        }

        if (length(result$info) > 0) {
            message("\nINFO:")
            for (i in result$info) message(sprintf("  ℹ %s", i))
        }

        if (result$valid) {
            message("\n✓ Configuration is valid")
        } else {
            message("\n✗ Configuration has errors - please fix before running pipeline")
        }
    }

    result
}

# ==============================================================================
# RNASEQ CONFIG VALIDATION
# ==============================================================================

validate_rnaseq_config <- function(config, base_path, result, verbose) {

    # Required fields
    required <- c("counts_file", "metadata_file")
    result <- check_required_fields(config, required, result)

    # File existence
    result <- check_file_exists(config, "counts_file", base_path, result, verbose)
    result <- check_file_exists(config, "metadata_file", base_path, result, verbose)
    result <- check_file_exists(config, "custom_gene_mapping", base_path, result, verbose, required = FALSE)
    result <- check_file_exists(config, "custom_gmt_file", base_path, result, verbose, required = FALSE)

    # Check organism specification
    result <- check_organism(config, result, verbose)

    # Column specifications
    result <- check_column_spec(config, "sample_id_column", result)
    result <- check_column_spec(config, "group_col", result)

    # Numeric parameters
    result <- check_numeric_range(config, "alpha", 0, 1, result)
    result <- check_numeric_range(config, "lfc_threshold", 0, 10, result)
    result <- check_numeric_range(config, "min_count", 0, 1000, result)

    # Contrasts
    result <- check_contrasts(config, result, verbose)

    # Pathway parameters
    result <- check_pathway_config(config, result, verbose)

    # Design formula
    result <- check_design_formula(config, result)

    result
}

# ==============================================================================
# PROTEOMICS CONFIG VALIDATION
# ==============================================================================

validate_proteomics_config <- function(config, base_path, result, verbose) {

    # Check nested input structure
    input <- config$input %||% config

    # Required fields
    result <- check_file_exists_nested(config, c("input", "quant_matrix"), base_path, result, verbose)
    result <- check_file_exists_nested(config, c("input", "metadata"), base_path, result, verbose)

    # Alternative flat structure
    if (is.null(config$input)) {
        result <- check_file_exists(config, "intensity_file", base_path, result, verbose, required = FALSE)
        result <- check_file_exists(config, "metadata_file", base_path, result, verbose, required = FALSE)
    }

    # Check organism
    result <- check_organism(config, result, verbose)

    # Processing parameters
    processing <- config$processing %||% config
    if (!is.null(processing$normalization_method)) {
        valid_methods <- c("median", "vsn", "quantile", "none")
        if (!(processing$normalization_method %in% valid_methods)) {
            result$warnings <- c(result$warnings,
                sprintf("Unrecognized normalization_method: %s (valid: %s)",
                        processing$normalization_method, paste(valid_methods, collapse = ", ")))
        }
    }

    # Imputation
    imputation <- config$imputation %||% config
    if (!is.null(imputation$method)) {
        valid_methods <- c("QRILC", "KNN", "half_min", "mean", "median", "none")
        if (!(imputation$method %in% valid_methods)) {
            result$warnings <- c(result$warnings,
                sprintf("Unrecognized imputation method: %s", imputation$method))
        }
    }

    # Filtering thresholds
    filtering <- config$filtering %||% config
    result <- check_numeric_range_nested(config, c("filtering", "global_min_presence"), 0, 1, result)
    result <- check_numeric_range_nested(config, c("filtering", "group_min_presence"), 0, 1, result)

    result
}

# ==============================================================================
# METABOLOMICS CONFIG VALIDATION
# ==============================================================================

validate_metabolomics_config <- function(config, base_path, result, verbose) {

    # Check nested input structure
    result <- check_file_exists_nested(config, c("input", "feature_matrix"), base_path, result, verbose)
    result <- check_file_exists_nested(config, c("input", "metadata"), base_path, result, verbose)

    # Alternative flat structure
    if (is.null(config$input)) {
        result <- check_file_exists(config, "feature_file", base_path, result, verbose, required = FALSE)
        result <- check_file_exists(config, "metadata_file", base_path, result, verbose, required = FALSE)
    }

    # Data type
    if (!is.null(config$input$data_type)) {
        valid_types <- c("untargeted", "targeted")
        if (!(config$input$data_type %in% valid_types)) {
            result$warnings <- c(result$warnings,
                sprintf("Unrecognized data_type: %s (valid: %s)",
                        config$input$data_type, paste(valid_types, collapse = ", ")))
        }
    }

    # Normalization
    processing <- config$processing %||% config
    if (!is.null(processing$normalization_method)) {
        valid_methods <- c("PQN", "median", "quantile", "none")
        if (!(processing$normalization_method %in% valid_methods)) {
            result$info <- c(result$info,
                sprintf("Non-standard normalization method: %s", processing$normalization_method))
        }
    }

    # Batch correction
    batch <- config$batch_correction %||% list()
    if (!is.null(batch$method) && batch$method == "loess") {
        # Loess requires QC samples
        result$info <- c(result$info,
            "Loess batch correction requires QC samples in the data")
    }

    result
}

# ==============================================================================
# MULTIOMICS CONFIG VALIDATION
# ==============================================================================

validate_multiomics_config <- function(config, base_path, result, verbose) {

    # Check omics_present
    omics_present <- config$omics_present
    if (is.null(omics_present) || length(omics_present) < 2) {
        result$errors <- c(result$errors,
            "multiomics pipeline requires at least 2 omics types in 'omics_present'")
    } else {
        valid_omics <- c("transcriptomics", "proteomics", "metabolomics", "rnaseq")
        invalid <- setdiff(omics_present, valid_omics)
        if (length(invalid) > 0) {
            result$warnings <- c(result$warnings,
                sprintf("Unrecognized omics types: %s", paste(invalid, collapse = ", ")))
        }
    }

    # Metadata file
    result <- check_file_exists(config, "metadata", base_path, result, verbose)

    # Per-omics data files
    for (omics in omics_present %||% character(0)) {
        omics_config <- config[[omics]]
        if (is.null(omics_config)) {
            result$errors <- c(result$errors,
                sprintf("Missing configuration section for omics type: %s", omics))
        } else {
            # Check for data file
            if (!is.null(omics_config$counts_file)) {
                result <- check_file_exists_in_section(config, omics, "counts_file", base_path, result, verbose)
            }
            if (!is.null(omics_config$intensity_file)) {
                result <- check_file_exists_in_section(config, omics, "intensity_file", base_path, result, verbose)
            }
            if (!is.null(omics_config$feature_matrix)) {
                result <- check_file_exists_in_section(config, omics, "feature_matrix", base_path, result, verbose)
            }
        }
    }

    # Integration methods
    integration <- config$integration %||% list()
    if (!is.null(integration$methods)) {
        valid_methods <- c("MOFA2", "DIABLO", "SNF", "SPLS", "correlation")
        invalid <- setdiff(integration$methods, valid_methods)
        if (length(invalid) > 0) {
            result$warnings <- c(result$warnings,
                sprintf("Unrecognized integration methods: %s", paste(invalid, collapse = ", ")))
        }
    }

    # MOFA2 parameters
    if (!is.null(integration$mofa2)) {
        result <- check_numeric_range_nested(config, c("integration", "mofa2", "num_factors"), 1, 50, result)
    }

    # DIABLO parameters
    if (!is.null(integration$diablo)) {
        result <- check_numeric_range_nested(config, c("integration", "diablo", "ncomp"), 1, 20, result)
    }

    result
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Check that required fields are present
check_required_fields <- function(config, required, result) {
    for (field in required) {
        if (is.null(config[[field]])) {
            result$errors <- c(result$errors, sprintf("Missing required field: %s", field))
        }
    }
    result
}

#' Check if a file exists
check_file_exists <- function(config, field, base_path, result, verbose, required = TRUE) {
    file_path <- config[[field]]

    if (is.null(file_path)) {
        if (required) {
            result$errors <- c(result$errors, sprintf("Missing required file path: %s", field))
        }
        return(result)
    }

    # Resolve relative path
    if (!file.exists(file_path)) {
        full_path <- file.path(base_path, file_path)
        if (!file.exists(full_path)) {
            if (required) {
                result$errors <- c(result$errors,
                    sprintf("File not found: %s (field: %s)", file_path, field))
            } else {
                result$warnings <- c(result$warnings,
                    sprintf("Optional file not found: %s (field: %s)", file_path, field))
            }
        } else if (verbose) {
            result$info <- c(result$info, sprintf("Found: %s", full_path))
        }
    } else if (verbose) {
        result$info <- c(result$info, sprintf("Found: %s", file_path))
    }

    result
}

#' Check file exists in nested config
check_file_exists_nested <- function(config, path, base_path, result, verbose, required = TRUE) {
    value <- config
    for (p in path) {
        value <- value[[p]]
        if (is.null(value)) break
    }

    field_name <- paste(path, collapse = ".")

    if (is.null(value)) {
        if (required && length(path) <= 2) {
            # Only error if it's a primary required field
            result$warnings <- c(result$warnings, sprintf("Missing nested field: %s", field_name))
        }
        return(result)
    }

    # Check file existence
    if (!file.exists(value)) {
        full_path <- file.path(base_path, value)
        if (!file.exists(full_path)) {
            if (required) {
                result$errors <- c(result$errors,
                    sprintf("File not found: %s (field: %s)", value, field_name))
            }
        }
    }

    result
}

#' Check file in omics section
check_file_exists_in_section <- function(config, section, field, base_path, result, verbose) {
    file_path <- config[[section]][[field]]
    if (!is.null(file_path) && !file.exists(file_path)) {
        full_path <- file.path(base_path, file_path)
        if (!file.exists(full_path)) {
            result$errors <- c(result$errors,
                sprintf("File not found in %s.%s: %s", section, field, file_path))
        }
    }
    result
}

#' Check organism specification
check_organism <- function(config, result, verbose) {
    organism <- config$organism

    if (is.null(organism)) {
        result$info <- c(result$info,
            "No organism specified - will attempt auto-detection or skip annotation")
        return(result)
    }

    # Check if it's a known organism
    known_organisms <- c(
        "Homo sapiens", "human",
        "Mus musculus", "mouse",
        "Rattus norvegicus", "rat",
        "Danio rerio", "zebrafish",
        "Drosophila melanogaster", "fly",
        "Caenorhabditis elegans", "worm",
        "Saccharomyces cerevisiae", "yeast",
        "Arabidopsis thaliana"
    )

    organism_lower <- tolower(organism)
    if (!(organism_lower %in% tolower(known_organisms))) {
        result$info <- c(result$info,
            sprintf("Non-model organism: %s - annotation will require custom mapping file or biomaRt", organism))

        # Check for custom mapping
        if (is.null(config$custom_gene_mapping) && is.null(config$annotation$custom_mapping_file)) {
            result$info <- c(result$info,
                "Consider providing custom_gene_mapping or annotation.custom_mapping_file")
        }
    }

    result
}

#' Check column specification
check_column_spec <- function(config, field, result) {
    value <- config[[field]]
    if (is.null(value)) {
        result$info <- c(result$info, sprintf("Using default for: %s", field))
    }
    result
}

#' Check numeric parameter is in valid range
check_numeric_range <- function(config, field, min_val, max_val, result) {
    value <- config[[field]]
    if (!is.null(value)) {
        if (!is.numeric(value)) {
            result$errors <- c(result$errors,
                sprintf("%s must be numeric (got: %s)", field, class(value)))
        } else if (value < min_val || value > max_val) {
            result$warnings <- c(result$warnings,
                sprintf("%s = %s is outside typical range [%s, %s]",
                        field, value, min_val, max_val))
        }
    }
    result
}

#' Check numeric range in nested config
check_numeric_range_nested <- function(config, path, min_val, max_val, result) {
    value <- config
    for (p in path) {
        value <- value[[p]]
        if (is.null(value)) return(result)
    }

    field_name <- paste(path, collapse = ".")
    if (!is.numeric(value)) {
        result$errors <- c(result$errors,
            sprintf("%s must be numeric (got: %s)", field_name, class(value)))
    } else if (value < min_val || value > max_val) {
        result$warnings <- c(result$warnings,
            sprintf("%s = %s is outside typical range [%s, %s]",
                    field_name, value, min_val, max_val))
    }

    result
}

#' Check contrasts specification
check_contrasts <- function(config, result, verbose) {
    contrasts <- config$contrasts

    if (is.null(contrasts) || length(contrasts) == 0) {
        result$info <- c(result$info,
            "No contrasts specified - will auto-generate pairwise comparisons")
        return(result)
    }

    for (i in seq_along(contrasts)) {
        contrast <- contrasts[[i]]

        # Check format: should be "A - B" or list with numerator/denominator
        if (is.character(contrast)) {
            if (!grepl("\\s*-\\s*", contrast)) {
                result$warnings <- c(result$warnings,
                    sprintf("Contrast '%s' should be in format 'group1 - group2'", contrast))
            }
        } else if (is.list(contrast)) {
            if (is.null(contrast$numerator) || is.null(contrast$denominator)) {
                result$warnings <- c(result$warnings,
                    sprintf("Contrast %d should have 'numerator' and 'denominator' fields", i))
            }
        }
    }

    result
}

#' Check pathway configuration
check_pathway_config <- function(config, result, verbose) {
    # Check for custom GMT
    has_custom_gmt <- !is.null(config$custom_gmt_file) ||
                      !is.null(config$pathway$custom_gmt_file) ||
                      !is.null(config$enrichment$custom_gmt_file)

    # Check pathway database
    pathway_db <- config$pathway_database %||%
                  config$pathway$database %||%
                  config$enrichment$databases

    if (is.null(pathway_db) && !has_custom_gmt) {
        result$info <- c(result$info,
            "No pathway database or custom GMT specified - enrichment will use defaults (GO, KEGG)")
    }

    # Check pathway sizes
    min_size <- config$pathway_min_size %||% config$pathway$min_size
    max_size <- config$pathway_max_size %||% config$pathway$max_size

    if (!is.null(min_size) && !is.null(max_size) && min_size >= max_size) {
        result$errors <- c(result$errors,
            "pathway_min_size must be less than pathway_max_size")
    }

    result
}

#' Check design formula
check_design_formula <- function(config, result) {
    formula <- config$design_formula

    if (!is.null(formula)) {
        # Basic syntax check
        if (!grepl("^~", formula)) {
            result$warnings <- c(result$warnings,
                sprintf("design_formula should start with '~' (got: %s)", formula))
        }

        # Check for common issues
        if (grepl("\\+\\s*$", formula) || grepl("^~\\s*$", formula)) {
            result$errors <- c(result$errors,
                sprintf("Invalid design formula: %s", formula))
        }
    }

    result
}

# ==============================================================================
# NULL COALESCING OPERATOR
# ==============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

# ==============================================================================
# CONVENIENCE FUNCTION
# ==============================================================================

#' Quick validation check - returns TRUE/FALSE
#'
#' @param config Configuration list
#' @param pipeline_type Pipeline type
#' @param base_path Base path for files
#'
#' @return Logical
#' @export
is_valid_config <- function(config, pipeline_type, base_path = ".") {
    result <- validate_config(config, pipeline_type, base_path, verbose = FALSE)
    result$valid
}
