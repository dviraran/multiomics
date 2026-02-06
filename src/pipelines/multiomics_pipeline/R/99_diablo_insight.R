#!/usr/bin/env Rscript

# =============================================================================
# DIABLO Insight Reporter
# =============================================================================
# This script loads the existing DIABLO results from the targets store or output
# files and generates a concise "Executive Summary" of the top driving features.
# It is designed to be run manually to get quick insights.
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(targets)
    library(mixOmics)
})

# Set working directory to project root if needed
# (Adjust this path if running interactively from elsewhere)
if (dir.exists("src/pipelines/multiomics_pipeline")) {
    setwd("src/pipelines/multiomics_pipeline")
}

cat("================================================================================\n")
cat("                       DIABLO INSIGHT REPORT                                    \n")
cat("================================================================================\n\n")

# 1. Load Data
# -----------------------------------------------------------------------------
cat("Loading DIABLO results...\n")

tryCatch(
    {
        # Try loading from targets store first
        targets::tar_load(integration_results)

        if (is.null(integration_results$diablo)) {
            stop("DIABLO results not found in integration_results")
        }

        diablo_res <- integration_results$diablo
        model <- diablo_res$model # This might differ depending on how it was saved

        # If 'model' is not directly in the list, we might need to rely on the saved CSVs
        # But let's check the structure based on 07_diablo.R
        # Result structure: list(model = diablo_result, results = diablo_results, ...)

        if ("result" %in% names(diablo_res)) {
            # Sometimes wrapped differently
            full_model <- diablo_res$result
        } else if ("model" %in% names(diablo_res)) {
            full_model <- diablo_res$model
        } else {
            # Fallback: if we only have the 'results' list with loadings/scores
            full_model <- NULL
        }

        loadings_list <- diablo_res$results$loadings
    },
    error = function(e) {
        cat("\nError loading targets: ", e$message, "\n")
        cat("Attempting to read from output CSVs...\n")
        # Fallback to CSV reading would go here, but let's assume targets works for now
        q(status = 1)
    }
)

if (is.null(loadings_list)) {
    cat("Could not retrieve feature loadings.\n")
    q(status = 1)
}

# 2. Extract Top Drivers (Component 1)
# -----------------------------------------------------------------------------
cat("\nAnalyzing Component 1 (Primary Separation)...\n")
cat("These are the features most responsible for distinguishing your conditions.\n\n")

omics_layers <- setdiff(names(loadings_list), "Y")

for (layer in omics_layers) {
    cat(sprintf("--- %s ---\n", toupper(layer)))

    # Get loadings for this layer
    # Loadings format: rows = features, cols = components
    layer_loadings <- loadings_list[[layer]]

    if (is.data.frame(layer_loadings)) {
        layer_loadings <- as.matrix(layer_loadings)
    }

    # Check if Component 1 exists
    if (ncol(layer_loadings) >= 1) {
        # Extract Comp 1
        comp1 <- layer_loadings[, 1]

        # Sort by absolute value
        top_indices <- order(abs(comp1), decreasing = TRUE)[1:10]
        top_features <- comp1[top_indices]

        # Print table
        print_df <- data.frame(
            Rank = 1:10,
            Feature = names(top_features),
            Importance = round(abs(top_features), 4),
            Direction = ifelse(top_features > 0, "(+)", "(-)")
        )

        print(print_df, row.names = FALSE)
        cat("\n")
    } else {
        cat("  No Component 1 found.\n")
    }
}

# 3. Stability / Performance Check
# -----------------------------------------------------------------------------
# Check if error rates exist
cat("--- PERFORMANCE METRICS ---\n")
if (!is.null(diablo_res$results$performance)) {
    perf <- diablo_res$results$performance

    if (!is.null(perf$optimal_ncomp)) {
        cat(sprintf("Optimal Number of Components: %s\n", perf$optimal_ncomp))
    }

    if (!is.null(perf$error_rates)) {
        if (is.list(perf$error_rates)) {
            # It is a list of matrices (one per distance metric)
            # We want the minimum error from the first column (overall error)
            min_vals <- sapply(perf$error_rates, function(x) {
                if (is.matrix(x)) min(x[, 1], na.rm = TRUE) else min(x, na.rm = TRUE)
            })
            min_err <- min(min_vals, na.rm = TRUE)
        } else if (is.data.frame(perf$error_rates) && "error_rate" %in% names(perf$error_rates)) {
            min_err <- min(perf$error_rates$error_rate, na.rm = TRUE)
        } else {
            min_err <- min(perf$error_rates, na.rm = TRUE)
        }

        if (is.finite(min_err)) {
            cat(sprintf("Best Overall Error Rate: %.2f%%\n", min_err * 100))
        }
    }
} else {
    cat("Performance metrics not available in the saved object.\n")
}

cat("\n================================================================================\n")
cat("Interpretation Tip: High importance features in differnet layers often\n")
cat("biologically interact. Check 'outputs/plots/diablo_circos.png' to see connections.\n")
cat("================================================================================\n")
