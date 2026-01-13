#!/usr/bin/env Rscript
# =============================================================================
# Run All Multi-Omics Pipelines
# =============================================================================
# Master script to execute all four pipelines in sequence:
# 1. RNA-seq pipeline
# 2. Proteomics pipeline
# 3. Metabolomics pipeline (optional)
# 4. Multi-omics integration pipeline
#
# Usage:
#   Rscript run_all_pipelines.R [--parallel] [--skip-metabolomics] [--clean]
#
# Options:
#   --parallel         Run single-omics pipelines in parallel (requires future)
#   --skip-metabolomics Skip the metabolomics pipeline
#   --clean            Clean all targets caches before running
#   --only=PIPELINE    Run only specified pipeline (rnaseq, proteomics, metabolomics, multiomics)
# =============================================================================

library(targets)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
run_parallel <- "--parallel" %in% args
skip_metabolomics <- "--skip-metabolomics" %in% args
clean_first <- "--clean" %in% args
only_pipeline <- NULL

for (arg in args) {
  if (grepl("^--only=", arg)) {
    only_pipeline <- sub("^--only=", "", arg)
  }
}

# Directories
examples_dir <- getwd()
multiomics_root <- dirname(examples_dir)

cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║           Multi-Omics Pipeline Runner                            ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

cat("Root directory:", multiomics_root, "\n")
cat("Parallel mode:", run_parallel, "\n")
cat("Skip metabolomics:", skip_metabolomics, "\n\n")

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

run_pipeline <- function(pipeline_name, pipeline_dir) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("  Running:", toupper(pipeline_name), "pipeline\n")
  cat("  Directory:", pipeline_dir, "\n")
  cat(rep("=", 70), "\n\n", sep = "")

  start_time <- Sys.time()

  # Change to pipeline directory
  original_dir <- getwd()
  setwd(pipeline_dir)

  tryCatch({
    # Clean if requested
    if (clean_first) {
      cat("Cleaning targets cache...\n")
      tar_destroy(ask = FALSE)
    }

    # Run the pipeline
    cat("Starting pipeline...\n\n")
    tar_make()

    end_time <- Sys.time()
    duration <- difftime(end_time, start_time, units = "mins")

    cat("\n✓", pipeline_name, "pipeline completed in", round(duration, 2), "minutes\n")

    # Return success
    return(list(
      pipeline = pipeline_name,
      status = "success",
      duration = as.numeric(duration),
      errors = NULL
    ))

  }, error = function(e) {
    cat("\n✗ Error in", pipeline_name, "pipeline:", e$message, "\n")
    return(list(
      pipeline = pipeline_name,
      status = "error",
      duration = NA,
      errors = e$message
    ))

  }, finally = {
    setwd(original_dir)
  })
}

# -----------------------------------------------------------------------------
# Define Pipeline Execution Order
# -----------------------------------------------------------------------------

pipelines <- list(
  rnaseq = file.path(multiomics_root, "rnaseq_pipeline"),
  proteomics = file.path(multiomics_root, "proteomics_pipeline"),
  metabolomics = file.path(multiomics_root, "metabolomics_pipeline"),
  multiomics = file.path(multiomics_root, "multiomics_pipeline")
)

# Filter pipelines based on options
if (!is.null(only_pipeline)) {
  if (!only_pipeline %in% names(pipelines)) {
    stop("Unknown pipeline: ", only_pipeline, "\nValid options: ", paste(names(pipelines), collapse = ", "))
  }
  pipelines <- pipelines[only_pipeline]
  cat("Running only:", only_pipeline, "pipeline\n\n")
}

if (skip_metabolomics && is.null(only_pipeline)) {
  pipelines$metabolomics <- NULL
  cat("Skipping metabolomics pipeline\n\n")
}

# -----------------------------------------------------------------------------
# Execute Pipelines
# -----------------------------------------------------------------------------

results <- list()
overall_start <- Sys.time()

if (run_parallel && length(pipelines) > 1 && is.null(only_pipeline)) {
  # Run single-omics pipelines in parallel, then multiomics
  cat("Running single-omics pipelines in parallel...\n")

  if (requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)) {

    library(future)
    library(future.apply)

    # Set up parallel workers
    plan(multisession, workers = min(3, length(pipelines) - 1))

    # Run single-omics pipelines in parallel
    single_omics <- pipelines[names(pipelines) != "multiomics"]

    single_results <- future_lapply(names(single_omics), function(name) {
      run_pipeline(name, single_omics[[name]])
    }, future.seed = TRUE)

    names(single_results) <- names(single_omics)
    results <- c(results, single_results)

    # Reset to sequential
    plan(sequential)

    # Run multiomics last (depends on single-omics outputs)
    if ("multiomics" %in% names(pipelines)) {
      results$multiomics <- run_pipeline("multiomics", pipelines$multiomics)
    }

  } else {
    cat("Warning: future package not available. Running sequentially.\n")
    run_parallel <- FALSE
  }
}

if (!run_parallel) {
  # Run pipelines sequentially
  for (name in names(pipelines)) {
    results[[name]] <- run_pipeline(name, pipelines[[name]])
  }
}

overall_end <- Sys.time()
overall_duration <- difftime(overall_end, overall_start, units = "mins")

# -----------------------------------------------------------------------------
# Summary Report
# -----------------------------------------------------------------------------

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                    PIPELINE EXECUTION SUMMARY                    ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

summary_df <- data.frame(
  Pipeline = character(),
  Status = character(),
  Duration_min = numeric(),
  stringsAsFactors = FALSE
)

for (name in names(results)) {
  r <- results[[name]]
  summary_df <- rbind(summary_df, data.frame(
    Pipeline = r$pipeline,
    Status = r$status,
    Duration_min = round(r$duration, 2),
    stringsAsFactors = FALSE
  ))

  status_icon <- if (r$status == "success") "✓" else "✗"
  cat(sprintf("  %s %-15s %s (%.1f min)\n",
              status_icon, r$pipeline, r$status,
              ifelse(is.na(r$duration), 0, r$duration)))
}

cat("\n")
cat("Total execution time:", round(overall_duration, 2), "minutes\n")

# Count successes/failures
n_success <- sum(sapply(results, function(x) x$status == "success"))
n_total <- length(results)

cat("Pipelines completed:", n_success, "/", n_total, "\n")

# List any errors
errors <- sapply(results, function(x) x$errors)
errors <- errors[!sapply(errors, is.null)]

if (length(errors) > 0) {
  cat("\nErrors encountered:\n")
  for (name in names(errors)) {
    cat("  -", name, ":", errors[[name]], "\n")
  }
}

# -----------------------------------------------------------------------------
# Output Locations
# -----------------------------------------------------------------------------

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("                        OUTPUT LOCATIONS                           \n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

for (name in names(pipelines)) {
  report_path <- file.path(pipelines[[name]], "outputs", "report")
  if (dir.exists(report_path)) {
    reports <- list.files(report_path, pattern = "\\.html$", full.names = TRUE)
    if (length(reports) > 0) {
      cat(name, "report:\n")
      cat("  ", reports[1], "\n\n")
    }
  }
}

# Save summary
summary_file <- file.path(examples_dir, "pipeline_run_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("Summary saved to:", summary_file, "\n")

cat("\n✓ All done!\n")
