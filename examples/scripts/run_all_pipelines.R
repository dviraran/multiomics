#!/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript
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
#   Rscript run_all_pipelines.R [--parallel] [--skip-metabolomics] [--clean] [--no-prompt]
#
# Options:
#   --parallel          Run single-omics pipelines in parallel (requires future)
#   --skip-metabolomics Skip the metabolomics pipeline
#   --clean             Clean all targets caches before running (skip prompt)
#   --no-prompt         Don't prompt for cache action, use cache if available
#   --only=PIPELINE     Run only specified pipeline (rnaseq, proteomics, metabolomics, multiomics)
#   --run-name=NAME     Specify custom run name (default: auto-generated)
#
# Outputs are saved to: {current_dir}/runs/{dataset}_{timestamp}/
# =============================================================================

library(targets)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
run_parallel <- "--parallel" %in% args
skip_metabolomics <- "--skip-metabolomics" %in% args
clean_first <- "--clean" %in% args
no_prompt <- "--no-prompt" %in% args
only_pipeline <- NULL
custom_run_name <- NULL

for (arg in args) {
  if (grepl("^--only=", arg)) {
    only_pipeline <- sub("^--only=", "", arg)
  }
  if (grepl("^--run-name=", arg)) {
    custom_run_name <- sub("^--run-name=", "", arg)
  }
}

# Directories
start_dir <- getwd()
multiomics_root <- normalizePath(file.path(dirname(start_dir), ".."), mustWork = FALSE)

# Try to find multiomics root by looking for rnaseq_pipeline
if (!dir.exists(file.path(multiomics_root, "rnaseq_pipeline"))) {
  # Maybe we're already in the root or examples
  if (dir.exists(file.path(start_dir, "rnaseq_pipeline"))) {
    multiomics_root <- start_dir
  } else if (dir.exists(file.path(dirname(start_dir), "rnaseq_pipeline"))) {
    multiomics_root <- dirname(start_dir)
  } else {
    stop("Cannot find multiomics pipelines. Please run from the examples directory or project root.")
  }
}

cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║           Multi-Omics Pipeline Runner                            ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

cat("Root directory:", multiomics_root, "\n")
cat("Working directory:", start_dir, "\n")
cat("Parallel mode:", run_parallel, "\n")
cat("Skip metabolomics:", skip_metabolomics, "\n\n")

# -----------------------------------------------------------------------------
# Helper Functions for Run Management
# -----------------------------------------------------------------------------

#' Detect dataset name from pipeline config
#' @param config_path Path to config.yml
#' @return Dataset name (e.g., "nci60", "stategra")
detect_dataset <- function(config_path) {
  if (!file.exists(config_path)) return("unknown")

  tryCatch({
    config <- yaml::read_yaml(config_path)

    # Try to extract from input file paths
    input_paths <- c(

      config$counts_file,
      config$metadata_file,
      config$input_file,
      config$intensity_file
    )

    for (path in input_paths) {
      if (!is.null(path)) {
        # Look for known dataset names in path
        if (grepl("nci60", path, ignore.case = TRUE)) return("nci60")
        if (grepl("stategra", path, ignore.case = TRUE)) return("stategra")
        # Extract from path pattern like "data/datasetname/"
        match <- regmatches(path, regexpr("data/([^/]+)/", path))
        if (length(match) > 0) {
          dataset <- gsub("data/|/", "", match)
          if (nchar(dataset) > 0) return(dataset)
        }
      }
    }

    return("analysis")
  }, error = function(e) {
    return("analysis")
  })
}

#' Prompt user for cache action
#' @param pipeline_name Name of the pipeline
#' @param targets_dir Path to _targets directory
#' @return Action: "reuse", "clean", or "skip"
prompt_cache_action <- function(pipeline_name, targets_dir) {
  has_cache <- dir.exists(targets_dir) && length(list.files(targets_dir)) > 0

  if (!has_cache) {
    return("run")
  }

  if (no_prompt) {
    return("reuse")
  }

  if (clean_first) {
    return("clean")
  }

  cat(sprintf("\n┌─ %s pipeline ─────────────────────────────────────\n", toupper(pipeline_name)))
  cat("│ Previous results found in cache.\n")
  cat("│ What would you like to do?\n")
  cat("│   1. Reuse cache (only run changed targets)\n")
  cat("│   2. Clean and rerun from scratch\n")
  cat("│   3. Skip this pipeline\n")
  cat("└────────────────────────────────────────────────────\n")

  repeat {
    choice <- readline("Enter choice (1/2/3): ")
    if (choice %in% c("1", "2", "3")) break
    cat("Invalid choice. Please enter 1, 2, or 3.\n")
  }

  switch(choice,
    "1" = "reuse",
    "2" = "clean",
    "3" = "skip"
  )
}

#' Collect outputs from pipeline to run directory
#' @param pipeline_dir Source pipeline directory
#' @param run_dir Destination run directory
#' @param pipeline_name Name of the pipeline
collect_outputs <- function(pipeline_dir, run_dir, pipeline_name) {
  source_outputs <- file.path(pipeline_dir, "outputs")
  dest_dir <- file.path(run_dir, pipeline_name)

  if (!dir.exists(source_outputs)) {
    cat("  No outputs to collect for", pipeline_name, "\n")
    return(invisible(NULL))
  }

  # Create destination directory

dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  # Copy all outputs
  files <- list.files(source_outputs, recursive = TRUE, full.names = TRUE)
  for (f in files) {
    rel_path <- sub(paste0(source_outputs, "/"), "", f)
    dest_path <- file.path(dest_dir, rel_path)
    dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
    file.copy(f, dest_path, overwrite = TRUE)
  }

  n_files <- length(files)
  cat("  Collected", n_files, "files from", pipeline_name, "\n")

  invisible(n_files)
}

#' Generate run info JSON
#' @param run_dir Run directory
#' @param dataset Dataset name
#' @param results Pipeline results list
generate_run_info <- function(run_dir, dataset, results, start_time, end_time) {
  run_info <- list(
    run_name = basename(run_dir),
    dataset = dataset,
    start_time = as.character(start_time),
    end_time = as.character(end_time),
    duration_minutes = as.numeric(difftime(end_time, start_time, units = "mins")),
    pipelines = lapply(results, function(r) {
      list(
        name = r$pipeline,
        status = r$status,
        duration_minutes = r$duration,
        error = r$errors
      )
    }),
    system_info = list(
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      platform = R.version$platform,
      hostname = Sys.info()["nodename"]
    )
  )

  json_path <- file.path(run_dir, "run_info.json")
  jsonlite::write_json(run_info, json_path, auto_unbox = TRUE, pretty = TRUE)
  cat("Run info saved to:", json_path, "\n")
}

#' Create run directory with timestamp
#' @param base_dir Base directory for runs
#' @param dataset Dataset name
#' @param custom_name Optional custom name
#' @return Path to created run directory
create_run_directory <- function(base_dir, dataset, custom_name = NULL) {
  runs_dir <- file.path(base_dir, "runs")
  dir.create(runs_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(custom_name)) {
    run_name <- custom_name
  } else {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
    run_name <- paste0(dataset, "_", timestamp)
  }

  run_dir <- file.path(runs_dir, run_name)

  # Handle existing directory
  if (dir.exists(run_dir)) {
    cat("\nRun directory already exists:", run_dir, "\n")
    cat("  1. Overwrite\n")
    cat("  2. Create new with suffix\n")
    choice <- readline("Enter choice (1/2): ")

    if (choice == "2") {
      i <- 2
      while (dir.exists(paste0(run_dir, "_v", i))) {
        i <- i + 1
      }
      run_dir <- paste0(run_dir, "_v", i)
    } else {
      unlink(run_dir, recursive = TRUE)
    }
  }

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nRun directory:", run_dir, "\n\n")

  return(run_dir)
}

# -----------------------------------------------------------------------------
# Pipeline Execution Function
# -----------------------------------------------------------------------------

run_pipeline <- function(pipeline_name, pipeline_dir, run_dir = NULL, cache_action = NULL) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("  Running:", toupper(pipeline_name), "pipeline\n")
  cat("  Directory:", pipeline_dir, "\n")
  cat(rep("=", 70), "\n\n", sep = "")

  # Determine cache action if not provided
  targets_dir <- file.path(pipeline_dir, "_targets")
  if (is.null(cache_action)) {
    cache_action <- prompt_cache_action(pipeline_name, targets_dir)
  }

  # Handle skip
  if (cache_action == "skip") {
    cat("Skipping", pipeline_name, "pipeline\n")
    return(list(
      pipeline = pipeline_name,
      status = "skipped",
      duration = 0,
      errors = NULL
    ))
  }

  start_time <- Sys.time()

  # Change to pipeline directory
  original_dir <- getwd()
  setwd(pipeline_dir)

  tryCatch({
    # Clean if requested
    if (cache_action == "clean") {
      cat("Cleaning targets cache...\n")
      tar_destroy(ask = FALSE)
    }

    # Run the pipeline
    cat("Starting pipeline...\n\n")
    tar_make()

    end_time <- Sys.time()
    duration <- difftime(end_time, start_time, units = "mins")

    cat("\n✓", pipeline_name, "pipeline completed in", round(duration, 2), "minutes\n")

    # Collect outputs to run directory if specified
    if (!is.null(run_dir)) {
      collect_outputs(pipeline_dir, run_dir, pipeline_name)
    }

    # Return success
    return(list(
      pipeline = pipeline_name,
      status = "success",
      duration = as.numeric(duration),
      errors = NULL
    ))

  }, error = function(e) {
    cat("\n✗ Error in", pipeline_name, "pipeline:", e$message, "\n")

    # Still collect partial outputs if run_dir specified
    if (!is.null(run_dir)) {
      tryCatch({
        collect_outputs(pipeline_dir, run_dir, pipeline_name)
      }, error = function(e2) NULL)
    }

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
# Detect Dataset and Create Run Directory
# -----------------------------------------------------------------------------

# Detect dataset from first available pipeline config
dataset_name <- "analysis"
for (pname in c("rnaseq", "proteomics", "metabolomics")) {
  if (pname %in% names(pipelines)) {
    config_path <- file.path(pipelines[[pname]], "config.yml")
    dataset_name <- detect_dataset(config_path)
    if (dataset_name != "unknown" && dataset_name != "analysis") break
  }
}

cat("Detected dataset:", dataset_name, "\n")

# Create run directory in current working directory
run_dir <- create_run_directory(start_dir, dataset_name, custom_run_name)

# -----------------------------------------------------------------------------
# Prompt for Cache Actions (before execution)
# -----------------------------------------------------------------------------

# If not in parallel mode, prompt for each pipeline upfront
cache_actions <- list()
if (!run_parallel) {
  cat("═══════════════════════════════════════════════════════════════════\n")
  cat("                     CACHE STATUS CHECK                            \n")
  cat("═══════════════════════════════════════════════════════════════════\n")

  for (name in names(pipelines)) {
    targets_dir <- file.path(pipelines[[name]], "_targets")
    cache_actions[[name]] <- prompt_cache_action(name, targets_dir)
  }
  cat("\n")
}

# -----------------------------------------------------------------------------
# Execute Pipelines
# -----------------------------------------------------------------------------

results <- list()
overall_start <- Sys.time()

if (run_parallel && length(pipelines) > 1 && is.null(only_pipeline)) {
  # Run single-omics pipelines in parallel, then multiomics
  cat("Running single-omics pipelines in parallel...\n")
  cat("Note: In parallel mode, cache will be reused if available.\n\n")

  if (requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)) {

    library(future)
    library(future.apply)

    # Set up parallel workers
    plan(multisession, workers = min(3, length(pipelines) - 1))

    # Run single-omics pipelines in parallel
    single_omics <- pipelines[names(pipelines) != "multiomics"]

    # Determine cache action for parallel runs
    parallel_cache_action <- if (clean_first) "clean" else "reuse"

    single_results <- future_lapply(names(single_omics), function(name) {
      run_pipeline(name, single_omics[[name]], run_dir, parallel_cache_action)
    }, future.seed = TRUE)

    names(single_results) <- names(single_omics)
    results <- c(results, single_results)

    # Reset to sequential
    plan(sequential)

    # Run multiomics last (depends on single-omics outputs)
    if ("multiomics" %in% names(pipelines)) {
      results$multiomics <- run_pipeline("multiomics", pipelines$multiomics, run_dir, parallel_cache_action)
    }

  } else {
    cat("Warning: future package not available. Running sequentially.\n")
    run_parallel <- FALSE
  }
}

if (!run_parallel) {
  # Run pipelines sequentially with pre-determined cache actions
  for (name in names(pipelines)) {
    results[[name]] <- run_pipeline(name, pipelines[[name]], run_dir, cache_actions[[name]])
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

  status_icon <- if (r$status == "success") "✓" else if (r$status == "skipped") "○" else "✗"
  cat(sprintf("  %s %-15s %s (%.1f min)\n",
              status_icon, r$pipeline, r$status,
              ifelse(is.na(r$duration), 0, r$duration)))
}

cat("\n")
cat("Total execution time:", round(overall_duration, 2), "minutes\n")

# Count successes/failures
n_success <- sum(sapply(results, function(x) x$status == "success"))
n_skipped <- sum(sapply(results, function(x) x$status == "skipped"))
n_total <- length(results)

cat("Pipelines completed:", n_success, "/", n_total)
if (n_skipped > 0) cat(" (", n_skipped, " skipped)", sep = "")
cat("\n")

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
# Save Run Info and Summary
# -----------------------------------------------------------------------------

# Generate run_info.json
generate_run_info(run_dir, dataset_name, results, overall_start, overall_end)

# Save summary CSV to run directory
summary_file <- file.path(run_dir, "pipeline_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)

# -----------------------------------------------------------------------------
# Output Locations
# -----------------------------------------------------------------------------

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("                        OUTPUT LOCATIONS                           \n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("All outputs saved to:\n")
cat("  ", run_dir, "\n\n")

# List available reports in run directory
for (name in names(pipelines)) {
  if (results[[name]]$status == "skipped") next

  report_dir <- file.path(run_dir, name)
  if (dir.exists(report_dir)) {
    # Find HTML reports
    reports <- list.files(report_dir, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
    if (length(reports) > 0) {
      cat(name, "report:\n")
      cat("  ", reports[1], "\n")
    }
  }
}

cat("\nRun summary:", summary_file, "\n")
cat("Run metadata:", file.path(run_dir, "run_info.json"), "\n")

cat("\n✓ All done!\n")
