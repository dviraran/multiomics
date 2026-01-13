#!/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript
# =============================================================================
# Run All Multi-Omics Pipelines
# =============================================================================
# Master script to execute all four pipelines in sequence using pipeline_runner.R
#
# Usage:
#   Rscript run_all_pipelines.R [options]
#
# Options:
#   --clean             Clean all targets caches before running
#   --no-prompt         Don't prompt for cache action, reuse cache if available
#   --skip-metabolomics Skip the metabolomics pipeline
#   --only=PIPELINE     Run only specified pipeline (rnaseq, proteomics, metabolomics, multiomics)
#   --run-name=NAME     Specify custom run name (default: auto-generated)
#   --config-dir=PATH   Directory containing config files (default: auto-detect)
#
# Outputs are saved to: {current_dir}/runs/{dataset}_{timestamp}/
# =============================================================================

# Find the multiomics root directory
find_root <- function() {
  # Try to get script path (works when run with Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) > 0) {
    script_path <- sub("^--file=", "", script_arg[1])
    script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
    candidate <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
    if (file.exists(file.path(candidate, "pipeline_runner.R"))) {
      return(candidate)
    }
  }

  # Try common locations relative to working directory
  candidates <- c(
    "..",           # If in examples/
    "../..",        # If in examples/scripts/
    ".",            # If in root
    Sys.getenv("MULTIOMICS_ROOT")
  )

  for (candidate in candidates) {
    if (candidate == "") next
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (file.exists(file.path(candidate, "pipeline_runner.R"))) {
      return(candidate)
    }
  }

  stop("Cannot find pipeline_runner.R. Run from the multiomics root or examples directory,\n",
       "or set MULTIOMICS_ROOT environment variable.")
}

multiomics_root <- find_root()

source(file.path(multiomics_root, "pipeline_runner.R"))

# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
clean_first <- "--clean" %in% args
no_prompt <- "--no-prompt" %in% args
skip_metabolomics <- "--skip-metabolomics" %in% args
only_pipeline <- NULL
custom_run_name <- NULL
config_dir <- NULL

for (arg in args) {
  if (grepl("^--only=", arg)) {
    only_pipeline <- sub("^--only=", "", arg)
  }
  if (grepl("^--run-name=", arg)) {
    custom_run_name <- sub("^--run-name=", "", arg)
  }
  if (grepl("^--config-dir=", arg)) {
    config_dir <- sub("^--config-dir=", "", arg)
  }
}

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Detect dataset name from config files
detect_dataset <- function(config_paths) {
  for (config_path in config_paths) {
    if (!file.exists(config_path)) next

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
    }, error = function(e) NULL)
  }
  return("analysis")
}

#' Prompt user for cache action
prompt_cache_action <- function(pipeline_name, pipeline_dir) {
  targets_dir <- file.path(pipeline_dir, "_targets")
  has_cache <- dir.exists(targets_dir) && length(list.files(targets_dir)) > 0

  if (!has_cache) return("run")
  if (no_prompt) return("reuse")
  if (clean_first) return("clean")

  cat(sprintf("\n--- %s pipeline ---\n", toupper(pipeline_name)))
  cat("Previous results found. What would you like to do?\n")
  cat("  1. Reuse cache (only run changed targets)\n")
  cat("  2. Clean and rerun from scratch\n")
  cat("  3. Skip this pipeline\n")

  repeat {
    choice <- readline("Enter choice (1/2/3): ")
    if (choice %in% c("1", "2", "3")) break
    cat("Invalid choice.\n")
  }

  switch(choice, "1" = "reuse", "2" = "clean", "3" = "skip")
}

#' Collect outputs from pipeline to run directory
collect_outputs <- function(pipeline_dir, run_dir, pipeline_name) {
  source_outputs <- file.path(pipeline_dir, "outputs")
  dest_dir <- file.path(run_dir, pipeline_name)

  if (!dir.exists(source_outputs)) {
    cat("  No outputs to collect for", pipeline_name, "\n")
    return(invisible(0))
  }

  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(source_outputs, recursive = TRUE, full.names = TRUE)
  for (f in files) {
    rel_path <- sub(paste0(source_outputs, "/"), "", f)
    dest_path <- file.path(dest_dir, rel_path)
    dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
    file.copy(f, dest_path, overwrite = TRUE)
  }

  cat("  Collected", length(files), "files from", pipeline_name, "\n")
  invisible(length(files))
}

#' Create run directory
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

  if (dir.exists(run_dir) && !no_prompt) {
    cat("\nRun directory exists:", run_dir, "\n")
    cat("  1. Overwrite\n")
    cat("  2. Create new with suffix\n")
    choice <- readline("Enter choice (1/2): ")

    if (choice == "2") {
      i <- 2
      while (dir.exists(paste0(run_dir, "_v", i))) i <- i + 1
      run_dir <- paste0(run_dir, "_v", i)
    } else {
      unlink(run_dir, recursive = TRUE)
    }
  }

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_dir
}

#' Generate run info JSON
generate_run_info <- function(run_dir, dataset, results, start_time, end_time) {
  run_info <- list(
    run_name = basename(run_dir),
    dataset = dataset,
    start_time = as.character(start_time),
    end_time = as.character(end_time),
    duration_minutes = as.numeric(difftime(end_time, start_time, units = "mins")),
    pipelines = lapply(results, function(r) {
      list(name = r$pipeline, status = r$status, duration_minutes = r$duration, error = r$error)
    }),
    system_info = list(
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      platform = R.version$platform
    )
  )

  jsonlite::write_json(run_info, file.path(run_dir, "run_info.json"),
                       auto_unbox = TRUE, pretty = TRUE)
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

cat("
================================================================================
                      Multi-Omics Pipeline Runner
================================================================================
")

start_dir <- getwd()
cat("Working directory:", start_dir, "\n")
cat("Multiomics root:", multiomics_root, "\n\n")

# Define pipelines and their config files
pipeline_names <- c("rnaseq", "proteomics", "metabolomics", "multiomics")

if (!is.null(only_pipeline)) {
  if (!only_pipeline %in% pipeline_names) {
    stop("Invalid pipeline: ", only_pipeline, ". Choose from: ", paste(pipeline_names, collapse = ", "))
  }
  pipeline_names <- only_pipeline
}

if (skip_metabolomics) {
  pipeline_names <- setdiff(pipeline_names, "metabolomics")
}

# Build pipeline info
pipelines <- list()
config_files <- list()

# Resolve config_dir to absolute path if provided
if (!is.null(config_dir)) {
  # If relative path, resolve relative to current working directory
  if (!grepl("^/", config_dir) && !grepl("^~", config_dir)) {
    config_dir <- normalizePath(file.path(start_dir, config_dir), mustWork = FALSE)
  } else {
    config_dir <- normalizePath(config_dir, mustWork = FALSE)
  }
  cat("Config directory:", config_dir, "\n")

  if (!dir.exists(config_dir)) {
    stop("Config directory not found: ", config_dir)
  }
}

for (name in pipeline_names) {
  pipeline_dir <- file.path(multiomics_root, paste0(name, "_pipeline"))
  pipelines[[name]] <- pipeline_dir

  # Try to find config file
  if (!is.null(config_dir)) {
    cfg <- file.path(config_dir, paste0(name, "_config.yml"))
    if (!file.exists(cfg)) {
      # Also try without _config suffix
      cfg <- file.path(config_dir, paste0(name, ".yml"))
    }
    if (!file.exists(cfg)) {
      cat("  Warning: No config found for", name, "in", config_dir, "\n")
      cfg <- file.path(pipeline_dir, "config.yml")
    }
  } else {
    cfg <- file.path(pipeline_dir, "config.yml")
  }
  config_files[[name]] <- if (file.exists(cfg)) normalizePath(cfg) else NULL
}

# Detect dataset
dataset_name <- detect_dataset(unlist(config_files))
cat("Detected dataset:", dataset_name, "\n")

# Create run directory
run_dir <- create_run_directory(start_dir, dataset_name, custom_run_name)
cat("Run directory:", run_dir, "\n\n")

# Prompt for cache actions (before execution)
cat("================================================================================\n")
cat("                           CACHE STATUS CHECK\n")
cat("================================================================================\n")

cache_actions <- list()
for (name in pipeline_names) {
  cache_actions[[name]] <- prompt_cache_action(name, pipelines[[name]])
}

# Execute pipelines
cat("\n================================================================================\n")
cat("                           EXECUTING PIPELINES\n")
cat("================================================================================\n")

results <- list()
overall_start <- Sys.time()

for (name in pipeline_names) {
  action <- cache_actions[[name]]

  if (action == "skip") {
    cat("\nSkipping", name, "pipeline\n")
    results[[name]] <- list(pipeline = name, status = "skipped", duration = 0, error = NULL)
    next
  }

  clean <- (action == "clean")
  config <- config_files[[name]]

  cat("\n")
  start_time <- Sys.time()

  # Use run_omics_pipeline from pipeline_runner.R
  result <- tryCatch({
    run_omics_pipeline(name, config, clean)
  }, error = function(e) {
    list(status = "error", error = e$message)
  })

  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "mins"))

  results[[name]] <- list(
    pipeline = name,
    status = result$status,
    duration = duration,
    error = result$error %||% NULL
  )

  # Collect outputs
  collect_outputs(pipelines[[name]], run_dir, name)
}

overall_end <- Sys.time()
overall_duration <- difftime(overall_end, overall_start, units = "mins")

# -----------------------------------------------------------------------------
# Summary Report
# -----------------------------------------------------------------------------

cat("\n================================================================================\n")
cat("                           PIPELINE SUMMARY\n")
cat("================================================================================\n\n")

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

  icon <- switch(r$status, "success" = "+", "skipped" = "o", "-")
  cat(sprintf("  [%s] %-15s %s (%.1f min)\n", icon, r$pipeline, r$status, r$duration))
}

cat("\nTotal execution time:", round(as.numeric(overall_duration), 2), "minutes\n")

n_success <- sum(sapply(results, function(x) x$status == "success"))
n_skipped <- sum(sapply(results, function(x) x$status == "skipped"))
cat("Pipelines completed:", n_success, "/", length(results))
if (n_skipped > 0) cat(" (", n_skipped, " skipped)", sep = "")
cat("\n")

# Show errors
errors <- Filter(function(x) !is.null(x$error), results)
if (length(errors) > 0) {
  cat("\nErrors:\n")
  for (r in errors) {
    cat("  -", r$pipeline, ":", r$error, "\n")
  }
}

# Save outputs
generate_run_info(run_dir, dataset_name, results, overall_start, overall_end)
write.csv(summary_df, file.path(run_dir, "pipeline_summary.csv"), row.names = FALSE)

# Output locations
cat("\n================================================================================\n")
cat("                           OUTPUT LOCATIONS\n")
cat("================================================================================\n\n")

cat("All outputs saved to:\n  ", run_dir, "\n\n")

for (name in pipeline_names) {
  if (results[[name]]$status == "skipped") next
  report_dir <- file.path(run_dir, name)
  if (dir.exists(report_dir)) {
    reports <- list.files(report_dir, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
    if (length(reports) > 0) {
      cat(name, "report:\n  ", reports[1], "\n")
    }
  }
}

cat("\nRun summary:", file.path(run_dir, "pipeline_summary.csv"), "\n")
cat("Run metadata:", file.path(run_dir, "run_info.json"), "\n")

cat("\nDone!\n")
