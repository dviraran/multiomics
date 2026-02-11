#' MOFA+ Analysis Wrapper for R
#'
#' This file provides R functions to run MOFA+ (Multi-Omics Factor Analysis)
#' using the Python mofapy2 package via a Python script.
#'
#' @author Generated for multiomics pipeline
#' @date 2026-02-06

#' Run MOFA+ analysis on multiple omics views
#'
#' @param views Named list of file paths to CSV files containing omics data.
#'   Each CSV should have features as rows and samples as columns.
#'   Example: list(rna = "data/rna.csv", protein = "data/prot.csv")
#' @param outfile Path to save the trained MOFA model (.hdf5 file)
#' @param outdir Directory to save CSV outputs (factors, weights, R2)
#' @param n_factors Number of factors to train (default: 10)
#' @param seed Random seed for reproducibility (default: 42)
#' @param max_iter Maximum number of training iterations (default: 1000)
#' @param convergence_mode Convergence mode: "fast", "medium", or "slow" (default: "fast")
#' @param scale_views Logical, whether to scale views to unit variance (default: FALSE)
#' @param verbose Logical, whether to print verbose output (default: TRUE)
#' @param python_exec Path to Python executable (default: "python3")
#'
#' @return List containing:
#'   \item{success}{Logical indicating if MOFA ran successfully}
#'   \item{model_file}{Path to saved HDF5 model file}
#'   \item{factors}{Data frame of latent factors (samples x factors)}
#'   \item{weights}{Named list of data frames with weights for each view}
#'   \item{variance_explained}{Data frame of variance explained (R2) per view}
#'
#' @examples
#' \dontrun{
#' # Run MOFA with two omics views
#' result <- run_mofa(
#'   views = list(
#'     rna = "data/abundance_RNA.csv",
#'     protein = "data/abundance_protein.csv"
#'   ),
#'   outfile = "results/mofa_model.hdf5",
#'   outdir = "results/mofa_outputs",
#'   n_factors = 15,
#'   seed = 42
#' )
#'
#' # Access results
#' factors <- result$factors
#' rna_weights <- result$weights$rna
#' }
#'
#' @export
run_mofa <- function(views,
                     outfile,
                     outdir,
                     n_factors = 10,
                     seed = 42,
                     max_iter = 1000,
                     convergence_mode = c("fast", "medium", "slow"),
                     scale_views = FALSE,
                     verbose = TRUE,
                     python_exec = "python3") {

  # Input validation
  if (!is.list(views) || length(views) == 0) {
    stop("views must be a non-empty named list")
  }

  if (is.null(names(views)) || any(names(views) == "")) {
    stop("views must be a named list (e.g., list(rna = 'file.csv', protein = 'file2.csv'))")
  }

  # Check that all view files exist
  for (view_name in names(views)) {
    if (!file.exists(views[[view_name]])) {
      stop("View file not found: ", views[[view_name]])
    }
  }

  # Check Python script exists
  script_path <- "scripts/run_mofa.py"
  if (!file.exists(script_path)) {
    stop("MOFA script not found: ", script_path,
         "\nPlease ensure run_mofa.py is in the scripts directory")
  }

  convergence_mode <- match.arg(convergence_mode)

  # Create output directory
  if (!dir.exists(dirname(outfile))) {
    dir.create(dirname(outfile), recursive = TRUE)
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # Build view arguments
  view_args <- sapply(names(views), function(name) {
    sprintf("%s=%s", name, views[[name]])
  })
  view_args_str <- paste(view_args, collapse = " ")

  # Build command
  cmd <- sprintf(
    "%s '%s' --views %s --outfile '%s' --outdir '%s' --factors %d --seed %d --max_iter %d --convergence_mode %s",
    python_exec,
    script_path,
    view_args_str,
    outfile,
    outdir,
    n_factors,
    seed,
    max_iter,
    convergence_mode
  )

  if (scale_views) {
    cmd <- paste(cmd, "--scale_views")
  }

  if (verbose) {
    cmd <- paste(cmd, "--verbose")
  }

  if (verbose) {
    message("Running MOFA+ analysis...")
    message("  Views: ", paste(names(views), collapse = ", "))
    message("  Factors: ", n_factors)
    message("  Output: ", outdir)
    message()
  }

  # Execute command
  exit_code <- system(cmd, intern = FALSE)

  # Check if successful
  if (exit_code != 0) {
    stop("MOFA script failed with exit code ", exit_code)
  }

  # Load outputs
  if (verbose) {
    message("\nLoading MOFA results...")
  }

  result <- list(
    success = TRUE,
    model_file = outfile
  )

  # Load factors
  factors_file <- file.path(outdir, "factors.csv")
  if (file.exists(factors_file)) {
    result$factors <- read.csv(factors_file, row.names = 1)
    if (verbose) {
      message("  Loaded factors: ", nrow(result$factors), " samples x ",
              ncol(result$factors), " factors")
    }
  } else {
    warning("Factors file not found: ", factors_file)
  }

  # Load weights for each view
  result$weights <- list()
  for (view_name in names(views)) {
    weights_file <- file.path(outdir, paste0("weights_", view_name, ".csv"))
    if (file.exists(weights_file)) {
      result$weights[[view_name]] <- read.csv(weights_file, row.names = 1)
      if (verbose) {
        message("  Loaded weights for ", view_name, ": ",
                nrow(result$weights[[view_name]]), " features x ",
                ncol(result$weights[[view_name]]), " factors")
      }
    } else {
      warning("Weights file not found for view ", view_name, ": ", weights_file)
    }
  }

  # Load variance explained
  r2_total_file <- file.path(outdir, "variance_explained_total.csv")
  r2_per_factor_file <- file.path(outdir, "variance_explained_per_factor.csv")

  if (file.exists(r2_total_file)) {
    result$variance_explained_total <- read.csv(r2_total_file, row.names = 1)
    if (verbose) {
      message("\n  Variance explained (R2):")
      for (view_name in names(views)) {
        if (view_name %in% colnames(result$variance_explained_total)) {
          r2 <- result$variance_explained_total["Total_R2", view_name]
          message("    ", view_name, ": ", round(r2, 2), "%")
        }
      }
    }
  }

  if (file.exists(r2_per_factor_file)) {
    result$variance_explained_per_factor <- read.csv(r2_per_factor_file, row.names = 1)
  }

  if (verbose) {
    message("\nMOFA analysis completed successfully!")
  }

  return(result)
}


#' Prepare abundance matrix for MOFA
#'
#' Helper function to prepare an abundance matrix (e.g., from DESeq2, edgeR)
#' for MOFA analysis by saving it as a CSV file.
#'
#' @param mat Matrix or data.frame with features as rows and samples as columns
#' @param outfile Path to save CSV file
#' @param feature_col Optional column name for feature IDs (if mat is a data.frame)
#'
#' @return Path to saved file (invisibly)
#'
#' @export
prepare_mofa_view <- function(mat, outfile, feature_col = NULL) {

  # Convert to data.frame if needed
  if (is.matrix(mat)) {
    df <- as.data.frame(mat)
  } else if (is.data.frame(mat)) {
    df <- mat
  } else {
    stop("mat must be a matrix or data.frame")
  }

  # Set feature IDs as rownames if specified
  if (!is.null(feature_col)) {
    if (!(feature_col %in% colnames(df))) {
      stop("feature_col '", feature_col, "' not found in data.frame")
    }
    rownames(df) <- df[[feature_col]]
    df[[feature_col]] <- NULL
  }

  # Ensure we have row names
  if (is.null(rownames(df)) || any(rownames(df) == "")) {
    stop("mat must have row names (feature IDs)")
  }

  # Create output directory if needed
  if (!dir.exists(dirname(outfile))) {
    dir.create(dirname(outfile), recursive = TRUE)
  }

  # Write CSV
  write.csv(df, outfile, row.names = TRUE)

  message("Saved MOFA view: ", outfile)
  message("  ", nrow(df), " features x ", ncol(df), " samples")

  invisible(outfile)
}


#' Plot MOFA factor variance explained
#'
#' Create a barplot showing variance explained by each factor per view
#'
#' @param mofa_result Result object from run_mofa()
#' @param top_n Number of top factors to plot (default: all)
#'
#' @return ggplot object
#'
#' @export
plot_mofa_variance <- function(mofa_result, top_n = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("tidyr package required for plotting")
  }

  if (is.null(mofa_result$variance_explained_per_factor)) {
    stop("Variance explained data not found in MOFA result")
  }

  df <- mofa_result$variance_explained_per_factor

  # Subset to top N factors if requested
  if (!is.null(top_n) && top_n < nrow(df)) {
    # Sum R2 across views for each factor
    factor_totals <- rowSums(df, na.rm = TRUE)
    top_factors <- names(sort(factor_totals, decreasing = TRUE)[1:top_n])
    df <- df[top_factors, , drop = FALSE]
  }

  # Convert to long format
  df$Factor <- rownames(df)
  df_long <- tidyr::pivot_longer(
    df,
    cols = -Factor,
    names_to = "View",
    values_to = "R2"
  )

  # Create plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Factor, y = R2, fill = View)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Variance Explained by MOFA Factors",
      x = "Factor",
      y = "Variance Explained (R²%)",
      fill = "View"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  return(p)
}
