# =============================================================================
# MetaboAnalystR integration wrappers (optional)
# =============================================================================

#' Initialize an mSet from a combined sample x feature table using MetaboAnalystR
#'
#' This will attempt to create a new MetaboAnalystR mSet object and populate it
#' by calling Read.TextData and SanityCheckData. Returns the mSet or NULL on
#' failure (with warnings logged).
init_mset_from_table <- function(combined_file, config) {
  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    log_message("MetaboAnalystR not available; skipping MetaboAnalyst-specific steps")
    return(NULL)
  }

  log_message("Initializing MetaboAnalyst mSet from: ", combined_file)

  if (!file.exists(combined_file)) {
    warning("Combined file not found: ", combined_file)
    return(NULL)
  }

  mSet <- NULL
  tryCatch({
    # Create data object for concentration data (conc)
    mSet <- MetaboAnalystR::InitDataObjects("conc", "stat")

    # Read text data (row-wise; discrete/continuous indicated by data_type)
    MetaboAnalystR::Read.TextData(mSet, combined_file, "rowu", config$input$data_type %||% "disc")

    # Basic sanity checks (structure, sample names, missing values)
    MetaboAnalystR::SanityCheckData(mSet)

    log_message("MetaboAnalyst mSet initialized and sanity-checked")
  }, error = function(e) {
    warning("MetaboAnalystR initialization failed: ", e$message)
    mSet <<- NULL
  })

  return(mSet)
}

#' Run MetaboAnalyst preprocessing: ReplaceMin, PreparePrenormData, Normalization
#'
#' These wrappers call MetaboAnalystR functions and save diagnostic plots if
#' available; they fall back to pipeline routines when MetaboAnalystR is not
#' installed or when the underlying functions fail.
metabo_preprocess <- function(mSet, config) {
  if (is.null(mSet)) return(NULL)

  tryCatch({
    if (!is.null(config$metaboanalyst$replace_min) && config$metaboanalyst$replace_min) {
      frac <- config$metaboanalyst$replace_min_fraction %||% 0.5
      MetaboAnalystR::ReplaceMin(mSet, value = frac)
    }

    MetaboAnalystR::PreparePrenormData(mSet)

    # Normalize (try to do log + mean centering if requested)
    norm_method <- config$metaboanalyst$normalization_method %||% "LogNorm"
    MetaboAnalystR::Normalization(mSet, method = norm_method)

    # Create summary plots if requested
    if (isTRUE(config$metaboanalyst$plot_norm_summary)) {
      out <- tryCatch({
        MetaboAnalystR::PlotNormSummary(mSet)
      }, error = function(e) {
        warning("PlotNormSummary failed: ", e$message); NULL
      })

      # Save sample-level normalization plots if available
      if (isTRUE(config$metaboanalyst$plot_sample_norm_summary)) {
        tryCatch({
          MetaboAnalystR::PlotSampleNormSummary(mSet)
        }, error = function(e) {
          warning("PlotSampleNormSummary failed: ", e$message)
        })
      }
    }

    return(mSet)
  }, error = function(e) {
    warning("MetaboAnalyst preprocessing failed: ", e$message)
    return(NULL)
  })
}

#' Extract normalized matrix from MetaboAnalyst mSet
#'
#' This tries MetaboAnalystR getters first, then falls back to pipeline saved file
extract_normalized_from_mset <- function(mSet, config) {
  # Try MetaboAnalystR internal data structures
  if (!is.null(mSet)) {
    safe_try <- tryCatch({
      if (!is.null(mSet$dataSet$norm)) {
        mat <- mSet$dataSet$norm
        # Ensure features x samples
        if (nrow(mat) == 0 || ncol(mat) == 0) stop("Empty normalized matrix in mSet")
        return(mat)
      } else {
        # Try typical getter
        if (exists("GetNormalizedData", where = asNamespace("MetaboAnalystR"), mode = "function")) {
          return(MetaboAnalystR::GetNormalizedData(mSet))
        }
      }
      stop("Normalized data not found in mSet")
    }, error = function(e) {
      warning("Could not extract normalized data from mSet: ", e$message)
      return(NULL)
    })

    if (!is.null(safe_try)) return(safe_try)
  }

  # Fallback: read normalized CSV from pipeline outputs
  norm_file <- file.path(config$output$output_dir, "tables", "normalized_matrix.csv")
  if (file.exists(norm_file)) {
    df <- readr::read_csv(norm_file)
    mat <- as.matrix(dplyr::select(df, -feature_id))
    rownames(mat) <- df$feature_id
    return(mat)
  }

  warning("No normalized matrix available from MetaboAnalyst or pipeline outputs")
  return(NULL)
}

#' Run MetaboAnalyst core analyses (Volcano, PLS, RF) if available
#'
run_metabo_core_analyses <- function(mSet, config) {
  if (is.null(mSet)) return(NULL)

  out <- list()

  # Volcano
  out$volcano <- tryCatch({
    MetaboAnalystR::Volcano.Anal(mSet)
  }, error = function(e) {
    warning("Volcano.Anal failed: ", e$message); NULL
  })

  # PLSR
  out$plsr <- tryCatch({
    MetaboAnalystR::PLSR.Anal(mSet)
  }, error = function(e) {
    warning("PLSR.Anal failed: ", e$message); NULL
  })

  # RF via MetaboAnalyst if available
  out$rf <- tryCatch({
    MetaboAnalystR::RF.Anal(mSet)
  }, error = function(e) {
    warning("RF.Anal failed: ", e$message); NULL
  })

  return(out)
}

#' Export normalized matrix as samples x features TSV
#'
#' @param normalized_data List returned by normalize_data() target
#' @param config Configuration list
#' @return Path to TSV file (or NULL)
export_normalized_samples_tsv <- function(normalized_data, config) {
  if (is.null(normalized_data) || is.null(normalized_data$matrix)) {
    warning("No normalized data to export")
    return(NULL)
  }

  mat <- normalized_data$matrix
  # Convert to samples x features
  df <- as.data.frame(t(mat))
  df <- tibble::rownames_to_column(df, var = config$input$sample_id_column %||% "sample_id")

  out_name <- config$export$normalized_samples_tsv_name %||% "normalized_matrix_samples_x_features.tsv"
  out_path <- file.path(config$output$output_dir, "tables", out_name)

  readr::write_tsv(df, out_path)
  log_message("Saved: ", out_path)
  return(out_path)
}


#' Perform Global Test Quantitative Enrichment Analysis (QEA) using MetaboAnalystR
#'
#' This function implements a self-contained QEA workflow using the Global Test
#' (QEA) implementation in MetaboAnalystR. It follows the recommended pipeline:
#' InitDataObjects -> Read.TextData -> CrossReferencing -> ReplaceMin ->
#' PreparePrenormData -> Normalization -> SetCurrentMsetLib ->
#' CalculateGlobalTestScore -> PlotQEA.Overview (bar & dot).
#'
#' @param conc_data A data.frame with metabolites in rows. First column must be
#'   compound name (e.g. "Compound"); remaining columns are samples.
#' @param group_vec A vector of group labels (length == number of samples).
#' @param norm_method Normalization method to pass to MetaboAnalystR::Normalization
#'   (default "SumNorm"). "ProbNorm" is also supported.
#' @param replace_min_frac Fraction used by ReplaceMin (default 0.5: replace with 1/2 min).
#' @param lib Library to use for pathway QEA (default "smpdb_pathway").
#' @param area Numeric area argument for SetCurrentMsetLib (mammalian = 0).
#' @return A list with elements: mSetObj (MetaboAnalystR mSet object),
#'   pathway_results_df (data.frame with p.values and topology), plots (list of
#'   file paths to barplot and dotplot PNGs and any recorded plot objects).
#' @export
performGlobalTestQEA <- function(conc_data, group_vec, norm_method = "SumNorm",
                                 replace_min_frac = 0.5,
                                 lib = "smpdb_pathway", area = 0,
                                 output_dir = NULL, save_outputs = TRUE) {
  # Basic checks
  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    stop("MetaboAnalystR is required. Install via devtools::install_github('xia-lab/MetaboAnalystR')")
  }

  if (!is.data.frame(conc_data)) stop("conc_data must be a data.frame with compounds in rows and samples in columns")
  if (ncol(conc_data) < 3) stop("conc_data must have at least one compound column and two samples")

  # Ensure first column is compound name
  first_col <- names(conc_data)[1]
  if (!tolower(first_col) %in% c("compound", "name", "metabolite", "feature")) {
    names(conc_data)[1] <- "Compound"
  }

  sample_names <- names(conc_data)[-1]
  if (length(group_vec) != length(sample_names)) stop("Length of group_vec must match number of samples (columns in conc_data excluding the first column)")

  # Report missingness
  missing_pct <- mean(is.na(as.matrix(conc_data[,-1]))) * 100
  if (missing_pct > 30) warning(sprintf("High missingness detected: %.1f%% of values are NA. ReplaceMin and imputation may be needed.", missing_pct))

  # Prepare temporary CSV for MetaboAnalyst Read.TextData
  tmpf <- tempfile(fileext = ".csv")
  utils::write.csv(conc_data, tmpf, row.names = FALSE, quote = FALSE)

  mSet <- NULL
  tryCatch({
    # 1. Init mSet for QEA
    mSet <- MetaboAnalystR::InitDataObjects("conc", "msetqea", FALSE)

    # 2. Read data (row-wise; concentration data)
    mSet <- MetaboAnalystR::Read.TextData(mSet, tmpf, "rowu", "conc")

    # Attach sample class labels (tries several common storage points)
    try({
      mSet$dataSet$cls <- as.character(group_vec)
    }, silent = TRUE)

    # 3. Cross reference compound names
    mSet <- MetaboAnalystR::CrossReferencing(mSet, "name")
    mSet <- MetaboAnalystR::CreateMappingResultTable(mSet)

    # 4. Sanity check and missing value replacement
    mSet <- MetaboAnalystR::SanityCheckData(mSet)
    mSet <- MetaboAnalystR::ReplaceMin(mSet, value = replace_min_frac)

    # 5. Normalization
    mSet <- MetaboAnalystR::PreparePrenormData(mSet)
    norm_succeeded <- FALSE
    try({
      mSet <- MetaboAnalystR::Normalization(mSet, norm_method)
      norm_succeeded <- TRUE
    }, silent = TRUE)
    if (!norm_succeeded) {
      try({
        mSet <- MetaboAnalystR::Normalization(mSet, norm_method, "NULL", "NULL", "PIF_178")
        norm_succeeded <- TRUE
      }, silent = TRUE)
    }
    if (!norm_succeeded) warning("Normalization via MetaboAnalystR::Normalization failed; proceeding with raw (but ReplaceMin applied) data")

    # 6. Select pathway library
    mSet <- tryCatch({
      MetaboAnalystR::SetCurrentMsetLib(mSet, lib, area)
    }, error = function(e) {
      warning(sprintf("SetCurrentMsetLib('%s', %s) failed: %s. Trying 'pathway' library as fallback.", lib, area, e$message))
      MetaboAnalystR::SetCurrentMsetLib(mSet, "pathway", area)
    })

    # 7. Calculate Global Test score
    mSet <- MetaboAnalystR::CalculateGlobalTestScore(mSet)

  }, error = function(e) {
    stop("MetaboAnalystR QEA pipeline failed: ", e$message)
  })

  # 8. Extract results
  pathway_df <- NULL
  try({
    if (!is.null(mSet$analSet$qea.global) && is.data.frame(mSet$analSet$qea.global)) {
      pathway_df <- mSet$analSet$qea.global
    } else if (!is.null(mSet$analSet$qea) && is.list(mSet$analSet$qea)) {
      df_candidates <- purrr::keep(mSet$analSet$qea, is.data.frame)
      if (length(df_candidates) >= 1) pathway_df <- df_candidates[[1]]
    }

    if (is.null(pathway_df)) {
      files <- list.files(pattern = "qea.*(csv|txt)$", ignore.case = TRUE)
      if (length(files) > 0) {
        f <- files[which.max(file.info(files)$mtime)]
        pathway_df <- tryCatch(readr::read_csv(f), error = function(e) NULL)
      }
    }
  }, silent = TRUE)

  if (is.null(pathway_df)) {
    warning("Could not automatically extract QEA results from 'mSet'. Returning mSet object only.")
    return(list(mSetObj = mSet, pathway_results_df = NULL, plots = list(), saved_paths = list()))
  }

  # Normalize result column names and add FDR
  names(pathway_df) <- tolower(names(pathway_df))
  pcol <- intersect(c("p.value", "pvalue", "p_val", "p"), names(pathway_df))
  if (length(pcol) == 0) {
    warning("No p-value column found in QEA results; returning raw table")
  } else {
    pcol <- pcol[1]
    pathway_df$p.adj <- stats::p.adjust(pathway_df[[pcol]], method = "fdr")
  }

  topo_col <- intersect(c("topology", "topo", "topology_score"), names(pathway_df))
  if (length(topo_col) == 0) topo_col <- intersect(c("impact", "centrality"), names(pathway_df))
  topo_col <- if (length(topo_col) == 0) NULL else topo_col[1]

  out_df <- pathway_df
  if (!is.null(topo_col)) names(out_df)[names(out_df) == topo_col] <- "topology"
  if (!("pathway" %in% names(out_df))) {
    pname <- intersect(c("pathway", "name", "id"), names(out_df))
    if (length(pname) >= 1) names(out_df)[names(out_df) == pname[1]] <- "pathway"
  }

  # 9. Generate plots
  plots <- list()
  oldwd <- getwd()
  tmp_plot_dir <- tempfile("qea_plots_")
  dir.create(tmp_plot_dir, showWarnings = FALSE)
  on.exit(setwd(oldwd), add = TRUE)
  setwd(tmp_plot_dir)

  try({
    MetaboAnalystR::PlotQEA.Overview(mSet, "qea_bar_", "bar", "png", 72, width = NA)
    bar_files <- list.files(pattern = "qea_bar_.*\\.png$")
    if (length(bar_files) > 0) plots$barplot_path <- file.path(tmp_plot_dir, bar_files[which.max(file.info(bar_files)$mtime)])
  }, silent = TRUE)

  try({
    MetaboAnalystR::PlotQEA.Overview(mSet, "qea_dot_", "dot", "png", 72, width = NA)
    dot_files <- list.files(pattern = "qea_dot_.*\\.png$")
    if (length(dot_files) > 0) plots$dotplot_path <- file.path(tmp_plot_dir, dot_files[which.max(file.info(dot_files)$mtime)])
  }, silent = TRUE)

  if (is.null(plots$barplot_path)) {
    try({
      MetaboAnalystR::PlotQEA.Overview(mSet, "qea_", "bar", "png", 72, width = NA)
      bar_files <- list.files(pattern = "qea_.*bar.*\\.png$|qea_.*bar.*png$")
      if (length(bar_files) > 0) plots$barplot_path <- file.path(tmp_plot_dir, bar_files[which.max(file.info(bar_files)$mtime)])
    }, silent = TRUE)
  }

  if (is.null(plots$dotplot_path)) {
    try({
      MetaboAnalystR::PlotQEA.Overview(mSet, "qea_", "dot", "png", 72, width = NA)
      dot_files <- list.files(pattern = "qea_.*dot.*\\.png$|qea_.*dot.*png$")
      if (length(dot_files) > 0) plots$dotplot_path <- file.path(tmp_plot_dir, dot_files[which.max(file.info(dot_files)$mtime)])
    }, silent = TRUE)
  }

  # Save outputs if requested
  saved_paths <- list()
  if (isTRUE(save_outputs)) {
    if (!is.null(output_dir)) {
      out_dir <- file.path(output_dir, "qea")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      try({
        saved_results_csv <- file.path(out_dir, "qea_pathway_results.csv")
        readr::write_csv(out_df, saved_results_csv)
        saved_paths$results_csv <- saved_results_csv
      }, silent = TRUE)

      if (!is.null(plots$barplot_path) && file.exists(plots$barplot_path)) {
        dest <- file.path(out_dir, basename(plots$barplot_path))
        file.copy(plots$barplot_path, dest, overwrite = TRUE)
        saved_paths$barplot_path <- dest
      }
      if (!is.null(plots$dotplot_path) && file.exists(plots$dotplot_path)) {
        dest <- file.path(out_dir, basename(plots$dotplot_path))
        file.copy(plots$dotplot_path, dest, overwrite = TRUE)
        saved_paths$dotplot_path <- dest
      }

      try({
        saved_mset <- file.path(out_dir, "qea_mset.rds")
        saveRDS(mSet, saved_mset)
        saved_paths$mset_rds <- saved_mset
      }, silent = TRUE)

    } else {
      warning("save_outputs=TRUE but output_dir is NULL; no files were saved.")
    }
  }

  # Return main outputs and saved paths
  return(list(
    mSetObj = mSet,
    pathway_results_df = out_df,
    plots = plots,
    saved_paths = saved_paths
  ))
}


# -------------------------
# Example usage (synthetic data)
# -------------------------
#' ## Example
#' # Generate synthetic dataset: 100 metabolites x 12 samples, two groups
#' # set.seed(42)
#' # n_met <- 100; n_samp <- 12
#' # met_names <- paste0("M", sprintf("%03d", 1:n_met))
#' # groups <- rep(c("Control", "Treatment"), each = n_samp/2)
' # # Create signal for Treatment group for a subset of metabolites
#' # mat <- matrix(rnorm(n_met * n_samp, mean = 5, sd = 2), nrow = n_met)
#' # mat[sample(1:(n_met*n_samp), size = floor(0.15*n_met*n_samp))] <- NA # inject missingness
#' # mat[1:10, groups == "Treatment"] <- mat[1:10, groups == "Treatment"] + 1.5
#' # df <- data.frame(Compound = met_names, mat)
' # # names(df)[-1] <- paste0("S", 1:n_samp)
#' # res <- performGlobalTestQEA(df, groups, norm_method = "SumNorm")
#' # head(res$pathway_results_df)
#' # print(res$plots)

