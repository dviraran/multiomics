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
