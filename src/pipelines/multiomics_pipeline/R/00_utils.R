# =============================================================================
# Utility Functions for Multi-Omics Pipeline
# =============================================================================

#' Load configuration file
load_config <- function(config_path = "config.yml") {
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  config <- yaml::read_yaml(config_path)
  config <- set_config_defaults(config)
  log_message("DEBUG_LOAD: Loaded config from ", config_path)
  return(config)
}

#' Set default values
set_config_defaults <- function(config) {
  if (is.null(config$global$sample_id_column)) {
    config$global$sample_id_column <- "sample_id"
  }
  if (is.null(config$harmonization$sample_mode)) {
    config$harmonization$sample_mode <- "intersection"
  }
  if (is.null(config$harmonization$duplicate_strategy)) {
    config$harmonization$duplicate_strategy <- "keep_max_mean"
  }
  if (is.null(config$output$output_dir)) {
    config$output$output_dir <- "outputs"
  }
  return(config)
}

#' Create output directories
create_output_dirs <- function(config) {
  output_dir <- config$output$output_dir
  dirs <- c(output_dir, file.path(output_dir, c("tables", "plots", "report", "logs")))
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  invisible(dirs)
}

#' Save table to CSV
save_table <- function(df, filename, config, subdir = "tables") {
  output_path <- file.path(config$output$output_dir, subdir, filename)
  readr::write_csv(df, output_path)
  log_message("Saved: ", output_path)
  return(output_path)
}

#' Save plot
save_plot <- function(plot, filename, config, width = NULL, height = NULL, subdir = "plots") {
  if (is.null(width)) width <- config$output$fig_width %||% 10
  if (is.null(height)) height <- config$output$fig_height %||% 8
  dpi <- config$output$plot_dpi %||% 300
  output_path <- file.path(config$output$output_dir, subdir, filename)
  ggplot2::ggsave(output_path, plot = plot, width = width, height = height, dpi = dpi)
  log_message("Saved: ", output_path)
  return(output_path)
}

#' Log message with timestamp
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(msg)
}

#' Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Sanitize sample names
sanitize_names <- function(names) {
  names <- trimws(names)
  names <- gsub("[^a-zA-Z0-9_.-]", "_", names)
  return(names)
}

#' Parse UniProt accession
parse_uniprot_accession <- function(ids) {
  uniprot_pattern <- "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?"
  extracted <- stringr::str_extract(ids, paste0("(sp|tr)\\|", uniprot_pattern, "\\|"))
  extracted <- stringr::str_replace_all(extracted, "(sp|tr)\\|", "")
  extracted <- stringr::str_replace_all(extracted, "\\|$", "")
  no_match <- is.na(extracted)
  if (any(no_match)) {
    direct_match <- stringr::str_extract(ids[no_match], uniprot_pattern)
    extracted[no_match] <- direct_match
  }
  return(extracted)
}

#' Strip Ensembl version suffix
strip_ensembl_version <- function(ids) {
  gsub("\\.[0-9]+$", "", ids)
}

#' Read GMT file
read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)
  sets <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      set_name <- parts[1]
      members <- parts[3:length(parts)]
      members <- members[members != ""]
      sets[[set_name]] <- members
    }
  }
  return(sets)
}

#' Calculate CV
calc_cv <- function(x, na.rm = TRUE) {
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}

#' Fisher's method for combining p-values
fisher_combine_pvalues <- function(pvals) {
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) == 0) {
    return(NA)
  }
  chi_sq <- -2 * sum(log(pvals))
  df <- 2 * length(pvals)
  combined_p <- pchisq(chi_sq, df, lower.tail = FALSE)
  return(combined_p)
}

#' Stouffer's method for combining p-values
stouffer_combine_pvalues <- function(pvals, weights = NULL) {
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals < 1]
  if (length(pvals) == 0) {
    return(NA)
  }
  if (is.null(weights)) weights <- rep(1, length(pvals))
  z_scores <- qnorm(1 - pvals)
  combined_z <- sum(weights * z_scores) / sqrt(sum(weights^2))
  combined_p <- 1 - pnorm(combined_z)
  return(combined_p)
}
