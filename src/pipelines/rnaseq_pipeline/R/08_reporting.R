# R/08_reporting.R
# Utility functions for saving outputs and creating summaries

#' Save data frame to CSV file
#' @param data Data frame to save
#' @param file_path Output file path
#' @return The file path (for targets tracking)
save_csv <- function(data, file_path) {
  # Ensure directory exists
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)

  # Convert to data frame if needed
  if (!is.data.frame(data)) {
    if (is.matrix(data)) {
      data <- as.data.frame(data)
    } else {
      data <- as.data.frame(as.list(data))
    }
  }

  write.csv(data, file_path, row.names = FALSE)
  message("Saved: ", file_path)

  file_path
}


#' Create pipeline summary
#' @param filtering_stats Filtering statistics
#' @param annotation_stats Annotation statistics
#' @param de_summary DE analysis summary
#' @param pathway_results Pathway analysis results
#' @param params Pipeline parameters
#' @return Summary data frame
create_pipeline_summary <- function(filtering_stats,
                                     annotation_stats,
                                     de_summary,
                                     pathway_results,
                                     params) {

  summary_rows <- list()

  # Basic info
  summary_rows$organism <- data.frame(
    category = "Input",
    metric = "Organism",
    value = params$organism
  )

  summary_rows$design <- data.frame(
    category = "Input",
    metric = "Design formula",
    value = params$design_formula
  )

  # Filtering stats
  if (is.data.frame(filtering_stats)) {
    for (i in 1:nrow(filtering_stats)) {
      summary_rows[[paste0("filter_", i)]] <- data.frame(
        category = "Filtering",
        metric = filtering_stats$metric[i],
        value = as.character(filtering_stats$value[i])
      )
    }
  }

  # Annotation stats
  if (is.list(annotation_stats)) {
    summary_rows$annot_source <- data.frame(
      category = "Annotation",
      metric = "Source",
      value = annotation_stats$source
    )
    summary_rows$annot_total <- data.frame(
      category = "Annotation",
      metric = "Total genes",
      value = as.character(annotation_stats$total_genes)
    )
    summary_rows$annot_success <- data.frame(
      category = "Annotation",
      metric = "Annotated genes",
      value = as.character(annotation_stats$annotated_genes)
    )
    summary_rows$annot_rate <- data.frame(
      category = "Annotation",
      metric = "Success rate",
      value = paste0(round(annotation_stats$success_rate * 100, 1), "%")
    )
  }

  # DE summary
  if (is.data.frame(de_summary) && nrow(de_summary) > 0) {
    for (i in 1:nrow(de_summary)) {
      contrast <- de_summary$contrast[i]
      summary_rows[[paste0("de_", i, "_sig")]] <- data.frame(
        category = "Differential Expression",
        metric = paste0(contrast, " - significant"),
        value = as.character(de_summary$significant[i])
      )
      summary_rows[[paste0("de_", i, "_up")]] <- data.frame(
        category = "Differential Expression",
        metric = paste0(contrast, " - upregulated"),
        value = as.character(de_summary$up[i])
      )
      summary_rows[[paste0("de_", i, "_down")]] <- data.frame(
        category = "Differential Expression",
        metric = paste0(contrast, " - downregulated"),
        value = as.character(de_summary$down[i])
      )
    }
  }

  # Pathway summary
  if (is.list(pathway_results) && length(pathway_results) > 0) {
    for (contrast_name in names(pathway_results)) {
      contrast_res <- pathway_results[[contrast_name]]
      for (analysis_name in names(contrast_res)) {
        res_df <- contrast_res[[analysis_name]]
        if (!is.null(res_df) && nrow(res_df) > 0 && "padj" %in% colnames(res_df)) {
          n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
          summary_rows[[paste0("pathway_", contrast_name, "_", analysis_name)]] <- data.frame(
            category = "Pathway Analysis",
            metric = paste0(contrast_name, " - ", analysis_name),
            value = paste0(n_sig, " significant (padj < 0.05)")
          )
        }
      }
    }
  }

  # Combine all rows
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  summary_df
}


#' Format numbers for display in reports
#' @param x Numeric value
#' @param digits Number of decimal places
#' @return Formatted string
format_number <- function(x, digits = 2) {
  if (is.na(x)) return("NA")
  if (abs(x) >= 1e6) {
    return(paste0(round(x / 1e6, digits), "M"))
  } else if (abs(x) >= 1e3) {
    return(paste0(round(x / 1e3, digits), "K"))
  } else {
    return(round(x, digits))
  }
}


#' Format p-value for display
#' @param p P-value
#' @param digits Significant digits
#' @return Formatted string
format_pvalue <- function(p, digits = 3) {
  if (is.na(p)) return("NA")
  if (p < 10^(-digits)) {
    return(paste0("< ", 10^(-digits)))
  } else {
    return(signif(p, digits))
  }
}
