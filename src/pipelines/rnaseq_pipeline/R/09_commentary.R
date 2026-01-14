# R/09_commentary.R
# Functions for AI-generated figure commentary
# Supports Claude vision API and data-driven fallback

#' Build table of all figures with metadata
#' @param qc_plots List of QC plot paths from generate_qc_plots()
#' @param de_plots List of DE plot paths from generate_de_plots()
#' @param pathway_plots List of pathway plot paths from generate_pathway_plots()
#' @param batch_analysis Batch analysis results (may contain plot paths)
#' @param gsva_results GSVA analysis results (may contain plot paths)
#' @param wgcna_results WGCNA analysis results (may contain plot paths)
#' @param deconvolution_results Deconvolution results (may contain plot paths)
#' @param config Pipeline configuration
#' @return Data frame with figure metadata
build_figures_table <- function(qc_plots = NULL,
                                 de_plots = NULL,
                                 pathway_plots = NULL,
                                 batch_analysis = NULL,
                                 gsva_results = NULL,
                                 wgcna_results = NULL,
                                 deconvolution_results = NULL,
                                 config = NULL) {

  figures <- list()

  # QC plots
  if (!is.null(qc_plots)) {
    if (!is.null(qc_plots$library_sizes)) {
      figures$library_sizes <- data.frame(
        figure_id = "library_sizes",
        filepath = qc_plots$library_sizes,
        plot_type = "barplot",
        section = "QC",
        title = "Library Sizes",
        description = "Total sequencing depth per sample",
        x_axis = "Sample",
        y_axis = "Total Counts (millions)",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(qc_plots$detected_genes)) {
      figures$detected_genes <- data.frame(
        figure_id = "detected_genes",
        filepath = qc_plots$detected_genes,
        plot_type = "barplot",
        section = "QC",
        title = "Detected Genes per Sample",
        description = "Number of genes with non-zero counts per sample",
        x_axis = "Sample",
        y_axis = "Number of Genes",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(qc_plots$correlation_heatmap)) {
      figures$correlation_heatmap <- data.frame(
        figure_id = "sample_correlation",
        filepath = qc_plots$correlation_heatmap,
        plot_type = "heatmap",
        section = "QC",
        title = "Sample-to-Sample Correlation",
        description = "Pearson correlation matrix between samples based on normalized expression",
        x_axis = "Samples",
        y_axis = "Samples",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(qc_plots$pca)) {
      figures$pca <- data.frame(
        figure_id = "pca_plot",
        filepath = qc_plots$pca,
        plot_type = "scatter",
        section = "QC",
        title = "PCA Plot",
        description = "Principal component analysis of samples based on top variable genes",
        x_axis = "PC1",
        y_axis = "PC2",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(qc_plots$pca_scree)) {
      figures$pca_scree <- data.frame(
        figure_id = "pca_scree",
        filepath = qc_plots$pca_scree,
        plot_type = "barplot",
        section = "QC",
        title = "PCA Variance Explained",
        description = "Variance explained by each principal component",
        x_axis = "Principal Component",
        y_axis = "% Variance",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  # DE plots (volcano, MA, heatmap per contrast)
  if (!is.null(de_plots)) {
    for (plot_name in names(de_plots)) {
      filepath <- de_plots[[plot_name]]

      if (grepl("^volcano_", plot_name)) {
        contrast <- sub("^volcano_", "", plot_name)
        figures[[plot_name]] <- data.frame(
          figure_id = plot_name,
          filepath = filepath,
          plot_type = "volcano",
          section = "Differential Expression",
          title = paste("Volcano Plot:", contrast),
          description = "Statistical significance vs fold change for all genes",
          x_axis = "log2 Fold Change",
          y_axis = "-log10(adjusted p-value)",
          contrast = contrast,
          stringsAsFactors = FALSE
        )
      } else if (grepl("^ma_", plot_name)) {
        contrast <- sub("^ma_", "", plot_name)
        figures[[plot_name]] <- data.frame(
          figure_id = plot_name,
          filepath = filepath,
          plot_type = "ma_plot",
          section = "Differential Expression",
          title = paste("MA Plot:", contrast),
          description = "Fold change vs mean expression for all genes",
          x_axis = "log10(Mean Expression)",
          y_axis = "log2 Fold Change",
          contrast = contrast,
          stringsAsFactors = FALSE
        )
      } else if (grepl("^heatmap_", plot_name)) {
        contrast <- sub("^heatmap_top_", "", plot_name)
        figures[[plot_name]] <- data.frame(
          figure_id = plot_name,
          filepath = filepath,
          plot_type = "heatmap",
          section = "Differential Expression",
          title = paste("Top DE Genes Heatmap:", contrast),
          description = "Expression patterns of top differentially expressed genes",
          x_axis = "Samples",
          y_axis = "Genes",
          contrast = contrast,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Pathway plots
  if (!is.null(pathway_plots)) {
    for (plot_name in names(pathway_plots)) {
      filepath <- pathway_plots[[plot_name]]
      # Extract contrast and analysis type from name
      # Format: pathway_[contrast]_[analysis].png
      parts <- strsplit(sub("^pathway_", "", plot_name), "_(?=[^_]+$)", perl = TRUE)[[1]]
      if (length(parts) >= 2) {
        contrast <- parts[1]
        analysis <- parts[2]
      } else {
        contrast <- plot_name
        analysis <- "enrichment"
      }

      figures[[plot_name]] <- data.frame(
        figure_id = plot_name,
        filepath = filepath,
        plot_type = "dotplot",
        section = "Pathway Analysis",
        title = paste("Pathway Enrichment:", contrast),
        description = paste("Enriched pathways from", analysis, "analysis"),
        x_axis = "Enrichment Score",
        y_axis = "Pathway",
        contrast = contrast,
        stringsAsFactors = FALSE
      )
    }
  }

  # Batch analysis plots
  if (!is.null(batch_analysis) && is.list(batch_analysis)) {
    if (!is.null(batch_analysis$plots)) {
      for (plot_name in names(batch_analysis$plots)) {
        filepath <- batch_analysis$plots[[plot_name]]
        if (!is.null(filepath) && is.character(filepath) && file.exists(filepath)) {
          figures[[paste0("batch_", plot_name)]] <- data.frame(
            figure_id = paste0("batch_", plot_name),
            filepath = filepath,
            plot_type = "batch",
            section = "Batch Effects",
            title = paste("Batch Analysis:", gsub("_", " ", plot_name)),
            description = "Visualization of batch effects in the data",
            x_axis = NA_character_,
            y_axis = NA_character_,
            contrast = NA_character_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  # GSVA results plots
  if (!is.null(gsva_results) && is.list(gsva_results)) {
    if (!is.null(gsva_results$plots)) {
      for (plot_name in names(gsva_results$plots)) {
        filepath <- gsva_results$plots[[plot_name]]
        if (!is.null(filepath) && is.character(filepath) && file.exists(filepath)) {
          figures[[paste0("gsva_", plot_name)]] <- data.frame(
            figure_id = paste0("gsva_", plot_name),
            filepath = filepath,
            plot_type = "heatmap",
            section = "GSVA Analysis",
            title = paste("GSVA:", gsub("_", " ", plot_name)),
            description = "Gene set variation analysis results",
            x_axis = "Samples",
            y_axis = "Gene Sets",
            contrast = NA_character_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  # WGCNA results plots
  if (!is.null(wgcna_results) && is.list(wgcna_results)) {
    if (!is.null(wgcna_results$plots)) {
      for (plot_name in names(wgcna_results$plots)) {
        filepath <- wgcna_results$plots[[plot_name]]
        if (!is.null(filepath) && is.character(filepath) && file.exists(filepath)) {
          figures[[paste0("wgcna_", plot_name)]] <- data.frame(
            figure_id = paste0("wgcna_", plot_name),
            filepath = filepath,
            plot_type = "network",
            section = "WGCNA Co-Expression",
            title = paste("WGCNA:", gsub("_", " ", plot_name)),
            description = "Weighted gene co-expression network analysis",
            x_axis = NA_character_,
            y_axis = NA_character_,
            contrast = NA_character_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  # Deconvolution results plots
  if (!is.null(deconvolution_results) && is.list(deconvolution_results)) {
    if (!is.null(deconvolution_results$plots)) {
      for (plot_name in names(deconvolution_results$plots)) {
        filepath <- deconvolution_results$plots[[plot_name]]
        if (!is.null(filepath) && is.character(filepath) && file.exists(filepath)) {
          figures[[paste0("deconv_", plot_name)]] <- data.frame(
            figure_id = paste0("deconv_", plot_name),
            filepath = filepath,
            plot_type = "barplot",
            section = "Cell Type Deconvolution",
            title = paste("Deconvolution:", gsub("_", " ", plot_name)),
            description = "Cell type composition estimates",
            x_axis = "Samples",
            y_axis = "Cell Type Proportion",
            contrast = NA_character_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  # Combine all figures
  if (length(figures) == 0) {
    message("No figures found to document")
    return(data.frame(
      figure_id = character(),
      filepath = character(),
      plot_type = character(),
      section = character(),
      title = character(),
      description = character(),
      x_axis = character(),
      y_axis = character(),
      contrast = character(),
      stringsAsFactors = FALSE
    ))
  }

  figures_tbl <- do.call(rbind, figures)
  rownames(figures_tbl) <- NULL

  # Add experiment context
  if (!is.null(config)) {
    figures_tbl$organism <- config$organism
    figures_tbl$group_col <- config$group_col
    figures_tbl$design_formula <- config$design_formula
  }

 message("Built figures table with ", nrow(figures_tbl), " figures")

  figures_tbl
}


#' Generate commentary for all figures
#' @param figures_tbl Figures metadata table
#' @param config Pipeline configuration
#' @param qc_metrics QC metrics for data-driven fallback
#' @param pca_result PCA results for data-driven fallback
#' @param sample_correlation Correlation matrix for data-driven fallback
#' @param de_summary DE summary for data-driven fallback
#' @param outlier_detection Outlier flags for data-driven fallback
#' @param batch_analysis Batch analysis results for data-driven fallback
#' @param gsva_results GSVA results for data-driven fallback
#' @param wgcna_results WGCNA results for data-driven fallback
#' @param deconvolution_results Deconvolution results for data-driven fallback
#' @param output_dir Directory for commentary outputs
#' @return Data frame with commentary for each figure
generate_all_commentary <- function(figures_tbl,
                                     config,
                                     qc_metrics = NULL,
                                     pca_result = NULL,
                                     sample_correlation = NULL,
                                     de_summary = NULL,
                                     outlier_detection = NULL,
                                     batch_analysis = NULL,
                                     gsva_results = NULL,
                                     wgcna_results = NULL,
                                     deconvolution_results = NULL,
                                     output_dir = "outputs/commentary") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Check if commentary is enabled
  commentary_enabled <- isTRUE(config$commentary_enabled)
  backend <- config$commentary_backend %||% "none"

  if (!commentary_enabled) {
    message("Commentary disabled in config. Using placeholder text.")
    backend <- "none"
  }

  commentary_list <- list()

  for (i in seq_len(nrow(figures_tbl))) {
    fig <- figures_tbl[i, ]

    message("Generating commentary for: ", fig$figure_id)

    # Check if file exists
    if (!file.exists(fig$filepath)) {
      message("  Warning: Figure file not found: ", fig$filepath)
      commentary_list[[fig$figure_id]] <- create_placeholder_commentary(
        fig$figure_id,
        reason = "Figure file not found"
      )
      next
    }

    # Build context for the figure
    context <- build_figure_context(
      fig = fig,
      config = config,
      qc_metrics = qc_metrics,
      pca_result = pca_result,
      sample_correlation = sample_correlation,
      de_summary = de_summary,
      outlier_detection = outlier_detection
    )

    # Generate commentary based on backend
    if (backend == "claude") {
      commentary <- run_claude_commentary(
        image_path = fig$filepath,
        figure_id = fig$figure_id,
        context = context,
        config = config,
        output_dir = output_dir
      )
    } else if (backend == "openai") {
      commentary <- run_openai_commentary(
        image_path = fig$filepath,
        figure_id = fig$figure_id,
        context = context,
        config = config,
        output_dir = output_dir
      )
    } else {
      # Data-driven fallback (no LLM)
      commentary <- generate_fallback_commentary(
        fig = fig,
        context = context,
        qc_metrics = qc_metrics,
        pca_result = pca_result,
        sample_correlation = sample_correlation,
        de_summary = de_summary,
        outlier_detection = outlier_detection
      )
    }

    commentary_list[[fig$figure_id]] <- commentary

    # Save individual JSON
    json_path <- file.path(output_dir, paste0(fig$figure_id, "_commentary.json"))
    jsonlite::write_json(commentary, json_path, auto_unbox = TRUE, pretty = TRUE)
  }

  # Combine into data frame
  commentary_tbl <- do.call(rbind, lapply(names(commentary_list), function(id) {
    c <- commentary_list[[id]]
    data.frame(
      figure_id = id,
      title = c$title %||% NA_character_,
      what_is_this = c$what_is_this %||% NA_character_,
      observations = paste(c$observations, collapse = " | ") %||% NA_character_,
      issues_checks = paste(c$issues_checks, collapse = " | ") %||% NA_character_,
      next_steps = paste(c$next_steps, collapse = " | ") %||% NA_character_,
      confidence = c$confidence %||% NA_character_,
      backend = c$backend %||% "unknown",
      generated_at = c$generated_at %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))

  # Save combined outputs
  csv_path <- file.path(output_dir, "commentary_summary.csv")
  write.csv(commentary_tbl, csv_path, row.names = FALSE)

  json_combined_path <- file.path(output_dir, "commentary_all.json")
  jsonlite::write_json(commentary_list, json_combined_path, auto_unbox = TRUE, pretty = TRUE)

  message("Saved commentary to ", output_dir)
  message("  - Individual JSONs: ", length(commentary_list), " files")
  message("  - Combined CSV: ", csv_path)
  message("  - Combined JSON: ", json_combined_path)

  # Return as list for easier access in Rmd
  attr(commentary_tbl, "commentary_list") <- commentary_list

  commentary_tbl
}


#' Build context information for a figure
#' @param fig Figure row from figures_tbl
#' @param config Pipeline config
#' @param ... Additional data objects
#' @return List with context information
build_figure_context <- function(fig,
                                  config,
                                  qc_metrics = NULL,
                                  pca_result = NULL,
                                  sample_correlation = NULL,
                                  de_summary = NULL,
                                  outlier_detection = NULL) {

  context <- list(
    figure_id = fig$figure_id,
    plot_type = fig$plot_type,
    section = fig$section,
    title = fig$title,
    description = fig$description,
    x_axis = fig$x_axis,
    y_axis = fig$y_axis,
    contrast = fig$contrast
  )

  # Add experiment summary
  if (!is.null(config)) {
    context$organism <- config$organism
    context$group_col <- config$group_col
    context$design_formula <- config$design_formula
  }

  # Add sample count
  if (!is.null(qc_metrics)) {
    context$n_samples <- nrow(qc_metrics)
    context$sample_names <- qc_metrics$sample
  }

  # Add group information
  if (!is.null(qc_metrics) && !is.null(config$group_col)) {
    # Groups would need to be passed in from metadata
  }

  # Add PCA variance info
 if (!is.null(pca_result)) {
    context$pca_var_pc1 <- round(pca_result$var_explained[1], 1)
    context$pca_var_pc2 <- round(pca_result$var_explained[2], 1)
    context$pca_n_genes <- pca_result$n_genes_used
  }

  # Add outlier info
  if (!is.null(outlier_detection)) {
    n_outliers <- sum(outlier_detection$is_outlier)
    context$n_outliers <- n_outliers
    if (n_outliers > 0) {
      context$outlier_samples <- outlier_detection$sample[outlier_detection$is_outlier]
    }
  }

  # Add DE summary for relevant plots
  if (!is.null(de_summary) && !is.na(fig$contrast)) {
    contrast_row <- de_summary[de_summary$contrast == fig$contrast, ]
    if (nrow(contrast_row) > 0) {
      context$de_significant <- contrast_row$significant[1]
      context$de_up <- contrast_row$up[1]
      context$de_down <- contrast_row$down[1]
      context$de_total <- contrast_row$total_genes[1]
    }
  }

  context
}


#' Run Claude vision API for commentary
#' @param image_path Path to the image file
#' @param figure_id Figure identifier
#' @param context Context information
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Commentary list
run_claude_commentary <- function(image_path,
                                   figure_id,
                                   context,
                                   config,
                                   output_dir) {

  # Build paths
  context_json <- file.path(output_dir, paste0(figure_id, "_context.json"))
  output_json <- file.path(output_dir, paste0(figure_id, "_commentary.json"))

  # Save context for the script
  jsonlite::write_json(context, context_json, auto_unbox = TRUE, pretty = TRUE)

  # Get script path
  script_path <- here::here("scripts", "figure_commentary_claude.py")

  if (!file.exists(script_path)) {
    message("  Claude commentary script not found: ", script_path)
    return(create_placeholder_commentary(figure_id, reason = "Script not found"))
  }

  # Build command
  model <- config$commentary_model %||% "claude-sonnet-4-20250514"
  max_tokens <- config$commentary_max_tokens %||% 1500

  cmd <- sprintf(
    "python3 '%s' --image '%s' --figure_id '%s' --context_json '%s' --out_json '%s' --model '%s' --max_tokens %d",
    script_path,
    image_path,
    figure_id,
    context_json,
    output_json,
    model,
    max_tokens
  )

  # Execute with timeout and retry
  max_retries <- config$commentary_max_retries %||% 2
  retry_delay <- config$commentary_retry_delay %||% 5

  for (attempt in seq_len(max_retries + 1)) {
    result <- tryCatch({
      system(cmd, intern = TRUE, ignore.stderr = FALSE, timeout = 120)

      # Check if output was created
      if (file.exists(output_json)) {
        commentary <- jsonlite::read_json(output_json)
        commentary$backend <- "claude"
        commentary$generated_at <- as.character(Sys.time())
        return(commentary)
      } else {
        stop("Output JSON not created")
      }
    }, error = function(e) {
      if (attempt <= max_retries) {
        message("  Attempt ", attempt, " failed: ", e$message, ". Retrying in ", retry_delay, "s...")
        Sys.sleep(retry_delay)
        return(NULL)
      } else {
        message("  Claude API call failed after ", max_retries + 1, " attempts: ", e$message)
        return(create_placeholder_commentary(figure_id, reason = paste("API error:", e$message)))
      }
    })

    if (!is.null(result) && is.list(result)) {
      return(result)
    }
  }

  create_placeholder_commentary(figure_id, reason = "API call failed")
}


#' Run OpenAI GPT-4 Vision API for commentary
#' @param image_path Path to the image file
#' @param figure_id Figure identifier
#' @param context Context information
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Commentary list
run_openai_commentary <- function(image_path,
                                   figure_id,
                                   context,
                                   config,
                                   output_dir) {

  # Build paths
  context_json <- file.path(output_dir, paste0(figure_id, "_context.json"))
  output_json <- file.path(output_dir, paste0(figure_id, "_commentary.json"))

  # Save context for the script
  jsonlite::write_json(context, context_json, auto_unbox = TRUE, pretty = TRUE)

  # Get script path
  script_path <- here::here("scripts", "figure_commentary_openai.py")

  if (!file.exists(script_path)) {
    message("  OpenAI commentary script not found: ", script_path)
    return(create_placeholder_commentary(figure_id, reason = "Script not found"))
  }

  # Build command - use openai_model if specified, otherwise default to gpt-4o
  model <- config$commentary_openai_model %||% "gpt-4o"
  max_tokens <- config$commentary_max_tokens %||% 1500

  cmd <- sprintf(
    "python3 '%s' --image '%s' --figure_id '%s' --context_json '%s' --out_json '%s' --model '%s' --max_tokens %d",
    script_path,
    image_path,
    figure_id,
    context_json,
    output_json,
    model,
    max_tokens
  )

  # Execute with timeout and retry
  max_retries <- config$commentary_max_retries %||% 2
  retry_delay <- config$commentary_retry_delay %||% 5

  for (attempt in seq_len(max_retries + 1)) {
    result <- tryCatch({
      system(cmd, intern = TRUE, ignore.stderr = FALSE, timeout = 120)

      # Check if output was created
      if (file.exists(output_json)) {
        commentary <- jsonlite::read_json(output_json)
        commentary$backend <- "openai"
        commentary$generated_at <- as.character(Sys.time())
        return(commentary)
      } else {
        stop("Output JSON not created")
      }
    }, error = function(e) {
      if (attempt <= max_retries) {
        message("  Attempt ", attempt, " failed: ", e$message, ". Retrying in ", retry_delay, "s...")
        Sys.sleep(retry_delay)
        return(NULL)
      } else {
        message("  OpenAI API call failed after ", max_retries + 1, " attempts: ", e$message)
        return(create_placeholder_commentary(figure_id, reason = paste("API error:", e$message)))
      }
    })

    if (!is.null(result) && is.list(result)) {
      return(result)
    }
  }

  create_placeholder_commentary(figure_id, reason = "API call failed")
}


#' Create placeholder commentary when LLM is unavailable
#' @param figure_id Figure identifier
#' @param reason Reason for placeholder
#' @return Commentary list with placeholder text
create_placeholder_commentary <- function(figure_id, reason = "Commentary disabled") {
  list(
    figure_id = figure_id,
    title = "Commentary Unavailable",
    what_is_this = paste("Automated commentary is not available.", reason),
    observations = list("Please refer to the figure and its caption for interpretation."),
    issues_checks = list(),
    next_steps = list(),
    confidence = "none",
    limitations = reason,
    backend = "placeholder",
    generated_at = as.character(Sys.time())
  )
}


#' Generate data-driven fallback commentary (no LLM)
#' @param fig Figure metadata row
#' @param context Figure context
#' @param ... Data objects for analysis
#' @return Commentary list
generate_fallback_commentary <- function(fig,
                                          context,
                                          qc_metrics = NULL,
                                          pca_result = NULL,
                                          sample_correlation = NULL,
                                          de_summary = NULL,
                                          outlier_detection = NULL) {

  commentary <- list(
    figure_id = fig$figure_id,
    title = fig$title,
    backend = "data_driven",
    generated_at = as.character(Sys.time()),
    confidence = "medium"
  )

  # Route to specific fallback function based on plot type
  if (fig$figure_id == "library_sizes" && !is.null(qc_metrics)) {
    commentary <- c(commentary, fallback_library_sizes(qc_metrics))
  } else if (fig$figure_id == "detected_genes" && !is.null(qc_metrics)) {
    commentary <- c(commentary, fallback_detected_genes(qc_metrics))
  } else if (fig$figure_id == "sample_correlation" && !is.null(sample_correlation)) {
    commentary <- c(commentary, fallback_correlation(sample_correlation, outlier_detection))
  } else if (fig$figure_id == "pca_plot" && !is.null(pca_result)) {
    commentary <- c(commentary, fallback_pca(pca_result, outlier_detection, context))
  } else if (fig$figure_id == "pca_scree" && !is.null(pca_result)) {
    commentary <- c(commentary, fallback_pca_scree(pca_result))
  } else if (grepl("^volcano_", fig$figure_id) && !is.null(de_summary)) {
    commentary <- c(commentary, fallback_volcano(de_summary, fig$contrast))
  } else if (grepl("^ma_", fig$figure_id) && !is.null(de_summary)) {
    commentary <- c(commentary, fallback_ma_plot(de_summary, fig$contrast))
  } else if (grepl("^heatmap_", fig$figure_id)) {
    commentary <- c(commentary, fallback_heatmap(de_summary, fig$contrast))
  } else if (grepl("^pathway_", fig$figure_id)) {
    commentary <- c(commentary, fallback_pathway())
  } else {
    # Generic fallback
    commentary$what_is_this <- fig$description
    commentary$observations <- list("Data-driven analysis not available for this plot type.")
    commentary$issues_checks <- list()
    commentary$next_steps <- list()
    commentary$confidence <- "low"
  }

  commentary
}


#' Fallback commentary for library sizes plot
#' @param qc_metrics QC metrics data frame
#' @return List with commentary fields
fallback_library_sizes <- function(qc_metrics) {
  total_counts <- qc_metrics$total_counts

  # Calculate statistics
  median_depth <- median(total_counts)
  min_depth <- min(total_counts)
  max_depth <- max(total_counts)
  cv <- sd(total_counts) / mean(total_counts)

  # Identify low samples (below 25th percentile - 1.5*IQR)
  q1 <- quantile(total_counts, 0.25)
  iqr <- IQR(total_counts)
  low_threshold <- q1 - 1.5 * iqr
  low_samples <- qc_metrics$sample[total_counts < max(low_threshold, median_depth * 0.5)]

  observations <- list(
    sprintf("Median library size: %.1fM reads", median_depth / 1e6),
    sprintf("Range: %.1fM to %.1fM reads", min_depth / 1e6, max_depth / 1e6),
    sprintf("Coefficient of variation: %.1f%%", cv * 100)
  )

  issues <- list()
  next_steps <- list()

  if (length(low_samples) > 0) {
    issues <- c(issues, sprintf("Samples with notably low depth: %s", paste(low_samples, collapse = ", ")))
    next_steps <- c(next_steps, "Consider investigating low-depth samples for technical issues")
  }

  if (cv > 0.5) {
    issues <- c(issues, "High variability in library sizes detected")
    next_steps <- c(next_steps, "Verify normalization adequately corrects for depth differences")
  }

  if (median_depth < 1e6) {
    issues <- c(issues, "Overall sequencing depth appears low")
    next_steps <- c(next_steps, "Consider whether depth is sufficient for detecting lowly expressed genes")
  }

  if (length(issues) == 0) {
    issues <- list("No major concerns identified")
  }

  list(
    what_is_this = "This barplot shows the total number of sequencing reads (library size) for each sample. Library size affects statistical power and should be reasonably consistent across samples.",
    observations = observations,
    issues_checks = issues,
    next_steps = next_steps
  )
}


#' Fallback commentary for detected genes plot
#' @param qc_metrics QC metrics data frame
#' @return List with commentary fields
fallback_detected_genes <- function(qc_metrics) {
  detected <- qc_metrics$detected_genes

  median_detected <- median(detected)
  min_detected <- min(detected)
  max_detected <- max(detected)
  cv <- sd(detected) / mean(detected)

  # Identify samples with notably fewer detected genes
  q1 <- quantile(detected, 0.25)
  iqr <- IQR(detected)
  low_threshold <- q1 - 1.5 * iqr
  low_samples <- qc_metrics$sample[detected < max(low_threshold, median_detected * 0.7)]

  observations <- list(
    sprintf("Median genes detected: %s", format(median_detected, big.mark = ",")),
    sprintf("Range: %s to %s genes", format(min_detected, big.mark = ","), format(max_detected, big.mark = ",")),
    sprintf("Detection rate variability (CV): %.1f%%", cv * 100)
  )

  issues <- list()
  next_steps <- list()

  if (length(low_samples) > 0) {
    issues <- c(issues, sprintf("Samples with fewer detected genes: %s", paste(low_samples, collapse = ", ")))
    next_steps <- c(next_steps, "Investigate if these samples have quality issues or represent distinct biology")
  }

  if (cv > 0.15) {
    issues <- c(issues, "Notable variation in gene detection across samples")
  }

  if (median_detected < 10000) {
    issues <- c(issues, "Overall gene detection appears relatively low")
    next_steps <- c(next_steps, "Check if filtering parameters are appropriate for this dataset")
  }

  if (length(issues) == 0) {
    issues <- list("Gene detection is consistent across samples")
  }

  list(
    what_is_this = "This barplot shows the number of genes detected (with at least one count) per sample. Samples with fewer detected genes may indicate technical issues or distinct biological states.",
    observations = observations,
    issues_checks = issues,
    next_steps = next_steps
  )
}


#' Fallback commentary for correlation heatmap
#' @param sample_correlation Correlation matrix
#' @param outlier_detection Outlier flags
#' @return List with commentary fields
fallback_correlation <- function(sample_correlation, outlier_detection = NULL) {
  # Get off-diagonal correlations
  cor_mat <- sample_correlation
  diag(cor_mat) <- NA

  # Calculate per-sample median correlation
  median_cors <- apply(cor_mat, 1, median, na.rm = TRUE)
  overall_median <- median(median_cors)
  min_cor <- min(median_cors)

  # Identify samples with low correlation
  low_cor_threshold <- overall_median - 2 * sd(median_cors)
  low_cor_samples <- names(median_cors)[median_cors < low_cor_threshold]

  observations <- list(
    sprintf("Overall median sample correlation: %.3f", overall_median),
    sprintf("Lowest median correlation: %.3f", min_cor)
  )

  if (overall_median > 0.9) {
    observations <- c(observations, "High correlation suggests consistent expression profiles")
  } else if (overall_median < 0.7) {
    observations <- c(observations, "Lower correlation may indicate batch effects or biological heterogeneity")
  }

  issues <- list()
  next_steps <- list()

  if (length(low_cor_samples) > 0) {
    issues <- c(issues, sprintf("Samples with low correlation: %s", paste(low_cor_samples, collapse = ", ")))
    next_steps <- c(next_steps, "Consider whether these samples should be excluded or investigated further")
  }

  if (!is.null(outlier_detection)) {
    cor_outliers <- outlier_detection$sample[outlier_detection$correlation_outlier]
    if (length(cor_outliers) > 0) {
      issues <- c(issues, sprintf("Correlation-based outliers flagged: %s", paste(cor_outliers, collapse = ", ")))
    }
  }

  if (length(issues) == 0) {
    issues <- list("Sample correlations appear consistent")
  }

  list(
    what_is_this = "This heatmap shows Pearson correlation coefficients between all pairs of samples based on normalized expression values. High correlation (red) indicates similar expression profiles.",
    observations = observations,
    issues_checks = issues,
    next_steps = next_steps
  )
}


#' Fallback commentary for PCA plot
#' @param pca_result PCA results list
#' @param outlier_detection Outlier flags
#' @param context Figure context with group info
#' @return List with commentary fields
fallback_pca <- function(pca_result, outlier_detection = NULL, context = NULL) {
  var_exp <- pca_result$var_explained
  pca_data <- pca_result$pca_data

  observations <- list(
    sprintf("PC1 explains %.1f%% of variance", var_exp[1]),
    sprintf("PC2 explains %.1f%% of variance", var_exp[2]),
    sprintf("Total variance in PC1-PC2: %.1f%%", var_exp[1] + var_exp[2])
  )

  # Check for group separation if group column exists
  group_col <- context$group_col
  if (!is.null(group_col) && group_col %in% colnames(pca_data)) {
    groups <- unique(pca_data[[group_col]])
    if (length(groups) > 1) {
      # Calculate group centroids
      centroids <- aggregate(cbind(PC1, PC2) ~ pca_data[[group_col]], data = pca_data, mean)

      # Calculate between-group distances
      if (nrow(centroids) >= 2) {
        distances <- dist(centroids[, c("PC1", "PC2")])
        mean_dist <- mean(distances)

        # Compare to within-group spread
        within_spread <- mean(sapply(groups, function(g) {
          group_data <- pca_data[pca_data[[group_col]] == g, ]
          if (nrow(group_data) > 1) {
            mean(dist(group_data[, c("PC1", "PC2")]))
          } else {
            0
          }
        }))

        if (mean_dist > within_spread * 1.5) {
          observations <- c(observations, "Groups show clear separation in PCA space")
        } else if (mean_dist > within_spread) {
          observations <- c(observations, "Groups show moderate separation in PCA space")
        } else {
          observations <- c(observations, "Groups show overlap in PCA space")
        }
      }
    }
  }

  issues <- list()
  next_steps <- list()

  # Check for outliers
  if (!is.null(outlier_detection)) {
    pca_outliers <- outlier_detection$sample[outlier_detection$pca_outlier]
    if (length(pca_outliers) > 0) {
      issues <- c(issues, sprintf("PCA-based outliers: %s", paste(pca_outliers, collapse = ", ")))
      next_steps <- c(next_steps, "Investigate outlier samples for technical or biological explanations")
    }
  }

  # Check variance explained
  if (var_exp[1] > 50) {
    issues <- c(issues, "PC1 explains >50% variance - may indicate a dominant technical effect")
    next_steps <- c(next_steps, "Check if PC1 correlates with batch, library size, or other technical factors")
  }

  if (var_exp[1] + var_exp[2] < 30) {
    observations <- c(observations, "Relatively low variance in first two PCs suggests complex data structure")
  }

  if (length(issues) == 0) {
    issues <- list("No major concerns from PCA visualization")
  }

  list(
    what_is_this = "This scatter plot shows samples projected onto the first two principal components, which capture the main axes of variation in gene expression. Sample proximity indicates expression similarity.",
    observations = observations,
    issues_checks = issues,
    next_steps = next_steps
  )
}


#' Fallback commentary for PCA scree plot
#' @param pca_result PCA results list
#' @return List with commentary fields
fallback_pca_scree <- function(pca_result) {
  var_exp <- pca_result$var_explained
  n_pcs <- min(10, length(var_exp))

  cumvar <- cumsum(var_exp[1:n_pcs])
  pcs_for_50 <- which(cumvar >= 50)[1]
  pcs_for_80 <- which(cumvar >= 80)[1]

  observations <- list(
    sprintf("First PC explains %.1f%% of variance", var_exp[1]),
    sprintf("First %d PCs explain %.1f%% of variance", n_pcs, cumvar[n_pcs])
  )

  if (!is.na(pcs_for_50)) {
    observations <- c(observations, sprintf("50%% variance captured in first %d PCs", pcs_for_50))
  }

  issues <- list()
  next_steps <- list()

  if (var_exp[1] > 50) {
    issues <- c(issues, "First PC dominates variance - investigate whether this reflects biology or a technical artifact")
  }

  if (!is.na(pcs_for_80) && pcs_for_80 <= 3) {
    observations <- c(observations, "Data structure is well-captured by few components")
  } else if (is.na(pcs_for_80) || pcs_for_80 > 6) {
    observations <- c(observations, "Variance is distributed across many components")
  }

  if (length(issues) == 0) {
    issues <- list("Variance distribution appears reasonable")
  }

  list(
    what_is_this = "This scree plot shows the percentage of variance explained by each principal component. It helps determine how many PCs are needed to capture the main structure in the data.",
    observations = observations,
    issues_checks = issues,
    next_steps = next_steps
  )
}


#' Fallback commentary for volcano plot
#' @param de_summary DE summary data frame
#' @param contrast Contrast name
#' @return List with commentary fields
fallback_volcano <- function(de_summary, contrast) {
  # Find matching contrast
  contrast_clean <- gsub("_", " ", contrast)
  row <- de_summary[grepl(contrast, de_summary$contrast, fixed = TRUE) |
                    grepl(contrast_clean, de_summary$contrast), ]

  if (nrow(row) == 0) {
    return(list(
      what_is_this = "This volcano plot shows statistical significance (-log10 adjusted p-value) versus effect size (log2 fold change) for all genes.",
      observations = list("Summary statistics not available for this contrast"),
      issues_checks = list(),
      next_steps = list()
    ))
  }

  row <- row[1, ]

  observations <- list(
    sprintf("Total significant genes: %d", row$significant),
    sprintf("Upregulated: %d genes", row$up),
    sprintf("Downregulated: %d genes", row$down),
    sprintf("Percent significant: %.1f%%", row$pct_significant)
  )

  issues <- list()
  next_steps <- list()

  if (row$significant == 0) {
    issues <- c(issues, "No significantly differentially expressed genes detected")
    next_steps <- c(next_steps, "Consider relaxing significance threshold or investigating sample variability")
  } else if (row$significant < 10) {
    issues <- c(issues, "Very few significant genes detected")
    next_steps <- c(next_steps, "Low statistical power may limit biological interpretation")
  }

  if (row$significant > 0) {
    ratio <- row$up / max(row$down, 1)
    if (ratio > 3) {
      observations <- c(observations, "Strong bias toward upregulation")
    } else if (ratio < 0.33) {
      observations <- c(observations, "Strong bias toward downregulation")
    } else {
      observations <- c(observations, "Balanced distribution of up/down regulation")
    }
  }

  if (row$pct_significant > 30) {
    issues <- c(issues, "High proportion of significant genes - verify experimental conditions")
  }

  if (length(issues) == 0) {
    issues <- list("Results appear reasonable for differential expression analysis")
  }

  list(
    what_is_this = "This volcano plot displays the relationship between fold change magnitude (x-axis) and statistical significance (y-axis) for all tested genes. Points colored red (up) or blue (down) are statistically significant.",
    observations = observations,
    issues_checks = issues,
    next_steps = next_steps
  )
}


#' Fallback commentary for MA plot
#' @param de_summary DE summary data frame
#' @param contrast Contrast name
#' @return List with commentary fields
fallback_ma_plot <- function(de_summary, contrast) {
  row <- de_summary[grepl(contrast, de_summary$contrast, fixed = TRUE), ]

  observations <- list(
    "MA plot shows fold change as a function of mean expression",
    "Lowly expressed genes typically show higher variance in fold change"
  )

  if (nrow(row) > 0) {
    observations <- c(observations, sprintf("This contrast has %d significant genes", row$significant[1]))
  }

  issues <- list()
  next_steps <- list()

  issues <- c(issues, "Check for asymmetry or trends that might indicate normalization issues")

  list(
    what_is_this = "The MA plot shows log2 fold change versus mean expression. It helps identify expression-dependent biases and verify that fold changes are symmetrically distributed around zero.",
    observations = observations,
    issues_checks = issues,
    next_steps = next_steps
  )
}


#' Fallback commentary for heatmap
#' @param de_summary DE summary (unused but for consistency)
#' @param contrast Contrast name
#' @return List with commentary fields
fallback_heatmap <- function(de_summary = NULL, contrast = NULL) {
  list(
    what_is_this = "This heatmap shows expression patterns of top differentially expressed genes across samples. Colors indicate row-scaled expression values (z-scores).",
    observations = list(
      "Samples are clustered by expression similarity",
      "Genes are clustered by co-expression patterns",
      "Look for clear separation between experimental groups"
    ),
    issues_checks = list(
      "Check if samples cluster by condition (expected) or by batch (potential issue)",
      "Outlier samples may cluster separately from their group"
    ),
    next_steps = list(
      "Examine specific gene clusters for biological themes",
      "Consider pathway analysis of co-expressed gene modules"
    )
  )
}


#' Fallback commentary for pathway plot
#' @return List with commentary fields
fallback_pathway <- function() {
  list(
    what_is_this = "This dot plot shows enrichment scores and significance for pathway/gene set analysis. Dot size indicates gene set size; color indicates significance.",
    observations = list(
      "Pathways are ranked by enrichment score or significance",
      "Positive/negative enrichment indicates direction of regulation"
    ),
    issues_checks = list(
      "Very broad pathways may be less biologically informative",
      "Check for redundant pathways with overlapping gene content"
    ),
    next_steps = list(
      "Focus on pathways most relevant to your biological question",
      "Consider leading-edge analysis to identify key driver genes"
    )
  )
}


#' Get commentary for a specific figure
#' @param commentary_tbl Commentary table
#' @param figure_id Figure identifier
#' @return Commentary list or NULL
get_figure_commentary <- function(commentary_tbl, figure_id) {
  commentary_list <- attr(commentary_tbl, "commentary_list")

  if (!is.null(commentary_list) && figure_id %in% names(commentary_list)) {
    return(commentary_list[[figure_id]])
  }

  # Fall back to table lookup
  row <- commentary_tbl[commentary_tbl$figure_id == figure_id, ]
  if (nrow(row) == 0) {
    return(NULL)
  }

  # Reconstruct from table
  list(
    figure_id = row$figure_id,
    title = row$title,
    what_is_this = row$what_is_this,
    observations = strsplit(row$observations, " \\| ")[[1]],
    issues_checks = strsplit(row$issues_checks, " \\| ")[[1]],
    next_steps = strsplit(row$next_steps, " \\| ")[[1]],
    confidence = row$confidence,
    backend = row$backend
  )
}


#' Format commentary as HTML for Rmd
#' @param commentary Commentary list
#' @param include_static Include static educational text
#' @return HTML string
format_commentary_html <- function(commentary, include_static = TRUE) {
  if (is.null(commentary)) {
    return("<div class='commentary-block commentary-unavailable'><p><em>Commentary unavailable</em></p></div>")
  }

  html <- "<div class='commentary-block'>\n"

  # What is this figure
  html <- paste0(html, "<div class='commentary-section'>\n")
  html <- paste0(html, "<h5>What This Figure Shows</h5>\n")
  html <- paste0(html, "<p>", htmltools::htmlEscape(commentary$what_is_this), "</p>\n")
  html <- paste0(html, "</div>\n")

  # Observations
  if (length(commentary$observations) > 0 && commentary$observations[1] != "") {
    html <- paste0(html, "<div class='commentary-section'>\n")
    html <- paste0(html, "<h5>Key Observations</h5>\n<ul>\n")
    for (obs in commentary$observations) {
      html <- paste0(html, "<li>", htmltools::htmlEscape(obs), "</li>\n")
    }
    html <- paste0(html, "</ul></div>\n")
  }

  # Issues and checks
  if (length(commentary$issues_checks) > 0 && commentary$issues_checks[1] != "") {
    html <- paste0(html, "<div class='commentary-section'>\n")
    html <- paste0(html, "<h5>Potential Issues / Checks</h5>\n<ul>\n")
    for (issue in commentary$issues_checks) {
      html <- paste0(html, "<li>", htmltools::htmlEscape(issue), "</li>\n")
    }
    html <- paste0(html, "</ul></div>\n")
  }

  # Next steps
  if (length(commentary$next_steps) > 0 && commentary$next_steps[1] != "") {
    html <- paste0(html, "<div class='commentary-section'>\n")
    html <- paste0(html, "<h5>Recommended Next Steps</h5>\n<ul>\n")
    for (step in commentary$next_steps) {
      html <- paste0(html, "<li>", htmltools::htmlEscape(step), "</li>\n")
    }
    html <- paste0(html, "</ul></div>\n")
  }

  # Footer with metadata
  html <- paste0(html, "<div class='commentary-footer'>\n")
  html <- paste0(html, "<small>Generated by: ", commentary$backend %||% "unknown")
  if (!is.null(commentary$confidence) && commentary$confidence != "none") {
    html <- paste0(html, " | Confidence: ", commentary$confidence)
  }
  html <- paste0(html, "</small>\n</div>\n")

  html <- paste0(html, "</div>\n")

  html
}


#' Get static educational text for a plot type
#' @param plot_type Type of plot
#' @return HTML string with educational context
get_static_explanation <- function(plot_type) {
  explanations <- list(
    library_sizes = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>Library size (sequencing depth) directly affects the statistical power to detect differential expression.
      Large differences in library sizes between samples require appropriate normalization. Samples with very low
      library sizes may need to be excluded as they can have unreliable expression estimates.</p>
      </details>
      </div>
    ",
    detected_genes = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>The number of detected genes reflects both sequencing depth and RNA quality. Samples with notably
      fewer detected genes may have degraded RNA or other quality issues. Consistent detection across samples
      supports reliable downstream analysis.</p>
      </details>
      </div>
    ",
    sample_correlation = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>Sample-to-sample correlation reveals the global similarity structure of your data. Replicates should
      cluster together with high correlation. Low correlation with other samples may indicate outliers,
      batch effects, or mislabeled samples.</p>
      </details>
      </div>
    ",
    pca_plot = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>PCA reduces the high-dimensional gene expression data to a 2D visualization. The first principal
      components capture the largest sources of variation. Ideally, samples should separate by biological
      condition (e.g., treatment vs control) rather than technical factors (e.g., batch).</p>
      </details>
      </div>
    ",
    pca_scree = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>The scree plot helps assess data complexity. If PC1 explains most variance, there's a dominant
      signal (biological or technical). Distributed variance across many PCs suggests complex, multifactorial
      variation in the data.</p>
      </details>
      </div>
    ",
    volcano = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>Volcano plots combine fold change (biological significance) with statistical significance. Genes
      in the upper corners are both strongly changed and highly significant, making them priority candidates.
      The symmetry around x=0 can indicate normalization quality.</p>
      </details>
      </div>
    ",
    ma_plot = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>MA plots (M = log ratio, A = mean average) show how fold changes depend on expression level. They
      help identify expression-dependent biases. A good MA plot shows fold changes centered around zero across
      all expression levels, with increased spread at low expression.</p>
      </details>
      </div>
    ",
    heatmap = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>Heatmaps of top DE genes visualize expression patterns across samples. Clustering can reveal sample
      relationships and gene modules. Clear separation between conditions and co-expression of functionally
      related genes supports biological interpretation.</p>
      </details>
      </div>
    ",
    dotplot = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>Pathway enrichment plots summarize over-representation or gene set enrichment analysis results.
      They help interpret DE results in biological context by identifying coordinated changes in known
      pathways or functional categories.</p>
      </details>
      </div>
    "
  )

  if (plot_type %in% names(explanations)) {
    return(explanations[[plot_type]])
  }

  ""
}


# Null-coalescing operator for R
`%||%` <- function(x, y) if (is.null(x)) y else x
