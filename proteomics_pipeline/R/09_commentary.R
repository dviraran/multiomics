# =============================================================================
# Figure Commentary Generation
# =============================================================================
# Generates AI-powered or data-driven commentary for pipeline figures
# Supports Claude Vision API, OpenAI GPT-4 Vision, and deterministic fallback

#' Build table of all figures with metadata
#'
#' @param qc_results QC results from run_qc_analysis()
#' @param de_results Differential analysis results
#' @param pathway_results Pathway analysis results
#' @param ppi_results PPI network results
#' @param advanced_stats Advanced statistics results
#' @param config Configuration list
#' @return Data frame with figure metadata
build_figures_table <- function(qc_results = NULL,
                                 de_results = NULL,
                                 pathway_results = NULL,
                                 ppi_results = NULL,
                                 advanced_stats = NULL,
                                 config = NULL) {
  log_message("Building figures metadata table...")

  figures <- list()
  output_dir <- config$output$output_dir %||% "outputs"

  # QC plots
  if (!is.null(qc_results)) {
    # PCA PC1 vs PC2
    pca_file <- file.path(output_dir, "qc", "pca_pc1_pc2.png")
    if (file.exists(pca_file)) {
      figures$pca_pc1_pc2 <- data.frame(
        figure_id = "pca_pc1_pc2",
        filepath = pca_file,
        plot_type = "scatter",
        section = "QC",
        title = "PCA: PC1 vs PC2",
        description = "Principal component analysis showing sample clustering",
        x_axis = "PC1",
        y_axis = "PC2",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    # PCA Scree plot
    scree_file <- file.path(output_dir, "qc", "pca_scree.png")
    if (file.exists(scree_file)) {
      figures$pca_scree <- data.frame(
        figure_id = "pca_scree",
        filepath = scree_file,
        plot_type = "barplot",
        section = "QC",
        title = "PCA Scree Plot",
        description = "Variance explained by each principal component",
        x_axis = "Principal Component",
        y_axis = "% Variance Explained",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    # Correlation heatmap
    cor_file <- file.path(output_dir, "qc", "correlation_heatmap.png")
    if (file.exists(cor_file)) {
      figures$correlation_heatmap <- data.frame(
        figure_id = "correlation_heatmap",
        filepath = cor_file,
        plot_type = "heatmap",
        section = "QC",
        title = "Sample Correlation Heatmap",
        description = "Pearson correlation matrix between samples",
        x_axis = "Samples",
        y_axis = "Samples",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    # Clustering dendrogram
    dend_file <- file.path(output_dir, "qc", "clustering_dendrogram.png")
    if (file.exists(dend_file)) {
      figures$clustering_dendrogram <- data.frame(
        figure_id = "clustering_dendrogram",
        filepath = dend_file,
        plot_type = "dendrogram",
        section = "QC",
        title = "Sample Clustering Dendrogram",
        description = "Hierarchical clustering of samples based on expression",
        x_axis = "Samples",
        y_axis = "Distance",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    # UMAP
    umap_file <- file.path(output_dir, "qc", "umap.png")
    if (file.exists(umap_file)) {
      figures$umap <- data.frame(
        figure_id = "umap",
        filepath = umap_file,
        plot_type = "scatter",
        section = "QC",
        title = "UMAP",
        description = "UMAP dimensionality reduction showing sample relationships",
        x_axis = "UMAP1",
        y_axis = "UMAP2",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    # Batch effect PCA (if exists)
    batch_file <- file.path(output_dir, "qc", "batch_effect_pca.png")
    if (file.exists(batch_file)) {
      figures$batch_effect_pca <- data.frame(
        figure_id = "batch_effect_pca",
        filepath = batch_file,
        plot_type = "scatter",
        section = "QC",
        title = "PCA Colored by Batch",
        description = "PCA visualization to assess batch effects",
        x_axis = "PC1",
        y_axis = "PC2",
        contrast = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  # DE plots (volcano, MA, p-value distribution per contrast)
  if (!is.null(de_results) && !is.null(de_results$results)) {
    for (contrast_name in names(de_results$results)) {
      clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))

      # Volcano plot
      volcano_file <- file.path(output_dir, "plots", paste0("volcano_", clean_name, ".png"))
      if (file.exists(volcano_file)) {
        figures[[paste0("volcano_", clean_name)]] <- data.frame(
          figure_id = paste0("volcano_", clean_name),
          filepath = volcano_file,
          plot_type = "volcano",
          section = "Differential",
          title = paste("Volcano Plot:", contrast_name),
          description = "Statistical significance vs fold change",
          x_axis = "Log2 Fold Change",
          y_axis = "-Log10(P-value)",
          contrast = contrast_name,
          stringsAsFactors = FALSE
        )
      }

      # MA plot
      ma_file <- file.path(output_dir, "plots", paste0("ma_plot_", clean_name, ".png"))
      if (file.exists(ma_file)) {
        figures[[paste0("ma_", clean_name)]] <- data.frame(
          figure_id = paste0("ma_", clean_name),
          filepath = ma_file,
          plot_type = "ma_plot",
          section = "Differential",
          title = paste("MA Plot:", contrast_name),
          description = "Fold change vs average intensity",
          x_axis = "Average Log2 Intensity",
          y_axis = "Log2 Fold Change",
          contrast = contrast_name,
          stringsAsFactors = FALSE
        )
      }

      # P-value distribution
      pval_file <- file.path(output_dir, "plots", paste0("pvalue_dist_", clean_name, ".png"))
      if (file.exists(pval_file)) {
        figures[[paste0("pvalue_", clean_name)]] <- data.frame(
          figure_id = paste0("pvalue_", clean_name),
          filepath = pval_file,
          plot_type = "histogram",
          section = "Differential",
          title = paste("P-value Distribution:", contrast_name),
          description = "Distribution of raw p-values",
          x_axis = "P-value",
          y_axis = "Count",
          contrast = contrast_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Pathway plots
  if (!is.null(pathway_results) && !is.null(pathway_results$results)) {
    for (contrast_name in names(pathway_results$results)) {
      clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))

      # ORA plot
      ora_file <- file.path(output_dir, "plots", paste0("pathway_ora_", clean_name, ".png"))
      if (file.exists(ora_file)) {
        figures[[paste0("pathway_ora_", clean_name)]] <- data.frame(
          figure_id = paste0("pathway_ora_", clean_name),
          filepath = ora_file,
          plot_type = "dotplot",
          section = "Pathway",
          title = paste("ORA Enrichment:", contrast_name),
          description = "Over-representation analysis results",
          x_axis = "-Log10(Adjusted P-value)",
          y_axis = "Pathway",
          contrast = contrast_name,
          stringsAsFactors = FALSE
        )
      }

      # GSEA plot
      gsea_file <- file.path(output_dir, "plots", paste0("pathway_gsea_", clean_name, ".png"))
      if (file.exists(gsea_file)) {
        figures[[paste0("pathway_gsea_", clean_name)]] <- data.frame(
          figure_id = paste0("pathway_gsea_", clean_name),
          filepath = gsea_file,
          plot_type = "barplot",
          section = "Pathway",
          title = paste("GSEA Enrichment:", contrast_name),
          description = "Gene set enrichment analysis results",
          x_axis = "Normalized Enrichment Score",
          y_axis = "Pathway",
          contrast = contrast_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Combine all figures
  if (length(figures) == 0) {
    log_message("No figures found to document")
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
    figures_tbl$organism <- config$input$organism
    figures_tbl$condition_col <- config$design$condition_column
  }

  log_message("Built figures table with ", nrow(figures_tbl), " figures")

  figures_tbl
}


#' Generate commentary for all figures
#'
#' @param figures_tbl Figures metadata table
#' @param config Pipeline configuration
#' @param qc_results QC results for data-driven fallback
#' @param de_results DE results for data-driven fallback
#' @param imputed_data Imputed data for context
#' @param ppi_results PPI network results
#' @param advanced_stats Advanced statistics results
#' @param output_dir Output directory for commentary files
#' @return Data frame with commentary for each figure
generate_all_commentary <- function(figures_tbl,
                                     config,
                                     qc_results = NULL,
                                     de_results = NULL,
                                     imputed_data = NULL,
                                     ppi_results = NULL,
                                     advanced_stats = NULL,
                                     output_dir = NULL) {

  if (is.null(output_dir)) {
    output_dir <- file.path(config$output$output_dir %||% "outputs", "commentary")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Check if commentary is enabled
  commentary_enabled <- isTRUE(config$commentary$enabled)
  backend <- config$commentary$backend %||% "none"

  if (!commentary_enabled) {
    log_message("Commentary disabled in config. Using placeholder text.")
    backend <- "none"
  }

  commentary_list <- list()

  for (i in seq_len(nrow(figures_tbl))) {
    fig <- figures_tbl[i, ]

    log_message("Generating commentary for: ", fig$figure_id)

    # Check if file exists
    if (!file.exists(fig$filepath)) {
      log_message("  Warning: Figure file not found: ", fig$filepath)
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
      qc_results = qc_results,
      de_results = de_results,
      imputed_data = imputed_data
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
        qc_results = qc_results,
        de_results = de_results
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

  log_message("Saved commentary to ", output_dir)
  log_message("  - Individual JSONs: ", length(commentary_list), " files")
  log_message("  - Combined CSV: ", csv_path)

  # Return as list for easier access in Rmd
  attr(commentary_tbl, "commentary_list") <- commentary_list

  commentary_tbl
}


#' Build context information for a figure
#'
#' @param fig Figure row from figures_tbl
#' @param config Pipeline config
#' @param ... Additional data objects
#' @return List with context information
build_figure_context <- function(fig,
                                  config,
                                  qc_results = NULL,
                                  de_results = NULL,
                                  imputed_data = NULL) {

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
    context$organism <- config$input$organism
    context$condition_col <- config$design$condition_column
    context$design_formula <- config$design$design_formula
    context$normalization <- config$processing$normalization_method
    context$imputation <- config$imputation$method
  }

  # Add sample count
  if (!is.null(imputed_data)) {
    context$n_samples <- ncol(imputed_data$matrix)
    context$n_features <- nrow(imputed_data$matrix)
  }

  # Add PCA variance info
  if (!is.null(qc_results) && !is.null(qc_results$pca)) {
    var_exp <- qc_results$pca$variance_explained
    context$pca_var_pc1 <- round(var_exp[1], 1)
    context$pca_var_pc2 <- round(var_exp[2], 1)
  }

  # Add outlier info
  if (!is.null(qc_results) && !is.null(qc_results$outliers)) {
    context$n_outliers <- qc_results$outliers$n_outliers
    if (context$n_outliers > 0) {
      flagged <- qc_results$outliers$flagged_samples
      context$outlier_samples <- flagged$sample[flagged$is_outlier]
    }
  }

  # Add DE summary for relevant plots
  if (!is.null(de_results) && !is.na(fig$contrast)) {
    contrast_name <- fig$contrast
    if (contrast_name %in% names(de_results$results)) {
      res <- de_results$results[[contrast_name]]
      context$de_significant <- res$n_sig
      context$de_up <- res$n_up
      context$de_down <- res$n_down
      context$de_total <- nrow(res$table)
    }
  }

  context
}


#' Run Claude Vision API for commentary
#'
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
  script_path <- file.path("scripts", "figure_commentary_claude.py")

  if (!file.exists(script_path)) {
    log_message("  Claude commentary script not found: ", script_path)
    return(create_placeholder_commentary(figure_id, reason = "Script not found"))
  }

  # Build command
  model <- config$commentary$claude_model %||% "claude-sonnet-4-20250514"
  max_tokens <- config$commentary$max_tokens %||% 1500

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
  max_retries <- config$commentary$max_retries %||% 2
  retry_delay <- config$commentary$retry_delay %||% 5

  for (attempt in seq_len(max_retries + 1)) {
    result <- tryCatch({
      system(cmd, intern = TRUE, ignore.stderr = FALSE, timeout = 120)

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
        log_message("  Attempt ", attempt, " failed: ", e$message, ". Retrying...")
        Sys.sleep(retry_delay)
        return(NULL)
      } else {
        log_message("  Claude API call failed: ", e$message)
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
#'
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
  script_path <- file.path("scripts", "figure_commentary_openai.py")

  if (!file.exists(script_path)) {
    log_message("  OpenAI commentary script not found: ", script_path)
    return(create_placeholder_commentary(figure_id, reason = "Script not found"))
  }

  # Build command
  model <- config$commentary$openai_model %||% "gpt-4o"
  max_tokens <- config$commentary$max_tokens %||% 1500

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
  max_retries <- config$commentary$max_retries %||% 2
  retry_delay <- config$commentary$retry_delay %||% 5

  for (attempt in seq_len(max_retries + 1)) {
    result <- tryCatch({
      system(cmd, intern = TRUE, ignore.stderr = FALSE, timeout = 120)

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
        log_message("  Attempt ", attempt, " failed: ", e$message, ". Retrying...")
        Sys.sleep(retry_delay)
        return(NULL)
      } else {
        log_message("  OpenAI API call failed: ", e$message)
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
#'
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
#'
#' @param fig Figure metadata row
#' @param context Figure context
#' @param qc_results QC results
#' @param de_results DE results
#' @return Commentary list
generate_fallback_commentary <- function(fig,
                                          context,
                                          qc_results = NULL,
                                          de_results = NULL) {

  commentary <- list(
    figure_id = fig$figure_id,
    title = fig$title,
    backend = "data_driven",
    generated_at = as.character(Sys.time()),
    confidence = "medium"
  )

  # Route to specific fallback function based on figure type
  if (fig$figure_id == "pca_pc1_pc2" && !is.null(qc_results)) {
    commentary <- c(commentary, fallback_pca(qc_results, context))
  } else if (fig$figure_id == "pca_scree" && !is.null(qc_results)) {
    commentary <- c(commentary, fallback_pca_scree(qc_results))
  } else if (fig$figure_id == "correlation_heatmap" && !is.null(qc_results)) {
    commentary <- c(commentary, fallback_correlation(qc_results))
  } else if (fig$figure_id == "clustering_dendrogram") {
    commentary <- c(commentary, fallback_dendrogram())
  } else if (fig$figure_id == "umap" && !is.null(qc_results)) {
    commentary <- c(commentary, fallback_umap(qc_results, context))
  } else if (fig$figure_id == "batch_effect_pca" && !is.null(qc_results)) {
    commentary <- c(commentary, fallback_batch_effect(qc_results))
  } else if (grepl("^volcano_", fig$figure_id) && !is.null(de_results)) {
    commentary <- c(commentary, fallback_volcano(de_results, fig$contrast))
  } else if (grepl("^ma_", fig$figure_id) && !is.null(de_results)) {
    commentary <- c(commentary, fallback_ma_plot(de_results, fig$contrast))
  } else if (grepl("^pvalue_", fig$figure_id)) {
    commentary <- c(commentary, fallback_pvalue_dist())
  } else if (grepl("^pathway_ora_", fig$figure_id)) {
    commentary <- c(commentary, fallback_pathway_ora())
  } else if (grepl("^pathway_gsea_", fig$figure_id)) {
    commentary <- c(commentary, fallback_pathway_gsea())
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


#' Fallback commentary for PCA plot
fallback_pca <- function(qc_results, context = NULL) {
  pca <- qc_results$pca
  var_exp <- pca$variance_explained
  outliers <- qc_results$outliers

  observations <- list(
    sprintf("PC1 explains %.1f%% of variance", var_exp[1]),
    sprintf("PC2 explains %.1f%% of variance", var_exp[2]),
    sprintf("Total variance in PC1-PC2: %.1f%%", var_exp[1] + var_exp[2])
  )

  # Check for group separation if available
  pca_coords <- pca$coordinates
  if ("condition" %in% colnames(pca_coords)) {
    conditions <- unique(pca_coords$condition)
    if (length(conditions) >= 2) {
      # Calculate group centroids
      centroids <- aggregate(cbind(PC1, PC2) ~ condition, data = pca_coords, mean)
      if (nrow(centroids) >= 2) {
        mean_dist <- mean(dist(centroids[, c("PC1", "PC2")]))
        within_spread <- mean(sapply(conditions, function(g) {
          group_data <- pca_coords[pca_coords$condition == g, ]
          if (nrow(group_data) > 1) {
            mean(dist(group_data[, c("PC1", "PC2")]))
          } else 0
        }))

        if (mean_dist > within_spread * 1.5) {
          observations <- c(observations, "Conditions show clear separation in PCA space")
        } else if (mean_dist > within_spread) {
          observations <- c(observations, "Conditions show moderate separation in PCA space")
        } else {
          observations <- c(observations, "Conditions show overlap in PCA space")
        }
      }
    }
  }

  issues <- list()
  next_steps <- list()

  # Check for outliers
  if (!is.null(outliers) && outliers$n_outliers > 0) {
    flagged <- outliers$flagged_samples
    pca_outliers <- flagged$sample[flagged$pca_outlier]
    if (length(pca_outliers) > 0) {
      issues <- c(issues, sprintf("PCA-based outliers flagged: %s", paste(pca_outliers, collapse = ", ")))
      next_steps <- c(next_steps, "Investigate outlier samples for technical or biological explanations")
    }
  }

  if (var_exp[1] > 50) {
    issues <- c(issues, "PC1 explains >50% variance - may indicate a dominant effect")
    next_steps <- c(next_steps, "Check if PC1 correlates with condition, batch, or technical factors")
  }

  if (length(issues) == 0) {
    issues <- list("No major concerns from PCA visualization")
  }

  list(
    what_is_this = "This scatter plot shows samples projected onto the first two principal components, which capture the main axes of variation in protein abundance.",
    observations = observations,
    issues_checks = issues,
    next_steps = next_steps
  )
}


#' Fallback commentary for PCA scree plot
fallback_pca_scree <- function(qc_results) {
  var_exp <- qc_results$pca$variance_explained
  n_pcs <- min(10, length(var_exp))
  cumvar <- cumsum(var_exp[1:n_pcs])

  observations <- list(
    sprintf("First PC explains %.1f%% of variance", var_exp[1]),
    sprintf("First %d PCs explain %.1f%% of variance", n_pcs, cumvar[n_pcs])
  )

  pcs_for_50 <- which(cumvar >= 50)[1]
  if (!is.na(pcs_for_50)) {
    observations <- c(observations, sprintf("50%% variance captured in first %d PCs", pcs_for_50))
  }

  issues <- list()
  if (var_exp[1] > 50) {
    issues <- c(issues, "First PC dominates variance - investigate for technical artifacts")
  } else {
    issues <- list("Variance distribution appears reasonable")
  }

  list(
    what_is_this = "This scree plot shows the percentage of variance explained by each principal component, helping determine data complexity.",
    observations = observations,
    issues_checks = issues,
    next_steps = list()
  )
}


#' Fallback commentary for correlation heatmap
fallback_correlation <- function(qc_results) {
  cor_results <- qc_results$correlation

  summary <- cor_results$summary
  min_cor <- as.numeric(summary$value[summary$metric == "Min correlation"])
  mean_cor <- as.numeric(summary$value[summary$metric == "Mean correlation"])

  observations <- list(
    sprintf("Mean sample correlation: %.3f", mean_cor),
    sprintf("Minimum sample correlation: %.3f", min_cor)
  )

  if (mean_cor > 0.9) {
    observations <- c(observations, "High overall correlation suggests consistent sample profiles")
  } else if (mean_cor < 0.7) {
    observations <- c(observations, "Lower correlation may indicate batch effects or biological heterogeneity")
  }

  issues <- list()
  if (min_cor < 0.7) {
    # Find low correlation samples
    per_sample <- cor_results$per_sample
    low_samples <- per_sample$sample[per_sample$mean_correlation < 0.7]
    if (length(low_samples) > 0) {
      issues <- c(issues, sprintf("Samples with low correlation: %s", paste(low_samples, collapse = ", ")))
    }
  }

  if (length(issues) == 0) {
    issues <- list("Sample correlations appear consistent")
  }

  list(
    what_is_this = "This heatmap shows Pearson correlation between all sample pairs. High correlation (warm colors) indicates similar abundance profiles.",
    observations = observations,
    issues_checks = issues,
    next_steps = list("Investigate any samples with notably low correlation")
  )
}


#' Fallback commentary for clustering dendrogram
fallback_dendrogram <- function() {
  list(
    what_is_this = "This dendrogram shows hierarchical clustering of samples based on protein abundance profiles using Ward's method.",
    observations = list(
      "Samples are clustered by expression similarity",
      "Branch height indicates relative distance between clusters",
      "Color coding shows condition assignments"
    ),
    issues_checks = list(
      "Check if samples cluster primarily by condition (expected) or by batch (potential issue)",
      "Samples that cluster away from their condition group may need investigation"
    ),
    next_steps = list()
  )
}


#' Fallback commentary for UMAP
fallback_umap <- function(qc_results, context = NULL) {
  observations <- list(
    "UMAP preserves local structure better than PCA",
    "Proximity indicates similar abundance profiles"
  )

  if (!is.null(qc_results$umap)) {
    umap_coords <- qc_results$umap$coordinates
    if ("condition" %in% colnames(umap_coords)) {
      observations <- c(observations, "Samples are colored by experimental condition")
    }
  }

  list(
    what_is_this = "This UMAP plot provides a non-linear dimensionality reduction view of sample relationships, emphasizing local structure.",
    observations = observations,
    issues_checks = list("Compare with PCA to ensure consistent patterns"),
    next_steps = list()
  )
}


#' Fallback commentary for batch effect PCA
fallback_batch_effect <- function(qc_results) {
  batch_results <- qc_results$batch_effects

  observations <- list("PCA colored by batch to assess technical variation")

  issues <- list()

  if (!is.null(batch_results)) {
    summary <- batch_results$summary
    pc1_var <- as.numeric(summary$value[summary$metric == "Batch variance in PC1"])
    pc2_var <- as.numeric(summary$value[summary$metric == "Batch variance in PC2"])

    observations <- c(observations,
      sprintf("Batch explains %.1f%% of PC1 variance", pc1_var),
      sprintf("Batch explains %.1f%% of PC2 variance", pc2_var)
    )

    if (pc1_var > 30 || pc2_var > 30) {
      issues <- c(issues, "Substantial batch effect detected - consider batch correction")
    }
  }

  if (length(issues) == 0) {
    issues <- list("Batch effects appear manageable")
  }

  list(
    what_is_this = "This PCA plot is colored by batch to assess whether technical batch effects drive sample clustering.",
    observations = observations,
    issues_checks = issues,
    next_steps = list("Consider including batch in the statistical model if effects are strong")
  )
}


#' Fallback commentary for volcano plot
fallback_volcano <- function(de_results, contrast_name) {
  if (is.null(contrast_name) || !contrast_name %in% names(de_results$results)) {
    return(list(
      what_is_this = "This volcano plot shows statistical significance vs fold change for all proteins.",
      observations = list("Summary not available"),
      issues_checks = list(),
      next_steps = list()
    ))
  }

  res <- de_results$results[[contrast_name]]

  observations <- list(
    sprintf("Total significant proteins: %d", res$n_sig),
    sprintf("Up-regulated: %d proteins", res$n_up),
    sprintf("Down-regulated: %d proteins", res$n_down),
    sprintf("Percent significant: %.1f%%", res$n_sig / nrow(res$table) * 100)
  )

  issues <- list()
  next_steps <- list()

  if (res$n_sig == 0) {
    issues <- c(issues, "No significantly differential proteins detected")
    next_steps <- c(next_steps, "Consider relaxing thresholds or checking sample variability")
  } else if (res$n_sig < 10) {
    issues <- c(issues, "Few significant proteins detected")
  }

  if (res$n_sig > 0) {
    ratio <- res$n_up / max(res$n_down, 1)
    if (ratio > 3) {
      observations <- c(observations, "Strong bias toward up-regulation")
    } else if (ratio < 0.33) {
      observations <- c(observations, "Strong bias toward down-regulation")
    }
  }

  if (length(issues) == 0) {
    issues <- list("Results appear reasonable for differential analysis")
  }

  list(
    what_is_this = "This volcano plot shows the relationship between fold change magnitude (x-axis) and statistical significance (y-axis) for all proteins.",
    observations = observations,
    issues_checks = issues,
    next_steps = next_steps
  )
}


#' Fallback commentary for MA plot
fallback_ma_plot <- function(de_results, contrast_name) {
  list(
    what_is_this = "This MA plot shows log2 fold change vs average intensity. It helps identify intensity-dependent biases.",
    observations = list(
      "Fold changes should be symmetrically distributed around zero",
      "Higher variance at low intensities is expected",
      "Asymmetry may indicate normalization issues"
    ),
    issues_checks = list("Check for systematic trends across intensity range"),
    next_steps = list()
  )
}


#' Fallback commentary for p-value distribution
fallback_pvalue_dist <- function() {
  list(
    what_is_this = "This histogram shows the distribution of raw p-values. A flat distribution indicates no global signal; a peak near 0 indicates true positives.",
    observations = list(
      "A uniform distribution suggests no systematic differential abundance",
      "An enrichment near 0 suggests true differential proteins",
      "A peak near 1 may indicate conservative testing or batch effects"
    ),
    issues_checks = list("Anti-conservative patterns (peak at 0, dip in middle) may indicate model issues"),
    next_steps = list()
  )
}


#' Fallback commentary for ORA pathway plot
fallback_pathway_ora <- function() {
  list(
    what_is_this = "This dot plot shows over-representation analysis (ORA) results. Dot size indicates gene count; color indicates significance.",
    observations = list(
      "Pathways are ranked by significance",
      "Larger dots indicate more proteins in the pathway"
    ),
    issues_checks = list(
      "Very broad pathways may be less informative",
      "Check for redundant pathways with overlapping genes"
    ),
    next_steps = list("Focus on pathways most relevant to your biological question")
  )
}


#' Fallback commentary for GSEA pathway plot
fallback_pathway_gsea <- function() {
  list(
    what_is_this = "This bar chart shows Gene Set Enrichment Analysis (GSEA) results. NES indicates enrichment strength and direction.",
    observations = list(
      "Positive NES indicates up-regulation in the condition",
      "Negative NES indicates down-regulation",
      "Larger absolute NES indicates stronger enrichment"
    ),
    issues_checks = list("Verify that top pathways are biologically plausible"),
    next_steps = list("Examine leading-edge genes driving the enrichment")
  )
}


#' Get commentary for a specific figure
#'
#' @param commentary_tbl Commentary table with list attribute
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
#'
#' @param commentary Commentary list
#' @return HTML string
format_commentary_html <- function(commentary) {
  if (is.null(commentary)) {
    return("<div class='commentary-block commentary-unavailable'><p><em>Commentary unavailable</em></p></div>")
  }

  html <- "<div class='commentary-block'>\n"

  # What is this figure
  html <- paste0(html, "<div class='commentary-section'>\n")
  html <- paste0(html, "<h5>What This Figure Shows</h5>\n")
  html <- paste0(html, "<p>", htmlEscape(commentary$what_is_this), "</p>\n")
  html <- paste0(html, "</div>\n")

  # Observations
  if (length(commentary$observations) > 0 && commentary$observations[1] != "") {
    html <- paste0(html, "<div class='commentary-section'>\n")
    html <- paste0(html, "<h5>Key Observations</h5>\n<ul>\n")
    for (obs in commentary$observations) {
      html <- paste0(html, "<li>", htmlEscape(obs), "</li>\n")
    }
    html <- paste0(html, "</ul></div>\n")
  }

  # Issues and checks
  if (length(commentary$issues_checks) > 0 && commentary$issues_checks[1] != "") {
    html <- paste0(html, "<div class='commentary-section'>\n")
    html <- paste0(html, "<h5>Potential Issues / Checks</h5>\n<ul>\n")
    for (issue in commentary$issues_checks) {
      html <- paste0(html, "<li>", htmlEscape(issue), "</li>\n")
    }
    html <- paste0(html, "</ul></div>\n")
  }

  # Next steps
  if (length(commentary$next_steps) > 0 && commentary$next_steps[1] != "") {
    html <- paste0(html, "<div class='commentary-section'>\n")
    html <- paste0(html, "<h5>Recommended Next Steps</h5>\n<ul>\n")
    for (step in commentary$next_steps) {
      html <- paste0(html, "<li>", htmlEscape(step), "</li>\n")
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


#' HTML escape helper
htmlEscape <- function(x) {
  x <- gsub("&", "&amp;", x)
  x <- gsub("<", "&lt;", x)
  x <- gsub(">", "&gt;", x)
  x <- gsub("\"", "&quot;", x)
  x
}


#' Get static educational text for a plot type
#'
#' @param plot_type Type of plot
#' @return HTML string
get_static_explanation <- function(plot_type) {
  explanations <- list(
    pca = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>PCA reveals the major sources of variation in proteomics data. Samples should ideally separate by
      biological condition rather than technical factors like batch. Outliers and unexpected clustering
      patterns can indicate sample quality issues.</p>
      </details>
      </div>
    ",
    correlation = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>Sample correlation matrices reveal global similarity patterns. Replicates should cluster together
      with high correlation. Low correlation may indicate technical problems, mislabeling, or distinct
      biological states.</p>
      </details>
      </div>
    ",
    volcano = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>Volcano plots combine fold change with statistical significance. Proteins in the upper corners
      are both biologically and statistically significant. The symmetry around zero can indicate
      normalization quality.</p>
      </details>
      </div>
    ",
    ma_plot = "
      <div class='static-explanation'>
      <details>
      <summary><strong>Why include this plot?</strong></summary>
      <p>MA plots show how fold changes relate to abundance level. They help identify intensity-dependent
      biases. Well-normalized data should show fold changes centered around zero across all intensity levels.</p>
      </details>
      </div>
    "
  )

  if (plot_type %in% names(explanations)) {
    return(explanations[[plot_type]])
  }

  ""
}
