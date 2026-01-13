# R/06_differential_expression.R
# Functions for differential expression analysis

#' Run differential expression analysis with DESeq2
#' @param dds Normalized DESeqDataSet object
#' @param contrasts List of contrasts to test
#' @param group_col Primary grouping column (used when no contrasts specified)
#' @param alpha Significance threshold
#' @param lfc_threshold Log2 fold change threshold
#' @return List of DE results tables
run_differential_expression <- function(dds,
                                        contrasts = list(),
                                        group_col = NULL,
                                        alpha = 0.05,
                                        lfc_threshold = 0) {

  # Run full DESeq2 analysis
  message("Running DESeq2 differential expression analysis...")
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  results_list <- list()

  # If no contrasts specified, get all pairwise comparisons for the main grouping variable
  if (length(contrasts) == 0) {
    design_vars <- all.vars(DESeq2::design(dds))

    # Use group_col if specified and valid, otherwise use first design variable
    if (!is.null(group_col) && group_col %in% design_vars) {
      main_var <- group_col
    } else {
      # Fall back to first variable in design (typically the main factor)
      main_var <- design_vars[1]
    }

    levels <- levels(SummarizedExperiment::colData(dds)[[main_var]])

    if (length(levels) == 2) {
      # Simple two-group comparison
      contrast_name <- paste0(main_var, "_", levels[2], "_vs_", levels[1])
      contrasts[[contrast_name]] <- list(
        factor = main_var,
        numerator = levels[2],
        denominator = levels[1]
      )
    } else {
      # Multiple groups - compare each to first level
      for (lvl in levels[-1]) {
        contrast_name <- paste0(main_var, "_", lvl, "_vs_", levels[1])
        contrasts[[contrast_name]] <- list(
          factor = main_var,
          numerator = lvl,
          denominator = levels[1]
        )
      }
    }

    message("Auto-generated ", length(contrasts), " contrast(s) for factor: ", main_var)
  }

  # Run each contrast
  for (contrast_name in names(contrasts)) {
    contrast <- contrasts[[contrast_name]]

    message("Processing contrast: ", contrast_name)

    tryCatch({
      # Get results
      res <- DESeq2::results(
        dds,
        contrast = c(contrast$factor, contrast$numerator, contrast$denominator),
        alpha = alpha,
        lfcThreshold = lfc_threshold
      )

      # Shrink log2 fold changes for visualization
      # Try ashr first, fall back to apeglm or normal if not available
      shrink_type <- if (requireNamespace("ashr", quietly = TRUE)) {
        "ashr"
      } else {
        "normal"
      }

      res_shrunk <- tryCatch({
        DESeq2::lfcShrink(
          dds,
          contrast = c(contrast$factor, contrast$numerator, contrast$denominator),
          type = shrink_type,
          quiet = TRUE
        )
      }, error = function(e) {
        # If shrinkage fails, use original LFC values
        warning("LFC shrinkage failed: ", e$message, ". Using unshrunk values.")
        res
      })

      # Combine results
      res_df <- as.data.frame(res)
      res_df$gene_id <- rownames(res_df)
      res_df$log2FoldChange_shrunk <- res_shrunk$log2FoldChange
      res_df$contrast <- contrast_name

      # Add significance flags
      res_df$significant <- !is.na(res_df$padj) & res_df$padj < alpha
      res_df$direction <- ifelse(res_df$log2FoldChange > 0, "up", "down")
      res_df$direction[!res_df$significant] <- "ns"

      # Reorder columns
      res_df <- res_df[, c(
        "gene_id", "contrast", "baseMean", "log2FoldChange",
        "log2FoldChange_shrunk", "lfcSE", "stat", "pvalue", "padj",
        "significant", "direction"
      )]

      # Sort by adjusted p-value
      res_df <- res_df[order(res_df$padj), ]

      results_list[[contrast_name]] <- res_df

      # Report summary
      n_sig <- sum(res_df$significant, na.rm = TRUE)
      n_up <- sum(res_df$direction == "up", na.rm = TRUE)
      n_down <- sum(res_df$direction == "down", na.rm = TRUE)
      message("  ", contrast_name, ": ", n_sig, " significant genes (",
              n_up, " up, ", n_down, " down) at padj < ", alpha)

    }, error = function(e) {
      warning("Failed to process contrast '", contrast_name, "': ", e$message)
      results_list[[contrast_name]] <- NULL
    })
  }

  results_list
}


#' Annotate DE results with gene information
#' @param de_results_list List of DE results
#' @param annotation Gene annotation data frame
#' @return List of annotated DE results
annotate_de_results <- function(de_results_list, annotation) {

  annotated_list <- lapply(names(de_results_list), function(contrast_name) {
    res <- de_results_list[[contrast_name]]

    if (is.null(res)) return(NULL)

    # Merge with annotation
    res_annotated <- merge(
      res,
      annotation,
      by = "gene_id",
      all.x = TRUE
    )

    # Reorder by adjusted p-value
    res_annotated <- res_annotated[order(res_annotated$padj), ]

    res_annotated
  })

  names(annotated_list) <- names(de_results_list)

  message("Annotated ", length(annotated_list), " DE result table(s)")

  annotated_list
}


#' Save DE results to CSV files
#' @param de_results_annotated List of annotated DE results
#' @param output_dir Output directory
#' @return Vector of output file paths
save_de_results <- function(de_results_annotated, output_dir = "outputs/de_results") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output_files <- c()

  for (contrast_name in names(de_results_annotated)) {
    res <- de_results_annotated[[contrast_name]]

    if (is.null(res)) next

    # Clean filename
    clean_name <- gsub("[^a-zA-Z0-9_-]", "_", contrast_name)
    output_file <- file.path(output_dir, paste0("de_results_", clean_name, ".csv"))

    write.csv(res, output_file, row.names = FALSE)
    output_files <- c(output_files, output_file)
  }

  message("Saved ", length(output_files), " DE result file(s) to ", output_dir)

  output_files
}


#' Summarize DE results across contrasts
#' @param de_results_annotated List of annotated DE results
#' @param alpha Significance threshold
#' @param lfc_threshold Log2 fold change threshold
#' @return Summary data frame
summarize_de_results <- function(de_results_annotated,
                                  alpha = 0.05,
                                  lfc_threshold = 0) {

  summary_list <- lapply(names(de_results_annotated), function(contrast_name) {
    res <- de_results_annotated[[contrast_name]]

    if (is.null(res)) {
      return(data.frame(
        contrast = contrast_name,
        total_genes = 0,
        significant = 0,
        up = 0,
        down = 0,
        pct_significant = 0,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      contrast = contrast_name,
      total_genes = nrow(res),
      significant = sum(res$significant, na.rm = TRUE),
      up = sum(res$direction == "up", na.rm = TRUE),
      down = sum(res$direction == "down", na.rm = TRUE),
      pct_significant = round(sum(res$significant, na.rm = TRUE) / nrow(res) * 100, 2),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_list)

  summary_df
}


#' Generate DE visualization plots
#' @param de_results_annotated List of annotated DE results
#' @param vst_counts VST-transformed counts for heatmaps
#' @param metadata Sample metadata
#' @param group_col Column for grouping
#' @param output_dir Output directory
#' @return List of plot file paths
generate_de_plots <- function(de_results_annotated,
                               vst_counts,
                               metadata,
                               group_col = "condition",
                               output_dir = "outputs/plots") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  plot_files <- list()

  for (contrast_name in names(de_results_annotated)) {
    res <- de_results_annotated[[contrast_name]]

    if (is.null(res) || nrow(res) == 0) next

    clean_name <- gsub("[^a-zA-Z0-9_-]", "_", contrast_name)

    # 1. Volcano plot
    res_plot <- res[!is.na(res$padj), ]

    p_volcano <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = direction), alpha = 0.5, size = 1) +
      scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey60")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
      labs(
        title = paste("Volcano Plot:", contrast_name),
        x = "log2 Fold Change",
        y = "-log10(adjusted p-value)"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")

    volcano_file <- file.path(output_dir, paste0("volcano_", clean_name, ".png"))
    ggsave(volcano_file, p_volcano, width = 8, height = 6)
    plot_files[[paste0("volcano_", clean_name)]] <- volcano_file

    # 2. MA plot
    p_ma <- ggplot(res_plot, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
      geom_point(aes(color = direction), alpha = 0.5, size = 1) +
      scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey60")) +
      geom_hline(yintercept = 0, color = "grey40") +
      labs(
        title = paste("MA Plot:", contrast_name),
        x = "log10(Mean Expression)",
        y = "log2 Fold Change"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")

    ma_file <- file.path(output_dir, paste0("ma_plot_", clean_name, ".png"))
    ggsave(ma_file, p_ma, width = 8, height = 6)
    plot_files[[paste0("ma_", clean_name)]] <- ma_file

    # 3. Top genes heatmap
    top_genes <- head(res$gene_id[res$significant], 50)

    if (length(top_genes) >= 5) {
      top_genes_present <- top_genes[top_genes %in% rownames(vst_counts)]

      if (length(top_genes_present) >= 5) {
        heatmap_mat <- vst_counts[top_genes_present, , drop = FALSE]

        # Scale rows
        heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

        heatmap_file <- file.path(output_dir, paste0("heatmap_top_", clean_name, ".png"))
        png(heatmap_file, width = 800, height = 1000)

        # Annotation for columns
        if (group_col %in% colnames(metadata)) {
          annotation_col <- data.frame(
            Group = metadata[[group_col]],
            row.names = rownames(metadata)
          )
        } else {
          annotation_col <- NA
        }

        pheatmap::pheatmap(
          heatmap_mat_scaled,
          main = paste("Top DE Genes:", contrast_name),
          annotation_col = annotation_col,
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          show_rownames = length(top_genes_present) <= 30,
          fontsize_row = 8,
          border_color = NA
        )
        dev.off()
        plot_files[[paste0("heatmap_", clean_name)]] <- heatmap_file
      }
    }
  }

  message("Generated DE plots in ", output_dir)

  plot_files
}
