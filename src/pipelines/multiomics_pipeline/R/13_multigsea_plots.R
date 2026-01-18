# =============================================================================
# MultiGSEA Correlation Plots
# =============================================================================

#' Run MultiGSEA Correlation Analysis
#'
#' @param enrichment_results List containing enrichment results from `run_multiomics_enrichment`.
#' @param config Pipeline configuration list.
#'
#' @return A list of ggplot objects or file paths.
#' @export
run_multigsea_plots <- function(enrichment_results, config) {
    log_message("=== Running MultiGSEA Correlation Analysis ===")

    if (is.null(enrichment_results) || is.null(enrichment_results$per_omics)) {
        log_message("No enrichment results available for MultiGSEA.")
        return(NULL)
    }

    mg_config <- config$enrichment$multigsea %||% list()
    if (!(mg_config$run_multigsea %||% TRUE)) {
        log_message("MultiGSEA analysis disabled in config.")
        return(NULL)
    }

    p_thresh <- mg_config$pvalue_threshold %||% 0.05
    corr_method <- mg_config$correlation_method %||% "pearson"

    # Output directory
    out_dir <- file.path(config$output$output_dir, "enrichment", "multigsea")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    # Extract results per omic
    per_omics <- enrichment_results$per_omics
    omics_names <- names(per_omics)

    if (length(omics_names) < 2) {
        log_message("Need at least 2 omics with enrichment results for MultiGSEA.")
        return(NULL)
    }

    # Identify pairs
    pairs <- combn(omics_names, 2, simplify = FALSE)
    plots <- list()

    for (pair in pairs) {
        omic1 <- pair[1]
        omic2 <- pair[2]

        res1 <- per_omics[[omic1]]$results
        res2 <- per_omics[[omic2]]$results

        if (is.null(res1) || is.null(res2)) next

        # Merge by term
        # Ensure column names are consistent or use 'term'
        common_terms <- intersect(res1$term, res2$term)

        if (length(common_terms) < 3) {
            log_message("Too few common terms between ", omic1, " and ", omic2)
            next
        }

        # Align data
        df1 <- res1[match(common_terms, res1$term), ]
        df2 <- res2[match(common_terms, res2$term), ]

        # Calculate signed -log10(FDR) ? Or just -log10?
        # User asked for "-log FDR".
        # Usually we want to know direction too? But enrichment results (ORA) don't always give direction unless we add it.
        # ORA results usually just have pvalue.
        # Let's use simple -log10(padj).

        get_score <- function(df) {
            padj <- df$padj
            # Handle zero p-values (replace with min non-zero or epsilon)
            min_nz <- min(padj[padj > 0], na.rm = TRUE)
            padj[padj == 0] <- min_nz / 10
            -log10(padj)
        }

        score1 <- get_score(df1)
        score2 <- get_score(df2)

        plot_df <- data.frame(
            term = common_terms,
            x = score1,
            y = score2,
            stringsAsFactors = FALSE
        )

        # Determine significance status for coloring
        # Sig if padj < threshold -> score > -log10(threshold)
        cut_score <- -log10(p_thresh)

        plot_df$status <- "Not Sig"
        plot_df$status[plot_df$x > cut_score & plot_df$y > cut_score] <- "Both Sig"
        plot_df$status[plot_df$x > cut_score & plot_df$y <= cut_score] <- paste0(omic1, " Sig")
        plot_df$status[plot_df$x <= cut_score & plot_df$y > cut_score] <- paste0(omic2, " Sig")

        # Calculate correlation
        cor_res <- cor.test(plot_df$x, plot_df$y, method = corr_method)
        cor_val <- round(cor_res$estimate, 3)
        p_val <- signif(cor_res$p.value, 3)

        # Plot
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = status, label = term)) +
            ggplot2::geom_point(alpha = 0.7, size = 2) +
            ggplot2::geom_hline(yintercept = cut_score, linetype = "dashed", color = "gray50") +
            ggplot2::geom_vline(xintercept = cut_score, linetype = "dashed", color = "gray50") +
            ggplot2::labs(
                title = paste0("MultiGSEA: ", omic1, " vs ", omic2),
                subtitle = paste0(corr_method, " cor = ", cor_val, ", p = ", p_val),
                x = paste0(omic1, " -log10(FDR)"),
                y = paste0(omic2, " -log10(FDR)")
            ) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = c(
                "Both Sig" = "red",
                "Not Sig" = "grey",
                setNames("blue", paste0(omic1, " Sig")),
                setNames("green", paste0(omic2, " Sig"))
            ))

        # Add labels for top points (e.g., both sig)
        if (requireNamespace("ggrepel", quietly = TRUE)) {
            label_df <- plot_df[plot_df$status == "Both Sig", ]
            if (nrow(label_df) > 0) {
                # Limit to top 10 by sum of scores
                label_df$sum_score <- label_df$x + label_df$y
                label_df <- label_df[order(label_df$sum_score, decreasing = TRUE), ]
                label_df <- head(label_df, 10)

                p <- p + ggrepel::geom_text_repel(
                    data = label_df,
                    ggplot2::aes(label = term),
                    size = 3,
                    max.overlaps = 10
                )
            }
        }

        # Save
        filename <- paste0("multigsea_", omic1, "_vs_", omic2)
        save_plot(p, filename, config, width = 8, height = 8, custom_path = out_dir)
        save_table(plot_df, paste0(filename, ".csv"), config) # Save underlying data

        plots[[filename]] <- p
    }

    log_message("MultiGSEA plots generated: ", length(plots))
    return(plots)
}
