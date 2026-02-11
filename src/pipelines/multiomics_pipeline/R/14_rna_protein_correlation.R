#' Compute RNA-Protein log2FC correlation
#'
#' This script defines a function that calculates Pearson correlation between
#' log2 fold‑change values of transcriptomics (RNA) and proteomics (protein)
#' for the common genes. It saves a CSV table of the correlations and a
#' histogram plot.
#'
#' @param mae_data MultiAssayExperiment object containing harmonized omics data
#' @param config List of configuration parameters (used for output paths)
#' @return Data frame with gene symbols and correlation values
#' @export
rna_protein_correlation <- function(mae_data, config) {
    # Load required libraries
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required for plotting")
    }

    # Helper to create output dirs
    create_dir <- function(path) {
        if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
    }

    # =========================================================================
    # Part 1: Expression Correlation (Abundance vs Abundance per gene)
    # =========================================================================

    # Extract normalized matrices
    rna_mat <- mae_data$harmonized_omics$transcriptomics$normalized_matrix
    prot_mat <- mae_data$harmonized_omics$proteomics$normalized_matrix

    # Gene mapping must be available
    gene_mapping <- mae_data$gene_mapping

    cor_df <- NULL
    summary_stats <- list(
        n_genes = 0,
        mean_cor = NA,
        median_cor = NA,
        pct_positive = NA,
        n_significant = 0
    )

    if (!is.null(gene_mapping) && !is.null(rna_mat) && !is.null(prot_mat)) {
        # Subset mapping for each omic
        rna_genes <- gene_mapping[gene_mapping$omics == "transcriptomics", ]
        prot_genes <- gene_mapping[gene_mapping$omics == "proteomics", ]

        # Find common gene symbols
        common_genes <- intersect(rna_genes$gene_symbol, prot_genes$gene_symbol)
        common_genes <- common_genes[!is.na(common_genes)]

        if (length(common_genes) >= 10) {
            # Map to feature IDs
            rna_features <- rna_genes$feature_id[match(common_genes, rna_genes$gene_symbol)]
            prot_features <- prot_genes$feature_id[match(common_genes, prot_genes$gene_symbol)]

            # Align samples
            common_samples <- intersect(colnames(rna_mat), colnames(prot_mat))

            if (length(common_samples) > 2) {
                # Compute per‑gene Pearson correlation
                correlations <- numeric(length(common_genes))
                names(correlations) <- common_genes

                for (i in seq_along(common_genes)) {
                    rna_expr <- rna_mat[rna_features[i], common_samples]
                    prot_expr <- prot_mat[prot_features[i], common_samples]
                    if (sd(rna_expr, na.rm = TRUE) > 0 && sd(prot_expr, na.rm = TRUE) > 0) {
                        correlations[i] <- cor(rna_expr, prot_expr, use = "pairwise.complete.obs")
                    } else {
                        correlations[i] <- NA
                    }
                }
                correlations <- correlations[!is.na(correlations)]

                if (length(correlations) > 0) {
                    # Save results
                    cor_df <- data.frame(
                        gene_symbol = names(correlations),
                        correlation = correlations,
                        stringsAsFactors = FALSE
                    )

                    csv_path <- file.path(config$output$output_dir, "tables", "rna_protein_correlations.csv")
                    create_dir(csv_path)
                    utils::write.csv(cor_df, csv_path, row.names = FALSE)

                    # Plot histogram
                    plot_path <- file.path(config$output$output_dir, "plots", "rna_protein_concordance.png")
                    create_dir(plot_path)

                    p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = correlation)) +
                        ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
                        ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
                        ggplot2::theme_minimal() +
                        ggplot2::labs(
                            title = "RNA‑Protein Expression Correlation Distribution",
                            x = "Pearson Correlation",
                            y = "Count"
                        )
                    ggplot2::ggsave(plot_path, p, width = 8, height = 6, dpi = 150)

                    # Update summary
                    summary_stats$n_genes <- length(correlations)
                    summary_stats$mean_cor <- mean(correlations)
                    summary_stats$median_cor <- median(correlations)
                    summary_stats$pct_positive <- 100 * mean(correlations > 0)
                    summary_stats$n_significant <- sum(abs(correlations) > 0.5)
                }
            }
        }
    }

    # =========================================================================
    # Part 2: Differential Concordance (Log2FC vs Log2FC)
    # =========================================================================

    de_concordance_df <- NULL

    # Extract DE/DA tables
    rna_obj <- mae_data$harmonized_omics$transcriptomics
    prot_obj <- mae_data$harmonized_omics$proteomics

    # We prefer the 'de_table' or 'da_table' attached to the objects
    rna_de <- if (!is.null(rna_obj$de_table)) rna_obj$de_table else rna_obj$da_table
    prot_de <- if (!is.null(prot_obj$da_table)) prot_obj$da_table else prot_obj$de_table

    # Helper to get padj col
    get_padj_col <- function(df) {
        if ("adj.P.Val" %in% colnames(df)) {
            return("adj.P.Val")
        }
        if ("padj" %in% colnames(df)) {
            return("padj")
        }
        if ("FDR" %in% colnames(df)) {
            return("FDR")
        }
        return(NULL)
    }

    if (!is.null(rna_de) && !is.null(prot_de) && !is.null(gene_mapping)) {
        # We need to map protein feature_ids to gene symbols, and RNA feature_ids to gene symbols
        # to join them.

        # Prepare RNA
        # RNA often has gene_id as feature_id. Check mapping.
        rna_map_sub <- gene_mapping[gene_mapping$omics == "transcriptomics", ]

        # Ensure 'log2FC' exists
        if ("log2FC" %in% colnames(rna_de)) {
            rna_de$gene_symbol <- rna_map_sub$gene_symbol[match(rna_de$feature_id, rna_map_sub$feature_id)]

            padj_col <- get_padj_col(rna_de)
            cols_to_keep <- c("gene_symbol", "log2FC")
            if (!is.null(padj_col)) cols_to_keep <- c(cols_to_keep, padj_col)

            rna_de_clean <- rna_de[!is.na(rna_de$gene_symbol), cols_to_keep, drop = FALSE]

            # Normalize names
            if (!is.null(padj_col)) {
                colnames(rna_de_clean) <- c("gene_symbol", "rna_log2FC", "rna_padj")
            } else {
                colnames(rna_de_clean) <- c("gene_symbol", "rna_log2FC")
                rna_de_clean$rna_padj <- NA
            }
        } else {
            rna_de_clean <- NULL
        }

        # Prepare Protein
        prot_map_sub <- gene_mapping[gene_mapping$omics == "proteomics", ]
        if ("log2FC" %in% colnames(prot_de)) {
            prot_de$gene_symbol <- prot_map_sub$gene_symbol[match(prot_de$feature_id, prot_map_sub$feature_id)]

            padj_col <- get_padj_col(prot_de)
            cols_to_keep <- c("gene_symbol", "log2FC")
            if (!is.null(padj_col)) cols_to_keep <- c(cols_to_keep, padj_col)

            prot_de_clean <- prot_de[!is.na(prot_de$gene_symbol), cols_to_keep, drop = FALSE]

            # Normalize names
            if (!is.null(padj_col)) {
                colnames(prot_de_clean) <- c("gene_symbol", "protein_log2FC", "protein_padj")
            } else {
                colnames(prot_de_clean) <- c("gene_symbol", "protein_log2FC")
                prot_de_clean$protein_padj <- NA
            }
        } else {
            prot_de_clean <- NULL
        }

        if (!is.null(rna_de_clean) && !is.null(prot_de_clean)) {
            # Merge
            de_merged <- merge(rna_de_clean, prot_de_clean, by = "gene_symbol")

            if (nrow(de_merged) > 10) {
                # 1. Define Significance Categories (FDR < 0.05)
                # Handle NAs in p-values by treating them as 1 (non-significant)
                rna_p <- de_merged$rna_padj
                rna_p[is.na(rna_p)] <- 1

                prot_p <- de_merged$protein_padj
                prot_p[is.na(prot_p)] <- 1

                is_sig_rna <- rna_p < 0.05
                is_sig_prot <- prot_p < 0.05

                de_merged$category <- "Non-sig"
                # Order of assignment matters for priority if needed, but mutually exclusive here except 'Both'
                de_merged$category[is_sig_rna & !is_sig_prot] <- "Sig RNA (Gold)"
                de_merged$category[!is_sig_rna & is_sig_prot] <- "Sig Protein (Purple)"
                de_merged$category[is_sig_rna & is_sig_prot] <- "Sig Both (Red)"

                # Set factor levels for plotting order (Non-sig on bottom, Both on top)
                de_merged$category <- factor(de_merged$category,
                    levels = c("Non-sig", "Sig RNA (Gold)", "Sig Protein (Purple)", "Sig Both (Red)")
                )

                de_merged$concordant <- sign(de_merged$rna_log2FC) == sign(de_merged$protein_log2FC)

                # 2. Calculate Correlations
                # Union of significant genes (Sig in RNA OR Sig in Protein OR Both)
                is_sig_union <- is_sig_rna | is_sig_prot

                cor_all <- cor(de_merged$rna_log2FC, de_merged$protein_log2FC, use = "complete.obs")

                if (sum(is_sig_union) > 3) {
                    cor_sig <- cor(de_merged$rna_log2FC[is_sig_union], de_merged$protein_log2FC[is_sig_union], use = "complete.obs")
                } else {
                    cor_sig <- NA
                }

                # Save table
                csv_path <- file.path(config$output$output_dir, "tables", "rna_protein_de_concordance.csv")
                create_dir(csv_path)
                utils::write.csv(de_merged, csv_path, row.names = FALSE)
                de_concordance_df <- de_merged

                # 3. Create Scatter Plot
                plot_path <- file.path(config$output$output_dir, "plots", "rna_protein_de_scatter.png")
                create_dir(plot_path)

                # Define requested colors
                # Gold, Purple, Red, Gray
                custom_colors <- c(
                    "Non-sig" = "gray80",
                    "Sig RNA (Gold)" = "gold3",
                    "Sig Protein (Purple)" = "purple",
                    "Sig Both (Red)" = "red"
                )

                subtitle_text <- paste0("All genes: r = ", round(cor_all, 3))
                if (!is.na(cor_sig)) {
                    subtitle_text <- paste0(subtitle_text, " | Sig. Union (FDR<0.05): r = ", round(cor_sig, 3))
                }

                p <- ggplot2::ggplot(de_merged, ggplot2::aes(x = rna_log2FC, y = protein_log2FC)) +
                    ggplot2::geom_vline(xintercept = 0, color = "gray90") +
                    ggplot2::geom_hline(yintercept = 0, color = "gray90") +
                    # Plot non-sig first (bottom layer), then others
                    ggplot2::geom_point(ggplot2::aes(color = category), alpha = 0.6, size = 1.5) +
                    ggplot2::geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE, linewidth = 0.5) +
                    ggplot2::scale_color_manual(values = custom_colors) +
                    ggplot2::theme_minimal() +
                    ggplot2::labs(
                        title = "Differential Concordance: RNA vs Protein",
                        subtitle = subtitle_text,
                        x = "RNA log2 Fold Change",
                        y = "Protein log2 Fold Change",
                        color = "Significance"
                    ) +
                    ggplot2::theme(legend.position = "bottom")

                ggplot2::ggsave(plot_path, p, width = 8, height = 7, dpi = 300)

                # 4. Translation Efficiency (TE) Analysis
                # TE ~ Protein / RNA => log2(TE) = log2(Protein) - log2(RNA)
                # We interpret log2FC differences as changes in TE
                de_merged$te_log2FC <- de_merged$protein_log2FC - de_merged$rna_log2FC

                # Save TE results
                te_path <- file.path(config$output$output_dir, "tables", "translation_efficiency.csv")
                utils::write.csv(de_merged[, c("gene_symbol", "rna_log2FC", "protein_log2FC", "te_log2FC", "category")],
                    te_path,
                    row.names = FALSE
                )

                # Plot 1: TE Distribution
                te_hist_path <- file.path(config$output$output_dir, "plots", "translation_efficiency_hist.png")
                p_hist <- ggplot2::ggplot(de_merged, ggplot2::aes(x = te_log2FC)) +
                    ggplot2::geom_histogram(bins = 40, fill = "darkcyan", color = "white", alpha = 0.8) +
                    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                    ggplot2::theme_minimal() +
                    ggplot2::labs(
                        title = "Distribution of Translation Efficiency Changes",
                        subtitle = "log2(TE) = log2(Protein FC) - log2(RNA FC)",
                        x = "TE log2 Fold Change",
                        y = "Count"
                    )
                ggplot2::ggsave(te_hist_path, p_hist, width = 8, height = 6)

                # Plot 2: TE Scatter Plot (RNA log2FC vs Protein log2FC, colored by TE)
                # Points above identity line = higher protein FC => positive TE (red)
                # Points below identity line = lower protein FC => negative TE (blue)
                te_scatter_path <- file.path(config$output$output_dir, "plots", "translation_efficiency_scatter.png")

                te_cor <- cor(de_merged$rna_log2FC, de_merged$protein_log2FC, use = "complete.obs")
                te_subtitle <- paste0(
                    "r = ", round(te_cor, 3),
                    " | Red = protein > RNA (high TE), Blue = protein < RNA (low TE)"
                )

                p_scatter <- ggplot2::ggplot(de_merged, ggplot2::aes(x = rna_log2FC, y = protein_log2FC, color = te_log2FC)) +
                    ggplot2::geom_hline(yintercept = 0, color = "gray90") +
                    ggplot2::geom_vline(xintercept = 0, color = "gray90") +
                    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
                    ggplot2::geom_point(alpha = 0.7, size = 1.5) +
                    ggplot2::scale_color_gradient2(
                        low = "blue", mid = "grey90", high = "red",
                        midpoint = 0, name = "log2(TE)\n(Prot FC - RNA FC)"
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::labs(
                        title = "Translation Efficiency: RNA vs Protein log2FC",
                        subtitle = te_subtitle,
                        x = "RNA log2 Fold Change",
                        y = "Protein log2 Fold Change"
                    )
                ggplot2::ggsave(te_scatter_path, p_scatter, width = 8, height = 7, dpi = 300)
            }
        }
    }

    # Return nested list structure required by the report
    list(
        rna_protein = list(
            correlations = cor_df,
            summary = summary_stats
        ),
        de_concordance = de_concordance_df
    )
}
