#' Executive Summary Generation for Pipeline Reports
#'
#' This module generates executive summaries for analysis reports,
#' providing key findings in a concise format at the top of reports.
#'
#' @description
#' Features:
#' - Automatic extraction of key findings from results
#' - Configurable summary length and focus areas
#' - Support for all pipeline types

# ==============================================================================
# MAIN SUMMARY FUNCTIONS
# ==============================================================================

#' Generate executive summary for RNA-seq analysis
#'
#' @param de_results Differential expression results
#' @param enrichment_results Pathway enrichment results
#' @param qc_results QC analysis results
#' @param config Configuration list
#' @return Character string with formatted summary
#' @export
generate_rnaseq_summary <- function(de_results,
                                     enrichment_results = NULL,
                                     qc_results = NULL,
                                     config = NULL) {

    summary_points <- character()

    # Sample information
    if (!is.null(qc_results)) {
        n_samples <- qc_results$n_samples %||% NA
        n_genes <- qc_results$n_genes %||% NA

        if (!is.na(n_samples) && !is.na(n_genes)) {
            summary_points <- c(summary_points,
                sprintf("Analyzed %d samples with %s genes after filtering",
                        n_samples, format(n_genes, big.mark = ",")))
        }

        # Outliers
        if (!is.null(qc_results$outliers)) {
            n_outliers <- sum(qc_results$outliers$is_outlier, na.rm = TRUE)
            if (n_outliers > 0) {
                summary_points <- c(summary_points,
                    sprintf("Detected %d potential outlier sample(s) based on PCA", n_outliers))
            }
        }
    }

    # DE results summary
    if (!is.null(de_results) && length(de_results) > 0) {
        for (contrast_name in names(de_results)) {
            de_table <- de_results[[contrast_name]]

            if (is.data.frame(de_table)) {
                # Count significant genes
                sig_col <- intersect(c("padj", "adj.P.Val", "FDR", "adj_pvalue"),
                                     names(de_table))[1]
                lfc_col <- intersect(c("log2FoldChange", "logFC", "log2FC"),
                                     names(de_table))[1]

                if (!is.na(sig_col) && !is.na(lfc_col)) {
                    alpha <- config$alpha %||% 0.05
                    lfc_thresh <- config$lfc_threshold %||% 0

                    sig_genes <- de_table[[sig_col]] < alpha & !is.na(de_table[[sig_col]])
                    if (lfc_thresh > 0) {
                        sig_genes <- sig_genes & abs(de_table[[lfc_col]]) > lfc_thresh
                    }

                    n_sig <- sum(sig_genes, na.rm = TRUE)
                    n_up <- sum(sig_genes & de_table[[lfc_col]] > 0, na.rm = TRUE)
                    n_down <- sum(sig_genes & de_table[[lfc_col]] < 0, na.rm = TRUE)

                    if (n_sig > 0) {
                        summary_points <- c(summary_points,
                            sprintf("%s: %d significant genes (%d up, %d down) at FDR < %.2f",
                                    contrast_name, n_sig, n_up, n_down, alpha))

                        # Top genes
                        top_genes <- get_top_genes(de_table, sig_col, lfc_col, n = 5)
                        if (length(top_genes) > 0) {
                            summary_points <- c(summary_points,
                                sprintf("  Top genes: %s", paste(top_genes, collapse = ", ")))
                        }
                    } else {
                        summary_points <- c(summary_points,
                            sprintf("%s: No significant genes at FDR < %.2f",
                                    contrast_name, alpha))
                    }
                }
            }
        }
    }

    # Enrichment summary
    if (!is.null(enrichment_results) && length(enrichment_results) > 0) {
        top_pathways <- get_top_pathways(enrichment_results, n = 3)
        if (length(top_pathways) > 0) {
            summary_points <- c(summary_points,
                sprintf("Top enriched pathways: %s", paste(top_pathways, collapse = "; ")))
        }
    }

    format_summary(summary_points, "RNA-seq Analysis Summary")
}


#' Generate executive summary for proteomics analysis
#'
#' @param da_results Differential abundance results
#' @param enrichment_results Enrichment results
#' @param qc_results QC results
#' @param config Configuration list
#' @return Character string with formatted summary
#' @export
generate_proteomics_summary <- function(da_results,
                                         enrichment_results = NULL,
                                         qc_results = NULL,
                                         config = NULL) {

    summary_points <- character()

    # Sample/protein info
    if (!is.null(qc_results)) {
        n_samples <- qc_results$n_samples %||% NA
        n_proteins <- qc_results$n_proteins %||% NA
        pct_missing <- qc_results$pct_missing %||% NA

        if (!is.na(n_samples) && !is.na(n_proteins)) {
            summary_points <- c(summary_points,
                sprintf("Analyzed %d samples with %s proteins",
                        n_samples, format(n_proteins, big.mark = ",")))
        }

        if (!is.na(pct_missing)) {
            summary_points <- c(summary_points,
                sprintf("Missing value rate: %.1f%% (imputed using %s)",
                        pct_missing, config$imputation$method %||% "default method"))
        }
    }

    # DA results
    if (!is.null(da_results) && length(da_results) > 0) {
        for (contrast_name in names(da_results)) {
            da_table <- da_results[[contrast_name]]

            if (is.data.frame(da_table)) {
                sig_col <- intersect(c("adj.P.Val", "padj", "FDR"), names(da_table))[1]
                lfc_col <- intersect(c("logFC", "log2FoldChange"), names(da_table))[1]

                if (!is.na(sig_col) && !is.na(lfc_col)) {
                    alpha <- config$differential$adj_pvalue_threshold %||% 0.05
                    lfc_thresh <- config$differential$log2fc_threshold %||% 0

                    sig <- da_table[[sig_col]] < alpha & !is.na(da_table[[sig_col]])
                    if (lfc_thresh > 0) {
                        sig <- sig & abs(da_table[[lfc_col]]) > lfc_thresh
                    }

                    n_sig <- sum(sig, na.rm = TRUE)
                    n_up <- sum(sig & da_table[[lfc_col]] > 0, na.rm = TRUE)
                    n_down <- sum(sig & da_table[[lfc_col]] < 0, na.rm = TRUE)

                    summary_points <- c(summary_points,
                        sprintf("%s: %d significant proteins (%d up, %d down)",
                                contrast_name, n_sig, n_up, n_down))
                }
            }
        }
    }

    format_summary(summary_points, "Proteomics Analysis Summary")
}


#' Generate executive summary for metabolomics analysis
#'
#' @param de_results Differential analysis results
#' @param enrichment_results Enrichment results
#' @param qc_results QC results
#' @param config Configuration list
#' @return Character string with formatted summary
#' @export
generate_metabolomics_summary <- function(de_results,
                                           enrichment_results = NULL,
                                           qc_results = NULL,
                                           config = NULL) {

    summary_points <- character()

    # Sample/feature info
    if (!is.null(qc_results)) {
        n_samples <- qc_results$n_samples %||% NA
        n_features <- qc_results$n_features %||% NA
        data_type <- config$input$data_type %||% "untargeted"

        if (!is.na(n_samples) && !is.na(n_features)) {
            summary_points <- c(summary_points,
                sprintf("Analyzed %d samples with %s %s metabolite features",
                        n_samples, format(n_features, big.mark = ","), data_type))
        }

        # Batch correction
        if (!is.null(qc_results$batch_correction_applied) && qc_results$batch_correction_applied) {
            summary_points <- c(summary_points,
                sprintf("Batch correction applied using %s method",
                        config$batch_correction$method %||% "standard"))
        }
    }

    # DE results
    if (!is.null(de_results) && !is.null(de_results$results)) {
        for (contrast_name in names(de_results$results)) {
            de_table <- de_results$results[[contrast_name]]$table

            if (is.data.frame(de_table)) {
                sig_col <- intersect(c("adj.P.Val", "padj", "FDR"), names(de_table))[1]

                if (!is.na(sig_col)) {
                    alpha <- 0.05
                    n_sig <- sum(de_table[[sig_col]] < alpha, na.rm = TRUE)

                    summary_points <- c(summary_points,
                        sprintf("%s: %d significantly altered metabolites (FDR < %.2f)",
                                contrast_name, n_sig, alpha))
                }
            }
        }
    }

    format_summary(summary_points, "Metabolomics Analysis Summary")
}


#' Generate executive summary for multi-omics analysis
#'
#' @param mofa_results MOFA2 results
#' @param diablo_results DIABLO results
#' @param concordance_results Concordance analysis results
#' @param config Configuration list
#' @return Character string with formatted summary
#' @export
generate_multiomics_summary <- function(mofa_results = NULL,
                                         diablo_results = NULL,
                                         concordance_results = NULL,
                                         config = NULL) {

    summary_points <- character()

    omics_present <- config$omics_present %||% c()
    summary_points <- c(summary_points,
        sprintf("Integrated %d omics types: %s",
                length(omics_present), paste(omics_present, collapse = ", ")))

    # MOFA results
    if (!is.null(mofa_results)) {
        n_factors <- mofa_results$n_factors %||% NA
        var_explained <- mofa_results$total_variance_explained %||% NA

        if (!is.na(n_factors)) {
            summary_points <- c(summary_points,
                sprintf("MOFA2: Identified %d latent factors", n_factors))

            if (!is.na(var_explained)) {
                summary_points <- c(summary_points,
                    sprintf("  Total variance explained: %.1f%%", var_explained * 100))
            }

            # Top contributing omics per factor
            if (!is.null(mofa_results$factor_contributions)) {
                top_factor <- mofa_results$factor_contributions[[1]]
                if (!is.null(top_factor)) {
                    summary_points <- c(summary_points,
                        sprintf("  Factor 1 top contributors: %s",
                                paste(names(head(sort(top_factor, decreasing = TRUE), 2)),
                                      collapse = ", ")))
                }
            }
        }
    }

    # DIABLO results
    if (!is.null(diablo_results)) {
        n_components <- diablo_results$n_components %||% NA
        error_rate <- diablo_results$error_rate %||% NA

        if (!is.na(n_components)) {
            summary_points <- c(summary_points,
                sprintf("DIABLO: %d components for supervised integration", n_components))

            if (!is.na(error_rate)) {
                summary_points <- c(summary_points,
                    sprintf("  Cross-validation error rate: %.1f%%", error_rate * 100))
            }
        }

        # Top features
        if (!is.null(diablo_results$top_features)) {
            for (omics in names(diablo_results$top_features)) {
                top_feats <- head(diablo_results$top_features[[omics]], 3)
                if (length(top_feats) > 0) {
                    summary_points <- c(summary_points,
                        sprintf("  Top %s features: %s", omics, paste(top_feats, collapse = ", ")))
                }
            }
        }
    }

    # Concordance
    if (!is.null(concordance_results)) {
        if (!is.null(concordance_results$method_agreement)) {
            agreement <- concordance_results$method_agreement
            summary_points <- c(summary_points,
                sprintf("Method concordance: MOFA-DIABLO agreement = %.1f%%",
                        agreement * 100))
        }
    }

    format_summary(summary_points, "Multi-Omics Integration Summary")
}


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Get top genes from DE results
get_top_genes <- function(de_table, sig_col, lfc_col, n = 5) {
    # Try to get gene symbol column
    symbol_col <- intersect(c("symbol", "gene_symbol", "gene_name", "SYMBOL"),
                            names(de_table))[1]

    if (is.na(symbol_col)) {
        # Fall back to row names or first column
        if (!is.null(rownames(de_table))) {
            gene_ids <- rownames(de_table)
        } else {
            gene_ids <- de_table[[1]]
        }
    } else {
        gene_ids <- de_table[[symbol_col]]
    }

    # Sort by significance and fold change
    sig_vals <- de_table[[sig_col]]
    lfc_vals <- abs(de_table[[lfc_col]])

    # Combined ranking
    ranks <- rank(sig_vals, na.last = TRUE) - rank(lfc_vals, na.last = TRUE)
    top_idx <- head(order(ranks), n)

    gene_ids[top_idx]
}


#' Get top pathways from enrichment results
get_top_pathways <- function(enrichment_results, n = 3) {
    pathways <- character()

    # Handle different result structures
    if (is.data.frame(enrichment_results)) {
        # Single data frame
        name_col <- intersect(c("pathway", "Description", "term", "name"),
                              names(enrichment_results))[1]
        if (!is.na(name_col)) {
            pathways <- head(enrichment_results[[name_col]], n)
        }
    } else if (is.list(enrichment_results)) {
        # List of results (per contrast or per database)
        for (result in enrichment_results) {
            if (is.data.frame(result)) {
                name_col <- intersect(c("pathway", "Description", "term", "name"),
                                      names(result))[1]
                if (!is.na(name_col)) {
                    pathways <- c(pathways, head(result[[name_col]], 2))
                }
            }
        }
        pathways <- head(unique(pathways), n)
    }

    # Clean up pathway names (remove prefix codes, truncate)
    pathways <- gsub("^GO:[0-9]+\\s*", "", pathways)
    pathways <- sapply(pathways, function(x) {
        if (nchar(x) > 50) paste0(substr(x, 1, 47), "...") else x
    })

    pathways
}


#' Format summary points into a nice output
format_summary <- function(summary_points, title = "Executive Summary") {
    if (length(summary_points) == 0) {
        return(paste0("## ", title, "\n\nNo summary available.\n"))
    }

    # Format as markdown
    formatted <- paste0("## ", title, "\n\n")

    for (point in summary_points) {
        if (startsWith(point, "  ")) {
            # Sub-point (indented)
            formatted <- paste0(formatted, "  - ", trimws(point), "\n")
        } else {
            formatted <- paste0(formatted, "- ", point, "\n")
        }
    }

    formatted
}


#' Generate summary for any pipeline type
#'
#' @param pipeline_type One of "rnaseq", "proteomics", "metabolomics", "multiomics"
#' @param results Named list of results from the pipeline
#' @param config Configuration list
#' @return Character string with formatted summary
#' @export
generate_summary <- function(pipeline_type, results, config = NULL) {
    switch(pipeline_type,
        "rnaseq" = generate_rnaseq_summary(
            de_results = results$de_results,
            enrichment_results = results$enrichment_results,
            qc_results = results$qc_results,
            config = config
        ),
        "proteomics" = generate_proteomics_summary(
            da_results = results$da_results,
            enrichment_results = results$enrichment_results,
            qc_results = results$qc_results,
            config = config
        ),
        "metabolomics" = generate_metabolomics_summary(
            de_results = results$de_results,
            enrichment_results = results$enrichment_results,
            qc_results = results$qc_results,
            config = config
        ),
        "multiomics" = generate_multiomics_summary(
            mofa_results = results$mofa_results,
            diablo_results = results$diablo_results,
            concordance_results = results$concordance_results,
            config = config
        ),
        paste("Unknown pipeline type:", pipeline_type)
    )
}


# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
