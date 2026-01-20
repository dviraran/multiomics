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

#' Run Pathview Visualization for Agreed Pathways
#'
#' @param enrichment_results List containing enrichment results.
#' @param mae_data MultiAssayExperiment data object (or list with harmonized_omics).
#' @param config Pipeline configuration list.
#'
#' @return List of generated plot paths.
#' @export
run_multigsea_pathview <- function(enrichment_results, mae_data, config) {
    log_message("=== Running MultiGSEA Pathview Visualization ===")

    if (!requireNamespace("pathview", quietly = TRUE)) {
        log_message("Package 'pathview' not installed. Skipping.")
        return(NULL)
    }

    mg_config <- config$enrichment$multigsea %||% list()
    if (!(mg_config$run_pathview %||% TRUE)) {
        log_message("Pathview analysis disabled in config.")
        return(NULL)
    }

    if (is.null(enrichment_results) || is.null(enrichment_results$per_omics)) {
        log_message("No enrichment results available.")
        return(NULL)
    }

    # Output directory
    out_dir <- file.path(config$output$output_dir, "enrichment", "multigsea", "pathview")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    # 1. Identify Agreed KEGG Pathways
    # We look for terms that are "KEGG" type and appear in >= 2 omics

    per_omics <- enrichment_results$per_omics
    omics_names <- names(per_omics)

    kegg_pathways <- list()

    for (omic in omics_names) {
        res <- per_omics[[omic]]$results
        if (!is.null(res) && nrow(res) > 0) {
            # Check for KEGG results (ID starts with hsa or purely numeric, or type column)
            # clusterProfiler KEGG IDs are usually "hsa12345"
            is_kegg <- FALSE
            if ("type" %in% colnames(res)) {
                kegg_res <- res[res$type == "KEGG", ]
            } else {
                # Heuristic: IDs start with "hsa" or "map" or numeric
                kegg_res <- res[grep("^hsa|^map|^[0-9]+$", res$ID), ]
            }

            if (nrow(kegg_res) > 0) {
                kegg_pathways[[omic]] <- kegg_res$ID
            }
        }
    }

    if (length(kegg_pathways) < 2) {
        log_message("Less than 2 omics have KEGG results. Skipping Pathview.")
        return(NULL)
    }

    # Find common pathways
    all_kegg <- unlist(kegg_pathways)
    if (length(all_kegg) == 0) {
        log_message("No KEGG pathways found.")
        return(NULL)
    }

    pathway_counts <- table(all_kegg)
    common_pathways <- names(pathway_counts)[pathway_counts >= 2]

    if (length(common_pathways) == 0) {
        log_message("No agreed KEGG pathways found between omics.")
        return(NULL)
    }

    log_message("Found ", length(common_pathways), " agreed KEGG pathways.")

    # 2. Prepare Data for Pathview
    # Gene Data (Transcriptomics + Proteomics)
    gene_data <- NULL

    # helper to get fold changes
    get_logfc <- function(omic_name, id_col = "entrez_id") {
        if (!omic_name %in% names(mae_data$harmonized_omics)) {
            return(NULL)
        }

        dat <- mae_data$harmonized_omics[[omic_name]]
        de <- dat$de_table %||% dat$da_table
        anno <- dat$feature_annotation

        if (is.null(de) || is.null(anno)) {
            return(NULL)
        }

        # Merge to get IDs
        # feature_id is common
        merged <- merge(de, anno, by = "feature_id")

        if (!id_col %in% colnames(merged)) {
            return(NULL)
        }

        # Get LogFC column
        fc_col <- grep("logFC|log2FoldChange", colnames(merged), ignore.case = TRUE, value = TRUE)[1]
        if (is.na(fc_col)) {
            return(NULL)
        }

        # Create named vector
        # Handle multiple features mapping to same Entrez ID: take mean
        vec <- tapply(merged[[fc_col]], merged[[id_col]], mean, na.rm = TRUE)
        return(vec)
    }

    # Transcriptomics
    rna_fc <- get_logfc("transcriptomics", "entrez_id")

    # Proteomics
    prot_fc <- get_logfc("proteomics", "entrez_id")

    # Combine Gene Data
    if (!is.null(rna_fc) && !is.null(prot_fc)) {
        # Create matrix
        all_genes <- unique(c(names(rna_fc), names(prot_fc)))
        gene_data <- matrix(NA,
            nrow = length(all_genes), ncol = 2,
            dimnames = list(all_genes, c("Transcriptomics", "Proteomics"))
        )

        idx_rna <- match(names(rna_fc), all_genes)
        gene_data[idx_rna, 1] <- rna_fc

        idx_prot <- match(names(prot_fc), all_genes)
        gene_data[idx_prot, 2] <- prot_fc
    } else if (!is.null(rna_fc)) {
        gene_data <- rna_fc
    } else if (!is.null(prot_fc)) {
        gene_data <- prot_fc
    }

    # Metabolomics Data (CPD Data)
    cpd_data <- NULL
    met_fc <- get_logfc("metabolomics", "kegg_id")
    if (!is.null(met_fc)) {
        cpd_data <- met_fc
    }

    if (is.null(gene_data) && is.null(cpd_data)) {
        log_message("No valid Entrez/KEGG IDs found in data for Pathview.")
        return(NULL)
    }

    # 3. Run Pathview
    generated_plots <- list()

    cwd <- getwd()
    setwd(out_dir)
    on.exit(setwd(cwd))

    for (pid in common_pathways) {
        # Clean ID
        clean_pid <- sub("^[a-z]+", "", pid)

        tryCatch(
            {
                pv.out <- pathview::pathview(
                    gene.data = gene_data,
                    cpd.data = cpd_data,
                    pathway.id = clean_pid,
                    species = "hsa",
                    out.suffix = "multiomics",
                    temp.file = TRUE,
                    kegg.dir = out_dir,
                    keys.align = "y",
                    kev.dir = NULL,
                    match.data = TRUE,
                    multi.state = !is.null(dim(gene_data)) && ncol(gene_data) > 1,
                    same.layer = FALSE
                )

                outfile <- paste0("hsa", clean_pid, ".multiomics.png")
                if (file.exists(outfile)) {
                    generated_plots[[pid]] <- file.path(out_dir, outfile)
                    log_message("Generated Pathview: ", outfile)
                }
            },
            error = function(e) {
                log_message("Pathview failed for ", pid, ": ", e$message)
            }
        )
    }

    return(generated_plots)
}
