# =============================================================================
# Plotting Functions for Multi-Omics Viewer
# =============================================================================
# Reusable ggplot2/plotly functions for common visualizations

# =============================================================================
# PCA Plots
# =============================================================================

#' Create interactive PCA scatter plot
#' @param pca_data Data frame with PC columns and sample info
#' @param pc_x PC for x-axis (e.g., "PC1")
#' @param pc_y PC for y-axis (e.g., "PC2")
#' @param color_var Variable to color by
#' @param var_explained Named vector of variance explained per PC
#' @param show_labels Show sample labels
#' @param show_ellipse Show confidence ellipses
#' @return plotly object
create_pca_plot <- function(pca_data, pc_x = "PC1", pc_y = "PC2",
                            color_var = NULL, var_explained = NULL,
                            show_labels = FALSE, show_ellipse = TRUE) {

    # Base plot
    p <- ggplot(pca_data, aes_string(x = pc_x, y = pc_y))

    if (!is.null(color_var) && color_var %in% colnames(pca_data)) {
        p <- p + aes_string(color = color_var, text = "sample")
    } else {
        p <- p + aes(text = sample)
    }

    p <- p + geom_point(size = 3, alpha = 0.8)

    # Add ellipses if requested and color variable is set
    if (show_ellipse && !is.null(color_var) && color_var %in% colnames(pca_data)) {
        p <- p + stat_ellipse(level = 0.95, linetype = "dashed", show.legend = FALSE)
    }

    # Add labels if requested
    if (show_labels && "sample" %in% colnames(pca_data)) {
        p <- p + geom_text(aes(label = sample), size = 3, vjust = 1.5,
                           show.legend = FALSE)
    }

    # Axis labels with variance explained
    x_label <- pc_x
    y_label <- pc_y

    if (!is.null(var_explained)) {
        pc_x_idx <- as.integer(gsub("PC", "", pc_x))
        pc_y_idx <- as.integer(gsub("PC", "", pc_y))

        if (pc_x_idx <= length(var_explained)) {
            x_label <- paste0(pc_x, " (", round(var_explained[pc_x_idx], 1), "%)")
        }
        if (pc_y_idx <= length(var_explained)) {
            y_label <- paste0(pc_y, " (", round(var_explained[pc_y_idx], 1), "%)")
        }
    }

    p <- p +
        labs(x = x_label, y = y_label) +
        theme_minimal() +
        theme(legend.position = "right")

    ggplotly(p, tooltip = c("text", "color"))
}

#' Create scree plot for PCA variance explained
#' @param var_explained Numeric vector of variance explained
#' @param n_pcs Number of PCs to show
#' @return ggplot object
create_scree_plot <- function(var_explained, n_pcs = 10) {
    n_pcs <- min(n_pcs, length(var_explained))
    var_exp <- var_explained[1:n_pcs]

    df <- data.frame(
        PC = factor(paste0("PC", seq_along(var_exp)),
                    levels = paste0("PC", seq_along(var_exp))),
        variance = var_exp,
        cumulative = cumsum(var_exp)
    )

    ggplot(df, aes(x = PC, y = variance)) +
        geom_col(fill = "#3498db", alpha = 0.8) +
        geom_line(aes(y = cumulative, group = 1), color = "#e74c3c", linewidth = 1) +
        geom_point(aes(y = cumulative), color = "#e74c3c", size = 2) +
        labs(title = "Variance Explained by PC",
             x = NULL,
             y = "% Variance") +
        theme_minimal()
}

# =============================================================================
# Volcano Plots
# =============================================================================

#' Create interactive volcano plot
#' @param de_data Data frame with log2FC, pvalue, padj columns
#' @param log2fc_col Column name for log2 fold change
#' @param pval_col Column name for p-value (or adjusted p-value)
#' @param label_col Column name for labels (gene names)
#' @param log2fc_thresh Threshold for significance
#' @param pval_thresh P-value threshold
#' @param title Plot title
#' @return plotly object
create_volcano_plot <- function(de_data, log2fc_col = "log2FoldChange",
                                pval_col = "padj", label_col = "symbol",
                                log2fc_thresh = 1, pval_thresh = 0.05,
                                title = "Volcano Plot") {

    # Determine significance
    de_data$neg_log10_pval <- -log10(de_data[[pval_col]] + 1e-300)
    de_data$significance <- "ns"
    de_data$significance[de_data[[pval_col]] < pval_thresh &
                          de_data[[log2fc_col]] > log2fc_thresh] <- "up"
    de_data$significance[de_data[[pval_col]] < pval_thresh &
                          de_data[[log2fc_col]] < -log2fc_thresh] <- "down"
    de_data$significance <- factor(de_data$significance, levels = c("up", "down", "ns"))

    # Create hover text
    if (label_col %in% colnames(de_data)) {
        de_data$hover_text <- paste0(
            de_data[[label_col]], "<br>",
            "log2FC: ", round(de_data[[log2fc_col]], 2), "<br>",
            "p-adj: ", format(de_data[[pval_col]], digits = 3)
        )
    } else {
        de_data$hover_text <- paste0(
            "log2FC: ", round(de_data[[log2fc_col]], 2), "<br>",
            "p-adj: ", format(de_data[[pval_col]], digits = 3)
        )
    }

    # Color palette
    colors <- c("up" = "#e74c3c", "down" = "#3498db", "ns" = "#95a5a6")

    p <- ggplot(de_data, aes_string(x = log2fc_col, y = "neg_log10_pval",
                                     color = "significance", text = "hover_text")) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_color_manual(values = colors,
                           labels = c("up" = "Up-regulated",
                                      "down" = "Down-regulated",
                                      "ns" = "Not significant")) +
        geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = c(-log2fc_thresh, log2fc_thresh),
                   linetype = "dashed", color = "gray50") +
        labs(title = title,
             x = "log2 Fold Change",
             y = "-log10(adjusted p-value)",
             color = "Significance") +
        theme_minimal() +
        theme(legend.position = "right")

    ggplotly(p, tooltip = "text")
}

# =============================================================================
# Heatmaps
# =============================================================================

#' Create interactive heatmap
#' @param matrix_data Matrix or data frame (features as rows, samples as columns)
#' @param row_annotation Optional row annotation data frame
#' @param col_annotation Optional column annotation data frame
#' @param scale_rows Z-score scale rows
#' @param cluster_rows Cluster rows
#' @param cluster_cols Cluster columns
#' @param title Plot title
#' @return heatmaply object
create_heatmap <- function(matrix_data, row_annotation = NULL, col_annotation = NULL,
                           scale_rows = TRUE, cluster_rows = TRUE, cluster_cols = TRUE,
                           title = "Heatmap") {

    # Ensure matrix format
    if (!is.matrix(matrix_data)) {
        matrix_data <- as.matrix(matrix_data)
    }

    # Scale if requested
    if (scale_rows) {
        matrix_data <- t(scale(t(matrix_data)))
    }

    # Create heatmap
    heatmaply(
        matrix_data,
        colors = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
        Rowv = cluster_rows,
        Colv = cluster_cols,
        row_side_colors = row_annotation,
        col_side_colors = col_annotation,
        main = title,
        fontsize_row = 8,
        fontsize_col = 8,
        showticklabels = c(TRUE, TRUE)
    )
}

# =============================================================================
# Pathway Enrichment Plots
# =============================================================================

#' Create pathway enrichment dot plot
#' @param pathway_data Data frame with pathway enrichment results
#' @param pathway_col Column with pathway names
#' @param nes_col Column with NES or enrichment score (NULL for ORA)
#' @param pval_col Column with p-value
#' @param size_col Column for dot size (e.g., gene set size)
#' @param n_pathways Number of pathways to show
#' @param title Plot title
#' @return plotly object
create_pathway_dotplot <- function(pathway_data, pathway_col = "pathway",
                                   nes_col = "NES", pval_col = "padj",
                                   size_col = "size", n_pathways = 20,
                                   title = "Pathway Enrichment") {

    # Select top pathways
    pathway_data <- pathway_data[order(pathway_data[[pval_col]]), ]
    pathway_data <- head(pathway_data, n_pathways)

    # Create plot
    pathway_data$neg_log10_pval <- -log10(pathway_data[[pval_col]] + 1e-300)

    # Truncate pathway names if too long
    pathway_data$pathway_short <- substr(pathway_data[[pathway_col]], 1, 50)

    if (!is.null(nes_col) && nes_col %in% colnames(pathway_data)) {
        # GSEA-style with NES
        p <- ggplot(pathway_data, aes_string(x = nes_col, y = "reorder(pathway_short, neg_log10_pval)")) +
            geom_point(aes_string(size = size_col, color = "neg_log10_pval")) +
            scale_color_gradient(low = "#f39c12", high = "#e74c3c",
                                 name = "-log10(p-adj)") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
            labs(title = title, x = "Normalized Enrichment Score", y = NULL,
                 size = "Gene Set Size") +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 9))
    } else {
        # ORA-style without NES
        p <- ggplot(pathway_data, aes(x = neg_log10_pval, y = reorder(pathway_short, neg_log10_pval))) +
            geom_point(aes_string(size = size_col), color = "#3498db", alpha = 0.7) +
            labs(title = title, x = "-log10(adjusted p-value)", y = NULL,
                 size = "Gene Set Size") +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 9))
    }

    ggplotly(p)
}

# =============================================================================
# Expression/Abundance Plots
# =============================================================================

#' Create expression boxplot for a single gene/feature
#' @param expr_data Data frame with expression values (long format)
#' @param feature_id Feature to plot
#' @param group_var Variable for grouping (x-axis)
#' @param expr_col Column with expression values
#' @param title Plot title
#' @return plotly object
create_expression_boxplot <- function(expr_data, feature_id, group_var = "condition",
                                      expr_col = "expression", title = NULL) {

    if (is.null(title)) {
        title <- paste("Expression of", feature_id)
    }

    p <- ggplot(expr_data, aes_string(x = group_var, y = expr_col, fill = group_var)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
        labs(title = title, x = NULL, y = "Expression") +
        theme_minimal() +
        theme(legend.position = "none")

    ggplotly(p)
}

# =============================================================================
# Correlation Plots
# =============================================================================
#' Create correlation scatter plot
#' @param data Data frame with x and y values
#' @param x_col Column for x-axis
#' @param y_col Column for y-axis
#' @param color_col Optional column for coloring
#' @param label_col Column for hover labels
#' @param title Plot title
#' @return plotly object
create_correlation_scatter <- function(data, x_col, y_col, color_col = NULL,
                                       label_col = NULL, title = "Correlation") {

    # Calculate correlation
    r <- cor(data[[x_col]], data[[y_col]], use = "pairwise.complete.obs")

    # Create hover text
    if (!is.null(label_col) && label_col %in% colnames(data)) {
        data$hover_text <- data[[label_col]]
    } else {
        data$hover_text <- rownames(data)
    }

    p <- ggplot(data, aes_string(x = x_col, y = y_col, text = "hover_text"))

    if (!is.null(color_col) && color_col %in% colnames(data)) {
        p <- p + geom_point(aes_string(color = color_col), alpha = 0.6)
    } else {
        p <- p + geom_point(color = "#3498db", alpha = 0.6)
    }

    p <- p +
        geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        labs(title = paste0(title, " (r = ", round(r, 3), ")"),
             x = x_col, y = y_col) +
        theme_minimal()

    ggplotly(p, tooltip = "text")
}

# =============================================================================
# Network Plots
# =============================================================================

#' Create interactive network visualization
#' @param edges Data frame with from, to, (optional) weight columns
#' @param nodes Data frame with id and optional attributes
#' @param color_var Variable for node coloring
#' @param size_var Variable for node sizing
#' @param title Plot title
#' @return visNetwork object
create_network_plot <- function(edges, nodes = NULL, color_var = NULL,
                                size_var = NULL, title = "Network") {

    # Ensure proper column names for visNetwork
    if (!"from" %in% colnames(edges)) {
        colnames(edges)[1] <- "from"
    }
    if (!"to" %in% colnames(edges)) {
        colnames(edges)[2] <- "to"
    }

    # Create nodes if not provided
    if (is.null(nodes)) {
        unique_nodes <- unique(c(edges$from, edges$to))
        nodes <- data.frame(id = unique_nodes, label = unique_nodes,
                            stringsAsFactors = FALSE)
    } else {
        if (!"id" %in% colnames(nodes)) {
            colnames(nodes)[1] <- "id"
        }
        if (!"label" %in% colnames(nodes)) {
            nodes$label <- nodes$id
        }
    }

    # Apply coloring
    if (!is.null(color_var) && color_var %in% colnames(nodes)) {
        # Map color variable to actual colors
        color_vals <- nodes[[color_var]]
        if (is.numeric(color_vals)) {
            # Gradient for numeric
            nodes$color <- colorRampPalette(c("#3498db", "#e74c3c"))(100)[
                as.integer(cut(color_vals, 100))
            ]
        } else {
            # Categorical colors
            unique_vals <- unique(color_vals)
            color_map <- setNames(
                categorical_colors[1:length(unique_vals)],
                unique_vals
            )
            nodes$color <- color_map[color_vals]
        }
    } else {
        nodes$color <- "#3498db"
    }

    # Apply sizing
    if (!is.null(size_var) && size_var %in% colnames(nodes)) {
        size_vals <- nodes[[size_var]]
        nodes$size <- scales::rescale(size_vals, to = c(10, 40))
    } else {
        nodes$size <- 20
    }

    # Create network
    visNetwork(nodes, edges, main = title) %>%
        visOptions(
            highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1),
            selectedBy = list(variable = color_var, multiple = FALSE),
            nodesIdSelection = TRUE
        ) %>%
        visPhysics(stabilization = list(iterations = 100)) %>%
        visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
        visLayout(randomSeed = 42)
}
