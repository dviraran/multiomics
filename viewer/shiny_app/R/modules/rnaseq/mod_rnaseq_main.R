# =============================================================================
# RNA-seq Main Module (LAZY LOADING)
# =============================================================================
# Orchestrates all RNA-seq visualizations in a tabbed interface
# Data is loaded only when this tab is first accessed

rnaseq_main_ui <- function(id) {
    ns <- NS(id)

    navset_card_tab(
        id = ns("rnaseq_tabs"),

        # QC Tab
        nav_panel(
            title = "QC",
            icon = icon("chart-bar"),
            rnaseq_qc_ui(ns("qc"))
        ),

        # PCA Tab
        nav_panel(
            title = "PCA",
            icon = icon("circle-dot"),
            rnaseq_pca_ui(ns("pca"))
        ),

        # Differential Expression Tab
        nav_panel(
            title = "Differential Expression",
            icon = icon("volcano"),
            rnaseq_de_ui(ns("de"))
        ),

        # Pathway Enrichment Tab
        nav_panel(
            title = "Pathway Enrichment",
            icon = icon("diagram-project"),
            rnaseq_pathway_ui(ns("pathway"))
        ),

        # Gene Browser Tab
        nav_panel(
            title = "Gene Browser",
            icon = icon("search"),
            rnaseq_browser_ui(ns("browser"))
        )
    )
}

rnaseq_main_server <- function(id, data_dir) {
    moduleServer(id, function(input, output, session) {

        # LAZY LOADING: Only load data when this module is accessed
        data <- reactiveVal(NULL)
        is_loaded <- reactiveVal(FALSE)

        # Load data on first access (triggered by any sub-module)
        load_data <- reactive({
            if (!is_loaded() && !is.null(data_dir)) {
                message("Loading RNA-seq data from: ", data_dir)
                loaded <- load_rnaseq_data(data_dir)
                data(loaded)
                is_loaded(TRUE)
            }
            data()
        })

        # Call sub-modules with lazy-loaded data
        rnaseq_qc_server("qc", load_data)
        rnaseq_pca_server("pca", load_data)
        rnaseq_de_server("de", load_data)
        rnaseq_pathway_server("pathway", load_data)
        rnaseq_browser_server("browser", load_data)
    })
}

# =============================================================================
# QC Sub-module
# =============================================================================

rnaseq_qc_ui <- function(id) {
    ns <- NS(id)

    tagList(
        # Library sizes - full width
        card(
            full_screen = TRUE,
            card_header("Library Sizes"),
            card_body(plotlyOutput(ns("library_sizes"), height = "400px"))
        ),

        # Detected genes - full width
        card(
            full_screen = TRUE,
            card_header("Detected Genes"),
            card_body(plotlyOutput(ns("detected_genes"), height = "400px"))
        ),

        # Sample correlation heatmap - full width
        card(
            full_screen = TRUE,
            card_header("Sample Correlation"),
            card_body(plotlyOutput(ns("correlation_heatmap"), height = "500px"))
        ),

        # QC metrics table - full width
        card(
            card_header("QC Metrics"),
            card_body(DTOutput(ns("qc_table")))
        )
    )
}

rnaseq_qc_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        # Library sizes bar plot
        output$library_sizes <- renderPlotly({
            data <- data_reactive()
            req(data, data$qc_metrics)

            qc <- data$qc_metrics
            count_col <- intersect(c("total_counts", "total_normalized", "library_size"),
                                   colnames(qc))[1]

            if (is.na(count_col)) return(NULL)

            qc <- qc[order(qc[[count_col]], decreasing = TRUE), ]

            p <- ggplot(qc, aes_string(x = "reorder(sample, -get(count_col))",
                                       y = count_col)) +
                geom_col(fill = "#3498db", alpha = 0.8) +
                labs(x = "Sample", y = "Total Counts") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

            ggplotly(p)
        })

        # Detected genes bar plot
        output$detected_genes <- renderPlotly({
            data <- data_reactive()
            req(data, data$qc_metrics)

            qc <- data$qc_metrics
            gene_col <- intersect(c("detected_genes", "n_genes", "genes_detected"),
                                  colnames(qc))[1]

            if (is.na(gene_col)) return(NULL)

            qc <- qc[order(qc[[gene_col]], decreasing = TRUE), ]

            p <- ggplot(qc, aes_string(x = "reorder(sample, -get(gene_col))",
                                       y = gene_col)) +
                geom_col(fill = "#2ecc71", alpha = 0.8) +
                labs(x = "Sample", y = "Detected Genes") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

            ggplotly(p)
        })

        # Correlation heatmap
        output$correlation_heatmap <- renderPlotly({
            data <- data_reactive()
            req(data, data$sample_correlation)

            cor_mat <- data$sample_correlation

            # Convert to matrix if needed
            if (!is.matrix(cor_mat)) {
                if (colnames(cor_mat)[1] %in% c("", "X", "sample")) {
                    rownames(cor_mat) <- cor_mat[[1]]
                    cor_mat <- cor_mat[, -1]
                }
                cor_mat <- as.matrix(cor_mat)
            }

            heatmaply(cor_mat,
                      colors = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
                      main = "",
                      fontsize_row = 8,
                      fontsize_col = 8)
        })

        # QC metrics table
        output$qc_table <- renderDT({
            data <- data_reactive()
            req(data, data$qc_metrics)

            datatable(
                data$qc_metrics,
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE
            ) %>%
                formatRound(columns = which(sapply(data$qc_metrics, is.numeric)),
                            digits = 2)
        })
    })
}

# =============================================================================
# PCA Sub-module
# =============================================================================

rnaseq_pca_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectInput(ns("pc_x"), "X-axis PC:", choices = paste0("PC", 1:10), selected = "PC1"),
            selectInput(ns("pc_y"), "Y-axis PC:", choices = paste0("PC", 1:10), selected = "PC2"),
            selectInput(ns("color_by"), "Color by:", choices = NULL),
            checkboxInput(ns("show_labels"), "Show sample labels", FALSE),
            checkboxInput(ns("show_ellipse"), "Show confidence ellipse", TRUE)
        ),

        tagList(
            # PCA scatter plot - full width
            card(
                full_screen = TRUE,
                card_header("PCA Plot"),
                card_body(plotlyOutput(ns("pca_plot"), height = "550px"))
            ),

            # Scree plot - full width
            card(
                full_screen = TRUE,
                card_header("Variance Explained"),
                card_body(plotOutput(ns("scree_plot"), height = "350px"))
            )
        )
    )
}

rnaseq_pca_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Update color choices when data loads
        observe({
            data <- data_reactive()
            req(data, data$pca_results)
            pca <- data$pca_results

            # Get non-PC columns
            pc_cols <- grep("^PC[0-9]", colnames(pca), value = TRUE)
            other_cols <- setdiff(colnames(pca), c(pc_cols, "sample"))

            # Add metadata columns if available
            if (!is.null(data$metadata)) {
                meta_cols <- get_metadata_columns(data$metadata)
                other_cols <- unique(c(other_cols, meta_cols))
            }

            updateSelectInput(session, "color_by", choices = other_cols,
                              selected = other_cols[1])
        })

        # Extract variance explained
        var_explained <- reactive({
            data <- data_reactive()
            req(data, data$pca_results)
            extract_variance_explained(data$pca_results)
        })

        # PCA plot
        output$pca_plot <- renderPlotly({
            data <- data_reactive()
            req(data, data$pca_results, input$pc_x, input$pc_y)

            pca_df <- data$pca_results

            # Merge with metadata if needed
            if (!is.null(input$color_by) &&
                !input$color_by %in% colnames(pca_df) &&
                !is.null(data$metadata)) {

                sample_col <- intersect(c("sample", "sample_id"), colnames(pca_df))[1]
                meta_sample_col <- intersect(c("sample", "sample_id"), colnames(data$metadata))[1]

                if (!is.na(sample_col) && !is.na(meta_sample_col)) {
                    pca_df <- merge(pca_df, data$metadata,
                                    by.x = sample_col, by.y = meta_sample_col,
                                    all.x = TRUE)
                }
            }

            create_pca_plot(
                pca_df,
                pc_x = input$pc_x,
                pc_y = input$pc_y,
                color_var = input$color_by,
                var_explained = var_explained(),
                show_labels = input$show_labels,
                show_ellipse = input$show_ellipse
            )
        })

        # Scree plot
        output$scree_plot <- renderPlot({
            req(var_explained())
            create_scree_plot(var_explained())
        })
    })
}

# =============================================================================
# Differential Expression Sub-module
# =============================================================================

rnaseq_de_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectInput(ns("contrast"), "Select Contrast:",
                        choices = NULL),
            hr(),
            sliderInput(ns("log2fc_thresh"), "log2FC threshold:",
                        min = 0, max = 3, value = 1, step = 0.25),
            sliderInput(ns("pval_thresh"), "P-adj threshold:",
                        min = 0.001, max = 0.1, value = 0.05, step = 0.005),
            hr(),
            uiOutput(ns("de_summary"))
        ),

        tagList(
            # Volcano plot - full width
            card(
                full_screen = TRUE,
                card_header("Volcano Plot"),
                card_body(plotlyOutput(ns("volcano_plot"), height = "550px"))
            ),

            # MA plot - full width
            card(
                full_screen = TRUE,
                card_header("MA Plot"),
                card_body(plotlyOutput(ns("ma_plot"), height = "450px"))
            ),

            # DE results table - full width
            card(
                card_header("Differential Expression Results"),
                card_body(DTOutput(ns("de_table")))
            )
        )
    )
}

rnaseq_de_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Update contrast choices
        observe({
            data <- data_reactive()
            req(data, data$de_results)
            contrasts <- get_contrasts(data$de_results)
            updateSelectInput(session, "contrast", choices = contrasts,
                              selected = contrasts[1])
        })

        # Get current DE results
        current_de <- reactive({
            data <- data_reactive()
            req(input$contrast, data, data$de_results)
            get_de_for_contrast(data$de_results, input$contrast)
        })

        # DE summary
        output$de_summary <- renderUI({
            req(current_de())

            de <- current_de()

            # Find column names
            pval_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(de))[1]
            fc_col <- intersect(c("log2FoldChange", "log2FC", "logFC"), colnames(de))[1]

            if (is.na(pval_col) || is.na(fc_col)) return(NULL)

            summary <- summarize_de_results(de, pval_col, fc_col,
                                            input$pval_thresh, input$log2fc_thresh)

            tagList(
                tags$div(class = "text-center",
                         tags$p(tags$strong("Total genes: "), format_number(summary$n_total)),
                         tags$p(tags$strong("Significant: "), format_number(summary$n_significant),
                                paste0(" (", summary$pct_significant, "%)")),
                         tags$p(style = "color: #e74c3c;",
                                tags$strong("Up: "), format_number(summary$n_up)),
                         tags$p(style = "color: #3498db;",
                                tags$strong("Down: "), format_number(summary$n_down))
                )
            )
        })

        # Volcano plot
        output$volcano_plot <- renderPlotly({
            req(current_de())

            de <- current_de()

            # Find column names
            pval_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(de))[1]
            fc_col <- intersect(c("log2FoldChange", "log2FC", "logFC"), colnames(de))[1]
            label_col <- intersect(c("symbol", "gene_symbol", "gene_name", "gene_id"),
                                   colnames(de))[1]

            if (is.na(pval_col) || is.na(fc_col)) return(NULL)

            create_volcano_plot(
                de,
                log2fc_col = fc_col,
                pval_col = pval_col,
                label_col = label_col,
                log2fc_thresh = input$log2fc_thresh,
                pval_thresh = input$pval_thresh,
                title = paste("Volcano:", input$contrast)
            )
        })

        # MA plot
        output$ma_plot <- renderPlotly({
            req(current_de())

            de <- current_de()

            # Find column names
            fc_col <- intersect(c("log2FoldChange", "log2FC", "logFC"), colnames(de))[1]
            mean_col <- intersect(c("baseMean", "AveExpr", "avg_intensity"), colnames(de))[1]
            pval_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(de))[1]

            if (is.na(fc_col) || is.na(mean_col)) return(NULL)

            de$significance <- "ns"
            de$significance[de[[pval_col]] < input$pval_thresh &
                              de[[fc_col]] > input$log2fc_thresh] <- "up"
            de$significance[de[[pval_col]] < input$pval_thresh &
                              de[[fc_col]] < -input$log2fc_thresh] <- "down"

            p <- ggplot(de, aes_string(x = paste0("log10(", mean_col, " + 1)"),
                                       y = fc_col, color = "significance")) +
                geom_point(alpha = 0.5, size = 1) +
                scale_color_manual(values = c("up" = "#e74c3c", "down" = "#3498db", "ns" = "#95a5a6")) +
                geom_hline(yintercept = 0, linetype = "dashed") +
                labs(x = "log10(Mean Expression)", y = "log2 Fold Change") +
                theme_minimal()

            ggplotly(p)
        })

        # DE table
        output$de_table <- renderDT({
            req(current_de())

            de <- current_de()

            # Select and order columns
            pval_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(de))[1]
            fc_col <- intersect(c("log2FoldChange", "log2FC", "logFC"), colnames(de))[1]

            # Round numeric columns
            numeric_cols <- which(sapply(de, is.numeric))

            datatable(
                de,
                options = list(
                    pageLength = 15,
                    scrollX = TRUE,
                    order = list(list(which(colnames(de) == pval_col) - 1, "asc"))
                ),
                rownames = FALSE,
                filter = "top"
            ) %>%
                formatRound(columns = numeric_cols, digits = 3) %>%
                formatSignif(columns = pval_col, digits = 3)
        })
    })
}

# =============================================================================
# Pathway Enrichment Sub-module
# =============================================================================

rnaseq_pathway_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectInput(ns("pathway_source"), "Pathway Database:",
                        choices = NULL),
            selectInput(ns("pathway_contrast"), "Contrast:",
                        choices = NULL),
            sliderInput(ns("n_pathways"), "Top pathways to show:",
                        min = 10, max = 50, value = 20)
        ),

        tagList(
            # Pathway dot plot - full width
            card(
                full_screen = TRUE,
                card_header("Pathway Enrichment"),
                card_body(plotlyOutput(ns("pathway_plot"), height = "600px"))
            ),

            # Pathway table - full width
            card(
                card_header("Enrichment Results"),
                card_body(DTOutput(ns("pathway_table")))
            )
        )
    )
}

rnaseq_pathway_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Update choices
        observe({
            data <- data_reactive()
            req(data, data$pathway_results)

            # Get available pathway sources
            sources <- names(data$pathway_results)
            updateSelectInput(session, "pathway_source", choices = sources,
                              selected = sources[1])

            # Get contrasts
            contrasts <- get_contrasts(data$de_results)
            updateSelectInput(session, "pathway_contrast", choices = contrasts,
                              selected = contrasts[1])
        })

        # Current pathway data
        current_pathways <- reactive({
            data <- data_reactive()
            req(input$pathway_source, data, data$pathway_results)

            pathways <- data$pathway_results[[input$pathway_source]]

            # Filter by contrast if column exists
            if (!is.null(input$pathway_contrast) && "contrast" %in% colnames(pathways)) {
                pathways <- pathways[pathways$contrast == input$pathway_contrast, ]
            }

            pathways
        })

        # Pathway dot plot
        output$pathway_plot <- renderPlotly({
            req(current_pathways())

            pathways <- current_pathways()

            # Find columns
            pathway_col <- intersect(c("pathway", "term", "Description", "pathway_name"),
                                     colnames(pathways))[1]
            pval_col <- intersect(c("padj", "p.adjust", "FDR"), colnames(pathways))[1]
            nes_col <- intersect(c("NES", "enrichmentScore"), colnames(pathways))
            nes_col <- if(length(nes_col) > 0) nes_col[1] else NULL
            size_col <- intersect(c("size", "setSize", "Count"), colnames(pathways))[1]

            if (is.na(pathway_col) || is.na(pval_col)) return(NULL)

            create_pathway_dotplot(
                pathways,
                pathway_col = pathway_col,
                nes_col = nes_col,
                pval_col = pval_col,
                size_col = size_col,
                n_pathways = input$n_pathways,
                title = paste("Pathway Enrichment:", input$pathway_source)
            )
        })

        # Pathway table
        output$pathway_table <- renderDT({
            req(current_pathways())

            pathways <- current_pathways()

            # Select key columns
            display_cols <- intersect(
                c("pathway", "term", "Description", "NES", "pval", "padj", "size", "leadingEdge"),
                colnames(pathways)
            )

            datatable(
                pathways[, display_cols, drop = FALSE],
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE,
                filter = "top"
            )
        })
    })
}

# =============================================================================
# Gene Browser Sub-module
# =============================================================================

rnaseq_browser_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectizeInput(ns("gene_select"), "Search Gene:",
                           choices = NULL,
                           options = list(
                               placeholder = "Start typing gene name...",
                               maxOptions = 100
                           )),
            selectInput(ns("group_var"), "Group by:",
                        choices = NULL)
        ),

        card(
            full_screen = TRUE,
            card_header(textOutput(ns("gene_title"))),
            card_body(plotlyOutput(ns("expression_plot"), height = "400px"))
        )
    )
}

rnaseq_browser_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Update gene choices
        observe({
            data <- data_reactive()
            req(data, data$normalized_counts)

            # Get gene IDs from row names or first column
            counts <- data$normalized_counts
            if ("gene_id" %in% colnames(counts)) {
                genes <- counts$gene_id
            } else if (is.character(rownames(counts))) {
                genes <- rownames(counts)
            } else {
                genes <- counts[[1]]
            }

            # Add gene symbols if available
            if (!is.null(data$gene_annotation)) {
                annot <- data$gene_annotation
                symbol_col <- intersect(c("symbol", "gene_symbol", "gene_name"),
                                        colnames(annot))[1]
                if (!is.na(symbol_col)) {
                    symbols <- annot[[symbol_col]]
                    genes <- unique(c(symbols[!is.na(symbols)], genes))
                }
            }

            updateSelectizeInput(session, "gene_select", choices = genes, server = TRUE)
        })

        # Update group choices
        observe({
            data <- data_reactive()
            if (!is.null(data) && !is.null(data$metadata)) {
                cols <- get_metadata_columns(data$metadata)
                updateSelectInput(session, "group_var", choices = cols,
                                  selected = cols[1])
            }
        })

        # Gene title
        output$gene_title <- renderText({
            req(input$gene_select)
            paste("Expression of", input$gene_select)
        })

        # Expression plot
        output$expression_plot <- renderPlotly({
            data <- data_reactive()
            req(input$gene_select, data, data$normalized_counts)

            counts <- data$normalized_counts
            gene <- input$gene_select

            # Prepare expression data
            expr_df <- prepare_expression_data(
                counts,
                gene,
                data$metadata,
                sample_col = "sample"
            )

            if (is.null(expr_df) || nrow(expr_df) == 0) {
                return(NULL)
            }

            # Create boxplot
            group_var <- if (!is.null(input$group_var) && input$group_var %in% colnames(expr_df)) {
                input$group_var
            } else {
                "sample"
            }

            create_expression_boxplot(
                expr_df,
                gene,
                group_var = group_var,
                expr_col = "expression",
                title = paste("Expression of", gene)
            )
        })
    })
}
