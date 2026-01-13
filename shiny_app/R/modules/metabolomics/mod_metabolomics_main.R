# =============================================================================
# Metabolomics Main Module
# =============================================================================

metabolomics_main_ui <- function(id) {
    ns <- NS(id)

    navset_card_tab(
        id = ns("metab_tabs"),

        # QC Tab
        nav_panel(
            title = "QC",
            icon = icon("chart-bar"),
            metabolomics_qc_ui(ns("qc"))
        ),

        # PCA Tab
        nav_panel(
            title = "PCA",
            icon = icon("chart-scatter"),
            metabolomics_pca_ui(ns("pca"))
        ),

        # Differential Analysis Tab
        nav_panel(
            title = "Differential Analysis",
            icon = icon("volcano"),
            metabolomics_da_ui(ns("da"))
        ),

        # Metabolite Browser Tab
        nav_panel(
            title = "Metabolite Browser",
            icon = icon("search"),
            metabolomics_browser_ui(ns("browser"))
        )
    )
}

metabolomics_main_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {
        metabolomics_qc_server("qc", data)
        metabolomics_pca_server("pca", data)
        metabolomics_da_server("da", data)
        metabolomics_browser_server("browser", data)
    })
}

# =============================================================================
# QC Sub-module
# =============================================================================

metabolomics_qc_ui <- function(id) {
    ns <- NS(id)

    layout_columns(
        col_widths = c(6, 6),

        # Missing values
        card(
            full_screen = TRUE,
            card_header("Missing Value Pattern"),
            card_body(plotlyOutput(ns("missing_plot"), height = "400px"))
        ),

        # Sample correlation
        card(
            full_screen = TRUE,
            card_header("Sample Correlation"),
            card_body(plotlyOutput(ns("correlation_heatmap"), height = "400px"))
        ),

        # Summary statistics
        card(
            card_header("Data Summary"),
            card_body(uiOutput(ns("data_summary")))
        )
    )
}

metabolomics_qc_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {

        output$missing_plot <- renderPlotly({
            req(data$normalized_matrix)

            mat <- data$normalized_matrix
            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "metabolite_id")) {
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            # Calculate missing per sample
            missing_per_sample <- colSums(is.na(mat))
            df <- data.frame(
                sample = names(missing_per_sample),
                missing = missing_per_sample
            )

            p <- ggplot(df, aes(x = reorder(sample, missing), y = missing)) +
                geom_col(fill = "#e74c3c", alpha = 0.7) +
                coord_flip() +
                labs(x = "Sample", y = "Missing Features") +
                theme_minimal()

            ggplotly(p)
        })

        output$correlation_heatmap <- renderPlotly({
            req(data$normalized_matrix)

            mat <- data$normalized_matrix
            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "metabolite_id")) {
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            cor_mat <- cor(mat, use = "pairwise.complete.obs")

            heatmaply(cor_mat,
                      colors = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
                      fontsize_row = 8,
                      fontsize_col = 8)
        })

        output$data_summary <- renderUI({
            req(data$normalized_matrix)

            mat <- data$normalized_matrix
            if (!is.matrix(mat)) {
                n_features <- nrow(mat) - ifelse(colnames(mat)[1] %in% c("", "X", "feature_id"), 0, 0)
                n_samples <- ncol(mat) - 1
                mat_numeric <- mat[, -1]
            } else {
                n_features <- nrow(mat)
                n_samples <- ncol(mat)
                mat_numeric <- mat
            }

            pct_missing <- round(100 * sum(is.na(mat_numeric)) / length(mat_numeric), 1)

            tagList(
                tags$p(tags$strong("Features: "), format_number(n_features)),
                tags$p(tags$strong("Samples: "), format_number(n_samples)),
                tags$p(tags$strong("Missing: "), paste0(pct_missing, "%"))
            )
        })
    })
}

# =============================================================================
# PCA Sub-module
# =============================================================================

metabolomics_pca_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectInput(ns("pc_x"), "X-axis PC:", choices = paste0("PC", 1:10), selected = "PC1"),
            selectInput(ns("pc_y"), "Y-axis PC:", choices = paste0("PC", 1:10), selected = "PC2"),
            selectInput(ns("color_by"), "Color by:", choices = NULL),
            checkboxInput(ns("show_labels"), "Show sample labels", FALSE)
        ),

        layout_columns(
            col_widths = c(8, 4),

            card(
                full_screen = TRUE,
                card_header("PCA Plot"),
                card_body(plotlyOutput(ns("pca_plot"), height = "500px"))
            ),

            card(
                card_header("Top Loadings"),
                card_body(DTOutput(ns("loadings_table")))
            )
        )
    )
}

metabolomics_pca_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Use provided PCA or compute
        pca_data <- reactive({
            if (!is.null(data$pca_scores)) {
                pca_df <- data$pca_scores
                if (!"sample" %in% colnames(pca_df)) {
                    pca_df$sample <- rownames(pca_df)
                }
                return(pca_df)
            }

            req(data$normalized_matrix)
            mat <- data$normalized_matrix

            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "metabolite_id")) {
                    rownames(mat) <- mat[[1]]
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            # Handle missing values
            row_na <- rowMeans(is.na(mat))
            mat <- mat[row_na < 0.5, ]
            for (i in 1:ncol(mat)) {
                mat[is.na(mat[, i]), i] <- mean(mat[, i], na.rm = TRUE)
            }

            pca <- prcomp(t(mat), scale. = TRUE)
            pca_df <- as.data.frame(pca$x[, 1:min(10, ncol(pca$x))])
            pca_df$sample <- rownames(pca_df)
            attr(pca_df, "var_explained") <- 100 * pca$sdev^2 / sum(pca$sdev^2)
            attr(pca_df, "loadings") <- pca$rotation

            pca_df
        })

        observe({
            req(pca_data())
            cols <- setdiff(colnames(pca_data()), c(paste0("PC", 1:20), "sample"))
            if (!is.null(data$metadata)) {
                meta_cols <- get_metadata_columns(data$metadata)
                cols <- unique(c(cols, meta_cols))
            }
            updateSelectInput(session, "color_by", choices = cols,
                              selected = if(length(cols) > 0) cols[1] else NULL)
        })

        output$pca_plot <- renderPlotly({
            req(pca_data(), input$pc_x, input$pc_y)

            pca_df <- pca_data()
            var_exp <- attr(pca_df, "var_explained")

            if (!is.null(input$color_by) && !input$color_by %in% colnames(pca_df) &&
                !is.null(data$metadata)) {
                meta_sample_col <- intersect(c("sample", "sample_id"), colnames(data$metadata))[1]
                if (!is.na(meta_sample_col)) {
                    pca_df <- merge(pca_df, data$metadata,
                                    by.x = "sample", by.y = meta_sample_col, all.x = TRUE)
                }
            }

            create_pca_plot(pca_df, input$pc_x, input$pc_y, input$color_by,
                            var_exp, input$show_labels, TRUE)
        })

        output$loadings_table <- renderDT({
            req(pca_data())

            if (!is.null(data$pca_loadings)) {
                loadings <- data$pca_loadings
            } else {
                loadings <- attr(pca_data(), "loadings")
                if (is.null(loadings)) return(NULL)
                loadings <- as.data.frame(loadings[, 1:min(5, ncol(loadings))])
                loadings$feature <- rownames(loadings)
            }

            # Get top loadings for PC1
            if ("PC1" %in% colnames(loadings)) {
                loadings <- loadings[order(abs(loadings$PC1), decreasing = TRUE), ]
            }

            datatable(
                head(loadings, 20),
                options = list(pageLength = 10, dom = "t"),
                rownames = FALSE
            ) %>% formatRound(columns = which(sapply(loadings, is.numeric)), digits = 3)
        })
    })
}

# =============================================================================
# Differential Analysis Sub-module
# =============================================================================

metabolomics_da_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectInput(ns("contrast"), "Select Contrast:", choices = NULL),
            hr(),
            sliderInput(ns("log2fc_thresh"), "log2FC threshold:",
                        min = 0, max = 3, value = 1, step = 0.25),
            sliderInput(ns("pval_thresh"), "P-adj threshold:",
                        min = 0.001, max = 0.1, value = 0.05, step = 0.005),
            hr(),
            uiOutput(ns("da_summary"))
        ),

        layout_columns(
            col_widths = c(6, 6),

            card(
                full_screen = TRUE,
                card_header("Volcano Plot"),
                card_body(plotlyOutput(ns("volcano_plot"), height = "450px"))
            ),

            card(
                card_header("Differential Results"),
                card_body(DTOutput(ns("da_table")))
            )
        )
    )
}

metabolomics_da_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        observe({
            req(data$de_results)
            contrasts <- get_contrasts(data$de_results)
            updateSelectInput(session, "contrast", choices = contrasts, selected = contrasts[1])
        })

        current_da <- reactive({
            req(input$contrast, data$de_results)
            get_de_for_contrast(data$de_results, input$contrast)
        })

        output$da_summary <- renderUI({
            req(current_da())

            da <- current_da()
            pval_col <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(da))[1]
            fc_col <- intersect(c("log2FC", "log2FoldChange", "logFC"), colnames(da))[1]

            if (is.na(pval_col) || is.na(fc_col)) return(NULL)

            summary <- summarize_de_results(da, pval_col, fc_col,
                                            input$pval_thresh, input$log2fc_thresh)

            tagList(
                tags$div(class = "text-center",
                         tags$p(tags$strong("Total: "), format_number(summary$n_total)),
                         tags$p(tags$strong("Significant: "), format_number(summary$n_significant)),
                         tags$p(style = "color: #e74c3c;", tags$strong("Up: "), format_number(summary$n_up)),
                         tags$p(style = "color: #3498db;", tags$strong("Down: "), format_number(summary$n_down))
                )
            )
        })

        output$volcano_plot <- renderPlotly({
            req(current_da())

            da <- current_da()
            pval_col <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(da))[1]
            fc_col <- intersect(c("log2FC", "log2FoldChange", "logFC"), colnames(da))[1]
            label_col <- intersect(c("metabolite_name", "feature_id", "compound"),
                                   colnames(da))[1]

            if (is.na(pval_col) || is.na(fc_col)) return(NULL)

            create_volcano_plot(da, fc_col, pval_col, label_col,
                                input$log2fc_thresh, input$pval_thresh,
                                paste("Volcano:", input$contrast))
        })

        output$da_table <- renderDT({
            req(current_da())
            datatable(current_da(),
                      options = list(pageLength = 15, scrollX = TRUE),
                      rownames = FALSE, filter = "top")
        })
    })
}

# =============================================================================
# Metabolite Browser Sub-module
# =============================================================================

metabolomics_browser_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectizeInput(ns("metabolite_select"), "Search Metabolite:",
                           choices = NULL,
                           options = list(placeholder = "Start typing...")),
            selectInput(ns("group_var"), "Group by:", choices = NULL)
        ),

        card(
            full_screen = TRUE,
            card_header(textOutput(ns("metabolite_title"))),
            card_body(plotlyOutput(ns("abundance_plot"), height = "400px"))
        )
    )
}

metabolomics_browser_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        observe({
            req(data$normalized_matrix)

            mat <- data$normalized_matrix
            features <- if ("feature_id" %in% colnames(mat)) mat$feature_id
            else if ("metabolite_id" %in% colnames(mat)) mat$metabolite_id
            else if (is.character(rownames(mat))) rownames(mat)
            else mat[[1]]

            updateSelectizeInput(session, "metabolite_select", choices = features, server = TRUE)
        })

        observe({
            if (!is.null(data$metadata)) {
                cols <- get_metadata_columns(data$metadata)
                updateSelectInput(session, "group_var", choices = cols,
                                  selected = if(length(cols) > 0) cols[1] else NULL)
            }
        })

        output$metabolite_title <- renderText({
            req(input$metabolite_select)
            paste("Abundance of", input$metabolite_select)
        })

        output$abundance_plot <- renderPlotly({
            req(input$metabolite_select, data$normalized_matrix)

            expr_df <- prepare_expression_data(
                data$normalized_matrix,
                input$metabolite_select,
                data$metadata,
                sample_col = "sample"
            )

            if (is.null(expr_df)) return(NULL)

            group_var <- if (!is.null(input$group_var) && input$group_var %in% colnames(expr_df)) {
                input$group_var
            } else {
                "sample"
            }

            create_expression_boxplot(expr_df, input$metabolite_select, group_var, "expression")
        })
    })
}
