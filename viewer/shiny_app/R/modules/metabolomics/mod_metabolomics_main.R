# =============================================================================
# Metabolomics Main Module (LAZY LOADING)
# =============================================================================
# Data is loaded only when this tab is first accessed

metabolomics_main_ui <- function(id) {
    ns <- NS(id)

    navset_card_tab(
        id = ns("metab_tabs"),

        nav_panel(
            title = "QC",
            icon = icon("chart-bar"),
            metabolomics_qc_ui(ns("qc"))
        ),

        nav_panel(
            title = "PCA",
            icon = icon("circle-dot"),
            metabolomics_pca_ui(ns("pca"))
        ),

        nav_panel(
            title = "Differential Analysis",
            icon = icon("volcano"),
            metabolomics_da_ui(ns("da"))
        ),

        nav_panel(
            title = "Metabolite Browser",
            icon = icon("search"),
            metabolomics_browser_ui(ns("browser"))
        )
    )
}

metabolomics_main_server <- function(id, data_dir) {
    moduleServer(id, function(input, output, session) {

        # LAZY LOADING
        data <- reactiveVal(NULL)
        is_loaded <- reactiveVal(FALSE)

        load_data <- reactive({
            if (!is_loaded() && !is.null(data_dir)) {
                message("Loading Metabolomics data from: ", data_dir)
                loaded <- load_metabolomics_data(data_dir)
                data(loaded)
                is_loaded(TRUE)
            }
            data()
        })

        metabolomics_qc_server("qc", load_data)
        metabolomics_pca_server("pca", load_data)
        metabolomics_da_server("da", load_data)
        metabolomics_browser_server("browser", load_data)
    })
}

# =============================================================================
# QC Sub-module
# =============================================================================

metabolomics_qc_ui <- function(id) {
    ns <- NS(id)

    layout_column_wrap(
        width = 0.5,  # Two columns side by side

        card(
            full_screen = TRUE,
            card_header("Missing Value Pattern"),
            card_body(plotlyOutput(ns("missing_heatmap"), height = "800px"))
        ),

        card(
            card_header("QC Summary"),
            card_body(uiOutput(ns("qc_summary")))
        )
    )
}

metabolomics_qc_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        output$missing_heatmap <- renderPlotly({
            data <- data_reactive()
            req(data, data$normalized_matrix)

            mat <- data$normalized_matrix
            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "metabolite_id")) {
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            missing_mat <- ifelse(is.na(mat), 1, 0)
            if (nrow(missing_mat) > 500) {
                missing_mat <- missing_mat[sample(1:nrow(missing_mat), 500), ]
            }

            heatmaply(missing_mat,
                      colors = c("white", "#e74c3c"),
                      main = "Missing Values",
                      fontsize_row = 6,
                      fontsize_col = 8,
                      showticklabels = c(FALSE, TRUE))
        })

        output$qc_summary <- renderUI({
            data <- data_reactive()
            req(data, data$normalized_matrix)

            mat <- data$normalized_matrix
            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "metabolite_id")) {
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            tagList(
                tags$p(tags$strong("Features: "), format_number(nrow(mat))),
                tags$p(tags$strong("Samples: "), format_number(ncol(mat))),
                tags$p(tags$strong("Missing: "), paste0(round(100 * sum(is.na(mat)) / length(mat), 1), "%"))
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
            checkboxInput(ns("show_labels"), "Show labels", FALSE)
        ),

        card(
            full_screen = TRUE,
            card_header("PCA Plot"),
            card_body(plotlyOutput(ns("pca_plot"), height = "800px"))
        )
    )
}

metabolomics_pca_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        pca_data <- reactive({
            data <- data_reactive()
            req(data)

            if (!is.null(data$pca_scores)) {
                return(data$pca_scores)
            }

            req(data$normalized_matrix)
            mat <- data$normalized_matrix
            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "metabolite_id")) {
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            row_na <- rowMeans(is.na(mat))
            mat <- mat[row_na < 0.5, ]

            for (i in 1:ncol(mat)) {
                mat[is.na(mat[, i]), i] <- mean(mat[, i], na.rm = TRUE)
            }

            pca <- prcomp(t(mat), scale. = TRUE)
            pca_df <- as.data.frame(pca$x[, 1:min(10, ncol(pca$x))])
            pca_df$sample <- rownames(pca_df)
            attr(pca_df, "var_explained") <- 100 * pca$sdev^2 / sum(pca$sdev^2)
            pca_df
        })

        observe({
            req(pca_data())
            data <- data_reactive()
            cols <- setdiff(colnames(pca_data()), c(paste0("PC", 1:20), "sample"))
            if (!is.null(data) && !is.null(data$metadata)) {
                cols <- unique(c(cols, get_metadata_columns(data$metadata)))
            }
            updateSelectInput(session, "color_by", choices = cols,
                              selected = if(length(cols) > 0) cols[1] else NULL)
        })

        output$pca_plot <- renderPlotly({
            data <- data_reactive()
            req(pca_data(), input$pc_x, input$pc_y)

            pca_df <- pca_data()

            if (!is.null(input$color_by) &&
                !input$color_by %in% colnames(pca_df) &&
                !is.null(data) && !is.null(data$metadata)) {
                meta_sample_col <- intersect(c("sample", "sample_id"), colnames(data$metadata))[1]
                if (!is.na(meta_sample_col)) {
                    pca_df <- merge(pca_df, data$metadata,
                                    by.x = "sample", by.y = meta_sample_col, all.x = TRUE)
                }
            }

            create_pca_plot(pca_df, input$pc_x, input$pc_y, input$color_by,
                            attr(pca_data(), "var_explained"),
                            input$show_labels, TRUE)
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
            selectInput(ns("contrast"), "Contrast:", choices = NULL),
            sliderInput(ns("log2fc_thresh"), "log2FC threshold:", min = 0, max = 3, value = 1),
            sliderInput(ns("pval_thresh"), "P-adj threshold:", min = 0.001, max = 0.1, value = 0.05),
            uiOutput(ns("da_summary"))
        ),

        layout_column_wrap(
            width = 0.5,  # Two columns side by side

            card(full_screen = TRUE, card_header("Volcano Plot"),
                 card_body(plotlyOutput(ns("volcano_plot"), height = "800px"))),
            card(card_header("Results Table"),
                 card_body(DTOutput(ns("da_table"))))
        )
    )
}

metabolomics_da_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        observe({
            data <- data_reactive()
            req(data, data$de_results)
            contrasts <- get_contrasts(data$de_results)
            updateSelectInput(session, "contrast", choices = contrasts, selected = contrasts[1])
        })

        current_da <- reactive({
            data <- data_reactive()
            req(input$contrast, data, data$de_results)
            get_de_for_contrast(data$de_results, input$contrast)
        })

        output$da_summary <- renderUI({
            req(current_da())
            da <- current_da()
            pval_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(da))[1]
            fc_col <- intersect(c("log2FC", "log2FoldChange", "logFC"), colnames(da))[1]
            if (is.na(pval_col) || is.na(fc_col)) return(NULL)

            summary <- summarize_de_results(da, pval_col, fc_col, input$pval_thresh, input$log2fc_thresh)
            tagList(
                tags$p(tags$strong("Total: "), format_number(summary$n_total)),
                tags$p(tags$strong("Significant: "), format_number(summary$n_significant)),
                tags$p(style = "color: #e74c3c;", tags$strong("Up: "), format_number(summary$n_up)),
                tags$p(style = "color: #3498db;", tags$strong("Down: "), format_number(summary$n_down))
            )
        })

        output$volcano_plot <- renderPlotly({
            req(current_da())
            da <- current_da()
            pval_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(da))[1]
            fc_col <- intersect(c("log2FC", "log2FoldChange", "logFC"), colnames(da))[1]
            label_col <- intersect(c("metabolite", "name", "feature_id"), colnames(da))[1]
            if (is.na(pval_col) || is.na(fc_col)) return(NULL)
            create_volcano_plot(da, fc_col, pval_col, label_col, input$log2fc_thresh, input$pval_thresh)
        })

        output$da_table <- renderDT({
            req(current_da())
            datatable(current_da(), options = list(pageLength = 15, scrollX = TRUE),
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
            selectizeInput(ns("metabolite_select"), "Search:", choices = NULL,
                           options = list(placeholder = "Start typing...")),
            selectInput(ns("group_var"), "Group by:", choices = NULL)
        ),

        card(
            full_screen = TRUE,
            card_header(textOutput(ns("metab_title"))),
            card_body(plotlyOutput(ns("abundance_plot"), height = "800px"))
        )
    )
}

metabolomics_browser_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        observe({
            data <- data_reactive()
            req(data, data$normalized_matrix)
            mat <- data$normalized_matrix
            features <- if ("feature_id" %in% colnames(mat)) mat$feature_id
            else if (is.character(rownames(mat))) rownames(mat)
            else mat[[1]]
            updateSelectizeInput(session, "metabolite_select", choices = features, server = TRUE)
        })

        observe({
            data <- data_reactive()
            if (!is.null(data) && !is.null(data$metadata)) {
                cols <- get_metadata_columns(data$metadata)
                updateSelectInput(session, "group_var", choices = cols,
                                  selected = if(length(cols) > 0) cols[1] else NULL)
            }
        })

        output$metab_title <- renderText({
            req(input$metabolite_select)
            paste("Abundance of", input$metabolite_select)
        })

        output$abundance_plot <- renderPlotly({
            data <- data_reactive()
            req(input$metabolite_select, data, data$normalized_matrix)

            expr_df <- prepare_expression_data(
                data$normalized_matrix, input$metabolite_select,
                data$metadata, sample_col = "sample"
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
