# =============================================================================
# Proteomics Main Module (LAZY LOADING)
# =============================================================================
# Data is loaded only when this tab is first accessed

proteomics_main_ui <- function(id) {
    ns <- NS(id)

    navset_card_tab(
        id = ns("prot_tabs"),

        # QC Tab
        nav_panel(
            title = "QC",
            icon = icon("chart-bar"),
            proteomics_qc_ui(ns("qc"))
        ),

        # PCA Tab
        nav_panel(
            title = "PCA",
            icon = icon("circle-dot"),
            proteomics_pca_ui(ns("pca"))
        ),

        # Differential Abundance Tab
        nav_panel(
            title = "Differential Abundance",
            icon = icon("volcano"),
            proteomics_da_ui(ns("da"))
        ),

        # PPI Network Tab
        nav_panel(
            title = "PPI Network",
            icon = icon("project-diagram"),
            proteomics_ppi_ui(ns("ppi"))
        ),

        # Protein Browser Tab
        nav_panel(
            title = "Protein Browser",
            icon = icon("search"),
            proteomics_browser_ui(ns("browser"))
        )
    )
}

proteomics_main_server <- function(id, data_dir) {
    moduleServer(id, function(input, output, session) {

        # LAZY LOADING: Only load data when this module is accessed
        data <- reactiveVal(NULL)
        is_loaded <- reactiveVal(FALSE)

        load_data <- reactive({
            if (!is_loaded() && !is.null(data_dir)) {
                message("Loading Proteomics data from: ", data_dir)
                loaded <- load_proteomics_data(data_dir)
                data(loaded)
                is_loaded(TRUE)
            }
            data()
        })

        proteomics_qc_server("qc", load_data)
        proteomics_pca_server("pca", load_data)
        proteomics_da_server("da", load_data)
        proteomics_ppi_server("ppi", load_data)
        proteomics_browser_server("browser", load_data)
    })
}

# =============================================================================
# QC Sub-module
# =============================================================================

proteomics_qc_ui <- function(id) {
    ns <- NS(id)

    tagList(
        # Missing values heatmap - full width
        card(
            full_screen = TRUE,
            card_header("Missing Value Pattern"),
            card_body(plotlyOutput(ns("missing_heatmap"), height = "450px"))
        ),

        # Sample correlation - full width
        card(
            full_screen = TRUE,
            card_header("Sample Correlation"),
            card_body(plotlyOutput(ns("correlation_heatmap"), height = "450px"))
        ),

        # QC metrics - full width
        card(
            card_header("QC Summary"),
            card_body(uiOutput(ns("qc_summary")))
        )
    )
}

proteomics_qc_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        # Missing values heatmap
        output$missing_heatmap <- renderPlotly({
            data <- data_reactive()
            req(data, data$normalized_matrix)

            mat <- data$normalized_matrix

            # Convert to matrix
            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "protein_id")) {
                    rownames(mat) <- mat[[1]]
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            # Create missing value indicator
            missing_mat <- ifelse(is.na(mat), 1, 0)

            # Subsample if too large
            if (nrow(missing_mat) > 500) {
                idx <- sample(1:nrow(missing_mat), 500)
                missing_mat <- missing_mat[idx, ]
            }

            heatmaply(missing_mat,
                      colors = c("white", "#e74c3c"),
                      main = "Missing Values (red = missing)",
                      fontsize_row = 6,
                      fontsize_col = 8,
                      showticklabels = c(FALSE, TRUE))
        })

        # Correlation heatmap
        output$correlation_heatmap <- renderPlotly({
            data <- data_reactive()
            req(data, data$normalized_matrix)

            mat <- data$normalized_matrix

            # Convert to matrix
            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "protein_id")) {
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            # Calculate correlation
            cor_mat <- cor(mat, use = "pairwise.complete.obs")

            heatmaply(cor_mat,
                      colors = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
                      main = "Sample Correlation",
                      fontsize_row = 8,
                      fontsize_col = 8)
        })

        # QC summary
        output$qc_summary <- renderUI({
            data <- data_reactive()
            req(data, data$normalized_matrix)

            mat <- data$normalized_matrix
            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "protein_id")) {
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            n_proteins <- nrow(mat)
            n_samples <- ncol(mat)
            pct_missing <- round(100 * sum(is.na(mat)) / length(mat), 1)

            tagList(
                tags$p(tags$strong("Proteins: "), format_number(n_proteins)),
                tags$p(tags$strong("Samples: "), format_number(n_samples)),
                tags$p(tags$strong("Missing values: "), paste0(pct_missing, "%"))
            )
        })
    })
}

# =============================================================================
# PCA Sub-module
# =============================================================================

proteomics_pca_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectInput(ns("pc_x"), "X-axis PC:", choices = paste0("PC", 1:10), selected = "PC1"),
            selectInput(ns("pc_y"), "Y-axis PC:", choices = paste0("PC", 1:10), selected = "PC2"),
            selectInput(ns("color_by"), "Color by:", choices = NULL),
            checkboxInput(ns("show_labels"), "Show sample labels", FALSE)
        ),

        card(
            full_screen = TRUE,
            card_header("PCA Plot"),
            card_body(plotlyOutput(ns("pca_plot"), height = "500px"))
        )
    )
}

proteomics_pca_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Compute PCA from normalized matrix if not provided
        pca_data <- reactive({
            data <- data_reactive()
            req(data)

            # Check if PCA results exist in QC data
            if (!is.null(data$qc) && !is.null(data$qc$pca_coordinates)) {
                return(data$qc$pca_coordinates)
            }

            # Otherwise compute from matrix
            req(data$normalized_matrix)

            mat <- data$normalized_matrix
            if (!is.matrix(mat)) {
                if (colnames(mat)[1] %in% c("", "X", "feature_id", "protein_id")) {
                    mat <- mat[, -1]
                }
                mat <- as.matrix(mat)
            }

            # Remove rows with too many NAs
            row_na <- rowMeans(is.na(mat))
            mat <- mat[row_na < 0.5, ]

            # Impute remaining NAs with column means
            for (i in 1:ncol(mat)) {
                mat[is.na(mat[, i]), i] <- mean(mat[, i], na.rm = TRUE)
            }

            # Compute PCA
            pca <- prcomp(t(mat), scale. = TRUE)

            # Create result data frame
            pca_df <- as.data.frame(pca$x[, 1:min(10, ncol(pca$x))])
            pca_df$sample <- rownames(pca_df)

            # Store variance explained
            attr(pca_df, "var_explained") <- 100 * pca$sdev^2 / sum(pca$sdev^2)

            pca_df
        })

        # Update color choices
        observe({
            req(pca_data())
            data <- data_reactive()

            cols <- setdiff(colnames(pca_data()), c(paste0("PC", 1:20), "sample"))

            if (!is.null(data) && !is.null(data$metadata)) {
                meta_cols <- get_metadata_columns(data$metadata)
                cols <- unique(c(cols, meta_cols))
            }

            updateSelectInput(session, "color_by", choices = cols,
                              selected = if(length(cols) > 0) cols[1] else NULL)
        })

        # PCA plot
        output$pca_plot <- renderPlotly({
            data <- data_reactive()
            req(pca_data(), input$pc_x, input$pc_y)

            pca_df <- pca_data()
            var_exp <- attr(pca_df, "var_explained")

            # Merge with metadata
            if (!is.null(input$color_by) &&
                !input$color_by %in% colnames(pca_df) &&
                !is.null(data) && !is.null(data$metadata)) {

                meta_sample_col <- intersect(c("sample", "sample_id"), colnames(data$metadata))[1]
                if (!is.na(meta_sample_col)) {
                    pca_df <- merge(pca_df, data$metadata,
                                    by.x = "sample", by.y = meta_sample_col,
                                    all.x = TRUE)
                }
            }

            create_pca_plot(
                pca_df,
                pc_x = input$pc_x,
                pc_y = input$pc_y,
                color_var = input$color_by,
                var_explained = var_exp,
                show_labels = input$show_labels,
                show_ellipse = TRUE
            )
        })
    })
}

# =============================================================================
# Differential Abundance Sub-module
# =============================================================================

proteomics_da_ui <- function(id) {
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
                card_header("MA Plot"),
                card_body(plotlyOutput(ns("ma_plot"), height = "450px"))
            ),

            card(
                card_header("Differential Abundance Results"),
                card_body(DTOutput(ns("da_table")))
            )
        )
    )
}

proteomics_da_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        observe({
            data <- data_reactive()
            req(data, data$de_results)
            contrasts <- get_contrasts(data$de_results)
            updateSelectInput(session, "contrast", choices = contrasts,
                              selected = contrasts[1])
        })

        current_da <- reactive({
            data <- data_reactive()
            req(input$contrast, data, data$de_results)
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
                         tags$p(tags$strong("Total proteins: "), format_number(summary$n_total)),
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
            label_col <- intersect(c("gene_symbol", "symbol", "protein_name", "feature_id"),
                                   colnames(da))[1]

            if (is.na(pval_col) || is.na(fc_col)) return(NULL)

            create_volcano_plot(da, fc_col, pval_col, label_col,
                                input$log2fc_thresh, input$pval_thresh,
                                paste("Volcano:", input$contrast))
        })

        output$ma_plot <- renderPlotly({
            req(current_da())

            da <- current_da()
            fc_col <- intersect(c("log2FC", "log2FoldChange", "logFC"), colnames(da))[1]
            mean_col <- intersect(c("AveExpr", "baseMean", "avg_intensity"), colnames(da))[1]
            pval_col <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(da))[1]

            if (is.na(fc_col) || is.na(mean_col)) return(NULL)

            da$significance <- "ns"
            da$significance[da[[pval_col]] < input$pval_thresh &
                              da[[fc_col]] > input$log2fc_thresh] <- "up"
            da$significance[da[[pval_col]] < input$pval_thresh &
                              da[[fc_col]] < -input$log2fc_thresh] <- "down"

            p <- ggplot(da, aes_string(x = mean_col, y = fc_col, color = "significance")) +
                geom_point(alpha = 0.5, size = 1) +
                scale_color_manual(values = c("up" = "#e74c3c", "down" = "#3498db", "ns" = "#95a5a6")) +
                geom_hline(yintercept = 0, linetype = "dashed") +
                labs(x = "Average Expression", y = "log2 Fold Change") +
                theme_minimal()

            ggplotly(p)
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
# PPI Network Sub-module
# =============================================================================

proteomics_ppi_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            sliderInput(ns("edge_weight"), "Min edge weight:",
                        min = 0, max = 1, value = 0.4, step = 0.1),
            selectInput(ns("node_color"), "Color nodes by:",
                        choices = c("Community", "Degree", "Significance")),
            selectInput(ns("node_size"), "Size nodes by:",
                        choices = c("Degree", "Betweenness", "Fixed")),
            hr(),
            uiOutput(ns("network_stats"))
        ),

        layout_columns(
            col_widths = c(8, 4),

            card(
                full_screen = TRUE,
                card_header("PPI Network"),
                card_body(visNetworkOutput(ns("network_plot"), height = "550px"))
            ),

            card(
                card_header("Hub Proteins"),
                card_body(DTOutput(ns("hub_table")))
            )
        )
    )
}

proteomics_ppi_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Network stats
        output$network_stats <- renderUI({
            data <- data_reactive()
            if (is.null(data) || is.null(data$ppi_edges)) {
                return(tags$p("No PPI network data available"))
            }

            edges <- data$ppi_edges
            n_nodes <- length(unique(c(edges[[1]], edges[[2]])))
            n_edges <- nrow(edges)

            tagList(
                tags$p(tags$strong("Nodes: "), format_number(n_nodes)),
                tags$p(tags$strong("Edges: "), format_number(n_edges))
            )
        })

        # Network visualization
        output$network_plot <- renderVisNetwork({
            data <- data_reactive()
            req(data, data$ppi_edges)

            edges <- data$ppi_edges

            # Rename columns
            if (!"from" %in% colnames(edges)) {
                colnames(edges)[1] <- "from"
                colnames(edges)[2] <- "to"
            }

            # Filter by weight if available
            if ("weight" %in% colnames(edges) || "combined_score" %in% colnames(edges)) {
                weight_col <- intersect(c("weight", "combined_score"), colnames(edges))[1]
                edges <- edges[edges[[weight_col]] >= input$edge_weight, ]
            }

            if (nrow(edges) == 0) {
                return(NULL)
            }

            # Create nodes
            all_nodes <- unique(c(edges$from, edges$to))
            nodes <- data.frame(id = all_nodes, label = all_nodes, stringsAsFactors = FALSE)

            # Add node attributes
            if (!is.null(data$ppi_nodes)) {
                node_data <- data$ppi_nodes
                id_col <- intersect(c("id", "protein_id", "node"), colnames(node_data))[1]

                if (!is.na(id_col)) {
                    nodes <- merge(nodes, node_data, by.x = "id", by.y = id_col, all.x = TRUE)
                }
            }

            # Add community info
            if (!is.null(data$communities)) {
                comm <- data$communities
                id_col <- intersect(c("protein_id", "id", "node"), colnames(comm))[1]
                comm_col <- intersect(c("community", "community_id", "module"), colnames(comm))[1]

                if (!is.na(id_col) && !is.na(comm_col)) {
                    nodes <- merge(nodes, comm[, c(id_col, comm_col)],
                                   by.x = "id", by.y = id_col, all.x = TRUE)
                }
            }

            # Apply coloring
            if (input$node_color == "Community" && "community" %in% colnames(nodes)) {
                n_comm <- length(unique(nodes$community))
                comm_colors <- categorical_colors[1:min(n_comm, 10)]
                nodes$color <- comm_colors[as.numeric(factor(nodes$community))]
            } else if (input$node_color == "Degree" && "degree" %in% colnames(nodes)) {
                nodes$color <- colorRampPalette(c("#3498db", "#e74c3c"))(100)[
                    as.integer(cut(nodes$degree, 100))
                ]
            } else {
                nodes$color <- "#3498db"
            }

            # Apply sizing
            if (input$node_size == "Degree" && "degree" %in% colnames(nodes)) {
                nodes$size <- scales::rescale(nodes$degree, to = c(10, 40))
            } else if (input$node_size == "Betweenness" && "betweenness" %in% colnames(nodes)) {
                nodes$size <- scales::rescale(nodes$betweenness, to = c(10, 40))
            } else {
                nodes$size <- 20
            }

            # Create network
            visNetwork(nodes, edges) %>%
                visOptions(
                    highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1),
                    nodesIdSelection = TRUE
                ) %>%
                visPhysics(stabilization = list(iterations = 100)) %>%
                visInteraction(navigationButtons = TRUE) %>%
                visLayout(randomSeed = 42)
        })

        # Hub proteins table
        output$hub_table <- renderDT({
            data <- data_reactive()
            req(data, data$hub_proteins)

            datatable(
                head(data$hub_proteins, 20),
                options = list(pageLength = 10, dom = "t"),
                rownames = FALSE
            )
        })
    })
}

# =============================================================================
# Protein Browser Sub-module
# =============================================================================

proteomics_browser_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectizeInput(ns("protein_select"), "Search Protein:",
                           choices = NULL,
                           options = list(placeholder = "Start typing...")),
            selectInput(ns("group_var"), "Group by:", choices = NULL)
        ),

        card(
            full_screen = TRUE,
            card_header(textOutput(ns("protein_title"))),
            card_body(plotlyOutput(ns("abundance_plot"), height = "400px"))
        )
    )
}

proteomics_browser_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        observe({
            data <- data_reactive()
            req(data, data$normalized_matrix)

            mat <- data$normalized_matrix
            proteins <- if ("feature_id" %in% colnames(mat)) mat$feature_id
            else if ("protein_id" %in% colnames(mat)) mat$protein_id
            else if (is.character(rownames(mat))) rownames(mat)
            else mat[[1]]

            if (!is.null(data$feature_annotations)) {
                symbol_col <- intersect(c("gene_symbol", "symbol"), colnames(data$feature_annotations))[1]
                if (!is.na(symbol_col)) {
                    symbols <- data$feature_annotations[[symbol_col]]
                    proteins <- unique(c(symbols[!is.na(symbols)], proteins))
                }
            }

            updateSelectizeInput(session, "protein_select", choices = proteins, server = TRUE)
        })

        observe({
            data <- data_reactive()
            if (!is.null(data) && !is.null(data$metadata)) {
                cols <- get_metadata_columns(data$metadata)
                updateSelectInput(session, "group_var", choices = cols,
                                  selected = if(length(cols) > 0) cols[1] else NULL)
            }
        })

        output$protein_title <- renderText({
            req(input$protein_select)
            paste("Abundance of", input$protein_select)
        })

        output$abundance_plot <- renderPlotly({
            data <- data_reactive()
            req(input$protein_select, data, data$normalized_matrix)

            expr_df <- prepare_expression_data(
                data$normalized_matrix,
                input$protein_select,
                data$metadata,
                sample_col = "sample"
            )

            if (is.null(expr_df)) return(NULL)

            group_var <- if (!is.null(input$group_var) && input$group_var %in% colnames(expr_df)) {
                input$group_var
            } else {
                "sample"
            }

            create_expression_boxplot(expr_df, input$protein_select, group_var, "expression")
        })
    })
}
