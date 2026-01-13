# =============================================================================
# Multi-omics Main Module
# =============================================================================

multiomics_main_ui <- function(id) {
    ns <- NS(id)

    navset_card_tab(
        id = ns("multi_tabs"),

        # Overview Tab
        nav_panel(
            title = "Overview",
            icon = icon("table"),
            multiomics_overview_ui(ns("overview"))
        ),

        # Cross-omics Correlations Tab
        nav_panel(
            title = "Cross-omics Correlations",
            icon = icon("chart-line"),
            multiomics_correlations_ui(ns("correlations"))
        ),

        # MOFA Tab
        nav_panel(
            title = "MOFA",
            icon = icon("layer-group"),
            multiomics_mofa_ui(ns("mofa"))
        ),

        # DIABLO Tab
        nav_panel(
            title = "DIABLO",
            icon = icon("bullseye"),
            multiomics_diablo_ui(ns("diablo"))
        ),

        # Gene-Protein Concordance Tab
        nav_panel(
            title = "RNA-Protein",
            icon = icon("dna"),
            multiomics_concordance_ui(ns("concordance"))
        )
    )
}

multiomics_main_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {
        multiomics_overview_server("overview", data)
        multiomics_correlations_server("correlations", data)
        multiomics_mofa_server("mofa", data)
        multiomics_diablo_server("diablo", data)
        multiomics_concordance_server("concordance", data)
    })
}

# =============================================================================
# Overview Sub-module
# =============================================================================

multiomics_overview_ui <- function(id) {
    ns <- NS(id)

    layout_columns(
        col_widths = c(4, 4, 4),

        # Summary cards
        card(
            card_header("Data Summary", class = "bg-primary text-white"),
            card_body(uiOutput(ns("data_summary")))
        ),

        card(
            card_header("Integration Methods", class = "bg-success text-white"),
            card_body(uiOutput(ns("methods_summary")))
        ),

        card(
            card_header("Omics Types", class = "bg-info text-white"),
            card_body(uiOutput(ns("omics_summary")))
        ),

        # MAE summary table
        card(
            card_header("MultiAssayExperiment Summary"),
            card_body(DTOutput(ns("mae_table")))
        )
    )
}

multiomics_overview_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {

        output$data_summary <- renderUI({
            tagList(
                if (!is.null(data$mae_summary)) {
                    n_samples <- if ("n_samples" %in% colnames(data$mae_summary)) {
                        sum(data$mae_summary$n_samples)
                    } else { "N/A" }

                    n_features <- if ("n_features" %in% colnames(data$mae_summary)) {
                        sum(data$mae_summary$n_features)
                    } else { "N/A" }

                    tagList(
                        tags$p(tags$strong("Total Samples: "), n_samples),
                        tags$p(tags$strong("Total Features: "), format_number(n_features))
                    )
                } else {
                    tags$p("No summary data available")
                }
            )
        })

        output$methods_summary <- renderUI({
            tagList(
                tags$p(icon(if(data$has_mofa) "check-circle" else "times-circle",
                            class = if(data$has_mofa) "text-success" else "text-muted"),
                       " MOFA2"),
                tags$p(icon(if(data$has_diablo) "check-circle" else "times-circle",
                            class = if(data$has_diablo) "text-success" else "text-muted"),
                       " DIABLO"),
                tags$p(icon(if(data$has_correlations) "check-circle" else "times-circle",
                            class = if(data$has_correlations) "text-success" else "text-muted"),
                       " Cross-omics Correlations")
            )
        })

        output$omics_summary <- renderUI({
            if (length(data$omics_types) > 0) {
                tagList(
                    lapply(data$omics_types, function(o) {
                        tags$p(icon("check-circle", class = "text-success"), " ", o)
                    })
                )
            } else {
                tags$p("No omics data loaded")
            }
        })

        output$mae_table <- renderDT({
            req(data$mae_summary)
            datatable(data$mae_summary,
                      options = list(pageLength = 10, dom = "t"),
                      rownames = FALSE)
        })
    })
}

# =============================================================================
# Cross-omics Correlations Sub-module
# =============================================================================

multiomics_correlations_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectInput(ns("omics_pair"), "Omics Pair:",
                        choices = NULL),
            sliderInput(ns("cor_threshold"), "Min |correlation|:",
                        min = 0, max = 1, value = 0.3, step = 0.05),
            hr(),
            uiOutput(ns("cor_summary"))
        ),

        layout_columns(
            col_widths = c(6, 6),

            # Correlation distribution
            card(
                full_screen = TRUE,
                card_header("Correlation Distribution"),
                card_body(plotlyOutput(ns("cor_distribution"), height = "400px"))
            ),

            # Top correlations table
            card(
                card_header("Top Correlations"),
                card_body(DTOutput(ns("cor_table")))
            ),

            # Correlation heatmap
            card(
                full_screen = TRUE,
                card_header("Correlation Heatmap (Top Features)"),
                card_body(plotlyOutput(ns("cor_heatmap"), height = "450px"))
            )
        )
    )
}

multiomics_correlations_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Update omics pair choices
        observe({
            pairs <- c()
            if (!is.null(data$crossomics_correlations)) {
                pairs <- names(data$crossomics_correlations)
            }
            if (!is.null(data$rna_protein_cors)) {
                pairs <- unique(c(pairs, "RNA_vs_Protein"))
            }
            if (length(pairs) == 0) {
                pairs <- "No correlation data"
            }
            updateSelectInput(session, "omics_pair", choices = pairs, selected = pairs[1])
        })

        # Get current correlation data
        current_cors <- reactive({
            req(input$omics_pair)

            if (input$omics_pair == "RNA_vs_Protein" && !is.null(data$rna_protein_cors)) {
                return(data$rna_protein_cors)
            }

            if (!is.null(data$crossomics_correlations) &&
                input$omics_pair %in% names(data$crossomics_correlations)) {
                return(data$crossomics_correlations[[input$omics_pair]])
            }

            NULL
        })

        # Correlation summary
        output$cor_summary <- renderUI({
            req(current_cors())

            cors <- current_cors()
            cor_col <- intersect(c("correlation", "rna_protein_correlation", "cor"),
                                 colnames(cors))[1]

            if (is.na(cor_col)) return(NULL)

            cor_vals <- cors[[cor_col]]
            cor_vals <- cor_vals[abs(cor_vals) >= input$cor_threshold]

            tagList(
                tags$p(tags$strong("Features: "), format_number(length(cor_vals))),
                tags$p(tags$strong("Mean r: "), round(mean(cor_vals, na.rm = TRUE), 3)),
                tags$p(tags$strong("Median r: "), round(median(cor_vals, na.rm = TRUE), 3)),
                tags$p(tags$strong("% Positive: "),
                       round(100 * mean(cor_vals > 0, na.rm = TRUE), 1), "%")
            )
        })

        # Correlation distribution
        output$cor_distribution <- renderPlotly({
            req(current_cors())

            cors <- current_cors()
            cor_col <- intersect(c("correlation", "rna_protein_correlation", "cor"),
                                 colnames(cors))[1]

            if (is.na(cor_col)) return(NULL)

            # Filter
            cors <- cors[abs(cors[[cor_col]]) >= input$cor_threshold, ]

            p <- ggplot(cors, aes_string(x = cor_col)) +
                geom_histogram(bins = 50, fill = "#3498db", color = "white", alpha = 0.7) +
                geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
                geom_vline(xintercept = median(cors[[cor_col]], na.rm = TRUE),
                           color = "#2ecc71", linetype = "dashed") +
                labs(title = "Feature-level Correlation Distribution",
                     x = "Correlation", y = "Count") +
                theme_minimal()

            ggplotly(p)
        })

        # Correlation table
        output$cor_table <- renderDT({
            req(current_cors())

            cors <- current_cors()
            cor_col <- intersect(c("correlation", "rna_protein_correlation", "cor"),
                                 colnames(cors))[1]

            if (is.na(cor_col)) return(NULL)

            # Filter and sort
            cors <- cors[abs(cors[[cor_col]]) >= input$cor_threshold, ]
            cors <- cors[order(abs(cors[[cor_col]]), decreasing = TRUE), ]

            datatable(head(cors, 100),
                      options = list(pageLength = 15, scrollX = TRUE),
                      rownames = FALSE)
        })

        # Correlation heatmap
        output$cor_heatmap <- renderPlotly({
            req(current_cors())

            cors <- current_cors()
            cor_col <- intersect(c("correlation", "rna_protein_correlation", "cor"),
                                 colnames(cors))[1]

            if (is.na(cor_col)) return(NULL)

            # Get top features
            cors <- cors[order(abs(cors[[cor_col]]), decreasing = TRUE), ]
            top_cors <- head(cors, 50)

            # Create simple heatmap of top correlations
            p <- ggplot(top_cors, aes(x = 1, y = reorder(rownames(top_cors), get(cor_col)),
                                      fill = get(cor_col))) +
                geom_tile() +
                scale_fill_gradient2(low = "#3498db", mid = "white", high = "#e74c3c",
                                     midpoint = 0, name = "Correlation") +
                labs(x = "", y = "Feature") +
                theme_minimal() +
                theme(axis.text.x = element_blank())

            ggplotly(p)
        })
    })
}

# =============================================================================
# MOFA Sub-module
# =============================================================================

multiomics_mofa_ui <- function(id) {
    ns <- NS(id)

    if_else_ui <- function(condition, true_ui, false_ui) {
        if (condition) true_ui else false_ui
    }

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectInput(ns("factor_x"), "X-axis Factor:", choices = paste0("Factor", 1:10), selected = "Factor1"),
            selectInput(ns("factor_y"), "Y-axis Factor:", choices = paste0("Factor", 1:10), selected = "Factor2"),
            selectInput(ns("color_by"), "Color by:", choices = NULL),
            hr(),
            selectInput(ns("weights_omics"), "Weights for omics:", choices = NULL),
            sliderInput(ns("n_weights"), "Top weights to show:", min = 10, max = 50, value = 20)
        ),

        layout_columns(
            col_widths = c(6, 6),

            # Variance explained heatmap
            card(
                card_header("Variance Explained per Factor"),
                card_body(plotlyOutput(ns("variance_heatmap"), height = "300px"))
            ),

            # Factor scatter plot
            card(
                full_screen = TRUE,
                card_header("Sample Factor Scores"),
                card_body(plotlyOutput(ns("factor_plot"), height = "400px"))
            ),

            # Top weights
            card(
                card_header("Top Feature Weights"),
                card_body(plotlyOutput(ns("weights_plot"), height = "400px"))
            ),

            # Weights table
            card(
                card_header("Feature Weights Table"),
                card_body(DTOutput(ns("weights_table")))
            )
        )
    )
}

multiomics_mofa_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Check if MOFA data is available
        has_mofa <- reactive({
            !is.null(data$mofa_factors) && !is.null(data$mofa_variance)
        })

        # Update choices
        observe({
            if (!is.null(data$metadata)) {
                cols <- get_metadata_columns(data$metadata)
                updateSelectInput(session, "color_by", choices = cols,
                                  selected = if(length(cols) > 0) cols[1] else NULL)
            }

            if (!is.null(data$mofa_weights)) {
                omics <- names(data$mofa_weights)
                updateSelectInput(session, "weights_omics", choices = omics,
                                  selected = omics[1])
            }
        })

        # Variance explained heatmap
        output$variance_heatmap <- renderPlotly({
            req(data$mofa_variance)

            var_df <- data$mofa_variance

            # Convert to matrix format if needed
            if (!is.matrix(var_df)) {
                # Assume first column is factor names
                if (colnames(var_df)[1] %in% c("", "X", "factor", "Factor")) {
                    rownames(var_df) <- var_df[[1]]
                    var_df <- var_df[, -1]
                }
                var_mat <- as.matrix(var_df)
            } else {
                var_mat <- var_df
            }

            heatmaply(var_mat,
                      colors = colorRampPalette(c("white", "#3498db", "#e74c3c"))(100),
                      main = "",
                      fontsize_row = 10,
                      fontsize_col = 10)
        })

        # Factor scatter plot
        output$factor_plot <- renderPlotly({
            req(data$mofa_factors, input$factor_x, input$factor_y)

            factors <- data$mofa_factors

            # Get sample column
            sample_col <- intersect(c("sample", "sample_id"), colnames(factors))
            if (length(sample_col) == 0 && !is.null(rownames(factors))) {
                factors$sample <- rownames(factors)
                sample_col <- "sample"
            } else {
                sample_col <- sample_col[1]
            }

            # Merge with metadata
            if (!is.null(input$color_by) && !input$color_by %in% colnames(factors) &&
                !is.null(data$metadata)) {
                meta_sample_col <- intersect(c("sample", "sample_id"), colnames(data$metadata))[1]
                if (!is.na(meta_sample_col)) {
                    factors <- merge(factors, data$metadata,
                                     by.x = sample_col, by.y = meta_sample_col, all.x = TRUE)
                }
            }

            # Create plot
            p <- ggplot(factors, aes_string(x = input$factor_x, y = input$factor_y,
                                            text = sample_col))

            if (!is.null(input$color_by) && input$color_by %in% colnames(factors)) {
                p <- p + geom_point(aes_string(color = input$color_by), size = 3, alpha = 0.7)
            } else {
                p <- p + geom_point(color = "#3498db", size = 3, alpha = 0.7)
            }

            p <- p +
                labs(x = input$factor_x, y = input$factor_y) +
                theme_minimal()

            ggplotly(p, tooltip = c("text", "color"))
        })

        # Top weights bar plot
        output$weights_plot <- renderPlotly({
            req(data$mofa_weights, input$weights_omics)

            weights <- data$mofa_weights[[input$weights_omics]]
            if (is.null(weights)) return(NULL)

            # Get feature and Factor1 columns
            feature_col <- intersect(c("feature", "feature_id", "gene"), colnames(weights))
            if (length(feature_col) == 0) {
                if (!is.null(rownames(weights))) {
                    weights$feature <- rownames(weights)
                    feature_col <- "feature"
                } else {
                    feature_col <- colnames(weights)[1]
                }
            } else {
                feature_col <- feature_col[1]
            }

            factor_col <- intersect(c("Factor1", "factor1", "V1"), colnames(weights))
            if (length(factor_col) == 0) {
                # Use first numeric column
                num_cols <- which(sapply(weights, is.numeric))
                factor_col <- colnames(weights)[num_cols[1]]
            } else {
                factor_col <- factor_col[1]
            }

            # Get top weights
            weights <- weights[order(abs(weights[[factor_col]]), decreasing = TRUE), ]
            top_weights <- head(weights, input$n_weights)

            top_weights$direction <- ifelse(top_weights[[factor_col]] > 0, "positive", "negative")

            p <- ggplot(top_weights, aes_string(x = paste0("reorder(", feature_col, ", ", factor_col, ")"),
                                                y = factor_col, fill = "direction")) +
                geom_col(alpha = 0.8) +
                scale_fill_manual(values = c("positive" = "#e74c3c", "negative" = "#3498db")) +
                coord_flip() +
                labs(x = "Feature", y = "Weight (Factor 1)") +
                theme_minimal() +
                theme(legend.position = "none")

            ggplotly(p)
        })

        # Weights table
        output$weights_table <- renderDT({
            req(data$mofa_weights, input$weights_omics)

            weights <- data$mofa_weights[[input$weights_omics]]
            if (is.null(weights)) return(NULL)

            datatable(weights,
                      options = list(pageLength = 15, scrollX = TRUE),
                      rownames = FALSE)
        })
    })
}

# =============================================================================
# DIABLO Sub-module
# =============================================================================

multiomics_diablo_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectInput(ns("comp_x"), "X-axis Component:", choices = paste0("comp", 1:5), selected = "comp1"),
            selectInput(ns("comp_y"), "Y-axis Component:", choices = paste0("comp", 1:5), selected = "comp2"),
            selectInput(ns("color_by"), "Color by:", choices = NULL),
            hr(),
            selectInput(ns("loadings_omics"), "Loadings for omics:", choices = NULL)
        ),

        layout_columns(
            col_widths = c(6, 6),

            # Sample scores plot
            card(
                full_screen = TRUE,
                card_header("Sample Scores"),
                card_body(plotlyOutput(ns("scores_plot"), height = "450px"))
            ),

            # CV error plot
            card(
                card_header("Cross-validation Error"),
                card_body(plotlyOutput(ns("cv_plot"), height = "350px"))
            ),

            # Feature loadings
            card(
                card_header("Feature Loadings"),
                card_body(plotlyOutput(ns("loadings_plot"), height = "400px"))
            ),

            # Selected features table
            card(
                card_header("Selected Features"),
                card_body(DTOutput(ns("features_table")))
            )
        )
    )
}

multiomics_diablo_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Update choices
        observe({
            if (!is.null(data$metadata)) {
                cols <- get_metadata_columns(data$metadata)
                updateSelectInput(session, "color_by", choices = cols,
                                  selected = if(length(cols) > 0) cols[1] else NULL)
            }

            if (!is.null(data$diablo_loadings)) {
                omics <- names(data$diablo_loadings)
                updateSelectInput(session, "loadings_omics", choices = omics,
                                  selected = omics[1])
            }
        })

        # Sample scores plot
        output$scores_plot <- renderPlotly({
            req(data$diablo_scores)

            # Get first omics type scores (they should be similar)
            scores <- data$diablo_scores[[1]]
            if (is.null(scores)) return(NULL)

            sample_col <- intersect(c("sample", "sample_id", "rowname"), colnames(scores))
            if (length(sample_col) == 0) {
                scores$sample <- rownames(scores)
                sample_col <- "sample"
            } else {
                sample_col <- sample_col[1]
            }

            # Merge with metadata
            if (!is.null(input$color_by) && !input$color_by %in% colnames(scores) &&
                !is.null(data$metadata)) {
                meta_sample_col <- intersect(c("sample", "sample_id"), colnames(data$metadata))[1]
                if (!is.na(meta_sample_col)) {
                    scores <- merge(scores, data$metadata,
                                    by.x = sample_col, by.y = meta_sample_col, all.x = TRUE)
                }
            }

            p <- ggplot(scores, aes_string(x = input$comp_x, y = input$comp_y,
                                           text = sample_col))

            if (!is.null(input$color_by) && input$color_by %in% colnames(scores)) {
                p <- p + geom_point(aes_string(color = input$color_by), size = 3, alpha = 0.7)
            } else {
                p <- p + geom_point(color = "#3498db", size = 3, alpha = 0.7)
            }

            p <- p +
                labs(x = input$comp_x, y = input$comp_y) +
                theme_minimal()

            ggplotly(p, tooltip = c("text", "color"))
        })

        # CV error plot
        output$cv_plot <- renderPlotly({
            req(data$diablo_cv)

            cv <- data$diablo_cv

            if ("ncomp" %in% colnames(cv) && "error_rate" %in% colnames(cv)) {
                p <- ggplot(cv, aes(x = ncomp, y = error_rate)) +
                    geom_line(color = "#3498db") +
                    geom_point(size = 3, color = "#3498db") +
                    labs(x = "Number of Components", y = "Error Rate") +
                    theme_minimal()

                ggplotly(p)
            } else {
                # Just show the table as a plot placeholder
                NULL
            }
        })

        # Loadings plot
        output$loadings_plot <- renderPlotly({
            req(data$diablo_loadings, input$loadings_omics)

            loadings <- data$diablo_loadings[[input$loadings_omics]]
            if (is.null(loadings)) return(NULL)

            # Get feature and comp1 columns
            feature_col <- intersect(c("feature", "feature_id", "gene", "protein"),
                                     colnames(loadings))
            if (length(feature_col) == 0) {
                if (!is.null(rownames(loadings))) {
                    loadings$feature <- rownames(loadings)
                    feature_col <- "feature"
                } else {
                    feature_col <- colnames(loadings)[1]
                }
            } else {
                feature_col <- feature_col[1]
            }

            comp_col <- intersect(c("comp1", "Comp1", "V1"), colnames(loadings))
            if (length(comp_col) == 0) {
                num_cols <- which(sapply(loadings, is.numeric))
                comp_col <- colnames(loadings)[num_cols[1]]
            } else {
                comp_col <- comp_col[1]
            }

            # Filter non-zero loadings and get top
            loadings <- loadings[loadings[[comp_col]] != 0, ]
            loadings <- loadings[order(abs(loadings[[comp_col]]), decreasing = TRUE), ]
            top_loadings <- head(loadings, 20)

            top_loadings$direction <- ifelse(top_loadings[[comp_col]] > 0, "positive", "negative")

            p <- ggplot(top_loadings, aes_string(x = paste0("reorder(", feature_col, ", ", comp_col, ")"),
                                                 y = comp_col, fill = "direction")) +
                geom_col(alpha = 0.8) +
                scale_fill_manual(values = c("positive" = "#e74c3c", "negative" = "#3498db")) +
                coord_flip() +
                labs(x = "Feature", y = "Loading (Component 1)") +
                theme_minimal() +
                theme(legend.position = "none")

            ggplotly(p)
        })

        # Features table
        output$features_table <- renderDT({
            req(data$diablo_loadings, input$loadings_omics)

            loadings <- data$diablo_loadings[[input$loadings_omics]]
            if (is.null(loadings)) return(NULL)

            datatable(loadings,
                      options = list(pageLength = 15, scrollX = TRUE),
                      rownames = FALSE)
        })
    })
}

# =============================================================================
# RNA-Protein Concordance Sub-module
# =============================================================================

multiomics_concordance_ui <- function(id) {
    ns <- NS(id)

    layout_sidebar(
        sidebar = sidebar(
            width = 250,
            selectizeInput(ns("gene_select"), "Search Gene:",
                           choices = NULL,
                           options = list(placeholder = "Start typing...")),
            hr(),
            uiOutput(ns("concordance_summary"))
        ),

        layout_columns(
            col_widths = c(6, 6),

            # Fold change scatter
            card(
                full_screen = TRUE,
                card_header("RNA vs Protein Fold Change"),
                card_body(plotlyOutput(ns("fc_scatter"), height = "450px"))
            ),

            # Correlation distribution
            card(
                card_header("Gene-level Correlation Distribution"),
                card_body(plotlyOutput(ns("cor_distribution"), height = "400px"))
            ),

            # Selected gene expression
            card(
                card_header("Gene Expression (RNA vs Protein)"),
                card_body(plotlyOutput(ns("gene_expression"), height = "350px"))
            ),

            # Concordance table
            card(
                card_header("Concordance Table"),
                card_body(DTOutput(ns("concordance_table")))
            )
        )
    )
}

multiomics_concordance_server <- function(id, data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Update gene choices
        observe({
            genes <- c()
            if (!is.null(data$rna_protein_de_concordance)) {
                gene_col <- intersect(c("gene_symbol", "gene", "symbol"),
                                      colnames(data$rna_protein_de_concordance))[1]
                if (!is.na(gene_col)) {
                    genes <- unique(data$rna_protein_de_concordance[[gene_col]])
                }
            }
            if (!is.null(data$rna_protein_cors)) {
                gene_col <- intersect(c("gene_symbol", "gene", "symbol"),
                                      colnames(data$rna_protein_cors))[1]
                if (!is.na(gene_col)) {
                    genes <- unique(c(genes, data$rna_protein_cors[[gene_col]]))
                }
            }
            updateSelectizeInput(session, "gene_select", choices = genes, server = TRUE)
        })

        # Concordance summary
        output$concordance_summary <- renderUI({
            if (!is.null(data$rna_protein_de_concordance)) {
                conc <- data$rna_protein_de_concordance

                # Calculate concordance stats
                if ("concordant" %in% colnames(conc)) {
                    n_concordant <- sum(conc$concordant == TRUE, na.rm = TRUE)
                    n_total <- nrow(conc)
                    pct_concordant <- round(100 * n_concordant / n_total, 1)

                    return(tagList(
                        tags$p(tags$strong("Genes analyzed: "), format_number(n_total)),
                        tags$p(tags$strong("Concordant: "), format_number(n_concordant),
                               paste0(" (", pct_concordant, "%)"))
                    ))
                }
            }

            tags$p("No concordance data")
        })

        # Fold change scatter
        output$fc_scatter <- renderPlotly({
            req(data$rna_protein_de_concordance)

            conc <- data$rna_protein_de_concordance

            # Find columns
            rna_fc <- intersect(c("rna_log2FC", "rna_logFC", "RNA_log2FC"),
                                colnames(conc))[1]
            prot_fc <- intersect(c("protein_log2FC", "prot_log2FC", "Protein_log2FC"),
                                 colnames(conc))[1]
            gene_col <- intersect(c("gene_symbol", "gene", "symbol"),
                                  colnames(conc))[1]

            if (is.na(rna_fc) || is.na(prot_fc)) return(NULL)

            # Calculate correlation
            r <- cor(conc[[rna_fc]], conc[[prot_fc]], use = "pairwise.complete.obs")

            # Add concordance coloring
            if ("concordant" %in% colnames(conc)) {
                conc$concordance <- ifelse(conc$concordant, "Concordant", "Discordant")
            } else {
                conc$concordance <- "Unknown"
            }

            p <- ggplot(conc, aes_string(x = rna_fc, y = prot_fc,
                                         color = "concordance",
                                         text = if(!is.na(gene_col)) gene_col else NULL)) +
                geom_point(alpha = 0.6) +
                geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed",
                            inherit.aes = FALSE, aes_string(x = rna_fc, y = prot_fc)) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
                scale_color_manual(values = c("Concordant" = "#2ecc71", "Discordant" = "#e74c3c",
                                              "Unknown" = "#95a5a6")) +
                labs(title = paste0("RNA vs Protein Fold Change (r = ", round(r, 3), ")"),
                     x = "RNA log2(Fold Change)",
                     y = "Protein log2(Fold Change)") +
                theme_minimal()

            ggplotly(p, tooltip = "text")
        })

        # Correlation distribution
        output$cor_distribution <- renderPlotly({
            req(data$rna_protein_cors)

            cors <- data$rna_protein_cors
            cor_col <- intersect(c("correlation", "rna_protein_correlation", "cor"),
                                 colnames(cors))[1]

            if (is.na(cor_col)) return(NULL)

            p <- ggplot(cors, aes_string(x = cor_col)) +
                geom_histogram(bins = 40, fill = "#3498db", color = "white", alpha = 0.7) +
                geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
                geom_vline(xintercept = median(cors[[cor_col]], na.rm = TRUE),
                           color = "#2ecc71", linetype = "dashed") +
                labs(title = "Per-gene RNA-Protein Correlation",
                     x = "Correlation", y = "Count") +
                theme_minimal()

            ggplotly(p)
        })

        # Gene expression plot
        output$gene_expression <- renderPlotly({
            req(input$gene_select, data$rna_matrix, data$prot_matrix)

            gene <- input$gene_select

            # Get RNA expression
            rna_expr <- NULL
            if (!is.null(data$rna_matrix)) {
                rna_mat <- data$rna_matrix
                if (gene %in% rownames(rna_mat)) {
                    rna_expr <- as.numeric(rna_mat[gene, ])
                    names(rna_expr) <- colnames(rna_mat)
                } else if (gene %in% rna_mat[[1]]) {
                    idx <- which(rna_mat[[1]] == gene)
                    rna_expr <- as.numeric(rna_mat[idx, -1])
                    names(rna_expr) <- colnames(rna_mat)[-1]
                }
            }

            # Get Protein expression
            prot_expr <- NULL
            if (!is.null(data$prot_matrix)) {
                prot_mat <- data$prot_matrix
                if (gene %in% rownames(prot_mat)) {
                    prot_expr <- as.numeric(prot_mat[gene, ])
                    names(prot_expr) <- colnames(prot_mat)
                } else if (gene %in% prot_mat[[1]]) {
                    idx <- which(prot_mat[[1]] == gene)
                    prot_expr <- as.numeric(prot_mat[idx, -1])
                    names(prot_expr) <- colnames(prot_mat)[-1]
                }
            }

            if (is.null(rna_expr) && is.null(prot_expr)) {
                return(NULL)
            }

            # Find common samples
            common_samples <- intersect(names(rna_expr), names(prot_expr))

            if (length(common_samples) == 0) {
                return(NULL)
            }

            df <- data.frame(
                sample = common_samples,
                RNA = rna_expr[common_samples],
                Protein = prot_expr[common_samples]
            )

            r <- cor(df$RNA, df$Protein, use = "complete.obs")

            p <- ggplot(df, aes(x = RNA, y = Protein, text = sample)) +
                geom_point(size = 3, alpha = 0.7, color = "#3498db") +
                geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", linetype = "dashed") +
                labs(title = paste0(gene, ": RNA vs Protein (r = ", round(r, 3), ")"),
                     x = "RNA Expression",
                     y = "Protein Abundance") +
                theme_minimal()

            ggplotly(p, tooltip = "text")
        })

        # Concordance table
        output$concordance_table <- renderDT({
            req(data$rna_protein_de_concordance)

            datatable(data$rna_protein_de_concordance,
                      options = list(pageLength = 15, scrollX = TRUE),
                      rownames = FALSE,
                      filter = "top")
        })
    })
}
