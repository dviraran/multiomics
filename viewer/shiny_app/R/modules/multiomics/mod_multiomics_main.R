# =============================================================================
# Multi-omics Main Module (LAZY LOADING)
# =============================================================================
# Data is loaded only when this tab is first accessed

multiomics_main_ui <- function(id) {
    ns <- NS(id)

    navset_card_tab(
        id = ns("multi_tabs"),

        nav_panel(
            title = "Overview",
            icon = icon("table"),
            multiomics_overview_ui(ns("overview"))
        ),

        nav_panel(
            title = "Cross-omics Correlations",
            icon = icon("chart-line"),
            multiomics_correlations_ui(ns("correlations"))
        ),

        nav_panel(
            title = "MOFA",
            icon = icon("layer-group"),
            multiomics_mofa_ui(ns("mofa"))
        ),

        nav_panel(
            title = "DIABLO",
            icon = icon("bullseye"),
            multiomics_diablo_ui(ns("diablo"))
        ),

        nav_panel(
            title = "RNA-Protein",
            icon = icon("dna"),
            multiomics_concordance_ui(ns("concordance"))
        )
    )
}

multiomics_main_server <- function(id, data_dir) {
    moduleServer(id, function(input, output, session) {

        # LAZY LOADING
        data <- reactiveVal(NULL)
        is_loaded <- reactiveVal(FALSE)

        load_data <- reactive({
            if (!is_loaded() && !is.null(data_dir)) {
                message("Loading Multi-omics data from: ", data_dir)
                loaded <- load_multiomics_data(data_dir)
                data(loaded)
                is_loaded(TRUE)
            }
            data()
        })

        multiomics_overview_server("overview", load_data)
        multiomics_correlations_server("correlations", load_data)
        multiomics_mofa_server("mofa", load_data)
        multiomics_diablo_server("diablo", load_data)
        multiomics_concordance_server("concordance", load_data)
    })
}

# =============================================================================
# Overview Sub-module
# =============================================================================

multiomics_overview_ui <- function(id) {
    ns <- NS(id)

    layout_columns(
        col_widths = c(4, 4, 4),

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

        card(
            card_header("MultiAssayExperiment Summary"),
            card_body(DTOutput(ns("mae_table")))
        )
    )
}

multiomics_overview_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        output$data_summary <- renderUI({
            data <- data_reactive()
            tagList(
                if (!is.null(data) && !is.null(data$mae_summary)) {
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
                    tags$p("No MAE summary available")
                }
            )
        })

        output$methods_summary <- renderUI({
            data <- data_reactive()
            if (is.null(data)) return(tags$p("Loading..."))
            tagList(
                tags$p(icon("check-circle", style = if(!is.null(data$has_mofa) && data$has_mofa) "color: green;" else "color: gray;"),
                       " MOFA2"),
                tags$p(icon("check-circle", style = if(!is.null(data$has_diablo) && data$has_diablo) "color: green;" else "color: gray;"),
                       " DIABLO"),
                tags$p(icon("check-circle", style = if(!is.null(data$has_correlations) && data$has_correlations) "color: green;" else "color: gray;"),
                       " Cross-omics Correlations")
            )
        })

        output$omics_summary <- renderUI({
            data <- data_reactive()
            if (is.null(data) || is.null(data$omics_types)) {
                return(tags$p("No omics types found"))
            }
            tagList(
                lapply(data$omics_types, function(x) tags$p(icon("circle", style = "color: #3498db;"), " ", x))
            )
        })

        output$mae_table <- renderDT({
            data <- data_reactive()
            req(data, data$mae_summary)
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
            selectInput(ns("omics_pair"), "Omics Pair:", choices = NULL),
            sliderInput(ns("cor_thresh"), "Min correlation:", min = 0, max = 1, value = 0.5, step = 0.1)
        ),

        layout_columns(
            col_widths = c(6, 6),
            card(full_screen = TRUE, card_header("Correlation Distribution"),
                 card_body(plotlyOutput(ns("cor_hist"), height = "400px"))),
            card(full_screen = TRUE, card_header("Top Correlations"),
                 card_body(DTOutput(ns("cor_table"))))
        )
    )
}

multiomics_correlations_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        observe({
            data <- data_reactive()
            if (!is.null(data) && !is.null(data$crossomics_correlations)) {
                pairs <- names(data$crossomics_correlations)
                updateSelectInput(session, "omics_pair", choices = pairs,
                                  selected = if(length(pairs) > 0) pairs[1] else NULL)
            }
        })

        current_cors <- reactive({
            data <- data_reactive()
            req(input$omics_pair, data, data$crossomics_correlations)
            data$crossomics_correlations[[input$omics_pair]]
        })

        output$cor_hist <- renderPlotly({
            req(current_cors())
            cors <- current_cors()
            cor_col <- intersect(c("correlation", "cor", "r", "spearman"), colnames(cors))[1]
            if (is.na(cor_col)) return(NULL)

            p <- ggplot(cors, aes_string(x = cor_col)) +
                geom_histogram(bins = 50, fill = "#3498db", alpha = 0.7) +
                geom_vline(xintercept = c(-input$cor_thresh, input$cor_thresh),
                           linetype = "dashed", color = "#e74c3c") +
                labs(x = "Correlation", y = "Count") +
                theme_minimal()
            ggplotly(p)
        })

        output$cor_table <- renderDT({
            req(current_cors())
            cors <- current_cors()
            cor_col <- intersect(c("correlation", "cor", "r", "spearman"), colnames(cors))[1]
            if (is.na(cor_col)) return(NULL)

            filtered <- cors[abs(cors[[cor_col]]) >= input$cor_thresh, ]
            filtered <- filtered[order(-abs(filtered[[cor_col]])), ]
            datatable(head(filtered, 100),
                      options = list(pageLength = 10, scrollX = TRUE),
                      rownames = FALSE)
        })
    })
}

# =============================================================================
# MOFA Sub-module
# =============================================================================

multiomics_mofa_ui <- function(id) {
    ns <- NS(id)

    layout_columns(
        col_widths = c(6, 6),

        card(full_screen = TRUE, card_header("Variance Explained by Factor"),
             card_body(plotlyOutput(ns("variance_heatmap"), height = "400px"))),

        card(full_screen = TRUE, card_header("Factor Scores"),
             card_body(
                 layout_sidebar(
                     sidebar = sidebar(
                         width = 150,
                         selectInput(ns("factor_x"), "X-axis:", choices = paste0("Factor", 1:10), selected = "Factor1"),
                         selectInput(ns("factor_y"), "Y-axis:", choices = paste0("Factor", 1:10), selected = "Factor2"),
                         selectInput(ns("color_by"), "Color:", choices = NULL)
                     ),
                     plotlyOutput(ns("factor_plot"), height = "350px")
                 )
             )),

        card(card_header("Top Feature Weights"),
             card_body(
                 selectInput(ns("weight_factor"), "Factor:", choices = NULL, width = "200px"),
                 DTOutput(ns("weights_table"))
             ))
    )
}

multiomics_mofa_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        observe({
            data <- data_reactive()
            if (!is.null(data)) {
                if (!is.null(data$mofa_factors)) {
                    factors <- grep("^Factor", colnames(data$mofa_factors), value = TRUE)
                    updateSelectInput(session, "factor_x", choices = factors, selected = factors[1])
                    updateSelectInput(session, "factor_y", choices = factors,
                                      selected = if(length(factors) > 1) factors[2] else factors[1])
                    updateSelectInput(session, "weight_factor", choices = factors, selected = factors[1])
                }
                if (!is.null(data$metadata)) {
                    cols <- get_metadata_columns(data$metadata)
                    updateSelectInput(session, "color_by", choices = cols,
                                      selected = if(length(cols) > 0) cols[1] else NULL)
                }
            }
        })

        output$variance_heatmap <- renderPlotly({
            data <- data_reactive()
            req(data, data$mofa_variance)

            var_df <- data$mofa_variance

            # Try to reshape if needed
            if ("factor" %in% colnames(var_df) && "view" %in% colnames(var_df)) {
                var_mat <- reshape2::dcast(var_df, factor ~ view, value.var = "variance_explained")
                rownames(var_mat) <- var_mat$factor
                var_mat <- as.matrix(var_mat[, -1])
            } else {
                var_mat <- as.matrix(var_df)
            }

            heatmaply(var_mat,
                      colors = viridis::viridis(100),
                      main = "% Variance Explained",
                      fontsize_row = 10,
                      fontsize_col = 10)
        })

        output$factor_plot <- renderPlotly({
            data <- data_reactive()
            req(data, data$mofa_factors, input$factor_x, input$factor_y)

            factors_df <- data$mofa_factors

            if (!is.null(input$color_by) &&
                !input$color_by %in% colnames(factors_df) &&
                !is.null(data$metadata)) {
                sample_col <- intersect(c("sample", "sample_id"), colnames(factors_df))[1]
                meta_col <- intersect(c("sample", "sample_id"), colnames(data$metadata))[1]
                if (!is.na(sample_col) && !is.na(meta_col)) {
                    factors_df <- merge(factors_df, data$metadata,
                                        by.x = sample_col, by.y = meta_col, all.x = TRUE)
                }
            }

            p <- ggplot(factors_df, aes_string(x = input$factor_x, y = input$factor_y))
            if (!is.null(input$color_by) && input$color_by %in% colnames(factors_df)) {
                p <- p + geom_point(aes_string(color = input$color_by), size = 3, alpha = 0.7)
            } else {
                p <- p + geom_point(size = 3, alpha = 0.7, color = "#3498db")
            }
            p <- p + theme_minimal()
            ggplotly(p)
        })

        output$weights_table <- renderDT({
            data <- data_reactive()
            req(data, data$mofa_weights, input$weight_factor)

            # Combine weights from all views
            all_weights <- do.call(rbind, lapply(names(data$mofa_weights), function(view) {
                w <- data$mofa_weights[[view]]
                if (!is.null(w)) {
                    w$view <- view
                    w
                } else { NULL }
            }))

            if (is.null(all_weights)) return(NULL)

            # Filter to selected factor
            weight_col <- input$weight_factor
            if (weight_col %in% colnames(all_weights)) {
                all_weights <- all_weights[order(-abs(all_weights[[weight_col]])), ]
            }

            datatable(head(all_weights, 50),
                      options = list(pageLength = 10, scrollX = TRUE),
                      rownames = FALSE)
        })
    })
}

# =============================================================================
# DIABLO Sub-module
# =============================================================================

multiomics_diablo_ui <- function(id) {
    ns <- NS(id)

    layout_columns(
        col_widths = c(6, 6),

        card(full_screen = TRUE, card_header("Sample Scores"),
             card_body(
                 selectInput(ns("score_view"), "View:", choices = NULL, width = "200px"),
                 plotlyOutput(ns("scores_plot"), height = "400px")
             )),

        card(full_screen = TRUE, card_header("Cross-validation Performance"),
             card_body(plotlyOutput(ns("cv_plot"), height = "400px"))),

        card(card_header("Selected Features"),
             card_body(DTOutput(ns("loadings_table"))))
    )
}

multiomics_diablo_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        observe({
            data <- data_reactive()
            if (!is.null(data) && !is.null(data$diablo_scores)) {
                views <- names(data$diablo_scores)
                updateSelectInput(session, "score_view", choices = views,
                                  selected = if(length(views) > 0) views[1] else NULL)
            }
        })

        output$scores_plot <- renderPlotly({
            data <- data_reactive()
            req(data, data$diablo_scores, input$score_view)

            scores <- data$diablo_scores[[input$score_view]]
            req(scores)

            comp_cols <- grep("^comp", colnames(scores), value = TRUE)
            if (length(comp_cols) < 2) return(NULL)

            p <- ggplot(scores, aes_string(x = comp_cols[1], y = comp_cols[2])) +
                geom_point(size = 3, alpha = 0.7, color = "#3498db") +
                labs(title = paste(input$score_view, "Scores")) +
                theme_minimal()
            ggplotly(p)
        })

        output$cv_plot <- renderPlotly({
            data <- data_reactive()
            req(data, data$diablo_cv)

            cv <- data$diablo_cv

            comp_col <- intersect(c("ncomp", "component", "n_components"), colnames(cv))[1]
            error_col <- intersect(c("error_rate", "error", "classification_error"), colnames(cv))[1]

            if (is.na(comp_col) || is.na(error_col)) return(NULL)

            p <- ggplot(cv, aes_string(x = comp_col, y = error_col)) +
                geom_line(color = "#3498db", linewidth = 1) +
                geom_point(color = "#e74c3c", size = 3) +
                labs(x = "Number of Components", y = "Error Rate") +
                theme_minimal()
            ggplotly(p)
        })

        output$loadings_table <- renderDT({
            data <- data_reactive()
            req(data, data$diablo_loadings)

            all_loadings <- do.call(rbind, lapply(names(data$diablo_loadings), function(view) {
                l <- data$diablo_loadings[[view]]
                if (!is.null(l)) {
                    l$view <- view
                    l
                } else { NULL }
            }))

            if (is.null(all_loadings)) return(NULL)

            datatable(head(all_loadings, 100),
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

    layout_columns(
        col_widths = c(6, 6),

        card(full_screen = TRUE, card_header("RNA vs Protein Fold Change"),
             card_body(plotlyOutput(ns("fc_scatter"), height = "450px"))),

        card(full_screen = TRUE, card_header("Correlation Distribution"),
             card_body(plotlyOutput(ns("cor_dist"), height = "450px"))),

        card(card_header("Concordance Summary"),
             card_body(DTOutput(ns("concordance_table"))))
    )
}

multiomics_concordance_server <- function(id, data_reactive) {
    moduleServer(id, function(input, output, session) {

        output$fc_scatter <- renderPlotly({
            data <- data_reactive()
            req(data, data$rna_protein_de_concordance)

            conc <- data$rna_protein_de_concordance

            rna_fc <- intersect(c("rna_log2FC", "rna_logFC", "log2FC_rna"), colnames(conc))[1]
            prot_fc <- intersect(c("prot_log2FC", "prot_logFC", "log2FC_prot"), colnames(conc))[1]

            if (is.na(rna_fc) || is.na(prot_fc)) return(NULL)

            p <- ggplot(conc, aes_string(x = rna_fc, y = prot_fc)) +
                geom_point(alpha = 0.5, color = "#3498db") +
                geom_smooth(method = "lm", se = TRUE, color = "#e74c3c") +
                geom_hline(yintercept = 0, linetype = "dashed") +
                geom_vline(xintercept = 0, linetype = "dashed") +
                labs(x = "RNA log2FC", y = "Protein log2FC") +
                theme_minimal()
            ggplotly(p)
        })

        output$cor_dist <- renderPlotly({
            data <- data_reactive()
            req(data, data$rna_protein_cors)

            cors <- data$rna_protein_cors
            cor_col <- intersect(c("correlation", "cor", "r", "spearman"), colnames(cors))[1]
            if (is.na(cor_col)) return(NULL)

            p <- ggplot(cors, aes_string(x = cor_col)) +
                geom_histogram(bins = 50, fill = "#3498db", alpha = 0.7) +
                geom_vline(xintercept = 0, linetype = "dashed") +
                labs(x = "RNA-Protein Correlation", y = "Count") +
                theme_minimal()
            ggplotly(p)
        })

        output$concordance_table <- renderDT({
            data <- data_reactive()
            req(data, data$rna_protein_de_concordance)
            datatable(head(data$rna_protein_de_concordance, 100),
                      options = list(pageLength = 15, scrollX = TRUE),
                      rownames = FALSE)
        })
    })
}
