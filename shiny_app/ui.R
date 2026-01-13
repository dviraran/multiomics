# =============================================================================
# Multi-Omics Viewer - User Interface
# =============================================================================

ui <- page_navbar(
    title = app_data$project_name,
    theme = app_theme,
    id = "main_nav",

    # Header with info
    header = tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),

    # ==========================================================================
    # Home/Overview Tab
    # ==========================================================================
    nav_panel(
        title = "Overview",
        icon = icon("house"),
        value = "overview",

        layout_columns(
            col_widths = c(4, 4, 4),

            # RNA-seq summary card
            card(
                card_header("RNA-seq Analysis", class = "bg-primary text-white"),
                card_body(
                    if (!is.null(app_data$rnaseq)) {
                        tagList(
                            tags$p(icon("check-circle", class = "text-success"),
                                   " Data loaded"),
                            tags$ul(
                                if (app_data$rnaseq$has_de) {
                                    tags$li(paste(length(app_data$rnaseq$contrasts), "contrast(s)"))
                                },
                                if (app_data$rnaseq$has_pathways) {
                                    tags$li("Pathway enrichment available")
                                },
                                if (app_data$rnaseq$has_pca) {
                                    tags$li("PCA results available")
                                }
                            )
                        )
                    } else {
                        tags$p(icon("times-circle", class = "text-muted"),
                               " No RNA-seq data")
                    }
                )
            ),

            # Proteomics summary card
            card(
                card_header("Proteomics Analysis", class = "bg-success text-white"),
                card_body(
                    if (!is.null(app_data$proteomics)) {
                        tagList(
                            tags$p(icon("check-circle", class = "text-success"),
                                   " Data loaded"),
                            tags$ul(
                                if (app_data$proteomics$has_de) {
                                    tags$li(paste(length(app_data$proteomics$contrasts), "contrast(s)"))
                                },
                                if (app_data$proteomics$has_ppi) {
                                    tags$li("PPI network available")
                                }
                            )
                        )
                    } else {
                        tags$p(icon("times-circle", class = "text-muted"),
                               " No Proteomics data")
                    }
                )
            ),

            # Multi-omics summary card
            card(
                card_header("Multi-omics Integration", class = "bg-info text-white"),
                card_body(
                    if (!is.null(app_data$multiomics)) {
                        tagList(
                            tags$p(icon("check-circle", class = "text-success"),
                                   " Data loaded"),
                            tags$ul(
                                if (length(app_data$multiomics$omics_types) > 0) {
                                    tags$li(paste("Omics:", paste(app_data$multiomics$omics_types, collapse = ", ")))
                                },
                                if (app_data$multiomics$has_mofa) {
                                    tags$li("MOFA results available")
                                },
                                if (app_data$multiomics$has_diablo) {
                                    tags$li("DIABLO results available")
                                },
                                if (app_data$multiomics$has_correlations) {
                                    tags$li("Cross-omics correlations available")
                                }
                            )
                        )
                    } else {
                        tags$p(icon("times-circle", class = "text-muted"),
                               " No Multi-omics data")
                    }
                )
            )
        ),

        # Instructions
        card(
            card_header("Getting Started"),
            card_body(
                tags$p("Use the navigation tabs above to explore your analysis results:"),
                tags$ul(
                    tags$li(tags$strong("RNA-seq:"), " Differential expression, pathway enrichment, PCA"),
                    tags$li(tags$strong("Proteomics:"), " Protein abundance, PPI networks"),
                    tags$li(tags$strong("Metabolomics:"), " Metabolite profiles, pathway context"),
                    tags$li(tags$strong("Multi-omics:"), " Cross-omics correlations, MOFA/DIABLO integration")
                ),
                tags$p(class = "text-muted",
                       "Click on points in plots to see details. Use dropdowns to select contrasts and color variables.")
            )
        )
    ),

    # ==========================================================================
    # RNA-seq Tab
    # ==========================================================================
    nav_panel(
        title = "RNA-seq",
        icon = icon("dna"),
        value = "rnaseq",

        if (!is.null(app_data$rnaseq)) {
            rnaseq_main_ui("rnaseq")
        } else {
            card(
                card_body(
                    tags$p(class = "text-center text-muted",
                           icon("info-circle"), " No RNA-seq data available")
                )
            )
        }
    ),

    # ==========================================================================
    # Proteomics Tab
    # ==========================================================================
    nav_panel(
        title = "Proteomics",
        icon = icon("cubes"),
        value = "proteomics",

        if (!is.null(app_data$proteomics)) {
            proteomics_main_ui("proteomics")
        } else {
            card(
                card_body(
                    tags$p(class = "text-center text-muted",
                           icon("info-circle"), " No Proteomics data available")
                )
            )
        }
    ),

    # ==========================================================================
    # Metabolomics Tab
    # ==========================================================================
    nav_panel(
        title = "Metabolomics",
        icon = icon("flask"),
        value = "metabolomics",

        if (!is.null(app_data$metabolomics)) {
            metabolomics_main_ui("metabolomics")
        } else {
            card(
                card_body(
                    tags$p(class = "text-center text-muted",
                           icon("info-circle"), " No Metabolomics data available")
                )
            )
        }
    ),

    # ==========================================================================
    # Multi-omics Tab
    # ==========================================================================
    nav_panel(
        title = "Multi-omics",
        icon = icon("project-diagram"),
        value = "multiomics",

        if (!is.null(app_data$multiomics)) {
            multiomics_main_ui("multiomics")
        } else {
            card(
                card_body(
                    tags$p(class = "text-center text-muted",
                           icon("info-circle"), " No Multi-omics integration data available")
                )
            )
        }
    ),

    # ==========================================================================
    # Footer
    # ==========================================================================
    nav_spacer(),
    nav_item(
        tags$a(
            href = "https://github.com/your-repo/multiomics",
            target = "_blank",
            icon("github"), " GitHub"
        )
    )
)
