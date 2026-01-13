# =============================================================================
# Multi-Omics Viewer - Server Logic
# =============================================================================

server <- function(input, output, session) {

    # ==========================================================================
    # RNA-seq Module
    # ==========================================================================
    if (!is.null(app_data$rnaseq)) {
        rnaseq_main_server("rnaseq", app_data$rnaseq)
    }

    # ==========================================================================
    # Proteomics Module
    # ==========================================================================
    if (!is.null(app_data$proteomics)) {
        proteomics_main_server("proteomics", app_data$proteomics)
    }

    # ==========================================================================
    # Metabolomics Module
    # ==========================================================================
    if (!is.null(app_data$metabolomics)) {
        metabolomics_main_server("metabolomics", app_data$metabolomics)
    }

    # ==========================================================================
    # Multi-omics Module
    # ==========================================================================
    if (!is.null(app_data$multiomics)) {
        multiomics_main_server("multiomics", app_data$multiomics)
    }
}
