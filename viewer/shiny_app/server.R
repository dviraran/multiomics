# =============================================================================
# Server - Multi-Omics Viewer (LAZY LOADING)
# =============================================================================
# Data is loaded lazily in each module when the tab is first accessed

server <- function(input, output, session) {

    # Pass data directories to modules (not the actual data)
    # Each module handles its own lazy loading

    # RNA-seq module - receives path to data directory
    if (app_data$has_rnaseq) {
        rnaseq_main_server("rnaseq", app_data$rnaseq_dir)
    }

    # Proteomics module
    if (app_data$has_proteomics) {
        proteomics_main_server("proteomics", app_data$proteomics_dir)
    }

    # Metabolomics module
    if (app_data$has_metabolomics) {
        metabolomics_main_server("metabolomics", app_data$metabolomics_dir)
    }

    # Multi-omics module
    if (app_data$has_multiomics) {
        multiomics_main_server("multiomics", app_data$multiomics_dir)
    }
}
