# =============================================================================
# Multi-Omics Viewer - Main Application
# =============================================================================
# Interactive Shiny application for exploring multi-omics analysis results
#
# To run:
#   1. source("../install.R")  # One-time setup
#   2. source("../run_viewer.R")  # Launch app
#
# Or directly:
#   shiny::runApp("shiny_app")

library(shiny)

# Source global settings and load data
source("global.R")

# Source UI and server
source("ui.R")
source("server.R")

# Run the application
shinyApp(ui = ui, server = server)
