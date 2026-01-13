# =============================================================================
# Create Distributable Viewer Package
# =============================================================================
# This script creates a standalone viewer package that can be shared with
# collaborators. It bundles the Shiny app code with pipeline output data.
#
# Usage:
#   source("create_viewer_package.R")
#
#   create_viewer_package(
#       output_dir = "my_project_viewer",
#       rnaseq_outputs = "rnaseq_pipeline/outputs",
#       proteomics_outputs = "proteomics_pipeline/outputs",
#       multiomics_outputs = "multiomics_pipeline/outputs",
#       project_name = "My Analysis Project"
#   )
# =============================================================================

#' Create a distributable viewer package
#'
#' @param output_dir Name of the output directory (will be created)
#' @param rnaseq_outputs Path to RNA-seq pipeline outputs (or NULL to skip)
#' @param proteomics_outputs Path to proteomics pipeline outputs (or NULL to skip)
#' @param metabolomics_outputs Path to metabolomics pipeline outputs (or NULL to skip)
#' @param multiomics_outputs Path to multi-omics pipeline outputs (or NULL to skip)
#' @param project_name Name to display in the viewer
#' @param create_zip Whether to create a zip file
#' @return Path to the created package
create_viewer_package <- function(
        output_dir = "multiomics_viewer_package",
        rnaseq_outputs = NULL,
        proteomics_outputs = NULL,
        metabolomics_outputs = NULL,
        multiomics_outputs = NULL,
        project_name = "Multi-Omics Analysis",
        create_zip = TRUE
) {

    cat("================================================================================\n")
    cat("Creating Multi-Omics Viewer Package\n")
    cat("================================================================================\n\n")

    # Get the directory containing this script (repo root)
    repo_root <- getwd()

    # Create output directory
    if (dir.exists(output_dir)) {
        cat("Warning: Output directory exists. Overwriting...\n")
        unlink(output_dir, recursive = TRUE)
    }
    dir.create(output_dir, recursive = TRUE)
    dir.create(file.path(output_dir, "data"), recursive = TRUE)

    cat("[1/5] Copying Shiny app code...\n")
    # Copy shiny_app directory
    shiny_source <- file.path(repo_root, "shiny_app")
    if (!dir.exists(shiny_source)) {
        stop("Cannot find shiny_app directory at: ", shiny_source)
    }
    file.copy(shiny_source, output_dir, recursive = TRUE)

    cat("[2/5] Copying launcher scripts...\n")
    # Copy install and run scripts
    file.copy(file.path(repo_root, "install_viewer.R"), output_dir)
    file.copy(file.path(repo_root, "run_viewer.R"), output_dir)

    cat("[3/5] Copying pipeline outputs...\n")
    # Copy RNA-seq outputs
    if (!is.null(rnaseq_outputs) && dir.exists(rnaseq_outputs)) {
        cat("  - Copying RNA-seq data...\n")
        rnaseq_dest <- file.path(output_dir, "data", "rnaseq")
        dir.create(rnaseq_dest, recursive = TRUE)
        copy_pipeline_outputs(rnaseq_outputs, rnaseq_dest, "rnaseq")
    }

    # Copy proteomics outputs
    if (!is.null(proteomics_outputs) && dir.exists(proteomics_outputs)) {
        cat("  - Copying Proteomics data...\n")
        prot_dest <- file.path(output_dir, "data", "proteomics")
        dir.create(prot_dest, recursive = TRUE)
        copy_pipeline_outputs(proteomics_outputs, prot_dest, "proteomics")
    }

    # Copy metabolomics outputs
    if (!is.null(metabolomics_outputs) && dir.exists(metabolomics_outputs)) {
        cat("  - Copying Metabolomics data...\n")
        metab_dest <- file.path(output_dir, "data", "metabolomics")
        dir.create(metab_dest, recursive = TRUE)
        copy_pipeline_outputs(metabolomics_outputs, metab_dest, "metabolomics")
    }

    # Copy multi-omics outputs
    if (!is.null(multiomics_outputs) && dir.exists(multiomics_outputs)) {
        cat("  - Copying Multi-omics data...\n")
        multi_dest <- file.path(output_dir, "data", "multiomics")
        dir.create(multi_dest, recursive = TRUE)
        copy_pipeline_outputs(multiomics_outputs, multi_dest, "multiomics")
    }

    cat("[4/5] Creating configuration file...\n")
    # Create config with project name
    config <- list(project_name = project_name)
    yaml::write_yaml(config, file.path(output_dir, "data", "config.yml"))

    # Create README
    create_readme(output_dir, project_name)

    cat("[5/5] Finalizing package...\n")
    # Create zip if requested
    zip_path <- NULL
    if (create_zip) {
        zip_path <- paste0(output_dir, ".zip")
        cat("  - Creating zip file:", zip_path, "\n")

        # Use zip command
        old_wd <- getwd()
        setwd(dirname(output_dir))
        zip(zip_path, basename(output_dir))
        setwd(old_wd)
    }

    cat("\n================================================================================\n")
    cat("Package created successfully!\n")
    cat("================================================================================\n\n")
    cat("Output directory:", normalizePath(output_dir), "\n")
    if (!is.null(zip_path)) {
        cat("Zip file:", normalizePath(zip_path), "\n")
    }
    cat("\nTo share with collaborators:\n")
    cat("1. Send them the zip file (or the folder)\n")
    cat("2. They unzip and open the folder in RStudio\n")
    cat("3. Run: source('install_viewer.R')  # one-time setup\n")
    cat("4. Run: source('run_viewer.R')      # to launch\n")

    invisible(output_dir)
}

#' Copy relevant outputs from a pipeline
#' @param source_dir Source directory
#' @param dest_dir Destination directory
#' @param pipeline_type Type of pipeline
copy_pipeline_outputs <- function(source_dir, dest_dir, pipeline_type) {

    # Define patterns of files to copy for each pipeline
    patterns <- switch(pipeline_type,
                       rnaseq = c(
                           "*.csv",
                           "de_results/*.csv",
                           "pathway_results/*.csv"
                       ),
                       proteomics = c(
                           "tables/*.csv",
                           "qc/*.csv",
                           "ppi_networks/*.csv"
                       ),
                       metabolomics = c(
                           "tables/*.csv",
                           "qc/*.csv",
                           "*.csv"
                       ),
                       multiomics = c(
                           "tables/*.csv",
                           "*.csv"
                       )
    )

    # Copy CSV files
    for (pattern in patterns) {
        # Handle subdirectories
        if (grepl("/", pattern)) {
            subdir <- dirname(pattern)
            file_pattern <- basename(pattern)

            src_subdir <- file.path(source_dir, subdir)
            if (dir.exists(src_subdir)) {
                dest_subdir <- file.path(dest_dir, subdir)
                dir.create(dest_subdir, recursive = TRUE, showWarnings = FALSE)

                files <- list.files(src_subdir, pattern = gsub("\\*", ".*", file_pattern),
                                    full.names = TRUE)
                for (f in files) {
                    file.copy(f, dest_subdir)
                }
            }
        } else {
            files <- list.files(source_dir, pattern = gsub("\\*", ".*", pattern),
                                full.names = TRUE)
            for (f in files) {
                if (file.info(f)$isdir == FALSE) {  # Don't copy directories
                    file.copy(f, dest_dir)
                }
            }
        }
    }
}

#' Create README file for the package
#' @param output_dir Output directory
#' @param project_name Project name
create_readme <- function(output_dir, project_name) {
    readme_content <- sprintf('
================================================================================
%s - Multi-Omics Viewer
================================================================================

This folder contains an interactive viewer for exploring your multi-omics
analysis results.

QUICK START
-----------
1. Open RStudio
2. Set this folder as your working directory:
   - Session -> Set Working Directory -> Choose Directory
   - Or: setwd("/path/to/this/folder")

3. Install required packages (one-time setup):
   source("install_viewer.R")

4. Launch the viewer:
   source("run_viewer.R")

5. The viewer will open in your web browser. Explore!


WHAT IS INCLUDED
-----------------
- shiny_app/     : The visualization application
- data/          : Your analysis results (CSV files)
- install_viewer.R : Package installer
- run_viewer.R     : Viewer launcher


TROUBLESHOOTING
---------------
- "Package not found" errors: Run source("install_viewer.R") again
- Blank plots: Check that data files exist in the data/ folder
- App won\'t start: Make sure you\'re in the correct directory

For more help, contact your bioinformatics collaborator.

Generated: %s
================================================================================
', project_name, Sys.time())

    writeLines(readme_content, file.path(output_dir, "README.txt"))
}

# Print usage info when sourced
cat("
================================================================================
Multi-Omics Viewer Package Creator
================================================================================

Usage:
  create_viewer_package(
      output_dir = \"my_project_viewer\",
      rnaseq_outputs = \"path/to/rnaseq/outputs\",
      proteomics_outputs = \"path/to/proteomics/outputs\",
      multiomics_outputs = \"path/to/multiomics/outputs\",
      project_name = \"My Project Name\"
  )

This will create a folder with everything needed to share the viewer.
================================================================================
")
