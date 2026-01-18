# =============================================================================
# MS-Helios Integration (JAR Wrapper)
# =============================================================================

#' Run MS-Helios Visualization
#'
#' @param mae_data A MultiAssayExperiment object containing harmonized omics data.
#' @param config A list containing pipeline configuration.
#'
#' @return A list containing the output directory and plot path.
#' @export
run_ms_helios <- function(mae_data, config) {
    # Check if MS-Helios is enabled
    mh_config <- config$ms_helios %||% list()
    if (!(mh_config$run_ms_helios %||% FALSE)) {
        log_message("MS-Helios analysis disabled in config. Skipping.")
        return(NULL)
    }

    # Dependencies
    jar_path <- file.path("src/tools", "MS-Helios.jar")
    if (!file.exists(jar_path)) {
        stop("MS-Helios.jar not found at: ", jar_path)
    }

    circos_bin <- Sys.which("circos")
    if (circos_bin == "") {
        warning("Circos executable not found in PATH. MS-Helios plot cannot be generated.")
        return(NULL)
    }

    log_message("Starting MS-Helios visualization...")

    # Output directory
    out_dir <- file.path(config$output$output_dir, "ms_helios")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    # Tracks Configuration (Radii)
    # Outer to Inner
    # Layout: Ideogram (0.8r) -> Track 1 (0.75-0.95) -> Track 2 (0.50-0.70) -> Track 3 (0.25-0.45)
    # This assumes 3 tracks max for now or dynamic calculation.

    matrices <- extract_matrices_for_integration(mae_data)
    omics_layers <- names(matrices)
    n_omics <- length(omics_layers)

    # Define geometry
    # Start from 0.95r, go down by width + gap
    # Gap = 0.05
    track_width <- 0.20
    track_gap <- 0.05

    start_r <- 0.95

    # Temp directory for runs
    temp_run_dir <- file.path(out_dir, "temp_runs")
    if (!dir.exists(temp_run_dir)) dir.create(temp_run_dir)

    plots_config_str <- ""
    data_files_copy <- list()
    base_config <- NULL
    base_samples <- NULL

    for (i in seq_along(omics_layers)) {
        omic <- omics_layers[i]
        log_message("Processing MS-Helios track for: ", omic)

        # 1. Prepare Data
        # Extract
        mat <- matrices[[omic]]

        # Check if log transform needed? Assuming data is already suitable (normalized)
        # Scale? MS-Helios can normalize, but we can do z-score here if configured
        if (mh_config$tracks[[omic]]$scale %||% "none" == "zscore") {
            mat <- t(scale(t(mat)))
        }

        # Format for MS-Helios: ACCESSION, Sample1, Sample2...
        # Ensure sample order matches across all omics (which MAE guarantees usually if harmonized)
        # We must enforce column order to be identical to first omic for ideogram consistency
        if (i == 1) {
            sample_order <- colnames(mat)
        } else {
            # Subset/Reorder to match first omic
            shared <- intersect(sample_order, colnames(mat))
            if (length(shared) < length(sample_order)) {
                warning("Omic ", omic, " has fewer samples than reference. Filling/Matching.")
            }
            mat <- mat[, sample_order, drop = FALSE]
        }

        # Write CSV
        df <- as.data.frame(mat)
        df <- tibble::rownames_to_column(df, "ACCESSION")

        input_csv <- file.path(temp_run_dir, paste0(omic, ".csv"))
        write.csv(df, input_csv, row.names = FALSE, quote = FALSE)

        # 2. Run MS-Helios
        # -a histogram (default)
        # -i input.csv
        # -o specific_out

        curr_out <- file.path(temp_run_dir, paste0(omic, "_out"))
        # Clean prev
        if (dir.exists(curr_out)) unlink(curr_out, recursive = TRUE)
        dir.create(curr_out)

        cmd <- paste(
            "java -jar", shQuote(jar_path),
            "-i", shQuote(input_csv),
            "-o", shQuote(curr_out),
            "-a histogram"
        )

        res <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

        if (res != 0) {
            warning("MS-Helios JAR failed for ", omic)
            next
        }

        # 3. Parse Output
        mh_out_dir <- file.path(curr_out, "ms-helios")
        config_file <- file.path(mh_out_dir, "MS-Helios.config")

        if (!file.exists(config_file)) {
            warning("Config not generated for ", omic)
            next
        }

        # Read config
        conf_lines <- readLines(config_file)

        # Capture base config from the first successful run (ideogram, etc.)
        if (is.null(base_config)) {
            base_config <- conf_lines
            # Also copy samples.txt
            file.copy(file.path(mh_out_dir, "samples.txt"), file.path(out_dir, "samples.txt"), overwrite = TRUE)
        }

        # Extract <plots> block
        start_idx <- grep("<plots>", conf_lines)
        end_idx <- grep("</plots>", conf_lines)

        if (length(start_idx) > 0 && length(end_idx) > 0) {
            plot_block <- conf_lines[(start_idx + 1):(end_idx - 1)]

            # Determine new radius
            r1 <- start_r - (i - 1) * (track_width + track_gap)
            r0 <- r1 - track_width

            # Convert to string format "0.95r"
            r1_str <- sprintf("%.2fr", r1)
            r0_str <- sprintf("%.2fr", r0)

            # Replace radius in plot block
            # Pattern: r1 = ...
            plot_block <- gsub("r1\\s*=\\s*[0-9.]+r", paste0("r1 = ", r1_str), plot_block)
            plot_block <- gsub("r0\\s*=\\s*[0-9.]+r", paste0("r0 = ", r0_str), plot_block)

            # Handle file paths
            # The config refers to "datatrack_0.txt". We need to rename it to unique name and copy
            # Find file lines
            file_lines_idx <- grep("file\\s*=", plot_block)

            for (idx in file_lines_idx) {
                line <- plot_block[idx]
                # Extract filename
                orig_file <- gsub(".*file\\s*=\\s*", "", line)
                orig_file <- trimws(orig_file)

                # New filename
                new_file <- paste0(omic, "_", orig_file)

                # Update line
                plot_block[idx] <- sub(orig_file, new_file, line)

                # Copy data file
                src_dat <- file.path(mh_out_dir, orig_file)
                dst_dat <- file.path(out_dir, new_file)
                if (file.exists(src_dat)) {
                    file.copy(src_dat, dst_dat, overwrite = TRUE)
                }
            }

            # Append to master plots
            plots_config_str <- paste0(plots_config_str, "\n", paste(plot_block, collapse = "\n"))
        }
    }

    # 4. Construct Master Config
    if (is.null(base_config)) {
        warning("No valid MS-Helios configs were generated.")
        return(NULL)
    }

    # We iterate base_config and replace the <plots> block
    final_config <- c()
    in_plots <- FALSE
    for (line in base_config) {
        if (grepl("<plots>", line)) {
            final_config <- c(final_config, "<plots>")
            final_config <- c(final_config, plots_config_str)
            in_plots <- TRUE
            next
        }
        if (grepl("</plots>", line)) {
            final_config <- c(final_config, "</plots>")
            in_plots <- FALSE
            next
        }
        if (!in_plots) {
            final_config <- c(final_config, line)
        }
    }

    # Write Master Config
    master_conf_path <- file.path(out_dir, "MS-Helios.config")
    writeLines(final_config, master_conf_path)

    # 5. Run Circos
    # Need to copy 'etc' configs?
    # MS-Helios output seems to just rely on "include etc/..."
    # Usually /etc is relative. MS-Helios didn't output an etc folder in test_output...
    # Wait, MS-Helios.config lines: `<<include etc/image.conf>>`
    # If I don't have etc/ folder in out_dir, Circos will fail unless it finds it in standard paths or I assume MS-Helios generates it?
    # In my list_dir test_output/ms-helios, I ONLY saw config and samples.txt.
    # AND datatrack.
    # NO etc folder.
    # This implies MS-Helios expects to be run in a context where 'etc' exists?
    # Or maybe MS-Helios logic: "1. Copy ms-helios/ to the circos directory."
    # "2. Run in ms-helios/: perl ../bin/circos ..."
    # It implies it relies on Circos's own configuration files being reachable.
    # Since I am using `circos` from conda, `<<include etc/image.conf>>` might refer to Circos system defaults IF configured right,
    # or I need to provide them?
    # Usually standard Circos config files are local.
    # If I run `circos -conf MS-Helios.config` and `etc/` is not there, it checks default locations.
    # Let's hope Bioconda Circos is configured correctly.

    # Let's try running it.
    cwd_old <- getwd()
    setwd(out_dir)
    on.exit(setwd(cwd_old))

    cmd_circos <- paste(shQuote(circos_bin), "-conf MS-Helios.config")
    res_circos <- system(cmd_circos)

    if (res_circos == 0) {
        log_message("Circos plot generated successfully.")
    } else {
        warning("Circos execution failed.")
    }

    return(list(
        out_dir = out_dir,
        config = master_conf_path
    ))
}
