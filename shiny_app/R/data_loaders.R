# =============================================================================
# Data Loading Functions for Multi-Omics Viewer
# =============================================================================
# Functions to load CSV files from each pipeline's output directory

#' Safe CSV reader with error handling
#' @param file_path Path to CSV file
#' @return Data frame or NULL if file doesn't exist
safe_read_csv <- function(file_path) {
    if (!file.exists(file_path)) {
        return(NULL)
    }
    tryCatch({
        read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
        warning("Failed to read: ", file_path, " - ", e$message)
        NULL
    })
}

#' Load all CSV files matching a pattern
#' @param dir_path Directory to search
#' @param pattern Filename pattern (regex)
#' @return Named list of data frames
load_csv_files <- function(dir_path, pattern = "\\.csv$") {
    if (!dir.exists(dir_path)) return(list())

    files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) return(list())

    result <- lapply(files, safe_read_csv)
    names(result) <- tools::file_path_sans_ext(basename(files))

    # Remove NULL entries
    result[!sapply(result, is.null)]
}

# =============================================================================
# RNA-seq Data Loader
# =============================================================================

#' Load RNA-seq pipeline outputs
#' @param data_dir Path to rnaseq data directory
#' @return List containing all RNA-seq data
load_rnaseq_data <- function(data_dir) {
    data <- list()

    # QC metrics
    data$qc_metrics <- safe_read_csv(file.path(data_dir, "qc_metrics.csv"))

    # PCA results
    data$pca_results <- safe_read_csv(file.path(data_dir, "pca_results.csv"))

    # Sample correlation
    data$sample_correlation <- safe_read_csv(file.path(data_dir, "sample_correlation.csv"))

    # Normalized counts
    data$normalized_counts <- safe_read_csv(file.path(data_dir, "normalized_counts.csv"))
    if (is.null(data$normalized_counts)) {
        data$normalized_counts <- safe_read_csv(file.path(data_dir, "vst_counts.csv"))
    }

    # Sample metadata
    data$metadata <- safe_read_csv(file.path(data_dir, "sample_metadata.csv"))

    # Gene annotation
    data$gene_annotation <- safe_read_csv(file.path(data_dir, "gene_annotation.csv"))

    # DE results (multiple files)
    de_dir <- file.path(data_dir, "de_results")
    if (dir.exists(de_dir)) {
        data$de_results <- load_csv_files(de_dir)
    } else {
        # Try loading de_results_*.csv from main directory
        de_files <- list.files(data_dir, pattern = "^de_results.*\\.csv$", full.names = TRUE)
        if (length(de_files) > 0) {
            data$de_results <- lapply(de_files, safe_read_csv)
            names(data$de_results) <- gsub("de_results_?", "", tools::file_path_sans_ext(basename(de_files)))
        }
    }

    # Pathway results (multiple files)
    pathway_dir <- file.path(data_dir, "pathway_results")
    if (dir.exists(pathway_dir)) {
        data$pathway_results <- load_csv_files(pathway_dir)
    } else {
        # Try loading pathway_*.csv from main directory
        pathway_files <- list.files(data_dir, pattern = "^pathway.*\\.csv$", full.names = TRUE)
        if (length(pathway_files) > 0) {
            data$pathway_results <- lapply(pathway_files, safe_read_csv)
            names(data$pathway_results) <- tools::file_path_sans_ext(basename(pathway_files))
        }
    }

    # Extract contrast names
    if (!is.null(data$de_results) && length(data$de_results) > 0) {
        data$contrasts <- names(data$de_results)
    } else {
        data$contrasts <- character(0)
    }

    # Check what data is available
    data$has_qc <- !is.null(data$qc_metrics)
    data$has_pca <- !is.null(data$pca_results)
    data$has_de <- !is.null(data$de_results) && length(data$de_results) > 0
    data$has_pathways <- !is.null(data$pathway_results) && length(data$pathway_results) > 0
    data$has_counts <- !is.null(data$normalized_counts)

    data
}

# =============================================================================
# Proteomics Data Loader
# =============================================================================

#' Load Proteomics pipeline outputs
#' @param data_dir Path to proteomics data directory
#' @return List containing all proteomics data
load_proteomics_data <- function(data_dir) {
    data <- list()

    # Check for tables subdirectory
    tables_dir <- file.path(data_dir, "tables")
    if (!dir.exists(tables_dir)) {
        tables_dir <- data_dir
    }

    # Normalized matrix
    data$normalized_matrix <- safe_read_csv(file.path(tables_dir, "normalized_matrix.csv"))

    # Feature annotations
    data$feature_annotations <- safe_read_csv(file.path(tables_dir, "feature_annotations.csv"))

    # Sample metadata
    data$metadata <- safe_read_csv(file.path(data_dir, "sample_metadata.csv"))
    if (is.null(data$metadata)) {
        data$metadata <- safe_read_csv(file.path(tables_dir, "sample_metadata.csv"))
    }

    # DE results
    de_files <- list.files(tables_dir, pattern = "^differential.*\\.csv$", full.names = TRUE)
    if (length(de_files) > 0) {
        data$de_results <- lapply(de_files, safe_read_csv)
        names(data$de_results) <- gsub("differential_?", "", tools::file_path_sans_ext(basename(de_files)))
    }

    # QC data
    qc_dir <- file.path(data_dir, "qc")
    if (dir.exists(qc_dir)) {
        data$qc <- load_csv_files(qc_dir)
    }

    # PPI networks
    ppi_dir <- file.path(data_dir, "ppi_networks")
    if (dir.exists(ppi_dir)) {
        data$ppi_edges <- safe_read_csv(file.path(ppi_dir, "ppi_network_edges.csv"))
        data$ppi_nodes <- safe_read_csv(file.path(ppi_dir, "ppi_network_nodes.csv"))
        data$hub_proteins <- safe_read_csv(file.path(ppi_dir, "hub_proteins.csv"))
        data$communities <- safe_read_csv(file.path(ppi_dir, "community_membership.csv"))
    }

    # Extract contrast names
    if (!is.null(data$de_results) && length(data$de_results) > 0) {
        data$contrasts <- names(data$de_results)
    } else {
        data$contrasts <- character(0)
    }

    # Check what data is available
    data$has_matrix <- !is.null(data$normalized_matrix)
    data$has_de <- !is.null(data$de_results) && length(data$de_results) > 0
    data$has_ppi <- !is.null(data$ppi_edges)

    data
}

# =============================================================================
# Metabolomics Data Loader
# =============================================================================

#' Load Metabolomics pipeline outputs
#' @param data_dir Path to metabolomics data directory
#' @return List containing all metabolomics data
load_metabolomics_data <- function(data_dir) {
    data <- list()

    # Check for tables subdirectory
    tables_dir <- file.path(data_dir, "tables")
    if (!dir.exists(tables_dir)) {
        tables_dir <- data_dir
    }

    # Normalized matrix
    data$normalized_matrix <- safe_read_csv(file.path(tables_dir, "normalized_matrix.csv"))

    # PCA
    data$pca_scores <- safe_read_csv(file.path(tables_dir, "pca_scores.csv"))
    if (is.null(data$pca_scores)) {
        data$pca_scores <- safe_read_csv(file.path(data_dir, "pca_scores.csv"))
    }
    data$pca_loadings <- safe_read_csv(file.path(tables_dir, "pca_loadings.csv"))

    # Sample metadata
    data$metadata <- safe_read_csv(file.path(data_dir, "sample_metadata.csv"))
    if (is.null(data$metadata)) {
        data$metadata <- safe_read_csv(file.path(tables_dir, "sample_metadata.csv"))
    }

    # DE results
    de_files <- list.files(tables_dir, pattern = "^differential.*\\.csv$", full.names = TRUE)
    if (length(de_files) > 0) {
        data$de_results <- lapply(de_files, safe_read_csv)
        names(data$de_results) <- gsub("differential_?", "", tools::file_path_sans_ext(basename(de_files)))
    }

    # QC data
    qc_dir <- file.path(data_dir, "qc")
    if (dir.exists(qc_dir)) {
        data$qc <- load_csv_files(qc_dir)
    }

    # Extract contrast names
    if (!is.null(data$de_results) && length(data$de_results) > 0) {
        data$contrasts <- names(data$de_results)
    } else {
        data$contrasts <- character(0)
    }

    # Check what data is available
    data$has_matrix <- !is.null(data$normalized_matrix)
    data$has_pca <- !is.null(data$pca_scores)
    data$has_de <- !is.null(data$de_results) && length(data$de_results) > 0

    data
}

# =============================================================================
# Multi-omics Data Loader
# =============================================================================

#' Load Multi-omics pipeline outputs
#' @param data_dir Path to multiomics data directory
#' @return List containing all multi-omics integration data
load_multiomics_data <- function(data_dir) {
    data <- list()

    # Check for tables subdirectory
    tables_dir <- file.path(data_dir, "tables")
    if (!dir.exists(tables_dir)) {
        tables_dir <- data_dir
    }

    # MAE summary
    data$mae_summary <- safe_read_csv(file.path(tables_dir, "mae_summary.csv"))

    # Sample alignment
    data$sample_alignment <- safe_read_csv(file.path(tables_dir, "sample_alignment.csv"))

    # Metadata
    data$metadata <- safe_read_csv(file.path(data_dir, "sample_metadata.csv"))
    if (is.null(data$metadata)) {
        data$metadata <- safe_read_csv(file.path(tables_dir, "sample_metadata.csv"))
    }

    # MOFA results
    data$mofa_factors <- safe_read_csv(file.path(tables_dir, "mofa_factors.csv"))
    data$mofa_variance <- safe_read_csv(file.path(tables_dir, "mofa_variance_explained.csv"))

    # MOFA weights per omics
    mofa_weight_files <- list.files(tables_dir, pattern = "^mofa_weights.*\\.csv$", full.names = TRUE)
    if (length(mofa_weight_files) > 0) {
        data$mofa_weights <- lapply(mofa_weight_files, safe_read_csv)
        names(data$mofa_weights) <- gsub("mofa_weights_?", "", tools::file_path_sans_ext(basename(mofa_weight_files)))
    }

    # DIABLO results
    diablo_score_files <- list.files(tables_dir, pattern = "^diablo_scores.*\\.csv$", full.names = TRUE)
    if (length(diablo_score_files) > 0) {
        data$diablo_scores <- lapply(diablo_score_files, safe_read_csv)
        names(data$diablo_scores) <- gsub("diablo_scores_?", "", tools::file_path_sans_ext(basename(diablo_score_files)))
    }

    diablo_loading_files <- list.files(tables_dir, pattern = "^diablo_loadings.*\\.csv$", full.names = TRUE)
    if (length(diablo_loading_files) > 0) {
        data$diablo_loadings <- lapply(diablo_loading_files, safe_read_csv)
        names(data$diablo_loadings) <- gsub("diablo_loadings_?", "", tools::file_path_sans_ext(basename(diablo_loading_files)))
    }

    data$diablo_cv <- safe_read_csv(file.path(tables_dir, "diablo_cv_error_rates.csv"))

    # Cross-omics correlations
    cor_files <- list.files(tables_dir, pattern = "^crossomics_correlations.*\\.csv$", full.names = TRUE)
    if (length(cor_files) > 0) {
        data$crossomics_correlations <- lapply(cor_files, safe_read_csv)
        names(data$crossomics_correlations) <- gsub("crossomics_correlations_?", "",
                                                     tools::file_path_sans_ext(basename(cor_files)))
    }

    # RNA-Protein specific correlations
    data$rna_protein_cors <- safe_read_csv(file.path(tables_dir, "rna_protein_correlations.csv"))
    data$rna_protein_de_concordance <- safe_read_csv(file.path(tables_dir, "rna_protein_de_concordance.csv"))

    # DE/DA results from individual omics
    data$rna_de_results <- safe_read_csv(file.path(tables_dir, "rna_de_results.csv"))
    data$prot_da_results <- safe_read_csv(file.path(tables_dir, "prot_da_results.csv"))

    # Pathway overlap
    data$pathway_overlap <- safe_read_csv(file.path(tables_dir, "pathway_overlap_jaccard.csv"))
    data$pathway_shared_specific <- safe_read_csv(file.path(tables_dir, "pathway_shared_specific.csv"))

    # Expression matrices for cross-omics visualization
    data$rna_matrix <- safe_read_csv(file.path(tables_dir, "rna_normalized_matrix.csv"))
    data$prot_matrix <- safe_read_csv(file.path(tables_dir, "prot_normalized_matrix.csv"))
    data$metab_matrix <- safe_read_csv(file.path(tables_dir, "metab_normalized_matrix.csv"))

    # Check what data is available
    data$has_mofa <- !is.null(data$mofa_factors)
    data$has_diablo <- !is.null(data$diablo_scores) && length(data$diablo_scores) > 0
    data$has_correlations <- !is.null(data$crossomics_correlations) ||
                              !is.null(data$rna_protein_cors)
    data$has_pathway_overlap <- !is.null(data$pathway_overlap)

    # Determine available omics types
    data$omics_types <- c()
    if (!is.null(data$rna_matrix)) data$omics_types <- c(data$omics_types, "RNA-seq")
    if (!is.null(data$prot_matrix)) data$omics_types <- c(data$omics_types, "Proteomics")
    if (!is.null(data$metab_matrix)) data$omics_types <- c(data$omics_types, "Metabolomics")

    data
}
