#' GMT File Utilities for Pathway Analysis
#'
#' This module provides functions for working with GMT (Gene Matrix Transposed) files,
#' which are the standard format for pathway/gene set definitions.
#'
#' @description
#' Key features:
#' - Read/write GMT files
#' - Validate GMT files and report coverage
#' - Generate GMT files from various sources (biomaRt, KEGG, GO)
#' - Merge and filter GMT files
#' - Support for non-model organisms via custom GMT

# ==============================================================================
# GMT FILE READING
# ==============================================================================

#' Read a GMT file
#'
#' @param gmt_file Path to GMT file
#' @param verbose Print progress messages
#'
#' @return Named list where names are pathway IDs and values are character vectors of gene IDs
#'
#' @details
#' GMT format: each line is tab-separated with:
#' pathway_name <TAB> description <TAB> gene1 <TAB> gene2 <TAB> ...
#'
#' @export
read_gmt <- function(gmt_file, verbose = TRUE) {

    if (!file.exists(gmt_file)) {
        stop(sprintf("GMT file not found: %s", gmt_file))
    }

    if (verbose) message(sprintf("Reading GMT file: %s", gmt_file))

    lines <- readLines(gmt_file, warn = FALSE)
    lines <- lines[lines != ""]  # Remove empty lines

    pathways <- list()
    descriptions <- character()

    for (line in lines) {
        parts <- strsplit(line, "\t")[[1]]

        if (length(parts) < 3) next  # Need at least name, description, 1 gene

        pathway_name <- parts[1]
        pathway_desc <- parts[2]
        genes <- parts[3:length(parts)]
        genes <- genes[genes != ""]  # Remove empty genes

        if (length(genes) > 0) {
            pathways[[pathway_name]] <- genes
            descriptions[pathway_name] <- pathway_desc
        }
    }

    # Add descriptions as attribute
    attr(pathways, "descriptions") <- descriptions

    if (verbose) {
        message(sprintf("  Loaded %d pathways", length(pathways)))
        message(sprintf("  Total unique genes: %d", length(unique(unlist(pathways)))))
        sizes <- sapply(pathways, length)
        message(sprintf("  Pathway sizes: min=%d, median=%d, max=%d",
                        min(sizes), median(sizes), max(sizes)))
    }

    pathways
}

#' Write a GMT file
#'
#' @param pathways Named list of pathways (names = pathway IDs, values = gene vectors)
#' @param output_file Path to output file
#' @param descriptions Optional named vector of pathway descriptions
#' @param verbose Print progress
#'
#' @export
write_gmt <- function(pathways, output_file, descriptions = NULL, verbose = TRUE) {

    if (is.null(descriptions)) {
        descriptions <- attr(pathways, "descriptions")
    }

    lines <- character(length(pathways))

    for (i in seq_along(pathways)) {
        name <- names(pathways)[i]
        genes <- pathways[[i]]
        desc <- descriptions[name] %||% ""

        lines[i] <- paste(c(name, desc, genes), collapse = "\t")
    }

    writeLines(lines, output_file)

    if (verbose) {
        message(sprintf("Wrote %d pathways to: %s", length(pathways), output_file))
    }

    invisible(output_file)
}

# ==============================================================================
# GMT VALIDATION
# ==============================================================================

#' Validate a GMT file against a set of feature IDs
#'
#' @param gmt_file Path to GMT file, or pre-loaded pathway list
#' @param feature_ids Character vector of feature IDs in your dataset
#' @param min_coverage Minimum proportion of genes that should match (0-1)
#' @param min_pathway_size Minimum genes per pathway after filtering
#' @param max_pathway_size Maximum genes per pathway
#' @param verbose Print detailed report
#'
#' @return List with validation results and filtered pathways
#'
#' @export
validate_gmt <- function(gmt_file,
                         feature_ids,
                         min_coverage = 0.1,
                         min_pathway_size = 5,
                         max_pathway_size = 500,
                         verbose = TRUE) {

    # Load GMT if path provided
    if (is.character(gmt_file) && length(gmt_file) == 1) {
        pathways <- read_gmt(gmt_file, verbose = FALSE)
    } else {
        pathways <- gmt_file
    }

    # Get all unique genes in GMT
    all_gmt_genes <- unique(unlist(pathways))

    # Calculate overlap
    matching_genes <- intersect(all_gmt_genes, feature_ids)
    coverage <- length(matching_genes) / length(feature_ids)
    gmt_coverage <- length(matching_genes) / length(all_gmt_genes)

    # Filter pathways to matching genes
    filtered_pathways <- lapply(pathways, function(genes) {
        intersect(genes, feature_ids)
    })

    # Apply size filters
    pathway_sizes <- sapply(filtered_pathways, length)
    keep <- pathway_sizes >= min_pathway_size & pathway_sizes <= max_pathway_size
    filtered_pathways <- filtered_pathways[keep]

    # Preserve descriptions
    if (!is.null(attr(pathways, "descriptions"))) {
        attr(filtered_pathways, "descriptions") <- attr(pathways, "descriptions")[names(filtered_pathways)]
    }

    result <- list(
        valid = coverage >= min_coverage,
        coverage = coverage,
        gmt_coverage = gmt_coverage,
        n_pathways_original = length(pathways),
        n_pathways_filtered = length(filtered_pathways),
        n_genes_in_gmt = length(all_gmt_genes),
        n_genes_matching = length(matching_genes),
        n_features_in_data = length(feature_ids),
        filtered_pathways = filtered_pathways,
        matching_genes = matching_genes
    )

    if (verbose) {
        message("\n=== GMT Validation Report ===")
        message(sprintf("GMT file: %d pathways, %d unique genes",
                        result$n_pathways_original, result$n_genes_in_gmt))
        message(sprintf("Your data: %d features", result$n_features_in_data))
        message(sprintf("\nOverlap:"))
        message(sprintf("  %.1f%% of your features found in GMT (%d genes)",
                        100 * result$coverage, result$n_genes_matching))
        message(sprintf("  %.1f%% of GMT genes found in your data",
                        100 * result$gmt_coverage))
        message(sprintf("\nAfter filtering (size %d-%d):",
                        min_pathway_size, max_pathway_size))
        message(sprintf("  %d pathways retained (from %d original)",
                        result$n_pathways_filtered, result$n_pathways_original))

        if (!result$valid) {
            message(sprintf("\nWARNING: Coverage (%.1f%%) is below threshold (%.1f%%)",
                            100 * result$coverage, 100 * min_coverage))
            message("Consider:")
            message("  1. Using a different GMT file with matching ID types")
            message("  2. Converting your feature IDs to match the GMT (e.g., to gene symbols)")
            message("  3. Creating a custom GMT for your organism")
        } else {
            message(sprintf("\nValidation PASSED: Coverage (%.1f%%) meets threshold",
                            100 * result$coverage))
        }
    }

    result
}

# ==============================================================================
# GMT GENERATION FROM DATABASES
# ==============================================================================

#' Generate GMT from biomaRt GO annotations
#'
#' @param organism Organism name (e.g., "Homo sapiens", "Mus musculus")
#' @param id_type Type of gene IDs to include ("ensembl", "symbol", "entrez")
#' @param ontology GO ontology: "BP" (biological process), "MF" (molecular function),
#'                 "CC" (cellular component), or "all"
#' @param output_file Optional path to save GMT file
#' @param verbose Print progress
#'
#' @return List of pathways (GMT format)
#' @export
generate_gmt_from_go <- function(organism,
                                  id_type = "symbol",
                                  ontology = "BP",
                                  output_file = NULL,
                                  verbose = TRUE) {

    if (!requireNamespace("biomaRt", quietly = TRUE)) {
        stop("biomaRt package required for GMT generation. Install with: BiocManager::install('biomaRt')")
    }

    # Source organism detection utilities
    source_file <- system.file("R", "organism_detection.R", package = "multiomics")
    if (source_file == "") {
        # Try local path
        shared_dir <- file.path(dirname(dirname(getwd())), "shared", "R")
        if (file.exists(file.path(shared_dir, "organism_detection.R"))) {
            source(file.path(shared_dir, "organism_detection.R"))
        }
    }

    # Get Ensembl dataset
    dataset <- get_ensembl_dataset(organism)
    if (is.null(dataset)) {
        stop(sprintf("Could not determine Ensembl dataset for: %s", organism))
    }

    if (verbose) message(sprintf("Connecting to Ensembl (dataset: %s)...", dataset))

    mart <- biomaRt::useMart("ensembl", dataset = dataset)

    # Determine gene attribute to retrieve
    gene_attr <- switch(id_type,
                        "ensembl" = "ensembl_gene_id",
                        "symbol" = "external_gene_name",
                        "entrez" = "entrezgene_id",
                        "external_gene_name")

    # Get GO annotations
    if (verbose) message("Retrieving GO annotations...")

    go_data <- biomaRt::getBM(
        attributes = c(gene_attr, "go_id", "name_1006", "namespace_1003"),
        mart = mart
    )

    # Filter by ontology
    if (ontology != "all") {
        namespace_map <- c(
            "BP" = "biological_process",
            "MF" = "molecular_function",
            "CC" = "cellular_component"
        )
        go_data <- go_data[go_data$namespace_1003 == namespace_map[ontology], ]
    }

    # Remove entries without GO ID or gene
    go_data <- go_data[!is.na(go_data$go_id) & go_data$go_id != "", ]
    go_data <- go_data[!is.na(go_data[[gene_attr]]) & go_data[[gene_attr]] != "", ]

    # Build pathway list
    if (verbose) message("Building pathway list...")

    pathways <- split(go_data[[gene_attr]], go_data$go_id)
    pathways <- lapply(pathways, unique)

    # Create descriptions
    desc_df <- unique(go_data[, c("go_id", "name_1006")])
    descriptions <- setNames(desc_df$name_1006, desc_df$go_id)

    # Filter empty pathways
    pathways <- pathways[sapply(pathways, length) > 0]

    attr(pathways, "descriptions") <- descriptions[names(pathways)]

    if (verbose) {
        message(sprintf("Generated %d GO pathways", length(pathways)))
        sizes <- sapply(pathways, length)
        message(sprintf("Pathway sizes: min=%d, median=%d, max=%d",
                        min(sizes), median(sizes), max(sizes)))
    }

    # Save if output file specified
    if (!is.null(output_file)) {
        write_gmt(pathways, output_file, verbose = verbose)
    }

    pathways
}

#' Generate GMT from KEGG pathways
#'
#' @param organism Organism code (e.g., "hsa" for human, "mmu" for mouse)
#' @param id_type Type of gene IDs: "entrez", "symbol", "ensembl"
#' @param output_file Optional path to save GMT file
#' @param verbose Print progress
#'
#' @return List of pathways (GMT format)
#' @export
generate_gmt_from_kegg <- function(organism,
                                    id_type = "entrez",
                                    output_file = NULL,
                                    verbose = TRUE) {

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("KEGGREST package required. Install with: BiocManager::install('KEGGREST')")
    }

    if (verbose) message(sprintf("Fetching KEGG pathways for organism: %s", organism))

    # Get list of pathways
    pathway_list <- KEGGREST::keggList("pathway", organism)

    if (length(pathway_list) == 0) {
        stop(sprintf("No KEGG pathways found for organism: %s", organism))
    }

    if (verbose) message(sprintf("Found %d pathways, retrieving genes...", length(pathway_list)))

    pathways <- list()
    descriptions <- character()

    # Fetch genes for each pathway (with rate limiting)
    for (i in seq_along(pathway_list)) {
        pathway_id <- names(pathway_list)[i]
        pathway_name <- pathway_list[i]

        # Remove organism prefix from name
        pathway_name <- gsub(paste0("^", organism), "", pathway_name)
        pathway_name <- gsub("^ - ", "", pathway_name)

        tryCatch({
            pathway_info <- KEGGREST::keggGet(pathway_id)[[1]]

            if (!is.null(pathway_info$GENE)) {
                # KEGG returns genes as named vector: "entrez_id" = "symbol; description"
                genes <- names(pathway_info$GENE)

                if (id_type == "symbol") {
                    # Extract symbols from values
                    symbols <- sapply(pathway_info$GENE, function(x) {
                        strsplit(x, ";")[[1]][1]
                    })
                    genes <- unname(symbols)
                }

                if (length(genes) > 0) {
                    # Use clean pathway ID
                    clean_id <- gsub("path:", "", pathway_id)
                    pathways[[clean_id]] <- genes
                    descriptions[clean_id] <- pathway_name
                }
            }
        }, error = function(e) {
            if (verbose) message(sprintf("  Warning: Failed to fetch %s", pathway_id))
        })

        # Rate limiting
        if (i %% 10 == 0) {
            Sys.sleep(0.5)
            if (verbose) message(sprintf("  Progress: %d/%d pathways", i, length(pathway_list)))
        }
    }

    attr(pathways, "descriptions") <- descriptions

    if (verbose) {
        message(sprintf("\nGenerated %d KEGG pathways", length(pathways)))
        if (length(pathways) > 0) {
            sizes <- sapply(pathways, length)
            message(sprintf("Pathway sizes: min=%d, median=%d, max=%d",
                            min(sizes), median(sizes), max(sizes)))
        }
    }

    # Save if output file specified
    if (!is.null(output_file)) {
        write_gmt(pathways, output_file, verbose = verbose)
    }

    pathways
}

# ==============================================================================
# GMT MANIPULATION
# ==============================================================================

#' Merge multiple GMT files
#'
#' @param gmt_files Character vector of GMT file paths, or list of pathway lists
#' @param prefixes Optional prefixes to add to pathway names (avoid collisions)
#' @param verbose Print progress
#'
#' @return Combined list of pathways
#' @export
merge_gmt <- function(gmt_files, prefixes = NULL, verbose = TRUE) {

    all_pathways <- list()
    all_descriptions <- character()

    for (i in seq_along(gmt_files)) {
        # Load GMT if path provided
        if (is.character(gmt_files[[i]]) && length(gmt_files[[i]]) == 1) {
            pathways <- read_gmt(gmt_files[[i]], verbose = FALSE)
        } else {
            pathways <- gmt_files[[i]]
        }

        # Add prefix if specified
        if (!is.null(prefixes) && length(prefixes) >= i) {
            names(pathways) <- paste0(prefixes[i], "_", names(pathways))
            if (!is.null(attr(pathways, "descriptions"))) {
                names(attr(pathways, "descriptions")) <- paste0(prefixes[i], "_",
                                                                 names(attr(pathways, "descriptions")))
            }
        }

        # Merge
        all_pathways <- c(all_pathways, pathways)
        if (!is.null(attr(pathways, "descriptions"))) {
            all_descriptions <- c(all_descriptions, attr(pathways, "descriptions"))
        }
    }

    attr(all_pathways, "descriptions") <- all_descriptions

    if (verbose) {
        message(sprintf("Merged %d GMT sources into %d pathways",
                        length(gmt_files), length(all_pathways)))
    }

    all_pathways
}

#' Filter GMT pathways by size
#'
#' @param pathways List of pathways (GMT format)
#' @param min_size Minimum number of genes
#' @param max_size Maximum number of genes
#'
#' @return Filtered list of pathways
#' @export
filter_gmt_by_size <- function(pathways, min_size = 10, max_size = 500) {
    sizes <- sapply(pathways, length)
    keep <- sizes >= min_size & sizes <= max_size

    filtered <- pathways[keep]

    # Preserve descriptions
    if (!is.null(attr(pathways, "descriptions"))) {
        attr(filtered, "descriptions") <- attr(pathways, "descriptions")[names(filtered)]
    }

    filtered
}

#' Convert GMT gene IDs using a mapping table
#'
#' @param pathways List of pathways (GMT format)
#' @param mapping Data frame with columns: from, to
#' @param verbose Print progress
#'
#' @return Pathways with converted gene IDs
#' @export
convert_gmt_ids <- function(pathways, mapping, verbose = TRUE) {

    # Create lookup
    id_map <- setNames(mapping$to, mapping$from)

    converted <- lapply(pathways, function(genes) {
        new_genes <- id_map[genes]
        new_genes <- new_genes[!is.na(new_genes)]
        unique(new_genes)
    })

    # Remove empty pathways
    converted <- converted[sapply(converted, length) > 0]

    # Preserve descriptions
    if (!is.null(attr(pathways, "descriptions"))) {
        attr(converted, "descriptions") <- attr(pathways, "descriptions")[names(converted)]
    }

    if (verbose) {
        message(sprintf("Converted GMT: %d -> %d pathways (removed %d empty)",
                        length(pathways), length(converted),
                        length(pathways) - length(converted)))
    }

    converted
}

# ==============================================================================
# GMT STATISTICS
# ==============================================================================

#' Get summary statistics for a GMT file
#'
#' @param pathways List of pathways (GMT format)
#' @return Data frame with pathway statistics
#' @export
gmt_summary <- function(pathways) {

    sizes <- sapply(pathways, length)
    descriptions <- attr(pathways, "descriptions")

    data.frame(
        pathway_id = names(pathways),
        description = if (!is.null(descriptions)) descriptions[names(pathways)] else NA_character_,
        n_genes = sizes,
        stringsAsFactors = FALSE
    )
}

# ==============================================================================
# HELPER
# ==============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

# Reference external function if needed
get_ensembl_dataset <- function(organism) {
    # Inline version for when annotation_utils.R not loaded
    organism <- tolower(gsub("\\s+", "_", organism))

    datasets <- list(
        "homo_sapiens" = "hsapiens_gene_ensembl",
        "human" = "hsapiens_gene_ensembl",
        "mus_musculus" = "mmusculus_gene_ensembl",
        "mouse" = "mmusculus_gene_ensembl",
        "rattus_norvegicus" = "rnorvegicus_gene_ensembl",
        "rat" = "rnorvegicus_gene_ensembl",
        "danio_rerio" = "drerio_gene_ensembl",
        "zebrafish" = "drerio_gene_ensembl",
        "drosophila_melanogaster" = "dmelanogaster_gene_ensembl",
        "fly" = "dmelanogaster_gene_ensembl",
        "caenorhabditis_elegans" = "celegans_gene_ensembl",
        "worm" = "celegans_gene_ensembl",
        "saccharomyces_cerevisiae" = "scerevisiae_gene_ensembl",
        "yeast" = "scerevisiae_gene_ensembl",
        "arabidopsis_thaliana" = "athaliana_gene_ensembl",
        "gallus_gallus" = "ggallus_gene_ensembl",
        "chicken" = "ggallus_gene_ensembl",
        "sus_scrofa" = "sscrofa_gene_ensembl",
        "pig" = "sscrofa_gene_ensembl",
        "bos_taurus" = "btaurus_gene_ensembl",
        "cow" = "btaurus_gene_ensembl"
    )

    datasets[[organism]]
}
