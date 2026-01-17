# =============================================================================
# Metabolite Mapping to Pathways (KEGG / SMPDB)
# =============================================================================

#' Main entry point for metabolite mapping
#'
#' @param ingested_data Data from ingest_data target
#' @param config Pipeline configuration
#' @return List with pathway mapping and summary
map_metabolites_to_pathways <- function(ingested_data, config) {
    log_message("=== Starting Metabolite Mapping ===")

    if (!isTRUE(config$enrichment$mapping$enabled)) {
        log_message("Mapping disabled in config. Skipping.")
        return(NULL)
    }

    # If manual mapping already exists, we might want to skip or merge.
    # For now, let's prioritize manual mapping if present in ingested_data.
    if (!is.null(ingested_data$annotations$pathway_mapping)) {
        log_message("Manual pathway mapping already exists. Skipping automatic mapping.")
        return(ingested_data$annotations$pathway_mapping)
    }

    database <- config$enrichment$mapping$database %||% "kegg"
    organism <- config$enrichment$mapping$organism %||% "hsa" # Default to human

    # 1. Gather names for mapping
    # Priority:
    # 1. ingested_data$annotations$annotation_table$metabolite_name
    # 2. ingested_data$feature_metadata$preliminary_id (or similar columns)
    # 3. feature_id (rownames of matrix)

    feature_ids <- rownames(ingested_data$matrix)
    mapping_input <- data.frame(
        feature_id = feature_ids,
        query_name = NA_character_,
        stringsAsFactors = FALSE
    )

    # Check annotations
    ann <- ingested_data$annotations$annotation_table
    if (!is.null(ann)) {
        name_col <- intersect(colnames(ann), c("metabolite_name", "compound_name", "name"))
        if (length(name_col) > 0) {
            mapping_input$query_name <- ann[[name_col[1]]][match(mapping_input$feature_id, ann$feature_id)]
        }
    }

    # Check feature_metadata if query_name is still NA
    feat_meta <- ingested_data$feature_metadata
    if (!is.null(feat_meta)) {
        name_col <- intersect(colnames(feat_meta), c("preliminary_id", "metabolite_name", "compound_name", "name"))
        if (length(name_col) > 0) {
            na_idx <- is.na(mapping_input$query_name) | mapping_input$query_name == "Unknown"
            mapping_input$query_name[na_idx] <- feat_meta[[name_col[1]]][match(mapping_input$feature_id[na_idx], feat_meta$feature_id)]
        }
    }

    # Fallback to feature_id for remaining NAs
    na_idx <- is.na(mapping_input$query_name) | mapping_input$query_name == "Unknown"
    mapping_input$query_name[na_idx] <- mapping_input$feature_id[na_idx]

    # Remove NAs and "Unknown"
    names_to_query <- unique(mapping_input$query_name[!is.na(mapping_input$query_name) & mapping_input$query_name != "Unknown" & mapping_input$query_name != ""])

    if (length(names_to_query) == 0) {
        log_message("No valid metabolite names found for mapping.")
        return(NULL)
    }

    pathway_mapping <- NULL

    if (database == "kegg") {
        pathway_mapping <- map_to_kegg(names_to_query, organism)
    } else if (database == "smpdb") {
        # SMPDB mapping often goes through HMDB IDs
        hmdb_ids <- NULL
        if (!is.null(feat_meta) && "hmdb_id" %in% colnames(feat_meta)) {
            hmdb_ids <- feat_meta$hmdb_id[match(mapping_input$feature_id, feat_meta$feature_id)]
        }
        pathway_mapping <- map_to_smpdb(names_to_query, hmdb_ids)
    }

    if (is.null(pathway_mapping) || nrow(pathway_mapping) == 0) {
        log_message("Mapping failed to produce results.")
        return(NULL)
    }

    # Merge back to feature_ids
    final_mapping <- merge(pathway_mapping, mapping_input, by = "query_name", all.x = TRUE)

    # Ensure columns expected by run_pathway_enrichment are present
    if (database == "kegg") {
        final_mapping$kegg_id <- final_mapping$compound_id
    }

    log_message(
        "Mapped ", length(unique(final_mapping$feature_id)), " features to ",
        length(unique(final_mapping$pathway_id)), " pathways."
    )

    # Save summary
    summary_df <- data.frame(
        database = database,
        total_queried = length(names_to_query),
        mapped_names = length(unique(final_mapping$query_name)),
        mapped_features = length(unique(final_mapping$feature_id)),
        unique_pathways = length(unique(final_mapping$pathway_id))
    )
    save_table(summary_df, "mapping_summary.csv", config, "qc")

    return(final_mapping)
}

#' Map names to KEGG pathways
#'
#' @param names Character vector of names
#' @param organism KEGG organism code
#' @return Pathway mapping data frame
map_to_kegg <- function(names, organism) {
    log_message("Querying KEGG for ", length(names), " names...")

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        log_message("KEGGREST package not found. Skipping KEGG mapping.")
        return(NULL)
    }

    mapping_list <- list()

    # Limit to 50 for trial run if needed, but here we'll try all (there are only 20 in example)
    if (length(names) > 100) {
        log_message("Limiting KEGG query to first 100 metabolites for performance.")
        names <- names[1:100]
    }

    for (nm in names) {
        try(
            {
                # Search for compound
                res <- suppressMessages(KEGGREST::keggFind("compound", nm))

                if (length(res) > 0) {
                    # Try to find an exact match in the names
                    # KEGG names often contain multiple aliases separated by ;
                    cpd_id <- names(res)[1]

                    # Check if any alias matches exactly
                    all_aliases <- tolower(unlist(strsplit(res, "; ")))
                    if (tolower(nm) %in% all_aliases) {
                        # Found an exact alias, but which compound ID does it belong to?
                        # Actually, keggFind returns a vector where names are IDs.
                        # Let's iterate over hits and see which one has the best match.
                        for (i in 1:length(res)) {
                            aliases <- tolower(unlist(strsplit(res[i], "; ")))
                            if (tolower(nm) %in% aliases) {
                                cpd_id <- names(res)[i]
                                break
                            }
                        }
                    }

                    cpd_id <- gsub("cpd:", "", cpd_id)

                    # Get Pathways for this compound
                    pws <- suppressMessages(KEGGREST::keggLink("pathway", paste0("cpd:", cpd_id)))

                    if (length(pws) > 0) {
                        # Map mapXXXXX to organismXXXXX
                        path_ids <- gsub("path:map", paste0(organism), pws)
                        # Keep only organism pathways
                        path_ids <- path_ids[grepl(paste0("^", organism), path_ids)]

                        if (length(path_ids) > 0) {
                            mapping_list[[nm]] <- data.frame(
                                query_name = nm,
                                compound_id = cpd_id,
                                pathway_id = path_ids,
                                pathway_name = path_ids, # Placeholder
                                stringsAsFactors = FALSE
                            )
                        }
                    }
                }
            },
            silent = TRUE
        )
    }

    if (length(mapping_list) == 0) {
        return(NULL)
    }

    df <- do.call(rbind, mapping_list)
    rownames(df) <- NULL

    # Retrieval of pathway names in bulk
    unique_pws <- unique(df$pathway_id)
    if (length(unique_pws) > 0) {
        try(
            {
                all_pws <- suppressMessages(KEGGREST::keggList("pathway", organism))
                names(all_pws) <- gsub("path:", "", names(all_pws))

                idx <- match(df$pathway_id, names(all_pws))
                valid_idx <- !is.na(idx)

                df$pathway_name[valid_idx] <- as.character(all_pws[idx[valid_idx]])
                df <- df[valid_idx, ] # Keep only those found in the organism
            },
            silent = TRUE
        )
    }

    return(df)
}

#' Map names to SMPDB pathways
#'
#' @param names Character vector of names
#' @param hmdb_ids Optional HMDB IDs
#' @return Pathway mapping data frame
map_to_smpdb <- function(names, hmdb_ids = NULL) {
    log_message("SMPDB mapping using HMDB IDs...")

    if (is.null(hmdb_ids) || length(hmdb_ids) == 0) {
        log_message("No HMDB IDs provided for SMPDB mapping.")
        return(NULL)
    }

    # MetaboAnalystR is usually the best way for SMPDB mapping if available
    # For now, we'll provide a simplified placeholder or instruction.
    log_message("SMPDB mapping is currently limited. Please use KEGG or provide a manual mapping.")
    return(NULL)
}
