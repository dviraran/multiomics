# =============================================================================
# Identifier Harmonization / Mapping Layer
# =============================================================================
#
# This file uses shared utilities from ../shared/R/ when available for:
# - Unified annotation with fallback chain
# - Organism and ID type detection
# - Support for non-model organisms

# Source shared utilities if available
.source_shared_utils <- function() {
  possible_paths <- c(
    file.path(dirname(dirname(getwd())), "shared", "R"),
    file.path(dirname(getwd()), "shared", "R"),
    file.path(getwd(), "..", "shared", "R"),
    file.path(getwd(), "..", "..", "shared", "R")
  )

  for (shared_dir in possible_paths) {
    if (dir.exists(shared_dir)) {
      for (util_file in c("annotation_utils.R", "organism_detection.R", "gmt_utils.R")) {
        full_path <- file.path(shared_dir, util_file)
        if (file.exists(full_path)) {
          source(full_path)
        }
      }
      return(TRUE)
    }
  }
  return(FALSE)
}

.source_shared_utils()

#' Generate centralized ID mapping file (Gene <-> Protein)
#'
#' @param config Configuration list
#' @param output_dir Directory to save the mapping file
#' @return Data frame with columns: gene_id, entrez_id, uniprot_id, gene_symbol
generate_id_mapping <- function(config, output_dir = "outputs/tables") {
  log_message("=== Generating Centralized ID Mapping ===")

  # Check if mapping file is provided in config (e.g., transcriptomics or proteomics mapping)
  # The user request implies checking global or specific configs.
  # If a specific file is provided in config$global$mapping_file (hypothetically) or we use the raw mapping files.
  # But here we focus on creating one if missing.

  # Check if we should use an existing file (if specified in a custom location)
  if (!is.null(config$global$gene_protein_mapping) && file.exists(config$global$gene_protein_mapping)) {
    log_message("Using user-provided mapping file: ", config$global$gene_protein_mapping)
    mapping_df <- read.csv(config$global$gene_protein_mapping, stringsAsFactors = FALSE)

    # Basic validation
    required <- c("gene_id", "entrez_id", "uniprot_id", "gene_symbol")
    if (!all(required %in% colnames(mapping_df))) {
      log_message("WARNING: Provided mapping file missing required columns: ", paste(setdiff(required, colnames(mapping_df)), collapse = ", "))
    }

    return(mapping_df)
  }

  organism <- config$global$organism %||% "human"
  mapping_df <- NULL

  log_message("Organism: ", organism)

  if (organism == "c_elegans") {
    if (requireNamespace("org.Ce.eg.db", quietly = TRUE)) {
      db <- org.Ce.eg.db::org.Ce.eg.db
      log_message("Using org.Ce.eg.db for mapping generation")

      # Get all WBGene IDs
      keys <- keys(db, keytype = "WORMBASE")

      tryCatch(
        {
          # Map WORMBASE -> ENTREZID
          entrez <- AnnotationDbi::mapIds(db, keys = keys, column = "ENTREZID", keytype = "WORMBASE", multiVals = "first")

          # Map WORMBASE -> SYMBOL
          symbol <- AnnotationDbi::mapIds(db, keys = keys, column = "SYMBOL", keytype = "WORMBASE", multiVals = "first")

          # Map WORMBASE -> UNIPROT (via Entrez or directly if supported)
          # org.Ce.eg.db often maps Entrez -> Uniprot.
          # Let's map Entrez -> Uniprot.

          # We need a dataframe of keys, entrez, symbol
          df <- data.frame(
            gene_id = keys,
            entrez_id = entrez,
            gene_symbol = symbol,
            stringsAsFactors = FALSE
          )

          # Filter for those with Entrez IDs to map to Uniprot
          valid_entrez <- df[!is.na(df$entrez_id), ]

          # Map ENTREZID -> UNIPROT
          # Note: input must be character
          uniprot <- AnnotationDbi::mapIds(db, keys = as.character(valid_entrez$entrez_id), column = "UNIPROT", keytype = "ENTREZID", multiVals = "first")

          # Add uniprot to df
          df$uniprot_id <- NA
          df$uniprot_id[!is.na(df$entrez_id)] <- uniprot

          mapping_df <- df
          log_message("Generated mapping for ", nrow(mapping_df), " genes")
          log_message("  - With Entrez ID: ", sum(!is.na(mapping_df$entrez_id)))
          log_message("  - With UniProt ID: ", sum(!is.na(mapping_df$uniprot_id)))
        },
        error = function(e) {
          log_message("Error generating mapping: ", e$message)
        }
      )
    }
  } else if (organism == "human") {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      db <- org.Hs.eg.db::org.Hs.eg.db
      # Similar logic for human (typically Ensembl -> Entrez -> Uniprot)
      keys <- keys(db, keytype = "ENSEMBL")
      entrez <- AnnotationDbi::mapIds(db, keys = keys, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
      symbol <- AnnotationDbi::mapIds(db, keys = keys, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

      df <- data.frame(gene_id = keys, entrez_id = entrez, gene_symbol = symbol, stringsAsFactors = FALSE)
      valid_entrez <- df[!is.na(df$entrez_id), ]
      uniprot <- AnnotationDbi::mapIds(db, keys = as.character(valid_entrez$entrez_id), column = "UNIPROT", keytype = "ENTREZID", multiVals = "first")
      df$uniprot_id <- NA
      df$uniprot_id[!is.na(df$entrez_id)] <- uniprot
      mapping_df <- df
    }
  }

  if (!is.null(mapping_df)) {
    # Ensure output dir exists
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out_file <- file.path(output_dir, "gene_protein_mapping.csv")
    write.csv(mapping_df, out_file, row.names = FALSE)
    log_message("Mapping saved to ", out_file)
    return(mapping_df)
  } else {
    log_message("Failed to generate mapping or organism not supported")
    return(NULL)
  }
}

#' Harmonize identifiers across all omics
#' @param gene_protein_mapping Centralized mapping table (optional)
harmonize_identifiers <- function(preprocessed_data, config, gene_protein_mapping = NULL) {
  log_message("=== Harmonizing Identifiers Across Omics ===")

  processed_omics <- preprocessed_data$processed_omics
  harmonized <- list()

  # Transcriptomics: Ensembl -> Gene Symbol
  if ("transcriptomics" %in% names(processed_omics)) {
    harmonized$transcriptomics <- harmonize_transcriptomics_ids(
      processed_omics$transcriptomics, config, gene_protein_mapping
    )
  }

  # Proteomics: UniProt/Accession -> Gene Symbol
  if ("proteomics" %in% names(processed_omics)) {
    harmonized$proteomics <- harmonize_proteomics_ids(
      processed_omics$proteomics, config, gene_protein_mapping
    )
  }

  # Metabolomics: Feature ID -> KEGG/HMDB/ChEBI
  if ("metabolomics" %in% names(processed_omics)) {
    harmonized$metabolomics <- harmonize_metabolomics_ids(
      processed_omics$metabolomics, config
    )
  }

  log_message("=== Identifier Harmonization Complete ===")

  list(
    harmonized_omics = harmonized,
    metadata = preprocessed_data$metadata,
    alignment = preprocessed_data$alignment
  )
}

#' Harmonize transcriptomics identifiers
harmonize_transcriptomics_ids <- function(rna_data, config, gene_protein_mapping = NULL) {
  log_message("Harmonizing transcriptomics identifiers...")

  mat <- rna_data$normalized_matrix
  mapping <- rna_data$mapping

  # Strip Ensembl version suffix (ENSG00000123456.1 -> ENSG00000123456)
  original_ids <- rownames(mat)
  stripped_ids <- strip_ensembl_version(original_ids)
  rownames(mat) <- stripped_ids

  # Create feature annotation
  feature_anno <- data.frame(
    feature_id = stripped_ids,
    original_id = original_ids,
    stringsAsFactors = FALSE
  )

  feature_anno$gene_symbol <- NA
  feature_anno$entrez_id <- NA
  feature_anno$uniprot_id <- NA

  # Apply mapping if provided
  if (!is.null(mapping) && nrow(mapping) > 0) {
    # ... (existing manual mapping logic) ...
    # Keep existing logic as override? Or prioritize centralized?
    # User said: "if the config already provide a mapping file use it" -> This refers to the file loaded into `mapping`.
    # So kept existing logic for `mapping`.

    # (Omitted changes to existing block for brevity, assuming standard priority)
    # Just pass through or careful not to break.
    # Wait, I need to insert the centralized mapping check AFTER or AS FALLBACK to the specific mapping.
  }

  # 1. Apply specific mapping if provided in config logic (loaded in ingestion)
  if (!is.null(mapping) && nrow(mapping) > 0) {
    mapping_cols <- colnames(mapping)

    ensembl_col <- grep("ensembl|gene_id", mapping_cols, ignore.case = TRUE, value = TRUE)[1]
    symbol_col <- grep("symbol|gene_name|hgnc", mapping_cols, ignore.case = TRUE, value = TRUE)[1]
    entrez_col <- grep("entrez|ncbi", mapping_cols, ignore.case = TRUE, value = TRUE)[1]

    if (!is.na(ensembl_col) && !is.na(symbol_col)) {
      # Strip versions from mapping too
      mapping[[ensembl_col]] <- strip_ensembl_version(mapping[[ensembl_col]])

      # Match
      idx <- match(stripped_ids, mapping[[ensembl_col]])

      matches <- !is.na(idx)
      if (sum(matches) > 0) {
        feature_anno$gene_symbol[matches] <- mapping[[symbol_col]][idx[matches]]
        if (!is.na(entrez_col)) {
          feature_anno$entrez_id[matches] <- mapping[[entrez_col]][idx[matches]]
        }
        n_mapped <- sum(!is.na(feature_anno$gene_symbol))
        log_message("Mapped ", n_mapped, "/", nrow(feature_anno), " genes to symbols via provided mapping")
      }
    }
  }

  # 2. Centralized Mapping Check (fill gaps)
  if (!is.null(gene_protein_mapping)) {
    # Helper to map if symbol missing
    missing <- is.na(feature_anno$gene_symbol)
    if (any(missing)) {
      log_message("Using centralized ID mapping for remaining ", sum(missing), " transcriptomics features")

      # Match stripped_ids (WBGene/Ensembl) to gene_id
      idx <- match(stripped_ids[missing], gene_protein_mapping$gene_id)

      matches_sub <- !is.na(idx)
      if (sum(matches_sub) > 0) {
        # Care with indexing: feature_anno[missing, ][matches_sub, ]
        # Simplest: iterate or map

        # Get mapped values
        mapped_symbol <- gene_protein_mapping$gene_symbol[idx[matches_sub]]
        mapped_entrez <- gene_protein_mapping$entrez_id[idx[matches_sub]]
        mapped_uniprot <- gene_protein_mapping$uniprot_id[idx[matches_sub]]

        # Assign back
        # Indices in full df where missing & matched
        which_missing <- which(missing)
        target_indices <- which_missing[matches_sub]

        feature_anno$gene_symbol[target_indices] <- mapped_symbol
        feature_anno$entrez_id[target_indices] <- mapped_entrez
        feature_anno$uniprot_id[target_indices] <- mapped_uniprot

        log_message("Mapped additional ", length(target_indices), " genes via centralized mapping")
      }
    }
  }

  # 3. Fallback: Heuristics based on organism and ID type
  organism <- config$global$organism %||% "human"

  if (all(is.na(feature_anno$gene_symbol))) {
    # If practically nothing mapped, check if IDs are themselves symbols
    if (!any(grepl("^ENS", stripped_ids)) && !any(grepl("^WBGene", stripped_ids))) {
      feature_anno$gene_symbol <- stripped_ids
      log_message("IDs appear to be gene symbols, using as-is")
    } else if (any(grepl("^WBGene", stripped_ids)) && organism == "c_elegans") {
      # WBGene IDs for C. elegans - try direct org.Ce.eg.db query
      log_message("WARNING: IDs are WBGene accessions but mapping failed. Attempting direct org.Ce.eg.db query...")

      if (requireNamespace("org.Ce.eg.db", quietly = TRUE)) {
        tryCatch(
          {
            db <- org.Ce.eg.db::org.Ce.eg.db
            mapped_symbols <- AnnotationDbi::mapIds(
              db,
              keys = stripped_ids,
              column = "SYMBOL",
              keytype = "WORMBASE",
              multiVals = "first"
            )
            feature_anno$gene_symbol <- mapped_symbols
            n_mapped <- sum(!is.na(feature_anno$gene_symbol))
            log_message("  Mapped ", n_mapped, "/", length(stripped_ids), " features via direct org.Ce.eg.db query")
          },
          error = function(e) {
            log_message("  ERROR: Direct org.Ce.eg.db query failed: ", e$message)
          }
        )
      }

      # If still unmapped, leave as NA and warn
      n_unmapped <- sum(is.na(feature_anno$gene_symbol))
      if (n_unmapped > 0) {
        log_message("  WARNING: ", n_unmapped, " WBGene IDs remain unmapped. Enrichment analysis will be limited.")
        log_message("  Consider providing a mapping file in config or checking org.Ce.eg.db installation.")
      }
    } else if (any(grepl("^ENS", stripped_ids))) {
      # Ensembl IDs - try direct query based on organism
      log_message("WARNING: IDs are Ensembl accessions but mapping failed. Attempting direct query...")

      org_pkg <- NULL
      keytype <- "ENSEMBL"
      if (organism == "human") {
        org_pkg <- "org.Hs.eg.db"
      } else if (organism == "mouse") {
        org_pkg <- "org.Mm.eg.db"
      } else if (organism == "c_elegans") {
        # Ensembl IDs for C. elegans are rare but possible
        org_pkg <- "org.Ce.eg.db"
        keytype <- "ENSEMBL"
      }

      if (!is.null(org_pkg) && requireNamespace(org_pkg, quietly = TRUE)) {
        tryCatch(
          {
            db <- get(org_pkg, envir = asNamespace(org_pkg))
            mapped_symbols <- AnnotationDbi::mapIds(
              db,
              keys = stripped_ids,
              column = "SYMBOL",
              keytype = keytype,
              multiVals = "first"
            )
            feature_anno$gene_symbol <- mapped_symbols
            n_mapped <- sum(!is.na(feature_anno$gene_symbol))
            log_message("  Mapped ", n_mapped, "/", length(stripped_ids), " features via direct ", org_pkg, " query")
          },
          error = function(e) {
            log_message("  ERROR: Direct ", org_pkg, " query failed: ", e$message)
          }
        )
      }

      n_unmapped <- sum(is.na(feature_anno$gene_symbol))
      if (n_unmapped > 0) {
        log_message("  WARNING: ", n_unmapped, " Ensembl IDs remain unmapped. Enrichment analysis will be limited.")
      }
    }
  } else {
    # Partial mapping handling - try to fill remaining gaps
    missing <- is.na(feature_anno$gene_symbol)
    if (any(missing)) {
      n_missing <- sum(missing)
      log_message("  ", n_missing, " features still unmapped after primary mapping")

      missing_ids <- stripped_ids[missing]

      # For WBGene IDs in C. elegans, attempt direct query
      if (organism == "c_elegans" && any(grepl("^WBGene", missing_ids))) {
        wb_idx <- grepl("^WBGene", missing_ids)
        if (requireNamespace("org.Ce.eg.db", quietly = TRUE)) {
          tryCatch(
            {
              db <- org.Ce.eg.db::org.Ce.eg.db
              mapped_symbols <- AnnotationDbi::mapIds(
                db,
                keys = missing_ids[wb_idx],
                column = "SYMBOL",
                keytype = "WORMBASE",
                multiVals = "first"
              )
              # Update feature_anno at correct indices
              missing_indices <- which(missing)
              feature_anno$gene_symbol[missing_indices[wb_idx]] <- mapped_symbols
              n_wb_mapped <- sum(!is.na(mapped_symbols))
              log_message("  Mapped ", n_wb_mapped, " additional WBGene IDs via org.Ce.eg.db")
            },
            error = function(e) {
              log_message("  ERROR: org.Ce.eg.db query failed: ", e$message)
            }
          )
        }
      }

      # Re-check what's still missing
      still_missing <- is.na(feature_anno$gene_symbol)
      if (any(still_missing)) {
        missing_ids_final <- stripped_ids[still_missing]
        are_accessions <- grepl("^(ENS|WBGene)", missing_ids_final)

        # Only use IDs as symbols if they don't look like database accessions
        if (any(!are_accessions)) {
          feature_anno$gene_symbol[still_missing][!are_accessions] <- missing_ids_final[!are_accessions]
          log_message("  Used ", sum(!are_accessions), " non-accession IDs as gene symbols")
        }

        # For remaining accessions, leave as NA and warn
        if (any(are_accessions)) {
          log_message("  WARNING: ", sum(are_accessions), " database accessions remain unmapped (left as NA)")
          log_message("  These features will be excluded from enrichment analysis")
        }
      }
    }
  }

  # Update DE table if present
  de_table <- rna_data$de_table
  if (!is.null(de_table)) {
    de_table$gene_id <- strip_ensembl_version(de_table$gene_id)
    if ("gene_symbol" %in% colnames(feature_anno)) {
      idx <- match(de_table$gene_id, feature_anno$feature_id)
      de_table$gene_symbol <- feature_anno$gene_symbol[idx]
    }
  }

  list(
    normalized_matrix = mat,
    de_table = de_table,
    feature_annotation = feature_anno,
    gmt = rna_data$gmt
  )
}

#' Harmonize proteomics identifiers
harmonize_proteomics_ids <- function(prot_data, config, gene_protein_mapping = NULL) {
  log_message("Harmonizing proteomics identifiers...")

  mat <- prot_data$normalized_matrix
  mapping <- prot_data$mapping

  # Parse protein accessions
  original_ids <- rownames(mat)
  parsed <- parse_protein_ids(original_ids)

  # Create feature annotation
  feature_anno <- data.frame(
    feature_id = parsed$clean_id,
    original_id = original_ids,
    accession_type = parsed$type,
    stringsAsFactors = FALSE
  )

  # Update matrix rownames
  rownames(mat) <- parsed$clean_id

  # Apply mapping if provided
  if (!is.null(mapping) && nrow(mapping) > 0) {
    mapping_cols <- colnames(mapping)

    # Find matching columns
    accession_col <- grep("accession|protein|uniprot", mapping_cols, ignore.case = TRUE, value = TRUE)[1]
    symbol_col <- grep("symbol|gene", mapping_cols, ignore.case = TRUE, value = TRUE)[1]

    if (!is.na(accession_col) && !is.na(symbol_col)) {
      # Try matching on parsed IDs
      idx <- match(parsed$clean_id, mapping[[accession_col]])

      # If poor matching, try original IDs
      if (sum(!is.na(idx)) < length(idx) * 0.1) {
        idx <- match(original_ids, mapping[[accession_col]])
      }

      feature_anno$gene_symbol <- mapping[[symbol_col]][idx]
      n_mapped <- sum(!is.na(feature_anno$gene_symbol))
      log_message("Mapped ", n_mapped, "/", nrow(feature_anno), " proteins to gene symbols via provided mapping")
    }
  }

  # Centralized Mapping Check
  if (!is.null(gene_protein_mapping)) {
    log_message("Using centralized ID mapping for Proteomics")
    # Match clean_id (UniProt) to uniprot_id
    idx <- match(parsed$clean_id, gene_protein_mapping$uniprot_id)

    matches <- !is.na(idx)
    if (sum(matches) > 0) {
      feature_anno$gene_symbol[matches] <- gene_protein_mapping$gene_symbol[idx[matches]]
      log_message("Mapped ", sum(matches), "/", nrow(feature_anno), " proteins to symbols via centralized mapping")
    }
  } else {
    # Fallback: Try to map using OrgDb if available
    organism <- config$global$organism %||% "human"
    org_pkg <- NULL

    if (organism == "c_elegans") {
      org_pkg <- "org.Ce.eg.db"
    } else if (organism == "human") {
      org_pkg <- "org.Hs.eg.db"
    }

    if (!is.null(org_pkg) && requireNamespace(org_pkg, quietly = TRUE)) {
      log_message("Attempting to map UniProt IDs using ", org_pkg)
      db <- get(org_pkg, envir = asNamespace(org_pkg))

      tryCatch(
        {
          # Try mapping UNIPROT -> SYMBOL
          mapped_symbols <- AnnotationDbi::mapIds(
            db,
            keys = parsed$clean_id,
            column = "SYMBOL",
            keytype = "UNIPROT",
            multiVals = "first"
          )

          feature_anno$gene_symbol <- mapped_symbols
          n_mapped <- sum(!is.na(feature_anno$gene_symbol))
          log_message("Mapped ", n_mapped, "/", nrow(feature_anno), " proteins to gene symbols via ", org_pkg)

          # If mapping is poor, try WORMBASE for C. elegans if ID looks like WormBase ID?
          # But here IDs are Uniprot.
        },
        error = function(e) {
          log_message("Database mapping failed: ", e$message)
        }
      )
    }
  }

  # Update DA table if present
  da_table <- prot_data$da_table
  if (!is.null(da_table)) {
    # Update feature IDs
    old_ids <- da_table$feature_id
    parsed_da <- parse_protein_ids(old_ids)
    da_table$feature_id <- parsed_da$clean_id

    if ("gene_symbol" %in% colnames(feature_anno)) {
      idx <- match(da_table$feature_id, feature_anno$feature_id)
      da_table$gene_symbol <- feature_anno$gene_symbol[idx]
    }
  }

  list(
    normalized_matrix = mat,
    da_table = da_table,
    feature_annotation = feature_anno
  )
}

#' Parse protein identifiers
parse_protein_ids <- function(ids) {
  clean_ids <- character(length(ids))
  types <- character(length(ids))

  for (i in seq_along(ids)) {
    id <- ids[i]

    # UniProt format: sp|P12345|GENE_HUMAN or tr|Q12345|...
    if (grepl("^(sp|tr)\\|", id)) {
      parts <- strsplit(id, "\\|")[[1]]
      clean_ids[i] <- parts[2]
      types[i] <- ifelse(parts[1] == "sp", "UniProt_SwissProt", "UniProt_TrEMBL")
    }
    # MaxQuant REV__ or CON__ prefixes
    else if (grepl("^(REV__|CON__)", id)) {
      clean_ids[i] <- sub("^(REV__|CON__)", "", id)
      types[i] <- "contaminant_reverse"
    }
    # UniProt accession (standard or extended)
    # Matches P12345 (6 char) or A0A061AKV1 (10 char)
    else if (grepl("^[A-Z][0-9][A-Z0-9]+", id)) {
      # Remove isoform suffix
      clean_ids[i] <- sub("-[0-9]+$", "", id)
      types[i] <- "UniProt"
    }
    # Ensembl protein
    else if (grepl("^ENSP", id)) {
      clean_ids[i] <- strip_ensembl_version(id)
      types[i] <- "Ensembl_protein"
    }
    # RefSeq protein
    else if (grepl("^[NXY]P_", id)) {
      clean_ids[i] <- sub("\\.[0-9]+$", "", id)
      types[i] <- "RefSeq"
    }
    # Default: use as-is
    else {
      clean_ids[i] <- id
      types[i] <- "unknown"
    }
  }

  list(clean_id = clean_ids, type = types)
}

#' Harmonize metabolomics identifiers
harmonize_metabolomics_ids <- function(metab_data, config) {
  log_message("Harmonizing metabolomics identifiers...")

  mat <- metab_data$normalized_matrix
  annotation <- metab_data$annotation
  feature_metadata <- metab_data$feature_metadata

  # Create feature annotation from feature IDs
  feature_ids <- rownames(mat)
  feature_anno <- data.frame(
    feature_id = feature_ids,
    stringsAsFactors = FALSE
  )

  # Add feature metadata if available (mz, rt, adduct)
  if (!is.null(feature_metadata) && nrow(feature_metadata) > 0) {
    meta_cols <- colnames(feature_metadata)
    id_col <- meta_cols[1] # Assume first column is ID

    idx <- match(feature_ids, feature_metadata[[id_col]])

    # Add mz
    mz_col <- grep("^mz$|^m/z$|^mass$", meta_cols, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(mz_col)) {
      feature_anno$mz <- feature_metadata[[mz_col]][idx]
    }

    # Add rt
    rt_col <- grep("^rt$|^retention", meta_cols, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(rt_col)) {
      feature_anno$rt <- feature_metadata[[rt_col]][idx]
    }

    # Add adduct
    adduct_col <- grep("adduct", meta_cols, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(adduct_col)) {
      feature_anno$adduct <- feature_metadata[[adduct_col]][idx]
    }
  }

  # Add annotation if available (KEGG, HMDB, ChEBI, compound name)
  if (!is.null(annotation) && nrow(annotation) > 0) {
    anno_cols <- colnames(annotation)
    id_col <- anno_cols[1]

    idx <- match(feature_ids, annotation[[id_col]])

    # KEGG
    kegg_col <- grep("kegg", anno_cols, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(kegg_col)) {
      feature_anno$kegg_id <- annotation[[kegg_col]][idx]
    }

    # HMDB
    hmdb_col <- grep("hmdb", anno_cols, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(hmdb_col)) {
      feature_anno$hmdb_id <- annotation[[hmdb_col]][idx]
    }

    # ChEBI
    chebi_col <- grep("chebi", anno_cols, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(chebi_col)) {
      feature_anno$chebi_id <- annotation[[chebi_col]][idx]
    }

    # Compound name
    name_col <- grep("name|compound|metabolite", anno_cols, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(name_col)) {
      feature_anno$compound_name <- annotation[[name_col]][idx]
    }

    # Super class / class
    class_col <- grep("class|category", anno_cols, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(class_col)) {
      feature_anno$compound_class <- annotation[[class_col]][idx]
    }

    n_annotated <- sum(!is.na(feature_anno$compound_name) |
      !is.na(feature_anno$kegg_id) |
      !is.na(feature_anno$hmdb_id))
    log_message("Annotated ", n_annotated, "/", nrow(feature_anno), " metabolite features")
  }

  # Update DA table if present
  da_table <- metab_data$da_table
  if (!is.null(da_table) && "compound_name" %in% colnames(feature_anno)) {
    idx <- match(da_table$feature_id, feature_anno$feature_id)
    da_table$compound_name <- feature_anno$compound_name[idx]
  }

  list(
    normalized_matrix = mat,
    da_table = da_table,
    feature_annotation = feature_anno,
    pathway_mapping = metab_data$pathway_mapping,
    gmt = metab_data$gmt
  )
}

#' Create unified gene-centric mapping for cross-omics comparison
create_gene_mapping <- function(harmonized_data) {
  log_message("Creating unified gene-centric mapping...")

  harmonized <- harmonized_data$harmonized_omics
  gene_map <- list()

  # RNA: feature_id -> gene_symbol

  if ("transcriptomics" %in% names(harmonized)) {
    rna_anno <- harmonized$transcriptomics$feature_annotation
    if ("gene_symbol" %in% colnames(rna_anno)) {
      gene_map$rna <- data.frame(
        feature_id = rna_anno$feature_id,
        gene_symbol = rna_anno$gene_symbol,
        omics = "transcriptomics",
        stringsAsFactors = FALSE
      )
    }
  }

  # Protein: feature_id -> gene_symbol
  if ("proteomics" %in% names(harmonized)) {
    prot_anno <- harmonized$proteomics$feature_annotation
    if ("gene_symbol" %in% colnames(prot_anno)) {
      gene_map$protein <- data.frame(
        feature_id = prot_anno$feature_id,
        gene_symbol = prot_anno$gene_symbol,
        omics = "proteomics",
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine
  if (length(gene_map) > 0) {
    unified <- do.call(rbind, gene_map)
    unified <- unified[!is.na(unified$gene_symbol), ]

    log_message("Unified gene mapping: ", nrow(unified), " feature-gene pairs")
    return(unified)
  }

  return(NULL)
}
