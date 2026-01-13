# =============================================================================
# Identifier Harmonization / Mapping Layer
# =============================================================================

#' Harmonize identifiers across all omics
harmonize_identifiers <- function(preprocessed_data, config) {
  log_message("=== Harmonizing Identifiers Across Omics ===")

  processed_omics <- preprocessed_data$processed_omics
  harmonized <- list()

  # Transcriptomics: Ensembl -> Gene Symbol
  if ("transcriptomics" %in% names(processed_omics)) {
    harmonized$transcriptomics <- harmonize_transcriptomics_ids(
      processed_omics$transcriptomics, config
    )
  }

  # Proteomics: UniProt/Accession -> Gene Symbol
  if ("proteomics" %in% names(processed_omics)) {
    harmonized$proteomics <- harmonize_proteomics_ids(
      processed_omics$proteomics, config
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
harmonize_transcriptomics_ids <- function(rna_data, config) {
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

  # Apply mapping if provided
  if (!is.null(mapping) && nrow(mapping) > 0) {
    # Expect mapping to have: ensembl_id, gene_symbol, entrez_id (optional)
    mapping_cols <- colnames(mapping)

    # Find matching columns
    ensembl_col <- grep("ensembl|gene_id", mapping_cols, ignore.case = TRUE, value = TRUE)[1]
    symbol_col <- grep("symbol|gene_name|hgnc", mapping_cols, ignore.case = TRUE, value = TRUE)[1]
    entrez_col <- grep("entrez|ncbi", mapping_cols, ignore.case = TRUE, value = TRUE)[1]

    if (!is.na(ensembl_col) && !is.na(symbol_col)) {
      # Strip versions from mapping too
      mapping[[ensembl_col]] <- strip_ensembl_version(mapping[[ensembl_col]])

      # Match
      idx <- match(stripped_ids, mapping[[ensembl_col]])
      feature_anno$gene_symbol <- mapping[[symbol_col]][idx]

      if (!is.na(entrez_col)) {
        feature_anno$entrez_id <- mapping[[entrez_col]][idx]
      }

      n_mapped <- sum(!is.na(feature_anno$gene_symbol))
      log_message("Mapped ", n_mapped, "/", nrow(feature_anno), " genes to symbols")
    }
  } else {
    # Try to infer if IDs are already gene symbols
    if (!any(grepl("^ENS", stripped_ids))) {
      feature_anno$gene_symbol <- stripped_ids
      log_message("IDs appear to be gene symbols, using as-is")
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
harmonize_proteomics_ids <- function(prot_data, config) {
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
      log_message("Mapped ", n_mapped, "/", nrow(feature_anno), " proteins to gene symbols")
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
    # Plain UniProt accession (P12345 or Q12345-1)
    else if (grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]", id)) {
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
    id_col <- meta_cols[1]  # Assume first column is ID

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
