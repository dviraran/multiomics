# =============================================================================
# Annotation and ID Mapping
# =============================================================================

#' Map feature IDs to gene symbols and other identifiers
#'
#' @param ingested_data Data from ingestion step
#' @param config Configuration list
#' @return Annotation data frame
annotate_features <- function(ingested_data, config) {
  log_message("=== Starting Feature Annotation ===")

  feature_ids <- rownames(ingested_data$matrix)
  feature_type <- config$input$feature_id_type
  organism <- config$input$organism

  # Initialize annotation data frame
  annotation_df <- data.frame(
    feature_id = feature_ids,
    original_id = feature_ids,
    stringsAsFactors = FALSE
  )

  # Parse IDs if UniProt format
  if (feature_type %in% c("uniprot", "protein_accession")) {
    annotation_df$uniprot_accession <- parse_uniprot_accession(feature_ids)
  } else {
    annotation_df$uniprot_accession <- feature_ids
  }

  # Check for custom mapping file first
  if (!is.null(config$optional_inputs$mapping_file) &&
      file.exists(config$optional_inputs$mapping_file)) {
    log_message("Using custom mapping file: ", config$optional_inputs$mapping_file)
    annotation_df <- apply_custom_mapping(annotation_df, config)
  } else {
    # Attempt database mapping
    annotation_df <- attempt_database_mapping(annotation_df, organism, config)
  }

  # Calculate mapping success
  mapping_summary <- create_mapping_summary(annotation_df)

  # Save annotation table
  save_table(annotation_df, "feature_annotations.csv", config, "tables")
  save_table(mapping_summary, "mapping_summary.csv", config, "qc")

  log_message("=== Feature Annotation Complete ===")

  list(
    annotations = annotation_df,
    mapping_summary = mapping_summary
  )
}

#' Apply custom mapping file
#'
#' @param annotation_df Current annotation data frame
#' @param config Configuration list
#' @return Updated annotation data frame
apply_custom_mapping <- function(annotation_df, config) {
  mapping <- readr::read_csv(config$optional_inputs$mapping_file, show_col_types = FALSE)

  # Find the ID column in mapping file
  id_cols <- c("feature_id", "uniprot", "accession", "protein_id")
  id_col <- intersect(tolower(colnames(mapping)), id_cols)[1]

  if (is.na(id_col)) {
    warning("Could not identify ID column in mapping file. Using first column.")
    id_col <- colnames(mapping)[1]
  }

  # Match IDs
  match_idx <- match(annotation_df$uniprot_accession, mapping[[id_col]])
  alt_match <- match(annotation_df$feature_id, mapping[[id_col]])
  match_idx[is.na(match_idx)] <- alt_match[is.na(match_idx)]

  # Add available columns
  available_cols <- c("gene_symbol", "entrez_id", "protein_name", "description",
                      "gene_name", "symbol", "entrez", "gene_id")

  for (col in intersect(tolower(colnames(mapping)), available_cols)) {
    original_col <- colnames(mapping)[tolower(colnames(mapping)) == col][1]
    new_col <- gsub("gene_name|symbol", "gene_symbol", col)
    new_col <- gsub("gene_id|entrez", "entrez_id", new_col)

    annotation_df[[new_col]] <- mapping[[original_col]][match_idx]
  }

  log_message("Applied custom mapping to ", sum(!is.na(match_idx)), " features")

  return(annotation_df)
}

#' Attempt to map IDs using online databases
#'
#' @param annotation_df Current annotation data frame
#' @param organism Organism name
#' @param config Configuration list
#' @return Updated annotation data frame
attempt_database_mapping <- function(annotation_df, organism, config) {
  mapping_service <- config$annotation$mapping_service %||% "UniProt.ws"

  # Get valid UniProt accessions
  valid_accessions <- annotation_df$uniprot_accession[!is.na(annotation_df$uniprot_accession)]

  if (length(valid_accessions) == 0) {
    log_message("No valid UniProt accessions found for mapping")
    annotation_df$gene_symbol <- NA
    annotation_df$entrez_id <- NA
    annotation_df$protein_name <- NA
    return(annotation_df)
  }

  log_message("Attempting database mapping for ", length(valid_accessions), " features")

  # Try UniProt.ws first
  if (mapping_service == "UniProt.ws") {
    tryCatch({
      annotation_df <- map_with_uniprot_ws(annotation_df, config)
    }, error = function(e) {
      log_message("UniProt.ws mapping failed: ", e$message)
      log_message("Falling back to biomaRt...")
      annotation_df <<- map_with_biomart(annotation_df, organism, config)
    })
  } else if (mapping_service == "biomaRt") {
    annotation_df <- map_with_biomart(annotation_df, organism, config)
  } else {
    log_message("Unknown mapping service. Skipping database mapping.")
    annotation_df$gene_symbol <- NA
    annotation_df$entrez_id <- NA
    annotation_df$protein_name <- NA
  }

  return(annotation_df)
}

#' Map IDs using UniProt.ws
#'
#' @param annotation_df Annotation data frame
#' @param config Configuration list
#' @return Updated annotation data frame
map_with_uniprot_ws <- function(annotation_df, config) {
  if (!requireNamespace("UniProt.ws", quietly = TRUE)) {
    stop("UniProt.ws package not available. Install with BiocManager::install('UniProt.ws')")
  }

  log_message("Querying UniProt database...")

  # Get valid accessions
  valid_idx <- !is.na(annotation_df$uniprot_accession)
  accessions <- unique(annotation_df$uniprot_accession[valid_idx])

  # Remove isoform suffixes for query
  base_accessions <- gsub("-[0-9]+$", "", accessions)

  # Chunk the queries to avoid timeout
  chunk_size <- config$annotation$query_chunk_size %||% 500
  chunks <- split(base_accessions, ceiling(seq_along(base_accessions) / chunk_size))

  results_list <- list()

  for (i in seq_along(chunks)) {
    log_message("Querying chunk ", i, " of ", length(chunks), "...")

    tryCatch({
      up <- UniProt.ws::UniProt.ws(taxId = 9606)  # Human, will work for cross-species queries

      result <- UniProt.ws::select(
        up,
        keys = chunks[[i]],
        columns = c("gene_primary", "protein_name", "xref_geneid"),
        keytype = "UniProtKB"
      )

      results_list[[i]] <- result
    }, error = function(e) {
      log_message("Chunk ", i, " failed: ", e$message)
    })
  }

  # Combine results
  if (length(results_list) > 0) {
    all_results <- do.call(rbind, results_list)

    # Map back to annotation_df
    base_acc <- gsub("-[0-9]+$", "", annotation_df$uniprot_accession)
    match_idx <- match(base_acc, all_results$From)

    annotation_df$gene_symbol <- all_results$gene_primary[match_idx]
    annotation_df$protein_name <- all_results$protein_name[match_idx]
    annotation_df$entrez_id <- all_results$xref_geneid[match_idx]

    log_message("UniProt mapping successful for ", sum(!is.na(annotation_df$gene_symbol)), " features")
  } else {
    log_message("UniProt mapping returned no results")
    annotation_df$gene_symbol <- NA
    annotation_df$entrez_id <- NA
    annotation_df$protein_name <- NA
  }

  return(annotation_df)
}

#' Map IDs using biomaRt
#'
#' @param annotation_df Annotation data frame
#' @param organism Organism name
#' @param config Configuration list
#' @return Updated annotation data frame
map_with_biomart <- function(annotation_df, organism, config) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    log_message("biomaRt package not available. Skipping database mapping.")
    annotation_df$gene_symbol <- NA
    annotation_df$entrez_id <- NA
    annotation_df$protein_name <- NA
    return(annotation_df)
  }

  log_message("Querying biomaRt...")

  # Determine dataset based on organism
  dataset <- get_biomart_dataset(organism, config)

  if (is.null(dataset)) {
    log_message("No biomaRt dataset found for organism: ", organism)
    annotation_df$gene_symbol <- NA
    annotation_df$entrez_id <- NA
    annotation_df$protein_name <- NA
    return(annotation_df)
  }

  tryCatch({
    mart <- biomaRt::useMart("ensembl", dataset = dataset)

    # Get valid accessions
    valid_idx <- !is.na(annotation_df$uniprot_accession)
    accessions <- unique(annotation_df$uniprot_accession[valid_idx])
    base_accessions <- gsub("-[0-9]+$", "", accessions)

    results <- biomaRt::getBM(
      attributes = c("uniprotswissprot", "external_gene_name", "entrezgene_id", "description"),
      filters = "uniprotswissprot",
      values = base_accessions,
      mart = mart
    )

    # Map back
    base_acc <- gsub("-[0-9]+$", "", annotation_df$uniprot_accession)
    match_idx <- match(base_acc, results$uniprotswissprot)

    annotation_df$gene_symbol <- results$external_gene_name[match_idx]
    annotation_df$entrez_id <- results$entrezgene_id[match_idx]
    annotation_df$protein_name <- results$description[match_idx]

    log_message("biomaRt mapping successful for ", sum(!is.na(annotation_df$gene_symbol)), " features")
  }, error = function(e) {
    log_message("biomaRt mapping failed: ", e$message)
    annotation_df$gene_symbol <- NA
    annotation_df$entrez_id <- NA
    annotation_df$protein_name <- NA
  })

  return(annotation_df)
}

#' Get biomaRt dataset name for organism
#'
#' @param organism Organism name
#' @param config Configuration list
#' @return Dataset name or NULL
get_biomart_dataset <- function(organism, config) {
  # Check if explicitly specified in config
  if (!is.null(config$annotation$biomart_dataset)) {
    return(config$annotation$biomart_dataset)
  }

  # Common mappings
  organism_map <- list(
    "homo sapiens" = "hsapiens_gene_ensembl",
    "human" = "hsapiens_gene_ensembl",
    "mus musculus" = "mmusculus_gene_ensembl",
    "mouse" = "mmusculus_gene_ensembl",
    "rattus norvegicus" = "rnorvegicus_gene_ensembl",
    "rat" = "rnorvegicus_gene_ensembl",
    "danio rerio" = "drerio_gene_ensembl",
    "zebrafish" = "drerio_gene_ensembl",
    "drosophila melanogaster" = "dmelanogaster_gene_ensembl",
    "saccharomyces cerevisiae" = "scerevisiae_gene_ensembl",
    "yeast" = "scerevisiae_gene_ensembl",
    "arabidopsis thaliana" = "athaliana_eg_gene",
    "caenorhabditis elegans" = "celegans_gene_ensembl"
  )

  return(organism_map[[tolower(organism)]])
}

#' Create mapping summary statistics
#'
#' @param annotation_df Annotation data frame
#' @return Summary data frame
create_mapping_summary <- function(annotation_df) {
  n_total <- nrow(annotation_df)

  summary_df <- data.frame(
    metric = c(
      "Total features",
      "Valid UniProt accessions",
      "Gene symbols mapped",
      "Entrez IDs mapped",
      "Protein names mapped",
      "Mapping success rate (%)"
    ),
    value = c(
      n_total,
      sum(!is.na(annotation_df$uniprot_accession)),
      sum(!is.na(annotation_df$gene_symbol)),
      sum(!is.na(annotation_df$entrez_id)),
      sum(!is.na(annotation_df$protein_name)),
      round(sum(!is.na(annotation_df$gene_symbol)) / n_total * 100, 1)
    ),
    stringsAsFactors = FALSE
  )

  return(summary_df)
}
