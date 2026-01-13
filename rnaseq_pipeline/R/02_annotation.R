# R/02_annotation.R
# Functions for gene ID processing and annotation

#' Process gene IDs - strip version suffixes and handle duplicates
#' @param counts Counts matrix with gene IDs as rownames
#' @param strip_version Whether to strip Ensembl version suffixes
#' @return List with processed counts and ID mapping
process_gene_ids <- function(counts, strip_version = TRUE) {

  original_ids <- rownames(counts)

  if (strip_version) {
    # Strip version suffix (e.g., ENSG00000123456.1 -> ENSG00000123456)
    stripped_ids <- sub("\\.[0-9]+$", "", original_ids)

    # Create mapping table
    id_mapping <- data.frame(
      original_id = original_ids,
      stripped_id = stripped_ids,
      stringsAsFactors = FALSE
    )

    # Check for duplicates after stripping
    dup_ids <- stripped_ids[duplicated(stripped_ids)]

    if (length(dup_ids) > 0) {
      message("Found ", length(unique(dup_ids)), " gene IDs with duplicates after version stripping")
      message("Aggregating counts by summing duplicates")

      # Aggregate by summing counts for duplicate IDs
      counts_df <- as.data.frame(counts)
      counts_df$gene_id <- stripped_ids

      counts_agg <- counts_df %>%
        group_by(gene_id) %>%
        summarise(across(everything(), sum)) %>%
        as.data.frame()

      rownames(counts_agg) <- counts_agg$gene_id
      counts_agg$gene_id <- NULL
      counts <- as.matrix(counts_agg)

      # Update mapping to note aggregated IDs
      id_mapping$was_aggregated <- duplicated(stripped_ids) |
        duplicated(stripped_ids, fromLast = TRUE)
    } else {
      rownames(counts) <- stripped_ids
      id_mapping$was_aggregated <- FALSE
    }

    message("Processed gene IDs: ", nrow(counts), " unique genes")

  } else {
    id_mapping <- data.frame(
      original_id = original_ids,
      stripped_id = original_ids,
      was_aggregated = FALSE,
      stringsAsFactors = FALSE
    )
  }

  list(
    counts = counts,
    id_mapping = id_mapping
  )
}


#' Get organism database name for biomaRt/AnnotationHub
#' @param organism User-provided organism string
#' @return List with dataset name and other organism info
get_organism_info <- function(organism) {

  # Common organism mappings
  organism_map <- list(
    # Human
    "human" = list(
      dataset = "hsapiens_gene_ensembl",
      orgdb = "org.Hs.eg.db",
      kegg = "hsa",
      taxid = 9606
    ),
    "homo sapiens" = list(
      dataset = "hsapiens_gene_ensembl",
      orgdb = "org.Hs.eg.db",
      kegg = "hsa",
      taxid = 9606
    ),
    # Mouse
    "mouse" = list(
      dataset = "mmusculus_gene_ensembl",
      orgdb = "org.Mm.eg.db",
      kegg = "mmu",
      taxid = 10090
    ),
    "mus musculus" = list(
      dataset = "mmusculus_gene_ensembl",
      orgdb = "org.Mm.eg.db",
      kegg = "mmu",
      taxid = 10090
    ),
    # Rat
    "rat" = list(
      dataset = "rnorvegicus_gene_ensembl",
      orgdb = "org.Rn.eg.db",
      kegg = "rno",
      taxid = 10116
    ),
    "rattus norvegicus" = list(
      dataset = "rnorvegicus_gene_ensembl",
      orgdb = "org.Rn.eg.db",
      kegg = "rno",
      taxid = 10116
    ),
    # Zebrafish
    "zebrafish" = list(
      dataset = "drerio_gene_ensembl",
      orgdb = "org.Dr.eg.db",
      kegg = "dre",
      taxid = 7955
    ),
    "danio rerio" = list(
      dataset = "drerio_gene_ensembl",
      orgdb = "org.Dr.eg.db",
      kegg = "dre",
      taxid = 7955
    ),
    # Fly
    "fly" = list(
      dataset = "dmelanogaster_gene_ensembl",
      orgdb = "org.Dm.eg.db",
      kegg = "dme",
      taxid = 7227
    ),
    "drosophila melanogaster" = list(
      dataset = "dmelanogaster_gene_ensembl",
      orgdb = "org.Dm.eg.db",
      kegg = "dme",
      taxid = 7227
    ),
    # Worm
    "worm" = list(
      dataset = "celegans_gene_ensembl",
      orgdb = "org.Ce.eg.db",
      kegg = "cel",
      taxid = 6239
    ),
    "caenorhabditis elegans" = list(
      dataset = "celegans_gene_ensembl",
      orgdb = "org.Ce.eg.db",
      kegg = "cel",
      taxid = 6239
    ),
    # Yeast
    "yeast" = list(
      dataset = "scerevisiae_gene_ensembl",
      orgdb = "org.Sc.sgd.db",
      kegg = "sce",
      taxid = 559292
    ),
    "saccharomyces cerevisiae" = list(
      dataset = "scerevisiae_gene_ensembl",
      orgdb = "org.Sc.sgd.db",
      kegg = "sce",
      taxid = 559292
    ),
    # Arabidopsis
    "arabidopsis" = list(
      dataset = "athaliana_eg_gene",
      orgdb = "org.At.tair.db",
      kegg = "ath",
      taxid = 3702
    ),
    "arabidopsis thaliana" = list(
      dataset = "athaliana_eg_gene",
      orgdb = "org.At.tair.db",
      kegg = "ath",
      taxid = 3702
    )
  )

  # Normalize organism name
  org_lower <- tolower(trimws(organism))

  if (org_lower %in% names(organism_map)) {
    info <- organism_map[[org_lower]]
    info$supported = TRUE
    return(info)
  }

  # Return unsupported organism info
  list(
    dataset = NA,
    orgdb = NA,
    kegg = NA,
    taxid = NA,
    supported = FALSE
  )
}


#' Annotate genes using various sources
#' @param gene_ids Vector of gene IDs
#' @param organism Organism name
#' @param gene_id_type Type of gene IDs (ensembl_gene_id, etc.)
#' @param mapping_file Optional path to custom mapping file
#' @param annotation_source Source to use: "auto", "biomart", "annotationhub", "file"
#' @return List with annotation data frame and statistics
annotate_genes <- function(gene_ids,
                           organism,
                           gene_id_type = "ensembl_gene_id",
                           mapping_file = NULL,
                           annotation_source = "auto") {

  annotation <- data.frame(
    gene_id = gene_ids,
    stringsAsFactors = FALSE
  )

  # If custom mapping file provided, use it
  if (!is.null(mapping_file) && file.exists(mapping_file)) {
    message("Loading annotation from custom mapping file: ", mapping_file)

    custom_map <- read.csv(mapping_file, stringsAsFactors = FALSE)

    # Expect first column to be gene IDs
    id_col <- colnames(custom_map)[1]

    # Merge with annotation
    annotation <- annotation %>%
      left_join(custom_map, by = setNames(id_col, "gene_id"))

    stats <- list(
      source = "custom_file",
      total_genes = length(gene_ids),
      annotated_genes = sum(!is.na(annotation[[colnames(custom_map)[2]]])),
      success_rate = mean(!is.na(annotation[[colnames(custom_map)[2]]]))
    )

    message("Custom annotation: ", stats$annotated_genes, "/", stats$total_genes,
            " genes mapped (", round(stats$success_rate * 100, 1), "%)")

    return(list(annotation = annotation, stats = stats))
  }

  # Get organism info
  org_info <- get_organism_info(organism)

  if (!org_info$supported) {
    warning("Organism '", organism, "' not in standard database. ",
            "Annotation will be limited. Consider providing a mapping_file.")

    stats <- list(
      source = "none",
      total_genes = length(gene_ids),
      annotated_genes = 0,
      success_rate = 0,
      message = "Unsupported organism - no annotation available"
    )

    return(list(annotation = annotation, stats = stats))
  }

  # Try biomaRt
  if (annotation_source %in% c("auto", "biomart")) {
    tryCatch({
      message("Attempting annotation via biomaRt...")

      ensembl <- biomaRt::useEnsembl(
        biomart = "genes",
        dataset = org_info$dataset
      )

      # Define attributes to retrieve
      attrs <- c(
        gene_id_type,
        "external_gene_name",
        "description",
        "entrezgene_id",
        "gene_biotype"
      )

      # Query biomaRt
      biomart_result <- biomaRt::getBM(
        attributes = attrs,
        filters = gene_id_type,
        values = gene_ids,
        mart = ensembl
      )

      # Rename columns
      colnames(biomart_result) <- c(
        "gene_id", "symbol", "description", "entrez_id", "gene_biotype"
      )

      # Handle duplicates (keep first)
      biomart_result <- biomart_result %>%
        distinct(gene_id, .keep_all = TRUE)

      # Merge with annotation
      annotation <- annotation %>%
        left_join(biomart_result, by = "gene_id")

      stats <- list(
        source = "biomart",
        total_genes = length(gene_ids),
        annotated_genes = sum(!is.na(annotation$symbol)),
        success_rate = mean(!is.na(annotation$symbol))
      )

      message("biomaRt annotation: ", stats$annotated_genes, "/", stats$total_genes,
              " genes mapped (", round(stats$success_rate * 100, 1), "%)")

      return(list(annotation = annotation, stats = stats))

    }, error = function(e) {
      warning("biomaRt annotation failed: ", e$message)
    })
  }

  # Try AnnotationHub/OrgDb
  if (annotation_source %in% c("auto", "annotationhub") && !is.na(org_info$orgdb)) {
    tryCatch({
      message("Attempting annotation via OrgDb (", org_info$orgdb, ")...")

      # Load OrgDb
      if (!requireNamespace(org_info$orgdb, quietly = TRUE)) {
        stop("OrgDb package not installed: ", org_info$orgdb)
      }

      orgdb <- get(org_info$orgdb)

      # Map Ensembl to Symbol and Entrez
      gene_info <- AnnotationDbi::select(
        orgdb,
        keys = gene_ids,
        columns = c("SYMBOL", "GENENAME", "ENTREZID"),
        keytype = "ENSEMBL"
      )

      colnames(gene_info) <- c("gene_id", "symbol", "description", "entrez_id")

      # Handle duplicates
      gene_info <- gene_info %>%
        distinct(gene_id, .keep_all = TRUE)

      # Merge
      annotation <- annotation %>%
        left_join(gene_info, by = "gene_id")

      stats <- list(
        source = "orgdb",
        total_genes = length(gene_ids),
        annotated_genes = sum(!is.na(annotation$symbol)),
        success_rate = mean(!is.na(annotation$symbol))
      )

      message("OrgDb annotation: ", stats$annotated_genes, "/", stats$total_genes,
              " genes mapped (", round(stats$success_rate * 100, 1), "%)")

      return(list(annotation = annotation, stats = stats))

    }, error = function(e) {
      warning("OrgDb annotation failed: ", e$message)
    })
  }

  # No annotation available
  warning("Could not obtain gene annotation. Pipeline will continue with gene IDs only.")

  stats <- list(
    source = "none",
    total_genes = length(gene_ids),
    annotated_genes = 0,
    success_rate = 0,
    message = "Annotation failed - proceeding with gene IDs only"
  )

  list(annotation = annotation, stats = stats)
}
