# =============================================================================
# 01_data_ingestion.R — Excel/CSV loading, metadata, sample alignment
# =============================================================================

ingest_data <- function(config) {
  log_message("=== Starting Data Ingestion ===")
  create_output_dirs(config)

  # Load feature matrix (+ split annotation columns)
  raw <- load_feature_matrix(config)
  mat              <- raw$matrix
  feature_metadata <- raw$feature_metadata
  group_row        <- raw$group_row          # may be NULL

  # Load metadata
  metadata <- load_metadata(config, group_row, colnames(mat))

  # Align samples
  aligned  <- align_samples(mat, metadata, config)
  mat      <- aligned$matrix
  metadata <- aligned$metadata

  log_message("=== Data Ingestion Complete: ", nrow(mat), " features x ",
              ncol(mat), " samples ===")

  list(
    matrix           = mat,
    metadata         = metadata,
    feature_metadata = feature_metadata
  )
}

# ---------- feature matrix ---------------------------------------------------

load_feature_matrix <- function(config) {
  path <- config$input$feature_matrix
  if (!file.exists(path)) stop("Feature matrix not found: ", path)

  ext <- tolower(tools::file_ext(path))
  has_group_row <- isTRUE(config$input$has_group_row)

  # --- read raw lines for group-row detection --------------------------------
  group_row <- NULL

  if (ext %in% c("xlsx", "xls")) {
    sheet <- config$input$feature_sheet %||% 1
    df <- as.data.frame(readxl::read_excel(path, sheet = sheet),
                        stringsAsFactors = FALSE)
  } else {
    # For TSV/CSV: if has_group_row, read header + group row separately
    if (has_group_row) {
      raw_lines <- readLines(path, n = 2)
      header    <- strsplit(raw_lines[1], "\t")[[1]]
      grp_vals  <- strsplit(raw_lines[2], "\t")[[1]]
      # Read data starting from line 3 (skip header + group row)
      df <- read.delim(path, header = FALSE, skip = 2,
                       stringsAsFactors = FALSE, check.names = FALSE)
      colnames(df) <- header

      # Build group_row: map sample columns → group labels
      annot_cols_cfg <- config$input$annotation_columns %||%
        c("Molecule", "HMDB")
      n_annot <- length(intersect(annot_cols_cfg, header))
      sample_headers <- header[(n_annot + 1):length(header)]
      grp_labels     <- grp_vals[(n_annot + 1):length(grp_vals)]
      group_row <- setNames(grp_labels, sample_headers)
    } else if (ext == "tsv" || ext == "txt") {
      df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      df <- as.data.frame(readr::read_csv(path, show_col_types = FALSE),
                          stringsAsFactors = FALSE)
    }
  }

  # Separate annotation columns from intensity data
  annot_cols <- config$input$annotation_columns %||%
    c("Molecule", "HMDB", "SMILES", "KEGG", "CAS")
  annot_present <- intersect(annot_cols, colnames(df))
  sample_cols   <- setdiff(colnames(df), annot_present)

  # Build feature IDs from HMDB|Molecule (matching notebook convention)
  if (all(c("HMDB", "Molecule") %in% annot_present)) {
    feature_ids <- paste0(df$HMDB, "|", df$Molecule)
  } else {
    feature_ids <- as.character(df[[1]])
    sample_cols <- sample_cols[-1]
  }

  # Feature metadata
  if (length(annot_present) > 0) {
    feature_metadata <- df[, annot_present, drop = FALSE]
    feature_metadata$feature_id <- feature_ids
  } else {
    feature_metadata <- data.frame(feature_id = feature_ids,
                                    stringsAsFactors = FALSE)
  }

  # Intensity matrix
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  rownames(mat) <- feature_ids
  storage.mode(mat) <- "numeric"

  log_message("Loaded matrix: ", nrow(mat), " features x ", ncol(mat),
              " columns (", length(annot_present), " annotation cols split)")

  list(matrix = mat, feature_metadata = feature_metadata, group_row = group_row)
}

# ---------- metadata ---------------------------------------------------------

load_metadata <- function(config, group_row = NULL, sample_names = NULL) {
  path <- config$input$metadata
  sample_col    <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  if (!is.null(path) && file.exists(path)) {
    # --- external metadata file ---
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx", "xls")) {
      sheet <- config$input$metadata_sheet %||% 1
      metadata <- as.data.frame(
        readxl::read_excel(path, sheet = sheet),
        stringsAsFactors = FALSE
      )
    } else if (ext %in% c("tsv", "txt")) {
      metadata <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      metadata <- as.data.frame(
        readr::read_csv(path, show_col_types = FALSE),
        stringsAsFactors = FALSE
      )
    }
  } else {
    # --- try embedded Excel sheet ---
    feat_path <- config$input$feature_matrix
    feat_ext  <- tolower(tools::file_ext(feat_path))
    if (feat_ext %in% c("xlsx", "xls")) {
      sheet <- config$input$metadata_sheet %||% "SampleSheet"
      log_message("Loading metadata from sheet '", sheet, "' of ", feat_path)
      metadata <- as.data.frame(
        readxl::read_excel(feat_path, sheet = sheet),
        stringsAsFactors = FALSE
      )
    } else if (!is.null(group_row)) {
      # --- auto-generate metadata from group row ---
      log_message("Auto-generating metadata from group row")
      metadata <- data.frame(
        SampleName = names(group_row),
        condition  = unname(group_row),
        stringsAsFactors = FALSE
      )
      colnames(metadata) <- c(sample_col, condition_col)
    } else if (!is.null(sample_names)) {
      # --- infer condition from sample name pattern (GT*_YEAR_*) ---
      log_message("Inferring metadata from sample names")
      cond <- sapply(sample_names, function(s) {
        strsplit(s, "_")[[1]][2]
      }, USE.NAMES = FALSE)
      metadata <- data.frame(
        SampleName = sample_names,
        condition  = cond,
        stringsAsFactors = FALSE
      )
      colnames(metadata) <- c(sample_col, condition_col)
    } else {
      stop("No metadata source available: set input$metadata, ",
           "use Excel with a SampleSheet, or enable has_group_row")
    }
  }

  if (!sample_col %in% colnames(metadata))
    stop("Sample ID column '", sample_col, "' not found in metadata. ",
         "Available: ", paste(colnames(metadata), collapse = ", "))

  if (!condition_col %in% colnames(metadata))
    stop("Condition column '", condition_col, "' not found in metadata. ",
         "Available: ", paste(colnames(metadata), collapse = ", "))

  metadata[[condition_col]] <- as.factor(metadata[[condition_col]])

  ref <- config$design$reference_level
  if (!is.null(ref) && ref %in% levels(metadata[[condition_col]])) {
    metadata[[condition_col]] <- relevel(metadata[[condition_col]], ref = ref)
  }

  log_message("Loaded metadata: ", nrow(metadata), " samples, ",
              nlevels(metadata[[condition_col]]), " condition levels")
  metadata
}

# ---------- sample alignment -------------------------------------------------

align_samples <- function(mat, metadata, config) {
  sample_col <- config$input$sample_id_column
  mat_samples  <- colnames(mat)
  meta_samples <- metadata[[sample_col]]

  common <- intersect(mat_samples, meta_samples)
  if (length(common) == 0)
    stop("No matching samples between matrix and metadata!\n",
         "Matrix: ", paste(head(mat_samples, 5), collapse = ", "), "\n",
         "Metadata: ", paste(head(meta_samples, 5), collapse = ", "))

  only_mat  <- setdiff(mat_samples, meta_samples)
  only_meta <- setdiff(meta_samples, mat_samples)
  if (length(only_mat) > 0)
    log_message("Dropping ", length(only_mat),
                " matrix columns not in metadata: ",
                paste(head(only_mat, 5), collapse = ", "))
  if (length(only_meta) > 0)
    log_message("Dropping ", length(only_meta),
                " metadata rows not in matrix: ",
                paste(head(only_meta, 5), collapse = ", "))

  mat      <- mat[, common, drop = FALSE]
  metadata <- metadata[match(common, metadata[[sample_col]]), ]

  log_message("Aligned ", length(common), " samples")
  list(matrix = mat, metadata = metadata)
}
