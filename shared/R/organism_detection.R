#' Organism and ID Type Detection Utilities
#'
#' This module provides functions to automatically detect:
#' - Organism from gene/protein IDs
#' - ID type (Ensembl, Entrez, UniProt, symbol, etc.)
#' - ID patterns for validation
#'
#' @description
#' Works with non-model organisms by detecting patterns rather than
#' requiring database lookups.

# ==============================================================================
# ID TYPE DETECTION
# ==============================================================================

#' Detect the type of gene/protein IDs
#'
#' @param ids Character vector of IDs to analyze
#' @param n_sample Number of IDs to sample for detection (default 100)
#'
#' @return Character string: "ensembl_gene_id", "ensembl_transcript_id",
#'         "entrez", "symbol", "uniprot", "refseq", or "unknown"
#'
#' @details
#' Detection is based on ID patterns:
#' - Ensembl gene: ENSG, ENSMUSG, ENSRNOG, etc.
#' - Ensembl transcript: ENST, ENSMUST, etc.
#' - Ensembl protein: ENSP, ENSMUSP, etc.
#' - Entrez: numeric only
#' - UniProt: [A-Z][0-9][A-Z0-9]{3}[0-9] pattern
#' - RefSeq: NM_, NR_, XM_, NP_, XP_ prefixes
#' - Symbol: typically short alphanumeric
#'
#' @export
detect_id_type <- function(ids, n_sample = 100) {

    if (length(ids) == 0) return("unknown")

    # Sample IDs for efficiency
    if (length(ids) > n_sample) {
        ids <- sample(ids, n_sample)
    }

    # Remove NA and empty
    ids <- ids[!is.na(ids) & ids != ""]
    if (length(ids) == 0) return("unknown")

    # Count matches for each pattern
    scores <- list(
        ensembl_gene_id = sum(grepl("^ENS[A-Z]*G[0-9]+", ids, ignore.case = FALSE)),
        ensembl_transcript_id = sum(grepl("^ENS[A-Z]*T[0-9]+", ids, ignore.case = FALSE)),
        ensembl_protein_id = sum(grepl("^ENS[A-Z]*P[0-9]+", ids, ignore.case = FALSE)),
        entrez = sum(grepl("^[0-9]+$", ids)),
        uniprot = sum(grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", ids)),
        refseq_mrna = sum(grepl("^[NX][MR]_[0-9]+", ids)),
        refseq_protein = sum(grepl("^[NXY]P_[0-9]+", ids)),
        wormbase = sum(grepl("^WBGene[0-9]+", ids)),
        flybase = sum(grepl("^FBgn[0-9]+", ids)),
        sgd = sum(grepl("^Y[A-P][LR][0-9]+[WC]", ids)),  # Yeast systematic names
        tair = sum(grepl("^AT[1-5MC]G[0-9]+", ids)),  # Arabidopsis
        mgi = sum(grepl("^MGI:[0-9]+", ids)),  # Mouse MGI
        rgd = sum(grepl("^RGD:[0-9]+", ids)),  # Rat RGD
        hgnc = sum(grepl("^HGNC:[0-9]+", ids)),  # Human HGNC
        ncbi_locus = sum(grepl("^LOC[0-9]+", ids))  # NCBI LOC IDs
    )

    # Find best match
    best <- names(which.max(unlist(scores)))
    best_score <- max(unlist(scores))

    # If best match covers <20% of IDs, might be gene symbols
    if (best_score / length(ids) < 0.2) {
        # Check if they look like gene symbols (short, alphanumeric, some uppercase)
        avg_len <- mean(nchar(ids))
        has_letters <- mean(grepl("[A-Za-z]", ids))
        mostly_upper <- mean(grepl("^[A-Z0-9-]+$", ids))

        if (avg_len < 15 && has_letters > 0.8 && mostly_upper > 0.5) {
            return("symbol")
        }
    }

    # If very high entrez score, but also could be other numeric IDs
    if (best == "entrez" && best_score / length(ids) > 0.9) {
        # Check if numbers are in typical Entrez range (1 to ~1000000)
        numeric_ids <- as.numeric(ids[grepl("^[0-9]+$", ids)])
        if (all(numeric_ids > 0 & numeric_ids < 10000000, na.rm = TRUE)) {
            return("entrez")
        }
    }

    # Return best match if it covers >50% of IDs
    if (best_score / length(ids) > 0.5) {
        return(best)
    }

    # Default to unknown if no clear pattern
    "unknown"
}

# ==============================================================================
# ORGANISM DETECTION FROM IDS
# ==============================================================================

#' Detect organism from gene IDs
#'
#' @param ids Character vector of gene IDs
#' @param n_sample Number of IDs to sample (default 100)
#'
#' @return List with:
#'   - organism: detected organism name or "unknown"
#'   - confidence: "high", "medium", "low", or "none"
#'   - evidence: description of how organism was detected
#'
#' @details
#' Detection based on species-specific ID prefixes:
#' - ENSG = Human
#' - ENSMUSG = Mouse
#' - ENSRNOG = Rat
#' - FBgn = Fly
#' - WBGene = Worm
#' - etc.
#'
#' @export
detect_organism_from_ids <- function(ids, n_sample = 100) {

    result <- list(
        organism = "unknown",
        confidence = "none",
        evidence = "No organism-specific patterns detected"
    )

    if (length(ids) == 0) return(result)

    # Sample for efficiency
    if (length(ids) > n_sample) {
        ids <- sample(ids, n_sample)
    }

    ids <- ids[!is.na(ids) & ids != ""]
    if (length(ids) == 0) return(result)

    # Species-specific Ensembl prefixes
    ensembl_patterns <- list(
        "Homo sapiens" = "^ENSG[0-9]",
        "Mus musculus" = "^ENSMUSG[0-9]",
        "Rattus norvegicus" = "^ENSRNOG[0-9]",
        "Danio rerio" = "^ENSDARG[0-9]",
        "Gallus gallus" = "^ENSGALG[0-9]",
        "Sus scrofa" = "^ENSSSCG[0-9]",
        "Bos taurus" = "^ENSBTAG[0-9]",
        "Ovis aries" = "^ENSOARG[0-9]",
        "Canis familiaris" = "^ENSCAFG[0-9]",
        "Felis catus" = "^ENSFCAG[0-9]",
        "Pan troglodytes" = "^ENSPTNG[0-9]",
        "Macaca mulatta" = "^ENSMMUG[0-9]",
        "Equus caballus" = "^ENSECAG[0-9]",
        "Oryctolagus cuniculus" = "^ENSOCUG[0-9]"
    )

    # Other species-specific patterns
    other_patterns <- list(
        "Drosophila melanogaster" = "^FBgn[0-9]",
        "Caenorhabditis elegans" = "^WBGene[0-9]",
        "Saccharomyces cerevisiae" = "^Y[A-P][LR][0-9]+[WC]",
        "Arabidopsis thaliana" = "^AT[1-5MC]G[0-9]"
    )

    # Check Ensembl patterns
    for (org in names(ensembl_patterns)) {
        pattern <- ensembl_patterns[[org]]
        n_match <- sum(grepl(pattern, ids))
        if (n_match / length(ids) > 0.5) {
            result$organism <- org
            result$confidence <- if (n_match / length(ids) > 0.8) "high" else "medium"
            result$evidence <- sprintf("%.0f%% of IDs match Ensembl %s pattern (%s)",
                                        100 * n_match / length(ids),
                                        gsub("_", " ", org),
                                        pattern)
            return(result)
        }
    }

    # Check other patterns
    for (org in names(other_patterns)) {
        pattern <- other_patterns[[org]]
        n_match <- sum(grepl(pattern, ids))
        if (n_match / length(ids) > 0.5) {
            result$organism <- org
            result$confidence <- if (n_match / length(ids) > 0.8) "high" else "medium"
            result$evidence <- sprintf("%.0f%% of IDs match %s pattern (%s)",
                                        100 * n_match / length(ids),
                                        gsub("_", " ", org),
                                        pattern)
            return(result)
        }
    }

    # Check for LOC IDs (NCBI - could be any organism)
    n_loc <- sum(grepl("^LOC[0-9]+", ids))
    if (n_loc / length(ids) > 0.5) {
        result$organism <- "unknown"
        result$confidence <- "low"
        result$evidence <- sprintf("%.0f%% of IDs are NCBI LOC format (organism-agnostic)",
                                    100 * n_loc / length(ids))
        return(result)
    }

    # Check for gene symbols (could be any organism)
    id_type <- detect_id_type(ids)
    if (id_type == "symbol") {
        # Try to guess from common human/mouse gene patterns
        human_markers <- c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "KRAS")
        mouse_markers <- c("Brca1", "Brca2", "Trp53", "Egfr", "Myc", "Kras")

        human_style <- sum(grepl("^[A-Z][A-Z0-9]+$", ids)) / length(ids)
        mouse_style <- sum(grepl("^[A-Z][a-z0-9]+$", ids)) / length(ids)

        if (human_style > 0.7) {
            result$organism <- "likely Homo sapiens"
            result$confidence <- "low"
            result$evidence <- "Gene symbols appear to be human-style (ALL CAPS)"
        } else if (mouse_style > 0.7) {
            result$organism <- "likely Mus musculus"
            result$confidence <- "low"
            result$evidence <- "Gene symbols appear to be mouse-style (Title Case)"
        } else {
            result$evidence <- "Gene symbols detected but organism unclear from capitalization"
        }
    }

    result
}

# ==============================================================================
# ID VALIDATION
# ==============================================================================

#' Validate IDs and report issues
#'
#' @param ids Character vector of IDs
#' @param expected_type Expected ID type (optional)
#' @param verbose Print messages
#'
#' @return List with validation results
#' @export
validate_ids <- function(ids, expected_type = NULL, verbose = TRUE) {

    result <- list(
        valid = TRUE,
        n_total = length(ids),
        n_unique = length(unique(ids)),
        n_na = sum(is.na(ids)),
        n_empty = sum(ids == "", na.rm = TRUE),
        n_duplicated = sum(duplicated(ids)),
        detected_type = detect_id_type(ids),
        issues = character(0),
        warnings = character(0)
    )

    # Check for NAs
    if (result$n_na > 0) {
        result$warnings <- c(result$warnings,
                             sprintf("%d NA values in IDs", result$n_na))
    }

    # Check for empty strings
    if (result$n_empty > 0) {
        result$warnings <- c(result$warnings,
                             sprintf("%d empty strings in IDs", result$n_empty))
    }

    # Check for duplicates
    if (result$n_duplicated > 0) {
        result$warnings <- c(result$warnings,
                             sprintf("%d duplicated IDs", result$n_duplicated))
    }

    # Check if type matches expected
    if (!is.null(expected_type) && expected_type != "auto") {
        if (result$detected_type != expected_type && result$detected_type != "unknown") {
            result$warnings <- c(result$warnings,
                                 sprintf("Expected %s IDs but detected %s",
                                         expected_type, result$detected_type))
        }
    }

    # Check for mixed ID types
    if (result$detected_type == "unknown") {
        # Sample different patterns
        patterns_found <- c(
            ensembl = sum(grepl("^ENS", ids)),
            entrez = sum(grepl("^[0-9]+$", ids)),
            symbol = sum(grepl("^[A-Za-z][A-Za-z0-9-]*$", ids) & nchar(ids) < 20),
            uniprot = sum(grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]", ids))
        )

        if (sum(patterns_found > length(ids) * 0.1) > 1) {
            result$warnings <- c(result$warnings,
                                 "Possible mixed ID types detected - consider standardizing")
        }
    }

    # Set valid flag
    if (length(result$issues) > 0) {
        result$valid <- FALSE
    }

    # Print summary
    if (verbose) {
        message(sprintf("\n=== ID Validation ==="))
        message(sprintf("Total IDs: %d", result$n_total))
        message(sprintf("Unique IDs: %d", result$n_unique))
        message(sprintf("Detected type: %s", result$detected_type))

        if (length(result$warnings) > 0) {
            message("\nWarnings:")
            for (w in result$warnings) {
                message(sprintf("  - %s", w))
            }
        }

        if (length(result$issues) > 0) {
            message("\nIssues (may cause problems):")
            for (i in result$issues) {
                message(sprintf("  - %s", i))
            }
        }
    }

    result
}

# ==============================================================================
# ID CONVERSION HELPERS
# ==============================================================================

#' Strip version suffix from Ensembl IDs
#'
#' @param ids Character vector of IDs
#' @return Character vector with version suffixes removed
#' @export
strip_ensembl_version <- function(ids) {
    gsub("\\.[0-9]+$", "", ids)
}

#' Check if IDs have Ensembl version suffixes
#'
#' @param ids Character vector of IDs
#' @return Logical
#' @export
has_ensembl_version <- function(ids) {
    # Sample for efficiency
    sample_ids <- if (length(ids) > 100) sample(ids, 100) else ids
    sample_ids <- sample_ids[!is.na(sample_ids)]

    # Check for pattern like ENSG00000123456.1
    versioned <- grepl("^ENS[A-Z]*[GPT][0-9]+\\.[0-9]+$", sample_ids)
    mean(versioned) > 0.5
}

# ==============================================================================
# ORGANISM NAME UTILITIES
# ==============================================================================

#' Normalize organism name to standard format
#'
#' @param organism Character, organism name in any format
#' @return Character, standardized organism name
#' @export
normalize_organism_name <- function(organism) {
    if (is.null(organism) || is.na(organism)) return(NULL)

    organism <- trimws(tolower(organism))

    # Common aliases
    aliases <- list(
        "human" = "Homo sapiens",
        "homo sapiens" = "Homo sapiens",
        "hsapiens" = "Homo sapiens",
        "hsa" = "Homo sapiens",
        "mouse" = "Mus musculus",
        "mus musculus" = "Mus musculus",
        "mmusculus" = "Mus musculus",
        "mmu" = "Mus musculus",
        "rat" = "Rattus norvegicus",
        "rattus norvegicus" = "Rattus norvegicus",
        "rnorvegicus" = "Rattus norvegicus",
        "rno" = "Rattus norvegicus",
        "zebrafish" = "Danio rerio",
        "danio rerio" = "Danio rerio",
        "drerio" = "Danio rerio",
        "dre" = "Danio rerio",
        "fly" = "Drosophila melanogaster",
        "drosophila melanogaster" = "Drosophila melanogaster",
        "dmelanogaster" = "Drosophila melanogaster",
        "dmel" = "Drosophila melanogaster",
        "worm" = "Caenorhabditis elegans",
        "caenorhabditis elegans" = "Caenorhabditis elegans",
        "celegans" = "Caenorhabditis elegans",
        "cele" = "Caenorhabditis elegans",
        "yeast" = "Saccharomyces cerevisiae",
        "saccharomyces cerevisiae" = "Saccharomyces cerevisiae",
        "scerevisiae" = "Saccharomyces cerevisiae",
        "sce" = "Saccharomyces cerevisiae",
        "arabidopsis" = "Arabidopsis thaliana",
        "arabidopsis thaliana" = "Arabidopsis thaliana",
        "athaliana" = "Arabidopsis thaliana",
        "ath" = "Arabidopsis thaliana",
        "chicken" = "Gallus gallus",
        "gallus gallus" = "Gallus gallus",
        "pig" = "Sus scrofa",
        "sus scrofa" = "Sus scrofa",
        "cow" = "Bos taurus",
        "bos taurus" = "Bos taurus",
        "cattle" = "Bos taurus",
        "sheep" = "Ovis aries",
        "ovis aries" = "Ovis aries",
        "dog" = "Canis familiaris",
        "canis familiaris" = "Canis familiaris",
        "canis lupus familiaris" = "Canis familiaris",
        "cat" = "Felis catus",
        "felis catus" = "Felis catus",
        "horse" = "Equus caballus",
        "equus caballus" = "Equus caballus",
        "rabbit" = "Oryctolagus cuniculus",
        "oryctolagus cuniculus" = "Oryctolagus cuniculus",
        "chimpanzee" = "Pan troglodytes",
        "pan troglodytes" = "Pan troglodytes",
        "chimp" = "Pan troglodytes",
        "rhesus" = "Macaca mulatta",
        "macaca mulatta" = "Macaca mulatta",
        "rhesus macaque" = "Macaca mulatta"
    )

    # Look up alias
    normalized <- aliases[[organism]]

    # If not found, try to capitalize properly (Genus species)
    if (is.null(normalized)) {
        parts <- strsplit(organism, "\\s+")[[1]]
        if (length(parts) >= 2) {
            # Capitalize first word, lowercase rest
            normalized <- paste(
                paste0(toupper(substr(parts[1], 1, 1)), substr(parts[1], 2, nchar(parts[1]))),
                paste(tolower(parts[-1]), collapse = " ")
            )
        } else {
            normalized <- organism
        }
    }

    normalized
}

#' Get organism short code (3-4 letters)
#'
#' @param organism Character, organism name
#' @return Character, short code like "hsa", "mmu", etc.
#' @export
get_organism_code <- function(organism) {
    organism <- normalize_organism_name(organism)
    if (is.null(organism)) return(NULL)

    codes <- list(
        "Homo sapiens" = "hsa",
        "Mus musculus" = "mmu",
        "Rattus norvegicus" = "rno",
        "Danio rerio" = "dre",
        "Drosophila melanogaster" = "dme",
        "Caenorhabditis elegans" = "cel",
        "Saccharomyces cerevisiae" = "sce",
        "Arabidopsis thaliana" = "ath",
        "Gallus gallus" = "gga",
        "Sus scrofa" = "ssc",
        "Bos taurus" = "bta",
        "Ovis aries" = "oar",
        "Canis familiaris" = "cfa",
        "Pan troglodytes" = "ptr",
        "Macaca mulatta" = "mcc"
    )

    code <- codes[[organism]]

    # Generate code if not in list
    if (is.null(code)) {
        parts <- strsplit(organism, "\\s+")[[1]]
        if (length(parts) >= 2) {
            code <- paste0(tolower(substr(parts[1], 1, 1)),
                           tolower(substr(parts[2], 1, 2)))
        }
    }

    code
}

#' Check if organism is a "model organism" with good annotation
#'
#' @param organism Character, organism name
#' @return Logical
#' @export
is_model_organism <- function(organism) {
    organism <- normalize_organism_name(organism)

    model_organisms <- c(
        "Homo sapiens",
        "Mus musculus",
        "Rattus norvegicus",
        "Danio rerio",
        "Drosophila melanogaster",
        "Caenorhabditis elegans",
        "Saccharomyces cerevisiae",
        "Arabidopsis thaliana"
    )

    organism %in% model_organisms
}
