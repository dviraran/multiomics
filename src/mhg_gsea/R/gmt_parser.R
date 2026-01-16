#' Parse GMT File
#'
#' Reads a Gene Matrix Transposed (GMT) file and converts it into a named list of gene sets.
#'
#' @param file A character string specifying the path to the GMT file.
#' @return A named list where each element is a character vector of gene symbols.
#' @export
parse_gmt <- function(file) {
    if (!file.exists(file)) stop("File not found: ", file)

    # Read all lines
    lines <- readLines(file)

    # Filter empty lines
    lines <- lines[trimws(lines) != ""]

    gene_sets <- lapply(lines, function(line) {
        parts <- strsplit(line, "\t")[[1]]

        # Must have at least name and description (and maybe genes)
        if (length(parts) < 2) {
            return(NULL)
        }

        name <- parts[1]
        # description <- parts[2] # skip description

        # Genes start from index 3 to end
        if (length(parts) >= 3) {
            genes <- parts[3:length(parts)]
            # Remove empty strings if any
            genes <- genes[genes != ""]
        } else {
            genes <- character(0)
        }

        # Return as list first to allow naming later
        return(list(name = name, genes = genes))
    })

    # Remove invalid entries
    gene_sets <- gene_sets[!sapply(gene_sets, is.null)]

    # Create named list
    result <- lapply(gene_sets, function(x) x$genes)
    names(result) <- sapply(gene_sets, function(x) x$name)

    return(result)
}
