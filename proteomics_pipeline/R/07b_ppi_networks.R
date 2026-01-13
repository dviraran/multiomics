# =============================================================================
# Protein-Protein Interaction Network Analysis
# =============================================================================
# PPI network construction and analysis using STRING database and CORUM
#
# Key analyses:
# - STRING database PPI retrieval
# - Network topology analysis (hubs, communities, centrality)
# - Active subnetwork identification
# - Protein complex analysis (CORUM)
# - Network-based enrichment
# - Subcellular localization enrichment
# =============================================================================

# -----------------------------------------------------------------------------
# Main Orchestrator Function
# -----------------------------------------------------------------------------

#' Run PPI Network Analysis
#'
#' Comprehensive protein-protein interaction network analysis
#'
#' @param da_results Differential abundance results data frame
#' @param normalized_data Normalized protein intensity matrix
#' @param config Pipeline configuration list
#' @return List containing PPI network results and plots
#' @export
run_ppi_network_analysis <- function(da_results, normalized_data, config) {
  log_message("=== Running PPI Network Analysis ===")

  # Get config settings
  ppi_config <- config$ppi_analysis %||% list()

  if (!(ppi_config$run_ppi %||% TRUE)) {
    log_message("PPI analysis disabled in config. Skipping.")
    return(NULL)
  }

  # Set up output directories
  output_dir <- file.path(config$output$output_dir, "ppi_networks")
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  # Get significant proteins for network construction
  sig_threshold <- ppi_config$significance_threshold %||% 0.05
  lfc_threshold <- ppi_config$lfc_threshold %||% 0

  # Handle case where da_results might be a list (from limma output)
  if (is.list(da_results) && !is.data.frame(da_results)) {
    # Try to extract the data frame from common structures
    if ("table" %in% names(da_results)) {
      da_results <- da_results$table
    } else if (length(da_results) > 0 && is.data.frame(da_results[[1]])) {
      da_results <- da_results[[1]]
    } else {
      log_message("WARNING: Could not extract data frame from da_results. Skipping PPI analysis.")
      return(NULL)
    }
  }

  # Check for required columns
  if (!is.data.frame(da_results)) {
    log_message("WARNING: da_results is not a data frame. Skipping PPI analysis.")
    return(NULL)
  }

  # Map column names (handle different naming conventions)
  if (!"padj" %in% colnames(da_results) && "adj.P.Val" %in% colnames(da_results)) {
    da_results$padj <- da_results$adj.P.Val
  }
  if (!"log2FoldChange" %in% colnames(da_results) && "logFC" %in% colnames(da_results)) {
    da_results$log2FoldChange <- da_results$logFC
  }
  if (!"protein_id" %in% colnames(da_results)) {
    # Use row names as protein IDs
    da_results$protein_id <- rownames(da_results)
  }

  sig_proteins <- da_results %>%
    dplyr::filter(padj < sig_threshold, abs(log2FoldChange) > lfc_threshold) %>%
    dplyr::pull(protein_id)

  log_message("  Found ", length(sig_proteins), " significant proteins for network analysis")

  if (length(sig_proteins) < 5) {
    log_message("WARNING: Too few significant proteins for network analysis. Skipping.")
    return(NULL)
  }

  # Build PPI network from STRING
  ppi_network <- build_string_network(
    proteins = sig_proteins,
    species = ppi_config$species %||% 9606,  # Human by default
    score_threshold = ppi_config$string_score_threshold %||% 400,
    config = config
  )

  if (is.null(ppi_network) || igraph::ecount(ppi_network$graph) == 0) {
    log_message("WARNING: No PPI edges found. Skipping network analysis.")
    return(NULL)
  }

  # Network topology analysis
  topology_results <- analyze_network_topology(
    graph = ppi_network$graph,
    da_results = da_results,
    output_dir = output_dir
  )

  # Community detection
  community_results <- detect_network_communities(
    graph = ppi_network$graph,
    da_results = da_results,
    output_dir = output_dir
  )

  # Active subnetwork identification
  active_subnet <- NULL
  if (ppi_config$active_subnetwork %||% TRUE) {
    active_subnet <- identify_active_subnetworks(
      graph = ppi_network$graph,
      da_results = da_results,
      config = config,
      output_dir = output_dir
    )
  }

  # Protein complex analysis (CORUM)
  complex_results <- NULL
  if (ppi_config$complex_analysis %||% TRUE) {
    complex_results <- analyze_protein_complexes(
      proteins = sig_proteins,
      da_results = da_results,
      species = ppi_config$species %||% 9606,
      output_dir = output_dir
    )
  }

  # Subcellular localization enrichment
  subcell_results <- analyze_subcellular_localization(
    proteins = sig_proteins,
    da_results = da_results,
    output_dir = output_dir
  )

  # Network-based pathway enrichment
  network_enrichment <- run_network_enrichment(
    graph = ppi_network$graph,
    community_results = community_results,
    config = config,
    output_dir = output_dir
  )

  # Generate visualizations
  figures <- generate_ppi_plots(
    ppi_network = ppi_network,
    topology_results = topology_results,
    community_results = community_results,
    complex_results = complex_results,
    subcell_results = subcell_results,
    da_results = da_results,
    config = config,
    plots_dir = plots_dir
  )

  # Compile results
  results <- list(
    network = ppi_network,
    topology = topology_results,
    communities = community_results,
    active_subnetwork = active_subnet,
    complexes = complex_results,
    subcellular = subcell_results,
    enrichment = network_enrichment,
    figures = figures,
    summary = create_ppi_summary(ppi_network, topology_results, community_results, complex_results)
  )

  # Save outputs
  save_ppi_outputs(results, output_dir)

  log_message("PPI network analysis complete!")
  return(results)
}

# -----------------------------------------------------------------------------
# STRING Network Construction
# -----------------------------------------------------------------------------

#' Build STRING PPI Network
#'
#' @param proteins Vector of protein identifiers
#' @param species NCBI taxonomy ID (9606 for human, 10090 for mouse)
#' @param score_threshold Minimum STRING score (0-1000)
#' @param config Pipeline configuration
#' @return List with igraph object and edge data
build_string_network <- function(proteins, species = 9606, score_threshold = 400, config = NULL) {
  log_message("Building STRING PPI network...")
  log_message("  Species: ", species, " (", get_species_name(species), ")")
  log_message("  Score threshold: ", score_threshold)

  # Check for STRINGdb package
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    log_message("WARNING: STRINGdb package not available. Attempting API fallback...")
    return(build_string_network_api(proteins, species, score_threshold))
  }

  tryCatch({
    # Initialize STRING database connection
    string_db <- STRINGdb::STRINGdb$new(
      version = "12.0",
      species = species,
      score_threshold = score_threshold,
      input_directory = ""
    )

    # Map proteins to STRING IDs
    proteins_df <- data.frame(protein = proteins, stringsAsFactors = FALSE)
    mapped <- string_db$map(proteins_df, "protein", removeUnmappedRows = TRUE)

    log_message("  Mapped ", nrow(mapped), " of ", length(proteins), " proteins to STRING")

    if (nrow(mapped) < 2) {
      log_message("WARNING: Too few proteins mapped to STRING.")
      return(NULL)
    }

    # Get interactions
    interactions <- string_db$get_interactions(mapped$STRING_id)

    log_message("  Found ", nrow(interactions), " interactions")

    if (nrow(interactions) == 0) {
      return(NULL)
    }

    # Create igraph network
    graph <- igraph::graph_from_data_frame(
      interactions[, c("from", "to", "combined_score")],
      directed = FALSE
    )

    # Add protein names as attributes
    string_to_protein <- setNames(mapped$protein, mapped$STRING_id)
    igraph::V(graph)$protein_name <- string_to_protein[igraph::V(graph)$name]

    # Remove self-loops and simplify
    graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)

    return(list(
      graph = graph,
      edges = interactions,
      mapping = mapped,
      n_nodes = igraph::vcount(graph),
      n_edges = igraph::ecount(graph)
    ))

  }, error = function(e) {
    log_message("ERROR in STRING network construction: ", e$message)
    log_message("Attempting API fallback...")
    return(build_string_network_api(proteins, species, score_threshold))
  })
}

#' Build STRING Network via API (Fallback)
build_string_network_api <- function(proteins, species, score_threshold) {
  log_message("  Using STRING API fallback...")

  # Limit to first 500 proteins (API limit)
  if (length(proteins) > 500) {
    log_message("WARNING: Limiting to first 500 proteins for API query")
    proteins <- proteins[1:500]
  }

  # Build API URL
  base_url <- "https://string-db.org/api/tsv/network"
  protein_list <- paste(proteins, collapse = "%0d")

  url <- paste0(
    base_url,
    "?identifiers=", protein_list,
    "&species=", species,
    "&required_score=", score_threshold,
    "&caller_identity=multiomics_pipeline"
  )

  tryCatch({
    # Fetch interactions
    response <- utils::read.delim(url(url), stringsAsFactors = FALSE)

    if (nrow(response) == 0) {
      log_message("  No interactions found via API")
      return(NULL)
    }

    # Create igraph network
    graph <- igraph::graph_from_data_frame(
      response[, c("preferredName_A", "preferredName_B", "score")],
      directed = FALSE
    )

    graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)

    log_message("  API returned ", igraph::ecount(graph), " edges")

    return(list(
      graph = graph,
      edges = response,
      mapping = NULL,
      n_nodes = igraph::vcount(graph),
      n_edges = igraph::ecount(graph)
    ))

  }, error = function(e) {
    log_message("ERROR: STRING API query failed: ", e$message)
    return(NULL)
  })
}

#' Get Species Name from Taxonomy ID
get_species_name <- function(taxid) {
  species_map <- c(
    "9606" = "Homo sapiens",
    "10090" = "Mus musculus",
    "10116" = "Rattus norvegicus",
    "7955" = "Danio rerio",
    "6239" = "Caenorhabditis elegans",
    "7227" = "Drosophila melanogaster"
  )
  species_map[as.character(taxid)] %||% paste0("Species ", taxid)
}

# -----------------------------------------------------------------------------
# Network Topology Analysis
# -----------------------------------------------------------------------------

#' Analyze Network Topology
#'
#' @param graph igraph network object
#' @param da_results DA results for node annotation
#' @param output_dir Output directory
#' @return List with topology metrics
analyze_network_topology <- function(graph, da_results, output_dir) {
  log_message("Analyzing network topology...")

  # Basic metrics
  n_nodes <- igraph::vcount(graph)
  n_edges <- igraph::ecount(graph)
  density <- igraph::edge_density(graph)
  diameter <- igraph::diameter(graph)
  avg_path_length <- igraph::mean_distance(graph)
  clustering_coef <- igraph::transitivity(graph, type = "global")

  log_message("  Nodes: ", n_nodes)
  log_message("  Edges: ", n_edges)
  log_message("  Density: ", round(density, 4))

  # Node centrality metrics
  degree <- igraph::degree(graph)
  betweenness <- igraph::betweenness(graph)
  closeness <- igraph::closeness(graph)
  eigenvector <- tryCatch({
    igraph::eigen_centrality(graph)$vector
  }, error = function(e) rep(NA, n_nodes))

  # Create node metrics data frame
  node_names <- igraph::V(graph)$name
  protein_names <- igraph::V(graph)$protein_name
  if (is.null(protein_names)) protein_names <- node_names

  node_metrics <- data.frame(
    node_id = node_names,
    protein_name = protein_names,
    degree = degree,
    betweenness = betweenness,
    closeness = closeness,
    eigenvector = eigenvector,
    stringsAsFactors = FALSE
  )

  # Add DA results if available
  if (!is.null(da_results)) {
    node_metrics <- merge(
      node_metrics,
      da_results[, c("protein_id", "log2FoldChange", "padj")],
      by.x = "protein_name",
      by.y = "protein_id",
      all.x = TRUE
    )
  }

  # Identify hub proteins (top 10% by degree)
  degree_threshold <- quantile(degree, 0.9)
  hub_proteins <- node_metrics %>%
    filter(degree >= degree_threshold) %>%
    arrange(desc(degree))

  log_message("  Hub proteins (top 10%): ", nrow(hub_proteins))

  # Identify bottleneck proteins (high betweenness)
  betweenness_threshold <- quantile(betweenness, 0.9)
  bottleneck_proteins <- node_metrics %>%
    filter(betweenness >= betweenness_threshold) %>%
    arrange(desc(betweenness))

  # Save node metrics
  write.csv(
    node_metrics,
    file.path(output_dir, "network_node_metrics.csv"),
    row.names = FALSE
  )

  write.csv(
    hub_proteins,
    file.path(output_dir, "hub_proteins.csv"),
    row.names = FALSE
  )

  return(list(
    global_metrics = list(
      n_nodes = n_nodes,
      n_edges = n_edges,
      density = density,
      diameter = diameter,
      avg_path_length = avg_path_length,
      clustering_coefficient = clustering_coef
    ),
    node_metrics = node_metrics,
    hub_proteins = hub_proteins,
    bottleneck_proteins = bottleneck_proteins
  ))
}

# -----------------------------------------------------------------------------
# Community Detection
# -----------------------------------------------------------------------------

#' Detect Network Communities
#'
#' @param graph igraph network object
#' @param da_results DA results data frame
#' @param output_dir Output directory
#' @return List with community detection results
detect_network_communities <- function(graph, da_results, output_dir) {
  log_message("Detecting network communities...")

  # Try multiple community detection algorithms
  communities <- list()

  # Louvain (fast, good for large networks)
  tryCatch({
    communities$louvain <- igraph::cluster_louvain(graph)
    log_message("  Louvain: ", length(unique(igraph::membership(communities$louvain))), " communities")
  }, error = function(e) log_message("  Louvain failed: ", e$message))

  # Walktrap
  tryCatch({
    communities$walktrap <- igraph::cluster_walktrap(graph)
    log_message("  Walktrap: ", length(unique(igraph::membership(communities$walktrap))), " communities")
  }, error = function(e) log_message("  Walktrap failed: ", e$message))

  # Fast greedy (modularity optimization)
  tryCatch({
    communities$fast_greedy <- igraph::cluster_fast_greedy(graph)
    log_message("  Fast greedy: ", length(unique(igraph::membership(communities$fast_greedy))), " communities")
  }, error = function(e) log_message("  Fast greedy failed: ", e$message))

  # Use Louvain as primary (best balance of speed and quality)
  if (!is.null(communities$louvain)) {
    primary_comm <- communities$louvain
  } else if (!is.null(communities$walktrap)) {
    primary_comm <- communities$walktrap
  } else {
    log_message("WARNING: No community detection succeeded")
    return(NULL)
  }

  # Get community membership
  membership <- igraph::membership(primary_comm)
  modularity <- igraph::modularity(primary_comm)

  log_message("  Primary modularity: ", round(modularity, 3))

  # Create community data frame
  node_names <- igraph::V(graph)$name
  protein_names <- igraph::V(graph)$protein_name
  if (is.null(protein_names)) protein_names <- node_names

  community_df <- data.frame(
    node_id = node_names,
    protein_name = protein_names,
    community = as.integer(membership),
    stringsAsFactors = FALSE
  )

  # Add DA results
  if (!is.null(da_results)) {
    community_df <- merge(
      community_df,
      da_results[, c("protein_id", "log2FoldChange", "padj")],
      by.x = "protein_name",
      by.y = "protein_id",
      all.x = TRUE
    )
  }

  # Summarize communities
  community_summary <- community_df %>%
    group_by(community) %>%
    summarise(
      n_proteins = n(),
      mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
      n_upregulated = sum(log2FoldChange > 0 & padj < 0.05, na.rm = TRUE),
      n_downregulated = sum(log2FoldChange < 0 & padj < 0.05, na.rm = TRUE),
      proteins = paste(head(protein_name, 10), collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(desc(n_proteins))

  # Save results
  write.csv(
    community_df,
    file.path(output_dir, "community_membership.csv"),
    row.names = FALSE
  )

  write.csv(
    community_summary,
    file.path(output_dir, "community_summary.csv"),
    row.names = FALSE
  )

  return(list(
    primary = primary_comm,
    all_methods = communities,
    membership = community_df,
    summary = community_summary,
    modularity = modularity,
    n_communities = max(membership)
  ))
}

# -----------------------------------------------------------------------------
# Active Subnetwork Identification
# -----------------------------------------------------------------------------

#' Identify Active Subnetworks
#'
#' @param graph igraph network object
#' @param da_results DA results with significance scores
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return List with active subnetwork results
identify_active_subnetworks <- function(graph, da_results, config, output_dir) {
  log_message("Identifying active subnetworks...")

  # Map node scores from DA results
  node_names <- igraph::V(graph)$protein_name
  if (is.null(node_names)) node_names <- igraph::V(graph)$name

  # Get significance scores (-log10 p-value, signed by direction)
  scores <- sapply(node_names, function(p) {
    idx <- which(da_results$protein_id == p)
    if (length(idx) == 0) return(0)

    pval <- da_results$pvalue[idx[1]]
    lfc <- da_results$log2FoldChange[idx[1]]

    if (is.na(pval) || is.na(lfc)) return(0)

    score <- -log10(pval + 1e-300) * sign(lfc)
    return(score)
  })

  igraph::V(graph)$score <- scores

  # Find connected subnetworks with high aggregate scores
  # Simple approach: identify densely connected high-score regions

  # Threshold for "active" nodes
  score_threshold <- quantile(abs(scores), 0.75, na.rm = TRUE)
  active_nodes <- which(abs(scores) >= score_threshold)

  if (length(active_nodes) < 3) {
    log_message("  Too few active nodes for subnetwork analysis")
    return(NULL)
  }

  # Extract subgraph of active nodes
  active_subgraph <- igraph::induced_subgraph(graph, active_nodes)

  # Find connected components in active subgraph
  components <- igraph::components(active_subgraph)

  # Get largest connected components
  comp_sizes <- table(components$membership)
  large_comps <- as.integer(names(comp_sizes)[comp_sizes >= 3])

  active_subnets <- lapply(large_comps, function(comp_id) {
    comp_nodes <- which(components$membership == comp_id)
    subnet <- igraph::induced_subgraph(active_subgraph, comp_nodes)

    node_names_sub <- igraph::V(subnet)$protein_name
    if (is.null(node_names_sub)) node_names_sub <- igraph::V(subnet)$name

    scores_sub <- igraph::V(subnet)$score

    list(
      graph = subnet,
      proteins = node_names_sub,
      scores = scores_sub,
      mean_score = mean(scores_sub, na.rm = TRUE),
      n_nodes = igraph::vcount(subnet),
      n_edges = igraph::ecount(subnet)
    )
  })

  # Sort by aggregate score
  active_subnets <- active_subnets[order(sapply(active_subnets, function(x) abs(x$mean_score)), decreasing = TRUE)]

  log_message("  Found ", length(active_subnets), " active subnetworks")

  # Save top subnetworks
  subnet_summary <- data.frame(
    subnet_id = seq_along(active_subnets),
    n_nodes = sapply(active_subnets, function(x) x$n_nodes),
    n_edges = sapply(active_subnets, function(x) x$n_edges),
    mean_score = sapply(active_subnets, function(x) x$mean_score),
    proteins = sapply(active_subnets, function(x) paste(head(x$proteins, 10), collapse = ", ")),
    stringsAsFactors = FALSE
  )

  write.csv(
    subnet_summary,
    file.path(output_dir, "active_subnetworks.csv"),
    row.names = FALSE
  )

  return(list(
    subnetworks = active_subnets,
    summary = subnet_summary
  ))
}

# -----------------------------------------------------------------------------
# Protein Complex Analysis (CORUM)
# -----------------------------------------------------------------------------

#' Analyze Protein Complexes
#'
#' @param proteins Vector of protein identifiers
#' @param da_results DA results data frame
#' @param species NCBI taxonomy ID
#' @param output_dir Output directory
#' @return List with complex analysis results
analyze_protein_complexes <- function(proteins, da_results, species, output_dir) {
  log_message("Analyzing protein complexes (CORUM)...")

  # Get CORUM database
  corum_db <- get_corum_database(species)

  if (is.null(corum_db)) {
    log_message("  CORUM database not available. Skipping complex analysis.")
    return(NULL)
  }

  # Find complexes containing our proteins
  proteins_upper <- toupper(proteins)

  complex_hits <- lapply(seq_len(nrow(corum_db)), function(i) {
    complex_proteins <- toupper(unlist(strsplit(corum_db$subunits_gene_name[i], ";")))
    overlap <- intersect(proteins_upper, complex_proteins)

    if (length(overlap) >= 2) {
      return(data.frame(
        complex_id = corum_db$complex_id[i],
        complex_name = corum_db$complex_name[i],
        n_total_subunits = length(complex_proteins),
        n_detected = length(overlap),
        detection_ratio = length(overlap) / length(complex_proteins),
        detected_proteins = paste(overlap, collapse = ";"),
        stringsAsFactors = FALSE
      ))
    }
    return(NULL)
  })

  complex_df <- do.call(rbind, complex_hits[!sapply(complex_hits, is.null)])

  if (is.null(complex_df) || nrow(complex_df) == 0) {
    log_message("  No protein complexes found with >= 2 detected subunits")
    return(NULL)
  }

  # Sort by detection ratio
  complex_df <- complex_df[order(complex_df$detection_ratio, decreasing = TRUE), ]

  log_message("  Found ", nrow(complex_df), " complexes with detected subunits")

  # Test for enrichment of complexes
  # Hypergeometric test for each complex
  total_proteins_in_corum <- length(unique(unlist(strsplit(corum_db$subunits_gene_name, ";"))))
  n_query_proteins <- length(proteins)

  complex_df$pvalue <- sapply(seq_len(nrow(complex_df)), function(i) {
    complex_size <- complex_df$n_total_subunits[i]
    n_overlap <- complex_df$n_detected[i]

    phyper(
      n_overlap - 1,
      complex_size,
      total_proteins_in_corum - complex_size,
      n_query_proteins,
      lower.tail = FALSE
    )
  })

  complex_df$padj <- p.adjust(complex_df$pvalue, method = "BH")

  # Add regulation direction for each complex
  complex_df$mean_log2FC <- sapply(complex_df$detected_proteins, function(prots) {
    prot_list <- unlist(strsplit(prots, ";"))
    matched <- da_results$protein_id %in% prot_list | toupper(da_results$protein_id) %in% prot_list

    if (any(matched)) {
      return(mean(da_results$log2FoldChange[matched], na.rm = TRUE))
    }
    return(NA)
  })

  # Save results
  write.csv(
    complex_df,
    file.path(output_dir, "protein_complex_analysis.csv"),
    row.names = FALSE
  )

  return(list(
    complexes = complex_df,
    n_total = nrow(complex_df),
    n_significant = sum(complex_df$padj < 0.05, na.rm = TRUE)
  ))
}

#' Get CORUM Database
get_corum_database <- function(species) {
  # Try to load from package or download
  corum_file <- system.file("extdata", "corum_core.csv", package = "multiomics")

  if (file.exists(corum_file)) {
    return(read.csv(corum_file, stringsAsFactors = FALSE))
  }

  # Try to download CORUM
  tryCatch({
    log_message("  Downloading CORUM database...")

    url <- "https://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt"
    corum <- read.delim(url(url), stringsAsFactors = FALSE)

    # Filter by species
    if (species == 9606) {
      corum <- corum[corum$Organism == "Human", ]
    } else if (species == 10090) {
      corum <- corum[corum$Organism == "Mouse", ]
    }

    # Standardize column names
    names(corum) <- gsub(" ", "_", tolower(names(corum)))

    # Required columns
    if (!all(c("complexid", "complexname", "subunits(gene_name)") %in% names(corum))) {
      # Try alternative column names
      corum$complex_id <- corum$complexid %||% corum$complex_id
      corum$complex_name <- corum$complexname %||% corum$complex_name
      corum$subunits_gene_name <- corum$`subunits(gene_name)` %||% corum$subunits_gene_name
    } else {
      corum$complex_id <- corum$complexid
      corum$complex_name <- corum$complexname
      corum$subunits_gene_name <- corum$`subunits(gene_name)`
    }

    return(corum)

  }, error = function(e) {
    log_message("WARNING: Could not load CORUM database: ", e$message)
    return(NULL)
  })
}

# -----------------------------------------------------------------------------
# Subcellular Localization Analysis
# -----------------------------------------------------------------------------

#' Analyze Subcellular Localization Enrichment
#'
#' @param proteins Vector of protein identifiers
#' @param da_results DA results data frame
#' @param output_dir Output directory
#' @return List with subcellular localization results
analyze_subcellular_localization <- function(proteins, da_results, output_dir) {
  log_message("Analyzing subcellular localization enrichment...")

  # Check for annotation package
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    log_message("  org.Hs.eg.db not available. Skipping localization analysis.")
    return(NULL)
  }

  tryCatch({
    # Get GO CC (Cellular Component) terms
    go_cc <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = proteins,
      keytype = "SYMBOL",
      columns = c("GOALL", "ONTOLOGYALL")
    )

    # Filter to CC terms only
    go_cc <- go_cc[go_cc$ONTOLOGYALL == "CC" & !is.na(go_cc$GOALL), ]

    if (nrow(go_cc) == 0) {
      log_message("  No GO CC annotations found.")
      return(NULL)
    }

    # Get GO term names
    go_terms <- AnnotationDbi::select(
      GO.db::GO.db,
      keys = unique(go_cc$GOALL),
      keytype = "GOID",
      columns = "TERM"
    )

    # Merge
    go_cc <- merge(go_cc, go_terms, by.x = "GOALL", by.y = "GOID")

    # Count proteins per location
    location_counts <- go_cc %>%
      group_by(TERM) %>%
      summarise(
        n_proteins = n_distinct(SYMBOL),
        proteins = paste(unique(SYMBOL), collapse = "; "),
        .groups = "drop"
      ) %>%
      arrange(desc(n_proteins))

    # Focus on major compartments
    major_locations <- c(
      "nucleus", "cytoplasm", "membrane", "mitochondrion",
      "endoplasmic reticulum", "Golgi", "lysosome", "ribosome",
      "extracellular", "cytoskeleton", "plasma membrane"
    )

    location_summary <- location_counts %>%
      filter(grepl(paste(major_locations, collapse = "|"), TERM, ignore.case = TRUE)) %>%
      head(20)

    # Add DA information
    location_summary$mean_log2FC <- sapply(location_summary$proteins, function(prots) {
      prot_list <- unlist(strsplit(prots, "; "))
      matched <- da_results$protein_id %in% prot_list

      if (any(matched)) {
        return(mean(da_results$log2FoldChange[matched], na.rm = TRUE))
      }
      return(NA)
    })

    # Save results
    write.csv(
      location_summary,
      file.path(output_dir, "subcellular_localization.csv"),
      row.names = FALSE
    )

    return(list(
      locations = location_summary,
      all_annotations = go_cc
    ))

  }, error = function(e) {
    log_message("  Subcellular analysis failed: ", e$message)
    return(NULL)
  })
}

# -----------------------------------------------------------------------------
# Network-Based Enrichment
# -----------------------------------------------------------------------------

#' Run Network-Based Pathway Enrichment
#'
#' @param graph igraph network object
#' @param community_results Community detection results
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return List with enrichment results per community
run_network_enrichment <- function(graph, community_results, config, output_dir) {
  log_message("Running network-based enrichment...")

  if (is.null(community_results)) {
    log_message("  No communities available for enrichment.")
    return(NULL)
  }

  # Check for enrichment packages
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    log_message("  clusterProfiler not available. Skipping enrichment.")
    return(NULL)
  }

  enrichment_results <- list()

  # Enrich each community
  community_df <- community_results$membership
  communities <- unique(community_df$community)

  for (comm in communities) {
    comm_proteins <- community_df$protein_name[community_df$community == comm]

    if (length(comm_proteins) < 3) next

    tryCatch({
      # Convert to Entrez IDs
      entrez <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = comm_proteins,
        keytype = "SYMBOL",
        column = "ENTREZID"
      )
      entrez <- entrez[!is.na(entrez)]

      if (length(entrez) < 3) next

      # GO enrichment
      go_result <- clusterProfiler::enrichGO(
        gene = entrez,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
      )

      if (!is.null(go_result) && nrow(go_result@result) > 0) {
        enrichment_results[[paste0("community_", comm)]] <- list(
          go = go_result@result,
          n_proteins = length(comm_proteins),
          proteins = comm_proteins
        )
      }

    }, error = function(e) {
      # Silently skip failed enrichments
    })
  }

  if (length(enrichment_results) == 0) {
    log_message("  No significant enrichments found in any community.")
    return(NULL)
  }

  log_message("  Enriched ", length(enrichment_results), " communities")

  # Save combined results
  all_enrichments <- lapply(names(enrichment_results), function(comm) {
    df <- enrichment_results[[comm]]$go
    df$community <- comm
    df
  })

  combined <- do.call(rbind, all_enrichments)

  write.csv(
    combined,
    file.path(output_dir, "community_enrichment.csv"),
    row.names = FALSE
  )

  return(enrichment_results)
}

# -----------------------------------------------------------------------------
# Visualization Functions
# -----------------------------------------------------------------------------

#' Generate PPI Network Plots
#'
#' @param ppi_network PPI network object
#' @param topology_results Topology analysis results
#' @param community_results Community detection results
#' @param complex_results Complex analysis results
#' @param subcell_results Subcellular localization results
#' @param da_results DA results
#' @param config Pipeline configuration
#' @param plots_dir Plots output directory
#' @return List of figure paths
generate_ppi_plots <- function(ppi_network, topology_results, community_results,
                                complex_results, subcell_results, da_results,
                                config, plots_dir) {
  log_message("Generating PPI network plots...")

  figures <- list()

  # 1. Main network visualization
  figures$network <- plot_ppi_network(
    graph = ppi_network$graph,
    community_results = community_results,
    da_results = da_results,
    plots_dir = plots_dir
  )

  # 2. Network topology summary
  figures$topology <- plot_network_topology(
    topology_results = topology_results,
    plots_dir = plots_dir
  )

  # 3. Community network
  if (!is.null(community_results)) {
    figures$communities <- plot_community_network(
      graph = ppi_network$graph,
      community_results = community_results,
      plots_dir = plots_dir
    )
  }

  # 4. Hub proteins
  figures$hubs <- plot_hub_proteins(
    topology_results = topology_results,
    da_results = da_results,
    plots_dir = plots_dir
  )

  # 5. Subcellular localization
  if (!is.null(subcell_results)) {
    figures$subcellular <- plot_subcellular_enrichment(
      subcell_results = subcell_results,
      plots_dir = plots_dir
    )
  }

  return(figures)
}

#' Plot PPI Network
plot_ppi_network <- function(graph, community_results, da_results, plots_dir) {
  log_message("  Creating PPI network visualization...")

  # Set node colors based on log2FC
  node_names <- igraph::V(graph)$protein_name
  if (is.null(node_names)) node_names <- igraph::V(graph)$name

  # Get log2FC for each node
  log2fc <- sapply(node_names, function(p) {
    idx <- which(da_results$protein_id == p)
    if (length(idx) > 0) return(da_results$log2FoldChange[idx[1]])
    return(0)
  })

  # Color scale
  col_fun <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("blue", "white", "red")
  )
  node_colors <- col_fun(pmax(pmin(log2fc, 2), -2))

  # Set node size based on degree
  degrees <- igraph::degree(graph)
  node_sizes <- 3 + 5 * (degrees / max(degrees))

  # Layout
  layout <- igraph::layout_with_fr(graph)

  # Plot
  fig_path <- file.path(plots_dir, "ppi_network.png")
  png(fig_path, width = 12, height = 10, units = "in", res = 150)

  par(mar = c(0, 0, 2, 0))
  plot(
    graph,
    layout = layout,
    vertex.color = node_colors,
    vertex.size = node_sizes,
    vertex.label = if (igraph::vcount(graph) <= 50) node_names else NA,
    vertex.label.cex = 0.6,
    vertex.label.color = "black",
    vertex.frame.color = NA,
    edge.width = 0.5,
    edge.color = adjustcolor("grey50", alpha.f = 0.5),
    main = "Protein-Protein Interaction Network"
  )

  # Add legend
  legend(
    "bottomleft",
    legend = c("Down-regulated", "No change", "Up-regulated"),
    fill = c("blue", "white", "red"),
    border = "black",
    title = "log2FC",
    cex = 0.8
  )

  dev.off()

  return(fig_path)
}

#' Plot Network Topology Summary
plot_network_topology <- function(topology_results, plots_dir) {
  log_message("  Creating topology summary plot...")

  node_metrics <- topology_results$node_metrics

  # Degree distribution
  p1 <- ggplot(node_metrics, aes(x = degree)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(title = "Degree Distribution", x = "Degree", y = "Count") +
    theme_minimal()

  # Degree vs betweenness
  p2 <- ggplot(node_metrics, aes(x = degree, y = betweenness)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    labs(title = "Degree vs Betweenness", x = "Degree", y = "Betweenness") +
    theme_minimal()

  # Centrality correlation heatmap
  centrality_cols <- c("degree", "betweenness", "closeness", "eigenvector")
  centrality_mat <- as.matrix(node_metrics[, centrality_cols])
  centrality_cor <- cor(centrality_mat, use = "complete.obs")

  p3 <- ggplot(reshape2::melt(centrality_cor), aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = "Centrality Correlations", x = "", y = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Global metrics text
  gm <- topology_results$global_metrics
  metrics_text <- paste0(
    "Nodes: ", gm$n_nodes, "\n",
    "Edges: ", gm$n_edges, "\n",
    "Density: ", round(gm$density, 4), "\n",
    "Diameter: ", gm$diameter, "\n",
    "Avg Path: ", round(gm$avg_path_length, 2), "\n",
    "Clustering: ", round(gm$clustering_coefficient, 3)
  )

  p4 <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = metrics_text, size = 4, hjust = 0.5) +
    labs(title = "Global Metrics") +
    theme_void()

  # Combine
  combined <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)

  fig_path <- file.path(plots_dir, "network_topology.png")
  ggsave(fig_path, combined, width = 12, height = 10)

  return(fig_path)
}

#' Plot Community Network
plot_community_network <- function(graph, community_results, plots_dir) {
  log_message("  Creating community network visualization...")

  # Get community membership
  membership <- community_results$membership

  # Map to graph nodes
  node_names <- igraph::V(graph)$protein_name
  if (is.null(node_names)) node_names <- igraph::V(graph)$name

  community_map <- setNames(membership$community, membership$protein_name)
  node_communities <- community_map[node_names]

  # Color by community
  n_comm <- max(node_communities, na.rm = TRUE)
  comm_colors <- RColorBrewer::brewer.pal(min(12, max(3, n_comm)), "Set3")
  node_colors <- comm_colors[(node_communities %% length(comm_colors)) + 1]

  # Layout
  layout <- igraph::layout_with_fr(graph)

  # Plot
  fig_path <- file.path(plots_dir, "network_communities.png")
  png(fig_path, width = 12, height = 10, units = "in", res = 150)

  par(mar = c(0, 0, 2, 0))
  plot(
    graph,
    layout = layout,
    vertex.color = node_colors,
    vertex.size = 5,
    vertex.label = NA,
    vertex.frame.color = NA,
    edge.width = 0.3,
    edge.color = adjustcolor("grey50", alpha.f = 0.3),
    main = paste0("Network Communities (n=", n_comm, ", modularity=",
                  round(community_results$modularity, 3), ")")
  )

  # Add legend for top communities
  top_comm <- head(unique(node_communities[!is.na(node_communities)]), 10)
  legend(
    "bottomleft",
    legend = paste0("Community ", top_comm),
    fill = comm_colors[(top_comm %% length(comm_colors)) + 1],
    cex = 0.7,
    title = "Top Communities"
  )

  dev.off()

  return(fig_path)
}

#' Plot Hub Proteins
plot_hub_proteins <- function(topology_results, da_results, plots_dir) {
  log_message("  Creating hub proteins plot...")

  hub_df <- topology_results$hub_proteins

  if (nrow(hub_df) == 0) {
    return(NULL)
  }

  # Top 20 hubs
  hub_df <- head(hub_df, 20)

  # Add color based on log2FC
  hub_df$color <- ifelse(is.na(hub_df$log2FoldChange), "grey",
                         ifelse(hub_df$log2FoldChange > 0, "red", "blue"))

  p <- ggplot(hub_df, aes(x = reorder(protein_name, degree), y = degree, fill = color)) +
    geom_bar(stat = "identity") +
    scale_fill_identity() +
    coord_flip() +
    labs(
      title = "Top Hub Proteins by Degree",
      x = "Protein",
      y = "Degree (Number of Interactions)"
    ) +
    theme_minimal()

  fig_path <- file.path(plots_dir, "hub_proteins.png")
  ggsave(fig_path, p, width = 8, height = 6)

  return(fig_path)
}

#' Plot Subcellular Enrichment
plot_subcellular_enrichment <- function(subcell_results, plots_dir) {
  log_message("  Creating subcellular localization plot...")

  locations <- subcell_results$locations

  if (is.null(locations) || nrow(locations) == 0) {
    return(NULL)
  }

  # Top 15 locations
  locations <- head(locations, 15)

  p <- ggplot(locations, aes(x = reorder(TERM, n_proteins), y = n_proteins, fill = mean_log2FC)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         na.value = "grey50") +
    coord_flip() +
    labs(
      title = "Subcellular Localization of Significant Proteins",
      x = "Cellular Compartment",
      y = "Number of Proteins",
      fill = "Mean log2FC"
    ) +
    theme_minimal()

  fig_path <- file.path(plots_dir, "subcellular_enrichment.png")
  ggsave(fig_path, p, width = 10, height = 8)

  return(fig_path)
}

# -----------------------------------------------------------------------------
# Output Functions
# -----------------------------------------------------------------------------

#' Create PPI Analysis Summary
create_ppi_summary <- function(ppi_network, topology_results, community_results, complex_results) {
  summary <- list(
    n_nodes = ppi_network$n_nodes,
    n_edges = ppi_network$n_edges,
    density = topology_results$global_metrics$density,
    n_hub_proteins = nrow(topology_results$hub_proteins),
    top_hubs = head(topology_results$hub_proteins$protein_name, 5)
  )

  if (!is.null(community_results)) {
    summary$n_communities <- community_results$n_communities
    summary$modularity <- community_results$modularity
  }

  if (!is.null(complex_results)) {
    summary$n_complexes_detected <- complex_results$n_total
    summary$n_significant_complexes <- complex_results$n_significant
  }

  return(summary)
}

#' Save PPI Analysis Outputs
save_ppi_outputs <- function(results, output_dir) {
  log_message("Saving PPI analysis outputs...")

  # Save network as graphml
  if (!is.null(results$network$graph)) {
    igraph::write_graph(
      results$network$graph,
      file.path(output_dir, "ppi_network.graphml"),
      format = "graphml"
    )
  }

  # Save edge list
  if (!is.null(results$network$edges)) {
    write.csv(
      results$network$edges,
      file.path(output_dir, "ppi_network_edges.csv"),
      row.names = FALSE
    )
  }

  # Save full results as RDS
  saveRDS(results, file.path(output_dir, "ppi_analysis_results.rds"))

  log_message("  Outputs saved to: ", output_dir)
}
