#' Prepare gene lists for enrichment analysis
#'
#' Extracts up- and down-regulated gene lists from differential accessibility
#' results and prepares them for pathway enrichment analysis.
#'
#' @param results_list Named list of differential accessibility results from
#'   run_differential_analysis()
#' @param organism Character string: "hsa" for human or "mmu" for mouse
#' @param sig_cutoff FDR significance threshold (default: 0.05)
#' @param logfc_cutoff Absolute log fold-change threshold (default: 1)
#'
#' @return List containing:
#'   \item{gene_lists}{Named list of Entrez ID vectors for up/down genes per contrast}
#'   \item{results_df}{Combined data frame of all results with Contrast column}
#'   \item{universe}{Character vector of all Entrez IDs (background)}
#'   \item{org_db}{OrgDb object for the organism}
#'
#' @importFrom dplyr bind_rows filter pull rename mutate
#' @importFrom AnnotationDbi mapIds
#' @export
prepare_gene_lists_for_enrichment <- function(results_list,
                                              organism,
                                              sig_cutoff = 0.05,
                                              logfc_cutoff = 1) {
  FDR <- logFC <- ENTREZID <- Entrez.ID <- ensembleID <- NULL

  # Load appropriate annotation database
  if (organism == "hsa") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package 'org.Hs.eg.db' is required for human analysis")
    }
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  } else if (organism == "mmu") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Package 'org.Mm.eg.db' is required for mouse analysis")
    }
    org_db <- org.Mm.eg.db::org.Mm.eg.db
  } else {
    stop("Unsupported organism: must be 'hsa' or 'mmu'")
  }

  cat("Loading annotation database for", organism, "\n")

  # Convert all results to named gene lists by direction
  gene_list_named <- list()

  for (contrast in names(results_list)) {
    df <- results_list[[contrast]]$table

    # Add Entrez IDs if not present
    df <- .add_entrez_ids(df, org_db)

    # Extract up-regulated genes
    up_genes <- df %>%
      dplyr::filter(FDR < sig_cutoff, logFC > logfc_cutoff) %>%
      dplyr::pull(ENTREZID) %>%
      na.omit() %>%
      unique() %>%
      as.character()

    # Extract down-regulated genes
    down_genes <- df %>%
      dplyr::filter(FDR < sig_cutoff, logFC < -logfc_cutoff) %>%
      dplyr::pull(ENTREZID) %>%
      na.omit() %>%
      unique() %>%
      as.character()

    # Clean contrast name for list names
    parts <- strsplit(contrast, "_vs_")[[1]]
    contrast_clean <- paste(.clean_group_name(parts[1]), "vs",
                            .clean_group_name(parts[2]), sep = "_")

    gene_list_named[[paste0(contrast_clean, ".up")]] <- up_genes
    gene_list_named[[paste0(contrast_clean, ".down")]] <- down_genes

    cat("  ", contrast_clean, ": ", length(up_genes), " up, ",
        length(down_genes), " down\n", sep = "")
  }

  # Combine all results into one dataframe
  de_results_df <- dplyr::bind_rows(
    lapply(names(results_list), function(contrast) {
      df <- results_list[[contrast]]$table
      df <- .add_entrez_ids(df, org_db)
      df$Contrast <- contrast
      return(df)
    })
  )

  # Create universe (all genes tested)
  universe_entrez <- de_results_df %>%
    dplyr::pull(ENTREZID) %>%
    na.omit() %>%
    unique() %>%
    as.character()

  cat("Universe contains", length(universe_entrez), "unique genes\n")

  return(list(
    gene_lists = gene_list_named,
    results_df = de_results_df,
    universe = universe_entrez,
    org_db = org_db
  ))
}

#' Add Entrez IDs to results dataframe
#' @keywords internal
.add_entrez_ids <- function(df, org_db) {

  Entrez.ID <- ensembleID <- ENTREZID <- NULL


  # Check if already has ENTREZID
  if ("ENTREZID" %in% colnames(df)) {
    return(df)
  }

  # Check if we need to map IDs or if they're already Entrez IDs
  if ("geneId" %in% colnames(df)) {
    # ChIPseeker format - geneId is already Entrez ID
    df$ENTREZID <- as.character(df$geneId)

  } else if ("Entrez.ID" %in% colnames(df)) {
    # Check if Entrez.ID column contains actual Entrez IDs or Ensembl IDs
    sample_id <- df$Entrez.ID[!is.na(df$Entrez.ID)][1]

    if (!is.na(sample_id) && grepl("^ENSMUSG|^ENSG", sample_id)) {
      # These are Ensembl IDs - need to map
      df <- df %>%
        dplyr::rename(ensembleID = Entrez.ID) %>%
        dplyr::mutate(
          ENTREZID = as.character(AnnotationDbi::mapIds(
            org_db,
            keys = ensembleID,
            column = "ENTREZID",
            keytype = "ENSEMBL",
            multiVals = "first"
          ))
        )
    } else {
      # These are already Entrez IDs - use directly
      df$ENTREZID <- as.character(df$Entrez.ID)
    }
  } else {
    stop("Cannot find gene ID column (geneId or Entrez.ID) in results")
  }

  return(df)
}

#' Generate GO enrichment plot for ATAC-seq results
#'
#' Performs Gene Ontology enrichment analysis with custom interactive visualization
#' including gene symbols in hover tooltips.
#'
#' @param gene_lists Named list of Entrez gene ID vectors (from prepare_gene_lists_for_enrichment)
#' @param de_results_df Data frame with differential accessibility results
#' @param universe_entrez Character vector of background gene universe
#' @param org_db OrgDb annotation database object
#' @param ont_category Character string: "BP", "MF", or "CC"
#' @param significance_threshold Numeric p-value cutoff (default 0.05)
#' @param top_n Integer number of top terms to plot (default 10)
#' @param n_cores Number of cores for parallelization (default: auto-detect)
#'
#' @return List containing interactive_plot, static_plot, and go_results data frame
#'
#' @importFrom clusterProfiler compareCluster enrichGO dotplot
#' @importFrom AnnotationDbi mapIds
#' @importFrom dplyr bind_rows filter mutate arrange group_by ungroup distinct slice_head pull case_when select
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient scale_size_manual guides guide_colorbar theme_minimal theme element_text ggtitle xlab ylab annotate theme_void
#' @importFrom plotly ggplotly layout
#' @export
generate_enrichment_plot_atac <- function(gene_lists,
                                          de_results_df,
                                          universe_entrez,
                                          org_db,
                                          ont_category,
                                          significance_threshold = 0.05,
                                          top_n = 10,
                                          n_cores = 1) {

  # Avoid NSE warnings
  Entrez <- Contrast <- Entrez.ID <- ENTREZID <- Cluster <- Description <- NULL
  GeneSymbols <- logFC <- FDR <- geneID <- p.adjust <- GeneRatio <- GeneRatioCategory <- NULL
  plot_label <- tooltip_text <- NULL

  # Auto-detect cores if not specified
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 2)
  }

  # Clean gene list names
  names(gene_lists) <- gsub("efit_|_results_df", "", names(gene_lists))

  # Ensure gene lists are named
  if (is.null(names(gene_lists))) {
    stop("Each gene list must be named!")
  }
  contrast_order <- names(gene_lists)

  # Prepare data for compareCluster
  data <- dplyr::bind_rows(lapply(contrast_order, function(contrast) {
    genes <- gene_lists[[contrast]]
    if (length(genes) == 0 || all(is.na(genes))) {
      return(NULL)
    }
    data.frame(
      Entrez = as.character(genes),
      Contrast = contrast,
      stringsAsFactors = FALSE
    )
  }))

  # Handle no genes case
  if (is.null(data) || nrow(data) == 0) {
    return(.create_empty_enrichment_plot(ont_category))
  }

  # Ensure no duplicates
  data <- data %>%
    dplyr::distinct(Entrez, Contrast, .keep_all = TRUE)

  # Run GO enrichment with parallelization
  cat("Running GO enrichment analysis (", ont_category, ") with ", n_cores, " cores...\n", sep = "")

  if (n_cores > 1) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("Package 'BiocParallel' is required for parallelization")
    }
    BiocParallel::register(BiocParallel::MulticoreParam(workers = n_cores))
  } else {
    # Serial processing
    BiocParallel::register(BiocParallel::SerialParam())
  }

  formula_res <- clusterProfiler::compareCluster(
    Entrez ~ Contrast,
    data = data,
    fun = "enrichGO",
    universe = na.omit(universe_entrez),
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = ont_category,
    pvalueCutoff = significance_threshold
  )

  # Handle no results
  if (is.null(formula_res) || nrow(formula_res@compareClusterResult) == 0) {
    return(.create_empty_enrichment_plot(ont_category))
  }

  # Ensure clusters are ordered correctly
  formula_res@compareClusterResult$Cluster <- factor(
    formula_res@compareClusterResult$Cluster,
    levels = contrast_order
  )

  # Filter by significance
  filtered_results <- subset(
    formula_res@compareClusterResult,
    p.adjust <= significance_threshold
  )

  if (nrow(filtered_results) == 0) {
    return(.create_empty_enrichment_plot(ont_category))
  }

  # Map Entrez IDs to gene symbols (parallelized)
  cat("Mapping gene symbols...\n")
  if (n_cores > 1) {
    filtered_results$GeneSymbols <- parallel::mclapply(seq_len(nrow(filtered_results)), function(i) {
      assign_gene_symbols(filtered_results$geneID[i], as.character(filtered_results$Cluster[i]))
    }, mc.cores = n_cores) %>% unlist()
  } else {
    # Serial processing
    filtered_results$GeneSymbols <- sapply(seq_len(nrow(filtered_results)), function(i) {
      assign_gene_symbols(filtered_results$geneID[i], as.character(filtered_results$Cluster[i]))
    })
  }

  # Save filtered GO results for download
  download_go_results <- filtered_results %>%
    dplyr::select(Cluster, Description, p.adjust, GeneSymbols, dplyr::everything())

  # Identify top GO terms for plotting
  top_GO_terms <- filtered_results %>%
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::pull(Description) %>%
    unique()

  # Filter for plotting
  formula_res@compareClusterResult <- filtered_results %>%
    dplyr::filter(Description %in% top_GO_terms)

  # Convert GeneRatio to numeric
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    dplyr::mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"),
                                     function(x) as.numeric(x[1]) / as.numeric(x[2])))

  # Reorder GO terms using hierarchical clustering
  formula_res@compareClusterResult <- .reorder_GO_terms(formula_res@compareClusterResult)

  # Create size categories for GeneRatio
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(formula_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)
  bin_labels <- c("<=0.01", "0.01 - 0.05", "0.05 - 0.10", ">=0.10")
  size_mapping <- c("<=0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, ">=0.10" = 8)

  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    dplyr::mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks,
                                          labels = bin_labels, include.lowest = TRUE, right = FALSE))

  # Create plot labels and tooltips
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    dplyr::mutate(
      p.adjust = as.numeric(as.character(p.adjust)),
      GeneRatio = as.numeric(as.character(GeneRatio)),
      plot_label = ifelse(
        sapply(strsplit(as.character(Description), " "), length) > 6,
        sapply(strsplit(as.character(Description), " "), function(words) {
          paste(c(words[1:3], "...", tail(words, 3)), collapse = " ")
        }),
        as.character(Description)
      )
    ) %>%
    dplyr::mutate(
      tooltip_text = paste(
        "Cluster: ", as.character(Cluster), "<br>",
        "GO Term: ", as.character(Description), "<br>",
        "p.adjust: ", signif(p.adjust, 3), "<br>",
        "GeneRatio: ", signif(GeneRatio, 3), "<br>",
        "Top Genes:<br>", as.character(GeneSymbols)
      )
    )

  # Create interactive plot
  p <- ggplot2::ggplot(formula_res@compareClusterResult, ggplot2::aes(
    x = Cluster,
    y = plot_label,
    size = GeneRatioCategory,
    color = p.adjust,
    text = tooltip_text
  )) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_size_manual(name = "Gene Ratio", values = size_mapping) +
    ggplot2::scale_color_gradient(low = "red", high = "blue",
                                  limits = c(min(formula_res@compareClusterResult$p.adjust, na.rm = TRUE),
                                             max(formula_res@compareClusterResult$p.adjust, na.rm = TRUE)),
                                  name = "p.adjust") +
    ggplot2::guides(color = ggplot2::guide_colorbar(title = "p.adjust")) +
    ggplot2::ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = "")) +
    ggplot2::xlab("Gene List") +
    ggplot2::ylab("GO Term") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )

  # Convert to interactive
  interactive_plot <- plotly::ggplotly(p, tooltip = "text") %>%
    plotly::layout(legend = list(title = list(text = "Gene Ratio")))

  # Static high-resolution plot
  static_plot <- clusterProfiler::dotplot(formula_res, showCategory = top_n) +
    ggplot2::ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = "")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )

  return(list(
    interactive_plot = interactive_plot,
    static_plot = static_plot,
    go_results = download_go_results
  ))
}

# Helper functions below...
#' @keywords internal
.create_empty_enrichment_plot <- function(ont_category) {
  message_plot <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 1, y = 1,
                      label = paste0("No significant GO enrichment found\n(", ont_category, ")"),
                      size = 6, hjust = 0.5) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = ""))

  return(list(
    interactive_plot = plotly::ggplotly(message_plot),
    static_plot = message_plot,
    go_results = NULL
  ))
}

#' @keywords internal
.add_gene_symbols_parallel <- function(filtered_results, org_db, de_results_df, n_cores) {
  Contrast <- ENTREZID <- FDR <- logFC <- NULL


  # Batch map all Entrez IDs to symbols
  all_entrez_ids <- unique(unlist(strsplit(filtered_results$geneID, "/")))
  gene_symbols_batch <- AnnotationDbi::mapIds(org_db,
                                              keys = all_entrez_ids,
                                              column = "SYMBOL",
                                              keytype = "ENTREZID",
                                              multiVals = "first") %>%
    na.omit()

  gene_symbols_vector <- gene_symbols_batch[all_entrez_ids]

  # Function to assign symbols for one row
  assign_gene_symbols <- function(peak_list, contrast_full) {
    contrast_base <- sub("\\.(up|down)$", "", contrast_full)
    direction <- sub("^.*\\.", "", contrast_full)

    entrez_ids <- unlist(strsplit(peak_list, "/"))

    de_sub <- de_results_df %>%
      dplyr::filter(
        grepl(contrast_base, Contrast),
        ENTREZID %in% entrez_ids,
        dplyr::case_when(
          direction == "up" ~ logFC > 0,
          direction == "down" ~ logFC < 0,
          TRUE ~ FALSE
        )
      ) %>%
      dplyr::arrange(FDR) %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::pull(ENTREZID)

    gene_symbols <- gene_symbols_vector[de_sub]

    if (length(gene_symbols) > 0) {
      paste(unique(gene_symbols), collapse = "<br>")
    } else {
      NA_character_
    }
  }

  # Parallelize
  filtered_results$GeneSymbols <- parallel::mclapply(seq_len(nrow(filtered_results)), function(i) {
    assign_gene_symbols(filtered_results$geneID[i], as.character(filtered_results$Cluster[i]))
  }, mc.cores = n_cores) %>% unlist()

  return(filtered_results)
}

#' @keywords internal
.reorder_GO_terms <- function(df) {
  term_matrix <- table(df$Description, df$Cluster)
  if (nrow(term_matrix) > 1 && length(unique(df$Description)) > 1) {
    term_dist <- dist(term_matrix, method = "binary")
    term_hclust <- hclust(term_dist, method = "ward.D2")
    term_order <- rownames(term_matrix)[term_hclust$order]
  } else {
    term_order <- unique(df$Description)
  }
  df$Description <- factor(df$Description, levels = term_order)
  return(df)
}

#' Generate KEGG enrichment plot for ATAC-seq results
#'
#' Performs KEGG pathway enrichment analysis with custom interactive visualization
#' including gene symbols in hover tooltips.
#'
#' @param gene_lists Named list of Entrez gene ID vectors (from prepare_gene_lists_for_enrichment)
#' @param de_results_df Data frame with differential accessibility results
#' @param universe_entrez Character vector of background gene universe
#' @param org_db OrgDb annotation database object
#' @param organism Character string: "hsa" for human or "mmu" for mouse (for KEGG)
#' @param significance_threshold Numeric p-value cutoff (default 0.05)
#' @param top_n Integer number of top terms to plot (default 10)
#' @param n_cores Number of cores for parallelization (default: 1)
#'
#' @return List containing interactive_plot, static_plot, and kegg_results data frame
#'
#' @importFrom clusterProfiler compareCluster enrichKEGG dotplot
#' @importFrom AnnotationDbi mapIds
#' @importFrom dplyr bind_rows filter mutate arrange group_by ungroup distinct slice_head pull case_when select
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient scale_size_manual guides guide_colorbar theme_minimal theme element_text ggtitle xlab ylab annotate theme_void
#' @importFrom plotly ggplotly layout
#' @export
generate_kegg_enrichment_plot_atac <- function(gene_lists,
                                               de_results_df,
                                               universe_entrez,
                                               org_db,
                                               organism,
                                               significance_threshold = 0.05,
                                               top_n = 10,
                                               n_cores = 1) {

  # Avoid NSE warnings
  Entrez <- Contrast <- ENTREZID <- Cluster <- Description <- NULL
  GeneSymbols <- logFC <- FDR <- geneID <- p.adjust <- GeneRatio <- GeneRatioCategory <- NULL
  plot_label <- tooltip_text <- NULL

  # Clean gene list names
  names(gene_lists) <- gsub("efit_|_results_df", "", names(gene_lists))

  # Ensure gene lists are named
  if (is.null(names(gene_lists))) {
    stop("Each gene list must be named!")
  }
  contrast_order <- names(gene_lists)

  # Prepare data for compareCluster
  data <- dplyr::bind_rows(lapply(contrast_order, function(contrast) {
    genes <- gene_lists[[contrast]]
    if (length(genes) == 0 || all(is.na(genes))) {
      return(NULL)
    }
    data.frame(
      Entrez = as.character(genes),
      Contrast = contrast,
      stringsAsFactors = FALSE
    )
  }))

  # Handle no genes case
  if (is.null(data) || nrow(data) == 0) {
    return(.create_empty_kegg_plot())
  }

  # Ensure no duplicates
  data <- data %>%
    dplyr::distinct(Entrez, Contrast, .keep_all = TRUE)

  # Run KEGG enrichment with parallelization
  cat("Running KEGG enrichment analysis with ", n_cores, " cores...\n", sep = "")

  if (n_cores > 1) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("Package 'BiocParallel' is required for parallelization")
    }
    BiocParallel::register(BiocParallel::MulticoreParam(workers = n_cores))
  } else {
    BiocParallel::register(BiocParallel::SerialParam())
  }

  kegg_res <- clusterProfiler::compareCluster(
    Entrez ~ Contrast,
    data = data,
    fun = "enrichKEGG",
    universe = na.omit(universe_entrez),
    organism = organism,
    keyType = "ncbi-geneid",
    pvalueCutoff = significance_threshold
  )

  # Handle no results
  if (is.null(kegg_res) || nrow(kegg_res@compareClusterResult) == 0) {
    return(.create_empty_kegg_plot())
  }

  # Ensure clusters are ordered correctly
  kegg_res@compareClusterResult$Cluster <- factor(
    kegg_res@compareClusterResult$Cluster,
    levels = contrast_order
  )

  # Filter by significance
  filtered_results <- subset(
    kegg_res@compareClusterResult,
    p.adjust <= significance_threshold
  )

  if (nrow(filtered_results) == 0) {
    return(.create_empty_kegg_plot())
  }

  # Map Entrez IDs to gene symbols (parallelized)
  cat("Mapping gene symbols...\n")
  filtered_results <- .add_gene_symbols_parallel(
    filtered_results,
    org_db,
    de_results_df,
    n_cores
  )

  # Save filtered KEGG results for download
  download_kegg_results <- filtered_results %>%
    dplyr::select(Cluster, Description, p.adjust, GeneSymbols, dplyr::everything())

  # Identify top KEGG pathways for plotting
  top_KEGG_terms <- filtered_results %>%
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::pull(Description) %>%
    unique()

  # Filter for plotting
  kegg_res@compareClusterResult <- filtered_results %>%
    dplyr::filter(Description %in% top_KEGG_terms)

  # Convert GeneRatio to numeric
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    dplyr::mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"),
                                     function(x) as.numeric(x[1]) / as.numeric(x[2])))

  # Reorder KEGG terms using hierarchical clustering
  kegg_res@compareClusterResult <- .reorder_GO_terms(kegg_res@compareClusterResult)

  # Create size categories for GeneRatio
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(kegg_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)
  bin_labels <- c("<=0.01", "0.01 - 0.05", "0.05 - 0.10", ">=0.10")
  size_mapping <- c("<=0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, ">=0.10" = 8)

  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    dplyr::mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks,
                                          labels = bin_labels, include.lowest = TRUE, right = FALSE))

  # Create plot labels and tooltips
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    dplyr::mutate(
      p.adjust = as.numeric(as.character(p.adjust)),
      GeneRatio = as.numeric(as.character(GeneRatio)),
      plot_label = ifelse(
        sapply(strsplit(as.character(Description), " "), length) > 6,
        sapply(strsplit(as.character(Description), " "), function(words) {
          paste(c(words[1:3], "...", tail(words, 3)), collapse = " ")
        }),
        as.character(Description)
      )
    ) %>%
    dplyr::mutate(
      tooltip_text = paste(
        "Cluster: ", as.character(Cluster), "<br>",
        "KEGG Pathway: ", as.character(Description), "<br>",
        "p.adjust: ", signif(p.adjust, 3), "<br>",
        "GeneRatio: ", signif(GeneRatio, 3), "<br>",
        "Top Genes:<br>", as.character(GeneSymbols)
      )
    )

  # Create interactive plot
  p <- ggplot2::ggplot(kegg_res@compareClusterResult, ggplot2::aes(
    x = Cluster,
    y = plot_label,
    size = GeneRatioCategory,
    color = p.adjust,
    text = tooltip_text
  )) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_size_manual(name = "Gene Ratio", values = size_mapping) +
    ggplot2::scale_color_gradient(low = "red", high = "blue",
                                  limits = c(min(kegg_res@compareClusterResult$p.adjust, na.rm = TRUE),
                                             max(kegg_res@compareClusterResult$p.adjust, na.rm = TRUE)),
                                  name = "p.adjust") +
    ggplot2::guides(color = ggplot2::guide_colorbar(title = "p.adjust")) +
    ggplot2::ggtitle("KEGG Pathway Enrichment") +
    ggplot2::xlab("Gene List") +
    ggplot2::ylab("KEGG Pathway") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )

  # Convert to interactive
  interactive_plot <- plotly::ggplotly(p, tooltip = "text") %>%
    plotly::layout(legend = list(title = list(text = "Gene Ratio")))

  # Static high-resolution plot
  static_plot <- clusterProfiler::dotplot(kegg_res, showCategory = top_n) +
    ggplot2::ggtitle("KEGG Pathway Enrichment") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )

  return(list(
    interactive_plot = interactive_plot,
    static_plot = static_plot,
    kegg_results = download_kegg_results
  ))
}

#' @keywords internal
.create_empty_kegg_plot <- function() {
  message_plot <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 1, y = 1,
                      label = "No significant KEGG enrichment found",
                      size = 6, hjust = 0.5) +
    ggplot2::theme_void() +
    ggplot2::ggtitle("KEGG Pathway Enrichment")

  return(list(
    interactive_plot = plotly::ggplotly(message_plot),
    static_plot = message_plot,
    kegg_results = NULL
  ))
}
