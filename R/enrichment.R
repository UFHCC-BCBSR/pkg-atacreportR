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
