#' @importFrom magrittr %>%
#' @importFrom stats na.omit prcomp setNames p.adjust
#' @importFrom utils head read.csv read.delim write.table
#' @importFrom grDevices chull
NULL

#' Parse parameter file for ATAC-seq analysis
#'
#' Reads a text file containing command-line style parameters and converts
#' them to a named list.
#'
#' @param filepath Character string path to parameter file
#'
#' @return Named list of parameters
#'
#' @export
parse_params <- function(filepath) {
  lines <- readLines(filepath)
  param_list <- list()

  for (line in lines) {
    if (grepl("^--", line)) {
      parts <- strsplit(line, " ", fixed = TRUE)[[1]]
      key <- gsub("^--", "", parts[1])
      value <- paste(parts[-1], collapse = " ")
      value <- gsub("^['\"]|['\"]$", "", value)

      if (!is.na(suppressWarnings(as.numeric(value)))) {
        value <- as.numeric(value)
      }
      param_list[[key]] <- value
    }
  }

  return(param_list)
}

#' Summarize ATAC-seq sample QC metrics
#'
#' Generates a comprehensive QC summary table including alignment statistics,
#' peak counts, FRiP scores, and normalization factors.
#'
#' @param dds DESeqDataSet object with sample metadata and QC metrics in colData
#' @param peak_files Optional named vector of per-sample peak files (for calculating
#'   open bases and per-sample peak counts). Names should match sample IDs.
#' @param flagstat_dir Optional directory containing flagstat files for alignment statistics
#' @param render_table Logical; if TRUE, returns a DT::datatable, otherwise a data frame
#'
#' @return Either a DT datatable or data frame with QC metrics
#'
#' @importFrom dplyr left_join mutate select distinct
#' @importFrom purrr map_dfr
#' @importFrom readr read_tsv cols_only col_character col_double
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyselect all_of
#' @importFrom DT datatable
#' @importFrom SummarizedExperiment colData
#' @importFrom stringr str_subset str_extract
#' @export
summarize_atac_sample_qc <- function(dds,
                                     peak_files = NULL,
                                     flagstat_dir = NULL,
                                     render_table = TRUE) {
  # Avoid NSE warnings
  Sample <- Aligned_reads <- Open_bases <- Factor <- Peaks <- FRIP <- NULL

  samples <- colnames(dds)
  coldata <- as.data.frame(SummarizedExperiment::colData(dds))

  # Identify extra metadata columns to include (exclude technical columns)
  exclude_cols <- c("sample", "sizeFactor", "fastq_1", "fastq_2", "replicate",
                    "group.lib.size", "norm.factors", "bigwig_path")
  extra_cols <- setdiff(colnames(coldata), exclude_cols)

  # === 1. SIZE FACTORS ===
  if ("sizeFactor" %in% colnames(coldata)) {
    factor_df <- tibble::tibble(
      Sample = samples,
      Factor = round(coldata[samples, "sizeFactor"], 3)
    )
  } else {
    factor_df <- tibble::tibble(Sample = samples, Factor = NA_real_)
  }

  # === 2. OPEN BASES & PER-SAMPLE PEAK COUNT ===
  if (!is.null(peak_files)) {
    open_bases_df <- purrr::map_dfr(samples, function(sample) {
      if (sample %in% names(peak_files)) {
        peak_file <- peak_files[[sample]]
        if (file.exists(peak_file)) {
          peaks <- readr::read_tsv(peak_file, col_names = FALSE,
                                   col_types = readr::cols_only(
                                     X1 = readr::col_character(),
                                     X2 = readr::col_double(),
                                     X3 = readr::col_double()
                                   ), show_col_types = FALSE)
          total_bases <- sum(peaks$X3 - peaks$X2)
          peak_count <- nrow(peaks)
          return(tibble::tibble(Sample = sample, Open_bases = total_bases, Peaks = peak_count))
        }
      }
      # Fallback: use consensus peak count from DDS
      tibble::tibble(Sample = sample, Open_bases = NA_real_, Peaks = nrow(dds))
    })
  } else {
    # No per-sample peak files - just report consensus count
    open_bases_df <- tibble::tibble(
      Sample = samples,
      Open_bases = NA_real_,
      Peaks = nrow(dds)
    )
  }

  # === 3. ALIGNED READS ===
  aligned_reads_df <- tibble::tibble(Sample = samples, Aligned_reads = NA_integer_)

  if (!is.null(flagstat_dir) && dir.exists(flagstat_dir)) {
    aligned_reads_df <- purrr::map_dfr(samples, function(sample) {
      # Look for flagstat file matching sample name
      flagstat_file <- list.files(flagstat_dir, pattern = paste0(sample, ".*\\.flagstat$"),
                                  full.names = TRUE)

      if (length(flagstat_file) > 0 && file.exists(flagstat_file[1])) {
        flagstat_lines <- readLines(flagstat_file[1])
        aligned_reads <- flagstat_lines %>%
          stringr::str_subset("mapped \\(") %>%
          stringr::str_extract("^\\d+") %>%
          as.numeric()
        return(tibble::tibble(Sample = sample, Aligned_reads = aligned_reads))
      }
      return(tibble::tibble(Sample = sample, Aligned_reads = NA_integer_))
    })
  }

  # === 4. FRIP SCORES ===
  # Extract FRiP from colData if available
  frip_df <- tibble::tibble(Sample = samples, FRIP = NA_real_)

  # Check for FRiP-related columns in colData
  frip_cols <- grep("frip|FRiP|Fraction", colnames(coldata), ignore.case = TRUE, value = TRUE)

  if (length(frip_cols) > 0) {
    # Use the first FRiP column found
    frip_values <- coldata[[frip_cols[1]]]

    # Convert to percentage if it looks like proportion (values < 1)
    if (all(frip_values[!is.na(frip_values)] <= 1)) {
      frip_values <- frip_values * 100
    }

    frip_df$FRIP <- round(frip_values, 1)
  }

  # === 5. JOIN ALL ===
  full_qc_df <- factor_df %>%
    dplyr::left_join(aligned_reads_df, by = "Sample") %>%
    dplyr::left_join(open_bases_df, by = "Sample") %>%
    dplyr::left_join(frip_df, by = "Sample") %>%
    dplyr::select(Sample, Aligned_reads, Open_bases, Factor, Peaks, FRIP)

  # Add extra metadata columns if available
  if (length(extra_cols) > 0) {
    metadata_df <- tibble::as_tibble(coldata[, extra_cols, drop = FALSE], rownames = "Sample")
    full_qc_df <- full_qc_df %>%
      dplyr::left_join(metadata_df, by = "Sample")
  }

  # Remove any duplicate rows
  full_qc_df <- full_qc_df %>% dplyr::distinct()

  # === 6. Optional Render ===
  if (render_table) {
    qc_caption <- "Sample quality control metrics. Some metrics may be unavailable depending on input data."
    return(DT::datatable(
      full_qc_df,
      caption = qc_caption,
      options = list(pageLength = 20, dom = 't'),
      rownames = FALSE
    ))
  }

  return(full_qc_df)
}

#' Generate PCA plot for ATAC-seq data
#'
#' Creates an interactive PCA plot with optional confidence ellipses for groups.
#'
#' @param dge DGEList object with counts
#' @param title Character string for plot title
#' @param show_legend Logical; whether to show legend
#' @param plot_var Character string specifying which variable to color by
#'
#' @return List with plot object and title
#'
#' @importFrom edgeR cpm
#' @importFrom plotly plot_ly add_trace add_polygons layout
#' @importFrom car dataEllipse
#' @export
plot_pca <- function(dge, title, show_legend = TRUE, plot_var) {
  # Avoid NSE warnings
  PC1 <- PC2 <- Group <- Sample <- NULL

  logcpm <- cpm(dge, log = TRUE)

  finite_genes <- apply(logcpm, 1, function(x) all(is.finite(x)))
  logcpm_clean <- logcpm[finite_genes, ]

  pca_res <- prcomp(t(logcpm_clean), center = TRUE, scale. = FALSE)

  scores <- as.data.frame(pca_res$x[, 1:2])
  percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)

  samples_df <- as.data.frame(dge$samples)
  scores$Sample <- rownames(scores)
  scores$Group <- samples_df[[plot_var]]

  fig <- plotly::plot_ly()

  # Add ellipses
  for (group in unique(scores$Group)) {
    group_data <- scores[scores$Group == group, ]
    if (nrow(group_data) >= 2) {
      if (all(is.finite(group_data$PC1)) && all(is.finite(group_data$PC2))) {
        tryCatch({
          if (nrow(group_data) >= 3) {
            ellipse_coords <- car::dataEllipse(group_data$PC1, group_data$PC2,
                                               levels = 0.68, plot.points = FALSE, draw = FALSE)
          } else {
            hull_indices <- chull(group_data$PC1, group_data$PC2)
            ellipse_coords <- as.matrix(group_data[hull_indices, c("PC1", "PC2")])
            center_x <- mean(ellipse_coords[,1])
            center_y <- mean(ellipse_coords[,2])
            ellipse_coords[,1] <- center_x + (ellipse_coords[,1] - center_x) * 1.5
            ellipse_coords[,2] <- center_y + (ellipse_coords[,2] - center_y) * 1.5
          }

          fig <- fig %>%
            plotly::add_polygons(x = ellipse_coords[,1], y = ellipse_coords[,2],
                                 name = paste(group, "CI"),
                                 opacity = 0.2,
                                 showlegend = FALSE,
                                 hoverinfo = "skip")
        }, error = function(e) {
          message("Skipping ellipse for group ", group, ": ", e$message)
        })
      }
    }
  }

  fig <- fig %>%
    plotly::add_trace(
      data = scores,
      x = ~PC1, y = ~PC2,
      color = ~Group,
      text = ~Sample,
      type = 'scatter', mode = 'markers',
      marker = list(size = 10),
      hovertemplate = "%{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>",
      showlegend = show_legend
    ) %>%
    plotly::layout(
      xaxis = list(title = paste0("PC1 (", percentVar[1], "%)")),
      yaxis = list(title = paste0("PC2 (", percentVar[2], "%)")),
      width = 500,
      height = 400
    )

  return(list(plot = fig, title = title))
}

#' Clean group names by removing leading X prefix
#'
#' @param name Character string to clean
#'
#' @return Cleaned character string
#'
#' @export
clean_group_name <- function(name) {
  if (grepl("^X\\d", name)) {
    return(substring(name, 2))
  }
  return(name)
}

#' Assign gene symbols to peak lists
#'
#' Internal helper function for enrichment analysis
#'
#' @param peak_list Character string of Entrez IDs
#' @param contrast_full Character string contrast name
#'
#' @return Character string of gene symbols
#'
#' @keywords internal
assign_gene_symbols <- function(peak_list, contrast_full) {
    # Avoid NSE warnings
    Contrast <- ENTREZID <- logFC <- FDR <- NULL
  # This function needs de_results_df and gene_symbols_vector from parent environment
  # Extract base contrast and direction
  contrast_base <- sub("\\.(up|down)$", "", contrast_full)
  direction <- sub("^.*\\.", "", contrast_full)

  # Split the peak list into Entrez IDs
  entrez_ids <- unlist(strsplit(peak_list, "/"))

  # Access parent environment variables
  de_results_df_filtered <- get("de_results_df_filtered", envir = parent.frame())
  gene_symbols_vector <- get("gene_symbols_vector", envir = parent.frame())

  # Filter de_results_df based on contrast and direction
  de_sub <- de_results_df_filtered %>%
    dplyr::filter(grepl(contrast_base, Contrast),
                  ENTREZID %in% entrez_ids,
                  dplyr::case_when(
                    direction == "up" ~ logFC > 0,
                    direction == "down" ~ logFC < 0
                  ))

  # Sort by FDR and get the top 20 peaks
  top_peaks <- de_sub %>%
    dplyr::arrange(FDR) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::pull(ENTREZID)

  top_peaks <- unlist(top_peaks)

  # Map Entrez IDs to gene symbols
  gene_symbols <- gene_symbols_vector[top_peaks]

  # Return gene symbols as a string (with line breaks)
  if (length(gene_symbols) > 0) {
    paste(unique(gene_symbols), collapse = "<br>")
  } else {
    NA_character_
  }
}

#' Generate GO enrichment plot for ATAC-seq results
#'
#' Performs Gene Ontology enrichment analysis with custom interactive visualization
#' including gene symbols in hover tooltips.
#'
#' @param gene_lists Named list of Entrez gene ID vectors
#' @param de_results_df Data frame with differential accessibility results
#' @param universe_entrez Character vector of background gene universe
#' @param ont_category Character string: "BP", "MF", or "CC"
#' @param significance_threshold Numeric p-value cutoff (default 0.05)
#' @param top_n Integer number of top terms to plot (default 10)
#' @param annotation_db Character string name of annotation database
#'
#' @return List containing interactive_plot, static_plot, and go_results data frame
#'
#' @importFrom clusterProfiler compareCluster enrichGO dotplot
#' @importFrom AnnotationDbi mapIds
#' @importFrom dplyr bind_rows filter mutate arrange group_by ungroup distinct slice_head pull case_when select
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient scale_size_manual guides guide_colorbar theme_minimal theme element_text ggtitle xlab ylab annotate theme_void
#' @importFrom plotly ggplotly layout
#' @export
generate_enrichment_plot_atac <- function(gene_lists, de_results_df, universe_entrez, ont_category, significance_threshold = 0.05, top_n = 10, annotation_db) {
  # Avoid NSE warnings
  Entrez <- Contrast <- Entrez.ID <- ENTREZID <- Cluster <- Description <- NULL
  GeneSymbols <- logFC <- FDR <- geneID <- p.adjust <- GeneRatio <- GeneRatioCategory <- NULL
  plot_label <- tooltip_text <- NULL

  # Load the annotation object
  annotation_obj <- get(annotation_db, envir = asNamespace(annotation_db))
  names(gene_lists) <- gsub("efit_|_results_df","",names(gene_lists))

  # Ensure gene lists are named and define contrast order
  if (is.null(names(gene_lists))) stop("Each gene list must be named!")
  contrast_order <- names(gene_lists)

  # FIX: Handle both Ensembl and Entrez ID formats
  sample_id <- de_results_df$Entrez.ID[!is.na(de_results_df$Entrez.ID)][1]
  if (grepl("^ENSMUSG|^ENSG", sample_id)) {
    cat("Mapping Ensembl IDs to Entrez IDs...\n")
    de_results_df <- de_results_df %>%
      dplyr::mutate(ENTREZID = mapIds(
        annotation_obj,
        keys = Entrez.ID,
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "first"
      ))
  } else {
    cat("Using existing Entrez IDs...\n")
    de_results_df$ENTREZID <- as.character(de_results_df$Entrez.ID)
  }

  # Prepare data for compareCluster
  data <- dplyr::bind_rows(lapply(contrast_order, function(contrast) {
    genes <- gene_lists[[contrast]]
    if (length(genes) == 0 || all(is.na(genes))) {
      return(NULL)
    }
    data.frame(
      Entrez = genes,
      Contrast = contrast,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(data) || nrow(data) == 0) {
    message_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = paste0("No significant genes found for enrichment\n(", ont_category, ")"), size = 6, hjust = 0.5) +
      ggplot2::theme_void() +
      ggplot2::ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = ""))

    interactive_plot <- plotly::ggplotly(message_plot)
    static_plot <- message_plot

    return(list(
      interactive_plot = interactive_plot,
      static_plot = static_plot,
      go_results = NULL
    ))
  }

  # Ensure no duplicates and valid formatting
  data <- data %>%
    dplyr::distinct(Entrez, Contrast, .keep_all = TRUE)
  data$Entrez <- as.character(data$Entrez)

  # PARALLELIZE: Run GO enrichment with parallel processing
  cat("Running GO enrichment analysis (parallelized)...\n")
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("Package 'BiocParallel' is required for parallelization")
  }
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 14))

  formula_res <- clusterProfiler::compareCluster(
    Entrez ~ Contrast,
    data = data,
    fun = "enrichGO",
    universe = na.omit(universe_entrez),
    OrgDb = annotation_obj,
    keyType = "ENTREZID",
    ont = ont_category,
    pvalueCutoff = significance_threshold
  )

  # Handle no results case
  if (is.null(formula_res) || nrow(formula_res@compareClusterResult) == 0) {
    formula_res <- methods::new("compareClusterResult",
                                compareClusterResult = data.frame(
                                  Cluster = factor(),
                                  ID = character(),
                                  Description = character(),
                                  GeneRatio = character(),
                                  BgRatio = character(),
                                  pvalue = numeric(),
                                  p.adjust = numeric(),
                                  qvalue = numeric(),
                                  geneID = character(),
                                  Count = integer(),
                                  stringsAsFactors = FALSE
                                ))
  }

  # Ensure clusters are ordered correctly
  formula_res@compareClusterResult$Cluster <- factor(
    formula_res@compareClusterResult$Cluster,
    levels = contrast_order
  )

  # Filter by significance threshold
  filtered_results <- subset(
    formula_res@compareClusterResult,
    p.adjust <= significance_threshold
  )

  # Handle the special case: no enrichment found
  if (nrow(filtered_results) == 0) {
    message_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = paste0("No significant GO enrichment found\n(", ont_category, ")"), size = 6, hjust = 0.5) +
      ggplot2::theme_void() +
      ggplot2::ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = ""))

    interactive_plot <- plotly::ggplotly(message_plot)
    static_plot <- message_plot

    return(list(
      interactive_plot = interactive_plot,
      static_plot = static_plot,
      go_results = NULL
    ))
  }

  # Step 1: Gather all Entrez IDs from filtered_results
  all_entrez_ids <- unique(unlist(strsplit(filtered_results$geneID, "/")))

  # Step 2: Map all Entrez IDs to gene symbols in one go
  gene_symbols_batch <- mapIds(annotation_obj,
                               keys = all_entrez_ids,
                               column = "SYMBOL",
                               keytype = "ENTREZID",
                               multiVals = "first") %>%
    na.omit()
  gene_symbols_vector <- gene_symbols_batch[all_entrez_ids]

  # PARALLELIZE: Gene symbol assignment
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for parallelization")
  }
  num_cores <- min(6, parallel::detectCores() - 2)

  # Create a function to assign gene symbols to a list of Entrez IDs
  assign_gene_symbols <- function(peak_list, contrast_full) {
    # Extract base contrast and direction
    contrast_base <- sub("\\.(up|down)$", "", contrast_full)
    direction <- sub("^.*\\.", "", contrast_full)

    # Split the peak list into Entrez IDs
    entrez_ids <- unlist(strsplit(peak_list, "/"))

    # Pre-filter de_results_df once for relevant contrasts
    de_results_df_filtered <- de_results_df %>%
      dplyr::filter(!is.na(ENTREZID) & ENTREZID != "")

    # Filter de_results_df based on contrast and direction
    de_sub <- de_results_df_filtered %>%
      dplyr::filter(grepl(contrast_base, Contrast),
                    ENTREZID %in% entrez_ids,
                    dplyr::case_when(
                      direction == "up" ~ logFC > 0,
                      direction == "down" ~ logFC < 0
                    ))

    # Sort by FDR and get the top 20 peaks
    top_peaks <- de_sub %>%
      dplyr::arrange(FDR) %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::pull(ENTREZID)

    top_peaks <- unlist(top_peaks)

    # Map Entrez IDs to gene symbols
    gene_symbols <- gene_symbols_vector[top_peaks]

    # Return gene symbols as a string (with line breaks)
    if (length(gene_symbols) > 0) {
      paste(unique(gene_symbols), collapse = "<br>")
    } else {
      NA_character_
    }
  }

  filtered_results$GeneSymbols <- parallel::mclapply(seq_len(nrow(filtered_results)), function(i) {
    peak_list <- filtered_results$geneID[i]
    contrast_full <- as.character(filtered_results$Cluster[i])
    assign_gene_symbols(peak_list, contrast_full)
  }, mc.cores = num_cores)

  # Convert list to character vector
  filtered_results$GeneSymbols <- unlist(filtered_results$GeneSymbols)

  # Save filtered GO results for download
  download_go_results <- filtered_results %>%
    dplyr::select(Cluster, Description, p.adjust, GeneSymbols, dplyr::everything())

  # Identify top `n` GO terms for plotting
  top_GO_terms <- filtered_results %>%
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::pull(Description) %>%
    unique()

  # Filter for plotting (only top GO terms)
  formula_res@compareClusterResult <- filtered_results %>%
    dplyr::filter(Description %in% top_GO_terms)

  # Convert GeneRatio to numeric
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    dplyr::mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))

  # Reorder GO terms using hierarchical clustering
  reorder_GO_terms <- function(df) {
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

  formula_res@compareClusterResult <- reorder_GO_terms(formula_res@compareClusterResult)

  # Define GeneRatio bins and size mapping
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(formula_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)
  bin_labels <- c("<=0.01", "0.01 - 0.05", "0.05 - 0.10", ">=0.10")

  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    dplyr::mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE))

  size_mapping <- c("<=0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, ">=0.10" = 8)

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
    ggplot2::xlab("DE Gene list") +
    ggplot2::ylab("GO Term") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )

  # Convert to interactive plot
  interactive_plot <- plotly::ggplotly(p, tooltip = "text") %>%
    plotly::layout(legend = list(title = list(text = "Gene Ratio")))

  # Static High-Resolution Plot
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

  # Return results
  return(list(
    interactive_plot = interactive_plot,
    static_plot = static_plot,
    go_results = download_go_results
  ))
}

#' Generate KEGG enrichment plot for ATAC-seq results
#'
#' Performs KEGG pathway enrichment analysis with custom interactive visualization
#' including gene symbols in hover tooltips.
#'
#' @param gene_lists Named list of Entrez gene ID vectors
#' @param de_results_df Data frame with differential accessibility results
#' @param universe_entrez Character vector of background gene universe
#' @param significance_threshold Numeric p-value cutoff (default 0.05)
#' @param top_n Integer number of top terms to plot (default 10)
#' @param annotation_db Character string name of annotation database
#'
#' @return List containing interactive_plot, static_plot, and kegg_results data frame
#'
#' @importFrom clusterProfiler compareCluster enrichKEGG dotplot
#' @importFrom AnnotationDbi mapIds
#' @importFrom dplyr bind_rows filter mutate arrange group_by ungroup distinct slice_head pull case_when select
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient scale_size_manual guides guide_colorbar theme_minimal theme element_text ggtitle xlab ylab annotate theme_void
#' @importFrom plotly ggplotly layout
#' @export
generate_kegg_enrichment_plot_atac <- function(gene_lists, de_results_df, universe_entrez, significance_threshold = 0.05, top_n = 10, annotation_db) {
  # Avoid NSE warnings
  Entrez <- Contrast <- Entrez.ID <- ENTREZID <- Cluster <- Description <- NULL
  GeneSymbols <- logFC <- FDR <- geneID <- p.adjust <- GeneRatio <- GeneRatioCategory <- NULL
  plot_label <- tooltip_text <- NULL

  # Load the annotation object
  annotation_obj <- get(annotation_db, envir = asNamespace(annotation_db))
  names(gene_lists) <- gsub("efit_|_results_df","",names(gene_lists))

  # Ensure gene lists are named and define contrast order
  if (is.null(names(gene_lists))) stop("Each gene list must be named!")
  contrast_order <- names(gene_lists)

  force(de_results_df)

  # FIX: Handle both Ensembl and Entrez ID formats
  sample_id <- de_results_df$Entrez.ID[!is.na(de_results_df$Entrez.ID)][1]
  if (grepl("^ENSMUSG|^ENSG", sample_id)) {
    cat("Mapping Ensembl IDs to Entrez IDs...\n")
    de_results_df <- de_results_df %>%
      dplyr::mutate(ENTREZID = mapIds(
        annotation_obj,
        keys = Entrez.ID,
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "first"
      ))
  } else {
    cat("Using existing Entrez IDs...\n")
    de_results_df$ENTREZID <- as.character(de_results_df$Entrez.ID)
  }

  # Prepare data for compareCluster
  data <- dplyr::bind_rows(lapply(contrast_order, function(contrast) {
    genes <- gene_lists[[contrast]]
    if (length(genes) == 0 || all(is.na(genes))) {
      return(NULL)
    }
    data.frame(
      Entrez = genes,
      Contrast = contrast,
      stringsAsFactors = FALSE
    )
  }))

  # Ensure no duplicates and valid formatting
  data <- data %>%
    dplyr::distinct(Entrez, Contrast, .keep_all = TRUE)
  data$Entrez <- as.character(data$Entrez)

  # PARALLELIZE: Run KEGG enrichment with parallel processing
  cat("Running KEGG enrichment analysis (parallelized)...\n")
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("Package 'BiocParallel' is required for parallelization")
  }
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 14))

  # Get organism from parent environment
  report_params <- get("report_params", envir = parent.frame())

  kegg_res <- clusterProfiler::compareCluster(
    Entrez ~ Contrast,
    data = data,
    fun = "enrichKEGG",
    universe = na.omit(universe_entrez),
    organism = report_params[["organism"]],
    keyType = "ncbi-geneid",
    pvalueCutoff = significance_threshold
  )

  # Handle no results case
  if (is.null(kegg_res) || nrow(kegg_res@compareClusterResult) == 0) {
    kegg_res <- methods::new("compareClusterResult",
                             compareClusterResult = data.frame(
                               Cluster = factor(),
                               ID = character(),
                               Description = character(),
                               GeneRatio = character(),
                               BgRatio = character(),
                               pvalue = numeric(),
                               p.adjust = numeric(),
                               qvalue = numeric(),
                               geneID = character(),
                               Count = integer(),
                               stringsAsFactors = FALSE
                             ))
  }

  # Ensure clusters are ordered correctly
  kegg_res@compareClusterResult$Cluster <- factor(
    kegg_res@compareClusterResult$Cluster,
    levels = contrast_order
  )

  # Filter by significance threshold
  filtered_results <- subset(
    kegg_res@compareClusterResult,
    p.adjust <= significance_threshold
  )

  # Handle the special case: no enrichment found
  if (nrow(filtered_results) == 0) {
    message_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = "No significant KEGG enrichment found", size = 6, hjust = 0.5) +
      ggplot2::theme_void() +
      ggplot2::ggtitle("KEGG Pathway Enrichment")

    interactive_plot <- plotly::ggplotly(message_plot)
    static_plot <- message_plot

    return(list(
      interactive_plot = interactive_plot,
      static_plot = static_plot,
      kegg_results = NULL
    ))
  }

  # Step 1: Gather all Entrez IDs from filtered_results
  all_entrez_ids <- unique(unlist(strsplit(filtered_results$geneID, "/")))

  # Step 2: Map all Entrez IDs to gene symbols in one go
  gene_symbols_batch <- mapIds(annotation_obj,
                               keys = all_entrez_ids,
                               column = "SYMBOL",
                               keytype = "ENTREZID",
                               multiVals = "first") %>%
    na.omit()
  gene_symbols_vector <- gene_symbols_batch[all_entrez_ids]

  # Pre-filter de_results_df once for relevant contrasts
  de_results_df_filtered <- de_results_df %>%
    dplyr::filter(!is.na(ENTREZID) & ENTREZID != "")

  # Create a function to assign gene symbols to a list of Entrez IDs
  assign_gene_symbols <- function(peak_list, contrast_full) {
    # Extract base contrast and direction
    contrast_base <- sub("\\.(up|down)$", "", contrast_full)
    direction <- sub("^.*\\.", "", contrast_full)

    # Split the peak list into Entrez IDs
    entrez_ids <- unlist(strsplit(peak_list, "/"))

    # Filter de_results_df based on contrast and direction
    de_sub <- de_results_df_filtered %>%
      dplyr::filter(grepl(contrast_base, Contrast),
                    ENTREZID %in% entrez_ids,
                    dplyr::case_when(
                      direction == "up" ~ logFC > 0,
                      direction == "down" ~ logFC < 0
                    ))

    # Sort by FDR and get the top 20 peaks
    top_peaks <- de_sub %>%
      dplyr::arrange(FDR) %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::pull(ENTREZID)

    top_peaks <- unlist(top_peaks)

    # Map Entrez IDs to gene symbols
    gene_symbols <- gene_symbols_vector[top_peaks]

    # Return gene symbols as a string (with line breaks)
    if (length(gene_symbols) > 0) {
      paste(unique(gene_symbols), collapse = "<br>")
    } else {
      NA_character_
    }
  }

  # PARALLELIZE: Gene symbol assignment
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for parallelization")
  }
  num_cores <- min(6, parallel::detectCores() - 2)

  filtered_results$GeneSymbols <- parallel::mclapply(seq_len(nrow(filtered_results)), function(i) {
    peak_list <- filtered_results$geneID[i]
    contrast_full <- as.character(filtered_results$Cluster[i])
    assign_gene_symbols(peak_list, contrast_full)
  }, mc.cores = num_cores)

  # Convert list to character vector
  filtered_results$GeneSymbols <- unlist(filtered_results$GeneSymbols)

  # Save filtered KEGG results for download
  download_kegg_results <- filtered_results %>%
    dplyr::select(Cluster, Description, p.adjust, GeneSymbols, dplyr::everything())

  # Identify top `n` KEGG pathways for plotting
  top_KEGG_terms <- filtered_results %>%
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::pull(Description) %>%
    unique()

  # Filter for plotting (only top KEGG pathways)
  kegg_res@compareClusterResult <- filtered_results %>%
    dplyr::filter(Description %in% top_KEGG_terms)

  # Convert GeneRatio to numeric
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    dplyr::mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))

  # Reorder KEGG terms using hierarchical clustering
  reorder_KEGG_terms <- function(df) {
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

  kegg_res@compareClusterResult <- reorder_KEGG_terms(kegg_res@compareClusterResult)

  # Define GeneRatio bins and size mapping
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(kegg_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)
  bin_labels <- c("<=0.01", "0.01 - 0.05", "0.05 - 0.10", ">=0.10")

  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    dplyr::mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE))

  size_mapping <- c("<=0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, ">=0.10" = 8)

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
    ggplot2::xlab("DE Gene list") +
    ggplot2::ylab("KEGG Pathway") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )

  # Convert to interactive plot
  interactive_plot <- plotly::ggplotly(p, tooltip = "text") %>%
    plotly::layout(legend = list(title = list(text = "Gene Ratio")))

  # Static High-Resolution Plot
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

  # Return results
  return(list(
    interactive_plot = interactive_plot,
    static_plot = static_plot,
    kegg_results = download_kegg_results
  ))
}

#' Create downloadable PNG button for plots
#'
#' @param plot_object ggplot2 plot object
#' @param output_name Character string for filename
#' @param width Numeric plot width in inches
#' @param height Numeric plot height in inches
#' @param dpi Numeric resolution
#'
#' @return HTML download button
#'
#' @importFrom ggplot2 ggsave
#' @importFrom base64enc dataURI
#' @importFrom htmltools HTML
#' @export
download_button_png <- function(plot_object, output_name = "plot",
                                width = 12, height = 6, dpi = 300) {
  temp_file <- tempfile(fileext = ".png")
  ggplot2::ggsave(temp_file, plot = plot_object, width = width, height = height, dpi = dpi)

  encoded_img <- base64enc::dataURI(file = temp_file, mime = "image/png")

  download_button <- htmltools::HTML(paste0(
    '<a href="', encoded_img, '" download="', output_name, '.png" ',
    'class="btn btn-primary" style="padding:10px; font-size:16px; text-decoration:none; ',
    'color:white; background-color:#007BFF; border-radius:5px;">',
    'Download ', output_name, '</a>'
  ))

  return(download_button)
}

#' Copy file to web directory
#'
#' @param source_path Source file path
#' @param dest_dir Destination directory
#' @param filename Destination filename
#'
#' @return Logical indicating success
#'
#' @export
copy_file_to_web <- function(source_path, dest_dir, filename) {
  if (file.exists(source_path)) {
    dest_path <- file.path(dest_dir, filename)
    file.copy(source_path, dest_path, overwrite = TRUE)
    system2("chmod", args = c("644", dest_path))
    return(TRUE)
  }
  return(FALSE)
}

#' Create BigBed file for differential accessibility peaks
#'
#' Filters significant DA peaks and creates a BigBed track file for genome browser visualization.
#' Includes directional peak names indicating which group has increased accessibility.
#'
#' @param results_df Data frame with DA results
#' @param contrast_name Character string contrast name
#' @param output_dir Output directory path
#' @param organism Character string organism code ("mmu" or "hsa")
#'
#' @return Character string BigBed filename or NULL if creation failed
#'
#' @importFrom dplyr filter mutate select
#' @export
create_da_bigbed <- function(results_df, contrast_name, output_dir, organism) {
  # Avoid NSE warnings
  FDR <- logFC <- logCPM <- Gene.Name <- Chr <- Start <- End <- Strand <- NULL
  base_name <- direction_suffix <- peak_name <- score <- strand <- interval <- NULL
  . <- NULL  # For dplyr's .[[column]] syntax

  # This function needs report_params from parent environment
  report_params <- get("report_params", envir = parent.frame())

  sig_peaks <- results_df %>%
    dplyr::filter(FDR < report_params$sig_cutoff,
                  abs(logFC) > report_params$logFC_cutoff,
                  logCPM > 2,
                  !is.na(Gene.Name))

  if (nrow(sig_peaks) == 0) {
    cat("No significant peaks found for", contrast_name, "\n")
    return(NULL)
  }

  group_names <- strsplit(contrast_name, "_vs_")[[1]]
  group1 <- group_names[1]
  group2 <- if (length(group_names) > 1) group_names[2] else "baseline"

  chr_col <- if ("Chr" %in% colnames(sig_peaks)) "Chr" else "seqnames"
  start_col <- if ("Start" %in% colnames(sig_peaks)) "Start" else "start"
  end_col <- if ("End" %in% colnames(sig_peaks)) "End" else "end"
  strand_col <- if ("Strand" %in% colnames(sig_peaks)) "Strand" else "strand"

  cat("Detected columns for", contrast_name, ":\n")
  cat("  Chr:", chr_col, "Start:", start_col, "End:", end_col, "Strand:", strand_col, "\n")

  standard_chroms <- if (organism == "mmu") {
    paste0("chr", c(1:19, "X", "Y", "M"))
  } else {
    paste0("chr", c(1:22, "X", "Y", "M"))
  }

  bed_data <- sig_peaks %>%
    dplyr::mutate(
      Chr = as.character(.[[chr_col]]),
      Start = .[[start_col]],
      End = .[[end_col]],
      Strand = if (strand_col %in% colnames(.)) .[[strand_col]] else ".",
      base_name = ifelse(!is.na(Gene.Name) & Gene.Name != "", Gene.Name, interval),
      direction_suffix = ifelse(logFC > report_params$logFC_cutoff,
                                paste0("_MORE_ACC_IN_", toupper(group1)),
                                paste0("_MORE_ACC_IN_", toupper(group2))),
      peak_name = paste0(base_name, direction_suffix)
    ) %>%
    dplyr::mutate(Chr = ifelse(grepl("^chr", Chr), Chr, paste0("chr", Chr))) %>%
    dplyr::filter(Chr %in% standard_chroms) %>%
    dplyr::select(Chr, Start, End, peak_name, logFC, FDR, Strand) %>%
    dplyr::mutate(
      score = as.integer(ifelse(logFC > report_params$logFC_cutoff,
                                pmin(1000, pmax(600, 600 + abs(logFC) * 100)),
                                pmin(500, pmax(200, 500 - abs(logFC) * 100)))),
      strand = Strand
    ) %>%
    dplyr::select(Chr, Start, End, peak_name, score, strand)

  if (nrow(bed_data) == 0) {
    cat("No peaks on standard chromosomes for", contrast_name, "\n")
    return(NULL)
  }

  cat("Sample directional peak names for", contrast_name, ":\n")
  sample_names <- head(bed_data$peak_name, 3)
  cat(paste(sample_names, collapse = "\n"), "\n")
  cat("Sample scores:", head(bed_data$score, 3), "\n")
  cat("Peaks retained after chromosome filtering:", nrow(bed_data), "\n\n")

  bed_file <- file.path(output_dir, paste0(contrast_name, "_DA_peaks.bed"))
  write.table(bed_data, bed_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  sorted_bed_file <- file.path(output_dir, paste0(contrast_name, "_DA_peaks_sorted.bed"))
  system2("sort", args = c("-k1,1", "-k2,2n", bed_file),
          stdout = sorted_bed_file, env = c("LC_COLLATE=C"))

  # Determine genome and get chrom.sizes file
  genome_alias <- if (organism == "mmu") "mm10" else "hg38"

  # Try to use bundled chrom.sizes file from package
  chrom_sizes_file <- system.file("extdata",
                                  paste0(genome_alias, ".chrom.sizes"),
                                  package = "atacreportR")

  # If not found in package, download to output directory
  if (chrom_sizes_file == "" || !file.exists(chrom_sizes_file)) {
    cat("Chrom.sizes not bundled with package, downloading...\n")
    chrom_sizes_file <- file.path(output_dir, paste0(genome_alias, ".chrom.sizes"))

    if (!file.exists(chrom_sizes_file)) {
      download.file(
        paste0("http://hgdownload.soe.ucsc.edu/goldenPath/",
               genome_alias, "/bigZips/", genome_alias, ".chrom.sizes"),
        chrom_sizes_file,
        quiet = TRUE
      )
      cat("Downloaded", basename(chrom_sizes_file), "\n")
    }
  } else {
    cat("Using bundled chrom.sizes file\n")
  }

  bedToBigBed_path <- Sys.which("bedToBigBed")
  if (bedToBigBed_path == "") {
    alt_path <- "/apps/ucsc/20210803/bedToBigBed"
    if (file.exists(alt_path)) {
      bedToBigBed_path <- alt_path
    } else {
      stop("bedToBigBed not found")
    }
  }

  bigbed_file <- file.path(output_dir, paste0(contrast_name, "_DA_peaks.bb"))
  result <- system2(bedToBigBed_path, args = c(sorted_bed_file, chrom_sizes_file, bigbed_file),
                    stdout = TRUE, stderr = TRUE)

  if (file.exists(bigbed_file) && file.size(bigbed_file) > 1000) {
    system2("chmod", args = c("644", bigbed_file))
    cat("Successfully created BigBed:", basename(bigbed_file), "\n")
    return(paste0(contrast_name, "_DA_peaks.bb"))
  } else {
    cat("BigBed conversion failed for", contrast_name, "\n")
    if (length(result) > 0) cat("Error:", paste(result, collapse = "\n"), "\n")
    return(NULL)
  }
}

#' Preprocess differential accessibility results for visualization
#'
#' Adds grouping and significance columns for plotting.
#'
#' @param df Data frame with DA results
#' @param contrast Character string contrast name
#'
#' @return Data frame with added visualization columns
#'
#' @importFrom dplyr case_when mutate
#' @importFrom stringr str_trim str_split
#' @export
preprocess_result <- function(df, contrast) {
  # Avoid NSE warnings
  FDR <- logFC <- NULL

  # Access report_params from parent environment
  report_params <- get("report_params", envir = parent.frame())
  # Avoid NSE warnings
  FDR <- logFC <- NULL

  # Access report_params from parent environment
  report_params <- get("report_params", envir = parent.frame())

  groups <- stringr::str_trim(stringr::str_split(gsub("\\.", "-", contrast), "_vs_", simplify = TRUE))
  df$group1 <- groups[1]
  df$group2 <- groups[2]

  df$Accessibility <- dplyr::case_when(
    df$FDR < report_params$sig_cutoff & df$logFC > report_params$logFC_cutoff ~
      paste0("More accessible in ", df$group1),
    df$FDR < report_params$sig_cutoff & df$logFC < -report_params$logFC_cutoff ~
      paste0("More accessible in ", df$group2),
    TRUE ~ "Not Significant"
  )

  df$Significant <- dplyr::case_when(
    df$FDR < report_params$sig_cutoff & df$logFC > report_params$logFC_cutoff ~ "more",
    df$FDR < report_params$sig_cutoff & df$logFC < -report_params$logFC_cutoff ~ "less",
    TRUE ~ "ns"
  )

  df$negLogFDR <- -log10(df$FDR)

  df
}

#' Convert DESeqDataSet to DGEList
#'
#' Creates an edgeR DGEList object from a DESeqDataSet, transferring
#' counts, sample metadata (colData), and peak annotations (rowData).
#' Useful for converting from DESeq2 workflow to edgeR workflow.
#'
#' @param dds DESeqDataSet object with counts, and optionally colData and rowData
#' @param norm_method Normalization method: "TMM" (default), "none", or "manual"
#' @param norm_factors Optional numeric vector of manual normalization factors
#'   (required if norm_method = "manual"). Must have length equal to ncol(dds).
#'
#' @return DGEList object with:
#'   \item{counts}{Count matrix from dds}
#'   \item{samples}{Sample metadata from colData(dds) plus library sizes and norm factors}
#'   \item{genes}{Peak annotations from rowData(dds), if available}
#'
#' @examples
#' \dontrun{
#' # TMM normalization (default)
#' dge <- dds_to_dgelist(my_dds)
#'
#' # No normalization
#' dge_raw <- dds_to_dgelist(my_dds, norm_method = "none")
#'
#' # Manual normalization factors (e.g., from spike-ins)
#' spikein_factors <- calculate_spikein_factors(spikein_dds)
#' dge <- dds_to_dgelist(my_dds, norm_method = "manual", norm_factors = spikein_factors)
#' }
#'
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom SummarizedExperiment assay colData rowData
#' @export
dds_to_dgelist <- function(dds, norm_method = "TMM", norm_factors = NULL) {

  # Validate input
  if (!inherits(dds, "DESeqDataSet")) {
    stop("dds must be a DESeqDataSet object")
  }

  # Extract counts
  counts <- SummarizedExperiment::assay(dds, "counts")

  # Create DGEList
  dge <- edgeR::DGEList(counts = counts)

  # Add sample metadata from colData
  coldata_df <- as.data.frame(SummarizedExperiment::colData(dds))
  if (ncol(coldata_df) > 0) {
    dge$samples <- cbind(dge$samples, coldata_df)
    cat("Added", ncol(coldata_df), "sample metadata columns to dge$samples\n")
  }

  # Add peak annotations from rowData if available
  rowdata_df <- as.data.frame(SummarizedExperiment::rowData(dds))
  if (ncol(rowdata_df) > 0) {
    # Ensure rownames match
    if ("interval" %in% colnames(rowdata_df)) {
      rownames(rowdata_df) <- rowdata_df$interval
    }
    dge$genes <- rowdata_df[rownames(dge), ]
    cat("Added", ncol(rowdata_df), "peak annotation columns to dge$genes\n")
  } else {
    cat("No peak annotations in rowData - dge$genes is empty\n")
  }

  # Apply normalization
  if (norm_method == "TMM") {
    cat("Calculating TMM normalization factors...\n")
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    cat("TMM normalization complete. Factors range:",
        round(min(dge$samples$norm.factors), 3), "-",
        round(max(dge$samples$norm.factors), 3), "\n")

  } else if (norm_method == "manual") {
    if (is.null(norm_factors)) {
      stop("norm_factors must be provided when norm_method = 'manual'")
    }
    if (length(norm_factors) != ncol(dge)) {
      stop("Length of norm_factors (", length(norm_factors),
           ") must match number of samples (", ncol(dge), ")")
    }
    cat("Applying manual normalization factors...\n")
    dge$samples$norm.factors <- norm_factors
    cat("Manual normalization applied. Factors range:",
        round(min(norm_factors), 3), "-",
        round(max(norm_factors), 3), "\n")

  } else if (norm_method == "none") {
    cat("No normalization applied (all norm.factors = 1)\n")
    dge$samples$norm.factors <- rep(1, ncol(dge))

  } else {
    stop("norm_method must be 'TMM', 'manual', or 'none'")
  }

  cat("Created DGEList:", nrow(dge), "peaks x", ncol(dge), "samples\n")

  return(dge)
}
