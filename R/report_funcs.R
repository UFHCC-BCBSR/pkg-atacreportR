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
#' @param results_df Data frame with DA results including FDR, logFC, logCPM, and annotation columns
#' @param contrast_name Character string contrast name (e.g., "Heat_vs_Control")
#' @param output_dir Output directory path for generated files
#' @param organism Character string organism code ("mmu" for mouse, "hsa" for human)
#' @param sig_cutoff FDR significance threshold (default: 0.05)
#' @param logfc_cutoff Absolute log fold-change threshold (default: 1)
#' @param min_logcpm Minimum average accessibility threshold (default: 2)
#' @param bedToBigBed_path Optional path to bedToBigBed executable (auto-detected if NULL)
#'
#' @return Character string BigBed filename or NULL if creation failed
#'
#' @importFrom dplyr filter mutate select
#' @export
create_da_bigbed <- function(results_df,
                             contrast_name,
                             output_dir,
                             organism,
                             sig_cutoff = 0.05,
                             logfc_cutoff = 1,
                             min_logcpm = 2,
                             bedToBigBed_path = NULL) {

  # Avoid NSE warnings
  FDR <- logFC <- logCPM <- Gene.Name <- Chr <- Start <- End <- Strand <- NULL
  base_name <- direction_suffix <- peak_name <- score <- strand <- interval <- NULL
  . <- NULL

  # Filter for significant peaks
  sig_peaks <- results_df %>%
    dplyr::filter(
      FDR < sig_cutoff,
      abs(logFC) > logfc_cutoff,
      logCPM > min_logcpm
    )

  # Filter for peaks with gene names if Gene.Name column exists and has values
  if ("Gene.Name" %in% colnames(sig_peaks)) {
    n_before <- nrow(sig_peaks)
    sig_peaks <- sig_peaks %>%
      dplyr::filter(!is.na(Gene.Name) & Gene.Name != "")
    cat("Filtered", n_before - nrow(sig_peaks), "peaks without gene names\n")
  }

  if (nrow(sig_peaks) == 0) {
    cat("No significant peaks found for", contrast_name, "\n")
    return(NULL)
  }

  # Extract group names from contrast
  group_names <- strsplit(contrast_name, "_vs_")[[1]]
  group1 <- group_names[1]
  group2 <- if (length(group_names) > 1) group_names[2] else "baseline"

  # Detect column name formats (HOMER vs ChIPseeker)
  chr_col <- if ("Chr" %in% colnames(sig_peaks)) "Chr" else "seqnames"
  start_col <- if ("Start" %in% colnames(sig_peaks)) "Start" else "start"
  end_col <- if ("End" %in% colnames(sig_peaks)) "End" else "end"
  strand_col <- if ("Strand" %in% colnames(sig_peaks)) "Strand" else "strand"

  cat("Processing", nrow(sig_peaks), "significant peaks for", contrast_name, "\n")
  cat("Using columns - Chr:", chr_col, "Start:", start_col, "End:", end_col, "\n")

  # Define standard chromosomes
  standard_chroms <- if (organism == "mmu") {
    paste0("chr", c(1:19, "X", "Y", "M"))
  } else {
    paste0("chr", c(1:22, "X", "Y", "M"))
  }

  # Determine peak name: use Gene.Name if available, otherwise interval
  has_gene_names <- "Gene.Name" %in% colnames(sig_peaks)
  has_interval <- "interval" %in% colnames(sig_peaks)

  # Create BED format data
  bed_data <- sig_peaks %>%
    dplyr::mutate(
      Chr = as.character(.[[chr_col]]),
      Start = .[[start_col]],
      End = .[[end_col]],
      Strand = if (strand_col %in% colnames(.)) .[[strand_col]] else "."
    )

  # Add peak names
  if (has_gene_names) {
    bed_data <- bed_data %>%
      dplyr::mutate(base_name = Gene.Name)
  } else if (has_interval) {
    bed_data <- bed_data %>%
      dplyr::mutate(base_name = interval)
  } else {
    # Fallback: create names from coordinates
    bed_data <- bed_data %>%
      dplyr::mutate(base_name = paste0(Chr, ":", Start, "-", End))
  }

  # Add directional suffix
  bed_data <- bed_data %>%
    dplyr::mutate(
      direction_suffix = ifelse(logFC > logfc_cutoff,
                                paste0("_MORE_ACC_IN_", toupper(group1)),
                                paste0("_MORE_ACC_IN_", toupper(group2))),
      peak_name = paste0(base_name, direction_suffix)
    ) %>%
    # Ensure chr prefix
    dplyr::mutate(Chr = ifelse(grepl("^chr", Chr), Chr, paste0("chr", Chr))) %>%
    # Filter to standard chromosomes
    dplyr::filter(Chr %in% standard_chroms) %>%
    # Create score column (higher for more significant)
    dplyr::mutate(
      score = as.integer(ifelse(logFC > logfc_cutoff,
                                pmin(1000, pmax(600, 600 + abs(logFC) * 100)),
                                pmin(500, pmax(200, 500 - abs(logFC) * 100)))),
      strand = Strand
    ) %>%
    dplyr::select(Chr, Start, End, peak_name, score, strand)

  if (nrow(bed_data) == 0) {
    cat("No peaks on standard chromosomes for", contrast_name, "\n")
    return(NULL)
  }

  cat("Final peak count:", nrow(bed_data), "\n")
  cat("Sample peak names:", paste(head(bed_data$peak_name, 3), collapse = ", "), "\n\n")

  # Write BED file
  bed_file <- file.path(output_dir, paste0(contrast_name, "_DA_peaks.bed"))
  write.table(bed_data, bed_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  # Sort BED file
  sorted_bed_file <- file.path(output_dir, paste0(contrast_name, "_DA_peaks_sorted.bed"))
  system2("sort", args = c("-k1,1", "-k2,2n", bed_file),
          stdout = sorted_bed_file, env = c("LC_COLLATE=C"))

  # Get chrom.sizes file
  genome_alias <- if (organism == "mmu") "mm10" else "hg38"
  chrom_sizes_file <- .get_chrom_sizes(genome_alias, output_dir)

  # Find bedToBigBed executable
  if (is.null(bedToBigBed_path)) {
    bedToBigBed_path <- Sys.which("bedToBigBed")
    if (bedToBigBed_path == "") {
      # Try common locations
      alt_paths <- c(
        "/apps/ucsc/20210803/bedToBigBed",
        "/usr/local/bin/bedToBigBed",
        "/usr/bin/bedToBigBed"
      )
      for (path in alt_paths) {
        if (file.exists(path)) {
          bedToBigBed_path <- path
          break
        }
      }
    }
  }

  if (bedToBigBed_path == "" || !file.exists(bedToBigBed_path)) {
    stop("bedToBigBed executable not found. Please install UCSC tools or provide path.")
  }

  # Convert to BigBed
  bigbed_file <- file.path(output_dir, paste0(contrast_name, "_DA_peaks.bb"))
  result <- system2(bedToBigBed_path,
                    args = c(sorted_bed_file, chrom_sizes_file, bigbed_file),
                    stdout = TRUE, stderr = TRUE)

  # Check if successful
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

#' Get chromosome sizes file
#' @keywords internal
.get_chrom_sizes <- function(genome_alias, output_dir) {
  # Try package-bundled file first
  chrom_sizes_file <- system.file("extdata",
                                  paste0(genome_alias, ".chrom.sizes"),
                                  package = "atacreportR")

  # If not found, download to output directory
  if (chrom_sizes_file == "" || !file.exists(chrom_sizes_file)) {
    cat("Downloading chrom.sizes file for", genome_alias, "...\n")
    chrom_sizes_file <- file.path(output_dir, paste0(genome_alias, ".chrom.sizes"))

    if (!file.exists(chrom_sizes_file)) {
      download.file(
        paste0("http://hgdownload.soe.ucsc.edu/goldenPath/",
               genome_alias, "/bigZips/", genome_alias, ".chrom.sizes"),
        chrom_sizes_file,
        quiet = TRUE
      )
      cat("Downloaded chrom.sizes file\n")
    }
  } else {
    cat("Using bundled chrom.sizes file\n")
  }

  return(chrom_sizes_file)
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

#' Plot differential accessibility results by genomic annotation
#'
#' Creates interactive plotly visualizations showing the distribution of
#' genomic annotations for differential and non-differential peaks.
#'
#' @param results_list Named list of differential accessibility results from
#'   run_differential_analysis()
#' @param sig_cutoff FDR significance threshold (default: 0.05)
#' @param logfc_cutoff Absolute log fold-change threshold (default: 1)
#' @param plot_height Height of each plot in pixels (default: 400)
#' @param plot_width Width of each plot in pixels (default: 550)
#'
#' @return HTML tagList containing interactive plotly plots for each contrast
#'
#' @importFrom dplyr bind_rows filter mutate case_when
#' @importFrom stringr str_split
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_manual labs theme_minimal theme element_text margin
#' @importFrom plotly ggplotly layout subplot
#' @importFrom htmltools tagList tags
#' @export
plot_da_by_annotation <- function(results_list,
                                  sig_cutoff = 0.05,
                                  logfc_cutoff = 1,
                                  plot_height = 400,
                                  plot_width = 550) {

  # Check required columns
  if (length(results_list) == 0) {
    stop("results_list is empty")
  }

  test_df <- results_list[[1]]$table
  required_cols <- c("FDR", "logFC", "Annotation_short")
  missing <- setdiff(required_cols, colnames(test_df))
  if (length(missing) > 0) {
    stop("Missing required columns in results: ", paste(missing, collapse = ", "))
  }

  # Define a better color palette with good contrast
  annotation_colors <- c(
    "promoter-TSS" = "#E41A1C",    # Red
    "TTS" = "#377EB8",              # Blue
    "exon" = "#4DAF4A",             # Green
    "intron" = "#984EA3",           # Purple
    "Intergenic" = "#FF7F00"        # Orange
  )

  # STEP 1: Combine all contrasts
  all_results <- bind_rows(
    lapply(names(results_list), function(contrast) {
      df <- results_list[[contrast]]$table
      df$contrast <- contrast

      parts <- stringr::str_split(contrast, "_vs_", simplify = TRUE)
      left <- .clean_group_name(parts[1])
      right <- .clean_group_name(parts[2])

      df <- df %>%
        dplyr::mutate(
          significance_group = dplyr::case_when(
            FDR < sig_cutoff & logFC > logfc_cutoff ~
              paste0("More Accessible\nin ", right),
            FDR < sig_cutoff & logFC < -logfc_cutoff ~
              paste0("More Accessible\nin ", left),
            TRUE ~ "Not Significant"
          ),
          is_significant = significance_group != "Not Significant"
        )
      return(df)
    })
  )

  # STEP 2: Create plots for each contrast
  all_plots <- lapply(unique(all_results$contrast), function(cn) {
    contrast_df <- all_results %>%
      dplyr::filter(contrast == cn, !is.na(Annotation_short))

    parts <- stringr::str_split(cn, "_vs_", simplify = TRUE)
    contrast_display <- paste(.clean_group_name(parts[1]), "vs",
                              .clean_group_name(parts[2]))

    # Split into significant and non-significant
    sig_df <- contrast_df %>% dplyr::filter(is_significant == TRUE)
    nonsig_df <- contrast_df %>% dplyr::filter(is_significant == FALSE)

    # --- PLOT 1: Significant peaks ---
    p1 <- ggplot2::ggplot(sig_df, ggplot2::aes(x = significance_group, fill = Annotation_short)) +
      ggplot2::geom_bar(stat = "count", position = "stack", width = 0.6) +
      ggplot2::geom_text(stat = "count", ggplot2::aes(label = ggplot2::after_stat(count)),
                         position = ggplot2::position_stack(vjust = 0.5), size = 3.5,
                         color = "white", fontface = "bold") +
      ggplot2::scale_fill_manual(values = annotation_colors, name = "Genomic Region") +
      ggplot2::labs(
        title = paste("Differential Peaks:", contrast_display),
        subtitle = paste("FDR <", sig_cutoff, ", |logFC| >", logfc_cutoff),
        x = "",
        y = "Number of Peaks"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, face = "bold", size = 10),
        plot.title = ggplot2::element_text(size = 12, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
        plot.margin = ggplot2::margin(10, 10, 10, 10),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 10, face = "bold")
      )

    # --- PLOT 2: Non-significant peaks ---
    p2 <- ggplot2::ggplot(nonsig_df, ggplot2::aes(x = "Not\nSignificant", fill = Annotation_short)) +
      ggplot2::geom_bar(stat = "count", position = "stack", width = 0.6) +
      ggplot2::geom_text(stat = "count", ggplot2::aes(label = ggplot2::after_stat(count)),
                         position = ggplot2::position_stack(vjust = 0.5), size = 3.5,
                         color = "white", fontface = "bold") +
      ggplot2::scale_fill_manual(values = annotation_colors, name = "Genomic Region") +
      ggplot2::labs(
        title = paste("Non-Differential Peaks:", contrast_display),
        subtitle = "Peaks not meeting significance criteria",
        x = "",
        y = "Number of Peaks"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(face = "bold", size = 10),
        plot.title = ggplot2::element_text(size = 12, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
        plot.margin = ggplot2::margin(10, 10, 10, 10),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 10, face = "bold")
      )

    # Convert to plotly
    plotly1 <- plotly::ggplotly(p1, tooltip = c("fill", "count"))
    plotly2 <- plotly::ggplotly(p2, tooltip = c("fill", "count"))

    # Remove legend entries from second plot traces
    for (i in seq_along(plotly2$x$data)) {
      plotly2$x$data[[i]]$showlegend <- FALSE
    }

    # Combine plots side-by-side
    combined_plot <- plotly::subplot(
      plotly1, plotly2,
      nrows = 1,
      shareY = FALSE,
      titleX = TRUE,
      titleY = TRUE,
      margin = 0.08
    ) %>%
      plotly::layout(
        title = list(
          text = paste0("<b>Peak Annotation by Differential Accessibility: ", contrast_display, "</b><br>",
                        "<sub>FDR < ", sig_cutoff, ", |logFC| > ", logfc_cutoff, "</sub>"),
          x = 0.5,
          xanchor = "center",
          font = list(size = 14)
        ),
        height = plot_height,
        showlegend = TRUE,
        legend = list(
          orientation = "h",
          x = 0.5,
          xanchor = "center",
          y = -0.15,
          title = list(text = "<b>Genomic Region</b>", font = list(size = 11))
        ),
        annotations = list(
          list(text = "<b>Differential Peaks</b>",
               x = 0.25, y = 1.05,
               xref = "paper", yref = "paper",
               xanchor = "center", showarrow = FALSE,
               font = list(size = 12)),
          list(text = "<b>Non-Differential Peaks</b>",
               x = 0.75, y = 1.05,
               xref = "paper", yref = "paper",
               xanchor = "center", showarrow = FALSE,
               font = list(size = 12))
        )
      )

    # Return plot with section break
    list(
      combined_plot,
      htmltools::tags$div(style = "margin: 40px 0;")
    )
  })

  # Return all plots
  htmltools::tagList(unlist(all_plots, recursive = FALSE))
}

#' Clean group names for display
#' @keywords internal
.clean_group_name <- function(name) {
  name <- sub("^X", "", name)
  name <- gsub("_", " ", name)
  return(name)
}

#' Create MA plots for differential accessibility results
#'
#' Generates interactive MA plots showing the relationship between average
#' accessibility (log2 CPM) and differential accessibility (log2 fold change).
#'
#' @param results_list Named list of differential accessibility results from
#'   run_differential_analysis()
#' @param sig_cutoff FDR significance threshold (default: 0.05)
#' @param logfc_cutoff Absolute log fold-change threshold (default: 1)
#' @param point_size Size of points in plot (default: 5)
#' @param point_opacity Opacity of points (default: 0.6)
#'
#' @return HTML tagList containing interactive plotly MA plots for each contrast
#'
#' @importFrom plotly plot_ly layout
#' @importFrom htmltools tagList
#' @export
plot_ma <- function(results_list,
                    sig_cutoff = 0.05,
                    logfc_cutoff = 1,
                    point_size = 5,
                    point_opacity = 0.6) {

  # Check if results_list is empty
  if (length(results_list) == 0) {
    stop("results_list is empty")
  }

  # Create MA plots for each contrast
  ma_plots <- lapply(names(results_list), function(contrast) {
    df <- results_list[[contrast]]$table

    # Parse contrast name
    parts <- strsplit(contrast, "_vs_")[[1]]
    group1 <- .clean_group_name(parts[1])
    group2 <- .clean_group_name(parts[2])
    contrast_display <- paste(group1, "vs", group2)

    # Add accessibility classification
    df$Accessibility <- ifelse(
      df$FDR < sig_cutoff & df$logFC > logfc_cutoff,
      paste0("More accessible in ", group2),
      ifelse(
        df$FDR < sig_cutoff & df$logFC < -logfc_cutoff,
        paste0("More accessible in ", group1),
        "Not Significant"
      )
    )

    # Define colors
    color_values <- setNames(
      c("#E41A1C", "#377EB8", "#999999"),  # Red, Blue, Gray
      c(
        paste0("More accessible in ", group1),
        paste0("More accessible in ", group2),
        "Not Significant"
      )
    )

    # Create hover text
    df$hover_text <- paste0(
      "Gene: ", ifelse(is.na(df$Gene.Name) | df$Gene.Name == "", "N/A", df$Gene.Name), "<br>",
      "logFC: ", round(df$logFC, 2), "<br>",
      "logCPM: ", round(df$logCPM, 2), "<br>",
      "FDR: ", format(df$FDR, scientific = TRUE, digits = 3)
    )

    # Create plotly MA plot
    p <- plotly::plot_ly(
      data = df,
      x = ~logCPM,
      y = ~logFC,
      type = 'scatter',
      mode = 'markers',
      text = ~hover_text,
      hoverinfo = 'text',
      color = ~Accessibility,
      colors = color_values,
      marker = list(size = point_size, opacity = point_opacity)
    ) %>%
      plotly::layout(
        title = list(
          text = paste0("<b>MA Plot: ", contrast_display, "</b><br>",
                        "<sub>Average Accessibility vs Fold Change</sub>"),
          font = list(size = 14)
        ),
        xaxis = list(title = "Average Expression (log2 CPM)"),
        yaxis = list(title = "log2 Fold Change"),
        legend = list(
          title = list(text = "<b>Accessibility</b>"),
          orientation = "v",
          x = 1.02,
          xanchor = "left"
        ),
        hovermode = "closest",
        height = 500,
        width = 700
      )

    return(p)
  })

  # Return all plots as tagList
  htmltools::tagList(ma_plots)
}

#' Create volcano plots for differential accessibility results
#'
#' Generates interactive volcano plots showing the relationship between
#' log2 fold change and statistical significance (-log10 FDR).
#'
#' @param results_list Named list of differential accessibility results from
#'   run_differential_analysis()
#' @param sig_cutoff FDR significance threshold (default: 0.05)
#' @param logfc_cutoff Absolute log fold-change threshold (default: 1)
#' @param point_size Size of points in plot (default: 5)
#' @param point_opacity Opacity of points (default: 0.6)
#'
#' @return HTML tagList containing interactive plotly volcano plots for each contrast
#'
#' @importFrom plotly plot_ly layout
#' @importFrom htmltools tagList
#' @export
plot_volcano <- function(results_list,
                         sig_cutoff = 0.05,
                         logfc_cutoff = 1,
                         point_size = 5,
                         point_opacity = 0.6) {

  # Check if results_list is empty
  if (length(results_list) == 0) {
    stop("results_list is empty")
  }

  # Create volcano plots for each contrast
  volcano_plots <- lapply(names(results_list), function(contrast) {
    df <- results_list[[contrast]]$table

    # Parse contrast name
    parts <- strsplit(contrast, "_vs_")[[1]]
    group1 <- .clean_group_name(parts[1])
    group2 <- .clean_group_name(parts[2])
    contrast_display <- paste(group1, "vs", group2)

    # Calculate -log10(FDR)
    df$negLogFDR <- -log10(df$FDR)

    # Add accessibility classification
    df$Accessibility <- ifelse(
      df$FDR < sig_cutoff & df$logFC > logfc_cutoff,
      paste0("More accessible in ", group2),
      ifelse(
        df$FDR < sig_cutoff & df$logFC < -logfc_cutoff,
        paste0("More accessible in ", group1),
        "Not Significant"
      )
    )

    # Define colors
    color_values <- setNames(
      c("#E41A1C", "#377EB8", "#999999"),  # Red, Blue, Gray
      c(
        paste0("More accessible in ", group1),
        paste0("More accessible in ", group2),
        "Not Significant"
      )
    )

    # Create hover text
    df$hover_text <- paste0(
      "Gene: ", ifelse(is.na(df$Gene.Name) | df$Gene.Name == "", "N/A", df$Gene.Name), "<br>",
      "logFC: ", round(df$logFC, 2), "<br>",
      "-log10(FDR): ", round(df$negLogFDR, 2), "<br>",
      "FDR: ", format(df$FDR, scientific = TRUE, digits = 3)
    )

    # Create plotly volcano plot
    p <- plotly::plot_ly(
      data = df,
      x = ~logFC,
      y = ~negLogFDR,
      type = 'scatter',
      mode = 'markers',
      text = ~hover_text,
      hoverinfo = 'text',
      color = ~Accessibility,
      colors = color_values,
      marker = list(size = point_size, opacity = point_opacity)
    ) %>%
      plotly::layout(
        title = list(
          text = paste0("<b>Volcano Plot: ", contrast_display, "</b><br>",
                        "<sub>Fold Change vs Statistical Significance</sub>"),
          font = list(size = 14)
        ),
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10(FDR)"),
        legend = list(
          title = list(text = "<b>Accessibility</b>"),
          orientation = "v",
          x = 1.02,
          xanchor = "left"
        ),
        hovermode = "closest",
        height = 500,
        width = 700
      )

    return(p)
  })

  # Return all plots as tagList
  htmltools::tagList(volcano_plots)
}
