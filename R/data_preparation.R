#' Prepare analysis data from various input sources
#'
#' Main data preparation function that loads or creates DDS objects, peak annotations,
#' and QC data from user-provided files or generates them on-the-fly.
#'
#' @param sample_sheet Path to sample metadata file (CSV/TSV)
#' @param organism Organism code ("hsa" or "mmu")
#' @param output_dir Directory where generated files should be saved
#' @param project_id Project identifier for naming output files (default: "project")
#' @param dds_file Optional path to existing DDS RData file
#' @param peak_annotation Optional path to existing peak annotation file
#' @param peak_files Optional named vector of peak files (names = sample IDs, values = paths)
#' @param bam_files Optional named vector of BAM files (names = sample IDs, values = paths)
#' @param bigwig_files Optional named vector of BigWig files (names = sample IDs, values = paths)
#' @param qc_flagstat_dir Optional directory containing flagstat files
#' @param qc_frip_file Optional path to FRiP summary file
#'
#' @return DESeqDataSet object with:
#'   - counts in assays
#'   - sample metadata, QC metrics, and bigwig paths in colData
#'   - peak annotations in rowData
#'
#' @export
prepare_analysis_data <- function(sample_sheet,
                                  organism,
                                  output_dir,
                                  project_id = "project",
                                  dds_file = NULL,
                                  peak_annotation = NULL,
                                  peak_files = NULL,
                                  bam_files = NULL,
                                  bigwig_files = NULL,
                                  qc_flagstat_dir = NULL,
                                  qc_frip_file = NULL) {

  cat("Loading and preparing analysis data...\n")

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Load sample metadata
  sample_info <- load_sample_metadata(sample_sheet)

  # Track what we create vs what was provided
  created_dds <- FALSE
  created_anno <- FALSE

  # Load or create DDS
  if (!is.null(dds_file) && file.exists(dds_file)) {
    cat("Loading existing DDS file...\n")
    dds <- load_dds_file(dds_file)
  } else {
    if (is.null(peak_files) || is.null(bam_files)) {
      stop("Either dds_file must be provided, or both peak_files and bam_files must be provided")
    }
    cat("Creating DDS from peak and BAM files...\n")
    dds <- create_dds_from_peaks(
      peak_files  = peak_files,
      bam_files   = bam_files,
      sample_info = sample_info
    )
    created_dds <- TRUE
  }

  # Load or generate peak annotation
  if (!is.null(peak_annotation) && file.exists(peak_annotation)) {
    cat("Loading existing peak annotation...\n")
    peaks_anno <- load_peak_annotation(peak_annotation)
  } else {
    cat("Generating peak annotation with ChIPseeker...\n")
    peaks_anno <- generate_peak_annotation(dds, organism)
    created_anno <- TRUE
  }

  # Load QC data if available
  qc_data <- load_qc_data(
    flagstat_dir = qc_flagstat_dir,
    frip_file = qc_frip_file
  )

  # Validate consistency before merging
  validate_data_consistency(dds, sample_info, peaks_anno)

  # === MERGE EVERYTHING INTO DDS ===

  # 1. Merge sample metadata into colData
  cat("Merging sample metadata into DDS colData...\n")
  dds <- .merge_sample_info_to_dds(dds, sample_info)

  # 2. Merge peak annotations into rowData
  cat("Merging peak annotations into DDS rowData...\n")
  dds <- .merge_peak_anno_to_dds(dds, peaks_anno)

  # 3. Add QC metrics to colData (if available)
  if (!is.null(qc_data)) {
    cat("Adding QC metrics to DDS colData...\n")
    dds <- .add_qc_to_dds(dds, qc_data)
  }

  # 4. Add BigWig file paths to colData (if provided)
  if (!is.null(bigwig_files)) {
    cat("Adding BigWig file paths to DDS colData...\n")
    dds <- .add_bigwig_paths_to_dds(dds, bigwig_files)
  }

  # 5. Add peak file paths to colData (if provided) ### NEW!
  if (!is.null(peak_files)) {
    cat("Adding peak file paths to DDS colData...\n")
    dds <- .add_peak_paths_to_dds(dds, peak_files)
  }

  # === SAVE GENERATED FILES ===

  if (created_dds) {
    dds_path <- file.path(output_dir, paste0(project_id, ".dds.RData"))
    cat("Saving DDS to:", dds_path, "\n")
    dds_to_save <- dds
    save(dds_to_save, file = dds_path)

    peaks_txt <- file.path(output_dir, paste0(project_id, ".consensus-peaks.txt"))
    cat("Writing consensus peaks to:", peaks_txt, "\n")
    .write_consensus_peaks(SummarizedExperiment::rowRanges(dds), peaks_txt)
  } else {
    cat("DDS was provided by user; not overwriting on disk.\n")
  }

  if (created_anno) {
    anno_path <- file.path(output_dir, paste0(project_id, ".annotated.consensus-peaks.txt"))
    cat("Saving peak annotation to:", anno_path, "\n")
    .write_peak_annotation(SummarizedExperiment::rowData(dds), anno_path)
  } else {
    cat("Peak annotation was provided by user; not overwriting on disk.\n")
  }

  cat("\n=== Data preparation complete! ===\n")
  cat("DDS dimensions:", nrow(dds), "peaks x", ncol(dds), "samples\n")
  cat("colData columns:", paste(colnames(SummarizedExperiment::colData(dds)), collapse = ", "), "\n")
  cat("rowData columns:", paste(colnames(SummarizedExperiment::rowData(dds)), collapse = ", "), "\n")

  return(dds)
}

#' Merge sample info into DDS colData
#' @keywords internal
.merge_sample_info_to_dds <- function(dds, sample_info) {
  sample_info_matched <- sample_info[match(colnames(dds), sample_info$sample), ]

  if (!identical(sample_info_matched$sample, colnames(dds))) {
    stop("Sample names do not match between DDS and sample sheet")
  }

  existing_cols <- colnames(SummarizedExperiment::colData(dds))
  new_cols <- setdiff(colnames(sample_info_matched), c("sample", existing_cols))

  if (length(new_cols) > 0) {
    SummarizedExperiment::colData(dds)[new_cols] <- sample_info_matched[new_cols]
    cat("  Added", length(new_cols), "sample metadata columns\n")
  }

  return(dds)
}

#' Merge peak annotations into DDS rowData
#' @keywords internal
.merge_peak_anno_to_dds <- function(dds, peaks_anno) {
  # Ensure rownames match
  peaks_anno_matched <- peaks_anno[match(rownames(dds), peaks_anno$interval), ]

  if (!identical(peaks_anno_matched$interval, rownames(dds))) {
    stop("Peak intervals do not match between DDS and annotation")
  }

  # Add all annotation columns to rowData
  existing_cols <- colnames(SummarizedExperiment::rowData(dds))
  new_cols <- setdiff(colnames(peaks_anno_matched), c("interval", existing_cols))

  if (length(new_cols) > 0) {
    SummarizedExperiment::rowData(dds)[new_cols] <- peaks_anno_matched[new_cols]
    cat("  Added", length(new_cols), "peak annotation columns\n")
  }

  return(dds)
}

#' Add QC metrics to DDS colData
#' @keywords internal
.add_qc_to_dds <- function(dds, qc_data) {
  n_metrics <- 0

  # Add FRiP scores if available
  if (!is.null(qc_data$frip)) {
    # Match FRiP data to samples
    frip_df <- qc_data$frip

    # Try to match by sample name (handle various column name formats)
    sample_col <- NULL
    for (col in c("sample", "Sample", "sample_id", "ID")) {
      if (col %in% colnames(frip_df)) {
        sample_col <- col
        break
      }
    }

    if (!is.null(sample_col)) {
      frip_matched <- frip_df[match(colnames(dds), frip_df[[sample_col]]), ]

      # Add FRiP columns (excluding the sample name column)
      frip_cols <- setdiff(colnames(frip_df), c(sample_col, colnames(SummarizedExperiment::colData(dds))))
      if (length(frip_cols) > 0) {
        SummarizedExperiment::colData(dds)[frip_cols] <- frip_matched[frip_cols]
        n_metrics <- n_metrics + length(frip_cols)
      }
    }
  }

  # Could add flagstat metrics here in the future
  # if (!is.null(qc_data$flagstat)) { ... }

  if (n_metrics > 0) {
    cat("  Added", n_metrics, "QC metric column(s)\n")
  }

  return(dds)
}

#' Add BigWig file paths to DDS colData
#' @keywords internal
.add_bigwig_paths_to_dds <- function(dds, bigwig_files) {
  # Match bigwig files to DDS samples
  bigwig_matched <- bigwig_files[match(colnames(dds), names(bigwig_files))]

  # Add as a new column
  SummarizedExperiment::colData(dds)$bigwig_path <- bigwig_matched

  cat("  Added BigWig file paths for", sum(!is.na(bigwig_matched)), "samples\n")

  return(dds)
}

#' Add peak file paths to DDS colData
#' @keywords internal
.add_peak_paths_to_dds <- function(dds, peak_files) {
  # Match peak files to DDS samples
  peak_matched <- peak_files[match(colnames(dds), names(peak_files))]

  # Add as a new column
  SummarizedExperiment::colData(dds)$peak_file <- peak_matched

  cat("  Added peak file paths for", sum(!is.na(peak_matched)), "samples\n")

  return(dds)
}

#' Write peak annotation to file
#' @keywords internal
.write_peak_annotation <- function(peak_anno, out_path) {
  peak_anno_df <- as.data.frame(peak_anno)
  peak_anno_clean <- peak_anno_df

  for (col in names(peak_anno_clean)) {
    if (is.list(peak_anno_clean[[col]])) {
      peak_anno_clean[[col]] <- sapply(peak_anno_clean[[col]], function(x) {
        if (is.null(x) || length(x) == 0) {
          return(NA_character_)
        } else if (length(x) == 1) {
          return(as.character(x))
        } else {
          return(paste(x, collapse = ";"))
        }
      })
    }
  }

  utils::write.table(peak_anno_clean, file = out_path, sep = "\t",
                     row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#' Write consensus peaks to BED-like text file
#' @keywords internal
.write_consensus_peaks <- function(peaks_gr, out_file) {
  df <- data.frame(
    chr   = as.character(GenomicRanges::seqnames(peaks_gr)),
    start = as.integer(GenomicRanges::start(peaks_gr) - 1),
    end   = as.integer(GenomicRanges::end(peaks_gr)),
    name  = if (!is.null(names(peaks_gr))) names(peaks_gr) else paste0("Peak_", seq_along(peaks_gr)),
    stringsAsFactors = FALSE
  )

  utils::write.table(df, file = out_file, sep = "\t", row.names = FALSE,
                     col.names = TRUE, quote = FALSE)
}
