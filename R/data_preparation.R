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
#' @return List containing:
#'   \item{dds}{DESeqDataSet object with count data}
#'   \item{sample_info}{Data frame with sample metadata}
#'   \item{peaks_anno}{Data frame with peak annotations}
#'   \item{qc_data}{List of QC metrics (optional)}
#'   \item{bigwig_files}{Named vector of BigWig file paths}
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
      sample_info = sample_info,
      organism    = organism
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

  # Validate consistency
  validate_data_consistency(dds, sample_info, peaks_anno)

  # Save generated files
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
    .write_peak_annotation(peaks_anno, anno_path)
  } else {
    cat("Peak annotation was provided by user; not overwriting on disk.\n")
  }

  cat("Data preparation complete!\n")

  return(list(
    dds          = dds,
    sample_info  = sample_info,
    peaks_anno   = peaks_anno,
    qc_data      = qc_data,
    bigwig_files = bigwig_files
  ))
}

#' Write peak annotation to file
#' @keywords internal
.write_peak_annotation <- function(peak_anno, out_path) {
  peak_anno_clean <- peak_anno
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
  utils::write.table(df, file = out_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
