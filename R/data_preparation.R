#' Prepare analysis data from various input sources
#'
#' Main data preparation function that loads or creates DDS objects, peak annotations,
#' and QC data from user-provided files or generates them on-the-fly.
#'
#' @param report_params List containing analysis parameters including file paths,
#'   organism specification, and output directory
#'
#' @return List containing:
#'   \item{dds}{DESeqDataSet object with count data}
#'   \item{sample_info}{Data frame with sample metadata}
#'   \item{peaks_anno}{Data frame with peak annotations}
#'   \item{file_specs}{List of input file specifications}
#'   \item{qc_data}{List of QC metrics (optional)}
#'   \item{bigwig_files}{Named vector of BigWig file paths}
#'
#' @export
prepare_analysis_data <- function(report_params) {
  cat("Loading and preparing analysis data...\n")

  outdir <- .resolve_output_dir(report_params)
  seqID  <- report_params$seqID %||% "project"

  file_specs <- parse_file_specifications(report_params)
  sample_info <- load_sample_metadata(report_params$sample_sheet)

  created_dds  <- FALSE
  created_anno <- FALSE

  if (!is.null(file_specs$dds_file) && file.exists(file_specs$dds_file)) {
    cat("Loading existing DDS file...\n")
    dds <- load_dds_file(file_specs$dds_file)
  } else {
    cat("Creating DDS from peak and BAM files...\n")
    dds <- create_dds_from_peaks(
      peak_files  = file_specs$peak_files,
      bam_files   = file_specs$bam_files,
      sample_info = sample_info,
      organism    = report_params$organism
    )
    created_dds <- TRUE
  }

  if (!is.null(file_specs$peak_annotation) && file.exists(file_specs$peak_annotation)) {
    cat("Loading existing peak annotation...\n")
    peak_anno <- load_peak_annotation(file_specs$peak_annotation)
  } else {
    cat("Generating peak annotation with ChIPseeker...\n")
    peak_anno <- generate_peak_annotation(dds, report_params$organism)
    created_anno <- TRUE
  }

  qc_data <- load_qc_data(file_specs$qc_files)
  validate_data_consistency(dds, sample_info, peak_anno)

  if (created_dds) {
    dds_path <- file.path(outdir, paste0(seqID, ".dds.RData"))
    cat("Saving DDS to:", dds_path, "\n")
    dds_to_save <- dds
    save(dds_to_save, file = dds_path)

    peaks_txt <- file.path(outdir, paste0(seqID, ".consensus-peaks.txt"))
    cat("Writing consensus peaks to:", peaks_txt, "\n")
    .write_consensus_peaks(SummarizedExperiment::rowRanges(dds), peaks_txt)
  } else {
    cat("DDS was provided by user; not overwriting on disk.\n")
  }

  if (created_anno) {
    anno_path <- file.path(outdir, paste0(seqID, ".annotated.consensus-peaks.txt"))
    cat("Saving peak annotation to:", anno_path, "\n")

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
    utils::write.table(peak_anno_clean, file = anno_path, sep = "\t",
                       row.names = FALSE, col.names = TRUE, quote = FALSE)
  } else {
    cat("Peak annotation was provided by user; not overwriting on disk.\n")
  }

  cat("Data preparation complete!\n")

  return(list(
    dds          = dds,
    sample_info  = sample_info,
    peaks_anno   = peak_anno,
    file_specs   = file_specs,
    qc_data      = qc_data,
    bigwig_files = file_specs$bigwig_files
  ))
}

#' Parse file specifications from report parameters
#'
#' Extracts and parses file paths for peaks, BAMs, BigWigs, and QC data
#' from comma-separated parameter strings.
#'
#' @param report_params List of analysis parameters
#'
#' @return List containing parsed file specifications
#'
#' @keywords internal
parse_file_specifications <- function(report_params) {
  peak_files <- NULL
  if (!is.null(report_params$peak_files) && report_params$peak_files != "") {
    pairs <- strsplit(report_params$peak_files, ",")[[1]]
    peak_files <- setNames(
      sapply(pairs, function(x) strsplit(x, ":")[[1]][2]),
      sapply(pairs, function(x) strsplit(x, ":")[[1]][1])
    )
  }

  bam_files <- NULL
  if (!is.null(report_params$bam_files) && report_params$bam_files != "") {
    pairs <- strsplit(report_params$bam_files, ",")[[1]]
    bam_files <- setNames(
      sapply(pairs, function(x) strsplit(x, ":")[[1]][2]),
      sapply(pairs, function(x) strsplit(x, ":")[[1]][1])
    )
  }

  bigwig_files <- NULL
  if (!is.null(report_params$bigwig_files) && report_params$bigwig_files != "") {
    pairs <- strsplit(report_params$bigwig_files, ",")[[1]]
    bigwig_files <- setNames(
      sapply(pairs, function(x) strsplit(x, ":")[[1]][2]),
      sapply(pairs, function(x) strsplit(x, ":")[[1]][1])
    )
  }

  list(
    dds_file = if (is.null(report_params$dds_file) || report_params$dds_file == "") NULL else report_params$dds_file,
    peak_files = peak_files,
    bam_files  = bam_files,
    peak_annotation = if (is.null(report_params$peak_annotation) || report_params$peak_annotation == "") NULL else report_params$peak_annotation,
    bigwig_files = bigwig_files,
    qc_files = list(
      flagstat_dir = report_params$qc_flagstat_dir %||% report_params[["qc_flagstat_dir"]] %||% NULL,
      frip_file    = report_params$qc_frip_file   %||% report_params[["qc_frip_file"]]   %||% NULL
    )
  )
}

#' Resolve output directory and create if needed
#' @keywords internal
.resolve_output_dir <- function(report_params) {
  outdir <- report_params$output_path %||% report_params[["output-path"]]
  if (is.null(outdir) || outdir == "") stop("output_path / --output-path is required in params.")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir
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

#' Null coalescing operator
#'
#' Returns y if x is NULL or empty string, otherwise returns x
#'
#' @param x First value to check
#' @param y Fallback value
#' @return Either x or y
#' @keywords internal
#' @name null-default
NULL

#' @rdname null-default
`%||%` <- function(x, y) if (is.null(x) || (is.character(x) && x == "")) y else x
