#' Create consensus peaks from multiple samples
#'
#' Merges peaks across samples and retains only reproducible peaks
#' found in multiple samples. Uses a reproducibility threshold requiring
#' peaks to be present in at least 30% of samples (minimum 2).
#'
#' @param peak_files Named character vector of peak file paths (BED format).
#'   Names should be sample identifiers.
#' @param min_samples Integer minimum number of samples required for consensus
#'   (default: NULL, auto-calculates as max(2, 30% of samples))
#'
#' @return GRanges object containing consensus peaks with names "Peak_1", "Peak_2", etc.
#'
#' @examples
#' \dontrun{
#' peak_files <- c(
#'   "sample1" = "sample1_peaks.bed",
#'   "sample2" = "sample2_peaks.bed"
#' )
#' consensus <- make_consensus_peaks(peak_files)
#' }
#'
#' @importFrom GenomicRanges GRanges reduce findOverlaps
#' @importFrom IRanges IRanges
#' @export
make_consensus_peaks <- function(peak_files, min_samples = NULL) {

  # Validate inputs
  if (length(peak_files) == 0) {
    stop("peak_files cannot be empty")
  }

  if (is.null(names(peak_files)) || any(names(peak_files) == "")) {
    stop("All peak files must be named with sample identifiers")
  }

  cat("Creating consensus peaks from", length(peak_files), "samples...\n")

  # Read all peak files
  all_peaks <- list()
  for (sample_name in names(peak_files)) {
    file_path <- peak_files[[sample_name]]
    cat("Reading peaks for", sample_name, "from", basename(file_path), "\n")
    peaks_gr <- .read_peak_file(file_path, sample_name)
    all_peaks[[sample_name]] <- peaks_gr
    cat("   ->", length(peaks_gr), "peaks found\n")
  }

  # Combine all peaks
  combined_peaks <- do.call(c, unname(all_peaks))
  cat("Total peaks before filtering:", length(combined_peaks), "\n")

  # Reduce overlapping peaks
  reduced_all <- GenomicRanges::reduce(combined_peaks)

  # Count how many samples each peak appears in
  overlaps <- GenomicRanges::findOverlaps(combined_peaks, reduced_all)
  peak_counts <- table(subjectHits(overlaps))

  # Determine reproducibility threshold
  if (is.null(min_samples)) {
    min_samples <- max(2, ceiling(length(peak_files) * 0.3))
  }

  # Filter for reproducible peaks
  reproducible_indices <- as.numeric(names(peak_counts)[peak_counts >= min_samples])
  consensus_peaks <- reduced_all[reproducible_indices]

  cat("After reproducibility filter (>=", min_samples, "samples):",
      length(consensus_peaks), "peaks\n")
  cat("Filtered out", length(reduced_all) - length(consensus_peaks),
      "singleton peaks\n")

  # Name the consensus peaks
  names(consensus_peaks) <- paste0("Peak_", seq_along(consensus_peaks))

  return(consensus_peaks)
}

#' Read peak file in BED format
#'
#' Reads a BED-format peak file (minimum 3 columns: chr, start, end)
#' and converts to GRanges object. Start coordinates are 0-based in BED
#' format and are converted to 1-based for GRanges.
#'
#' @param file_path Character string path to peak file (BED format)
#' @param sample_name Character string sample identifier
#'
#' @return GRanges object with peaks and sample metadata
#'
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom IRanges IRanges
#' @keywords internal
.read_peak_file <- function(file_path, sample_name) {

  if (!file.exists(file_path)) {
    stop("Peak file not found: ", file_path)
  }

  # Try to read the file
  df <- tryCatch({
    read.delim(file_path, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  }, error = function(e) {
    stop("Failed to read peak file ", file_path, ": ", e$message)
  })

  if (ncol(df) < 3) {
    stop("Peak file must have at least 3 columns (chr, start, end): ", file_path)
  }

  if (nrow(df) == 0) {
    warning("Peak file is empty: ", file_path)
    return(GRanges())
  }

  # Create GRanges (convert from 0-based BED to 1-based GRanges)
  peaks_gr <- GRanges(
    seqnames = df[, 1],
    ranges = IRanges(start = df[, 2] + 1, end = df[, 3])
  )

  GenomicRanges::mcols(peaks_gr)$sample <- sample_name

  return(peaks_gr)
}

#' Create DESeqDataSet from peak and BAM files
#'
#' Integrates consensus peak calling, read counting, and sample metadata
#' to create a DESeqDataSet object ready for differential accessibility analysis.
#'
#' @param peak_files Named character vector of peak file paths (BED format).
#'   Names should match sample identifiers in sample_info.
#' @param bam_files Named character vector of BAM file paths.
#'   Names should match sample identifiers in sample_info.
#' @param sample_info Data frame with sample metadata. Must have a "sample" column
#'   matching the names of peak_files and bam_files.
#' @param min_samples Integer minimum number of samples for consensus peaks
#'   (default: NULL, auto-calculates as max(2, 30% of samples))
#'
#' @return DESeqDataSet object with:
#'   - Count matrix as assay
#'   - Sample metadata in colData
#'   - Consensus peaks as rowRanges
#'
#' @examples
#' \dontrun{
#' dds <- create_dds_from_peaks(
#'   peak_files = c("s1" = "s1.bed", "s2" = "s2.bed"),
#'   bam_files = c("s1" = "s1.bam", "s2" = "s2.bam"),
#'   sample_info = data.frame(sample = c("s1", "s2"), condition = c("A", "B"))
#' )
#' }
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom SummarizedExperiment rowRanges<-
#' @export
create_dds_from_peaks <- function(peak_files,
                                  bam_files,
                                  sample_info,
                                  min_samples = NULL) {

  # Validate inputs
  if (is.null(names(peak_files)) || any(names(peak_files) == "")) {
    stop("peak_files must be a named vector with sample identifiers")
  }

  if (is.null(names(bam_files)) || any(names(bam_files) == "")) {
    stop("bam_files must be a named vector with sample identifiers")
  }

  if (!identical(sort(names(peak_files)), sort(names(bam_files)))) {
    stop("Sample names in peak_files and bam_files must match")
  }

  cat("Creating DDS object from peak and BAM files...\n")

  # Create consensus peaks
  consensus_peaks <- make_consensus_peaks(peak_files, min_samples = min_samples)

  if (length(consensus_peaks) == 0) {
    stop("No consensus peaks found. Check peak files and reproducibility threshold.")
  }

  # Count reads in peaks
  count_matrix <- count_bam_reads_in_peaks(bam_files, consensus_peaks)

  # Match sample info to count matrix
  sample_info_matched <- .match_sample_info_to_counts(sample_info, count_matrix)

  # Create DESeq2 object
  cat("Creating DESeq2 object...\n")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_info_matched,
    design = ~ 1
  )

  # Add peak coordinates as rowRanges
  SummarizedExperiment::rowRanges(dds) <- consensus_peaks

  cat("DDS object created:", nrow(dds), "peaks x", ncol(dds), "samples\n")

  return(dds)
}

#' Match sample info to count matrix columns
#'
#' Ensures sample metadata is aligned with count matrix columns.
#' Creates basic sample info if none provided.
#'
#' @param sample_info Data frame with sample metadata
#' @param count_matrix Matrix of counts
#'
#' @return Data frame with sample info matched to count matrix columns
#'
#' @keywords internal
.match_sample_info_to_counts <- function(sample_info, count_matrix) {

  count_sample_names <- colnames(count_matrix)

  if (is.null(sample_info) || nrow(sample_info) == 0) {
    cat("No sample info provided, creating basic sample info...\n")
    sample_info_matched <- data.frame(
      sample = count_sample_names,
      row.names = count_sample_names,
      stringsAsFactors = FALSE
    )
  } else {
    # Match sample info to count matrix
    missing_samples <- setdiff(count_sample_names, sample_info$sample)
    if (length(missing_samples) > 0) {
      warning("Samples in count matrix not found in sample_info: ",
              paste(missing_samples, collapse = ", "))
    }

    sample_info_matched <- sample_info[match(count_sample_names, sample_info$sample), ]
    rownames(sample_info_matched) <- count_sample_names
  }

  return(sample_info_matched)
}
