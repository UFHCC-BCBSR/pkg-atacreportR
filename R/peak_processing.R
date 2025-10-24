#' Create consensus peaks from multiple samples
#'
#' Merges peaks across samples and retains only reproducible peaks
#' found in multiple samples.
#'
#' @param peak_files Named character vector of peak file paths
#'
#' @return GRanges object containing consensus peaks
#'
#' @importFrom GenomicRanges GRanges reduce findOverlaps
#' @importFrom IRanges IRanges
#' @export
make_consensus_peaks <- function(peak_files) {
  # Avoid NSE warnings
  sample <- NULL

  cat("Creating consensus peaks from", length(peak_files), "samples...\n")

  all_peaks <- list()

  for (sample_name in names(peak_files)) {
    file_path <- peak_files[[sample_name]]
    cat("Reading peaks for", sample_name, "from", basename(file_path), "\n")
    peaks_gr <- read_peak_file(file_path, sample_name)
    all_peaks[[sample_name]] <- peaks_gr
    cat("   ->", length(peaks_gr), "peaks found\n")
  }

  combined_peaks <- do.call(c, unname(all_peaks))
  cat("Total peaks before filtering:", length(combined_peaks), "\n")

  reduced_all <- GenomicRanges::reduce(combined_peaks)
  overlaps <- findOverlaps(combined_peaks, reduced_all)
  peak_counts <- table(subjectHits(overlaps))

  min_samples <- max(2, ceiling(length(peak_files) * 0.3))
  reproducible_indices <- as.numeric(names(peak_counts)[peak_counts >= min_samples])
  consensus_peaks <- reduced_all[reproducible_indices]

  cat("After reproducibility filter (>=", min_samples, "samples):", length(consensus_peaks), "peaks\n")
  cat("Filtered out", length(reduced_all) - length(consensus_peaks), "singleton peaks\n")

  names(consensus_peaks) <- paste0("Peak_", seq_along(consensus_peaks))

  return(consensus_peaks)
}

#' Read peak file in BED format
#'
#' @param file_path Character string path to peak file
#' @param sample_name Character string sample identifier
#'
#' @return GRanges object with peaks and sample metadata
#'
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom IRanges IRanges
#' @keywords internal
read_peak_file <- function(file_path, sample_name) {
  if (!file.exists(file_path)) {
    stop("Peak file not found: ", file_path)
  }

  df <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)

  if (ncol(df) < 3) {
    stop("Peak file must have at least 3 columns (chr, start, end): ", file_path)
  }

  peaks_gr <- GRanges(
    seqnames = df[,1],
    ranges = IRanges(start = df[,2] + 1, end = df[,3])
  )

  mcols(peaks_gr)$sample <- sample_name

  return(peaks_gr)
}

#' Create DESeqDataSet from peak and BAM files
#'
#' Integrates consensus peak calling, read counting, and sample metadata
#' to create a DESeqDataSet object.
#'
#' @param peak_files Named character vector of peak file paths
#' @param bam_files Named character vector of BAM file paths
#' @param sample_info Data frame with sample metadata
#' @param organism Character string: "mmu" or "hsa"
#'
#' @return DESeqDataSet object with counts and metadata
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom SummarizedExperiment rowRanges<-
#' @export
create_dds_from_peaks <- function(peak_files, bam_files, sample_info, organism) {
  cat("Creating DDS object from peak and BAM files...\n")

  consensus_peaks <- make_consensus_peaks(peak_files)
  count_matrix <- count_bam_reads_in_peaks(bam_files, consensus_peaks)
  sample_info_matched <- match_sample_info_to_counts(sample_info, count_matrix)

  cat("Creating DESeq2 object...\n")
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_info_matched,
    design = ~ 1
  )

  rowRanges(dds) <- consensus_peaks

  cat("DDS object created:", nrow(dds), "peaks x", ncol(dds), "samples\n")

  return(dds)
}

#' Match sample info to count matrix columns
#'
#' @param sample_info Data frame with sample metadata
#' @param count_matrix Matrix of counts
#'
#' @return Data frame with sample info matched to count matrix
#'
#' @keywords internal
match_sample_info_to_counts <- function(sample_info, count_matrix) {
  count_sample_names <- colnames(count_matrix)

  if (is.null(sample_info) || nrow(sample_info) == 0) {
    cat("No sample info provided, creating basic sample info...\n")
    sample_info_matched <- data.frame(
      sample = count_sample_names,
      row.names = count_sample_names,
      stringsAsFactors = FALSE
    )
  } else {
    sample_info_matched <- sample_info[match(count_sample_names, sample_info$sample), ]
    rownames(sample_info_matched) <- count_sample_names
  }

  return(sample_info_matched)
}
