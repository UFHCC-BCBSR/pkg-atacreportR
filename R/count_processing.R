#' Count BAM reads in consensus peaks using featureCounts
#'
#' Uses Rsubread::featureCounts to quantify reads overlapping consensus peaks
#' across multiple BAM files.
#'
#' @param bam_files Named character vector of BAM file paths
#' @param consensus_peaks GRanges object containing consensus peak regions
#'
#' @return Matrix of read counts with peaks as rows and samples as columns
#'
#' @importFrom Rsubread featureCounts
#' @importFrom GenomicRanges seqnames start end
#' @export
count_bam_reads_in_peaks <- function(bam_files, consensus_peaks) {
  cat("Counting reads in", length(consensus_peaks), "consensus peaks across",
      length(bam_files), "samples...\n")

  saf_file <- create_saf_file(consensus_peaks)

  bam_vector <- as.character(bam_files)
  names(bam_vector) <- names(bam_files)

  missing_bams <- bam_vector[!file.exists(bam_vector)]
  if (length(missing_bams) > 0) {
    stop("BAM files not found: ", paste(names(missing_bams), collapse = ", "))
  }

  cat("Running featureCounts...\n")

  fc_result <- featureCounts(
    files = bam_vector,
    annot.ext = saf_file,
    isGTFAnnotationFile = FALSE,
    isPairedEnd = TRUE,
    nthreads = 14,
    verbose = TRUE
  )

  unlink(saf_file)

  count_matrix <- fc_result$counts
  colnames(count_matrix) <- names(bam_files)
  rownames(count_matrix) <- names(consensus_peaks)

  cat("Counting complete:", nrow(count_matrix), "peaks x", ncol(count_matrix), "samples\n")

  return(count_matrix)
}

#' Create SAF format annotation file for featureCounts
#'
#' Converts a GRanges object to SAF (Simplified Annotation Format) for use
#' with Rsubread::featureCounts.
#'
#' @param consensus_peaks GRanges object containing genomic regions
#'
#' @return Character string path to temporary SAF file
#'
#' @importFrom GenomicRanges seqnames start end
#' @keywords internal
create_saf_file <- function(consensus_peaks) {
  saf_df <- data.frame(
    GeneID = names(consensus_peaks),
    Chr = as.character(seqnames(consensus_peaks)),
    Start = start(consensus_peaks),
    End = end(consensus_peaks),
    Strand = "+",
    stringsAsFactors = FALSE
  )

  temp_saf <- tempfile(fileext = ".saf")
  write.table(saf_df, file = temp_saf, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)

  cat("Created SAF file with", nrow(saf_df), "features\n")

  return(temp_saf)
}
