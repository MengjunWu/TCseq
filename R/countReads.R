#' count mapped reads overlap genomic intervals
#'
#' This function counts mapped reads from multiple BAM files that overlap genomic intervals
#' in \code{genomicFeature} in a \code{TCA} object. The counting result is stored in '\code{count}'
#' slot of the \code{TCA} object.
#'
#' @param object a \code{TCA} object
#'
#' @param dir character string giving the directory where BAM files are stored.
#'
#' @param method character string giving the counting method. Options are '\code{summarizeOverlaps}'
#' and '\code{featureCounts}'. For Windows system, only '\code{summarizeOverlaps}' can be used, For Linux
#' system, both methods can be used.
#'
#' @param ... additional arguments passed to \code{\link{summarizeOverlaps}} and \code{featureCounts} in \code{Rsubread} package
#'
#' @details
#' This function provides two options '\code{summarizeOverlaps}' from GenomicAlignments package and '\code{featureCounts}'
#' from Rsubread package to count the aligned reads. As Rsubread package is only avaible for linux systems,
#' Windows users can only use '\code{summarizeOverlaps}'. The user could specify counting details by passing additional
#' arguments (...), otherwise the default settings of the two methods are used. For counting details, see
#' \code{\link{summarizeOverlaps}}, \code{featureCounts} in \code{Rsubread} package
#'
#'
#' @return
#' A TCA object with updated '\code{count}' slot.
#'
#' @author
#' Mengjun Wu
#'
#' @seealso
#' \code{\link{summarizeOverlaps}}, \code{featureCounts} in \code{Rsubread} package
#'
#'
#' @export
countReads <- function(object, dir, method = "summarizeoverlaps", ...) {
  if (!"BAMfile" %in% colnames(object@design)) {
    err <- paste0("Can not find 'BAMfile' in design, please check if the correspoinding field is missing or uses a different column name.")
    stop(err)
  }
  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  setwd(dir)
  bamfiles <- as.vector(object@design$BAMfile)
  features <- object@genomicFeature
  ignore.strand <- NULL
  if (is.null(features$strand)) {
    warning("No strand information is provided, strand is ignored in reads counting")
    ignore.strand <- TRUE
  }
  gi <- GenomicRanges::GRanges(seqnames = features$chr, strand = features$strand, id = features$id, ranges = IRanges::IRanges(features$start +
                                                                                                                                1, width = features$end - features$start))
  method1 <- tolower(method)
  if (method1 == "featureCounts" && .Platform$OS.type == "windows") {
    stop(" 'featureCounts' is only available in Linux/Mac OS system.")
  }
  count.table <- switch(method1, summarizeoverlaps = {
    bamfl <- Rsamtools::BamFileList(bamfiles, yieldSize = 1e+06)
    c <- GenomicAlignments::summarizeOverlaps(gi, bamfl, ignore.strand = ignore.strand, ...)
    count.table <- SummarizedExperiment::assays(c)$counts
    row.names(count.table) <- features$id
    count.table
  }, featureCounts = {
    warning("To use the featureCounts, you need first load 'Rsubread' pacakge.")
    gi_rsubread <- createAnnotationFile(gi)
    count.table <- NULL
    stra <- 0
    if (!ignore.strand) {
      stra <- 1
    }
    for (i in bamfiles) {
      m <- paste0("Mapping bamfile ", i, " to the reference peakset")
      message(m)
      o <- capture.output(x <- featureCounts(i, annot.ext = gi_rsubread, strandSpecific = stra, ...))
      count.table <- cbind(count.table, x$counts)
    }
    rm(o)
    count.table
  })
  count.table <- as.matrix(count.table)
  counts(object) <- count.table
  object
}
