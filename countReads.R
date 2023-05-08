#' count mapped reads overlap genomic intervals
#'
#' This function counts mapped reads from multiple BAM files 
#' overlapping genomic intervals in \code{genomicFeature} in a 
#' \code{TCA} object. The resulted count table is stored in 
#' \code{count} slot of the \code{TCA} object.
#'
#' @param object a \code{TCA} object.
#'
#' @param dir character string giving the directory of BAM files.
#'
#' @param method character string giving the counting method. Options
#' are "\code{summarizeOverlaps}" and "\code{featureCounts}". For
#' Windows system, only "\code{summarizeOverlaps}" can be used, For
#' Linux system, both methods can be used.
#'
#' @param ... additional arguments passed to
#' \code{\link{summarizeOverlaps}} in GenomicAlignments package 
#' or \code{\link{featureCounts}} in Rsubread package.
#'
#' @param zero.based Logical. If TRUE, the start positions of the
#' genomic intervals are \emph{0-based}, if FALSE, the start positions
#' will be \emph{1-based}.
#'
#' @details
#' This function provides two options to count the mapped reads: 
#' "\code{summarizeOverlaps}" in the GenomicAlignments package and 
#' "\code{featureCounts}" in the Rsubread package. As Rsubread package 
#' is only avaible for linux systems, Windows users can only choose
#' "\code{summarizeOverlaps}". The user could further customize the 
#' counting paramters by passing additional arguments (...), otherwise 
#' the default settings of the two methods will be used. For details 
#' of the counting parameters, see \code{\link{summarizeOverlaps}}, 
#' \code{\link{featureCounts}}.
#'
#'
#' @return
#' A TCA object with updated \code{count} slot.
#'
#' @author
#' Mengjun Wu
#'
#' @seealso
#' \code{\link{summarizeOverlaps}}, \code{\link{featureCounts}} 
#'
#'
#' @export
countReads <- function(object, dir, method = "summarizeoverlaps",
                       zero.based = TRUE,...) {
  name.col.tmp <- colnames(object@design)
  name.col.tmp <- tolower(name.col.tmp)
  colnames(object@design) <- name.col.tmp
  if (!"bamfile" %in% colnames(object@design)) {
    err <- paste0("Can not find information of bam files in design, please check whether the correspoinding field is missing or the column name is the same as required.")
    stop(err)
  }
  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  setwd(dir)
  bamfiles <- as.vector(object@design$bamfile)
  features <- object@genomicFeature
  ignore.strand <- NULL
  if (is.null(features$strand)) {
    warning("No strand information is provided, strand is ignored in reads counting")
    ignore.strand <- TRUE
  }
  gi <- makeGRangesFromDataFrame(features, keep.extra.columns = TRUE,
                                 starts.in.df.are.0based = zero.based)
  method1 <- tolower(method)
  if (method1 == "featureCounts" && .Platform$OS.type == "windows") {
    stop(" 'featureCounts' is only available in Linux/Mac OS system.")
  }
  count.table <- switch(method1, summarizeoverlaps = {
    bamfl <- Rsamtools::BamFileList(bamfiles, yieldSize = 1e+06)
    c <- GenomicAlignments::summarizeOverlaps(gi, bamfl,
                                              ignore.strand = ignore.strand, ...)
    count.table <- SummarizedExperiment::assays(c)$counts
    row.names(count.table) <- features$id
    count.table
  }, featureCounts = {
    warning("To use the featureCounts, you need to load 'Rsubread' package first")
    gi_rsubread <- createAnnotationFile(gi)
    stra <- 0
    if (!ignore.strand) {
      stra <- 1
    }
    for (i in bamfiles) {
      m <- paste0("counting reads in bamfile ", i)
      message(m)
      o <- capture.output(x <- featureCounts(i, annot.ext = gi_rsubread,
                                             strandSpecific = stra, ...))
      count.table <- x$counts
    }
    rm(o)
    count.table
  })
  count.table <- as.matrix(count.table)
  counts(object) <- count.table
  object
}
