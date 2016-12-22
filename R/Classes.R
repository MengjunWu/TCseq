#' clust class
#'
#'\code{clust} is a S4 class for storing results of a clustering
#'analysis for time course data.
#'
#'@section Slots:
#'Oject of this class contains the following slots:
#'\describe{
#'  \item{\code{method}}{clustering method that has been used}
#'  \item{\code{dist}}{distance method that has been used}
#'  \item{\code{data}}{a matrix of original or standardized data that
#'  has been used for the analysis}
#'  \item{\code{centers}}{a matrix of class centers}
#'  \item{\code{cluster}}{an integer vector of length \eqn{n} (\eqn{n}
#'  is the number of data points each integer indicates the cluster a
#'  data point belongs to. For the fuzzy cmeans clustering method, a
#'  data point is assigned to the closest cluster to which the data
#'  point has highest membership value.}
#'  \item{\code{membership}}{a matrix with membership values of the
#'  data points to all the clusters}
#'}
#'@details:
#'The clust objects are returned from \code{\link{timeclust}} and have
#'a show method printing a compact summary of their contents
#'
#'@author Mengjun Wu
#'
#'@seealso \code{\link{timeclust}}, \code{\link{@}}
#'@exportClass clust

clust <- setClass("clust", slots = c(method = "character",
                                     dist = "character",
                                     data = "matrix",
                                     centers = "matrix",
                                     cluster = "integer",
                                     membership = "matrix"))

#'@rdname TCA
#'@export
setClass("TCA", slots = c(design = "data.frame", counts = "matrix",
                          genomicFeature = "data.frame",
                          DBfit = "DGEGLM",contrasts = "matrix",
                          tcTable = "matrix", clusterRes = "clust"),
         prototype = list(counts = matrix(0L, 0L, 0L),
                          design <- data.frame()))

setValidity("TCA", function(object) {
  counts <- object@counts
  design <- object@design
  genomicFeature <- object@genomicFeature
  if (!is.numeric(counts)) {
    stop("All counts must be numeric.")
  }
  if (any(is.na(counts))) {
    stop("NA values are not allowed in counts. ")
  }
  if (any(counts < 0)) {
    stop("counts contain negative number(s), all counts must be positive")
  }
  if (!is.integer(counts)) {
    if (any(round(counts) != counts)) {
      stop("All counts must be intergers.")
    } else {
      mode(counts) <- "integer"
      warning("All counts are coerced to integers.")
    }
  }
  if (!identical(matrix(0L, 0L, 0L), counts)) {
    if (ncol(counts) != nrow(design)) {
      stop("Number of columns in 'counts' must equal to number of rows in 'design'.")
    }
    if (nrow(counts) != nrow(genomicFeature)) {
      stop("Number of rows in 'counts' must equal to number of rows in 'genomicFeature'")
    }
  }
  if (!sum(c("sampleid", "timepoint", "group") %in%
           tolower(colnames(design))) == 3) {
    err <- paste0("One or more following required fields in 'design' are missing: 'sampleid', 'timepoint', 'group', check if the columns are correctly named or if the corresponding information is provided.")
    stop(err)
  }
  if (!sum(c("id", "chr", "start", "end") %in%
           tolower(colnames(genomicFeature))) == 4) {
    err <- paste0("One or more following required fields in 'genomicFeature' are missing: 'id', 'chr', 'start','end', check if the columns are correctly named or if the corresponding information is provided.")
    stop(err)
  }
  TRUE
})

#'TCA class and constructor
#'
#'\code{TCA} is a S4 class for storing input data, results of
#'differential binding and clustering analysis. A \code{TCA} object
#'can be created by the constructor function from a table of sample
#'information, a table genomic coordinates of features, read
#'counts(optional).
#'
#'@param design a data frame containing information about
#'samples/libraries, For time course analysis, design should contain
#'at least three columns (case insensitive): \code{sampleid},
#'\code{timepoint} and \code{group} providing time point and group
#'information of each sample/library. If \code{counts} is not provided
#'when creating \code{TCA} object, the column \code{BAMfile} can be
#'included in the \code{design} data frame, providing corresponding
#'BAM filename of each sample/library, this information can be used
#'for generating count table by using
#'\code{\link{countReads}} function.
#'
#'@param counts an integer matrix containing read counts. Rows
#'correspond to genomic features and columns to samples/libraries.
#'
#'@param genomicFeature a data frame or a GRanges object containing
#'genomic coordinates of features of interest (e.g. genes in RNA-seq,
#'binding regions in ChIP-seq). If genomicFeature is a data frame,
#'four columns are required in \code{genomicFeature}: \code{id},
#'\code{chr}, \code{start}, \code{end}; if genomicFeature is a Granges
#'object, the metadata column "\code{id}" is required. For
#'\code{TCAFromSummarizedExperiment}, genomicFeature must be
#'provided if se is a SummarizedExperiment object.
#'
#'
#'@param se A SummarizedExperiment or a RangedSummarizedExperiment
#'object. The object might contain multiple assays (count table)
#'in the assay list, only the first one will be taken to construct
#'TCA object. For SummarizedExperiment object, \code{genomicFeature}
#'must be provided while for RangedSummarizedExperiment object,
#'the genomic features will be extracted directly from the object.
#'
#'@param zero.based Logical. If TRUE, the start positions of the
#'genomic ranges in the returned \code{TCA} object are \emph{0-based},
#'if FALSE, the start positions will be \emph{1-based}.
#'
#'@return A TCA object
#'
#'@details A TCA object can be created without providing read counts,
#'read counts can be provided by \code{\link{counts}} or generated by
#'\code{\link{countReads}}, the number of rows should equal to that in
#'\code{genomicFeature} and the number of columns should equal to number
#'of rows in \code{design}. Input data and analysis results in a TCA
#'object can be accessed by using corresponding accessors and functions.
#'The TCA objects also have a show method printing a compact summary of
#'their contents see \code{\link{counts}}, \code{\link{TCA.accessors}},
#'\code{\link{DBresult}}, \code{\link{tcTable}}, \code{\link{timeclust}}.
#'\code{clust}
#'@author Mengjun Wu
#'@seealso \code{\link{counts}}, \code{\link{TCA.accessors}},
#'\code{\link{DBresult}}, \code{\link{timeclust}}, \code{\link{clust}}
#'
#'@author Mengjun Wu
#'
#'@examples
#'#create data frame of experiment design: 4 time points and 2 replicates for each time point.
#'d <- data.frame(sampleID = 1:8, group = rep(c(1, 2, 3, 4), 2),
#'                timepoint = rep(c('0h', '24h', '48h', '72h'), 2))
#'
#'
#'#create data frame of genomic intervals of interest
#'gf <- data.frame(chr = c(rep('chr1', 3), rep('chr2', 2), rep('chr4', 2)),
#'                 start = seq(100, 2000, by = 300),
#'                 end = seq(100, 2000, by = 300) + 150,
#'                 id = paste0('peak', 1:7))
#'tca <- TCA(design = d, genomicFeature = gf)
#'genomicFeature(tca)
#'
#'#if count table is available
#'c <- matrix(sample(1000, 56), nrow = 7, dimnames = list(paste0('peak', 1:7), 1:8))
#'tca <- TCA(design = d, counts = c, genomicFeature = gf)
#'# replace the count table of a \code{TCA} object
#'c2 <- matrix(sample(500, 56), nrow = 7, dimnames = list(paste0('peak', 1:7), 1:8))
#'counts(tca) <- c2
#'
#'
#'@export
TCA <- function(design, counts = matrix(0L, 0L, 0L), genomicFeature,
                zero.based = TRUE) {

  if (!is.numeric(counts)) {
    stop("All counts must be numeric.")
  }
  if (any(is.na(counts))) {
    stop("NA values are not allowed in counts.")
  }
  if (any(counts < 0)) {
    stop("counts contain negative number(s), all counts must be positive")
  }
  if (!is.integer(counts)) {
    if (any(round(counts) != counts)) {
      stop("All counts must be intergers.")
    } else {
      mode(counts) <- "integer"
      warning("All counts are coerced to integers.")
    }
  }
  if (!is.data.frame(design)) {
    stop("design must be 'data.frame'.")
  }
  if (!is.data.frame(genomicFeature) &&
      !is(genomicFeature, "GRanges")) {
    stop("genomicFeature must be a data frame or a GRanges object.")
  }
  if (is.data.frame(genomicFeature)) {
    if (sum(c("id", "chr", "start", "end") %in%
            tolower(colnames(genomicFeature))) != 4) {
      err <- paste0("One or more following required fields in genomicFeature are missing: 'id', 'chr', 'start','end', check if the columns are correctly named or if the corresponding information is provided. ")
      stop(err)
    }
    if (sum(c("id", "chr", "start", "end") %in%
            colnames(genomicFeature)) != 4) {
      colnames(genomicFeature) <- tolower(colnames(genomicFeature))
      warning("colnames of genomicFeature are all forced to lowercase.")
    }
    if (!zero.based) {
      genomicFeature$start <- genomicFeature$start + 1
    }
  }
  if (is(genomicFeature, "GRanges")) {
    if(!"id" %in% tolower(colnames(elementMetadata(genomicFeature)))) {
      stop("Required metadata of genomicFeature is mising: 'id'.")
    }
    if(!"id" %in% colnames(elementMetadata(genomicFeature))) {
      colnames(elementMetadata(genomicFeature)) <- tolower(colnames(elementMetadata(genomicFeature)))
      warning("Names of metadata of genomicFeature are all forced to lowercase.")
    }
    genomicFeature <- as(genomicFeature, "data.frame")
    if (zero.based) {
      enomicFeature$start <- genomicFeature$start - 1
    }
  }
  if (!identical(matrix(0L, 0L, 0L), counts)) {
    if (ncol(counts) != nrow(design)) {
      stop("number of columns in 'counts' must equal to number of rows in 'design'.")
    }
    if (nrow(counts) != nrow(genomicFeature)) {
      stop("Number of rows in 'counts' must equal to number of rows in 'genomicFeature'")
    }
  }

  if (!sum(c("sampleid", "timepoint", "group") %in%
           tolower(colnames(design))) == 3) {
    err <- paste0("One or more following required fields in design are missing: 'sampleid', 'timepoint', 'group', check if the columns are correctly named or if the corresponding information is provided.")
    stop(err)
  }
  if (!sum(c("sampleid", "timepoint", "group") %in%
           colnames(design)) == 3) {
    colnames(design) <- tolower(colnames(design))
    warning("colnames of design are all forced to lowercase.")
  }

  object <- new("TCA", design = design, counts = counts,
                genomicFeature = genomicFeature)
  object
}

#'@rdname TCA
#'@export
TCAFromSummarizedExperiment <-function(se, genomicFeature=NULL){
  if (!is(se, "SummarizedExperiment") &&
      !is(se, "RangedSummarizedExperiment")) {
    stop("se must be a SummarizedExperiment or a RangedSummarizedExperiment object.")
  }
  if (is(se, "SummarizedExperiment")) {
    if (is.null(genomicFeature)) {
      stop("genomicFeature must be provided")
    }
    design <- as(colData(se), "data.frame")
    counts <- assay(se,1)
  }
  if (is(se, "RangedSummarizedExperiment")) {
    design <- as(colData(se), "data.frame")
    counts <- assay(se,1)
    genomicFeature <- rowRanges(se)
  }
  object <- TCA(design = design, counts = counts,
                genomicFeature = genomicFeature)
  object
}

#Set inheritance
#The Class LargeDataObject from limma has \code{show} method for objects of the class.

setIs("clust", "LargeDataObject")
setIs("TCA", "LargeDataObject")
