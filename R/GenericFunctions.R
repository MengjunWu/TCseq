#' @import ggplot2
#' @import edgeR
#' @import e1071
#' @import cluster
#' @import clusterSim
#' @import reshape2
#' @import grid
#' @import locfit
#' @import Rcpp
#' @importFrom BiocGenerics counts counts<- design
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment assays
#' @importFrom Rsamtools BamFile BamFileList
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom grDevices rainbow
#' @importFrom stats as.dist complete.cases cor cutree hclust kmeans model.matrix sd time
#' @importFrom utils capture.output read.table
#' @rawNamespace if(.Platform$OS.type != 'windows') {importFrom(Rsubread, featureCounts)}
NULL

counts.TCA <- function(object, normalization = "none", lib.norm = TRUE, log = FALSE, ...) {
    if (!normalization %in% c("none", "rpkm", "cpm")) {
        stop("'normalization method should one of 'none', 'rpkm', 'cpm'.")
    }
    if (normalization == "none") {
        t <- object@counts
    }
    if (normalization != "none") {
        genomicFeature <- object@genomicFeature
        group <- object@design$group
        y <- DGEList(counts = object@counts, group = group)
        if (lib.norm) {
            y <- calcNormFactors(y)
        }
        c <- switch(normalization, rpkm = {
            giwidth <- genomicFeature$end - genomicFeature$start
            t <- rpkm(y, normalized.lib.sizes = lib.norm, gene.length = giwidth, log = log, ...)
            t
        }, cpm = {
            t <- cpm(y, normalized.lib.sizes = lib.norm, log = log, ...)
            t
        })
    }
    t
}

#' Extracts counts of a TCA object.
#'
#' \code{counts} extract raw read counts stored in a \code{TCA} object or compute normalized counts.
#'
#' @name counts
#' @aliases counts counts,TCA-method counts<-,TCA-method
#' @param object a \code{TCA} object
#'
#' @param normalization character string giving the normalization method. Options
#' are '\code{none}' (original raw counts), '\code{cpm}' (counts per million),
#' '\code{rpkm}' (reads per kilobase per million).
#'
#' @param lib.norm logical indicating whether or not use effective library size
#' (see 'Details' below) when \code{normalization} is '\code{cpm}' or '\code{rpkm}'.
#'
#' @param log logical if \code{TRUE}, the returned value will be on a log2 scale.
#'
#' @param value an integer matrix
#'
#' @param ... additional arguments passed to \code{\link{cpm}} or \code{\link{rpkm}}
#'
#' @details when calculating normalized counts, library size can be rescaled to
#' minimize the log-fold changes between samples for most genomic features
#' (e.g. genes, binding sites) by multiplying a scale factor. The rescaled
#' library size is called effective library size. In this function, the scale
#' factor is calculated using the weighted trimmed mean of M-values (TMM, Robinson
#' et al (2010))
#'
#' If log2 values are computed, a small count would be added to avoid logarithm of zero.
#' a small count is set proportional to the library size, the average value of such small
#' counts of all libraries counts is set to 0.25 by default.
#'
#' @references
#' Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome biology, 11(3), 1.
#'
#' @return
#' An integer matrix
#'
#' @author
#' Mengjun Wu
#'
#' @examples
#' data(tca_ATAC)
#' c <- counts(tca_ATAC)
#' # normalized counts table
#' c_norm <- counts(tca_ATAC, normalization='rpkm')
#' @export
setMethod("counts", "TCA", counts.TCA)

#' @rdname counts
#' @exportMethod 'counts<-'
setMethod("counts<-", "TCA", function(object, value) {
    object@counts <- value
    validObject(object)
    object
})


#' Accessors to extract slots of a TCA class.
#'
#' The \code{design} slot stores experimental information of samples/libraries,
#' the \code{genomicFeature} slot stores genomic coordinates of features, the \code{tcTable}
#' slot stores  time couse data as a matrix, where rows are genomic features and columns time points.
#' the \code{clustResults} slot stores results of clustering analysis as a \code{clust} object.
#'
#' @name Accessors
#' @aliases design design,TCA-method genomicFeature,TCA-method tcTable,TCA-method clustResults,TCA-method
#' @param object \code{TCA} object object
#' @return
#' \code{design} returns a data frame. \code{genomicFeature} returns a data frame. \code{tcTable}
#' returns a numeric matrix. \code{clustResults} returns a \code{clust} object, see \code{\link{clust}}
#' for details.
#'
#' @author
#' Mengjun Wu
#'
#' @seealso
#' \code{\link{clust}}
#'
#' @examples
#' data(tca_ATAC)
#' genomicFeature(tca_ATAC)
#' tcTable(tca_ATAC)

#' @rdname Accessors
#' @export
setMethod("design", "TCA", function(object) {
    object@design
})

#' @rdname Accessors
#' @export

setGeneric("genomicFeature", function(object) standardGeneric("genomicFeature"))
setMethod("genomicFeature", "TCA", function(object) {
    object@genomicFeature
})

#' @rdname Accessors
#' @export

setGeneric("tcTable", function(object) standardGeneric("tcTable"))

#' @rdname Accessors
#' @export

setMethod("tcTable", "TCA", function(object) {
    object@tcTable
})

#' @rdname Accessors
#' @export

setGeneric("clustResults", function(object) standardGeneric("clustResults"))

#' @rdname Accessors
#' @export

setMethod("clustResults", "TCA", function(object) {
    object@clusterRes
})
