#' @import ggplot2
#' @import edgeR
#' @import e1071
#' @import cluster
#' @import reshape2
#' @import grid
#' @import locfit
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom BiocGenerics counts counts<- design
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools BamFile BamFileList
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom grDevices rainbow
#' @importFrom stats as.dist complete.cases cor cutree hclust kmeans model.matrix sd time
#' @importFrom utils capture.output read.table
#' @importFrom methods new validObject
NULL

counts.TCA <- function(object, normalization = "none", lib.norm = TRUE,
                       log = FALSE, ...) {
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
      t <- rpkm(y, normalized.lib.sizes = lib.norm, gene.length = giwidth,
                log = log, ...)
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
#' \code{counts} extract raw read counts  stored in a \code{TCA} object 
#' or compute normalized counts from the raw counts.
#'
#' @name counts
#' @aliases counts counts,TCA-method counts<-,TCA-method
#' @param object a \code{TCA} object.
#'
#' @param normalization character string giving the normalization method.
#' Options are "\code{none}" (original raw counts), "\code{cpm}" (counts
#' per million),
#' "\code{rpkm}" (reads per kilobase per million).
#'
#' @param lib.norm logical indicating whether or not use effective library
#' size (see Details below) when \code{normalization} is "\code{cpm}" or
#' "\code{rpkm}".
#'
#' @param log logical if \code{TRUE}, the returned value will be on a log2
#' scale.
#'
#' @param value an integer matrix.
#'
#' @param ... additional arguments passed to \code{\link{cpm}} or
#' \code{\link{rpkm}} in the edgeR package.
#'
#' @details when calculating normalized counts, library size can be rescaled
#' to minimize the log-fold changes between samples for most genomic features
#' (e.g. genes, binding sites) by multiplying a scale factor. The rescaled
#' library size is called effective library size. In this function, the scale
#' factor is calculated using the weighted trimmed mean of M-values (TMM,
#' Robinson et al (2010))
#'
#' If log2 values are computed, a small count would be added to avoid logarithm 
#' of zero. The actual added count will be scaled according to the library size,
#' for details see \code{\link{addPriorCount}} in the edgeR package
#' when not specified, the prior count is set to 0.25 by default.
#'
#' @references
#' Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for
#' differential expression analysis of RNA-seq data. Genome biology, 11(3), 1.
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
#' Accessors are provided to extract \code{design}, \code{genomicFeature},
#' \code{tcTable}, \code{clustResults} slots of a TCA class. The \code{design}
#' slot stores experimental information of samples/libraries, the
#' \code{genomicFeature} slot stores genomic coordinates of features, the
#' \code{tcTable} slot stores time couse data as a matrix, where rows are
#' genomic features and columns time points. The \code{clustResults} slot
#' stores results of clustering analysis as a \code{clust} object.
#'
#' @name TCA.accessors
#' @aliases design design,TCA-method genomicFeature,TCA-method
#' tcTable,TCA-method clustResults,TCA-method
#'
#' @param object \code{TCA} object object
#' @return
#' \code{design} returns a data frame. \code{genomicFeature} returns a data frame.
#' \code{tcTable} returns a numeric matrix. \code{clustResults} returns a
#' \code{clust} object, see \code{\link{clust}} for details.
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

#' @rdname TCA.accessors
#' @export
setMethod("design", "TCA", function(object) {
  object@design
})

#' @rdname TCA.accessors
#' @export

setGeneric("genomicFeature", function(object) standardGeneric("genomicFeature"))
setMethod("genomicFeature", "TCA", function(object) {
  object@genomicFeature
})

#' @rdname TCA.accessors
#' @export

setGeneric("tcTable", function(object) standardGeneric("tcTable"))

#' @rdname TCA.accessors
#' @export

setMethod("tcTable", "TCA", function(object) {
  object@tcTable
})

#' @rdname TCA.accessors
#' @export

setGeneric("clustResults", function(object) standardGeneric("clustResults"))

#' @rdname TCA.accessors
#' @export

setMethod("clustResults", "TCA", function(object) {
  object@clusterRes
})

#' Accessors to extract slots of a clust class.
#'
#' Accessors are provided to extract \code{data}, \code{centers}, \code{cluster}, 
#' and \code{membership} slots stored in a clust class.
#' @name clust.accessors
#' @aliases clustData clustData,clust-method clustCenters,clust-method
#' clustCluster,clust-method clustMembership,clust-method
#'
#' @param object \code{clust} object object
#' @return
#' \code{clustData} returns a data matrix. \code{clustCenters} returns a matrix of
#' centers. \code{clustCluster} returns an integer vector. \code{clustMembership}
#' returns a matrix of membership, see \code{\link{clust}} for details.
#'
#' @author
#' Mengjun Wu
#'
#' @seealso
#' \code{\link{clust}}

#' @rdname clust.accessors
#' @export
setGeneric("clustData", function(object) standardGeneric("clustData"))

#' @rdname clust.accessors
#' @export
setMethod("clustData", "clust", function(object) {
  object@data
})

#' @rdname clust.accessors
#' @export
setGeneric("clustCenters", function(object) standardGeneric("clustCenters"))

#' @rdname clust.accessors
#' @export

setMethod("clustCenters", "clust", function(object) {
  object@centers
})

#' @rdname clust.accessors
#' @export
setGeneric("clustCluster", function(object) standardGeneric("clustCluster"))

#' @rdname clust.accessors
#' @export

setMethod("clustCluster", "clust", function(object) {
  object@cluster
})

#' @rdname clust.accessors
#' @export
setGeneric("clustMembership", function(object) standardGeneric("clustMembership"))

#' @rdname clust.accessors
#' @export

setMethod("clustMembership", "clust", function(object) {
  object@membership
})
