#' Perform differential binding analysis
#'
#' This function performs differetial analysis by fitting read counts
#' to a negative binomial generalized linear model.
#'
#' @param object a \code{TCA} object.
#'
#' @param categories character string indicating levels of which factor (column
#' in the \code{design} slot) are compared in the differential analysis. For time
#' course analysis, the default factor is '\code{timepoint}'.
#'
#' @param norm.lib logical indicating whether or not use effective library size when perform
#' normalization. See 'Details' of \code{\link{counts}}
#'
#' @param filter.type character string indicating which type of count (raw or normalized) is used
#' when doing filtering. Options are '\code{raw}', '\code{cpm}', '\code{rpkm}', '\code{NULL}'.
#' '\code{NULL}' means no filtering will be performed.
#'
#' @param filter.value A numberic value; if values of selected \code{filter.type} ('\code{raw}',
#' '\code{cpm}', '\code{rpkm}') of a genomic feature are larger than the \code{filter.value} in
#' at least a certain number (\code{samplePassfilter}) of samples/libraries for any of the conditions,
#' such genomic feature will be kept; otherwise the genomic feature will be dropped
#'
#' @param samplePassfilter numberic value indicating the least number of samples/libraries a genomic
#' feature with counts (raw or normalized) more than \code{filter.value} for all conditions if such
#' genomic feature will be kept.
#'
#' @param ... additional arguments passed to \code{\link{glmFit}} from edgeR package.
#'
#' @details The differetial event is detected by using the generalized linear model (GLM) methods
#' (McCarthy et al, 2012). This function fits the read counts of each genes to a negative binomial
#' glms by using \code{\link{glmFit}} function from edgeR. To further test the significance of changes,
#' see \code{DBresult}, \code{TopDBresult}
#'

#' @return
#' A \code{TCA} object
#'
#' @author
#' Mengjun Wu, Lei Gu
#'
#' @references McCarthy,D.J.,Chen, Y., & Smyth, G. K.(2012). Differential expression analysis of multifactor RNA-Seq
#' experiments with respect to biological variation. Nucleic acids research 40, 4288-4297.
#'
#' @seealso \code{DBresult}, \code{TopDBresult}
#'
#' @examples
#' data(tca_ATAC)
#' tca_ATAC <- DBanalysis(tca_ATAC)
#'
#' @export
DBanalysis <- function(object, categories = "timepoint", norm.lib = TRUE, filter.type = NULL, filter.value = NULL,
                       samplePassfilter = 2, ...) {
  if (!categories %in% colnames(object@design)) {
    err <- paste0("Can not find ", categories, " in design, please check if the correspoinding field is missing or a different name is used.")
    stop(err)
  }

  object@contrasts <- contrastMatrix(object, categories)

  # require(edgeR)
  group <- object@design[[categories]]
  y <- DGEList(counts = object@counts, group = group)
  if (norm.lib) {
    y <- calcNormFactors(y)
  }
  if (!is.null(filter.type)) {
    if (is.null(filter.value)) {
      err <- paste0("\"filter.value\" is required to be specified for the chosen filter.type ", filter.type,
                    ".")
      stop("\"filter.value\" is required to be specified for chosen .")
    } else {
      y <- switch(filter.type, raw = {
        ind <- rowSums(y$counts > filter.value) >= samplePassfilter
        y <- y[ind, , keep.lib.sizes = FALSE]
        y
      }, cpm = {
        ind <- rowSums(cpm(y, ...) > filter.value) >= samplePassfilter
        y <- y[ind, , keep.lib.sizes = FALSE]
        y

      }, rpkm = {
        giwidth <- object@genomicFeature$end - object@genomicFeature$start
        ind <- rowSums(rpkm(y, gene.length = giwidth, ...) > filter.value) >= samplePassfilter
        y <- y[ind, , keep.lib.sizes = FALSE]
        y
      })
    }
  }

  design <- model.matrix(~0 + group, data = y$samples)
  colnames(design) <- levels(y$samples$group)
  design <- design[, unique(group)]
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design, ...)
  object@DBfit <- fit
  object
}

# initialize a contrast table with all possible comibinations of group in defined categories
contrastMatrix <- function(object, categories) {
  ca <- unique(object@design[[categories]])
  a <- length(ca)
  b <- 2 * choose(a, 2)
  contrastM <- matrix(0, a, b)
  name <- c()
  count <- 1
  count.col <- -1
  count.col2 <- 0
  for (i in 1:(a - 1)) {
    count = count + 1
    for (j in count:a) {
      count.col <- count.col + 2
      count.col2 <- count.col2 + 2
      n <- paste0(ca[j], "vs", ca[i])
      n1 <- paste0(ca[i], "vs", ca[j])
      name <- c(name, n, n1)
      contrastM[i, count.col] = -1
      contrastM[j, count.col] = 1
      contrastM[j, count.col2] = -1
      contrastM[i, count.col2] = 1
    }
  }
  dimnames(contrastM) <- list(ca, name)
  contrastM
}

