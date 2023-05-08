#' Perform differential expression analysis
#'
#' This function is a wrapper for the \code{\link{glmFit}} in edgeR package.
#'
#' @param object a \code{TCA} object.
#'
#' @param categories character string giving which column in \code{design} 
#' will be used for differential analysis. For time course analysis, the default
#' column is "\code{timepoint}".
#'
#' @param norm.lib logical indicating whether or not use effective
#' library size when perform normalization. See \code{\link{counts}} for more 
#' details.
#'
#' @param filter.type character string indicating which type of count
#' (raw or normalized) is used when performing filtering. Options are
#' "\code{raw}", "\code{cpm}", "\code{rpkm}", "\code{NULL}". No filtering will 
#' be performed when using "\code{NULL}'.
#'
#' @param filter.value a numberic value; minimum values of selected
#' \code{filter.type} ("\code{raw}", "\code{cpm}", "\code{rpkm}"). It is used in 
#' combination with \code{samplePassfilter}.
#'
#' @param samplePassfilter a numberic value indicating the minimum number
#' of samples/libraries in which a genomic feature has counts value 
#' (raw or normalized) more than \code{filter.value}. Smaller than this number, 
#' the genomic feature will be filtered out.
#'
#' @param ... additional arguments passed to \code{\link{glmFit}} from
#' \code{edgeR} package.
#'
#' @details The differetial event is detected by using the generalized
#' linear model (GLM) methods (McCarthy et al, 2012). This function
#' fits the read counts of each genes to a negative binomial glm by
#' using \code{\link{glmFit}} function from edgeR. To further test the
#' significance of changes, see \code{DBresult}, \code{TopDBresult}
#'

#' @return
#' A \code{TCA} object
#'
#' @author
#' Mengjun Wu, Lei Gu
#'
#' @references McCarthy,D.J.,Chen, Y., & Smyth, G. K.(2012). Differential
#' expression analysis of multifactor RNA-Seq experiments with respect to
#' biological variation. Nucleic acids research 40, 4288-4297.
#'
#' @seealso \code{DBresult}, \code{TopDBresult}
#'
#' @examples
#' data(tca_ATAC)
#' tca_ATAC <- DBanalysis(tca_ATAC)
#'
#' @export
DBanalysis <- function(object, categories = "timepoint", norm.lib = TRUE,
                       filter.type = NULL, filter.value = NULL,
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
      err <- paste0("\"filter.value\" is required to be specified for the chosen filter.type ",
                    filter.type, ".")
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
  name <- vector(mode = "character", length = b)
  count <- 1
  count.col <- -1
  count.col2 <- 0
  for (i in seq_len((a - 1))) {
    count = count + 1
    for (j in count:a) {
      count.col <- count.col + 2
      count.col2 <- count.col2 + 2
      n <- paste0(ca[j], "vs", ca[i])
      n1 <- paste0(ca[i], "vs", ca[j])
      name[count.col] <- n
      name[count.col2] <- n1
      contrastM[i, count.col] = -1
      contrastM[j, count.col] = 1
      contrastM[j, count.col2] = -1
      contrastM[i, count.col2] = 1
    }
  }
  dimnames(contrastM) <- list(ca, name)
  contrastM
}

