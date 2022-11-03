#' constructs time course table for clustering analysis
#'
#' This function constructs a time course table of which rows
#' corrsponding to genomic features and columns the timepoint.
#' values can be mean normalized read counts or log2-fold
#' changes compared to the first timepoint. The time course
#' table is used for clustering analysis.
#'
#' @param object a \code{TCA} object returned by \code{DBanalysis}.
#'
#' @param value character string, either '\code{expression}' or
#' '\code{FC}'. '\code{expression}' is the mean normalized read
#' counts of replicates, \code{FC}' is the log2-fold changes
#' compared to the first time point.
#'
#' @param lib.norm logical indicating whether or not use effective
#' library size (see 'Details' in \code{\link{counts}}).
#'
#' @param norm.method character string specifying the normalization
#' method if \code{value} is '\code{expression}'
#'
#' @param subset optinal character vector giving a subset of
#' genomic features, if not NULL, time course table is generated
#' for only this subset of genomic features.
#'
#' @param filter logical, whether to drop the genomic features
#' shows no significant changes (defined by \code{pvalue},
#' \code{pvalue.threshold},\code{abs.fold} and \code{direction})
#' between any two time points.
#'
#' @param pvalue character string specify the type of p-values
#' ('\code{PValue}' or adjusted p-value 'paj')
#'
#' @param pvalue.threshold a numeric value giving threshold of
#' selected p-value, only features with higher (ajusted) p-values
#' than the threshold are kept.
#'
#' @param abs.fold a numeric value, the least absolute log2-fold
#' changes
#'
#' @param direction character string specify the direction of fold
#' changes ('\code{up}' (positive fold changes), \code{down}'
#' (negative fold changes), \code{both}'(both positive and negative
#' fold changes)), features changes more than \code{abs.fold} in
#' the defined direction are kept.
#'
#' @param ... additional arguments passing to \code{\link{rpkm}},
#' \code{\link{cpm}}
#' @note
#' If '\code{expression}' in \code{value} is chosen, for replicates ,
#' the normalized expression value is first calculated for each
#' replicate, then mean value is taken to represent the normalized
#' expression value.
#'
#' @return
#' A \code{TCA} object
#'
#' @author
#' Mengjun Wu
#'
#' @examples
#' data(tca_ATAC)
#' tca_ATAC <- DBanalysis(tca_ATAC)
#' tca_ATAC <- timecourseTable(tca_ATAC, value = 'expression',
#'                             lib.norm = TRUE, norm.method = 'rpkm')
#'
#' @export
#'
#'
timecourseTable <- function(object, value = "expression", lib.norm = TRUE,
                            norm.method = "rpkm", subset = NULL,
                            filter = FALSE, pvalue = "fdr",
                            pvalue.threshold = 0.05, abs.fold = 2,
                            direction = "both", ...) {
  if (!value %in% c("expression", "FC")) {
    err <- paste0("The value of time course table should be either normalized expression table (value=\"expression\") or logarithm of fold changes (value=\"FC\")")
    stop(err)
  }
  group <- unique(object@design$timepoint)
  genointerval <- object@genomicFeature[object@genomicFeature$id %in%
                                          row.names(object@DBfit$counts), ]
  if (value == "expression") {
    count <- object@DBfit$counts
    if (lib.norm) {
      y <- DGEList(counts = count, group = object@design$timepoint)
      y <- calcNormFactors(y)
    } else {
      y <- DGEList(counts = count, group = object@design$timepoint)
    }
    if (!norm.method %in% c("rpkm", "cpm")) {
      err <- paste0("norm.method should be one of \"rpkm\" or \"cpm\".")
      stop(err)
    }
    tc <- switch(norm.method, rpkm = {
      giwidth <- genointerval$end - genointerval$start
      t <- rpkm(y, normalized.lib.size = lib.norm, gene.length = giwidth, ...)
      t
    }, cpm = {
      t <- cpm(y, normalized.lib.size = lib.norm, ...)
      t
    })
    tc <- data.frame(tc, stringsAsFactors = FALSE)
    colnames(tc) <- object@design$timepoint
    tc <- as.data.frame(sapply(unique(names(tc)), function(col) rowMeans(tc[names(tc) == col])))
  }
  if (value == "FC") {
    group1 <- group[1]
    group2 <- group[group != group1]
    tc <- matrix(0, nrow = dim(genointerval)[1])
    t <- DBresult(object, group1 = group1, group2 = group2,
                  top.sig = FALSE, result.type = "list")
    t <- as(t, "list")
    for (i in t) {
      tc <- cbind(tc, i$logFC)
    }
    colnames(tc) <- group
    rownames(tc) <- genointerval$id
  }
  tc <- as.matrix(tc)

  if (filter) {
    contrasts <- colnames(object@contrasts)
    if (pvalue == "PValue") {
      p <- "none"
      p2 <- "PValue"
    } else {
      p <- pvalue
      p2 <- "paj"
    }
    DBtmpfilter <- DBresult(object, contrasts = contrasts,
                            p.adjust = p, result.type = "list", 
                            pvalue.threshold = pvalue.threshold, 
                            abs.fold = abs.fold,
                            top.sig = TRUE)
    feature.filter <- c()
    for (i in DBtmpfilter) {
      feature.filter <- c(feature.filter, rownames(i))
    }
    tc <- tc[unique(feature.filter), ]
  }

  if (!is.null(subset)) {
    tc <- tc[row.names(tc) %in% subset, ]
  }

  object@tcTable <- tc
  object
}


