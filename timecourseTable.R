#' constructs time course table for clustering analysis
#'
#' This function constructs a time course table of which rows
#' are genomic features and columns time points.
#' values can be normalized expression levels or log2-fold
#' changes compared to a control time point. The time course
#' table is used for clustering analysis.
#'
#' @param object a \code{TCA} object returned by \code{DBanalysis}.
#'
#' @param value a character string, either "\code{expression}" or
#' "\code{FC}". "\code{expression}" is the mean normalized read
#' counts of replicates, "\code{FC}" is the log2-fold changes
#' compared to the first time point.
#' 
#' @param control.group a character string giving the time point to 
#' be compared with, i.e., the denominator in the fold changes. It 
#' should match one of the time points in the \code{design} table 
#' in the \code{TCA} object. 
#' 
#' @param lib.norm logical indicating whether or not use effective
#' library size (see "Details" in \code{\link{counts}}).
#'
#' @param norm.method a character string specifying the normalization
#' method if \code{value} is "\code{expression}"
#'
#' @param subset an optinal character vector giving a subset of
#' genomic features, if not NULL, time course table is generated
#' for only this subset of genomic features.
#'
#' @param filter logical, whether to drop the genomic features
#' shows no significant changes (defined by \code{pvalue},
#' \code{pvalue.threshold},\code{abs.fold} and \code{direction})
#' between any two time points.
#'
#' @param pvalue character string specify the type of p-values:
#' "\code{none}" is unadjusted p-value or one of adjusted p-value 
#' "\code{holm}", "\code{hochberg}", "\code{hommel}", "\code{bonferroni}", 
#' "\code{BH}", "\code{BY}", "\code{fdr}". 
#'
#' @param pvalue.threshold a numeric value giving threshold of
#' selected p-value, significant changes have lower
#' (adjusted) p-values than the threshold.
#'
#' @param abs.fold a numeric value, the least minimum log2-fold
#' changes. The returned genomic regions have changes 
#' with absolute log2-fold changes exceeding \code{abs.fold}.
#'
#' @param direction character string specify the direction of fold
#' changes. "\code{up}": positive fold changes; "\code{down}":
#' negative fold changes; "\code{both}": both positive and
#' negative fold changes.
#'
#' @param ... additional arguments passing to \code{\link{rpkm}},
#' \code{\link{cpm}}
#' @note
#' If "\code{expression}" in \code{value} is chosen, the average 
#' normalized expression values of replicates for each group will 
#' be calculated and returned.
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
timecourseTable <- function(object, value = "expression", control.group=NULL,
                            lib.norm = TRUE, norm.method = "rpkm", 
                            subset = NULL, filter = FALSE, pvalue = "fdr",
                            pvalue.threshold = 0.05, abs.fold = 2,
                            direction = "both", ...) {
  if (!value %in% c("expression", "FC")) {
    err <- paste0("The value of time course table should be either 
                  normalized expression table (value=\"expression\") or 
                  logarithm of fold changes (value=\"FC\")")
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
    tc <- NULL
    if(is.null(control.group)){
      err <- paste0("control group needs to be specified.")
      stop(err)
    }
    group1 <- control.group
    tc <- cbind(tc, rep(0, length(genointerval[, 1])))
    group2 <- group[group != group1]
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


