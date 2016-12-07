#' Extract significant differential events
#' This function extracts genomic regions that have significant differential events by log2-fold changes,p-values or adjusted p-values.
#' Features passing certain thresholds will be kept while other features will be dropped.
#'
#' @param x list of data frames of differential binding results (such as returned value of \code{\link{DBresult}}),
#'  required columns of the data frames are
#' \code{logFC}, \code{PValue}, \code{paj}.
#'
#' @param pvalue character string specify the type of p-values ('\code{PValue}' or adjusted p-value 'paj')
#'
#' @param pvalue.threshold a numeric value giving threshold of selected p-value, only features with higher
#' (ajusted) p-values than the threshold are kept.
#'
#' @param abs.fold a numeric value, the least absolute log2-fold changes
#'
#' @param direction character string specify the direction of fold changes ('\code{up}' (positive fold changes),
#''\code{down}' (negative fold changes), \code{both}'(both positive and negative fold changes)), features changes
#' more than \code{abs.fold} in the defined direction are kept.
#'
#' @return
#'A list of data frames
#'
#' @author
#' Mengjun Wu, Lei Gu
#'
#' @examples
#' data(tca_ATAC)
#' tca_ATAC <- DBanalysis(tca_ATAC)
#' res <- DBresult(tca_ATAC, group1 = '0h', group2 = c('24h', '72h'))
#' top_res <- DBresult.filter(res)
#'
#' @seealso
#' \code{\link{DBresult}}
#'
#' @export
DBresult.filter <- function(x, pvalue = "paj", pvalue.threshold = 0.05, abs.fold = 2, direction = "both") {
  d <- x
  if (abs.fold < 0) {
    err <- paste0("\"abs.fold\" should be postive number.")
    stop(err)
  }
  for (i in 1:length(d)) {
    dt <- d[[i]]

    if (direction == "up") {
      dt <- dt[which(dt$logFC >= abs.fold), ]
    }

    if (direction == "down") {
      dt <- dt[which(dt$logFC <= -abs.fold), ]
    }
    if (direction == "both") {
      if (abs.fold == 0) {
        dt <- rbind(dt[which(dt$logFC >= abs.fold), ], dt[which(dt$logFC < -abs.fold), ])
      } else {
        dt <- rbind(dt[which(dt$logFC >= abs.fold), ], dt[which(dt$logFC <= -abs.fold), ])
      }
    }
    dt <- dt[which(dt[, pvalue] < pvalue.threshold), ]
    d[[i]] <- dt
  }
  d
}
