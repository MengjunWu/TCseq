#' combine and merge multiple BED files
#'
#' This function merges overlapping genomic regions into a single feature. 
#' The merged single feature represents the widest genomic interval 
#' that covers all overlapping regions.
#'
#' @param data a data frame containg coordinates information of peaks
#' to be merged. Columns of the data frame should be consistent with
#' the BED format where the first column contains chromosome information,
#' the second column the starting position, and the third column 
#' the ending position.
#'
#' @param dir a character string giving the directory where BED files
#' are stored. If \code{data} is not given, the function will reads
#' in the BED files under \code{code}.
#'
#' @param pattern an \code{\link{regular expression}}, only files that
#' have names match the regular expression will be read in.
#'
#' @param merge logical indicating whether to merge overlapped regions
#' or not. If False, regions are simply combined.
#'
#' @param overlap a numberic value giving the least number of base(s)
#' two regions should overlap when merging them.
#'
#' @param ratio a numberic value giving the thresold of overlapping
#' ratio between two regions to merge them. See '\code{Details}' below
#' for the definition of the overlapping ratio.
#'
#' @return a data frame with four columns: \code{chr}, \code{start},
#' \code{stop}, \code{id}
#'
#' @details
#' The overlapping ratio (OR) is defined as:
#'
#' \deqn{ OR = \frac{n}{\min(length(a), length(b)}}
#'
#' \eqn{a}, \eqn{b} are two genomic regions, \eqn{n} is the number of
#' overlapping bases between region \eqn{a} and region \eqn{b}.
#'
#' @author
#' Mengjun Wu, Lei Gu
#'
#' @examples
#' peaks <- data.frame(chr = c(rep('chr1',4),rep('chr2', 3), rep('chr3',2)),
#'                     start = c(100,148,230,300,330,480,1000,700,801),
#'                     end = c(150,220,500,450,600,900,1050,760,900))
#'
#' merged_peaks <- peakreference(data = peaks, merge = TRUE, overlap = 1)
#'
#' @export

peakreference <- function(data = NULL, dir = NULL, pattern = NULL,
                          merge = TRUE, overlap = 1, ratio = NULL) {
  if (is.null(data) && is.null(dir)) {
    stop("Either a data.frame of genomic coordinates or a directory 
         for the BED files should be given")
  }
  if (!is.null(data)) {
    checkBEDformat(data)
    data[, 1] <- factor(data[, 1])
    peakset <- data
  }
  if (is.null(data) && !is.null(dir)) {
    old <- setwd(tempdir())
    on.exit(setwd(old), add = TRUE)
    setwd(dir)
    filenames <- list.files(pattern = pattern)
    if (length(filenames) == 0) {
      err <- paste0("Can not find file names containing '",
                    pattern, "'.")
      stop(err)
    }
    datalist <- lapply(filenames, function(x) {
      read.table(file = x, header = FALSE)
    })
    peakset <- do.call(rbind, datalist)
    checkBEDformat(peakset)
  }
  peakset <- peakset[order(peakset[, 1], peakset[, 2]), ]
  if (merge) {
    if (overlap <= 0 || round(overlap) != overlap) {
      stop("\"overlap\" must be integer and greater than 0.")
    }
    peakset.sub <- split(peakset, peakset[, 1],
                         drop = TRUE)
    level <- names(peakset.sub)
    mergedpeak <- c()
    for (i in seq_len(length(peakset.sub))) {
      temp <- peakset.sub[[i]]
      if (is.null(ratio)) {
        submerge <- intervalmerge(temp[, 2], temp[, 3],
                                  overlap = overlap)
      } else {
        submerge <- intervalmerge(temp[, 2], temp[, 3],
                                  ratio = ratio)
      }

      chr <- rep(level[i], length(submerge[, 1]))
      submerge1 <- data.frame(chr, submerge)
      mergedpeak <- rbind(mergedpeak, submerge1)
    }
    name <- paste0("peak", seq_len(length(mergedpeak[, 1])))
    mergedpeak <- data.frame(mergedpeak, name)
    colnames(mergedpeak) <- c("chr", "start", "end", "id")
    mergedpeak
  } else {
    peakset
  }

}


checkBEDformat <- function(data) {
  if (ncol(data) < 3) {
    stop("At least three columns should be provided. The first column contains chromosome name,
         the second column contains starting position, the third column contains ending position.")
  }
  if (class(as.vector(data[, 1])) != "character") {
    stop("The first column contains chromosome name and must be character.")
  }
  if (any(round(data[,2]) != data[,2]) &&
      any(round(data[,3]) != data[,3])) {
    stop("the second and third column contain starting and ending positions, must be numeric.")
  }
  }

intervalmerge <- function(a0, b0, overlap = NULL,
                          ratio = NULL) {
  if (length(a0) > 1) {
    a1 <- c(a0[1])
    b1 <- c(b0[1])
    merge <- NULL
    for (i in seq_len(length(a0) - 1)) {
      if (is.null(ratio)) {
        if (b1[length(b1)] - a0[i + 1] < overlap) {
          a1 <- append(a1, a0[i + 1])
          b1 <- append(b1, b0[i + 1])
        } else {
          b1[length(b1)] <- max(b1[length(b1)], b0[i + 1])
        }
      }
      if (is.null(overlap)) {
        len <- min((b1[length(b1)] - a1[length(b1)]),
                   (b0[i + 1] - a0[i + 1]))
        rt <- (b1[length(b1)] - a0[i + 1])/len
        if (rt < ratio) {
          a1 <- append(a1, a0[i + 1])
          b1 <- append(b1, b0[i + 1])
        } else {
          b1[length(b1)] <- max(b1[length(b1)], b0[i + 1])
        }
      }
    }
    merge <- cbind(a1, b1)
  }
  if (length(a0) <= 1) {
    a1 <- c(a0[1])
    b1 <- c(b0[1])
    merge <- cbind(a1, b1)
  }
  merge
}








