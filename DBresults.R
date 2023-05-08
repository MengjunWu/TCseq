#' This function tests for differential expression 
#'
#' This function is a wrapper for \code{\link{glmLRT}} in edgeR package. 
#' It performs likelihood ratio tests for given coefficinets contrasts 
#' after fitting read counts to a negative binomial glm by
#' \code{\link{DBanalysis}}. \code{DBresult} also extracts the
#' diffential analysis results of given contrasts at a chosen significance level. 
#' \code{DBresult.cluster} returns similar results but only 
#' contain genomic features belong to a given cluster.
#'
#' @name DBresult
#'
#' @param object a \code{TCA} object, for \code{DBresult},
#' \code{DBanalysis} should already be called on the object;
#' for \code{DBresult.cluster}, both \code{DBanalysis} and
#' \code{timeclust} should be already called.
#'
#' @param group1 character string giving the group to be compared with,
#' i.e., the denominator in the fold changes. group1 can be set NULL and 
#' will be ignored if the comparisons are passed to \code{contrasts}
#'
#' @param group2 a character vetor giving the other groups to 
#' compare with \code{group1}, i.e., the numerator in the fold changes.
#' group2 can be set NULL and will be ignored if the comparisons are 
#' passed to \code{contrasts}
#'
#' @param contrasts a character vector, each string in
#' the vector gives a contrast of two groups with the format
#' "group2vsgroup1", group1 is the denominator level in the fold
#' changes and group2 is the numerator
#' level in the fold changes.
#'
#' @param p.adjust character string specifying a correction method
#' for p-values. Options are "\code{holm}", "\code{hochberg}", 
#' "\code{hommel}", "\code{bonferroni}", "\code{BH}", "\code{BY}", 
#' "\code{fdr}", and "\code{none}". 
#'
#' @param top.sig logical if TRUE, only genomic regions with
#' given log2-fold changes and significance levels (p-value) 
#' will be returned. Log2-fold changes are defined by \code{abs.fold}
#' and \code{direction}; significance levels are defined by \code{pvalue} 
#' and \code{pvalue.threshold}
#'  
#' @param pvalue character string specify the type of p-values
#' used for defining the significance level(\code{PValue}
#' or adjusted p-value \code{paj})
#'
#' @param pvalue.threshold a numeric value giving threshold of
#' selected p-value, Significant changes have lower
#' (adjusted) p-values than the threshold.
#'
#' @param abs.fold a numeric value, the minimum absolute log2-fold
#' changes. The returned genomic regions have changes 
#' with absolute log2-fold changes exceeding \code{abs.fold}.
#'
#' @param direction character string specify the direction of fold
#' changes. "\code{up}": positive fold changes; "\code{down}":
#' negative fold changes; "\code{both}": both positive and
#' negative fold changes.  
#'
#' @param cluster an integer giving the number of cluster from which 
#' genomic features are extracted.
#'
#' @param  cmthreshold a numeric value, this argument is applicable
#' only if \code{cmeans}' clustering method is selected when calling
#' \code{\link{timeclust}} function. if not NULL, the result table of
#' genomic features that belong to the defined \code{cluster} and
#' the membership values to this cluster exceed \code{cmthreshold}
#' are extracted.
#'
#' @param result.type character string giving the data type of return
#' value. Options are "GRangesList" and "list".
#'
#' @details This function uses \code{\link{glmLRT}} from edgeR which
#' perform likelihood ratio tests for the significance of changes.
#' For more deatils,
#' see \code{\link{glmLRT}}
#'
#' @note If not NULL \code{group1}, \code{group2} and \code{contrasts},
#' result tables are extracted from comparisons in \code{constrasts}.
#'
#' @return
#' A list or a GRangesList.
#' If \code{result.type} is "GRangesList", a GRangesList is returned containing
#' the differential analysis results for all provided contrasts. Each GRanges 
#' object of the list is one contrast, the analysis results are contained in 4 
#' metadata columns:
#'
#' @return \code{logFC} log2-fold changes between two groups.
#'
#' @return \code{PValue} p-values.
#'
#' @return \code{paj} adjusted p-values
#'
#' @return \code{id} name of genomic features 
#'
#' If \code{result.type} is "list", a list of data frames is returned.
#' Each data frame contains one contrast with the following columns:
#'
#' @return \code{logFC} log2-fold changes between two groups.
#'
#' @return \code{PValue} p-values.
#'
#' @return \code{paj} adjusted p-values
#'
#' @return \code{chr}  name of chromosomes 
#'
#' @return \code{start} starting positions of features in the 
#' chromosomes
#'
#' @return \code{end} ending postitions of features in the chromosomes
#'
#' @return \code{id} name of genomic features
#'
#' @author
#' Mengjun Wu, Lei Gu
#'
#' @seealso
#'
#' \code{\link{glmLRT}}
#'
#' @examples
#' data(tca_ATAC)
#' tca_ATAC <- DBanalysis(tca_ATAC)
#' ### extract differntial analysis of 24h, 72h to 0h
#' # set the contrasts using the 'group1' and 'group2' paramters
#' res1 <- DBresult(tca_ATAC, group1 = '0h', group2 = c('24h', '72h'))
#' # one can get the same result by setting the contrasts using hte 'contrasts' parameter
#' res2 <- DBresult(tca_ATAC, contrasts = c('24hvs0h', '72hvs0h'))
#' # extract significant diffential events
#' res.sig <- DBresult(tca_ATAC, contrasts = c('24hvs0h', '72hvs0h'),
#'                    top.sig = TRUE)
#'
#' # extract differntial analysis of 24h, 72h to 0h of a given cluster
#' tca_ATAC <- timecourseTable(tca_ATAC, filter = TRUE)
#' tca_ATAC <- timeclust(tca_ATAC, algo = 'cm', k = 6)
#' res_cluster1 <- DBresult.cluster(tca_ATAC, group1 = '0h',
#'                                  group2 = c('24h', '72h'),
#'                                  cluster = 1)
#'
#'
#'
#' @export
DBresult <- function(object, group1 = NULL, group2 = NULL,
                     contrasts = NULL, p.adjust = "fdr",
                     top.sig = FALSE, pvalue = "paj",
                     pvalue.threshold = 0.05, abs.fold = 2,
                     direction = "both", result.type = "GRangesList") {
  if (is.null(group1) && is.null(group2) && is.null(contrasts)) {
    stop("Either information of groups to compare or \"contrasts\" should be provided")
  }
  if (!is.null(contrasts)){
    contrasts <- contrasts
  }else{
    if (sum(group1 %in% group2) > 0) {
      warning("Members in group1 are also found in group2, overlapped members are removed from group2.")
      group2 <- group2[-which(group2 %in% group1)]
    }
    contrasts <- contrastNames(group1, group2)
  }
  if (!p.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
    stop("Method for adjusting P-values should be one of following methods: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'. Character string is case sensitive.")
  }
  fit <- object@DBfit
  contrast.table <- object@contrasts
  gi <- object@genomicFeature[object@genomicFeature$id %in%
                                row.names(fit$coefficients), ]
  gi <- gi[, c("chr", "start", "end", "id")]
  res <- list()
  for (i in contrasts) {
    tmp <- glmLRT(fit, contrast = contrast.table[, i])
    restmp <- tmp$table[, c(1, 4)]
    adjustp <- p.adjust(restmp$PValue, method = p.adjust)
    restmp <- cbind(restmp, adjustp)
    colnames(restmp)[length(restmp[1, ])] <- "paj"
    restmp <- cbind(restmp, gi, stringsAsFactors = FALSE)
    res[[i]] <- restmp
  }
  if (top.sig) {
    res <- DBresult.filter(x = res, pvalue = pvalue,
                           pvalue.threshold = pvalue.threshold,
                           abs.fold = abs.fold,
                           direction = direction)
  }
  if (tolower(result.type) == "grangeslist") {
    gr <- as(do.call(rbind, unname(res)), "GRanges")
    res <- suppressWarnings(split(gr, rep(names(res), lengths(res))))
  }

  res
}

#' @rdname DBresult
#' @export
DBresult.cluster <- function(object, group1 = NULL, group2 = NULL,
                             contrasts = NULL, p.adjust = "fdr",
                             top.sig = FALSE, pvalue = "paj",
                             pvalue.threshold = 0.05, abs.fold = 2,
                             direction = "both",cluster, cmthreshold = NULL,
                             result.type = "GRangesList") {
  if (length(object@clusterRes@cluster) == 0) {
    stop("No cluster information provided, clustering analysis must be performed first")
  }
  DBres <- DBresult(object, group1 = group1, group2 = group2,
                    contrasts = contrasts, p.adjust = p.adjust,
                    top.sig = top.sig, pvalue = pvalue,
                    pvalue.threshold = pvalue.threshold,
                    abs.fold = abs.fold,
                    direction = direction, result.type = "list")
  names <- names(object@clusterRes@cluster)
  res <- list()
  contrast_name <- names(DBres)
  counter <- 0
  for (i in DBres) {
    restmp <- i
    counter <- counter + 1
    clusters <- object@clusterRes@cluster
    clusternames <- names[which(clusters == cluster)]
    if (!is.null(cmthreshold)) {
      membership <- object@clusterRes@membership[clusters == cluster,
                                                 cluster]
      if (is.null(membership)) {
        stop("No membership matrix found. To get membership matrix, please choose 'cmeans' clustering method when calling timeclust")
      } else {
        clusternames <- clusternames[which(membership > cmthreshold)]
      }
    }
    restmp <- restmp[clusternames, ]
    contrast <- contrast_name[counter]
    res[[contrast]] <- restmp
  }
  if (tolower(result.type) == "grangeslist") {
    gr <- as(do.call(rbind, unname(res)), "GRanges")
    res <- suppressWarnings(split(gr, rep(names(res), lengths(res))))
  }
  res
}

# contrast contrast by given strings, group1 is a string, group2 can be a string or a vector of strings
contrastNames <- function(group1, group2) {
  b <- length(group2)
  name <- vector(mode = "character", length = b)
  for (i in seq_len(b)) {
    n <- paste0(group2[i], "vs", group1)
    name[i] <- n
  }
  name
}

DBresult.filter <- function(x, pvalue = "paj", pvalue.threshold = 0.05,
                            abs.fold = 2, direction = "both") {
  if (abs.fold < 0) {
    err <- paste0("\"abs.fold\" should be postive number.")
    stop(err)
  }
  d <- x
  for (i in seq_len(length(d))) {
    dt <- d[[i]]
    if (direction == "up") {
      dt <- dt[which(dt$logFC >= abs.fold), ]
    }

    if (direction == "down") {
      dt <- dt[which(dt$logFC <= -abs.fold), ]
    }
    if (direction == "both") {
      if (abs.fold == 0) {
        dt <- rbind(dt[which(dt$logFC >= abs.fold), ],
                    dt[which(dt$logFC < -abs.fold), ])
      } else {
        dt <- rbind(dt[which(dt$logFC >= abs.fold), ],
                    dt[which(dt$logFC <= -abs.fold), ])
      }
    }
    dt <- dt[which(dt[, pvalue] < pvalue.threshold), ]
    d[[i]] <- dt
  }
d
}
