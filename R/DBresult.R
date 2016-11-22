#' Perform likelihood ratio tests and extract the differential analysis results
#'
#' This function performs likelihood ratio tests for given coefficinets contrasts
#' after fitting read counts to GLM by \code{\link{DBanalysis}}. \code{DBresult} extracts
#' the result tables of differential analysis of given contrasts and the result tables contain
#' all genomic features. \code{DBresult.cluster} returns similar result tables while the result
#' tables only contains genomic features belong to a given cluster.
#'
#' @name DBresult
#'
#' @param object a \code{TCA} object, for \code{DBresult}, \code{DBanalysis} should already be called
#' on the object; for \code{DBresult.cluster}, both \code{DBanalysis} and \code{timeclust} should be
#' already called
#'
#' @param group1 character string giving the level to be compared, that is the denominator in the fold changes.
#'
#' @param group2 a character vetor giving other levels to compared with \code{group1}. that are numerator
#' in the fold changes.
#'
#' @param contrasts a character vector, each charcter string in the vector gives a contrast of two groups with
#' the format 'group2vsgroup1', group1 is the denominator level in the fold changes and group2 is the numerator
#' level in the fold changes.
#'
#' @param p.adjust character string specifying a correction method for p-values. Options are 'holm', 'hochberg',
#' 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'.
#'
#' @param cluster an integer, the result tables of genomic features belong to
#' the \code{cluster} are extracted.
#'
#' @param  cmthreshold a numeric value, this argument is applicable only if '\code{cmeans}' clustering method
#' is selected when calling \code{\link{timeclust}} function. if not NULL, the result table of genomic features
#' that belong to the defined \code{cluster} and the membership values to this cluster exceed \code{cmthreshold}
#' are extracted.
#'
#' @details This function uses \code{\link{glmLRT}} from edgeR which perform likelihood ratio tests for testing
#' significance of changes. For more deatils, see \code{\link{glmLRT}}
#'
#' @note If not NULL \code{group1}, \code{group2} and \code{contrasts}, result tables are extracted from comparisons
#' in \code{constrasts}.
#'
#' @return
#' A named list of result tables for all provided contrasts, the names of the list are the contrasts. Each result table
#' is a data frame containing the following columns:
#'
#' @return \code{logFC} log2-fold change of differential event between two tested.
#'
#' @return \code{PValue} p-values.
#'
#' @return \code{paj} adjusted p-values
#'
#' @return \code{chr}  name of the chromosomes
#'
#' @return \code{start} starting position of the feature in the chromosome
#'
#' @return \code{end} ending postition of the feature in the chromosome
#'
#' @return \code{id} genomic feature name
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
#'
#' # extract differntial analysis of 24h, 72h to 0h of a given cluster
#' tca_ATAC <- timecourseTable(tca_ATAC, filter = TRUE)
#' tca_ATAC <- timeclust(tca_ATAC, algo = 'cm', k = 6)
# res_cluster1 <- DBresult.cluster(tca_ATAC, group1 = '0h', group2 = c('24h', '72h'), cluster = 1)
#'
#' @export
DBresult <- function(object, group1 = NULL, group2 = NULL, contrasts = NULL, p.adjust = "fdr") {

    res <- list()
    names <- NULL
    if (is.null(group1) && is.null(group2) && is.null(contrasts)) {
        stop("Either information of groups to compare or \"contrasts\" should be provided")
    }
    if (is.null(contrasts)) {
        if (is.null(group1) || is.null(group2)) {
            stop("One group is missing, two groups should be provided for differential event analysis.")
        }
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
    gi <- object@genomicFeature[object@genomicFeature$id %in% row.names(fit$coefficients), ]
    gi <- gi[, c("chr", "start", "end", "id")]
    for (i in contrasts) {
        tmp <- glmLRT(fit, contrast = contrast.table[, i])
        restmp <- tmp$table[, c(1, 4)]
        if (!p.adjust == "none") {
            adjustp <- p.adjust(restmp$PValue, method = p.adjust)
            restmp <- cbind(restmp, adjustp)
            colnames(restmp)[length(restmp[1, ])] <- "paj"
        }
        restmp <- cbind(restmp, gi, stringsAsFactors = FALSE)
        res[[i]] <- restmp
    }
    res
}

#' @rdname DBresult
#' @export
DBresult.cluster <- function(object, group1 = NULL, group2 = NULL, contrasts = NULL, p.adjust = "fdr", cluster,
    cmthreshold = NULL) {
    if (length(object@clusterRes@cluster) == 0) {
        stop("No cluster information provided, clustering analysis must be performed first")
    }
    DBres <- DBresult(object, group1 = group1, group2 = group2, contrasts = contrasts, p.adjust = p.adjust)
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
            membership <- object@clusterRes@membership[clusters == cluster, cluster]
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
    res
}

# contrast contrast by given strings, group1 is a string, group2 can be a string or a vector of strings
contrastNames <- function(group1, group2) {
    name <- c()
    for (i in group2) {
        n <- paste0(i, "vs", group1)
        name <- c(name, n)
    }
    name
}
