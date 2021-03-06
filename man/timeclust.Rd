% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/timeclust.R
\name{timeclust}
\alias{timeclust}
\title{time couse data clustering}
\usage{
timeclust(x, algo, k, dist = "euclidean", centers = NULL,
  standardize = TRUE, ...)
}
\arguments{
\item{x}{a \code{TCA} object returned from
\code{\link{timecourseTable}} or a matrix}

\item{algo}{character string giving a clustering method. Options
are \code{km}' (kmeans), '\code{pam}' (partitioning around medoids),
'\code{hc}' (hierachical clustering), '\code{cm}' (cmeans).}

\item{k}{numeric value between \eqn{1} and \eqn{n - 1} ( \eqn{n}
is the number of data points to be clustered ).}

\item{dist}{character string specifying method for
distance(dissimilarity) calculation. It should be one of
'\code{correlation}' or one of the distance measure method in
\code{\link{dist}} function (for example '\code{euclidean}',
'\code{manhattan}')}

\item{centers}{a numeric matrix giving intial centers for kmeams,
pam or cmeans. If given, Number of rows of centers must be equal
to k.}

\item{standardize}{logical, if TRUE, z-score transformation will
performed on the data before clustering. See 'Details' below.}

\item{...}{additional arguments passing to \code{\link{kmeans}},
\code{\link{pam}}, \code{\link{hclust}}, \code{\link{cmeans}}}
}
\value{
If x is a \code{TCA} object, a \code{TCA} object will be returned.
If x is a matrix, a \code{clust} object will be returned
}
\description{
This function performs clustering analysis of time course data.
}
\details{
two types of clustering methods are provided: hard clustering
(\code{\link{kmeans}}, \code{\link{pam}}, \code{\link{hclust}})
and soft clustering(\code{\link{cmeans}}). In Hard clustering,
a data point can only be allocated to exactly one cluster
(for \code{\link{hclust}}, \code{\link{cutree}} is used to cut
a tree into clusters), while in soft clustering (also known as
fuzzy clustering), a data point can be assigned to multiple
clusters, membership values are used to indicate to what
degree a data point belongs to each cluster. For more details,
see the help() page of each function.

To avoid the influence of expression level to the clustering
analysis, z-score transformation can be applied to covert the
expression values to z-scores by performing the following formula:

\deqn{z = \frac{x - \mu}{\sigma}}

\eqn{x} is value to be converted (e.g., a expression value of a
genomic feature in one condition), \eqn{\mu} is the population
mean (e.g., average expression value of a genomic feature in
different conditions), \eqn{\sigma} is the standard deviation
(e.g., standard deviation of expression of a genomic feature
in different conditions).
}
\examples{

x <- matrix(rnorm(1600,sd=0.3), nrow = 200,
            dimnames = list(paste0('peak', 1:200), 1:8))
clust_res <- timeclust(x, algo = 'cm', k = 4) # return a clust object

}
\author{
Mengjun Wu
}
\seealso{
\code{\link{clust}}, \code{\link{kmeans}},
\code{\link{pam}}, \code{\link{hclust}}, \code{\link{cutree}}
}

