#' Plot clustering results for time course data.
#'
#' This function plots the clusters generated from
#' \code{\link{timeclust}}. For fuzzy cmeans clustering, data points
#' are color-coded according to membership values, the color palettes
#' can be customized.
#'
#' @param object a \code{TCA} object or a \code{clust} object
#'
#' @param categories character string giving the x-axis label
#'
#' @param value character string giving the y-axis label
#'
#' @param cols integer value specifying number of columns in the final
#' layout.
#'
#' @param cl.color  character string specifying a color for hard
#' clustering.
#'
#' @param membership.color  color palettes, a character vector of
#' n colors
#'
#' @param title.size numeric value specifying the font size of title
#' of each
#' plot in the layout
#'
#' @param axis.line.size numeric value specifying the size of both
#' axis lines
#'
#' @param axis.title.size numeric value specifying the font size of
#' titles of both axis
#'
#' @param axis.text.size numeric value specifying the font size of
#' labels of both axis
#'
#' @param legend.title.size numeric value specifying the font size
#' of legend title
#'
#' @param legend.text.size numeric value specifying the font size of
#' legend text
#'
#' @return
#' Plot all clusters in one plot and return a list of ggplot objects,
#' each object is for one cluster. The ggplot object can be drawed by
#' calling \code{\link{print.ggplot}}
#'
#' @examples
#' x <- matrix(sample(500, 1600, replace = TRUE), nrow = 200,
#'             dimnames = list(paste0('peak', 1:200), 1:8))
#' clust_res <- timeclust(x, algo = 'cm', k = 4, standardize = TRUE)
#' p <- timeclustplot(clust_res, cols =2)
#' # to plot a individual cluster
#' print (p[[2]]) # plot cluster 2
#' print (p[[3]]) # plot cluster 3
#'
#' @author
#' Mengjun Wu
#' @export

timeclustplot <- function(object = NULL, categories = "timepoint",
                          value = "expression", cols = NULL,
                          cl.color = "gray50",
                          membership.color = rainbow(30, s = 3/4, v = 1, start = 1/6),
                          title.size = 18, axis.line.size = 0.6,
                          axis.title.size = 18,
                          axis.text.size = 16, legend.title.size = 14,
                          legend.text.size = 14) {

  if (class(object) != "clust" && class(object) != "TCA") {
    stop("object should be a 'timeclust' object or a 'TCA' object")
  }
  if (class(object) == "clust") {
    data <- object@data
    cluster <- object@cluster
    membership <- object@membership
  }
  if (class(object) == "TCA") {
    data <- object@clusterRes@data
    cluster <- object@clusterRes@cluster
    membership <- object@clusterRes@membership
  }
  ncl <- max(cluster)
  membercolor <- vector(length = length(cluster))
  membervalue <- list()
  counter <- 0
  if (!sum(dim(membership) == 0) == 2) {
    color <- membership.color
    colorseq <- seq(0, 1, length = length(color))
    for (i in seq_len(ncl)) {
      mtmp <- membership[cluster == i, i]
      membervalue[[i]] <- mtmp
      for (j in seq_len(length(mtmp))) {
        counter <- counter + 1
        ind <- which(abs(colorseq - mtmp[j]) == min(abs(colorseq - mtmp[j])))
        membercolor[counter] <- color[ind]
      }
    }
    membervalue <- unlist(membervalue)
    names(membercolor) <- membervalue
  }

  plotlist <- list()
  for (i in seq_len(ncl)) {
    title <- paste0("Cluster ", i)
    dtmp <- data[cluster == i, ]
    a <- which(cluster == i)
    if (length(a) == 1) {
      dtmp <- data.frame(time = 1:length(dtmp), value = dtmp)
      if (!sum(dim(membership) == 0) == 2) {
        m <- membership[cluster == i, i]
        colorname = toString(m)
        plotlist[[i]] <- ggplot(dtmp, aes(x = time, y = value)) +
          geom_line(colour = membercolor[colorname]) + theme_bw() +
          ggtitle(title) +
          scale_x_continuous(breaks = dtmp$time,
                             labels = row.names(dtmp)) +
          labs(x = categories, y = value) +
          theme(plot.title = element_text(size = title.size),
                axis.line.x = element_line(color = "black",
                                           size = axis.line.size),
                axis.line.y = element_line(color = "black",
                                           size = axis.line.size),
                axis.title = element_text(size = axis.title.size),
                axis.text = element_text(size = axis.text.size),
                legend.position = "none", panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
      } else {
        plotlist[[i]] <- ggplot(dtmp, aes(x = time, y = value)) +
          geom_line(colour = cl.color) + theme_bw() + ggtitle(title) +
          scale_x_continuous(breaks = dtmp$time,
                             labels = row.names(dtmp)) +
          labs(x = categories, y = value) +
          theme(plot.title = element_text(size = title.size),
                axis.line.x = element_line(color = "black",
                                           size = axis.line.size),
                axis.line.y = element_line(color = "black",
                                           size = axis.line.size),
                axis.title = element_text(size = axis.title.size),
                axis.text = element_text(size = axis.text.size),
                legend.position = "none", panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
      }
    } else {
      dtmp_m <- melt(dtmp)
      colnames(dtmp_m) <- c("group", "time", "value")
      if (sum(dim(membership) == 0) == 2) {
        plotlist[[i]] <- ggplot(dtmp_m, aes(x = time, y = value)) +
          geom_line(aes(group = group), colour = cl.color) +
          theme_bw() + ggtitle(title) +
          labs(x = categories, y = value) +
          theme(plot.title = element_text(size = title.size),
                axis.line.x = element_line(color = "black",
                                           size = axis.line.size),
                axis.line.y = element_line(color = "black",
                                           size = axis.line.size),
                axis.title = element_text(size = axis.title.size),
                axis.text = element_text(size = axis.text.size),
                legend.position = "none", panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
      }
      if (!sum(dim(membership) == 0) == 2) {
        mem <- membership[cluster == i, i]
        mem1 <- data.frame(group = names(mem), member = mem)
        dtmp_m1 <- merge(dtmp_m, mem1, by = "group")
        colnames(dtmp_m1) <- c("group", "time", "value", "membership")
        dtmp_m1 <- dtmp_m1[order(dtmp_m1[, 4]), ]
        new.factor <- unique(as.vector(dtmp_m1$group))
        dtmp_m1$group <- factor(dtmp_m1$group, levels = new.factor)

        plotlist[[i]] <- ggplot(dtmp_m1, aes(x = time, y = value,
                                             colour = membership)) +
          geom_line(aes(group = group)) +
          scale_colour_gradientn(colours = membership.color) +
          guides(colour = guide_colourbar()) + theme_bw() +
          ggtitle(title) + labs(x = categories, y = value) +
          theme(plot.title = element_text(size = title.size),
                axis.line.x = element_line(color = "black",
                                           size = axis.line.size),
                axis.line.y = element_line(color = "black",
                                           size = axis.line.size),
                axis.title = element_text(size = axis.title.size),
                axis.text = element_text(size = axis.text.size),
                legend.title = element_text(size = legend.title.size),
                legend.text = element_text(size = legend.title.size),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())


      }
    }

  }
  suppressWarnings(multiplot(plotlist = plotlist, cols = cols))
  plotlist
}

multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])

  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout),
                                               ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
