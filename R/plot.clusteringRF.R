#' S3 plot method for clusteringRF.summary objects
#'
#' @param x the cumulative similarity matrix
#' @param ... other arguments to be passed to and from other methods.

#' @return Provides the associated heatmap of the similarity matrix
#' @import dplyr gplots ggplot2 GGally
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.hclust complete.cases cutree
#' @export

plot.clusteringRF.summary <- function(x, ...){

  my_palette <- colorRampPalette(c("lightyellow1", "yellow",'orange2','red3', "red4"))(n = 299)

  map <- heatmap.2(1-x$distance_matrix,
                   trace = "none",
                   symkey = F,
                   symm = T,
                   dendrogram = "column",
                   scale = "none",
                   col = my_palette,
                   ...)

  return(invisible(map))
}

