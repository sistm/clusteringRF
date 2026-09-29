#' S3 print method for clusteringRF objects
#'
#' @param x a \code{clusteringRF}.
#' @param ... other arguments to be passed to and from other methods.
#' @export

print.clusteringRF <- function(x, ...){

  print(paste0("A clustering random forest with ", length(x), " trees"))

}

