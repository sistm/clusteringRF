#' Unsupervised Random Forest for clustering
#'
#' main wrapper
#'
#' @param X a data.frame to be clustered
#' @param ntrees number of \code{divlust} trees. Default is \code{500}.
#' @param mtry number of variables selected at each \code{divlust} tree node. 
#' Default is \code{ncol(X)}.
#' @param distance a character string, either "co-clustering" or "inertia".
#' @param pruning Logical (TRUE or FALSE). If TRUE, inertia trees are pruned. 
#' @param weighting Logical (TRUE or FALSE). If TRUE, inertia nodes are weighted 
#' in the calculation of variable importance.
#' @param ncores number of cpus to parallelize over. Default is
#' \code{parallel::detectCores()-1}.
#'
#' @return a list of length \code{n_trees}, each elements containing 3 matrices.
#'
#' @importFrom pbapply pblapply
#' @importFrom parallel detectCores
#' @importFrom tibble is_tibble
#'
#' @export
#'
#' @examples
#' library(palmerpenguins)
#' mypeng <- as.data.frame(penguins)
#' mypeng$year <- factor(as.character(mypeng$year),
#'                          levels=c("2007", "2008", "2009"),
#'                          ordered=TRUE)
#' set.seed(123)
#' forest_clust <- clusteringRF(na.omit(mypeng[mypeng$sex=="male", -c(1, 7)]), ntrees = 50, ncores = 1)
#'
#' if(interactive()){
#'  resume <- summary(forest_clust)
#'
#'  plot(resume$Matrix)
#' }
#'
#' if(interactive()){
#'    ncol <- 15
#'    nobs_pg <- 43
#'    diff_mean <- 5
#'    X <- matrix(rnorm(n = ncol*nobs_pg, mean = 0, sd = 1), ncol=ncol)
#'    Y <- matrix(rnorm(n = ncol*nobs_pg, mean = diff_mean, sd = 1), ncol=ncol)
#'    mysim_data <- rbind.data.frame(X,Y)
#'    res <- clusteringRF(mysim_data)
#'    s <- summary(res)
#'    plot(s)
#' }

clusteringRF <- function(X, ntrees = 500, mtry = ncol(X), distance = "co-clustering", pruning = TRUE, weighting = FALSE, ncores = parallel::detectCores()-1){

  stopifnot(is.data.frame(X))
  stopifnot(distance %in% c("co-clustering", "inertia"))

  if(tibble::is_tibble(X)){
    X <- as.data.frame(X)
    warning("X is a tibble. It is converted to a simple `data.frame` instead.")
  }

  if(is.null(rownames(X))){
    rownames(X) <- as.character(1:nrow(X))
    message("No rownames for `X`. using row number instead.")
  }

  forest <- pblapply(1:ntrees, function(i){
    tree(X, mtry, distance, pruning, weighting)
  }, cl = ncores)

  class(forest) <- "clusteringRF"
  #message("We advise you do use the `summary()` on this object to agregate the result of this forest, before plotting the summary itself.")

  return(forest)
}
