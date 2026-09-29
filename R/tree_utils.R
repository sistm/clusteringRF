#' @title Constructs the path rules for each cluster
#' @description Function which builds a data frame describing, for every cluster, 
#' the associated variable intervals.
#' @param description monothetic description of the leaves
#' @return a data frame with six columns : the cluster label, the variable name, the lower 
#' and upper bounds of the variable's range (for quantitative variables), 
#' the categories (for qualitative variables), and the variable type
#' @importFrom stringr str_split str_split_fixed str_trim str_remove_all str_extract_all
#' @importFrom utils tail
#' @export

make_path_cluster <- function(description){
  path_cluster <- data.frame(cluster = character(), 
                             variable = character(), 
                             min = numeric(), 
                             max = numeric(), 
                             categories = character(), 
                             type = character())
  
  for (cluster in names(description)){
    conds <- str_split(description[[cluster]], " , ")[[1]] |> str_trim()
    parts <- str_split_fixed(conds, "=", 2)
    vars <- str_remove_all(parts[,1], " ")
    exprs <- str_trim(parts[,2])
    for (v in unique(vars)){
      idx <- tail(which(vars == v),1)
      expr <- exprs[idx]
      if(str_detect(expr, "^\\{")){
        cats <- expr |> str_remove_all("[{}]") |> str_split(",") |> unlist() |> str_trim()
        path_cluster <- rbind(path_cluster, 
                              data.frame(cluster = cluster,
                                         variable = v,
                                         min = NA_real_, 
                                         max = NA_real_, categories = I(list(cats)), 
                                         type = "quali"))
      }else{
        limits <- as.numeric(str_extract_all(exprs[idx],"-?\\d+\\.?\\d*|-Inf|Inf")[[1]])
        path_cluster <- rbind(path_cluster, 
                              data.frame(cluster = cluster,
                                         variable = v,
                                         min = limits[1], 
                                         max = limits[2], categories = I(list(NA_character_)), 
                                         type = "quanti")) 
      }
    }
  }
  return(path_cluster)}


#' @title Predicts cluster membership
#' @description Function which assigns each observation to the right cluster.
#' @param X a data frame of observations
#' @param path_cluster a data frame produced by \code{make_path_cluster()}
#' @return a data frame with the predicted cluster label for each observation
#' @export
make_prediction <- function(X, path_cluster) {
  n        <- nrow(X)
  clusters <- unique(path_cluster$cluster)
  nK       <- length(clusters)
  
  rules_list <- split(path_cluster, path_cluster$cluster)[clusters]
  
  membership <- matrix(TRUE, n, nK, dimnames = list(NULL, clusters))
  
  for (k in seq_len(nK)) {
    rules <- rules_list[[k]]
    for (r in seq_len(nrow(rules))) {
      v <- rules$variable[r]
      t <- rules$type[r]
      if (t == "quali"){
        membership[, k] <- membership[, k] &
          (X[[v]] %in% rules$categories[r][[1]])
      }else{
        membership[, k] <- membership[, k] &
          (X[[v]] >= rules$min[r]) & (X[[v]] < rules$max[r])
      }
    }
  }
  
  pred_idx   <- max.col(membership, ties.method = "first")
  data.frame(prediction = clusters[pred_idx], stringsAsFactors = FALSE)
}


#' @title Computes Clustering Stability (CS) Metrics
#' @description Function which estimates the importance of each variable on the clustering
#' using a permuation approach.
#' @param X_oob a data frame of out-of-bag observations
#' @param path_cluster a data frame produced by \code{make_path_cluster()}
#' @param prediction a data frame produced by \code{make_prediction()}
#' @return a list with 2 elements : Out-of-bag Reassignment Rate (CS_ORR) which counts 
#' the proportion of observations whose cluster assignment changed after
#' permutation and ARI Stability (CS_ARI) which computes adjusted rand index between the
#' original and permuted cluster assignments. 
#' @importFrom mclust adjustedRandIndex
#' @export
make_CS_metrics <- function(X_oob, path_cluster, prediction){
  p_oob <- ncol(X_oob)
  CS_ORR <- numeric(p_oob)
  CS_ARI <- numeric(p_oob)
  names(CS_ORR) <- names(CS_ARI) <- names(X_oob)
  for (j in seq_len(p_oob)){
    X_perm <- X_oob
    X_perm[[j]] <- sample(X_perm[[j]])
    prediction_j <- make_prediction(X_perm, path_cluster)
    CS_ORR[j] <- sum(prediction$prediction != 
                        prediction_j$prediction)/nrow(X_oob)
    CS_ARI[j] <- 1 - adjustedRandIndex(prediction$prediction, 
                                    prediction_j$prediction)
  }
  return(list(CS_ORR = CS_ORR, CS_ARI = CS_ARI))
}


#' @title Computes Clustering Stability (CS) Metrics for inertia distance
#' @description Function which estimates the importance of each variable on the clustering
#' using an approach based on the inertia path.
#' @param X_oob a data frame of out-of-bag observations
#' @param path_cluster a data frame produced by \code{make_path_cluster()}
#' @param prediction a data frame produced by \code{make_prediction()}
#' @param dist_clusters a distance matrix 
#' @return a list with 1 element : Inertia Path (CS_IP) which computes the sum  
#' of inertia between the original and permuted clusters assignments. 
#' @importFrom mclust adjustedRandIndex
#' @export
make_CS_metrics_inertia <- function(X_oob, path_cluster, prediction, dist_clusters){
  p_oob <- ncol(X_oob)
  CS_IP <- numeric(p_oob)
  names(CS_IP) <- names(X_oob)
  for (j in seq_len(p_oob)){
    X_perm <- X_oob
    X_perm[[j]] <- sample(X_perm[[j]])
    prediction_j <- make_prediction(X_perm, path_cluster)
    for (i in 1:nrow(prediction)){
      CS_IP[j] <- CS_IP[j] + dist_clusters[prediction[[1]][[i]], prediction_j[[1]][[i]]]
    }
  }
  return(list(CS_IP = CS_IP))

}
