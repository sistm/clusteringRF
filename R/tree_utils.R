#' @title Constructs the path rules for each cluster
#' @description Function which builds a data frame describing, for every cluster, 
#' the associated variable intervals.
#' @param description monothetic description of the leaves
#' @return a data frame with 4 columns : cluster label, variable name, lower 
#' bound of the variable range and upper bound of the variable range
#' @importFrom stringr str_split str_split_fixed str_trim str_remove_all str_extract_all
#' @importFrom utils tail
#' @export

make_path_cluster <- function(description){
  path_cluster <- data.frame(cluster = character(), 
                             variable = character(), 
                             min = numeric(), 
                             max = numeric())
  
  for (cluster in names(description)){
    elems <- str_split(description[[cluster]], ",")[[1]] |> str_trim()
    parsed <- str_split_fixed(elems, "=", 2)
    vars <- str_remove_all(parsed[,1], " ")
    exprs <- str_trim(parsed[,2])
    for (v in unique(vars)){
      idx <- tail(which(vars == v),1)
      limits <- as.numeric(str_extract_all(exprs[idx],"-?\\d+\\.?\\d*|-Inf|Inf")[[1]])
      path_cluster <- rbind(path_cluster, 
                            data.frame(cluster = cluster,
                                       variable = v,
                                       min = limits[1], 
                                       max = limits[2]))
    }
  }
  return(path_cluster)}


#' @title Predicts cluster membership
#' @description Function which assigns each observation to the right cluster.
#' @param X a data frame of observations
#' @param path_cluster a data frame produced by \code{make_path_cluster()}
#' @return a data frame with the predicted cluster label for each observation
#' @export
make_prediction <- function(X, path_cluster){
  n <- nrow(X)
  prediction <- character(n)
  clusters <- unique(path_cluster$cluster)
  for (i in seq_len(n)){
    X_i <- X[i, , drop = FALSE]
    scores <- numeric(length(clusters))
    names(scores) <- clusters
    for (cluster in clusters){
      rules <- path_cluster[path_cluster$cluster == cluster, ]
      score <- all(sapply(rules$variable, function(v){
        X_i[[v]] >= rules[rules$variable == v, "min"] && 
          X_i[[v]] < rules[rules$variable == v, "max"]
      }))
      scores[cluster] <- score
    }
    prediction[i] <- names(which.max(scores))
  }
  return(data.frame(prediction = prediction, stringsAsFactors = FALSE))
}


#' @title Computes Clustering Stability Importance (CSI) Metrics for co-clustering distance
#' @description Function which estimates the importance of each variable on the clustering
#' using a permuation approach.
#' @param X_oob a data frame of out-of-bag observations
#' @param path_cluster a data frame produced by \code{make_path_cluster()}
#' @param prediction a data frame produced by \code{make_prediction()}
#' @return a list with 2 elements : Out-of-bag Reassignment Rate (CSI_ORR) which counts 
#' the proportion of observations whose cluster assignment changed after
#' permutation and ARI Stability Importance (CSI_ARI) which computes adjusted rand index between the
#' original and permuted cluster assignments. 
#' @importFrom mclust adjustedRandIndex
#' @export
make_CSI_metrics_coclust <- function(X_oob, path_cluster, prediction){
  p_oob <- ncol(X_oob)
  CSI_ORR <- numeric(p_oob)
  CSI_ARI <- numeric(p_oob)
  names(CSI_ORR) <- names(CSI_ARI) <- names(X_oob)
  for (j in seq_len(p_oob)){
    X_perm <- X_oob
    X_perm[[j]] <- sample(X_perm[[j]])
    prediction_j <- make_prediction(X_perm, path_cluster)
    CSI_ORR[j] <- sum(prediction$prediction != 
                        prediction_j$prediction)/nrow(X_oob)
    CSI_ARI[j] <- adjustedRandIndex(prediction$prediction, 
                                    prediction_j$prediction)
  }
  return(list(CSI_ORR = CSI_ORR, CSI_ARI = CSI_ARI))
}


#' @title Computes Clustering Stability Importance (CSI) Metrics for inertia distance
#' @description Function which estimates the importance of each variable on the clustering
#' using a permutation approach.
#' @param X_oob a data frame of out-of-bag observations
#' @param path_cluster a data frame produced by \code{make_path_cluster()}
#' @param prediction a data frame produced by \code{make_prediction()}
#' @param dist_clusters a distance matrix 
#' @return a list with 1 element : Inertia Permutation (CSI_IP) which computes the sum  
#' of inertia between the original and permuted clusters assignments. 
#' @importFrom mclust adjustedRandIndex
#' @export
make_CSI_metrics_inertia <- function(X_oob, path_cluster, prediction, dist_clusters){
  p_oob <- ncol(X_oob)
  CSI_IP <- numeric(p_oob)
  names(CSI_IP) <- names(X_oob)
  for (j in seq_len(p_oob)){
    X_perm <- X_oob
    X_perm[[j]] <- sample(X_perm[[j]])
    prediction_j <- make_prediction(X_perm, path_cluster)
    for (i in 1:nrow(prediction)){
      CSI_IP[j] <- CSI_IP[j] + dist_clusters[prediction[[1]][[i]], prediction_j[[1]][[i]]]
    }
  }
  return(list(CSI_IP = CSI_IP))
}