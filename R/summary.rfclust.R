#' summary S3 method for rfclust objects
#'
#' Merge all matrices and some analysis
#'
#' @param object the results of the apply function on the RF function (tree)
#' @param ... more parameters if required
#' @return Provides the cumulative similarity matrix and some analysis
#' @import dplyr gplots ggplot2 GGally
#' @export


summary.rfclust <- function(object, ...){

  ntrees <- length(object)

  matrices_sim <- lapply(object,'[[',1)
  sum_sim <- Reduce('+', matrices_sim)

  matrices_dist <- lapply(object,'[[',2)
  sum_dist <- Reduce('+', matrices_dist)

  matrices_absent <- lapply(object,'[[',3)
  sum_absent <- Reduce('+', matrices_absent)
  
  MDI_importance <- lapply(object, '[[', 5)
  
  MDA_importance_permutation_count <- lapply(object, '[[', 6)
  
  MDA_importance_ARI_permutation <- lapply(object, '[[', 7)

  if(object[[1]]$distance == "co-clustering"){

    pair_appearances <- sum_dist + sum_sim                                          #Somme des occurence des paires dans les forêts.
    similarity_matrix <- sum_sim / pair_appearances
    distance_matrix <- 1 - similarity_matrix

  }else if(object[[1]]$distance == "inertia"){
    sum_present <- ntrees - sum_absent
    distance_matrix <- sum_dist / sum_present
  }

  # Diagonale nulle
  diag(distance_matrix) <- 0

  # Si NA, mettre à 0
  distance_matrix[is.na(distance_matrix)] <- 0
  
  #Importance des variables MDI 
  MDI_importance_var <- list()
  for (j in 1:length(MDI_importance)){
    for (var in names(MDI_importance[[j]])){
      if (is.null(MDI_importance_var[[var]])){
        MDI_importance_var[[var]] <- 0
      }
      MDI_importance_var[[var]]  <- MDI_importance_var[[var]] + MDI_importance[[j]][[var]]
    }
  }
  
  MDI_importance_moy_var <- lapply(MDI_importance_var, function(x) x/ntrees)
  
  #Importance des variables MDA 
  MDA_importance_permutation_count_var <- list()
  for (j in 1:length(MDA_importance_permutation_count)){
    for (var in names(MDA_importance_permutation_count[[j]])){
      if (is.null(MDA_importance_permutation_count_var[[var]])){
        MDA_importance_permutation_count_var[[var]] <- 0
      }
      MDA_importance_permutation_count_var[[var]]  <- MDA_importance_permutation_count_var[[var]] + MDA_importance_permutation_count[[j]][[var]]
    }
  }
  
  MDA_importance_permutation_count_moy_var <- lapply(MDA_importance_permutation_count_var, function(x) x/ntrees)
  
  #Importance des variables MDA 
  MDA_importance_ARI_permutation_var <- list()
  for (j in 1:length(MDA_importance_ARI_permutation)){
    for (var in names(MDA_importance_ARI_permutation[[j]])){
      if (is.null(MDA_importance_ARI_permutation_var[[var]])){
        MDA_importance_ARI_permutation_count_var[[var]] <- 0
      }
      MDA_importance_ARI_permutation_var[[var]]  <- MDA_importance_ARI_permutation_var[[var]] + MDA_importance_ARI_permutation[[j]][[var]]
    }
  }
  
  MDA_importance_ARI_permutation_moy_var <- lapply(MDA_importance_ARI_permutation_var, function(x) x/ntrees)
  
  output_summary <- list("distance_matrix" = distance_matrix, "MDI_importance_variables" = MDI_importance_moy_var, 
                         "MDA_importance_permutation_counts_variables" = MDA_importance_permutation_count_moy_var,
                         "MDA_importance_ MDA_importance_ARI_permutation_variables" =  MDA_importance_ARI_permutation_moy_var)
  class(output_summary) <- "rfclust.summary"
  return(output_summary)

}

