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
  
  
  
  if(object[[1]]$distance == "co-clustering"){
    
    pair_appearances <- sum_dist + sum_sim                                          #Somme des occurence des paires dans les forêts.
    similarity_matrix <- sum_sim / pair_appearances
    distance_matrix <- 1 - similarity_matrix
    diag(distance_matrix) <- 0
    distance_matrix[is.na(distance_matrix)] <- 0
    
    
    CSI_ORR <- lapply(object, '[[', 6)
    
    CSI_ARI <- lapply(object, '[[', 7)
    
    #Importance des variables CSI coclust
    CSI_ORR_var <- list()
    for (j in 1:length(CSI_ORR)){
      for (var in names(CSI_ORR[[j]])){
        if (is.null(CSI_ORR_var[[var]])){
          CSI_ORR_var[[var]] <- 0
        }
        CSI_ORR_var[[var]]  <- CSI_ORR_var[[var]] + CSI_ORR[[j]][[var]]
      }
    }
    
    CSI_ORR_moy_var <- lapply(CSI_ORR_var, function(x) x/ntrees)
    
    
    CSI_ARI_var <- list()
    for (j in 1:length(CSI_ARI)){
      for (var in names(CSI_ARI[[j]])){
        if (is.null(CSI_ARI_var[[var]])){
          CSI_ARI_var[[var]] <- 0
        }
        CSI_ARI_var[[var]]  <- CSI_ARI_var[[var]] + CSI_ARI[[j]][[var]]
      }
    }
    
    CSI_ARI_moy_var <- lapply(CSI_ARI_var, function(x) x/ntrees)
    
    output_summary <- list("distance_matrix" = distance_matrix, "MDI_importance_variables" = MDI_importance_moy_var, 
                           "CSI_ORR_variables" = CSI_ORR_moy_var,
                           "CSI_ARI_variables" =  CSI_ARI_moy_var)
    
    
  }else if(object[[1]]$distance == "inertia"){
    sum_present <- ntrees - sum_absent
    distance_matrix <- sum_dist / sum_present
    diag(distance_matrix) <- 0
    distance_matrix[is.na(distance_matrix)] <- 0
    
    CSI_IP <- lapply(object, '[[', 6)
    
    #Importance des variables CSI inertia
    CSI_IP_var <- list()
    for (j in 1:length(CSI_IP)){
      for (var in names(CSI_IP[[j]])){
        if (is.null(CSI_IP_var[[var]])){
          CSI_IP_var[[var]] <- 0
        }
        CSI_IP_var[[var]]  <- CSI_IP_var[[var]] + CSI_IP[[j]][[var]]
      }
    }
    
    CSI_IP_moy_var <- lapply(CSI_IP_var, function(x) x/ntrees)
    
    output_summary <- list("distance_matrix" = distance_matrix, "MDI_importance_variables" = MDI_importance_moy_var, 
                           "CSI_IP_variables" = CSI_IP_moy_var)
  }
  class(output_summary) <- "rfclust.summary"
  return(output_summary)
  
}

