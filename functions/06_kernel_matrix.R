step_6_kernel <- function(y1, y2, gamma_c) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: define the kernel in step 6
  #
  #       - for each replicate in strata y_d, obtain its continuous vector and take pairwise kernels
  #
  #
  # input:
  #
  # - y1               (q_c-dim vector)
  # - y2               (q_c-dim vector)
  # - gamma_c          (scalar)               kernel bandwidth parameter
  #
  #
  # output:
  #
  # - kernel_value     (scalar)               exp(-gamma_c||y1 - y2||^2)
  #
  # ----------------------------------------------------------------------------
  
  diff <- y1 - y2              # q_c x 1
  return(exp(-gamma_c * sum(diff^2)))  # scalar
}

KDE_weights <- function(Y_c_k, y_c_query){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: get normalized KDE weights for each Y_c_k value wrt y_c_query
  #
  #
  # input:
  #
  # - Y_c_k            (n x q_c matrix)       continuous covariates of each replicate
  # - y_c_query        (q_c-dim vector)       query covariate value
  #
  #
  # output:
  #
  # - weights2         (n-dim vector)         kernel weights that sum to 1
  #
  # ----------------------------------------------------------------------------
  
  # 1) get gamma 
  
  gamma_c <- select_gamma_c_bandwidth_v2(Y_c_k)
  
  weights <- apply(Y_c_k, 1, function(row) {
    step_6_kernel(as.numeric(row), y_c_query, gamma_c) 
  })    
  weights2 <- weights / sum(weights) # normalize  
  
  return(weights2)
}



select_gamma_c_bandwidth_v2 <- function(Y_continuous_stratum) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: heuristic for selecting gamma_c based on median pairwise distance
  #
  #       - gamma_c = 1 / median_dist^2, or 1 if there are fewer than 2 replicates
  #         or the median distance is 0
  #
  #
  # input:
  #
  # - Y_continuous_stratum  (n_stratum x q_c matrix)
  #
  #
  # output:
  #
  # - gamma_c               (scalar)               kernel bandwidth parameter
  #
  # ----------------------------------------------------------------------------
  
  q_c <- ncol(Y_continuous_stratum)
  n_stratum <- nrow(Y_continuous_stratum)
  if (n_stratum < 2) {
    return(1.0)
  }
  
  pairs <- combn(n_stratum, 2)  # All unique index pairs (2 x K matrix)
  diffs <- Y_continuous_stratum[pairs[1, ], ] - Y_continuous_stratum[pairs[2, ], ]
  if(q_c == 1){
    distances <- abs(diffs)
  } else{
    distances <- sqrt(rowSums(diffs^2))
  }

  
  median_dist <- median(distances)
  
  if(median_dist > 0){
    return(1.0 / (median_dist^2))
  } else{
    return(1.0)
  }
}
