step_6_kernel <- function(y1, y2, gamma_c) {
  
  # ---------------------------------------------------------------------
  #
  # GOAL: defining the kernel in step 6
  #
  # - for each replicate in strata y_d, obtain its continuous vector and take pairwise kernels
  #
  # Input: 
  #
  # - y1 (q_c x 1)
  # - y2 (q_c x 1)
  # - gamma_c (scalar)
  #
  # 
  # Output: 
  # - kernel value (scalar)
  #
  # ---------------------------------------------------------------------
  
  diff <- y1 - y2              # q_c x 1
  return(exp(-gamma_c * sum(diff^2)))  # scalar
}

KDE_weights <- function(Y_c_k, y_c_query){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: get KDE weights for each Y_c_k value wrt y_c_query
  #
  # inputs:
  # 
  # - Y_c_k (n x q_c dim matrix)
  # - y_c_query (q_c-dim vector)
  #
  #
  # outputs:
  #
  # - weights2 (n-dim vector)
  #
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
  
  
  # Heuristic for selecting gamma_c based on median pairwise distance
  #
  # - faster?
  #
  #
  # Iinput: 
  #
  # - Y_continuous_stratum    (n_stratum x q_c matrix)
  # 
  #
  # Output:
  #
  # - gamma_c
  #
  # 
  # -------------------------
  
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
