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

construct_kernel_matrix_step_6 <- function(Y_continuous_stratum, gamma_c) {
  
  # 
  # GOAL: construct K_c^y_d matrix
  #
  # Input:
  # 
  # - Y_continuous_stratum    (n_stratum x q_c matrix)
  # - gamma_c                 (scalar)
  #
  #
  # Output: 
  #
  # - K_c (n_stratum x n_stratum matrix)
  #
  #
  
  n_stratum <- nrow(Y_continuous_stratum)
  K_c <- matrix(0, nrow=n_stratum, ncol=n_stratum)
  
  for (i in 1:n_stratum) {
    for (j in i:n_stratum) {
      
      # Compute RBF kernel between two q_c-dimensional vectors
      y_i <- Y_continuous_stratum[i, ]  # q_c x 1 vector
      y_j <- Y_continuous_stratum[j, ]  # q_c x 1 vector
      
      K_c[i, j] <- step_6_kernel(y_i, y_j, gamma_c)
      K_c[j, i] <- K_c[i, j]
    }
  }
  
  return(K_c)
}

select_gamma_c_bandwidth <- function(Y_continuous_stratum) {
  
  
  # Heuristic for selecting gamma_c based on median pairwise distance
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
  
  n_stratum <- nrow(Y_continuous_stratum)
  if (n_stratum < 2) {
    return(1.0)
  }
  
  distances <- c()
  for (i in 1:(n_stratum-1)) { # for each unique pair
    for (j in (i+1):n_stratum) {
      diff <- Y_continuous_stratum[i, ] - Y_continuous_stratum[j, ]
      distances <- c(distances, sqrt(sum(diff^2)))  # append their euclidean distance
    }
  }
  
  median_dist <- median(distances)
  
  if(median_dist > 0){
    return(1.0 / (median_dist^2))
  } else{
    return(1.0)
  }
}


