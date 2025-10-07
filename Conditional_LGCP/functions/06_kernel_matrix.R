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

construct_kernel_matrix_step_6 <- function(Y_continuous_stratum, y_query, gamma_c) {
  
  # ----------------------------------------------------------------------------
  # 
  # 
  # GOAL: construct K_c^y_d matrix
  #
  #       recall that it is k(. , y_c_k) * k(. , y_c_k)
  #       so we need to input y_query with every combination of y_c_k 
  # 
  # Input:
  # 
  # - Y_continuous_stratum    (n_stratum x q_c matrix)
  # - y_query                 (q_c dim vector)
  # - gamma_c                 (scalar)
  #
  #
  # Output: 
  #
  # - K_c (n_stratum x n_stratum matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  n_stratum <- nrow(Y_continuous_stratum)
  K_c <- matrix(0, nrow=n_stratum, ncol=n_stratum)
  
  for (i in 1:n_stratum) {
    for (j in i:n_stratum) {
      
      # Compute RBF kernel between two q_c-dimensional vectors
      y_i <- Y_continuous_stratum[i, ]  # q_c x 1 vector
      y_j <- Y_continuous_stratum[j, ]  # q_c x 1 vector
      
      K_c[i, j] <- step_6_kernel(y_i, y_query, gamma_c) * step_6_kernel(y_j, y_query, gamma_c)
      K_c[j, i] <- K_c[i, j]
    }
  }
  
  return(K_c)
}

construct_kernel_matrix_mice <- function(Y_continuous_stratum, y_query, gamma_time = 1/500, gamma_week = 1){
  
  # ----------------------------------------------------------------------------
  # 
  # 
  # GOAL: construct K_c^y_d matrix
  #
  #       recall that it is k(. , y_c_k) * k(. , y_c_k)
  #       so we need to input y_query with every combination of y_c_k 
  #       Y_continuous_stratum is [week #, minutes_elapsed]
  # 
  # Input:
  # 
  # - Y_continuous_stratum    (n_stratum x 2 matrix)
  #   - first column = week #
  #   - second column = minutes_elapsed
  # - y_query                 (q_c = 2-dim vector)
  # - gamma_time              (scalar)
  # - gamma_week              (scalar)
  #
  #
  # Output: 
  #
  # - K_c (n_stratum x n_stratum matrix)
  #
  #
  # ----------------------------------------------------------------------------  
  
  n_stratum <- nrow(Y_continuous_stratum)
  K_c <- matrix(0, nrow=n_stratum, ncol=n_stratum)
  
  for (i in 1:n_stratum) {
    for (j in i:n_stratum) {
      
      # Compute RBF kernel between two q_c-dimensional vectors
      y_i <- Y_continuous_stratum[i, ]  # 2 x 1 vector
      y_j <- Y_continuous_stratum[j, ]  # 2 x 1 vector
      
      # if the weeks are the same, focus on their time-disparity
      if(y_i[1] == y_j[1]){
        K_c[i, j] <- step_6_kernel(y_i[2], y_query[2], gamma_time) * step_6_kernel(y_j[2], y_query[2], gamma_time)
        K_c[j, i] <- K_c[i, j]        
      } else{ # if they are not the same week, focus on their week
        K_c[i, j] <- step_6_kernel(y_i[1], y_query[1], gamma_week) * step_6_kernel(y_j[1], y_query[1], gamma_week)
        K_c[j, i] <- K_c[i, j]        
      }
    }
  }
  
  return(K_c)
  
}

construct_kernel_matrix_v2 <- function(alpha_hat, Y_c_stratum, gamma_c) {
  
  #
  # GOAL: get K_inv for step 6 of JASA implementation (not used)
  #
  # 
  # input:
  # 
  # - alpha_hat
  # - Y_c_stratum
  # - gamma_c
  #
  #
  #
  
  # Input dimensions
  n_stratum <- dim(alpha_hat)[1]  # n_y_d
  p <- dim(alpha_hat)[2]  # p
  d <- dim(alpha_hat)[3]  # d
  m <- dim(eta_hat)[2]  # m (from eta_hat)
  q_c <- ncol(Y_c_stratum)  # q_c
  
  cat("Step 6 input dimensions:\n")
  cat("  n_stratum:", n_stratum, "\n")
  cat("  p:", p, "\n")
  cat("  d:", d, "\n")
  cat("  m:", m, "\n")
  cat("  q_c:", q_c, "\n")
  
  # Compute kernel matrix K_c: n_y_d x n_y_d
  K_c <- matrix(0, n_stratum, n_stratum)
  for (i in 1:n_stratum) {
    for (j in 1:n_stratum) {
      diff <- Y_c_stratum[l, ] - Y_c_stratum[k, ]
      K_c[i, j] <- exp(-gamma_c * sum(diff^2))
    }
  }
  K_hat <- K_c / n_stratum  # n_y_d x n_y_d
  K_reg <- K_hat + gamma_c * diag(n_stratum)  # n_y_d x n_y_d
  K_inv <- solve(K_reg)  # n_y_d x n_y_d
  
  cat("Step 6 output: K_inv dimension", paste(dim(K_inv), collapse = " x "), "\n")
  
  return(K_inv)
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
