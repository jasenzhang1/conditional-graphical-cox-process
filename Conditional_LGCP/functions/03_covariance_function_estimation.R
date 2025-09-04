# estimate G_{ij}(s, t)

base_covariance_function <- function(rho_ij, rho_i, rho_j, regularization=1e-10){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: fundamental function that computes log ratio of intensities
  #
  #
  # input:
  #
  # - rho_ij (m x m matrix)
  # - rho_i  (m-dim vector)
  # - rho_j  (m-dim vector)
  # - regularization (number)
  #
  # 
  # output:
  # 
  # - G_hat (m x m matrix)  covariance estimate
  #
  # 
  # ----------------------------------------------------------------------------
  

  # Avoid log(0) by adding regularization
  numerator <- pmax(rho_ij, regularization)
  denominator <- pmax(tcrossprod(rho_i, rho_j), regularization)
  
  return(log(numerator / denominator))
  
}

estimate_covariance_functions <- function(rho_hat, rho_hat_pairs, regularization=1e-10) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Want to estimate G_{ij}(s, t)  across all p x p pairs of processes
  #
  # 
  # Input: 
  # 
  # - rho_hat         (p x m matrix)
  # - rho_hat_pairs   (p x p x m x m array)
  #
  #
  # Output: 
  #
  # - G_hat           (p x p x m x m array)
  #
  # ----------------------------------------------------------------------------
  
  p <- nrow(rho_hat)
  n_time <- ncol(rho_hat)
  G_hat <- array(0, dim=c(p, p, n_time, n_time))
  
  for (i in 1:p) {
    for (j in 1:p) {
      for (s_idx in 1:n_time) {
        for (t_idx in 1:n_time) {
          # Extract scalar values at grid points
          numerator <- rho_hat_pairs[i, j, s_idx, t_idx]      # scalar
          denominator <- rho_hat[i, s_idx] * rho_hat[j, t_idx] # scalar * scalar
          
          # Avoid log(0) by adding regularization
          numerator <- max(numerator, regularization)
          denominator <- max(denominator, regularization)
          
          G_hat[i, j, s_idx, t_idx] <- log(numerator / denominator)
        }
      }
    }
  }
  
  return(G_hat)
}


estimate_covariance_functions_ii <- function(rho_i_list, regularization=1e-10) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Want to estimate G_{ij}(s, t)  across all p x p pairs of processes
  #
  # - but we only care about G_{i,i}(s,t)
  #
  # 
  # Input: 
  # 
  # - rho_i_list (list of p entries)
  #   - each entry is a list of 2 matrices
  #     - rho_i      (m-dim vector)            univariate intensity
  #     - rho_ii_mat (m x m dim matrix)        bivariate intensity
  #
  #
  # Output: 
  #
  # - G_hat           (m x m x p array)
  #
  # ----------------------------------------------------------------------------
  
  p <- length(rho_i_list)
  n_time <- length(rho_i_list[[1]]$rho_i)
  G_hat <- array(0, dim=c(n_time, n_time, p))
  
  for (i in 1:p) {

    numerator <- rho_i_list[[i]]$rho_ii_mat                             # m x m
    
    denominator <- outer(rho_i_list[[i]]$rho_i, rho_i_list[[i]]$rho_i)  # m x m
        
    # Avoid log(0) by adding regularization
    numerator <- pmax(numerator, regularization)
    denominator <- pmax(denominator, regularization)
    
    G_hat[, , i] <- log(numerator / denominator)

    
  }
  
  return(G_hat)
}

# section 2: functions to get cross-informed G_ii

# - start with rho's to get G_ii's and G_ij's 

compute_cross_dependency_weights <- function(G_hat) {
  
  #----------------------------------------------------------------------------
  #
  #
  # GOAL: helper function for estimate_covariance_functions_ii_cross_informed
  #       calculate w_ij
  #
  #
  # Input: 
  #
  # - G_hat (p x p x m x m array)
  #
  #
  # Output: 
  #
  # - w_matrix (p x p matrix of dependency weights)
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(G_hat)[1]
  
  w_matrix <- matrix(0, nrow=p, ncol=p)

  
  # calculate hs norm for each i =/= j entry
  for (i in 1:(p-1)) {
    for (j in (1+i):p) {
      G_ij <- G_hat[i, j, , ]
      hs_norm <- hilbert_schmidt_norm(G_ij)
      w_matrix[i, j] <- hs_norm
      
    }
  }
  
  # then normalize by dividing everything by the max
  w_matrix <- w_matrix / max(w_matrix)
  
  return(w_matrix)
}

compute_cross_informed_covariance <- function(G_hat, w_matrix, lambda_cross, p) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: calculate G_tilde_{i,i}(s,t) with information from G_{i,j}(s,t)
  #
  #       recall that G_{i,j}(s,t) was the ratio of bivariate rho and univariate rhos
  #
  # 
  # Input: 
  #
  # - G_hat          (p x p x m x m matrix)   log ratio of bivariate and univariate rho
  # - w_matrix       (p x p matrix)           w_ij values calculated in another function
  # - lambda_cross   (scalar)                 hyperparameter
  # 
  #
  # Output: 
  #
  # - G_tilde (p x m x m matrix) - enhanced marginal covariances
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(G_hat)[1]
  m <- dim(G_hat)[3]  
  
  G_tilde <- array(0, dim=c(p, m, m))
  
  # i = base process, j = other processes that can update it 
  for (i in 1:p) { 
    G_tilde[i, , ] <- G_hat[i, i, , ]
    for (j in 1:p) {
      if (i != j) {
        G_ij <- G_hat[i, j, , ]
        weight <- w_matrix[i, j]
        # Add weighted absolute value of cross-covariance
        G_tilde[i, , ] <- G_tilde[i, , ] + lambda_cross * weight * abs(G_ij)
      }
    }
  }
  
  return(G_tilde)
}

estimate_covariance_functions_ii_cross_informed <- function(rho_list, lambda_cross = 0.1, regularization=1e-10) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Want to estimate G_{i,i}(s,t) for all p processes
  #
  # - 8/27/25 
  #   - we want to use information in G_{i,j}(s,t) to get a better estimate about G_{i,i}(s,t)
  #
  # 
  # Input: 
  # 
  # - rho_list (list of two lists, one for univariate, one for bivariate)
  #     - first list   (values are m-dim vector)           univariate intensity, items are called 'i'
  #     - second list  (values are m x m dim matrix)        bivariate intensity, items are called 'i_j'
  # 
  # - regularization (number)
  #
  #
  # Output: 
  #
  # - G_hat_tilde         (m x m x p array)
  #
  # ----------------------------------------------------------------------------
  
  # 1) parameters
  p <- length(rho_list[[1]])
  n_time <- length(rho_list[[1]][['1']])
  G_hat <- array(0, dim=c(n_time, n_time, p))
  
  w_matrix <- matrix(0, p, p)
  
  correction_matrix <- array(0, dim=c(n_time, n_time, p))
  
  # 2) for process i, find its rho_ii and all rho_ij's 
  for (i in 1:p) {
    ii_key <- paste0(i, '_', i)
    rho_ii <- rho_list[[2]][[ii_key]] # m x m  
    rho_i <- rho_list[[1]][[as.character(i)]] # m-dim vec
    
    # 2.1) calculate base G_ii
    numerator <- pmax(rho_ii, regularization)
    denominator <- pmax(outer(rho_i, rho_i), regularization)
    
    G_hat[,,i] <- log(numerator / denominator) %>% symmetrize()
    
    for (j in i:p) {
      
      if (j != i) {
        ij_key <- paste0(i, '_', j)
        rho_ij <- rho_list[[2]][[ij_key]] # m x m  
        rho_j <- rho_list[[1]][[as.character(j)]] # m-dim vec
        
        # calculate G_ij and G_ji
        numerator <- pmax(rho_ij, regularization)
        denominator <- pmax(outer(rho_i, rho_j), regularization)
        G_ij <- log(numerator / denominator)

        # update correction matrix and w_matrix
        w_matrix[i,j] <- hilbert_schmidt_norm(G_ij)
        
        correction_matrix[,,i] <- correction_matrix[,,i] + w_matrix[i,j] * abs(G_ij)
        correction_matrix[,,j] <- correction_matrix[,,j] + w_matrix[i,j] * abs(t(G_ij)) # G_ji
      }
    }
  }
  
  # finally, get the maximum HS norm and normalize everything by it 
  w_max <- max(w_matrix)
  
  correction_matrix <- (lambda_cross / w_max) * correction_matrix 
  
  G_hat_tilde <- G_hat + correction_matrix
  
  return(G_hat_tilde)
}
