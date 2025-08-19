# estimate G_{ij}(s, t)

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
  # - but we only carea bout G_{i,i}(s,t)
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
