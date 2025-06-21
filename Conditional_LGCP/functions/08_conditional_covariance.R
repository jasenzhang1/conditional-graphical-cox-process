evaluate_kernel_weights_at_query <- function(Y_continuous_stratum, query_y_c, gamma_c) {
  
  #
  #
  # Compute kernel weights for a query point
  n_stratum <- nrow(Y_continuous_stratum)
  kernel_weights <- numeric(n_stratum)
  
  for (k in 1:n_stratum) {
    y_k <- Y_continuous_stratum[k, ]  # q_c x 1
    kernel_weights[k] <- rbf_kernel(query_y_c, y_k, gamma_c)  # scalar
  }
  
  return(kernel_weights)
}

evaluate_regression_at_query <- function(M_hat, Y_continuous_stratum, query_y_c, 
                                         eigenfunctions, gamma_c, p) {
  # Input: M_hat (list), Y_continuous_stratum (n_stratum x q_c)
  #        query_y_c (q_c x 1), eigenfunctions (list of p matrices m x d)
  # Output: V_conditional (list of length p^2, each element is m x m matrix)
  
  max_components <- ncol(eigenfunctions[[1]])  # d
  n_time <- nrow(eigenfunctions[[1]])          # m
  n_stratum <- nrow(Y_continuous_stratum)
  
  # Compute kernel weights: n_stratum x 1 vector
  kernel_weights <- evaluate_kernel_weights_at_query(Y_continuous_stratum, query_y_c, gamma_c)
  
  V_conditional <- list()
  
  for (i in 1:p) {
    for (j in 1:p) {
      key <- paste(i, j, sep="_")
      
      # Construct conditional covariance operator: m x m matrix
      V_cond_ij <- matrix(0, nrow=n_time, ncol=n_time)
      
      for (a in 1:max_components) {
        for (b in 1:max_components) {
          # Evaluate regression coefficient at query point: scalar
          M_ab_at_query <- sum(kernel_weights * M_hat[[key]][, a, b])
          
          eta_ia <- eigenfunctions[[i]][, a]  # m x 1 vector
          eta_jb <- eigenfunctions[[j]][, b]  # m x 1 vector
          # Tensor product: (m x 1) %*% (1 x m) = (m x m)
          tensor_prod <- tensor_product(eta_ia, eta_jb)
          # Weighted sum: scalar * (m x m) + (m x m) = (m x m)
          V_cond_ij <- V_cond_ij + M_ab_at_query * tensor_prod
        }
      }
      
      V_conditional[[key]] <- V_cond_ij
    }
  }
  
  return(V_conditional)
}