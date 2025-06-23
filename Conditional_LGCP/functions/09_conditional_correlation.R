matrix_inv_sqrt <- function(A, regularization=1e-6) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: robustly take the -1/2 power of a matrix
  #
  # Input: 
  #
  # - A (m x m symmetric positive definite matrix)
  # 
  # 
  # Output: 
  # 
  # - A^{-1/2} (m x m matrix)
  #
  #
  # ---------------------------------------------------------------------
  
  eigen_result <- eigen(A, symmetric=TRUE)
  eigenvals <- eigen_result$values  # m x 1 vector
  eigenvecs <- eigen_result$vectors # m x m matrix
  
  # Regularize small eigenvalues
  eigenvals <- pmax(eigenvals, regularization)  # m x 1
  inv_sqrt_eigenvals <- 1.0 / sqrt(eigenvals)  # m x 1
  
  # Reconstruct: (m x m) %*% diag(m x 1) %*% (m x m) = (m x m)
  return(eigenvecs %*% diag(inv_sqrt_eigenvals) %*% t(eigenvecs))
}


estimate_conditional_correlation <- function(V_conditional, gamma1, p) {
  
  # ---------------------------------------------------------------------
  #
  # GOAL: estimate conditional correlation
  #
  #
  # Input: 
  #
  # - V_conditional      (list of length p^2, each element m x m)
  # - gamma1             (scalar)
  # - p                  (scalar)
  #
  #
  # Output: 
  #
  # - C_conditional      (list of length p^2, each element m x m matrix)
  #
  #
  #------------------------------------------------------------------------
  
  n_time <- nrow(V_conditional[[1]])  # m
  C_conditional <- list()
  
  for (i in 1:p) {
    for (j in 1:p) {
      key <- paste(i, j, sep="_")
      
      if (i == j) {
        # Diagonal elements are identity: m x m
        C_conditional[[key]] <- diag(n_time)
      } else {
        # Off-diagonal correlation
        key_ii <- paste(i, i, sep="_")
        key_jj <- paste(j, j, sep="_")
        
        V_ii <- V_conditional[[key_ii]] + gamma1 * diag(n_time)   # m x m
        V_jj <- V_conditional[[key_jj]] + gamma1 * diag(n_time)   # m x m  
        V_ij <- V_conditional[[key]]                              # m x m
        
        # Matrix operations: (m x m) %*% (m x m) %*% (m x m) = (m x m)
        V_ii_inv_sqrt <- matrix_inv_sqrt(V_ii)  # m x m
        V_jj_inv_sqrt <- matrix_inv_sqrt(V_jj)  # m x m
        
        C_conditional[[key]] <- V_ii_inv_sqrt %*% V_ij %*% V_jj_inv_sqrt
      }
    }
  }
  
  return(C_conditional)
}