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

estimate_conditional_correlation_v2 <- function(V_conditional, gamma1, p) {
  
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
    for (j in i:p) {
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

estimate_conditional_correlation_v3 <- function(V_conditional, p, pinv_eps = 1e-6) {
  
  # ---------------------------------------------------------------------
  #
  # GOAL: estimate conditional correlation
  #
  # - v3: we adaptively estimate gamma depending on the negative eigenvalues
  # -     use psd_jitter in 00a_matrix_massaging
  # -     when i = j, we get the identity matrix
  #
  # Input: 
  #
  # - V_conditional      (list of length p^2, each element m x m)
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
    for (j in i:p) {
      key <- paste(i, j, sep="_")
      
      if (i == j) {
        # Diagonal elements are identity: m x m
        C_conditional[[key]] <- diag(n_time)
      } else {
        # Off-diagonal correlation
        key_ii <- paste(i, i, sep="_")
        key_jj <- paste(j, j, sep="_")
        
        V_ii <- psd_jitter(V_conditional[[key_ii]])  # m x m
        V_jj <- psd_jitter(V_conditional[[key_jj]])  # m x m  
        V_ij <- V_conditional[[key]]                 # m x m
        
        # Matrix operations: (m x m) %*% (m x m) %*% (m x m) = (m x m)
        V_ii_inv_sqrt <- matrix_inv_sqrt(V_ii)  # m x m
        V_jj_inv_sqrt <- matrix_inv_sqrt(V_jj)  # m x m
        
        C_conditional[[key]] <- V_ii_inv_sqrt %*% V_ij %*% V_jj_inv_sqrt
      }
    }
  }
  
  return(C_conditional)
}

estimate_conditional_correlation_v4 <- function(V_conditional, p, pinv_eps = 1e-6) {
  
  # ---------------------------------------------------------------------
  #
  # GOAL: estimate conditional correlation
  #
  # - v3: we adaptively estimate gamma depending on the negative eigenvalues
  # -     use psd_jitter in 00a_matrix_massaging
  #
  # - v4: 8/14/2025
  # -     we now let the i = j term be calculated like the rest 
  #
  # Input: 
  #
  # - V_conditional      (list of length p^2, each element m x m)
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
    for (j in i:p) {
      key <- paste(i, j, sep="_")
      

      key_ii <- paste(i, i, sep="_")
      key_jj <- paste(j, j, sep="_")
      
      V_ii <- psd_jitter(V_conditional[[key_ii]])  # m x m
      V_jj <- psd_jitter(V_conditional[[key_jj]])  # m x m  
      V_ij <- V_conditional[[key]]                 # m x m
      
      # Matrix operations: (m x m) %*% (m x m) %*% (m x m) = (m x m)
      V_ii_inv_sqrt <- matrix_inv_sqrt(V_ii)  # m x m
      V_jj_inv_sqrt <- matrix_inv_sqrt(V_jj)  # m x m
      
      C_conditional[[key]] <- V_ii_inv_sqrt %*% V_ij %*% V_jj_inv_sqrt

    }
  }
  
  return(C_conditional)
}



correlation_estimation_KL_cov <- function(eigendecomp, KL_cov){
  
  # ----------------------------------------------------------------------------
  # 
  #
  # GOAL: construct the correlation operator from the KL covariance method
  #
  # input:
  #
  # - eigendecomp (list of 3 entries)
  #
  #   - eigenvalues
  #   - eigenfunctions
  #   - n_dims
  #
  # - KL_cov (list of i_j entries, each is d x d matrix) each value represents covariance between KL coeffs of components a and b in process i and j
  #
  #
  # output:
  #
  # - C_cond   (list of i_j entries)
  #
  #
  # ----------------------------------------------------------------------------
  
  
  p <- length(eigendecomp[[1]])
  d <- dim(KL_cov[[1]])[1]
  m <- dim(eigendecomp[[2]][[1]])[1]
  Delta <- 1/m
  
  C_cond <- list()
  
  for(i in 1:p){
    for(j in i:p){
      key <- paste0(i, '_', j)
      
      eval_i <- eigendecomp[[1]][[i]]
      eval_j <- eigendecomp[[1]][[j]]
      
      evec_i <- eigendecomp[[2]][[i]]
      evec_j <- eigendecomp[[2]][[j]]
      
      d_i <- dim(evec_i)[2]
      d_j <- dim(evec_j)[2]
      
      cov_ij <- KL_cov[[key]]
      
      
      # cov / sqrt(var_i * var_j)
      coeffs <- cov_ij / sqrt(tcrossprod(eval_i, eval_j))
      
      if (any(is.nan(coeffs))) {
        stop("correlation construction process contains NaN values!")
      }      
      
      
      C_ij <- matrix(0, nrow = m, ncol = m)
      
      for(a in 1:d_i){
        for(b in 1:d_j){
          
          # REMEMBER, we have eigenvectors defined as Delta * (v^\top v) = 1. 
          # to keep them in the 
          C_ij <- C_ij + coeffs[a,b] * Delta * tcrossprod(evec_i[, a], evec_j[, b])
          
        }
      }
      
      C_cond[[key]] <- C_ij
      
    }
  }
  
  return(C_cond)
  
  
}
