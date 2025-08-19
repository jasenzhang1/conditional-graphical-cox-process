library(abind)

construct_cross_covariance_matrix <- function(alpha_hat_stratum) {
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: construct cross-covariance matrix
  #
  #
  # Input: 
  #
  # - alpha_hat_stratum (n_stratum x p x d matrix)
  #
  #
  # Output: 
  #
  # - V_YcXij (list of length p^2, each element is n_stratum x d x d array)
  #
  # ----------------------------------------------------------------------------
  
  n_stratum <- dim(alpha_hat_stratum)[1] # 89
  p <- dim(alpha_hat_stratum)[2]         # 5
  max_components <- dim(alpha_hat_stratum)[3]  # d
  
  V_YcXij <- list()
  
  for (i in 1:p) {
    for (j in 1:p) {
      
      # (n x d) times (n x d) gives us (n x d x d)
      # - outer product among the d dimensions
      # - elementwise concatenation along the n dimension
      
      if(max_components == 1){
        A <- alpha_hat_stratum[ , i ,]
        A <- matrix(A, nrow = length(A))
        B <- alpha_hat_stratum[ , j ,]   
        B <- matrix(B, nrow = length(B))
      } else{
        A <- alpha_hat_stratum[ , i ,]
        B <- alpha_hat_stratum[ , j ,]        
      }

      
      V_matrix <- abind(
        lapply(1:n_stratum, function(i) A[i, ] %o% B[i, ]),
        along = 0
        )   
  
      
      V_YcXij[[paste(i, j, sep="_")]] <- V_matrix
    }
  }
  
  return(V_YcXij)
}

construct_cross_covariance_matrix_v2 <- function(alpha_hat_stratum) {
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: construct cross-covariance matrix
  #
  # - only calculate i_j entries where j >= i to save time
  #
  # Input: 
  #
  # - alpha_hat_stratum (n_stratum x p x d matrix) tensor of KL coefficients
  #
  #
  # Output: 
  #
  # - V_YcXij (list of length p^2, each element is n_stratum x d x d array)
  #
  # ----------------------------------------------------------------------------
  
  n_stratum      <- dim(alpha_hat_stratum)[1] # 89
  p              <- dim(alpha_hat_stratum)[2] # 5
  max_components <- dim(alpha_hat_stratum)[3] # d
  
  V_YcXij <- list()
  
  for (i in 1:p) {
    for (j in i:p) {
      
      # (n x d) times (n x d) gives us (n x d x d)
      # - outer product among the d dimensions
      # - elementwise concatenation along the n dimension
      
      if(max_components == 1){
        A <- alpha_hat_stratum[ , i ,]   # becomes a vector
        A <- matrix(A, nrow = length(A)) # make it a n x 1 matrix
        B <- alpha_hat_stratum[ , j ,]   
        B <- matrix(B, nrow = length(B))
      } else{
        A <- alpha_hat_stratum[ , i ,]  # n x d matrix 
        B <- alpha_hat_stratum[ , j ,]        
      }
      
      
      # note that if i and j are switched, all (d x d) matrices are transposed
      V_matrix <- abind(
        lapply(1:n_stratum, function(i) A[i, ] %o% B[i, ]),
        along = 0
      )   
      
      
      V_YcXij[[paste(i, j, sep="_")]] <- V_matrix
    }
  }
  
  return(V_YcXij)
}

construct_cross_covariance_matrix_v3 <- function(alpha_hat_stratum, ncores) {
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: construct cross-covariance matrix
  #
  # - only calculate i_j entries where j >= i to save time
  # - 8/12/2025 and parallelize!
  #
  # Input: 
  #
  # - alpha_hat_stratum (n_stratum x p x d matrix) tensor of KL coefficients
  #
  #
  # Output: 
  #
  # - V_YcXij (list of length p^2, each element is n_stratum x d x d array)
  #
  # ----------------------------------------------------------------------------
  
  n_stratum      <- dim(alpha_hat_stratum)[1] # 89
  p              <- dim(alpha_hat_stratum)[2] # 5
  max_components <- dim(alpha_hat_stratum)[3] # d
  
  # Generate (i, j) index pairs where j >= i
  index_pairs <- do.call(rbind, lapply(1:p, function(i) cbind(i, i:p)))
  
  # Parallel computation for each (i, j) pair
  results <- mclapply(
    1:nrow(index_pairs),
    function(idx) {
      i <- index_pairs[idx, 1]
      j <- index_pairs[idx, 2]
      
      # Extract A and B matrices
      if (max_components == 1) {
        A <- matrix(alpha_hat_stratum[, i, ], nrow = n_stratum)
        B <- matrix(alpha_hat_stratum[, j, ], nrow = n_stratum)
      } else {
        A <- alpha_hat_stratum[, i, ] # n x d
        B <- alpha_hat_stratum[, j, ] # n x d
      }
      
      # Create n_stratum x d x d array of outer products
      V_matrix <- abind(
        lapply(1:n_stratum, function(k) A[k, ] %o% B[k, ]),
        along = 0
      )
      
      list(name = paste(i, j, sep = "_"), value = V_matrix)
    },
    mc.cores = ncores
  )
  
  # Combine into a named list
  V_YcXij <- setNames(lapply(results, `[[`, "value"),
                      sapply(results, `[[`, "name"))
  
  return(V_YcXij)
}

estimate_regression_operators <- function(K_c, V_YcXij, gamma_c, p) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate regression operator M_{X_{ij}}
  #
  #
  # Input: 
  # 
  # - K_c           (n_stratum x n_stratum)
  # - V_YcXij       (list)
  # - gamma_c       (scalar)
  # - p             (scalar)
  #
  #
  # Output: 
  # 
  # - M_hat (list of length p^2, each element is n_stratum x d x d array)
  #
  # 
  # ----------------------------------------------------------------------------
  
  n_stratum <- nrow(K_c)
  max_components <- dim(V_YcXij[[1]])[2]  # d
  
  # Regularized inverse: (n_stratum x n_stratum)^{-1} = (n_stratum x n_stratum)
  K_c_reg_inv <- solve(K_c + gamma_c * diag(n_stratum))
  
  M_hat <- list()
  
  for (i in 1:p) {
    for (j in 1:p) {
      key <- paste(i, j, sep="_")
      V_matrix <- V_YcXij[[key]]  # n_stratum x d x d
      
      # Initialize regression operator: n_stratum x d x d
      M_hat[[key]] <- array(0, dim=c(n_stratum, max_components, max_components))
      
      # Matrix multiplication for each (a,b) component
      for (a in 1:max_components) {
        for (b in 1:max_components) {
          # (n_stratum x n_stratum) %*% (n_stratum x 1) = (n_stratum x 1)
          M_hat[[key]][, a, b] <- K_c_reg_inv %*% V_matrix[, a, b]
        }
      }
    }
  }
  
  return(M_hat)
}

estimate_regression_operators_v2 <- function(K_c, V_YcXij, gamma_c, p) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate regression operator M_{X_{ij}}
  #
  # - utilize the V_YcXij list where j >= i
  # - (a, b) cannot be truncated, because the elements of the outer product are different
  #
  # Input: 
  # 
  # - K_c           (n_stratum x n_stratum)
  # - V_YcXij       (list)
  # - gamma_c       (scalar)
  # - p             (scalar)
  #
  #
  # Output: 
  # 
  # - M_hat (list of length p^2, each element is n_stratum x d x d array)
  #
  # 
  # ----------------------------------------------------------------------------
  
  n_stratum <- nrow(K_c)
  max_components <- dim(V_YcXij[[1]])[2]  # d
  
  # Regularized inverse: (n_stratum x n_stratum)^{-1} = (n_stratum x n_stratum)
  K_c_reg_inv <- solve(K_c + gamma_c * diag(n_stratum))
  
  M_hat <- list()
  
  for (i in 1:p) {
    for (j in i:p) {
      key_ij <- paste(i, j, sep="_")
      
      V_matrix_ij <- V_YcXij[[key_ij]]  # n_stratum x d x d
      
      # Initialize regression operator: n_stratum x d x d
      M_hat[[key_ij]] <- array(0, dim=c(n_stratum, max_components, max_components))
      
      # Matrix multiplication for each (a,b) component
      for (a in 1:max_components) {
        for (b in 1:max_components) {
          # (n_stratum x n_stratum) %*% (n_stratum x 1) = (n_stratum x 1)
          M_hat[[key_ij]][, a, b] <- K_c_reg_inv %*% V_matrix_ij[, a, b]
        }
      }
    }
  }
  
  return(M_hat)
}

estimate_regression_operators_v3 <- function(K_c, V_YcXij, p) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate regression operator M_{X_{ij}}
  #
  # - utilize the V_YcXij list where j >= i
  # - (a, b) cannot be truncated, because the elements of the outer product are different
  # - v3 = doing (a, b) calculations in parallel
  #
  # - 8/12/2025 - changed gamma_c to be the minimum gamma necessary to be PD (using psd_jitter)
  #
  # Input: 
  # 
  # - K_c           (n_stratum x n_stratum)
  # - V_YcXij       (list)
  # - gamma_c       (scalar)
  # - p             (scalar)
  #
  #
  # Output: 
  # 
  # - M_hat (list of length p^2, each element is n_stratum x d x d array)
  #
  # 
  # ----------------------------------------------------------------------------
  
  n_stratum <- nrow(K_c)
  max_components <- dim(V_YcXij[[1]])[2]  # d
  
  # Regularized inverse: (n_stratum x n_stratum)^{-1} = (n_stratum x n_stratum)
  # K_c_reg_inv <- solve(K_c + gamma_c * diag(n_stratum))
  K_c_reg_inv <- solve(psd_jitter(K_c, pinv_eps = 1e-4))
  
  M_hat <- list()
  
  for (i in 1:p) {
    for (j in i:p) {
      key_ij <- paste(i, j, sep="_")
      
      V_matrix_ij <- V_YcXij[[key_ij]]  # n_stratum x d x d
      
      # Initialize regression operator: n_stratum x d x d
      M_hat[[key_ij]] <- array(0, dim=c(n_stratum, max_components, max_components))
      
      # Matrix multiplication for each (a,b) component (flattened for speed)
      
      V_matrix_ij_flat <- matrix(V_matrix_ij, nrow = n_stratum, ncol = max_components^2) # n x (d^2)
      result <- K_c_reg_inv %*% V_matrix_ij_flat # (n x n) %*% (n x d^2)
      
      M_hat[[key_ij]] <- array(result, dim = c(n_stratum, max_components, max_components))
    }
  }
  
  return(M_hat)
}
