library(abind)

construct_cross_covariance_matrix <- function(alpha_hat_stratum) {
  
  #
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

estimate_regression_operators <- function(K_c, V_YcXij, gamma_c, p) {
  
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