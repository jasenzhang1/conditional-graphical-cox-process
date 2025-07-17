generate_covariance_matrix <- function(time_grid, lengthscale = 1.0, variance = 1.0) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Generate Matérn covariance matrix for GP
  #
  # - covariance matrix uses nu = 1/2 to get an exponential kernel
  # 
  # Input: 
  #
  # - time_grid      (m x 1 vector of times)
  # - lengthscale    (scalar)
  # - variance       (scalar)
  #
  # 
  # Output: 
  #
  # - K (m x m covariance matrix)
  #
  #
  # - --------------------------------------------------------------------------
  
  m <- length(time_grid)
  K <- matrix(0, m, m)
  
  # Exponential covariance: K(s,t) = variance * exp(-|s-t|/lengthscale)
  for (i in 1:m) {
    for (j in i:m) {
      r <- abs(time_grid[i] - time_grid[j])
      K[i, j] <- variance * exp(-r / lengthscale)
      K[j, i] <- K[i, j]
    }
  }
  
  # Add small regularization for numerical stability
  return(K + 1e-6 * diag(m))
}



sample_conditional_precision <- function(precision_spec, y_c, time_grid) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Sample precision operator for given conditioning values
  #
  # - y_d's effect is from the underlying adjacency matrix
  #
  # Input: 
  #
  # - precision_spec (list; from generate_precision_operators in code 23)
  #   - adjacency = adj_matrix,              # E_{y_d}, 0's on diag
  #   - beta_coefficients = beta_coeffs,     (p x p x q_c matrix)
  #   - signal_strength = theta,             theta
  #   -  dependence_type = dependence_type,   # Type of h_ij(y_c)
  #   - p = p,
  #   - q_c = q_c
  # 
  # - y_c            (q_c x 1 vector of continuous covariates)
  # - time_grid (m x 1 vector)
  #
  #
  # Output: 
  #
  # - P_block (p x p x m x m array) - precision operator
  #
  #
  # ----------------------------------------------------------------------------
  
  p <- precision_spec$p
  m <- length(time_grid)
  theta <- precision_spec$signal_strength
  
  # Base covariance and precision for temporal structure
  base_cov <- generate_covariance_matrix(time_grid, lengthscale = 2.0, variance = 1.0) # (m x m matrix)
  base_precision <- solve(base_cov)  # K_base^{-1}
  
  
  # Initialize precision operator: P (p x p x m x m)
  P_block <- array(0, dim = c(p, p, m, m))
  
  # Diagonal blocks: [P]_{i,i} = base_precision
  for (i in 1:p) {
    P_block[i, i, , ] <- base_precision
  }
  
  # Off-diagonal blocks for edges: [P]_{i,j} = theta * h_ij(y_c) * base_precision
  
  # for each i < j edge: 
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      if (precision_spec$adjacency[i, j] == 1) {
        # Extract beta_ij coefficients
        beta_ij <- precision_spec$beta_coefficients[i, j, ]
        
        # Compute conditional dependence strength h_ij(y_c)
        h_val <- conditional_dependence_function(
          y_c, beta_ij, precision_spec$dependence_type
        )
        
        # Off-diagonal precision elements
        off_diag_strength <- theta * h_val # scalar x scalar 
        off_diag_precision <- off_diag_strength * base_precision  # scalar * (m x m matrix)
        
        P_block[i, j, , ] <- off_diag_precision
        P_block[j, i, , ] <- off_diag_precision  # Symmetry
      }
    }
  }
  
  return(P_block)
}

assemble_block_matrix <- function(precision_operators){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: assemble block matrix
  #
  #
  # input:
  # 
  # - precision_operators (p x p x m x m)
  #
  # 
  # output:
  # 
  # - precision_matrix (pm x pm)
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(precision_operators)[1]
  m <- dim(precision_operators)[3]
  
  precision_matrix <- matrix(0, p*m, p*m)
  
  for(i in 1:p){
    i_start <- (i-1)*m + 1
    i_end <- i*m 
    for(j in i:p){
      j_start <- (j-1)*m + 1
      j_end <- j*m
      
      precision_matrix[i_start:i_end, j_start:j_end] <- precision_operators[i, j, ,]
      precision_matrix[j_start:j_end, i_start:i_end] <- precision_operators[i, j, ,]
    }
  }
  
  return(precision_matrix)
  
}

generate_log_intensity_functions <- function(precision_operators, time_grid, sample_mode = 'simple',
                                             baseline_mean = 0, seed = NULL) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Generate log-intensity functions from conditional Gaussian process
  #
  # 
  # Input: 
  #
  # - precision_operators (p x p x m x m)
  # - time_grid (m x 1)
  # - baseline_mean
  # - seed
  #
  #
  # Output:
  # 
  # - X_functions (p x m matrix) - log-intensity functions
  #
  #
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  p <- dim(precision_operators)[1]
  m <- length(time_grid)
  
  
  
  # For simplicity, generate independent marginals
  # In full implementation, would use block matrix inversion
  
  X_functions <- array(0, dim = c(p, m))
  
  if(sample_mode == 'simple'){
    
    for (i in 1:p) {
      # Marginal precision for process i: [P]_{i,i}
      marginal_precision <- precision_operators[i, i, , ]
      marginal_cov <- solve(marginal_precision)  # Convert to covariance
      
      # Sample from GP: X_i ~ GP(baseline_mean, marginal_cov)
      X_functions[i, ] <- rmvnorm(1, 
                                  mean = rep(baseline_mean, m), 
                                  sigma = marginal_cov)
    }
  } else{
    precision_matrix <- assemble_block_matrix(precision_operators) #pm x pm 
    
    L_mat <- chol(precision_mat)
    z <- rnorm(p * m)
    mult_GP <- backsolve(L_mat, z)
    
    X_functions <- matrix(mult_GP, nrow = p, ncol = m)
    
  }
  
  return(X_functions)
}

