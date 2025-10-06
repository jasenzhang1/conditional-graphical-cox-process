library(MASS)

exp_kernel <- function(x, y, gamma, variance = 1) variance * exp(- gamma * abs(x - y))
rbf_kernel <- function(x, y, gamma, variance = 1) variance * exp(- gamma * (x - y)^2 )

generate_covariance_matrix <- function(time_grid, kernel = 'rbf', gamma = 1.0, variance = 1.0, nugget = 0) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Generate covariance matrix for GP
  #
  # - 'rbf' = Gaussian
  # - 'rbf_pd' = Gaussian with gamma large enough to be pd
  # - 'exponential' = Matern
  # - 'polynomial'
  # 
  # Input: 
  #
  # - time_grid      (m x 1 vector of times)
  # - gamma          (scalar)
  # - variance       (scalar)
  # - nugget         (scalar)   noise term
  #
  # 
  # Output: 
  #
  # - K (m x m covariance matrix)
  #
  #
  # - --------------------------------------------------------------------------
  
  m <- length(time_grid)
  
  delta_t <- time_grid[2] - time_grid[1]
  
  # gamma <- log(2 * m) * 2 * m  # ideal gamma to generate RBF kernel that is PD
  
  K <- matrix(0, m, m)
  
  
  # Exponential covariance: K(s,t) = variance * exp(-|s-t|/lengthscale)
  for (i in 1:m) {
    for (j in i:m) {
      
      if(kernel == 'rbf'){
        r <- (time_grid[i] - time_grid[j])^2
        val <- variance * exp(- gamma * r)
      } 
      
      else if(kernel == 'rbf_pd'){
        r <- (time_grid[i] - time_grid[j])^2
        val <- variance * exp(- gamma * r)        
      }
      
      else if(kernel == 'exponential'){
        r <- abs(time_grid[i] - time_grid[j])
        val <- variance * exp(-gamma * r)
      } 
      
      else if (kernel == 'polynomial'){
        val <- variance * (1 + time_grid[i] * time_grid[j])^5
      }
      else{
        stop('unsupported kernel')
      }
      
      K[i, j] <- val
      K[j, i] <- val
    }
  }
  
  # nugget term
  
  diag(K) <- diag(K) + nugget
  
  return(K)
}

generate_truncated_covariance_matrix <- function(time_grid, kernel_name, gamma, variance, var_explained){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Generate covariance matrix for GP
  #
  # - but we want to truncate the rbf so that its inverse has a meaningful HS norm 
  #
  # - 'rbf' = Gaussian
  # - 'exponential' = Matern
  # 
  # Input: 
  #
  # - time_grid      (m x 1 vector of times)
  # - gamma          (scalar)
  # - variance       (scalar)
  #
  # 
  # Output: 
  #
  # - K (m x m covariance matrix)
  #
  #
  # - --------------------------------------------------------------------------
  
  kernel_fn <- if(kernel_name=="exp") exp_kernel else rbf_kernel
  
  n <- length(time_grid)
  
  delta_t <- 1/n
  K <- outer(time_grid, time_grid, function(a,b) kernel_fn(a,b, gamma))
  A <- delta_t * K
  ev <- eigen(A, symmetric = TRUE)
  
  cum_frac <- cumsum(ev$values) / sum(ev$values)
  m <- which(cum_frac >= var_explained)[1]
  
  vals <- ev$values[1:m]
  vecs <- ev$vectors[,1:m]
  # truncated inverse of A
  inv_trunc <- vecs %*% diag(1/vals, nrow=m, ncol=m) %*% t(vecs)
  hs <- sqrt(sum(inv_trunc^2))   # Frobenius norm
  
  


  return(list(cov_mat = K, m=n, dt=dt, HS=hs, vals=vals))
}

prec_mat_massager <- function(prec_mat, manual_thresh = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: if we manually generate a precision matrix, we need to:
  #
  # 1) invert it to be a covariance matrix
  # 2) normalize it to be a correlation matrix
  # 3) invert it again to be a standardized precision matrix
  # 4) make it a partial correlation matrix
  #
  #
  #
  # input: 
  #
  # - prec_mat (p x p matrix)  un-normalized precision matrix, where 0's mean no adjacency
  #
  # 
  # output:
  #
  # - list of precision matrix features
  #   - adj_mat         (p x p matrix)    adjacency matrix with 0's on the diagonal. 1 = dependent, 0 = independent
  #   - prec_mat        (p x p matrix)    precision matrix from the massaged correlation matrix
  #   - cor_mat         (p x p matrix)    correlation matrix 
  #   - partial_cor_mat (p x p matrix)    partial correlation matrix
  #   - threshold_p    (number)           smallest non-zero off-diagonal value. Anything less than this threshold will be assumed to be 0
  #
  # ----------------------------------------------------------------------------
  

  
  # 1) start with prec_mat, make it cov_mat, then normalize it to be cor_mat
  cov_mat <- solve_sym(prec_mat)
  cor_mat <- cov_mat %>% cov2cor() %>% sym()
  
  # 2) prec_mat2 = precision matrix from correlation matrix
  prec_mat2 <- solve_sym(cor_mat)
  
  
  # 3) partial correlation 
  D <- diag(1 / sqrt(diag(prec_mat))) 
  partial_cor_mat <- -D %*% prec_mat %*% D  
  diag(partial_cor_mat) <- 1  
  
  # 4) define the threshold as the min value that should not be nonzero
  threshold_p <- min(abs(prec_mat2[prec_mat != 0]))  
  
  # 5) adjacency_matrix takes on a value of 1 if original prec_mat is nonzero
  #    or, we define 1 with an indicator function
  
  if(!is.null(manual_thresh)){
    adj_mat <- matrix(0, p, p)
    adj_mat[prec_mat >= manual_thresh] <- 1
    diag(adj_mat) <- 0       
  } else{
    adj_mat <- matrix(0, p, p)
    adj_mat[prec_mat != 0] <- 1
    diag(adj_mat) <- 0      
  }
  

  
  return(list(adj_mat = adj_mat, 
              prec_mat_og = prec_mat,
              prec_mat = prec_mat2, 
              cor_mat = cor_mat, 
              cov_mat = cov_mat,
              partial_cor_mat = partial_cor_mat, 
              simu_mat = prec_mat2,
              threshold_p = threshold_p))
}

generate_sparse_precision_matrix <- function(y_c_k, p, adj_type, adj_params){
  
  # ----------------------------------------------------------------------------
  # 
  # 7/24/2025
  #
  # GOAL: generate a ground truth sparse precision matrix and correlation matrix 
  #
  # - randomly setting off-diagonal entries to 0 will not guarantee PD!!
  #
  # 
  # input:
  # 
  # - p          (number)             dimension of precision matrix
  # - y_c_k      (q_c dim vector)     continuous covariate vector
  # - adj_type   (string)             type of precision matrix (e.g. banded)
  #                                   check `main_simulation_results_notes.txt` for details
  # - adj_params (vector)             vector of misc parameters
  #
  # 
  # output:
  #
  # - list of precision matrix features
  #   - adj_mat      (p x p matrix)    adjacency matrix with 0's on the diagonal. 1 = dependent, 0 = independent
  #   - prec_mat     (p x p matrix)    precision matrix from the massaged correlation matrix
  #   - cor_mat      (p x p matrix)    correlation matrix 
  #   - threshold_p  (number)          smallest non-zero off-diagonal value. Anything less than this threshold will be assumed to be 0
  # 
  # 
  # ----------------------------------------------------------------------------
  
  prec_mat <- diag(1, p)
  
  if(adj_type == 'banded_c0'){
    
    # adj_params = adj_params = [value = 1, rho = 0.3]
    # 0.3's on off diagonals - constant over time
    # nothing else
    
    rho <- adj_params[2]
    
    prec_mat[row(prec_mat) == col(prec_mat) - 1] <- rho
    prec_mat[row(prec_mat) == col(prec_mat) + 1] <- rho 
    
    result <- prec_mat_massager(prec_mat) # helper function above
  }  
  
  if(adj_type %in% c('banded_c1', 'banded_c2')){
    
    # adj_params = adj_params = [y_min = 0, y_max = 1, rho = 0.3]
    # 0.3's on off diagonals - constant over time
    # nothing else
    
    rho <- adj_params[3]
    
    prec_mat[row(prec_mat) == col(prec_mat) - 1] <- rho
    prec_mat[row(prec_mat) == col(prec_mat) + 1] <- rho 
    
    result <- prec_mat_massager(prec_mat) # helper function above
  }
  
  if(adj_type == 'banded_v1'){
    
    # adj_params = [rho = 0.3]
    
    rho <- adj_params[1]
    
    # 0.3's on off diagonals - constant over time
    # for y_c_k = 1, add nothing
    # for y_c_k = 2, make next set of off-diagonals 0.3
    # for y_c_k = 2, make next set of off-diagonals 0.3
    # etc...
    
    for(i in 1:y_c_k){
      prec_mat[row(prec_mat) == col(prec_mat) - i] <- rho
      prec_mat[row(prec_mat) == col(prec_mat) + i] <- rho     
    }
    
    result <- prec_mat_massager(prec_mat) 
  }
  
  if(adj_type == 'banded_v2'){
    
    # adj_params = [k = 3, cs = 0.4, epsilon = 0.05]
    
    k <- adj_params[1]
    cs <- adj_params[2]
    epsilon <- adj_params[3]
    
    banded_params <- list(k = k)
    alpha_funcs_banded <- create_banded_alpha(p, k = k, covariate_strength = cs)
    
    result_banded <- construct_gershgorin_precision_matrix(
      p, y_c_k, alpha_funcs_banded, 
      structure_type = "banded",
      structure_params = banded_params,
      epsilon = 0.05
    )
    
    result <- prec_mat_massager(result_banded$precision_matrix)
  } 
  
  if(adj_type %in% c('banded_trig', 'banded_trig2')){
    
    # adj_params = [y_min = 0, y_max = 1, rho_max = 0.9]
    #
    # AR(1) precision matrix, where rho(y) = rho_max * cos(2 * pi * y)
    #
    # - this allows us to start with a very positive network, then no network, then very negative, then none, then positive again
    
    rho_max <- adj_params[3]
    rho <- rho_max * cos(2 * pi * y_c_k)

    mat <- matrix(0, p, p)
    for(i in 1:p){
      for(j in i:p){
        mat[i,j] <- rho^(abs(i-j))
        mat[j,i] <- mat[i, j]
      } 
    }
    
    result <- prec_mat_massager(mat, rho^3)
  } 
  
  if(adj_type == 'sparse_v1'){

    
    # utilizes y_c_k 
    
    # adj_params = [s, connection_prob, covariate_strength, epsilon]
    
    s <- adj_params[1] # 10
    cp <- adj_params[2] # 0.01
    cs <- adj_params[3] # 0.5
    epsilon <- adj_params[4] # 0.04
    
    # Sparse structure with s=10 connections per node
    sparse_params <- list(s = s)
    alpha_funcs_sparse <- create_sparse_alpha(p, s = s, connection_prob = cp, covariate_strength = cs)
    
    result_sparse <- construct_gershgorin_precision_matrix(
      p, y_c_k, alpha_funcs_sparse,
      structure_type = "sparse",
      structure_params = sparse_params, 
      epsilon = epsilon
    )
    
    result <- prec_mat_massager(result_sparse$precision_matrix)
  }

    
  return(result)
  

}



sample_conditional_precision_v3 <- function(time_grid, time_grid_est,
                                            base_kernel_params, 
                                            prec_mat_truth, 
                                            y_c) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Sample precision operator for given conditioning values
  #
  #
  # Input: 
  #
  # - time_grid
  # - time_grid_est
  #
  # - base_kernel_params
  #   - base_gamma
  #   - base_kernel
  #   - base_variance
  #   - base_GP_mean
  #
  #
  # - prec_mat_truth (list of data about p x p matrix)
  #   - adj_mat
  #   - prec_mat
  #   - cov_mat
  #   - partial_cor_mat
  #   - threshold_p
  #   - simu_mat (THE PRECISION MATRIX THAT WE WILL BE USING)
  # 
  # - y_c            (q_c x 1 vector of continuous covariates)
  #
  #
  # Output: 
  #
  # - P_block (p x p x m x m array) - precision operator
  #
  #
  # ----------------------------------------------------------------------------
  

  p <- dim(prec_mat_truth$cor_mat)[1]
  m <- length(time_grid)

  
  
  # 1) merge time_grids
  time_grid_both <- union(time_grid, time_grid_est) %>% sort()
  time_grid_idx <- match(time_grid, time_grid_both)
  time_grid_est_idx <- match(time_grid_est, time_grid_both)
  
  # 2.1) get RBF kernel for mxm, m_est x m_est, and m2 x m2
  base_cov_both <- generate_covariance_matrix(time_grid_both, 
                                         kernel = base_kernel_params$base_kernel, 
                                         gamma = base_kernel_params$base_gamma, 
                                         variance = base_kernel_params$base_variance) # (m_est x m_est matrix)
  

  base_cov_both <- psd_jitter(base_cov_both)
  
  base_cov <- base_cov_both[time_grid_idx, time_grid_idx]
  
  base_cov_est <- base_cov_both[time_grid_est_idx, time_grid_est_idx]
  
  # means
  base_mean      <- rep(base_kernel_params$base_GP_mean, length(time_grid))
  base_mean_est  <- rep(base_kernel_params$base_GP_mean, length(time_grid_est))
  base_mean_both <- rep(base_kernel_params$base_GP_mean, length(time_grid_both))
  
  # 2.2) inverse - K_base^{-1}
  base_precision <- solve_sym(base_cov)  # use solve_sym() which takes the inverse then ensures it's symmetric
  base_precision_both <- solve_sym(base_cov_both)
  base_precision_est <- solve_sym(base_cov_est)
  
  # 3) get pm x pm matrices for variance and precision
  #    also get pm-dim mean vector
  
  # GP_simu_var_both  <- kronecker(prec_mat_truth$cor_mat, base_cov_both)
  # GP_simu_prec_both <- kronecker(prec_mat_truth$simu_mat, base_precision_both)
  # GP_simu_mean_both <- rep(base_kernel_params$base_GP_mean, length(time_grid_both))
  # 
  # GP_simu_var  <- kronecker(prec_mat_truth$cor_mat, base_cov)
  # GP_simu_prec <- kronecker(prec_mat_truth$simu_mat, base_precision)
  # GP_simu_mean <- rep(base_kernel_params$base_GP_mean, length(time_grid))
  # 
  # GP_simu_var_est  <- kronecker(prec_mat_truth$cor_mat, base_cov_est)
  # GP_simu_prec_est <- kronecker(prec_mat_truth$simu_mat, base_precision_est)
  # GP_simu_mean_est <- rep(base_kernel_params$base_GP_mean, length(time_grid_est))
  
  # verify if cor_mat and simu_mat are inverses
  #           base_cov_both and base_precision_both are inverses
  # summary(as.numeric(prec_mat_truth$cor_mat - solve(prec_mat_truth$simu_mat)))
  # summary(as.numeric(base_cov_both - solve(base_precision_both)))
  
  
  # 4) statistics to report
  
  num_edges <- sum(prec_mat_truth$adj_mat) / 2
  total_possible_edges <- p * (p-1)/2
  obs_sparsity <- num_edges / total_possible_edges
  
  return(list(P_block_kronecker = list(base_cov = base_cov, # m x m
                                       base_precision = base_precision,  
                                       base_mean = base_mean,
                                       
                                       base_cov_est = base_cov_est, # m_est x m_est
                                       base_precision_est = base_precision_est,
                                       base_mean_est = base_mean_est,
                                       
                                       base_cov_both = base_cov_both, # m2 x m2
                                       base_precision_both = base_precision_both,
                                       base_mean_both = base_mean_both,
                                       
                                       # lists
                                       base_kernel_params = base_kernel_params, 
                                       prec_mat_truth = prec_mat_truth),
 
                                       
                                       # mean, var, prec, for all 3 time_grids
                                       # GP_simu_mean = GP_simu_mean,
                                       # GP_simu_var = GP_simu_var,
                                       # GP_simu_prec = GP_simu_prec,
                                       # GP_simu_mean_est = GP_simu_mean_est,
                                       # GP_simu_var_est = GP_simu_var_est,
                                       # GP_simu_prec_est = GP_simu_prec_est,
                                       # GP_simu_mean_both = GP_simu_mean_both,        
                                       # GP_simu_var_both = GP_simu_var_both,
                                       # GP_simu_prec_both = GP_simu_prec_both),  
              
         
              
              adj_mat = prec_mat_truth$adj_mat, 
              num_edges = num_edges, 
              obs_sparsity = obs_sparsity,
              threshold_p = prec_mat_truth$threshold_p
              ))
}

conditional_precision_before_graph <- function(simu_settings, y_c, time_grid) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Sample precision operator for given conditioning values
  #
  #
  # Input: 
  #
  # - simu_settings (list of the following)
  #
  #   - beta_coefficients
  #   - signal_strength
  #   - dependence_type
  #   - y_c_borders
  #   - threshold, p, q_c
  # 
  # - y_c         (q_c x 1 vector of continuous covariates)
  # - time_grid   (m x 1 vector)
  #
  # Output: 
  #
  # - P_block (p x p x m x m array) - precision operator
  #
  #
  # ----------------------------------------------------------------------------
  
  p <- simu_settings$p
  theta <- simu_settings$signal_strength
  beta_coefficients <- simu_settings$beta_coefficients
  m <- length(time_grid)
  
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
      
      beta_ij <- beta_coefficients[i, j, ]
      
      h_val <- conditional_dependence_function(
        y_c, beta_ij, precision_spec$dependence_type
      )      
      
      prec_ij <- base_precision * theta * h_val
      prec_ij_norm <- sum(prec_ij^2)
      
      
      if (precision_spec$adjacency[[region_id]][i, j] == 1) {
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

assemble_block_matrix_24 <- function(precision_operators){
  
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
      precision_matrix[j_start:j_end, i_start:i_end] <- t(precision_operators[i, j, ,])
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
    precision_matrix <- assemble_block_matrix_24(precision_operators) #pm x pm 
    
    # method 1: full inverse into MVNorm ---------------------------------------
    cov_mat_full <- solve(precision_matrix)
    wow <- cov_mat_full[1:50, 1:50]
    wow2 <- cov_mat_full[51:100, 51:100]
    
    cov_mat_full <- (cov_mat_full + t(cov_mat_full)) / 2
    cov_mat_full <- psd_jitter(cov_mat_full)
    
    sample_vec <- mvrnorm(1, mu = rep(baseline_mean, p*m), Sigma = cov_mat_full)
    X_functions <- matrix(sample_vec, nrow = p, ncol = m, byrow = TRUE)
    
    
    # method 2: cholesky -------------------------------------------------------
    precision_matrix_sym <- psd_jitter(precision_matrix)
    
    # case when min eigenvalue is negative
    
    L_mat <- chol(precision_matrix_sym)
    z <- rnorm(p * m)
    mult_GP <- backsolve(L_mat, z)
    
    X_functions <- matrix(mult_GP, nrow = p, ncol = m)
    
  }
  
  return(X_functions)
}

generate_log_intensity_functions_full_mat <- function(cov_matrix_both, time_grid, time_grid_est, 
                                                      sample_mode = 'simple',
                                                      baseline_mean, seed = NULL) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Generate log-intensity functions from conditional Gaussian process
  #
  # 
  # Input: 
  #
  # - cov_matrix_both (pm2 x pm2 matrix)
  # - time_grid       (m x 1)
  # - time_grid_est   (m_est vector)
  # - baseline_mean   (m x 1 vector)
  # - seed            (number)
  #
  #
  # Output:
  # 
  # - X_functions      (p x m matrix) - log-intensity functions
  #
  #
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  time_grid_both <- sort(union(time_grid, time_grid_est))
  time_grid_idx <- match(time_grid, time_grid_both)
  time_grid_est_idx <- match(time_grid_est, time_grid_both)
  
  m <- length(time_grid_both)
  
  
  
  p <- dim(cov_matrix_both)[1] / m
  
  
  
  
  # For simplicity, generate independent marginals
  # In full implementation, would use block matrix inversion
  
  X_functions <- array(0, dim = c(p, m))
  
  if(sample_mode == 'simple'){
    
    for (i in 1:p) {
      
      start_ind <- (i-1) * m + 1
      end_ind <- i*m
      
      # Marginal precision for process i: [P]_{i,i}
      marginal_cov <- cov_matrix_both[start_ind:end_ind, start_ind:end_ind]

      
      # Sample from GP: X_i ~ GP(baseline_mean, marginal_cov)
      X_functions[i, ] <- rmvnorm(1, 
                                  mean = baseline_mean, 
                                  sigma = marginal_cov)
    }
  } else{
    
    # method 1: full inverse into MVNorm ---------------------------------------

    sample_vec <- mvrnorm(1, mu = rep(baseline_mean, p), Sigma = cov_matrix_both)
    X_functions <- matrix(sample_vec, nrow = p, ncol = m, byrow = TRUE)
    
    X_functions_base <- X_functions[, time_grid_idx]
    X_functions_est <- X_functions[, time_grid_est_idx]
  }
  
  return(list(X_functions, X_functions_base, X_functions_est))
}

