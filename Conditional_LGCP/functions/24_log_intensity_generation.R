library(MASS)

generate_covariance_matrix <- function(time_grid, kernel = 'rbf', gamma = 1.0, variance = 1.0) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Generate covariance matrix for GP
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
  
  return(K)
}

generate_sparse_precision_matrix <- function(p, graph_type = 'banded'){
  
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
  # - p          (number)  dimension of precision matrix
  # - graph_type (string)  type of precision matrix (e.g. banded)
  #
  # 
  # output:
  #
  # - list of precision matrix features
  #   - adj_mat      (p x p matrix)    adjacency matrix with 0's on the diagonal. 1 = dependent, 0 = independent
  #   - prec_mat     (p x p matrix)    precision matrix from the massaged correlation matrix
  #   - cor_at       (p x p matrix)    correlation matrix 
  #   - threshold_p  (number)          smallest non-zero off-diagonal value. Anything less than this threshold will be assumed to be 0
  # 
  # 
  # ----------------------------------------------------------------------------
  
  prec_mat <- diag(1, p)
  adj_mat <- matrix(0, p, p)
  
  if(graph_type == 'banded'){
    prec_mat[row(prec_mat) == col(prec_mat) - 1] <- 0.3
    prec_mat[row(prec_mat) == col(prec_mat) + 1] <- 0.3 
    
    cor_mat <- solve(prec_mat) %>% cov2cor()
    
    prec_mat2 <- solve(cor_mat)
    
    # define the threshold as the min value that should not be nonzero
    threshold_p <- min(abs(prec_mat2[prec_mat != 0]))  
    
    # adjacency_matrix
    adj_mat[prec_mat != 0] <- 1
    diag(adj_mat) <- 0
    
  }
  
  return(list(adj_mat = adj_mat, prec_mat = prec_mat2, cor_mat = cor_mat, threshold_p = threshold_p))
  

}

sample_conditional_precision <- function(precision_spec, y_c, region_id) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Sample precision operator for given conditioning values
  #
  #
  # Input: 
  #
  # - precision_spec (list; from generate_precision_operators in code 23)
  # 
  #   - beta_coefficients   (p x p x q_c matrix)   matrix of coefficients to model continuous covariate effects
  #   - signal_strength     (number)               theta value for signal strength
  #   - dependence_type     (string)               how do continuous covariates affect the model
  #   - time_grid           (m dim vector)         discretized time points
  #   - base_kernel_params  (list)                 m x m kernel parameters
  # 
  #   - y_c_borders 
  #   - sparsity 
  #   - p = p,
  #   - q_c = q_c
  # 
  # - y_c            (q_c x 1 vector of continuous covariates)
  # - region_id (string) name of region
  #
  #
  # Output: 
  #
  # - P_block (p x p x m x m array) - precision operator
  #
  #
  # ----------------------------------------------------------------------------
  
  p <- precision_spec$p
  m <- length(precision_spec$time_grid)
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


sample_conditional_precision_v3 <- function(simu_settings, 
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
  # - simu_settings (list; from generate_precision_operators in code 23)
  #   - beta_coefficients = beta_coeffs,     (p x p x q_c matrix)
  #   - signal_strength = theta,             theta
  #   - dependence_type = dependence_type,   # Type of h_ij(y_c)
  #   - y_c_borders
  #   - sparsity
  #   - p = p,
  #   - q_c = q_c
  #   - time_grid
  #   - base_kernel
  #     - base_gamma
  #     - base_kernel
  #     - base_variance
  #     - base_GP_mean
  #
  # - prec_mat_truth (list of data about p x p matrix)
  #   - adj_mat
  #   - prec_mat
  #   - cov_mat
  #   - threshold_p
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
  
  # load parameters 
  
  threshold <- NA
  threshold_p <- prec_mat_truth$threshold_p
  adj_mat <- prec_mat_truth$adj_mat
  time_grid <- simu_settings$time_grid
  
  p <- simu_settings$p
  m <- length(time_grid)
  theta <- simu_settings$signal_strength
  sparsity <- simu_settings$sparsity
  
  
  # Base covariance and precision for temporal structure (m x m)
  base_cov <- generate_covariance_matrix(time_grid, 
                                         kernel = simu_settings$base_kernel_params$base_kernel, 
                                         gamma = simu_settings$base_kernel_params$base_gamma, 
                                         variance = simu_settings$base_kernel_params$base_variance) # (m x m matrix)
  base_cov <- psd_jitter(base_cov)
  
  base_precision <- solve(base_cov)  # K_base^{-1}
  
  base_precision <- (base_precision + t(base_precision))/2
  
  delta_t <- time_grid[2] - time_grid[1]
  
  
  # create matrix of theta * h_{ij} --------------------------------------------
  
  # H_mat <- matrix(0, p, p)
  # 
  # for(i in 1:p){
  #   for(j in i:p){
  #     
  #     beta_ij <- simu_settings$beta_coefficients[i, j, ]
  #     
  #     h_val <- conditional_dependence_function(
  #       y_c, beta_ij, simu_settings$dependence_type
  #     )
  #     
  #     H_mat[i, j] <- theta * h_val
  #     H_mat[j, i] <- theta * h_val
  #   }
  # }
  # 
  # H_mat_cor <- cov2cor(crossprod(H_mat))        # random correlation matrix
  # H_mat_cor <- (H_mat_cor + t(H_mat_cor))/2
  # H_mat_prec <- solve(H_mat_cor)                # random precision matrix
  # H_mat_prec <- (H_mat_prec + t(H_mat_prec))/2
  # 
  # 
  # # HS norms to get adjacency matrix 
  # adj_mat <- matrix(0, p, p)
  # 
  # HS_norm_mat <- matrix(0, p, p)
  # HS_norms <- c()
  # 
  # base_precision_HS <- sqrt(sum(base_precision^2)) * delta_t  # HS norm = sqrt() * delta_t
  # 
  # 
  # # calculate HS norms for each (i \neq j)
  # for (i in 1:(p-1)) {
  #   for (j in (i+1):p) {
  #     
  #     HS_norm <- base_precision_HS * abs(H_mat_prec[i, j])  # HS norm of c * A = |c| * ||A||_HS
  #     
  #     HS_norms <- c(HS_norms, HS_norm)
  #     
  #     HS_norm_mat[i,j] <- HS_norm
  #     HS_norm_mat[j,i] <- HS_norm
  #   }
  # }
  # 
  # # now, we select which edges are kept with the sparsity parameter
  # 
  # if(is.na(threshold)){
  #   threshold <- quantile(HS_norms, 1 - sparsity) %>% unname() # shouldn't be used during estimation
  #   
  # }
  # 
  # 
  # 
  # 
  # adj_mat <- HS_norm_mat
  # adj_mat[adj_mat < threshold] <- 0
  # adj_mat[adj_mat > 0] <- 1
  # 
  # H_mat_prec_thresh <- H_mat_prec
  # H_mat_prec_thresh[! adj_mat + diag(p)] <- 0
  # 
  # H_mat_cor_thresh <- solve(H_mat_prec_thresh) %>% cov2cor()
  # H_mat_cor_thresh <- (H_mat_cor_thresh + t(H_mat_cor_thresh))/2
  
  GP_simu_var = kronecker(prec_mat_truth$cor_mat, base_cov)
  GP_simu_mean <- rep(base_kernel_params$base_GP_mean, m)
  
  # with our ground truth precision matrix H_mat_prec_thresh, find threshold for which below it are just 0's
  
  # if(is.na(threshold_p)){
  #   threshold_p <- min(abs(H_mat_prec_thresh[H_mat_prec_thresh != 0])) 
  # }
  
  
  # statistics to report
  
  num_edges <- sum(adj_mat) / 2
  total_possible_edges <- p * (p-1)/2
  obs_sparsity <- num_edges / total_possible_edges
  
  return(list(P_block_kronecker = list(base_cov = base_cov,
                                       base_precision = base_precision,  # m x m
                                       base_gamma = base_kernel_params$base_gamma,
                                       base_kernel = base_kernel_params$base_kernel,
                                       base_variance = base_kernel_params$base_variance,
                                       base_GP_mean = base_kernel_params$base_GP_mean,
                                       # H_mat_cor = H_mat_cor,            # p x p 
                                       # H_mat_prec = H_mat_prec,
                                       # H_mat_cor_thresh = H_mat_cor_thresh,
                                       # H_mat_prec_thresh = H_mat_prec_thresh,
                                      
                                       prec_mat = prec_mat_truth$prec_mat,  # p x p
                                       cor_mat  = prec_mat_truth$cor_mat,   # p x p 
                                       GP_simu_mean = GP_simu_mean,         # pm x 1 vec
                                       GP_simu_var = GP_simu_var),          # pm x pm mat
              
         
              
              adj_mat = adj_mat, 
              num_edges = num_edges, 
              target_sparsity = sparsity,
              obs_sparsity = obs_sparsity,
              # threshold = threshold,
              threshold_p = threshold_p
              #HS_norms = HS_norm_mat
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

generate_log_intensity_functions_full_mat <- function(cov_matrix, time_grid, sample_mode = 'simple',
                                                      baseline_mean, seed = NULL) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Generate log-intensity functions from conditional Gaussian process
  #
  # 
  # Input: 
  #
  # - cov_matrix (pm x pm matrix)
  # - time_grid (m x 1)
  # - baseline_mean (m x 1 vector)
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
  
  m <- length(time_grid)
  p <- dim(cov_matrix)[1] / m
  
  
  
  
  # For simplicity, generate independent marginals
  # In full implementation, would use block matrix inversion
  
  X_functions <- array(0, dim = c(p, m))
  
  if(sample_mode == 'simple'){
    
    for (i in 1:p) {
      
      start_ind <- (i-1) * m + 1
      end_ind <- i*m
      
      # Marginal precision for process i: [P]_{i,i}
      marginal_cov <- cov_matrix[start_ind:end_ind, start_ind:end_ind]

      
      # Sample from GP: X_i ~ GP(baseline_mean, marginal_cov)
      X_functions[i, ] <- rmvnorm(1, 
                                  mean = baseline_mean, 
                                  sigma = marginal_cov)
    }
  } else{
    
    # method 1: full inverse into MVNorm ---------------------------------------

    sample_vec <- mvrnorm(1, mu = rep(baseline_mean, p), Sigma = cov_matrix)
    X_functions <- matrix(sample_vec, nrow = p, ncol = m, byrow = TRUE)
    
    
  }
  
  return(X_functions)
}

