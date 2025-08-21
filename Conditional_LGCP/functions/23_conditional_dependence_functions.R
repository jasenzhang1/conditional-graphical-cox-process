

conditional_dependence_function <- function(y_c, beta_ij, type = "linear") {
  
  
  # ----------------------------------------------------------------------------
  #
  # 
  # GOAL: Compute conditional dependence h_ij(y_c)
  #
  # - we offer linear, quadratic, and non-linear dependencies
  # 
  # 
  # Input: 
  #
  # - y_c        (q_c x 1 - dim vector)      vector of continuous covariates
  # - beta_ij    (q_c x 1 - dim vector)      vector of beta coefficients
  # - type       (string)                    either "linear", "quadratic", or "nonlinear"
  #
  # 
  #  Output: 
  #
  # - scalar dependence strength
  #
  #
  # ----------------------------------------------------------------------------
  
  switch(type,
         "constant" = {
           # h_ij(y_c) = beta_ij[1]
           beta_ij[1]
         },
         "linear" = {
           # h_ij(y_c) = 1 + beta_ij^T * y_c
           1 + sum(beta_ij * y_c)
         },
         "quadratic" = {
           # h_ij(y_c) = 1 + beta_ij^T * y_c + 0.1 * ||y_c||^2
           1 + sum(beta_ij * y_c) + 0.1 * sum(y_c^2)
         },
         "nonlinear" = {
           # h_ij(y_c) = 1 + sin(beta_ij^T * y_c) * exp(-||y_c||^2/2)
           linear_term <- sum(beta_ij * y_c)
           1 + sin(linear_term) * exp(-sum(y_c^2)/2)
         },
         1  # Default: no dependence
  )
}


generate_precision_operators <- function(adj_matrix, theta, q_c, 
                                         dependence_type = "linear", seed = NULL) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Generate true precision operators for all discrete strata
  #
  # - here, we just sample beta coefficients uniformly 
  # 
  # Input: 
  #
  # - adj_matrix  (p x p matrix)    matrix of 1's and 0's, with 0 on diags
  # - theta       (number)          signal strength
  # - q_c         (integer)         dimension of continuous covariates
  #
  # 
  # Output: 
  #
  # - specification for precision operators
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  p <- nrow(adj_matrix)
  
  # Generate regression coefficients beta_ij for continuous dependence
  # beta_ij ~ Uniform([-0.5, 0.5]^{q_c}) for each edge (i,j)
  beta_coeffs <- array(0, dim = c(p, p, q_c))
  for (i in 1:p) {
    for (j in 1:p) {
      if (adj_matrix[i, j] == 1) {
        beta_coeffs[i, j, ] <- runif(q_c, -0.5, 0.5)  # why are beta_ij and beta_ji different?
      }
    }
  }
  
  return(list(
    adjacency = adj_matrix,              # E_{y_d}             (p x p adj mat)
    beta_coefficients = beta_coeffs,     # {beta_ij}           (p x p x q_c)
    signal_strength = theta,             # theta               (scalar)
    dependence_type = dependence_type,   # Type of h_ij(y_c)   (string)
    p = p,                               # number of processes (scalar)
    q_c = q_c                            # dim cont. var       (scalar)
  ))
}

permute_indices <- function(sizes) {
  
  # if I supply c(2,3), return a list of strings:
  #
  # - "1_1", "1_2", "1_3", "2_1", "2_2", "2_3"  
  #
  #
  #
  # sizes (vector) 
  #
  # 
  
  
  # Create a list of 1:x, 1:y, ... for each element
  index_lists <- lapply(sizes, function(n) seq_len(n))
  
  # Create all combinations
  combos <- expand.grid(index_lists)
  
  # Collapse each row into a string
  apply(combos, 1, function(row) paste(row, collapse = "_"))
}

generate_precision_operators_and_matrix <- function(p, theta, q_c, 
                                                    y_c_borders,
                                                    threshold,
                                                    dependence_type = "constant", seed = NULL) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Generate true precision operators for all discrete strata
  #
  # - here, we just sample beta coefficients uniformly from [-0.5, 0.5]
  # 
  # Input: 
  #
  # - p           (number)          dimension of process
  # - theta       (number)          signal strength
  # - q_c         (integer)         dimension of continuous covariates
  # - y_c_borders (list of q_c vectors)   each list tells us the bordering values to create new ground truth graphs
  # - threshold   (number)          value of h_ij for which we decide if (i, j) is an edge
  # 
  # Output: 
  #
  # - specification for precision operators
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  
  # Generate regression coefficients beta_ij for continuous dependence
  # beta_ij ~ Uniform([-0.5, 0.5]^{q_c}) for each edge (i,j)
  beta_coeffs <- array(0, dim = c(p, p, q_c))
  for (i in 1:p) {
    for (j in 1:p) {
      beta_coeffs[i, j, ] <- runif(q_c, -0.5, 0.5)  # why are beta_ij and beta_ji different?
    }
  }
  
  if(dependence_type == 'constant'){ # used for the base case of E_{y_c, y_d} is constant over time
    
    adj_mat <- matrix(0, p, p)
    
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        c_ij <- beta_coeffs[i, j, 1] # between -0.5 ,0.5
        
        if(abs(c_ij) > threshold){
          adj_mat[i,j] <- 1
        }
      }
    }
    
    adj_mat <- adj_mat + t(adj_mat)
    
    adj_mat_list <- list()
    
    num_regions <-  sapply(y_c_borders, length) + 1 # 
    
    regions <- permute_indices(num_regions)
    
    for(i in regions){
      adj_mat_list[[i]] <- adj_mat
    }
    
    
  }

  
  return(list(
    adjacency = adj_mat_list,            # list of E_{y_d}     (list of p x p adj mat)
    beta_coefficients = beta_coeffs,     # {beta_ij}           (p x p x q_c)
    signal_strength = theta,             # theta               (scalar)
    dependence_type = dependence_type,   # Type of h_ij(y_c)   (string)
    p = p,                               # number of processes (scalar)
    q_c = q_c                            # dim cont. var       (scalar)
  ))
}


collect_beta_and_parameters <- function(p, time_grid, time_grid_est, theta, q_c, 
                                        y_c_borders,
                                        sparsity,
                                        base_kernel_params,
                                        dependence_type, 
                                        adj_type, 
                                        seed = NULL) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Generate beta coefficients and collect all model params in a list
  #
  # - here, we just sample beta coefficients uniformly from [-0.5, 0.5]
  # 
  # Input: 
  #
  # - p                  (number)          dimension of process
  # - time_grid          (number)          time-grid points
  # - time_grid_est      (vector) 
  # - theta              (number)          signal strength
  # - q_c                (integer)         dimension of continuous covariates
  # - y_c_borders        (list of q_c vectors)   each list tells us the bordering values to create new ground truth graphs
  # - sparsity           (number)           proportion of edges
  # 
  # - base_kernel_params (list)            list of mxm kernel params
  #   - base_gamma       (number)          gamma param for RBF kernel
  #   - base_kernel      (string)          name of kernel such as 'RBF'
  #   - base_variance    (number)          variance parameter in kernel
  #   - base_GP_mean     (number)          constant mean value in GP
  #
  # - dependence_type    (string)          how do the continuous covariates affect ground truth?
  # - adj_type           (string)          pxp adjacency matrix and how it changes via covariates
  # - seed               (number)          simulation seed number
  # 
  # 
  # Output: 
  #
  # - specification for precision operators
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  m <- length(time_grid)
  
  # Generate regression coefficients beta_ij for continuous dependence
  # beta_ij ~ Uniform([-0.5, 0.5]^{q_c}) for each edge (i,j)
  beta_coeffs <- array(0, dim = c(p, p, q_c))
  for (i in 1:p) {
    for (j in 1:p) {
      beta_coeffs[i, j, ] <- runif(q_c, -0.5, 0.5)  # why are beta_ij and beta_ji different?
    }
  }
  
  
  return(list(
    beta_coefficients = beta_coeffs,     # {beta_ij}           (p x p x q_c)
    signal_strength = theta,             # theta               (scalar)
    dependence_type = dependence_type,   # Type of h_ij(y_c)   (string)
    adj_type = adj_type,                 #                     (string)
    time_grid = time_grid,               # time grid
    time_grid_est = time_grid_est,
    base_kernel_params = base_kernel_params,
    y_c_borders = y_c_borders,
    sparsity = sparsity,
    p = p,                               # number of processes (scalar)
    q_c = q_c                            # dim cont. var       (scalar)
  ))
}
