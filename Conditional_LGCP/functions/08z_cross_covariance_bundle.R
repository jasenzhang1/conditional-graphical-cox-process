

steps_78_OG <- function(kl_coeffs, y_c_strata, query_y_c, eigenfunctions, ncores){
  
  #
  # v3, v3, v2
  # 
  p <- dim(kl_coeffs)[2]
  
  # step 6
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  K_c <- construct_kernel_matrix_step_6(y_c_strata, query_y_c, gamma_c)
  
  # steps 7 and 8
  V_YcXij <- construct_cross_covariance_matrix_v3(kl_coeffs, ncores)
  M_hat <- estimate_regression_operators_v3(K_c, V_YcXij, p)
  V_cond <- evaluate_regression_at_query_v2(M_hat, y_c_strata, query_y_c, eigenfunctions, gamma_c, p)  
  
  return(V_cond)
  
}

steps_78_JASA <- function(kl_coeffs, y_c_strata, query_y_c, eigenfunctions, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # 9/9/2025
  #
  # GOAL: obtain conditional covariance immediately, with the JASA routine
  #
  # v4, v3, v3
  #
  # input:
  # 
  # - kl_coeffs             (n_stratum x p x d matrix) tensor of KL coefficients
  # - y_c_strata            (n_stratum x q_c matrix)   matrix of continuous covariates
  # - query_y_c             (q_c dim vector)           vector of query covariates
  # - eigenfunctions        (list of p matrices, each of which is m x d_i)
  # - ncores                (integer)
  #
  #
  # output:
  #
  # - V_cond                (list of length p^2, each element is m x m matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(kl_coeffs)[2]
  
  # step 6
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  K_c <- construct_kernel_matrix_step_6(y_c_strata, query_y_c, gamma_c)
  
  V_YcXij <- construct_cross_covariance_matrix_v4(kl_coeffs, y_c_strata, query_y_c, gamma_c, ncores)
  M_hat <- estimate_regression_operators_v3(K_c, V_YcXij, p)
  V_cond <- evaluate_regression_at_query_v3(M_hat, y_c_strata, query_y_c, eigenfunctions, gamma_c, p)  
  
  return(V_cond)
}



# basis coefficient bundling

basis_coefficient_method_JASA <- function(kl_coeffs, y_c_strata, query_y_c, eigenfunctions, ncores){
  
  # ----------------------------------------------------------------------------
  #
  #
  # 9/11/2025
  #
  # GOAL: obtain conditional precision immediately, with the JASA routine
  #
  #       we use the basis coefficient method to reduce computational complexity
  #
  #
  # input:
  # 
  # - kl_coeffs             (n_stratum x p x d matrix) tensor of KL coefficients
  # - y_c_strata            (n_stratum x q_c matrix)   matrix of continuous covariates
  # - query_y_c             (q_c dim vector)           vector of query covariates
  # - eigenfunctions        (list of p matrices, each of which is m x d_i)
  # - ncores                (integer)
  #
  #
  # output:
  #
  # - V_cond                (list of length p^2, each element is m x m matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  # 1) cross-covariance = w(y_c) * alpha * alpha
  

  n_stratum      <- dim(kl_coeffs)[1] 
  p              <- dim(kl_coeffs)[2] 
  max_components <- dim(kl_coeffs)[3] 
  
  # step 6
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  K_c <- construct_kernel_matrix_step_6(y_c_strata, query_y_c, gamma_c) 
  K_c_reg_inv <- solve_sym(psd_jitter(K_c, pinv_eps = 1e-4))
  
  
  # evaluate the kernel of y_c with every y_c_k (n-dim vec)
  weights <- apply(y_c_strata, 1, function(row) {
    step_6_kernel(as.numeric(row), query_y_c, gamma_c) 
  })  
  
  # Generate (i, j) index pairs where j >= i
  index_pairs <- do.call(rbind, lapply(1:p, function(i) cbind(i, i:p)))
  
  # Parallel computation for each (i, j) pair
  results <- pbmclapply(
    1:nrow(index_pairs),
    function(idx) {
      i <- index_pairs[idx, 1]
      j <- index_pairs[idx, 2]
      
      # Extract A and B matrices
      if (max_components == 1) {
        A <- matrix(kl_coeffs[, i, ], nrow = n_stratum)
        B <- matrix(kl_coeffs[, j, ], nrow = n_stratum)
      } else {
        A <- kl_coeffs[, i, ] # n x d
        B <- kl_coeffs[, j, ] # n x d
      }
      
      # Create n_stratum x d x d array of outer products NOTICE, ITS WEIGHTS * ALPHA_I * ALPHA_J
      V_matrix_ij <- abind(
        lapply(1:n_stratum, function(k) weights[k] * A[k, ] %o% B[k, ]),
        along = 0
      )
      
      # Matrix multiplication with K_inverse to get regression operator
      
      V_matrix_ij_flat <- matrix(V_matrix_ij, nrow = n_stratum, ncol = max_components^2) # n x (d^2)
      result <- K_c_reg_inv %*% V_matrix_ij_flat # (n x n) %*% (n x d^2)
      
      
      # expanding and taking the mean across samples
      M_ij <- array(result, dim = c(n_stratum, max_components, max_components))  
      M_ij_mean <- apply(M_ij, c(2, 3), mean)
      
      
      M_ij_mean
      
    },
    mc.cores = ncores)
  
  
  names(results) <- paste0(index_pairs[,1], '_', index_pairs[,2])
  return(results)
  
  
  
  
  
}