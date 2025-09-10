

steps_78_OG <- function(kl_coeffs, y_c_strata, query_y_c, eigenfunctions, ncores){
  
  #
  # v3, v3, v2
  # 
  p <- dim(kl_coeffs)[2]
  
  # step 6
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  K_c <- construct_kernel_matrix_step_6(y_c_strata, gamma_c)
  
  # steps 7 and 8
  V_YcXij <- construct_cross_covariance_matrix_v3(kl_coeffs, y_c_strata, query_y_c, ncores)
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
  K_c <- construct_kernel_matrix_step_6(y_c_strata, gamma_c)
  
  V_YcXij <- construct_cross_covariance_matrix_v4(kl_coeffs, y_c_strata, query_y_c, ncores)
  M_hat <- estimate_regression_operators_v3(K_c, V_YcXij, p)
  V_cond <- evaluate_regression_at_query_v3(M_hat, y_c_strata, query_y_c, eigenfunctions, gamma_c, p)  
  
  return(V_cond)
}
