

full_conditional_estimation <- function(data_df, y_c_strata, query_y_cs, Tseq, threshold, ncores){
  
  #
  #
  #
  #
  # - data_df         (data.table with 'feature_id', 'time', and 'subject_num')
  # - y_c_strata      (n x q_c matrix of continuous values for each subject)
  # - query_y_cs      (matrix, each row is a y_c query to do estimation on)
  # - Tseq            (values from 0 to Tmax to approximate at)
  # - threshold       (number, HS norm threshold value to recover adj_mat)
  # - ncores

  estimated_graphs <- list()
  
  n <- dim(y_c_strata)[1]

    
  # prepare times
  t6_total <- 0
  t7_total <- 0
  t8_total <- 0 
  
  t0 <- Sys.time()
  
  
  subject_nums <- 1:n
  s_yd <- 1:n
  
  data_df4 <- data_df 

  
  
  
  # size of dataset
  print(paste0('number of subjects: ', length(s_yd)))
  print(paste0('number of spikes: ', nrow(data_df4)))
  
  
  # 6) retroactively get parameters for subjects (ntrain, patient_sel) and neurons (p, feature_sel)
  
  ntrain <- length(unique(data_df4$subject_num))
  p <- length(unique(data_df4$feature_id))
  
  
  patient_sel = unique(data_df4$subject_num) %>% sort()
  feature_sel = unique(data_df4$feature_id) %>% sort()
  
  
  # 7) estimation 
  
  t1 <- Sys.time()
  
  rho_list <- estimate_intensities_stratum_parallel_v3(data_df4, patient_sel, feature_sel, Tseq, ncores)
  t2 <- Sys.time()
  
  recovered_intensities <- do.call(rbind, lapply(rho_list, function(x) x$rho_i)) # p x m matrix of intensities

  time_elapsed <- round(as.numeric(difftime(t2, t1, units = 'mins')), 2)
  print(paste0('checkpoint 1: ', time_elapsed, ' mins'))
  
  # part 3
  g_ij_st <- estimate_covariance_functions_ii(rho_list)
  t3 <- Sys.time()
  
  # part 4
  eigen_decomp <- compute_eigendecomposition_ii(g_ij_st)
  t4 <- Sys.time()
  
  # part 5
  # redundant alphas are set to 0
  kl_coeffs <- estimate_kl_coefficients_parallel(data_df4, eigen_decomp$eigenfunctions, eigen_decomp$n_dims, 
                                                 patient_sel, feature_sel, Tseq, ncores)
  t5 <- Sys.time()
  
  time_elapsed <- round(as.numeric(difftime(t5, t2, units = 'mins')), 2)
  print(paste0('checkpoint 2: ', time_elapsed, ' mins'))    
  
  
  # here, we fit to y_c_strata -----------------------------------------------

    
  t6a <- Sys.time()
  
  # part 6
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  K_c <- construct_kernel_matrix_step_6(y_c_strata, gamma_c) # (n_stratum x n_stratum)
  
  
  t6 <- Sys.time()
  t6_total <- t6_total + round(as.numeric(difftime(t6, t6a, units = 'mins')), 2)
  
  # part 7
  
  # V_YcXij_og <- construct_cross_covariance_matrix(kl_coeffs)
  V_YcXij <- construct_cross_covariance_matrix_v2(kl_coeffs)
  # M_hat_og <- estimate_regression_operators(K_c, V_YcXij_og, gamma_c, p)
  M_hat <- estimate_regression_operators_v3(K_c, V_YcXij, gamma_c, p)
  
  t7 <- Sys.time()
  t7_total <- t7_total + round(as.numeric(difftime(t7, t6, units = 'mins')), 2)
  
  time_elapsed <- round(as.numeric(difftime(t7, t6a, units = 'mins')), 2)
  print(paste0('checkpoint 3: ', time_elapsed, ' mins'))  
  
  # part 8 - parallelize this
  
  cont_inds <- 1:nrow(query_y_cs)
  estimated_graphs_v2 <- pbmclapply(cont_inds, function(cont_ind) {
    
    
    query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()
    
    V_cond <- evaluate_regression_at_query_v2(M_hat, y_c_strata, query_y_c, eigen_decomp$eigenfunctions, gamma_c, p)
    
    
    # part 9
    gamma1 <- 0.01
    C_cond <- estimate_conditional_correlation_v3(V_cond, p)
    
  
    # part 10
    gamma2 <- 0.01
    P_cond <- estimate_precision_operator_v3(C_cond, p)
    
    
    
    # part 11
    #graph_threshold <- select_threshold_by_stability(P_cond, p)
    final_graph_estimates <- estimate_graph(P_cond, C_cond, V_cond, threshold, p)
    
    
    # 8) save graphs (recall that we have ground truths of the form a_b_c)
    
    final_graph_estimates
    
    
  }, mc.cores = ncores) # done with all y_c levels
  
  
  
  names(estimated_graphs_v2) <- apply(query_y_cs, 1, function(x){paste(x, collapse = '_')})
  
  t8 <- Sys.time()
  
  
  time_elapsed <- round(as.numeric(difftime(t8, t7, units = 'mins')), 2)
  print(paste0('checkpoint 4: ', time_elapsed, ' mins'))   
  t8_total <- t8_total + time_elapsed
    

  
  t9 <- Sys.time()
  
  # 9) save runtimes ---------------------------------------------------------
  
  times_end   <- c(t1, t2, t3, t4, t5, t9)
  times_start <- c(t0, t1, t2, t3, t4, t0)    
  times_diff <- difftime(times_end, times_start, units = 'mins') %>% as.numeric() %>% round(2)
  
  times_diff <- c(times_diff[1:5], t6_total, t7_total, t8_total, times_diff[6])
  
  
  
  
  return(list(rho_list = rho_list,                   # list of intensities (rho_i) and bivariate intensities (rho_ii)
              estimates = estimated_graphs_v2, 
              intensities = recovered_intensities,
              run_time = times_diff))

}

conditional_estimation_from_log_intensities <- function(X_k, GP_cov, threshold, ncores){
  
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL:
  # 
  # - we start with ground truth log-intensities and try to recover the precision matrix
  #
  #
  # - X_k    (p x m matrix)     log-intensities of p processes at m timepoints
  # - GP_cov (m x m matrix)     base covariance of the GP
  # - threshold       (number, HS norm threshold value to recover adj_mat)
  # - ncores
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(X_k)[1]
  patient_sel <- 1:p
  
  estimated_graphs <- list()
  
  
  recovered_intensities <- exp(X_k) # p x m matrix of intensities
  
  bivariate_intensities <- lapply(1:nrow(recovered_intensities), function(i) { # rho_ii = rho_i %o% rho_i
    v <- recovered_intensities[i, ]
    v %o% v  * exp(GP_cov) # outer product + GP covariance
  })
  
  # merge rho_i and rho_ii together to rho_list
  rho_list <- lapply(1:p, function(i) {
    list(
      rho_i  = recovered_intensities[i, ],       # m-dim vector
      rho_ii_mat = bivariate_intensities[[i]]   # m x m matrix
    )
  })

  
  # part 3
  g_ij_st <- estimate_covariance_functions_ii(rho_list)

  # part 4
  eigen_decomp <- compute_eigendecomposition_ii(g_ij_st)
  
  # part 4 validation
  g_ij_st_approx <- validate_eigendecomposition_ii(g_ij_st, eigen_decomp)
  validate_eigendecomposition_ii_visualization(g_ij_st, g_ij_st_approx, 1)
  summary(as.numeric(g_ij_st[1,,] - g_ij_st_approx[1,,]))
  
  # part 5
  # redundant alphas are set to 0
  kl_coeffs <- estimate_kl_coefficients_parallel(data_df4, eigen_decomp$eigenfunctions, eigen_decomp$n_dims, 
                                                 patient_sel, feature_sel, Tseq, ncores)
  t5 <- Sys.time()
  
  time_elapsed <- round(as.numeric(difftime(t5, t2, units = 'mins')), 2)
  print(paste0('checkpoint 2: ', time_elapsed, ' mins'))    
  
  
  # here, we fit to y_c_strata -----------------------------------------------
  
  
  t6a <- Sys.time()
  
  # part 6
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  K_c <- construct_kernel_matrix_step_6(y_c_strata, gamma_c) # (n_stratum x n_stratum)
  
  
  t6 <- Sys.time()
  t6_total <- t6_total + round(as.numeric(difftime(t6, t6a, units = 'mins')), 2)
  
  # part 7
  
  # V_YcXij_og <- construct_cross_covariance_matrix(kl_coeffs)
  V_YcXij <- construct_cross_covariance_matrix_v2(kl_coeffs)
  # M_hat_og <- estimate_regression_operators(K_c, V_YcXij_og, gamma_c, p)
  M_hat <- estimate_regression_operators_v3(K_c, V_YcXij, gamma_c, p)
  
  t7 <- Sys.time()
  t7_total <- t7_total + round(as.numeric(difftime(t7, t6, units = 'mins')), 2)
  
  time_elapsed <- round(as.numeric(difftime(t7, t6a, units = 'mins')), 2)
  print(paste0('checkpoint 3: ', time_elapsed, ' mins'))  
  
  # part 8 - parallelize this
  
  cont_inds <- 1:nrow(query_y_cs)
  estimated_graphs_v2 <- pbmclapply(cont_inds, function(cont_ind) {
    
    
    query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()
    
    V_cond <- evaluate_regression_at_query_v2(M_hat, y_c_strata, query_y_c, eigen_decomp$eigenfunctions, gamma_c, p)
    
    
    # part 9
    gamma1 <- 0.01
    C_cond <- estimate_conditional_correlation_v3(V_cond, p)
    
    
    # part 10
    gamma2 <- 0.01
    P_cond <- estimate_precision_operator_v3(C_cond, p)
    
    
    
    # part 11
    #graph_threshold <- select_threshold_by_stability(P_cond, p)
    final_graph_estimates <- estimate_graph(P_cond, C_cond, V_cond, threshold, p)
    
    
    # 8) save graphs (recall that we have ground truths of the form a_b_c)
    
    final_graph_estimates
    
    
  }, mc.cores = ncores) # done with all y_c levels
  
  
  
  names(estimated_graphs_v2) <- apply(query_y_cs, 1, function(x){paste(x, collapse = '_')})
  
  t8 <- Sys.time()
  
  
  time_elapsed <- round(as.numeric(difftime(t8, t7, units = 'mins')), 2)
  print(paste0('checkpoint 4: ', time_elapsed, ' mins'))   
  t8_total <- t8_total + time_elapsed
  
  
  
  t9 <- Sys.time()
  
  # 9) save runtimes ---------------------------------------------------------
  
  times_end   <- c(t1, t2, t3, t4, t5, t9)
  times_start <- c(t0, t1, t2, t3, t4, t0)    
  times_diff <- difftime(times_end, times_start, units = 'mins') %>% as.numeric() %>% round(2)
  
  times_diff <- c(times_diff[1:5], t6_total, t7_total, t8_total, times_diff[6])
  
  
  
  
  return(list(rho_list = rho_list,                   # list of intensities (rho_i) and bivariate intensities (rho_ii)
              estimates = estimated_graphs_v2, 
              intensities = recovered_intensities,
              run_time = times_diff))
  
}