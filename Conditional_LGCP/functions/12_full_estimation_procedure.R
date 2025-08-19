library(grid)

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

full_conditional_estimation_with_truths <- function(dataset, ncores){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: troubleshoot the full conditional estimation procedure with ground truths between steps
  #
  #
  # Input:
  #
  # - data_df         (data.table with 'feature_id', 'time', and 'subject_num')
  # - y_c_strata      (n x q_c matrix of continuous values for each subject)
  # - query_y_cs      (matrix, each row is a y_c query to do estimation on)
  # - Tseq            (values from 0 to Tmax to approximate at)
  # - threshold       (number, HS norm threshold value to recover adj_mat)
  # - ncores
  #
  # 
  #
  # ----------------------------------------------------------------------------
  
  # 0) load  
  
  time_grid <- dataset$time_grid
  time_grid_est <- dataset$time_grid_est
  time_grid_both <- dataset$time_grid_both
  true_graphs <- dataset$true_graphs
  
  data_df <- convert_data_for_estimation(dataset$subject_data) %>% as.data.table()
  data_df$time <- data_df$time / dataset$simulation_params$T_max # normalize to [0, 1]
  
  X_k_truth <- dataset$X_k_truth
  X_k_coarse_truth <- dataset$X_k_coarse_truth
  X_k_both_truth <- dataset$X_k_both_truth
  
  threshold <- true_graphs[[1]]$threshold_p
  cor_mat_pm_truth <- true_graphs[[1]]$P_block_kronecker$GP_simu_var
  base_GP_mean <- true_graphs[[1]]$P_block_kronecker$base_GP_mean
  
  # 1.1) visualize true log-intensities - good
  
  g_11 <- grid.arrange(visualize_log_intensity(X_k_truth[1:5,,1], time_grid, 'Finer Grid'),
                       visualize_log_intensity(X_k_coarse_truth[1:5,,1], time_grid_est, 'Coarser Grid'),
                       visualize_log_intensity(X_k_both_truth[1:5,,1], time_grid_both, 'Combined'),
                       nrow = 1
                       )  
  
  # 1.2) does a kronecker extractor work? - yes, if we condition on mxm matrix having diagonal of 1's 
  
  wow <- dataset$true_graphs[[1]]$P_block_kronecker
  
  pm_mat <- wow$GP_simu_var
  m_mat <- wow$base_cov
  p_mat <- wow$cor_mat
  pm_mat2 <- kronecker(p_mat, m_mat)
  
  table(pm_mat2 == pm_mat) # sanity check
  
  
  pm_decomp <- kronecker_decomp(pm_mat, p, m)         # decomp with no constraints
  pm_decomp2 <- kronecker_decomp_diag1(pm_mat, p, m)  # decomp with mxm matrix having 1's on diag
  pm_decomp3 <- kronecker_psd_factor(pm_mat, p, m, enforce_trace = FALSE)  # decomp with mxm matrix being PSD
  
  # original decomp, better decomp, and truth
  arr_mat <- matrix(1:6, nrow = 2, byrow = F)
  g_12 <- grid.arrange(visualize_matrix_heatmap(pm_decomp$p_mat, 'Default Decomp P'),
                       visualize_matrix_heatmap(pm_decomp$m_mat, 'Default Decomp M'),
                       visualize_matrix_heatmap(pm_decomp2$p_mat, 'Second Decomp P'),
                       visualize_matrix_heatmap(pm_decomp2$m_mat, 'Second Decomp M'),
                       visualize_matrix_heatmap(p_mat, 'True P'),
                       visualize_matrix_heatmap(m_mat, 'True M'), 
                       layout_matrix = arr_mat) 
  
  
  # more sanity checks
  pm_decomp_reconstruct <- kronecker(pm_decomp$p_mat, pm_decomp$m_mat)
  
  summary(as.numeric(pm_decomp_reconstruct - pm_mat))
  summary(as.numeric(pm_decomp2$p_mat - p_mat))
  
  # 2) visualize rho from X_truth ----------------------------------------------
  
  # 2.1) rho_i_est (mean intensity) from X_truths
  mat_list <- lapply(dataset$subject_data, function(x) exp(x$X_functions))
  rho_simu <- Reduce("+", mat_list) / length(mat_list)
  
  mat_list_coarse <- lapply(dataset$subject_data, function(x) exp(x$X_functions_coarse))
  rho_simu_coarse <- Reduce("+", mat_list_coarse) / length(mat_list_coarse)  
  
  
  
  # 2.2) rho_i_truth (ground truth of GP mean) theoretical truth
  
  kernel_params = dataset$true_graphs[[1]]$P_block_kronecker
  
  rho_i_mean <- exp(kernel_params$base_GP_mean + 0.5 * kernel_params$base_variance)
  GP_means_m <- kernel_params$GP_simu_mean  # m-dim vec
  GP_means_pm <- rep(GP_means_m, p)
  
  GP_kernel_pm <- kernel_params$GP_simu_var # pm x pm matrix
  
  GP_means_exp_pm <- exp(GP_means_pm + 0.5 * diag(GP_kernel_pm)) # ground truth of GP mean
  GP_mean_ground_truth <- matrix(GP_means_exp_pm, nrow = p, ncol = m) # put it in a matrix
  
  g_GP_baseline <- visualize_log_intensity(GP_mean_ground_truth[1:5,], time_grid, 'Rho i Truth from Theory') # graph
  
  # g_data_2_one_subject <- visualize_log_intensity(exp(dataset$X_k_truth[1:5,,1]), time_grid, 'Rho i Truth from X_Truth') + geom_hline(yintercept = rho_i_mean)
  
  g_data_2 <- visualize_log_intensity(rho_simu[1:5,], time_grid, 'Rho i Truth from X_Truth') + geom_hline(yintercept = rho_i_mean) # average intensities across all subjects
  
  g_22 <- grid.arrange(g_data_2, g_GP_baseline, nrow = 1) # compare graphs
  
  
  
  # 3) pre-processing before estimation ------------------------------------------
  
  simu_threshold <- dataset$simulation_params$threshold_p
  
  # 3.1) first, convert dataset into a format that can be used for estimation
  
  df_estimate <- convert_data_for_estimation(dataset$subject_data) %>% as.data.table()
  df_estimate$time <- df_estimate$time / T_max # normalize to [0, 1]
  
  query_yd <- rep(1, n) %>% as.data.frame()  # assume they are all from the same discrete strata
  y_c_strata <- dataset$Y_continuous
  query_y_cs <- dataset$Y_continuous %>% unique()
  discrete_strata <- query_yd %>% unique()
  Tseq_est <- time_grid_est
  
  

  # 3.2) more prep 
  
  estimated_graphs <- list()
  m_est <- length(Tseq_est)
  subject_nums <- 1:n
  s_yd <- 1:n
  data_df4 <- df_estimate 
  
  
  # size of dataset
  print(paste0('number of subjects: ', length(s_yd)))
  print(paste0('number of spikes: ', nrow(data_df4)))
  
  
  # 3.3) retroactively get parameters for subjects (ntrain, patient_sel) and neurons (p, feature_sel)
  
  ntrain <- length(unique(data_df4$subject_num))
  p <- length(unique(data_df4$feature_id))
  
  
  patient_sel = unique(data_df4$subject_num) %>% sort()
  feature_sel = unique(data_df4$feature_id) %>% sort()
  
  
  # 4) estimation --------------------------------------------------------------
  
  
  # part 1 - get X_k_est (p x m_est x n)
  
  result_list <- data_df4 %>%
    group_by(feature_id, subject_num) %>%
    summarise(
      vec = list(estimate_log_intensity_function(time, Tseq_est)),  
      .groups = "drop"
    )
  
  
  X_k_est <- array(
    data = unlist(result_list$vec),
    dim = c(m_est, n, p)
  ) %>% aperm(c(3, 1, 2))
  

  
  # part 2 - rho_i_est from X_k_hat
  rho_list <- estimate_intensities_stratum_parallel_v3(data_df4, patient_sel, feature_sel, Tseq_est, ncores)
  recovered_intensities <- do.call(rbind, lapply(rho_list, function(x) x$rho_i)) # p x m matrix of estimated mean intensities
  
  
  # rho_i_truth = exp(mu(t) + 0.5 * diag(GP_cov))
  rho_i_truth_value <- exp(kernel_params$base_GP_mean + 0.5 * kernel_params$base_variance)
  rho_i_truth <- replicate(p, exp(kernel_params$GP_simu_mean + 0.5 * diag(kernel_params$base_cov))) %>% t()
  rho_i_coarse_truth <- replicate(p, exp(kernel_params$GP_simu_mean_est + 0.5 * diag(kernel_params$base_cov_est))) %>% t()
  
  # rho_ii_truth = rho_i_true * t(rho_i_true) * exp(GP_cov)
  rho_ii_truth <- array(
    unlist(apply(rho_i_truth, 1, function(v) tcrossprod(v))),
    dim = c(m, m, p)
  )
  rho_ii_truth <- sweep(rho_ii_truth, 1:2, exp(kernel_params$base_cov), `*`) # each outer product --> hadamart with exp(K)
  
  rho_ii_coarse_truth <- array(
    unlist(apply(rho_i_coarse_truth, 1, function(v) tcrossprod(v))),
    dim = c(m_est, m_est, p)
  )
  rho_ii_coarse_truth <- sweep(rho_ii_coarse_truth, 1:2, exp(kernel_params$base_cov_est), `*`) # each outer product --> hadamart with exp(K)

  
  # calculate ground truth rho_ii from ground truth X_k's 

  Lambda_k_truth <- exp(X_k_truth)
  Lambda_k_coarse_truth <- exp(X_k_coarse_truth)
  
  rho_ii_truth_from_X <- array(0, dim = c(m, m, p))
  rho_ii_coarse_truth_from_X <- array(0, dim = c(m_est, m_est, p))
  
  for (i in 1:p) {
    mats <- array(0, dim = c(m, m, n))
    coarse_mats <- array(0, dim = c(m_est, m_est, n))
    for (k in 1:n) {
      v <- Lambda_k_truth[i, , k]
      mats[, , k] <- tcrossprod(v)  # faster than outer(v, v)
      v <- Lambda_k_coarse_truth[i, , k]
      coarse_mats[, , k] <- tcrossprod(v)  # faster than outer(v, v)      
    }
    
    # Average across replicates
    rho_ii_truth_from_X[ , , i] <- apply(mats, c(1, 2), mean)
    rho_ii_coarse_truth_from_X[ , , i] <- apply(coarse_mats, c(1, 2), mean)
  }  
  

  
  rho_ii_est <- simplify2array(lapply(rho_list, function(x) x$rho_ii_mat)) # m x m x p matrix of estimated mean intensities
  
  
  # how far do they deviate?
  summary(as.numeric(recovered_intensities - rho_i_coarse_truth)) # need more subjects, still a lot of variability (-51, 77)
                                                                  # 1000 n, 15 mins,    variability = (-20, 23)
  
  
  compute_operator_error_v2(rho_ii_truth[,,1], rho_ii_est[,,1], time_grid, Tseq_est)
  compute_operator_error_v2(rho_ii_est[,,1], rho_ii_est[,,2], Tseq_est, Tseq_est)
  
  
  # visualization of part 2
  
  g_est_21 <- grid.arrange(visualize_log_intensity(recovered_intensities[1:5,], Tseq_est, 'Rho i est from X_hat') + geom_hline(yintercept = rho_i_truth), # rho_i_truth vs rho_i_est
                           g_data_2, 
                           g_GP_baseline, nrow = 1)
  
  g_est_22 <- grid.arrange(visualize_matrix_heatmap(rho_ii_truth[,,1], 'Truth Theory', 50000, 170000),
                           visualize_matrix_heatmap(rho_ii_coarse_truth[,,1], 'Coarse Truth Theory', 50000, 170000),
                           visualize_matrix_heatmap(rho_ii_truth_from_X[,,1], 'Truth X 1', 50000, 170000),
                           visualize_matrix_heatmap(rho_ii_est[,,1], 'Estimate 1', 50000, 170000),
                           visualize_matrix_heatmap(rho_ii_truth_from_X[,,2], 'Truth X 2', 50000, 170000),
                           visualize_matrix_heatmap(rho_ii_est[,,2], 'Estimate 2', 50000, 170000),
                           layout_matrix = arr_mat)
  

  
  # part 3 - GP covariance estimation ------------------------------------------
  
  rho_list_X_truth <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_simu[i, ],     # 50-dim vector
      rho_ii_mat = rho_ii_truth_from_X[, , i]    # 50x50 matrix
    )
  })
  
  rho_list_X_coarse_truth <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_simu_coarse[i, ],     # 50-dim vector
      rho_ii_mat = rho_ii_coarse_truth_from_X[, , i]    # 50x50 matrix
    )
  }) 
  
  rho_list_coarse_truth <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_i_coarse_truth[i, ],     # 50-dim vector
      rho_ii_mat = rho_ii_coarse_truth[, , i]    # 50x50 matrix
    )
  })  
  
  rho_list_truth <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_i_truth[i, ],     # 50-dim vector
      rho_ii_mat = rho_ii_truth[, , i]    # 50x50 matrix
    )
  })     
  
  g_ij_st_est            <- estimate_covariance_functions_ii(rho_list)
  g_ij_st_X_coarse_truth <- estimate_covariance_functions_ii(rho_list_X_coarse_truth)
  g_ij_st_X_truth        <- estimate_covariance_functions_ii(rho_list_X_truth)
  g_ij_st_coarse_truth   <- estimate_covariance_functions_ii(rho_list_coarse_truth)
  g_ij_st_truth          <- estimate_covariance_functions_ii(rho_list_truth)
  
  g_ij_st_truth2 <- kernel_params$base_cov # ground truth should be the GP covariance matrix 
  g_ij_st_coarse_truth2 <- generate_covariance_matrix(Tseq_est, 
                                                     kernel = base_kernel_params$base_kernel, 
                                                     gamma = base_kernel_params$base_gamma, 
                                                     variance = base_kernel_params$base_variance) %>% psd_jitter()
  
  summary(as.numeric(g_ij_st_truth2 - g_ij_st_truth[,,1]))
  summary(as.numeric(g_ij_st_coarse_truth2 - g_ij_st_coarse_truth[,,1]))


  # ground truth covariance matrix
  g_est_31 <- grid.arrange(visualize_matrix_heatmap(g_ij_st_truth[,,1], 'Truth Theory', -1, 1),  
                           visualize_matrix_heatmap(g_ij_st_coarse_truth[,,1], 'Coarse Truth Theory', -1, 1),
                           visualize_matrix_heatmap(g_ij_st_X_truth[,,1], 'Truth X', -1, 1),    
                           visualize_matrix_heatmap(g_ij_st_X_coarse_truth[,,2], 'Coarse Truth X', -1, 1),
                           textGrob("3. Covariance Function\nEstimation", gp = gpar(fontsize = 14)),
                           visualize_matrix_heatmap(g_ij_st_est[,,1], 'Estimate', -1, 1),
                           layout_matrix = arr_mat)
  

  # distance between covariance matrices
  
  
  
  compute_operator_error_v2(g_ij_st_truth[,,1], g_ij_st_est[,,1], time_grid, Tseq_est)
  compute_operator_error_v2(g_ij_st_coarse_truth[,,1], g_ij_st_est[,,1], Tseq_est, Tseq_est)
  compute_operator_error_v2(g_ij_st_est[,,1], g_ij_st_est[,,2], Tseq_est, Tseq_est)

  
  
  # part 4 - eigendecomposition of GP covariance -------------------------------
  
  # perform eigendecomposition
  eigen_decomp_est            <- compute_eigendecomposition_ii(g_ij_st_est)
  eigen_decomp_X_coarse_truth <- compute_eigendecomposition_ii(g_ij_st_X_coarse_truth)
  eigen_decomp_X_truth        <- compute_eigendecomposition_ii(g_ij_st_X_truth)  
  eigen_decomp_coarse_truth   <- compute_eigendecomposition_ii(g_ij_st_coarse_truth)
  eigen_decomp_truth          <- compute_eigendecomposition_ii(g_ij_st_truth)
  
  
  # validate
  g_ij_st_est_decomp             <- validate_eigendecomposition_ii(g_ij_st_est, eigen_decomp_est)
  g_ij_st_X_coarse_truth_decomp  <- validate_eigendecomposition_ii(g_ij_st_X_coarse_truth, eigen_decomp_X_coarse_truth)  
  g_ij_st_X_truth_decomp         <- validate_eigendecomposition_ii(g_ij_st_X_truth, eigen_decomp_X_truth)  
  g_ij_st_coarse_truth_decomp    <- validate_eigendecomposition_ii(g_ij_st_coarse_truth, eigen_decomp_coarse_truth)  
  g_ij_st_truth_decomp           <- validate_eigendecomposition_ii(g_ij_st_truth, eigen_decomp_truth)

  
  
  # 4.1) do the eigenfunctions between coarse and regular truths differ? - yes - the g_ij_st is too different from the truth
  
  g_est_41 <- grid.arrange(visualize_log_intensity(t(eigen_decomp_truth$eigenfunctions[[1]]),           time_grid, 'Truth Theory'),
                           visualize_log_intensity(t(eigen_decomp_coarse_truth$eigenfunctions[[1]]),    Tseq_est,  'Coarse Truth Theory'),
                           visualize_log_intensity(t(eigen_decomp_X_truth$eigenfunctions[[1]]),         time_grid, 'Truth X'),
                           visualize_log_intensity(t(eigen_decomp_X_coarse_truth$eigenfunctions[[1]]),  Tseq_est,  'Coarse Truth X'),
                           textGrob("4. Eigenfunction\nBasis of\nProcess 1", gp = gpar(fontsize = 14)),
                           visualize_log_intensity(t(eigen_decomp_est$eigenfunctions[[1]]),             Tseq_est,  'Estimate'),
                           layout_matrix = arr_mat)
  
  # 4.2) are all eigenfunctions orthogonal? - yes
  
  t(eigen_decomp_est$eigenfunctions[[1]]) %*% eigen_decomp_est$eigenfunctions[[1]]
  t(eigen_decomp_truth$eigenfunctions[[1]]) %*% eigen_decomp_truth$eigenfunctions[[1]]
  t(eigen_decomp_coarse_truth$eigenfunctions[[1]]) %*% eigen_decomp_coarse_truth$eigenfunctions[[1]]
  
  # did eigendecomposition reconstruct the g matrix to the best of its ability? yes
  
  summary(as.numeric(g_ij_st_est[,,1] - g_ij_st_est_decomp[,,1])) # good (-3e-3, 3e-3)
  summary(as.numeric(g_ij_st_truth[,,1] - g_ij_st_truth_decomp[,,1])) # good (-3e-3, 5e-3)
  summary(as.numeric(g_ij_st_coarse_truth[,,1] - g_ij_st_coarse_truth_decomp[,,1])) # very good (-7e-4, 9e-4)
  
  
  # visualize g_ij_st vs its eigendecomposition reconstruction


  g_est_42 <- grid.arrange(visualize_matrix_heatmap(g_ij_st_truth[,,1], 'Truth Theory (TT)', -1, 1),  # ground truth 50 x 50 covariance
                           visualize_matrix_heatmap(g_ij_st_truth_decomp[,,1], 'TT Reconstruct', -1, 1), # ground truth reconstructed
                           visualize_matrix_heatmap(g_ij_st_coarse_truth[,,1], 'Coarse Truth Theory (CTT)', -1, 1), # ground truth 19 x 19 covariance
                           visualize_matrix_heatmap(g_ij_st_coarse_truth_decomp[,,1], 'CTT Reconstruct', -1, 1), # ground truth 19 x 19 covariance
                           visualize_matrix_heatmap(g_ij_st_est[,,1], 'Estimate (E)', -1, 1),    # estimate 19 x 19 covariance
                           visualize_matrix_heatmap(g_ij_st_est_decomp[,,1], 'E Reconstruct', -1, 1), # reconstructed estimate
                           layout_matrix = arr_mat)    
  
  # part 5 - KL expansion ------------------------------------------------------
  
  # 5.1) center X_k_est and X_k_truth

  
  mean_mat_est <- apply(X_k_est, c(1, 2), mean)  # result is p x m matrix
  X_k_est_center <- sweep(X_k_est, c(1, 2), mean_mat_est, FUN = "-")
  
  
  mean_mat_truth <- apply(X_k_truth, c(1, 2), mean)  # result is p x m matrix
  X_k_truth_center_est <- sweep(X_k_truth, c(1, 2), mean_mat_truth, FUN = "-") 
  X_k_truth_center <- X_k_truth - base_GP_mean
  
  
  mean_mat_coarse_truth <- apply(X_k_coarse_truth, c(1, 2), mean)  # result is p x m matrix
  X_k_coarse_truth_center_est <- sweep(X_k_coarse_truth, c(1, 2), mean_mat_coarse_truth, FUN = "-")   
  X_k_coarse_truth_center <- X_k_coarse_truth - base_GP_mean
  
  
   
  # 5.2) get KL coeffs   (n x p x d)
  

  
  kl_coeffs_v2                <- estimate_kl_coefficients_parallel_v2(X_k_est_center,               eigen_decomp_est$eigenfunctions,            Tseq_est,  ncores)
  kl_coeffs_X_coarse_truth_v2 <- estimate_kl_coefficients_parallel_v2(X_k_coarse_truth_center_est,  eigen_decomp_X_coarse_truth$eigenfunctions, Tseq_est,  ncores)
  kl_coeffs_X_truth_v2        <- estimate_kl_coefficients_parallel_v2(X_k_truth_center_est,         eigen_decomp_X_truth$eigenfunctions,        time_grid, ncores)  
  kl_coeffs_coarse_truth_v2   <- estimate_kl_coefficients_parallel_v2(X_k_coarse_truth_center_est,  eigen_decomp_coarse_truth$eigenfunctions,   Tseq_est,  ncores)
  kl_coeffs_truth_v2          <- estimate_kl_coefficients_parallel_v2(X_k_truth_center_est,         eigen_decomp_truth$eigenfunctions,          time_grid, ncores)
    
  
  
  # 5.3) Graph X_k (true, true_reconstruct, estimate, estimate_reconstruct)  for the first subject
  
  df1 <- validate_kl_full(X_k_est,          eigen_decomp_est$eigenfunctions,            Tseq_est,  'Estimate', ncores)
  df2 <- validate_kl_full(X_k_coarse_truth, eigen_decomp_X_coarse_truth$eigenfunctions, Tseq_est,  'Coarse Truth X', ncores)
  df3 <- validate_kl_full(X_k_truth,        eigen_decomp_X_truth$eigenfunctions,        time_grid, 'Truth X', ncores)
  df4 <- validate_kl_full(X_k_coarse_truth, eigen_decomp_coarse_truth$eigenfunctions,   Tseq_est,  'Coarse Truth Theory', ncores)
  df5 <- validate_kl_full(X_k_truth,        eigen_decomp_truth$eigenfunctions,          time_grid, 'Truth Theory', ncores)
  
  df_both <- rbind(df1, df5)
  df_both$cat <- factor(df_both$cat)
  
  df_all <- rbind(df1, df2, df3, df4, df5)
  df_all$cat <- factor(df_all$cat) 
  
  df_no_reconstruct <- df_all %>% filter(cat %in% c('Estimate', 'Coarse Truth X', 'Truth X', 'Coarse Truth Theory', 'Truth Theory'))
  
  # compare log intensities for subject 1, KL coefficients do a fine job at approximating
  # but the estimated intensity is just awful 
  
  g_est_51 <- ggplot() + 
    geom_line(data = df_both, aes(x = Var2, y = value, color = cat)) + 
    facet_wrap(~ Var1) + 
    ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
    ylab('Log Intensity') + 
    xlab('Time') + 
    labs(color = "Truth or Estimate") + 
    theme_bw()
    
  g_est_52 <- ggplot() + 
    geom_line(data = df_all, aes(x = Var2, y = value, color = cat)) + 
    facet_wrap(~ Var1) + 
    ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
    ylab('Log Intensity') + 
    xlab('Time') + 
    labs(color = "Truth or Estimate") + 
    theme_bw()   
  
  g_est_53 <- ggplot() + 
    geom_line(data = df_no_reconstruct, aes(x = Var2, y = value, color = cat)) + 
    facet_wrap(~ Var1) + 
    ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
    ylab('Log Intensity') + 
    xlab('Time') + 
    labs(color = "Truth or Estimate") + 
    theme_bw()  
  
  # here, we fit to y_c_strata -----------------------------------------------
  
  
  # part 6
  
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  K_c <- construct_kernel_matrix_step_6(y_c_strata, gamma_c) # (n_stratum x n_stratum)
  
  table(K_c) # adjacent weeks have a correlation of 0.3678
  
  
  # part 7 ---------------------------------------------------------------------

  V_YcXij_est             <- construct_cross_covariance_matrix_v3(kl_coeffs_v2, ncores)
  V_YcXij_X_coarse_truth  <- construct_cross_covariance_matrix_v3(kl_coeffs_X_coarse_truth_v2, ncores)
  V_YcXij_X_truth         <- construct_cross_covariance_matrix_v3(kl_coeffs_X_truth_v2, ncores)  
  V_YcXij_coarse_truth    <- construct_cross_covariance_matrix_v3(kl_coeffs_coarse_truth_v2, ncores)
  V_YcXij_truth           <- construct_cross_covariance_matrix_v3(kl_coeffs_truth_v2, ncores)
  
  
  M_hat_est            <- estimate_regression_operators_v3(K_c, V_YcXij_est, p)
  M_hat_X_coarse_truth <- estimate_regression_operators_v3(K_c, V_YcXij_X_coarse_truth, p)
  M_hat_X_truth        <- estimate_regression_operators_v3(K_c, V_YcXij_X_truth, p)  
  M_hat_coarse_truth   <- estimate_regression_operators_v3(K_c, V_YcXij_coarse_truth, p)
  M_hat_truth          <- estimate_regression_operators_v3(K_c, V_YcXij_truth, p)

  
  # 7.1) do these regression coefficients match? - idk, truth matches coarse truth
  
  g_est_71 <- grid.arrange(visualize_matrix_heatmap(V_YcXij_truth[['1_2']][1,,], 'Truth Theory'),
                           visualize_matrix_heatmap(V_YcXij_coarse_truth[['1_2']][1,,], 'Coarse Truth Theory'),
                           visualize_matrix_heatmap(V_YcXij_X_truth[['1_2']][1,,], 'Truth X'),               
                           visualize_matrix_heatmap(V_YcXij_X_coarse_truth[['1_2']][1,,], 'Coarse Truth X'),
                           textGrob("7. Cross Covariance\n of KL Coeffs", gp = gpar(fontsize = 14)),
                           visualize_matrix_heatmap(V_YcXij_est[['1_2']][1,,], 'Estimate'),
                           layout_matrix = arr_mat)
  
  g_est_72 <- grid.arrange(visualize_matrix_heatmap(M_hat_truth[['1_2']][1,,], 'Truth Theory'),
                           visualize_matrix_heatmap(M_hat_coarse_truth[['1_2']][1,,], 'Coarse Truth Theory'),
                           visualize_matrix_heatmap(M_hat_X_truth[['1_2']][1,,], 'Truth X'),
                           visualize_matrix_heatmap(M_hat_X_coarse_truth[['1_2']][1,,], 'Coarse Truth X'),
                           textGrob("7. Regression Operator\n of KL Coeffs", gp = gpar(fontsize = 14)),
                           visualize_matrix_heatmap(M_hat_est[['1_2']][1,,], 'Estimate'),
                           layout_matrix = arr_mat)  
  
  # part 8 - regress on y_c ----------------------------------------------------
  
  print('at covariate loop')
  
  cont_inds <- 1:nrow(query_y_cs)
  estimated_graphs_v2 <- pbmclapply(cont_inds, function(cont_ind) {
    
    
    query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()
    
    V_cond_est            <- evaluate_regression_at_query_v2(M_hat_est, y_c_strata, query_y_c, eigen_decomp_est$eigenfunctions, gamma_c, p)
    V_cond_X_coarse_truth <- evaluate_regression_at_query_v2(M_hat_X_coarse_truth, y_c_strata, query_y_c, eigen_decomp_X_coarse_truth$eigenfunctions, gamma_c, p)
    V_cond_X_truth        <- evaluate_regression_at_query_v2(M_hat_X_truth, y_c_strata, query_y_c, eigen_decomp_X_truth$eigenfunctions, gamma_c, p)    
    V_cond_coarse_truth   <- evaluate_regression_at_query_v2(M_hat_coarse_truth, y_c_strata, query_y_c, eigen_decomp_coarse_truth$eigenfunctions, gamma_c, p)
    V_cond_truth          <- evaluate_regression_at_query_v2(M_hat_truth, y_c_strata, query_y_c, eigen_decomp_truth$eigenfunctions, gamma_c, p)

    
    
    # use KL coefficients to obtain V_cond
    # V_cond_truth_12 <- validate_V_ij_v2(kl_coeffs_truth_v2, eigen_decomp_truth$eigenfunctions, 1, 2, y_c_strata, query_y_c, gamma_c)
    # V_cond_truth_11 <- validate_V_ij_v2(kl_coeffs_truth_v2, eigen_decomp_truth$eigenfunctions, 1, 1, y_c_strata, query_y_c, gamma_c)
    
    
    # 8.2) visualizations 
    
    # for pair 1_2, plot estimate and truths - this one should look correlated
    g_est_81 <- grid.arrange(visualize_matrix_heatmap(V_cond_truth[['1_2']], 'Truth Theory'),
                             visualize_matrix_heatmap(V_cond_coarse_truth[['1_2']], 'Coarse Truth Theory'),
                             visualize_matrix_heatmap(V_cond_X_truth[['1_2']], 'Truth X'), 
                             visualize_matrix_heatmap(V_cond_X_coarse_truth[['1_2']], 'Coarse Truth X'), 
                             textGrob("8. Conditional\nCovariance Operator\nProcess 1 with 2", gp = gpar(fontsize = 14)),
                             visualize_matrix_heatmap(V_cond_est[['1_2']], 'Estimate'), 
                             layout_matrix = arr_mat)    
    
    # but these two versions of the truth aren't exactly the same
    # grid.arrange(visualize_pm_block_matrix_heatmap(V_cond_truth[['1_2']]), 
    #              visualize_pm_block_matrix_heatmap(V_cond_truth_12), 
    #              layout_matrix = arr_mat)
    

    # visualize again for process 3 with 6 - they should not look correlated
    g_est_82 <- grid.arrange(visualize_matrix_heatmap(V_cond_truth[['3_6']], 'Truth Theory'),
                             visualize_matrix_heatmap(V_cond_coarse_truth[['3_6']], 'Coarse Truth Theory'),
                             visualize_matrix_heatmap(V_cond_X_truth[['3_6']], 'Truth X'), 
                             visualize_matrix_heatmap(V_cond_X_coarse_truth[['3_6']], 'Coarse Truth X'), 
                             textGrob("8. Conditional\nCovariance Operator\nProcess 3 with 6", gp = gpar(fontsize = 14)),
                             visualize_matrix_heatmap(V_cond_est[['3_6']], 'Estimate'), 
                             layout_matrix = arr_mat) 
    

    
    # 8.3) visualize the entire pm x pm matrix
    
    V_cond_est_full            <- assemble_block_matrix_v2(V_cond_est, p, m_est)
    V_cond_X_coarse_truth_full <- assemble_block_matrix_v2(V_cond_X_coarse_truth, p, m_est)
    V_cond_X_truth_full        <- assemble_block_matrix_v2(V_cond_X_truth, p, m)    
    V_cond_coarse_truth_full   <- assemble_block_matrix_v2(V_cond_coarse_truth, p, m_est)
    V_cond_truth_full          <- assemble_block_matrix_v2(V_cond_truth, p, m)
    
    

    
    # est vs truth vs ground truth
    g_est_83 <- grid.arrange(visualize_pm_block_matrix_heatmap(V_cond_est_full, 'Estimate'),
                             visualize_pm_block_matrix_heatmap(V_cond_truth_full, 'Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(kernel_params$GP_simu_var, 'Ground Truth'), 
                             nrow = 1)
    
    # all 5 + ground truth
    g_est_84 <- grid.arrange(visualize_pm_block_matrix_heatmap(kernel_params$GP_simu_var, 'Ground Truth'), 
                             visualize_pm_block_matrix_heatmap(V_cond_est_full, 'Estimate'),
                             visualize_pm_block_matrix_heatmap(V_cond_truth_full, 'Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(V_cond_coarse_truth_full, 'Coarse Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(V_cond_X_truth_full, 'Truth X'), 
                             visualize_pm_block_matrix_heatmap(V_cond_X_coarse_truth_full, 'Coarse Truth X'), 
                             layout_matrix = arr_mat)    
    
        
    
    # part 9 - correlation operator --------------------------------------------
    
    C_cond_est            <- estimate_conditional_correlation_v4(V_cond_est, p)
    C_cond_X_coarse_truth <- estimate_conditional_correlation_v4(V_cond_X_coarse_truth, p)
    C_cond_X_truth        <- estimate_conditional_correlation_v4(V_cond_X_truth, p)    
    C_cond_coarse_truth   <- estimate_conditional_correlation_v4(V_cond_coarse_truth, p)
    C_cond_truth          <- estimate_conditional_correlation_v4(V_cond_truth, p)
    
    
    # 9.1) visualize i_j matrix
    
    g_est_91 <- grid.arrange(visualize_pm_block_matrix_heatmap(C_cond_truth[['1_2']], 'Truth Theory'),
                             visualize_pm_block_matrix_heatmap(C_cond_coarse_truth[['1_2']], 'Coarse Truth Theory'),
                             visualize_pm_block_matrix_heatmap(C_cond_X_truth[['1_2']], 'Truth X'),
                             visualize_pm_block_matrix_heatmap(C_cond_X_coarse_truth[['1_2']], 'Coarse Truth X'),
                             textGrob("9. Conditional\nCorrelation Operator\nProcess 1 with 2", gp = gpar(fontsize = 14)),
                             visualize_pm_block_matrix_heatmap(C_cond_est[['1_2']], 'Estimate'),
                 layout_matrix = arr_mat)  
    
    g_est_92 <- grid.arrange(visualize_pm_block_matrix_heatmap(C_cond_truth[['3_6']], 'Truth Theory'),
                             visualize_pm_block_matrix_heatmap(C_cond_coarse_truth[['3_6']], 'Coarse Truth Theory'),
                             visualize_pm_block_matrix_heatmap(C_cond_X_truth[['3_6']], 'Truth X'),
                             visualize_pm_block_matrix_heatmap(C_cond_X_coarse_truth[['3_6']], 'Coarse Truth X'),
                             textGrob("9. Conditional\nCorrelation Operator\nProcess 3 with 6", gp = gpar(fontsize = 14)),
                             visualize_pm_block_matrix_heatmap(C_cond_est[['3_6']], 'Estimate'),
                             layout_matrix = arr_mat)      
    
    # 9.2) visualize the entire pm x pm block 
    
    C_cond_est_full            <- assemble_block_matrix_v2(C_cond_est, p, m_est)
    C_cond_X_coarse_truth_full <- assemble_block_matrix_v2(C_cond_X_coarse_truth, p, m_est)
    C_cond_X_truth_full        <- assemble_block_matrix_v2(C_cond_X_truth, p, m)    
    C_cond_coarse_truth_full   <- assemble_block_matrix_v2(C_cond_coarse_truth, p, m_est)
    C_cond_truth_full          <- assemble_block_matrix_v2(C_cond_truth, p, m)
    
    
    g_est_93 <- grid.arrange(visualize_pm_block_matrix_heatmap(C_cond_est_full, 'Estimate'),
                             visualize_pm_block_matrix_heatmap(C_cond_truth_full, 'Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(kernel_params$GP_simu_var, 'Ground Truth'), 
                             nrow = 1)    
    
    
    # all 5 + ground truth
    g_est_94 <- grid.arrange(visualize_pm_block_matrix_heatmap(kernel_params$GP_simu_var, 'Ground Truth'), 
                             visualize_pm_block_matrix_heatmap(C_cond_est_full, 'Estimate'),
                             visualize_pm_block_matrix_heatmap(C_cond_truth_full, 'Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(C_cond_coarse_truth_full, 'Coarse Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(C_cond_X_truth_full, 'Truth X'), 
                             visualize_pm_block_matrix_heatmap(C_cond_X_coarse_truth_full, 'Coarse Truth X'), 
                             layout_matrix = arr_mat)       
    
    
    # part 10 - precision operator ---------------------------------------------
    
    P_cond_est            <- estimate_precision_operator_v3(C_cond_est, p)
    P_cond_X_coarse_truth <- estimate_precision_operator_v3(C_cond_X_coarse_truth, p)
    P_cond_X_truth        <- estimate_precision_operator_v3(C_cond_X_truth, p)    
    P_cond_coarse_truth   <- estimate_precision_operator_v3(C_cond_coarse_truth, p)
    P_cond_truth          <- estimate_precision_operator_v3(C_cond_truth, p)
    
    prec_truth <- kronecker(true_graphs[[1]]$P_block_kronecker$prec_mat, true_graphs[[1]]$P_block_kronecker$base_precision)
    prec_coarse_truth <- kronecker(true_graphs[[1]]$P_block_kronecker$prec_mat, true_graphs[[1]]$P_block_kronecker$base_precision_est)
    
    
    P_ground_truth_12 <- extract_block_structure_ij(prec_truth, p, m, 1, 2)
    P_ground_truth_36 <- extract_block_structure_ij(prec_truth, p, m, 3, 6)
    
    # 10.1) visualize
    
    g_est_101 <- grid.arrange(visualize_pm_block_matrix_heatmap(P_ground_truth_12, 'Truth Ground'),
                              visualize_pm_block_matrix_heatmap(P_cond_est[['1_2']], 'Estimate'),
                              visualize_pm_block_matrix_heatmap(P_cond_truth[['1_2']], 'Truth Theory'),
                              visualize_pm_block_matrix_heatmap(P_cond_coarse_truth[['1_2']], 'Coarse Truth Theory'),
                              visualize_pm_block_matrix_heatmap(P_cond_X_truth[['1_2']], 'Truth X'),
                              visualize_pm_block_matrix_heatmap(P_cond_X_coarse_truth[['1_2']], 'Coarse Truth X'),
                              layout_matrix = arr_mat)  
    
    g_est_102 <- grid.arrange(visualize_pm_block_matrix_heatmap(P_ground_truth_36, 'Truth Ground'),
                              visualize_pm_block_matrix_heatmap(P_cond_est[['3_6']], 'Estimate'),
                              visualize_pm_block_matrix_heatmap(P_cond_truth[['3_6']], 'Truth Theory'),
                              visualize_pm_block_matrix_heatmap(P_cond_coarse_truth[['3_6']], 'Coarse Truth Theory'),
                              visualize_pm_block_matrix_heatmap(P_cond_X_truth[['3_6']], 'Truth X'),
                              visualize_pm_block_matrix_heatmap(P_cond_X_coarse_truth[['3_6']], 'Coarse Truth X'),
                              layout_matrix = arr_mat)     
    
    
    # 10.2) visualize the entire pm x pm block
    
    P_cond_est_full            <- assemble_block_matrix_v2(P_cond_est, p, m_est)
    P_cond_X_coarse_truth_full <- assemble_block_matrix_v2(P_cond_X_coarse_truth, p, m_est)
    P_cond_X_truth_full        <- assemble_block_matrix_v2(P_cond_X_truth, p, m)    
    P_cond_coarse_truth_full   <- assemble_block_matrix_v2(P_cond_coarse_truth, p, m_est)
    P_cond_truth_full          <- assemble_block_matrix_v2(P_cond_truth, p, m)
    
    
    
    

    
    g_est_103 <- grid.arrange(visualize_pm_block_matrix_heatmap(P_cond_est_full, 'Estimate'),
                              visualize_pm_block_matrix_heatmap(P_cond_truth_full, 'Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(prec_coarse_truth, 'Ground Truth'), 
                              nrow = 1)        
    

    # all 5 + ground truth
    g_est_104 <- grid.arrange(visualize_pm_block_matrix_heatmap(prec_coarse_truth, 'Ground Truth'), 
                              visualize_pm_block_matrix_heatmap(P_cond_est_full, 'Estimate'),
                              visualize_pm_block_matrix_heatmap(P_cond_truth_full, 'Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(P_cond_coarse_truth_full, 'Coarse Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(P_cond_X_truth_full, 'Truth X'), 
                              visualize_pm_block_matrix_heatmap(P_cond_X_coarse_truth_full, 'Coarse Truth X'), 
                              layout_matrix = arr_mat)     
    
    
    # part 11
    #graph_threshold <- select_threshold_by_stability(P_cond, p)
    final_graph_estimates                  <- estimate_graph(P_cond_est,            C_cond_est,            V_cond_est, threshold, p)
    final_graph_estimates_X_coarse_truth   <- estimate_graph(P_cond_X_coarse_truth, C_cond_X_coarse_truth, V_cond_X_coarse_truth, threshold, p)
    final_graph_estimates_X_truth          <- estimate_graph(P_cond_X_truth,        C_cond_X_truth,        V_cond_X_truth, threshold, p)    
    final_graph_estimates_coarse_truth     <- estimate_graph(P_cond_coarse_truth,   C_cond_coarse_truth,   V_cond_coarse_truth, threshold, p)
    final_graph_estimates_truth            <- estimate_graph(P_cond_truth,          C_cond_truth,          V_cond_truth, threshold, p)
    
    g_est_111 <- grid.arrange(visualize_pm_block_matrix_heatmap(true_graphs[[1]]$adj_mat, 'Ground Truth'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates$w_mat, 'Estimate'),
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_coarse_truth$w_mat, 'Coarse Truth X'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_truth$w_mat, 'Truth X'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_coarse_truth$w_mat, 'Coarse Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_truth$w_mat, 'Truth Theory'), 
                              layout_matrix = arr_mat)  
    

    
    # ROC curve 
    
    w_mat_ground_truth        <- hilbert_schmidt_norm_pm(prec_truth, p, m)
    w_mat_coarse_ground_truth <- hilbert_schmidt_norm_pm(prec_coarse_truth, p, m_est)
    w_mat_X_coarse_truth      <- hilbert_schmidt_norm_pm(P_cond_X_coarse_truth_full, p, m_est)
    w_mat_X_truth             <- hilbert_schmidt_norm_pm(P_cond_X_truth_full, p, m)    
    w_mat_coarse_truth        <- hilbert_schmidt_norm_pm(P_cond_coarse_truth_full, p, m_est)
    w_mat_truth               <- hilbert_schmidt_norm_pm(P_cond_truth_full, p, m)
    w_mat_est                 <- hilbert_schmidt_norm_pm(P_cond_est_full, p, m_est)
    
    
    arr_mat2 <- matrix(1:8, nrow = 2, byrow = F)
    g_est_112 <- grid.arrange(visualize_pm_block_matrix_heatmap(w_mat_ground_truth, 'Ground Truth'), 
                              visualize_pm_block_matrix_heatmap(w_mat_coarse_ground_truth, 'Coarse Ground Truth'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_truth$w_mat, 'Truth Theory'),   
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_coarse_truth$w_mat, 'Coarse Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_truth$w_mat, 'Truth X'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_coarse_truth$w_mat, 'Coarse Truth X'), 
                              textGrob("11. Hilbert Schmidt\n Norm", gp = gpar(fontsize = 14)),
                              visualize_pm_block_matrix_heatmap(final_graph_estimates$w_mat, 'Estimate'),
                              layout_matrix = arr_mat2)     

    
    roc_ground_truth        <- roc_with_threshold(w_mat_ground_truth,        true_graphs[[1]]$adj_mat, 'Ground Truth')
    roc_coarse_ground_truth <- roc_with_threshold(w_mat_coarse_ground_truth, true_graphs[[1]]$adj_mat, 'Coarse Ground Truth')
    roc_truth               <- roc_with_threshold(w_mat_truth,               true_graphs[[1]]$adj_mat, 'Truth Theory')   
    roc_coarse_truth        <- roc_with_threshold(w_mat_coarse_truth,        true_graphs[[1]]$adj_mat, 'Coarse Truth Theory')
    roc_X_truth             <- roc_with_threshold(w_mat_X_truth,             true_graphs[[1]]$adj_mat, 'Truth X')
    roc_X_coarse_truth      <- roc_with_threshold(w_mat_X_coarse_truth,      true_graphs[[1]]$adj_mat, 'Coarse Truth X')
    roc_est                 <- roc_with_threshold(w_mat_est,                 true_graphs[[1]]$adj_mat, 'Estimate')
    
    g_est_113 <- grid.arrange(roc_ground_truth$plot,
                              roc_coarse_ground_truth$plot,
                              roc_truth$plot,
                              roc_coarse_truth$plot,
                              roc_X_truth$plot,
                              roc_X_coarse_truth$plot,
                              textGrob("11. ROC Curve", gp = gpar(fontsize = 14)),
                              roc_est$plot,
                              layout_matrix = arr_mat2)
    
  
    
    
    # 8) save graphs (recall that we have ground truths of the form a_b_c)
    
    list(g_81  = g_est_81,  g_82  = g_est_82,  g_83  = g_est_83,  g_84  = g_est_84,
         g_91  = g_est_91,  g_92  = g_est_92,  g_93  = g_est_93,  g_94  = g_est_94,
         g_101 = g_est_101, g_102 = g_est_102, g_103 = g_est_103, g_104 = g_est_104,
         g_111 = g_est_111, g_112 = g_est_112, g_113 = g_est_113)
    
    
  }, mc.cores = ncores) # done with all y_c levels
  
  
  
  names(estimated_graphs_v2) <- apply(query_y_cs, 1, function(x){paste(x, collapse = '_')})

  
  # original graphs
  
  estimated_graphs_part_1 <- list(g_21 = g_est_21, g_22 = g_est_22,
                                  g_31 = g_est_31,
                                  g_41 = g_est_41, g_42 = g_est_42,
                                  g_51 = g_est_51, g_52 = g_est_52, g_53 = g_est_53,
                                  g_71 = g_est_71, g_72 = g_est_72)
  
  
  
  return(list(estimated_graphs_part_1 = estimated_graphs_part_1,      
              estimated_graphs_part_2 = estimated_graphs_v2))
  
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