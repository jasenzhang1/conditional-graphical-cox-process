# truncated from "with truths" for mice data


# DO NOT KEEP TRACK OF TRUTHS THROUGHOUT THE PROCESS

full_conditional_estimation_with_no_truth <- function(dataset, method, ncores, dir){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: regular estimation without keeping track of truths
  #
  # 
  #
  #
  # Supports:
  #   - "OG", "JASA"  → covariate loop after step_5 (v2-style)
  #   - "CPGM"        → covariate loop after step_1 (v3-style)
  #
  # Input:
  #
  # - dataset  (list of the following features)
  #
  #   - event_times         (n*p-dim list)      each item is a vector of timestamps for process i in subject k
  #   - X_k_truth           (p x m x n matrix)
  #   - X_k_coarse_truth    (p x m_est x n matrix)
  #   - X_k_both_truth      (p x m_both x n matrix)
  #   - Y_continuous        (n x q_c matrix)
  #   - simulation_params   (list of params such as n, p, query_y_cs, adj_type, adj_params, time_grids, seed)
  #   - true_graphs         (list of true settings for each query_y_cs)
  #   - summary_stats       (list of summary stats)
  #
  # 
  # - method   ('CPGM')  
  # - ncores   (integer)
  # - dir      (string)     folder name, such as "simu_results/block_banded_v2"
  # 
  #
  # ----------------------------------------------------------------------------
  
  if (!(method %in% c("CPGM"))) {
    stop("Error 12b: estimation method must be 'CPGM'")
  }
  
  # load  
  
  time_grid_est <- dataset$simulation_params$time_grid_est
  time_grid <- dataset$simulation_params$time_grid
  time_grid_both <- dataset$simulation_params$time_grid_both
  n <- dim(dataset$Y_continuous)[1]
  m_est <- length(time_grid_est)
  


  # ----------------------------------------------------------------------------
  # Step 0 - Check and preprocess data
  # ----------------------------------------------------------------------------
  
  step_0_events <- step_0_keep_events(dataset, k = 1, i_vec = 1:5)

  processed_data <- step_0_preprocess(dataset)
  data_df4     <- processed_data[[1]]
  y_c_strata   <- processed_data[[2]]
  query_y_cs   <- processed_data[[3]]
  patient_sel  <- processed_data[[4]]
  feature_sel  <- processed_data[[5]]
  rm(processed_data)
  
  p <- length(feature_sel)
  full <- F
  
  # ----------------------------------------------------------------------------
  # Step 1 - Log intensities
  # ----------------------------------------------------------------------------
  
  step_1 <- step_1_log_intensities(dataset, data_df4, time_grid_est, NA, NA, full)
  
  # start covariate loop 
  # ----------------------------------------------------------------------------
  
  print('at covariate loop')
  cont_inds <- 1:nrow(query_y_cs)
  estimated_graphs <- lapply(cont_inds, function(cont_ind) {
    
    print(paste0(cont_ind, ' out of ', length(cont_inds)))
    
    query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()
    
    if(method %in% c('OG', 'JASA')){
      rho_kernel <- F
      i_neq_j <- F
      query_y_c_step_2 <- 'hold'
      y_c_strata_step_2 <- 'hold'
    } else{
      rho_kernel <- T
      i_neq_j <- T
      query_y_c_step_2 <- query_y_c
      y_c_strata_step_2 <- y_c_strata     
    }
    
    # Step 2: rho_i estimation (per subject)
    
    
    step_2 <- step_2_rho_i(dataset, data_df4, kernel_params_i, rho_kernel,
                           patient_sel, feature_sel, time_grid, time_grid_est,
                           query_y_c_step_2, y_c_strata_step_2, ncores, full)           
    
    
    # Step 2b onward
    step_2b <- step_2_rho_ij(step_1, step_2, kernel_params_i, i_neq_j, full)
    step_3  <- step_3_g_ij(step_2, step_2b, kernel_params_i, i_neq_j, full)
    
    
    norm_G <- T # try normalizing G_ii = G_ii * diag(delta_t)
    norm_vec <- F # normalize eta so that Delta * \eta^\top * \eta = 1
    
    step_4 <- tryCatch({
      step_4_eigendecomp(step_3, p, time_grid, time_grid_est, norm_G, norm_vec, full)
    }, error = function(e) {
      cat("Error in step_4, saving dataset...\n")
      save(dataset, file = file.path(dir, "dataset.RData"))
      stop(e)
    })
    
    # step 5 to 9 split:
    
    if(method %in% c('OG', 'JASA')){
      step_5 <- step_5_KL_expansion(step_1, step_4, kernel_params_i, time_grid, time_grid_est, ncores)
      step_8 <- steps_78(step_4, step_5, kernel_params_i, y_c_strata, query_y_c, method, ncores)
      step_9 <- step_9_C_cond_from_V_cond(step_8, kernel_params_i)
    } else{
      step_5 <- step_5_KL_covariance(step_3, step_4, norm_G, full)
      step_9 <- step_9_C_cond_from_KL_cov(step_4, step_5, kernel_params_i, full)
    }
    
    # reunite at step 10 onwards
    
    block <- F
    MP <- F
    
    # step_10 <- step_10_P_cond(step_9, kernel_params_i, p, block, MP)
    result <- tryCatch({
      step_10 <- step_10_P_cond(step_9, kernel_params_i, p, block, MP, full)
    }, error = function(e) {
      cat("Error occurred in step_10, saving dataset...\n")
      save(dataset, file = paste0(dir, "/dataset.RData"))
      cat("Dataset saved to dataset.RData\n")
      stop(e)  # Re-throw the error
    })
    
    
    step_11 <- step_11_HS_norms(step_10, adj_mat_i, p, full)

    
    list(step_2 = step_2, step_2b = step_2b, step_3 = step_3,
         step_4 = step_4, step_5 = step_5, step_9 = step_9,
         step_10 = step_10, step_11 = step_11)
  })
  
  # ----------------------------------------------------------------------------
  # Optional: reorganize by steps instead of by y_c_query
  # ----------------------------------------------------------------------------
  
  steps <- unique(unlist(lapply(estimated_graphs, names)))
  reorganized <- setNames(lapply(steps, function(step) {
    sapply(estimated_graphs, `[[`, step, simplify = FALSE)
  }), steps)
  
  # Return step_0, step_1 (shared) + reorganized per-subject steps
  estimated_graphs_part_1 <- list(step_0_events = step_0_events,
                                  step_1 = step_1)
  
  all_results <- c(estimated_graphs_part_1, reorganized)
  all_results$y_c_query <- query_y_cs
  all_results$time_grid <- time_grid
  all_results$time_grid_est <- time_grid_est
  all_results$time_grid_both <- time_grid_both
  
  return(all_results)
  
}

full_conditional_estimation_mice <- function(data_df4, patient_sel, feature_sel, Tseq_est, terse, ncores){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: regular estimation without keeping track of truths
  #
  #
  # Input:
  #
  # - data_df4         (data.frame of 3 columns)
  #   - 'feature_id'
  #   - 'time'
  #   - 'subject_num'
  #
  # - patient_sel      (n-dim vector)   patient ID's 
  # - feature_sel      (p-dim vector)   process ID's
  #
  # - terse            (boolean)        if true, return much less
  # - ncores           (integer)        number of cores
  #
  # 
  #
  # ----------------------------------------------------------------------------
  
  # 0) load  
  
  time_grid <- dataset$time_grid
  time_grid_est <- dataset$time_grid_est
  time_grid_both <- dataset$time_grid_both
  true_graphs <- dataset$true_graphs
  true_graph_indices <- sort(names(true_graphs)) 
  
  data_df <- convert_data_for_estimation(dataset$subject_data) %>% as.data.table()
  data_df$time <- data_df$time / dataset$simulation_params$T_max # normalize to [0, 1]
  
  X_k_truth <- dataset$X_k_truth
  X_k_coarse_truth <- dataset$X_k_coarse_truth
  X_k_both_truth <- dataset$X_k_both_truth
  
  n <- dim(dataset$Y_continuous)[1]
  
  p <- length(feature_sel)
  n <- length(patient_sel)
  m_est <- length(Tseq_est)
  # estimation --------------------------------------------------------------
  
  
  
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
  
  rho_ii_truth <- sweep(rho_ii_truth, 1:2, exp(kernel_params$base_cov), `*`) # each outer product --> hadamart with exp(K) (base cov's are the same for everyone)
  
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
  
  g_est_21 <- grid.arrange(visualize_log_intensity(recovered_intensities[1:5,], Tseq_est, 'Rho i est from X_hat') + geom_hline(yintercept = rho_i_mean), # rho_i_truth vs rho_i_est
                           g_data_2, 
                           g_GP_baseline, nrow = 1)
  
  arr_mat <- matrix(1:6, nrow = 2, byrow = F)
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
  X_k_truth_center <- X_k_truth - kernel_params$base_GP_mean
  
  
  mean_mat_coarse_truth <- apply(X_k_coarse_truth, c(1, 2), mean)  # result is p x m matrix
  X_k_coarse_truth_center_est <- sweep(X_k_coarse_truth, c(1, 2), mean_mat_coarse_truth, FUN = "-")   
  X_k_coarse_truth_center <- X_k_coarse_truth - kernel_params$base_GP_mean
  
  
  
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
  
  g_est_61 <- visualize_matrix_heatmap(log(K_c), 'K_c Matrix') # log values for better comparison
  
  
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
    
    kernel_params_i = dataset$true_graphs[[cont_ind]]$P_block_kronecker
    
    adj_mat_i <- true_graphs[[cont_ind]]$adj_mat
    
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
                             visualize_pm_block_matrix_heatmap(kernel_params_i$GP_simu_var, 'Ground Truth'), 
                             nrow = 1)
    
    # all 5 + ground truth
    g_est_84 <- grid.arrange(visualize_pm_block_matrix_heatmap(kernel_params_i$GP_simu_var, 'Ground Truth'), 
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
                             visualize_pm_block_matrix_heatmap(kernel_params_i$GP_simu_var, 'Ground Truth'), 
                             nrow = 1)    
    
    
    # all 5 + ground truth
    g_est_94 <- grid.arrange(visualize_pm_block_matrix_heatmap(kernel_params_i$GP_simu_var, 'Ground Truth'), 
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
    
    prec_truth <- kronecker(kernel_params_i$prec_mat, kernel_params_i$base_precision)
    prec_coarse_truth <- kronecker(kernel_params_i$prec_mat, kernel_params_i$base_precision_est)
    
    
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
    final_graph_estimates                  <- estimate_graph(P_cond_est,            C_cond_est,            V_cond_est, p)
    final_graph_estimates_X_coarse_truth   <- estimate_graph(P_cond_X_coarse_truth, C_cond_X_coarse_truth, V_cond_X_coarse_truth, p)
    final_graph_estimates_X_truth          <- estimate_graph(P_cond_X_truth,        C_cond_X_truth,        V_cond_X_truth, p)    
    final_graph_estimates_coarse_truth     <- estimate_graph(P_cond_coarse_truth,   C_cond_coarse_truth,   V_cond_coarse_truth, p)
    final_graph_estimates_truth            <- estimate_graph(P_cond_truth,          C_cond_truth,          V_cond_truth, p)
    
    g_est_111 <- grid.arrange(visualize_pm_block_matrix_heatmap(adj_mat_i, 'Ground Truth'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates$w_mat, 'Estimate'),
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_coarse_truth$w_mat, 'Coarse Truth X'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_truth$w_mat, 'Truth X'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_coarse_truth$w_mat, 'Coarse Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_truth$w_mat, 'Truth Theory'), 
                              layout_matrix = arr_mat)  
    
    
    
    # ROC curve 
    
    w_mat_ground_truth        <- hilbert_schmidt_norm_pm(prec_truth, p, m)
    w_mat_coarse_ground_truth <- hilbert_schmidt_norm_pm(prec_coarse_truth, p, m_est)
    
    diag(w_mat_ground_truth) <- 0
    diag(w_mat_coarse_ground_truth) <- 0
    
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
    
    
    roc_ground_truth        <- roc_with_threshold(w_mat_ground_truth,        adj_mat_i, 'Ground Truth')
    roc_coarse_ground_truth <- roc_with_threshold(w_mat_coarse_ground_truth, adj_mat_i, 'Coarse Ground Truth')
    roc_truth               <- roc_with_threshold(w_mat_truth,               adj_mat_i, 'Truth Theory')   
    roc_coarse_truth        <- roc_with_threshold(w_mat_coarse_truth,        adj_mat_i, 'Coarse Truth Theory')
    roc_X_truth             <- roc_with_threshold(w_mat_X_truth,             adj_mat_i, 'Truth X')
    roc_X_coarse_truth      <- roc_with_threshold(w_mat_X_coarse_truth,      adj_mat_i, 'Coarse Truth X')
    roc_est                 <- roc_with_threshold(w_mat_est,                 adj_mat_i, 'Estimate')
    
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
    
    if(terse){
      list(g_84  = g_est_84, 
           g_94  = g_est_94, 
           g_104 = g_est_104, 
           g_112 = g_est_112,  
           g_113 = g_est_113)
    } else{
      list(g_81  = g_est_81,  g_82  = g_est_82,  g_83  = g_est_83,  g_84  = g_est_84,
           g_91  = g_est_91,  g_92  = g_est_92,  g_93  = g_est_93,  g_94  = g_est_94,
           g_101 = g_est_101, g_102 = g_est_102, g_103 = g_est_103, g_104 = g_est_104,
           g_111 = g_est_111, g_112 = g_est_112, g_113 = g_est_113)      
    }
    
    
    
  }, mc.cores = ncores) # done with all y_c levels
  
  
  
  names(estimated_graphs_v2) <- apply(query_y_cs, 1, function(x){paste(x, collapse = '_')})
  
  
  # what to output
  
  if(terse){
    estimated_graphs_part_1 <- list(g_01 = g_01, g_21 = g_est_21, g_22 = g_est_22)
    
    return(list(estimated_graphs_part_1 = estimated_graphs_part_1,      
                estimated_graphs_part_2 = estimated_graphs_v2))      
    
  } else{
    estimated_graphs_part_1 <- list(g_01 = g_01, g_21 = g_est_21, g_22 = g_est_22,
                                    g_31 = g_est_31,
                                    g_41 = g_est_41, g_42 = g_est_42,
                                    g_51 = g_est_51, g_52 = g_est_52, g_53 = g_est_53,
                                    g_71 = g_est_71, g_72 = g_est_72)
    
    return(list(estimated_graphs_part_1 = estimated_graphs_part_1,      
                estimated_graphs_part_2 = estimated_graphs_v2))        
    
  }
  
  
  
}

full_conditional_estimation_CPGM <- function(processed_data, time_grid_est, method, terse, ncores, dir){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: troubleshoot the full conditional estimation procedure with ground truths between steps
  #
  # - instead of previously, where we only used G_ii and discarded G_ij
  # - here we actively look to estimate G_ij and use it
  # - 9/3/2025 - use G_ij^k, for each subject
  #
  # - CPGM method
  #
  # Input:
  #
  # - processed_data
  # - time_grid_est (m-dim vec) time discretization
  # - method   ('CPGM')
  # - terse    (boolean) if true, return much less
  # - ncores
  # - dir      (string) folder name, such as "simu_results_banded_trig2_2"
  # 
  #
  # ----------------------------------------------------------------------------
  
  # check for errors 
  if (!(method %in% c("CPGM"))) {
    stop("Error: method must be CPGM")
  }    
  
  
  arr_mat <- matrix(1:6, nrow = 2, byrow = F)
  arr_mat_8 <- matrix(1:8, nrow = 2, byrow = F)
  arr_mat_10 <- matrix(1:10, nrow = 2, byrow = F)  
  
  # part 0 - check ground truth precisions and process data --------------------
  
  # process data
  processed_data <- step_0_preprocess(dataset)
  data_df4 = processed_data[[1]]
  y_c_strata = processed_data[[2]]
  query_y_cs = processed_data[[3]]
  patient_sel = processed_data[[4]]
  feature_sel = processed_data[[5]]
  rm(processed_data)

  # load parameters
  
  n <- dim(y_c_strata)[1]
  m_est <- length(time_grid_est) 
    
  
  i_neq_j <- T
  rho_kernel <- T
  
  # part 1 - log intensities ---------------------------------------------------
  
  step_1 <- step_1_log_intensities(dataset, data_df4, time_grid_est, time_grid, time_grid_both)
  
  
  # now, we regress on y_c  ----------------------------------------------------
  
  print('at covariate loop')
  cont_inds <- 1:nrow(query_y_cs)
  estimated_graphs_v2 <- lapply(cont_inds, function(cont_ind) {  
    
    
    print(paste0(cont_ind, ' out of ', nrow(query_y_cs)))
    
    # obtain this query's parameters
    kernel_params_i = dataset$true_graphs[[cont_ind]]$P_block_kronecker
    adj_mat_i <- dataset$true_graphs[[cont_ind]]$adj_mat
    query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()    
    
    # part 2 - intensities -----------------------------------------------------
    
    step_2 <- step_2_rho_i(dataset, data_df4, kernel_params_i, rho_kernel, patient_sel, feature_sel,
                           time_grid, time_grid_est, query_y_c, y_c_strata, ncores)
    
    
    # part 2.1 - bivariate intensities -----------------------------------------
    
    step_2b <- step_2_rho_ij(step_1, step_2, kernel_params_i, i_neq_j)
    
    
    
    # part 3 - GP covariance estimation ----------------------------------------
    
    step_3 <- step_3_g_ij(step_2, step_2b, kernel_params_i, i_neq_j)
    
    
    # part 4 - eigendecomposition of GP covariance -----------------------------
    
    # careful about some processes being empty
    
    # step_4 <- step_4_eigendecomp(step_3, p, time_grid, time_grid_est)
    result <- tryCatch({
      step_4 <- step_4_eigendecomp(step_3, p, time_grid, time_grid_est)
    }, error = function(e) {
      cat("Error occurred in step_4, saving dataset...\n")
      save(dataset, file = paste0(dir, "/dataset.RData"))
      cat("Dataset saved to dataset.RData\n")
      stop(e)  # Re-throw the error
    })    
    
    
    # part 5 - covariance of KL coefficients -----------------------------------
    
    step_5 <- step_5_KL_covariance(step_3, step_4)
    
    # part 9 - conditional_correlation from KL covariance ----------------------
    
    step_9 <- step_9_C_cond_from_KL_cov(step_4, step_5, kernel_params_i)
    
    # part 10 - precision operator ---------------------------------------------
    
    block <- F
    MP <- F
    
    # step_10 <- step_10_P_cond(step_9, kernel_params_i, p, block, MP)
    
    result <- tryCatch({
      step_10 <- step_10_P_cond(step_9, kernel_params_i, p, block, MP)
    }, error = function(e) {
      cat("Error occurred in step_10, saving dataset...\n")
      save(dataset, file = paste0(dir, "/dataset.RData"))
      cat("Dataset saved to dataset.RData\n")
      stop(e)  # Re-throw the error
    })       
    
    # part 11 - HS norms of precision matrix -----------------------------------
    
    step_11 <- step_11_HS_norms(step_10, adj_mat_i, p)
    
    # part 12 - ROC curve ------------------------------------------------------
    
    step_12 <- step_12_ROC(step_11, adj_mat_i)
    
    # part 13 - metrics
    
    metrics <- get_metrics(list(step_10$P_cond_est_full,   # P_cond_est_full
                                step_9$C_cond_est_full,    # C_cond_est_full
                                NULL,
                                step_10$P_cond_coarse_ground_truth_full,   # P_cond_coarse_ground_truth_full
                                step_9$C_cond_corase_ground_truth_full,    # C_cond_coarse_ground_truth_full
                                NULL,
                                step_11$w_mat_est,   # w_mat_est
                                adj_mat_i,
                                step_12$roc_est,   # roc_est
                                step_2$rho_i_est,    # rho_i_est
                                step_2$rho_i_coarse_truth,    # rho_i_coarse_truth
                                step_2b$rho_ii_est,   # rho_ii_est
                                step_2b$rho_ii_coarse_truth,   # rho_ii_coarse_truth
                                step_3$g_ij_est,    # g_ij_est
                                step_3$g_ij_coarse_truth))   # g_ij_coarse_truth     
    
    
    # 8) save graphs (recall that we have ground truths of the form a_b_c)
    
    list(step_2 = step_2,
         step_2b = step_2b,
         step_3 = step_3,
         step_4 = step_4,
         step_5 = step_5,
         step_9 = step_9,
         step_10 = step_10,
         step_11 = step_11,
         step_12 = step_12,
         metrics = metrics,
         kernel_params_i = kernel_params_i)
    
    
    
  })# , mc.cores = ncores) # done with all y_c levels
  
  
  
  names(estimated_graphs_v2) <- apply(query_y_cs, 1, function(x){paste(x, collapse = '_')})
  
  
  # what to output
  
  
  estimated_graphs_part_1 <- list(g_01 = g_01, step_1 = step_1)
  
  return(list(estimated_graphs_part_1 = estimated_graphs_part_1,      
              estimated_graphs_part_2 = estimated_graphs_v2))      
  
  
  
}

