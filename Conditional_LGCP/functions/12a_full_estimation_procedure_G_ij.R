# use the G_ij route

full_conditional_estimation_with_truths_G_ij <- function(dataset, terse, ncores){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: troubleshoot the full conditional estimation procedure with ground truths between steps
  #
  # - instead of previously, where we only used G_ii and discarded G_ij
  # - here we actively look to estimate G_ij and use it
  # - 9/3/2025 - use G_ij^k, for each subject
  #
  # Input:
  #
  # - dataset
  # - terse    (boolean) if true, return much less
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
  true_graph_indices <- sort(names(true_graphs)) 
  
  data_df <- convert_data_for_estimation(dataset$subject_data) %>% as.data.table()
  data_df$time <- data_df$time / dataset$simulation_params$T_max # normalize to [0, 1]
  
  n <- dim(dataset$Y_continuous)[1]
  m_est <- length(time_grid_est)
  m     <- length(time_grid)
  

  
  # 0.1) check ground truth precision pxp matrices
  
  prec_ground_truths <- lapply(dataset$true_graphs, function(item) item$P_block_kronecker$prec_mat)
  
  lay_mat <- matrix(c(1:5, NA), nrow = 2)
  
  g_01 <- grid.arrange(visualize_nonneg_matrix_heatmap(prec_ground_truths[['1']], '1', -1, 1),
                       visualize_nonneg_matrix_heatmap(prec_ground_truths[['2']], '2', -1, 1),
                       visualize_nonneg_matrix_heatmap(prec_ground_truths[['3']], '3', -1, 1),
                       visualize_nonneg_matrix_heatmap(prec_ground_truths[['4']], '4', -1, 1),
                       visualize_nonneg_matrix_heatmap(prec_ground_truths[['5']], '5', -1, 1),
                       layout_matrix = lay_mat
  )  
  
  # check if any prec mat is too nonnegative!
  check_psd <- sapply(prec_ground_truths, function(M) {
    min(eigen(M, symmetric = TRUE, only.values = TRUE)$values) > -1e-10
  })  
  
  if(! any(check_psd)){
    warning("code 01: some ground truth prec mats are not psd, proceeding anyway")
  }
  

  
  # 3) pre-processing before estimation ------------------------------------------
  
  
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
  
  # part 1 - log intensities
  #
  # - X_k_truth                 (p x m_truth x n)
  # - X_k_coarse_truth          (p x m_est   x n)
  # - X_k_both_truth            (p x m_both  x n)
  # - X_k_est                   (p x m_est   x n)
  #
  # - Lambda_k_truth            (p x m_truth x n)
  # - Lambda_k_coarse_truth     (p x m_est   x n)

  
  
  X_k_est <- subject_specific_log_intensity(data_df4, Tseq_est)
  X_k_truth <- dataset$X_k_truth
  X_k_coarse_truth <- dataset$X_k_coarse_truth
  X_k_both_truth <- dataset$X_k_both_truth
  
  Lambda_k_truth <- exp(X_k_truth)
  Lambda_k_coarse_truth <- exp(X_k_coarse_truth)
  
  
  # 1.1) visualize true log-intensities - good
  
  arr_mat <- matrix(1:6, nrow = 2, byrow = F)
  arr_mat_8 <- matrix(1:8, nrow = 2, byrow = F)
  arr_mat_10 <- matrix(1:10, nrow = 2, byrow = F)
  g_11 <- grid.arrange(visualize_log_intensity(X_k_est[1:5,,1], time_grid_est, 'Estimate'),
                       visualize_log_intensity(X_k_coarse_truth[1:5,,1], time_grid_est, 'Coarser Truth'),
                       visualize_log_intensity(X_k_truth[1:5,,1], time_grid, 'Finer Truth'),
                       visualize_log_intensity(X_k_both_truth[1:5,,1], time_grid_both, 'Combined Truth'),
                       textGrob("0. Log Intensity\n of first replicate", gp = gpar(fontsize = 14)),
                       layout_matrix = arr_mat)   
  
  # part 2 - intensities
  #
  # - rho_i_truth               (p x m_truth)
  # - rho_i_coarse_truth        (p x m_est)
  # - rho_i_X_truth             (p x m_truth)
  # - rho_i_X_coarse_truth      (p x m_est)
  # - rho_i_est                 (p x m_est)
  # 
  # - rho_i_k_est               (p x m_est x n)
  
  
  # subject-level intensity
  rho_i_k_est <- exp(X_k_est)
  
  
  
  # rho_i_truth = exp(mu(t) + 0.5 * diag(GP_cov))
  kernel_params       <- true_graphs$`1`$P_block_kronecker
  rho_i_truth_value   <- exp(kernel_params$base_GP_mean + 0.5 * kernel_params$base_variance)
  rho_i_truth         <- replicate(p, exp(kernel_params$GP_simu_mean + 0.5 * diag(kernel_params$base_cov))) %>% t() 
  rho_i_coarse_truth  <- replicate(p, exp(kernel_params$GP_simu_mean_est + 0.5 * diag(kernel_params$base_cov_est))) %>% t()  

  # rho_i (mean intensity) from X_truths
  mat_list <- lapply(dataset$subject_data, function(x) exp(x$X_functions))
  rho_i_X_truth <- Reduce("+", mat_list) / length(mat_list)
  
  mat_list_coarse <- lapply(dataset$subject_data, function(x) exp(x$X_functions_coarse))
  rho_i_X_coarse_truth <- Reduce("+", mat_list_coarse) / length(mat_list_coarse)    
  
  # rho_i_est from data
  rho_list <- estimate_intensities_stratum_parallel_v4(data_df4, patient_sel, feature_sel, Tseq_est, F, ncores)
  rho_i_est <- do.call(rbind, lapply(rho_list[[1]], function(v) as.numeric(v))) # p x m matrix of estimated mean intensities  
  
  
  
  # part 2 visualization 


  g_22 <- grid.arrange(visualize_log_intensity(rho_i_truth[1:5,], time_grid, 'Truth Theory'),
                       visualize_log_intensity(rho_i_coarse_truth[1:5,], time_grid_est, 'Coarse Truth Theory'),
                       visualize_log_intensity(rho_i_X_truth[1:5,], time_grid, 'Truth X'),
                       visualize_log_intensity(rho_i_X_coarse_truth[1:5,], time_grid_est, 'Coarse Truth X'),
                       textGrob("1. Rho_i\nEstimation for \n first 5 processes", gp = gpar(fontsize = 14)),
                       visualize_log_intensity(rho_i_est[1:5,], time_grid_est, 'Estimate'),
                       layout_matrix = arr_mat)
  
  
  g_23 <- visualize_log_intensity(t(rho_i_k_est[1,,1:5]), time_grid_est, 'Rho_i_k', 'Replicate')
  
  
  # checking
  rho_i_est_check <- apply(rho_i_k_est, c(1,2), mean)
  summary(as.numeric(rho_i_est - rho_i_est_check))    # rho_i_k_est is correct

  # part 2.1 - bivariate intensities -------------------------------------------
  #
  # - rho_i_truth               (p x m_truth)
  # - rho_i_coarse_truth        (p x m_est)
  # - rho_i_X_truth             (p x m_truth)
  # - rho_i_X_coarse_truth      (p x m_est)
  # - rho_ii_est                (p-dim list, m_est x m_est)
  # 
  # - rho_ij_k_est               (p x m_est x n)  

  
  # subject-level bivariate intensity
  rho_ij_k_est <- subject_specific_bivariate_intensity(data_df4, Tseq_est, Tseq_est, 'i', ncores)  

  
  rho_ii_est <- rho_list[[2]]
  rho_ii_truth <- list()
  rho_ii_coarse_truth <- list()
  rho_ii_X_truth <- list()
  rho_ii_X_coarse_truth <- list()
  for(i in 1:p){
    for(j in i:i){
      key <- paste0(i, '_', j)
      
      # theory
      rho_ii_truth[[key]]        <- tcrossprod(rho_i_truth[i,],        rho_i_truth[j,])        * exp(extract_block_structure_ij(kernel_params$GP_simu_var, m, i, j))
      rho_ii_coarse_truth[[key]] <- tcrossprod(rho_i_coarse_truth[i,], rho_i_coarse_truth[j,]) * exp(extract_block_structure_ij(kernel_params$GP_simu_var_est, m_est, i, j))
      
      # X_k portion
      mats <- array(0, dim = c(m, m))
      coarse_mats <- array(0, dim = c(m_est, m_est))  
      
      for (k in 1:n) {
        vi <- Lambda_k_truth[i, , k]
        vj <- Lambda_k_truth[j, , k]
        mats <- mats + tcrossprod(vi, vj)  # faster than outer(v, v)
        
        vi <- Lambda_k_coarse_truth[i, , k]
        vj <- Lambda_k_coarse_truth[j, , k]
        coarse_mats <- coarse_mats + tcrossprod(vi, vj)  # faster than outer(v, v) 
        
        # Average across replicates
        rho_ii_X_truth[[key]] <- mats / n
        rho_ii_X_coarse_truth[[key]] <- coarse_mats / n      
      }
    }
  }
  

  g_est_22 <- grid.arrange(visualize_matrix_heatmap(rho_ii_truth[['1_1']], 'Truth Theory', 40000, 200000),
                           visualize_matrix_heatmap(rho_ii_coarse_truth[['1_1']], 'Coarse Truth Theory', 40000, 200000),
                           visualize_matrix_heatmap(rho_ii_X_truth[['1_1']], 'Truth X', 40000, 200000),
                           visualize_matrix_heatmap(rho_ii_X_coarse_truth[['1_1']], 'Coarse Truth X', 40000, 200000),
                           textGrob("2. Rho ii\nEstimation \n for process 1", gp = gpar(fontsize = 14)),
                           visualize_matrix_heatmap(rho_ii_est[['1_1']], 'Estimate', 40000, 200000),
                           layout_matrix = arr_mat)
  
  
  # verify - take the sample mean of rho_ij_k_est where i = j
  
  df_pair <- rho_ij_k_est$pair_ID
  keep_idx <- which(df_pair[[1]] == df_pair[[2]])
  dt_selected <- rho_ij_k_est$bivariate_intensities[, ..keep_idx] 
  
  rho_ii_est_2 <- lapply(dt_selected, function(col_list) {
    # col_list is a list of 30 matrices (one per row)
    arr <- array(unlist(col_list), dim = c(m_est, m_est, length(col_list)))
    apply(arr, c(1, 2), mean)  # average over 3rd dimension
  })
  
  summary(as.numeric(rho_ii_est_2[[1]] - rho_ii_est[[1]])) # looks good!
  

  # part 3 - GP covariance estimation ------------------------------------------
  
  
  
  
  # bundle rho's to prepare for covariance function estimation
  # - rho_list_truth        (p-dim list, with rho_i and rho_ii_mat)
  # - rho_list_coarse_truth (p-dim list, with rho_i and rho_ii_mat)
  # - rho_list_X_truth             
  # - rho_list_X_coarse_truth      
  # - rho_list_est                 
  
  
  rho_list_truth <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_i_truth[i, ],     
      rho_ii_mat = rho_ii_truth[[i]]  
    )
  }) 
  
  rho_list_coarse_truth <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_i_coarse_truth[i, ],     
      rho_ii_mat = rho_ii_coarse_truth[[i]] 
    )
  })   
  
  rho_list_X_truth <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_i_X_truth[i, ],     
      rho_ii_mat = rho_ii_X_truth[[i]]  
    )
  })
  
  rho_list_X_coarse_truth <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_i_X_coarse_truth[i, ],   
      rho_ii_mat = rho_ii_X_coarse_truth[[i]]   
    )
  }) 
  
  rho_list_est <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_i_est[i, ],   
      rho_ii_mat = rho_ii_est[[i]]   
    )
  }) 
  
  #  covariance function estimation
  # - g_ii_truth                (m x m x p)
  # - g_ii_coarse_truth         (m_est x m_est x p)
  # - g_ii_X_truth              (m x m x p)
  # - g_ii_X_coarse_truth       (m_est x m_est x p)
  # - g_ii_est                  (m_est x m_est x p)
  #
  # - g_ij_k_est      (list of all i_j)   (m_est x m_est x n)
  # - g_ij_est_later  (list of all i_j)   (m_est x m_est)
  # - g_ii_est_later            (m_est x m_est x p)
  
  g_ii_est            <- estimate_covariance_functions_ii(rho_list_est)
  g_ii_X_coarse_truth <- estimate_covariance_functions_ii(rho_list_X_coarse_truth)
  g_ii_X_truth        <- estimate_covariance_functions_ii(rho_list_X_truth)
  g_ii_coarse_truth   <- estimate_covariance_functions_ii(rho_list_coarse_truth)
  g_ii_truth          <- estimate_covariance_functions_ii(rho_list_truth)
  

  # subject-specific G
  g_ij_k_est <- estimate_subject_specific_covariance(rho_i_k_est, rho_ij_k_est$bivariate_intensities, rho_ij_k_est$pair_ID)
  
  # take the average across all subjects
  g_ij_est_later <- lapply(g_ij_k_est, function(mat) {
    apply(mat, c(1, 2), mean)
  })
  

  g_ii_est_later <- array(
    unlist(g_ij_est_later[keep_idx]),
    dim = c(m_est, m_est, p)
  )
  
  # does jensen's hold?
  
  summary(as.numeric(g_ii_est[,,1] - g_ii_est_later[,,1]))
  
  i_val <- 1
  j_val <- 1
  k_val <- 3
  pair_val <- 1 # 2 and 8
  g_ij_k_val <- g_ij_k_est[[pair_val]][,,k_val]
  rho_i_k_val <- rho_i_k_est[i_val,,k_val]
  rho_j_k_val <- rho_i_k_est[j_val,,k_val]
  rho_ij_k_val <- rho_ij_k_est$bivariate_intensities[[pair_val]][[k_val]]
  
  g_ij_k_val_2 <- log(rho_ij_k_val / outer(rho_i_k_val, rho_j_k_val))
  
  summary(as.numeric(g_ij_k_val - g_ij_k_val_2))
  
  # visualization
  g_est_31 <- grid.arrange(visualize_matrix_heatmap(g_ii_truth[,,1], 'Truth Theory', -1, 1),  
                           visualize_matrix_heatmap(g_ii_coarse_truth[,,1], 'Coarse Truth Theory', -1, 1),
                           visualize_matrix_heatmap(g_ii_X_truth[,,1], 'Truth X', -1, 1),    
                           visualize_matrix_heatmap(g_ii_X_coarse_truth[,,1], 'Coarse Truth X', -1, 1),
                           visualize_matrix_heatmap(g_ii_est[,,1], 'Estimate', -1, 1),
                           visualize_matrix_heatmap(g_ii_est_later[,,1], 'Estimate Later Avg', -1, 1),
                           textGrob("3. Covariance Function\nEstimation (G_ii)", gp = gpar(fontsize = 14)),
                           visualize_matrix_heatmap(g_ij_k_est[[1]][,,1], 'Single Subject Estimate', -1, 1),
                           layout_matrix = arr_mat_8)
  
  
 
  
  
  # part 4 - eigendecomposition of GP covariance -------------------------------
  
  # perform eigendecomposition
  
  eigen_decomp_est            <- compute_eigendecomposition_ii(g_ii_est)
  eigen_decomp_est_later      <- compute_eigendecomposition_ii(g_ii_est_later)
  eigen_decomp_X_coarse_truth <- compute_eigendecomposition_ii(g_ii_X_coarse_truth)
  eigen_decomp_X_truth        <- compute_eigendecomposition_ii(g_ii_X_truth)  
  eigen_decomp_coarse_truth   <- compute_eigendecomposition_ii(g_ii_coarse_truth)
  eigen_decomp_truth          <- compute_eigendecomposition_ii(g_ii_truth)
  
  
  # validate
  g_ii_est_decomp             <- validate_eigendecomposition_ii(g_ii_est, eigen_decomp_est)
  g_ii_est_later_decomp       <- validate_eigendecomposition_ii(g_ii_est_later, eigen_decomp_est_later)
  g_ii_X_coarse_truth_decomp  <- validate_eigendecomposition_ii(g_ii_X_coarse_truth, eigen_decomp_X_coarse_truth)  
  g_ii_X_truth_decomp         <- validate_eigendecomposition_ii(g_ii_X_truth, eigen_decomp_X_truth)  
  g_ii_coarse_truth_decomp    <- validate_eigendecomposition_ii(g_ii_coarse_truth, eigen_decomp_coarse_truth)  
  g_ii_truth_decomp           <- validate_eigendecomposition_ii(g_ii_truth, eigen_decomp_truth)
  
  
  
  # 4.1) do the eigenfunctions between coarse and regular truths differ? - yes - the g_ij is too different from the truth
  
  g_est_41 <- grid.arrange(visualize_log_intensity(t(eigen_decomp_truth$eigenfunctions[[1]]),           time_grid, 'Truth Theory'),
                           visualize_log_intensity(t(eigen_decomp_coarse_truth$eigenfunctions[[1]]),    Tseq_est,  'Coarse Truth Theory'),
                           visualize_log_intensity(t(eigen_decomp_X_truth$eigenfunctions[[1]]),         time_grid, 'Truth X'),
                           visualize_log_intensity(t(eigen_decomp_X_coarse_truth$eigenfunctions[[1]]),  Tseq_est,  'Coarse Truth X'),
                           visualize_log_intensity(t(eigen_decomp_est$eigenfunctions[[1]]),             Tseq_est,  'Estimate'),
                           visualize_log_intensity(t(eigen_decomp_est_later$eigenfunctions[[1]]),       Tseq_est,  'Estimate Later'),
                           textGrob("4. Eigenfunction\nBasis of\nProcess 1", gp = gpar(fontsize = 14)),
                           layout_matrix = arr_mat_8)
  
  # 4.2) are all eigenfunctions orthogonal? - yes
  
  t(eigen_decomp_est$eigenfunctions[[1]]) %*% eigen_decomp_est$eigenfunctions[[1]]
  t(eigen_decomp_truth$eigenfunctions[[1]]) %*% eigen_decomp_truth$eigenfunctions[[1]]
  t(eigen_decomp_coarse_truth$eigenfunctions[[1]]) %*% eigen_decomp_coarse_truth$eigenfunctions[[1]]
  
  # did eigendecomposition reconstruct the g matrix to the best of its ability? yes
  
  summary(as.numeric(g_ii_est[,,1] - g_ii_est_decomp[,,1])) # good (-3e-3, 3e-3)
  summary(as.numeric(g_ii_truth[,,1] - g_ii_truth_decomp[,,1])) # good (-3e-3, 5e-3)
  summary(as.numeric(g_ii_coarse_truth[,,1] - g_ii_coarse_truth_decomp[,,1])) # very good (-7e-4, 9e-4)
  
  
  # visualize g_ij vs its eigendecomposition reconstruction
  
  
  g_est_42 <- grid.arrange(visualize_matrix_heatmap(g_ii_truth[,,1], 'Truth Theory (TT)', -1, 1),  # ground truth 50 x 50 covariance
                           visualize_matrix_heatmap(g_ii_truth_decomp[,,1], 'TT Reconstruct', -1, 1), # ground truth reconstructed
                           visualize_matrix_heatmap(g_ii_coarse_truth[,,1], 'Coarse Truth Theory (CTT)', -1, 1), # ground truth 19 x 19 covariance
                           visualize_matrix_heatmap(g_ii_coarse_truth_decomp[,,1], 'CTT Reconstruct', -1, 1), # ground truth 19 x 19 covariance
                           visualize_matrix_heatmap(g_ii_est[,,1], 'Estimate (E)', -1, 1),    # estimate 19 x 19 covariance
                           visualize_matrix_heatmap(g_ii_est_decomp[,,1], 'E Reconstruct', -1, 1), # reconstructed estimate
                           visualize_matrix_heatmap(g_ii_est_later[,,1], 'Estimate Later (EL)', -1, 1),    # estimate 19 x 19 covariance
                           visualize_matrix_heatmap(g_ii_est_later_decomp[,,1], 'EL Reconstruct', -1, 1), # reconstructed estimate   
                           textGrob("4. Eigenfunction\nReconstruction", gp = gpar(fontsize = 14)),
                           layout_matrix = arr_mat_10)    
  
  # part 5 - KL expansion ------------------------------------------------------
  
  
  # 5.1) center X_k_est and X_k_truth
  
  # 
  # - X_k_est_center                (p x m_est x n)
  # - X_k_truth_center              (p x m x n)       X_k_truth - acutal mean
  # - X_k_truth_center_est          (p x m x n)       X_k_truth - empirical mean
  # - X_k_coarse_truth_center       (p x m_est x n)   X_k_coarse_truth - actual mean
  # - X_k_coarse_truth_center_est   (p x m_est x n)   X_k_coarse_truth - empirical mean
  
  mean_mat_est <- apply(X_k_est, c(1, 2), mean)  # result is p x m matrix
  X_k_est_center <- sweep(X_k_est, c(1, 2), mean_mat_est, FUN = "-")
  
  
  mean_mat_truth <- apply(X_k_truth, c(1, 2), mean)  # result is p x m matrix
  X_k_truth_center_est <- sweep(X_k_truth, c(1, 2), mean_mat_truth, FUN = "-") 
  X_k_truth_center <- X_k_truth - kernel_params$base_GP_mean
  
  
  mean_mat_coarse_truth <- apply(X_k_coarse_truth, c(1, 2), mean)  # result is p x m matrix
  X_k_coarse_truth_center_est <- sweep(X_k_coarse_truth, c(1, 2), mean_mat_coarse_truth, FUN = "-")   
  X_k_coarse_truth_center <- X_k_coarse_truth - kernel_params$base_GP_mean
  
  
  
  # 5.2) get KL coeffs  
  # 
  # - kl_coeffs_v2                     (n x p x d)   coeffs from X_k_est
  # - kl_coeffs_X_coarse_truth_v2      (n x p x d)   coeffs from X_k_coarse_truth
  # - kl_coeffs_X_truth_v2             (n x p x d)   coeffs from X_k_truth
  # - kl_coeffs_coarse_truth_v2        (n x p x d)   coeffs from X_k_coarse_truth but eigendecomp from theory
  # - kl_coeffs_truth_v2               (n x p x d)   coeffs from X_k_truth        but eigendecomp from theory
  
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
  
  # estimate and truth theory + reconstruct
  g_est_51 <- ggplot() + 
    geom_line(data = df_both, aes(x = Var2, y = value, color = cat)) + 
    facet_wrap(~ Var1) + 
    ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
    ylab('Log Intensity') + 
    xlab('Time') + 
    labs(color = "Truth or Estimate") + 
    theme_bw()
  
  # estimate, truth x and truth theory + reconstruct 
  g_est_52 <- ggplot() + 
    geom_line(data = df_all, aes(x = Var2, y = value, color = cat)) + 
    facet_wrap(~ Var1) + 
    ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
    ylab('Log Intensity') + 
    xlab('Time') + 
    labs(color = "Truth or Estimate") + 
    theme_bw()   
  
  # estimate, truth x and truth theory + no reconstruct 
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
  
  # get c_ij

  c_ij_ab_k       <- estimate_cross_kl_coefficients_parallel(g_ij_k_est, df_pair, eigen_decomp_est$eigenfunctions,       Tseq_est, ncores)
  c_ij_ab_k_later <- estimate_cross_kl_coefficients_parallel(g_ij_k_est, df_pair, eigen_decomp_est_later$eigenfunctions, Tseq_est, ncores)  
  
  # get M_hat 
  M_hat_g_ij <- estimate_rkhs_regression_operators(c_ij_ab_k, df_pair, y_c_strata, gamma_c, ncores)
  
  # first model ----------------------------------------------------------------
  
  # list of length (p+1) choose 2
  
  # - V_YcXij_est                (list of i_j entries) --> (n x d x d matrix)
  # - V_YcXij_X_coarse_truth     (list of i_j entries) --> (n x d x d matrix)
  # - V_YcXij_X_truth            (list of i_j entries) --> (n x d x d matrix)
  # - V_YcXij_coarse_truth       (list of i_j entries) --> (n x d x d matrix)
  # - V_YcXij_truth              (list of i_j entries) --> (n x d x d matrix)
  
  
  V_YcXij_est             <- construct_cross_covariance_matrix_v3(kl_coeffs_v2,                ncores)
  V_YcXij_X_coarse_truth  <- construct_cross_covariance_matrix_v3(kl_coeffs_X_coarse_truth_v2, ncores)
  V_YcXij_X_truth         <- construct_cross_covariance_matrix_v3(kl_coeffs_X_truth_v2,        ncores)  
  V_YcXij_coarse_truth    <- construct_cross_covariance_matrix_v3(kl_coeffs_coarse_truth_v2,   ncores)
  V_YcXij_truth           <- construct_cross_covariance_matrix_v3(kl_coeffs_truth_v2,          ncores)
  
  
  # - M_hat_est                (list of i_j entries) --> (n x d_i x d_j matrix)
  # - M_hat_X_coarse_truth     (list of i_j entries) --> (n x d_i x d_j matrix)
  # - M_hat_X_truth            (list of i_j entries) --> (n x d_i x d_j matrix)
  # - M_hat_coarse_truth       (list of i_j entries) --> (n x d_i x d_j matrix)
  # - M_hat_truth              (list of i_j entries) --> (n x d_i x d_j matrix) 
  
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
    
    # - V_cond_g_ij             (list of i_j entries) --> (m_est x m_est matrix)
    # - V_cond_est              (list of i_j entries) --> (m_est x m_est matrix)
    # - V_cond_X_coarse_truth   (list of i_j entries) --> (m_est x m_est matrix)
    # - V_cond_X_truth          (list of i_j entries) --> (m x m matrix)
    # - V_cond_coarse_truth     (list of i_j entries) --> (m_est x m_est matrix)
    # - V_cond_truth            (list of i_j entries) --> (m x m matrix)
    
    V_cond_g_ij           <- evaluate_regression_at_query_g_ij(M_hat_g_ij, df_pair, y_c_strata, query_y_c, eigen_decomp_est$eigenfunctions,            gamma_c)
    V_cond_est            <- evaluate_regression_at_query_v2(M_hat_est,             y_c_strata, query_y_c, eigen_decomp_est$eigenfunctions,            gamma_c, p)
    V_cond_X_coarse_truth <- evaluate_regression_at_query_v2(M_hat_X_coarse_truth,  y_c_strata, query_y_c, eigen_decomp_X_coarse_truth$eigenfunctions, gamma_c, p)
    V_cond_X_truth        <- evaluate_regression_at_query_v2(M_hat_X_truth,         y_c_strata, query_y_c, eigen_decomp_X_truth$eigenfunctions,        gamma_c, p)    
    V_cond_coarse_truth   <- evaluate_regression_at_query_v2(M_hat_coarse_truth,    y_c_strata, query_y_c, eigen_decomp_coarse_truth$eigenfunctions,   gamma_c, p)
    V_cond_truth          <- evaluate_regression_at_query_v2(M_hat_truth,           y_c_strata, query_y_c, eigen_decomp_truth$eigenfunctions,          gamma_c, p)
    
    
    
    # use KL coefficients to obtain V_cond
    # V_cond_truth_12 <- validate_V_ij_v2(kl_coeffs_truth_v2, eigen_decomp_truth$eigenfunctions, 1, 2, y_c_strata, query_y_c, gamma_c)
    # V_cond_truth_11 <- validate_V_ij_v2(kl_coeffs_truth_v2, eigen_decomp_truth$eigenfunctions, 1, 1, y_c_strata, query_y_c, gamma_c)
    
    
    # 8.2) visualizations 
    
    # for pair 1_2, plot estimate and truths - this one should look correlated
    g_est_81 <- grid.arrange(visualize_matrix_heatmap(V_cond_truth[['1_2']], 'Truth Theory'),
                             visualize_matrix_heatmap(V_cond_coarse_truth[['1_2']], 'Coarse Truth Theory'),
                             visualize_matrix_heatmap(V_cond_X_truth[['1_2']], 'Truth X'), 
                             visualize_matrix_heatmap(V_cond_X_coarse_truth[['1_2']], 'Coarse Truth X'), 
                             visualize_matrix_heatmap(V_cond_est[['1_2']], 'Estimate'), 
                             visualize_matrix_heatmap(V_cond_g_ij[['1_2']], 'Estimate G_ij'), 
                             textGrob("8. Conditional\nCovariance Operator\nProcess 1 with 2", gp = gpar(fontsize = 14)),
                             layout_matrix = arr_mat_8)    
    
    # but these two versions of the truth aren't exactly the same
    # grid.arrange(visualize_pm_block_matrix_heatmap(V_cond_truth[['1_2']]), 
    #              visualize_pm_block_matrix_heatmap(V_cond_truth_12), 
    #              layout_matrix = arr_mat)
    
    
    # visualize again for process 3 with 6 - they should not look correlated
    g_est_82 <- grid.arrange(visualize_matrix_heatmap(V_cond_truth[['3_6']], 'Truth Theory'),
                             visualize_matrix_heatmap(V_cond_coarse_truth[['3_6']], 'Coarse Truth Theory'),
                             visualize_matrix_heatmap(V_cond_X_truth[['3_6']], 'Truth X'), 
                             visualize_matrix_heatmap(V_cond_X_coarse_truth[['3_6']], 'Coarse Truth X'), 
                             visualize_matrix_heatmap(V_cond_est[['3_6']], 'Estimate'), 
                             visualize_matrix_heatmap(V_cond_g_ij[['3_6']], 'Estimate G_ij'), 
                             textGrob("8. Conditional\nCovariance Operator\nProcess 3 with 6", gp = gpar(fontsize = 14)),
                             layout_matrix = arr_mat_8) 
    
    
    
    
    # 8.3) visualize the entire pm x pm matrix
     
    V_cond_g_ij_full           <- assemble_block_matrix_v2(V_cond_g_ij,           p, m_est)
    V_cond_est_full            <- assemble_block_matrix_v2(V_cond_est,            p, m_est)
    V_cond_X_coarse_truth_full <- assemble_block_matrix_v2(V_cond_X_coarse_truth, p, m_est)
    V_cond_X_truth_full        <- assemble_block_matrix_v2(V_cond_X_truth,        p, m)    
    V_cond_coarse_truth_full   <- assemble_block_matrix_v2(V_cond_coarse_truth,   p, m_est)
    V_cond_truth_full          <- assemble_block_matrix_v2(V_cond_truth,          p, m)
    
    
    
    
    # est vs truth vs ground truth
    g_est_83 <- grid.arrange(visualize_pm_block_matrix_heatmap(V_cond_est_full, 'Estimate'),
                             visualize_pm_block_matrix_heatmap(V_cond_truth_full, 'Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(kernel_params_i$GP_simu_var, 'Ground Truth'), 
                             nrow = 1)
    
    # all 5 + ground truth
    g_est_84 <- grid.arrange(visualize_pm_block_matrix_heatmap(kernel_params_i$GP_simu_var, 'Ground Truth'), 
                             textGrob("8. Conditional\nCovariance Operator\nAll Processes", gp = gpar(fontsize = 14)),
                             visualize_pm_block_matrix_heatmap(V_cond_est_full, 'Estimate'),
                             visualize_pm_block_matrix_heatmap(V_cond_g_ij_full, 'Estimate G_ij'),
                             visualize_pm_block_matrix_heatmap(V_cond_truth_full, 'Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(V_cond_coarse_truth_full, 'Coarse Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(V_cond_X_truth_full, 'Truth X'), 
                             visualize_pm_block_matrix_heatmap(V_cond_X_coarse_truth_full, 'Coarse Truth X'), 
                             layout_matrix = arr_mat_8)    
    
    
    
    # part 9 - correlation operator --------------------------------------------
    
    C_cond_g_ij           <- estimate_conditional_correlation_v4(V_cond_g_ij, p)
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
                             visualize_pm_block_matrix_heatmap(C_cond_est[['1_2']], 'Estimate'),
                             visualize_pm_block_matrix_heatmap(C_cond_g_ij[['1_2']], 'Estimate G_ij'),
                             textGrob("9. Conditional\nCorrelation Operator\nProcess 1 with 2", gp = gpar(fontsize = 14)),                             
                             layout_matrix = arr_mat_8)  
    
    g_est_92 <- grid.arrange(visualize_pm_block_matrix_heatmap(C_cond_truth[['3_6']], 'Truth Theory'),
                             visualize_pm_block_matrix_heatmap(C_cond_coarse_truth[['3_6']], 'Coarse Truth Theory'),
                             visualize_pm_block_matrix_heatmap(C_cond_X_truth[['3_6']], 'Truth X'),
                             visualize_pm_block_matrix_heatmap(C_cond_X_coarse_truth[['3_6']], 'Coarse Truth X'),
                             visualize_pm_block_matrix_heatmap(C_cond_est[['3_6']], 'Estimate'),
                             visualize_pm_block_matrix_heatmap(C_cond_g_ij[['3_6']], 'Estimate G_ij'),
                             textGrob("9. Conditional\nCorrelation Operator\nProcess 3 with 6", gp = gpar(fontsize = 14)),                             
                             layout_matrix = arr_mat_8)      
    
    # 9.2) visualize the entire pm x pm block 
    
    C_cond_g_ij_full           <- assemble_block_matrix_v2(C_cond_g_ij,           p, m_est)
    C_cond_est_full            <- assemble_block_matrix_v2(C_cond_est,            p, m_est)
    C_cond_X_coarse_truth_full <- assemble_block_matrix_v2(C_cond_X_coarse_truth, p, m_est)
    C_cond_X_truth_full        <- assemble_block_matrix_v2(C_cond_X_truth,        p, m)    
    C_cond_coarse_truth_full   <- assemble_block_matrix_v2(C_cond_coarse_truth,   p, m_est)
    C_cond_truth_full          <- assemble_block_matrix_v2(C_cond_truth,          p, m)
    
    
    g_est_93 <- grid.arrange(visualize_pm_block_matrix_heatmap(C_cond_est_full, 'Estimate'),
                             visualize_pm_block_matrix_heatmap(C_cond_truth_full, 'Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(kernel_params_i$GP_simu_var, 'Ground Truth'), 
                             nrow = 1)    
    
    
    # all 5 + ground truth
    g_est_94 <- grid.arrange(visualize_pm_block_matrix_heatmap(kernel_params_i$GP_simu_var, 'Ground Truth'), 
                             textGrob("9. Conditional\nCorrelation Operator\nAll Processes", gp = gpar(fontsize = 14)),
                             visualize_pm_block_matrix_heatmap(C_cond_est_full, 'Estimate'),
                             visualize_pm_block_matrix_heatmap(C_cond_g_ij_full, 'Estimate G_ij'),
                             visualize_pm_block_matrix_heatmap(C_cond_truth_full, 'Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(C_cond_coarse_truth_full, 'Coarse Truth Theory'), 
                             visualize_pm_block_matrix_heatmap(C_cond_X_truth_full, 'Truth X'), 
                             visualize_pm_block_matrix_heatmap(C_cond_X_coarse_truth_full, 'Coarse Truth X'), 
                             layout_matrix = arr_mat_8)       
    
    
    # part 10 - precision operator ---------------------------------------------
    
    P_cond_g_ij           <- estimate_precision_operator_v3(C_cond_g_ij, p)
    P_cond_est            <- estimate_precision_operator_v3(C_cond_est, p)
    P_cond_X_coarse_truth <- estimate_precision_operator_v3(C_cond_X_coarse_truth, p)
    P_cond_X_truth        <- estimate_precision_operator_v3(C_cond_X_truth, p)    
    P_cond_coarse_truth   <- estimate_precision_operator_v3(C_cond_coarse_truth, p)
    P_cond_truth          <- estimate_precision_operator_v3(C_cond_truth, p)
    
    prec_truth <- kronecker(kernel_params_i$prec_mat, kernel_params_i$base_precision)
    prec_coarse_truth <- kronecker(kernel_params_i$prec_mat, kernel_params_i$base_precision_est)
    
    
    P_ground_truth_12 <- extract_block_structure_ij(prec_truth, m, 1, 2)
    P_ground_truth_36 <- extract_block_structure_ij(prec_truth, m, 3, 6)
    
    # 10.1) visualize
    
    g_est_101 <- grid.arrange(visualize_pm_block_matrix_heatmap(P_ground_truth_12, 'Truth Ground'),
                              textGrob("10. Conditional\nPrecision Operator\n Process 1 and 2", gp = gpar(fontsize = 14)),
                              visualize_pm_block_matrix_heatmap(P_cond_est[['1_2']], 'Estimate'),
                              visualize_pm_block_matrix_heatmap(P_cond_g_ij[['1_2']], 'Estimate G_ij'),
                              visualize_pm_block_matrix_heatmap(P_cond_truth[['1_2']], 'Truth Theory'),
                              visualize_pm_block_matrix_heatmap(P_cond_coarse_truth[['1_2']], 'Coarse Truth Theory'),
                              visualize_pm_block_matrix_heatmap(P_cond_X_truth[['1_2']], 'Truth X'),
                              visualize_pm_block_matrix_heatmap(P_cond_X_coarse_truth[['1_2']], 'Coarse Truth X'),
                              layout_matrix = arr_mat_8)  
    
    g_est_102 <- grid.arrange(visualize_pm_block_matrix_heatmap(P_ground_truth_36, 'Truth Ground'),
                              textGrob("10. Conditional\nPrecision Operator\n Process 3 and 6", gp = gpar(fontsize = 14)),
                              visualize_pm_block_matrix_heatmap(P_cond_est[['3_6']], 'Estimate'),
                              visualize_pm_block_matrix_heatmap(P_cond_g_ij[['3_6']], 'Estimate G_ij'),
                              visualize_pm_block_matrix_heatmap(P_cond_truth[['3_6']], 'Truth Theory'),
                              visualize_pm_block_matrix_heatmap(P_cond_coarse_truth[['3_6']], 'Coarse Truth Theory'),
                              visualize_pm_block_matrix_heatmap(P_cond_X_truth[['3_6']], 'Truth X'),
                              visualize_pm_block_matrix_heatmap(P_cond_X_coarse_truth[['3_6']], 'Coarse Truth X'),
                              layout_matrix = arr_mat_8)     
    
    
    # 10.2) visualize the entire pm x pm block
    
    P_cond_g_ij_full           <- assemble_block_matrix_v2(P_cond_g_ij, p, m_est)
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
                              textGrob("10. Conditional\nPrecision Operator\n All Processes", gp = gpar(fontsize = 14)),
                              visualize_pm_block_matrix_heatmap(P_cond_est_full, 'Estimate'),
                              visualize_pm_block_matrix_heatmap(P_cond_g_ij_full, 'Estimate G_ij'),
                              visualize_pm_block_matrix_heatmap(P_cond_truth_full, 'Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(P_cond_coarse_truth_full, 'Coarse Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(P_cond_X_truth_full, 'Truth X'), 
                              visualize_pm_block_matrix_heatmap(P_cond_X_coarse_truth_full, 'Coarse Truth X'), 
                              layout_matrix = arr_mat_8)     
    
    
    # part 11
    #graph_threshold <- select_threshold_by_stability(P_cond, p)
    
    final_graph_estimates                  <- estimate_graph(P_cond_est,            C_cond_est,            V_cond_est, p)
    final_graph_estimates_g_ij             <- estimate_graph(P_cond_g_ij,           C_cond_g_ij,           V_cond_g_ij, p)
    final_graph_estimates_X_coarse_truth   <- estimate_graph(P_cond_X_coarse_truth, C_cond_X_coarse_truth, V_cond_X_coarse_truth, p)
    final_graph_estimates_X_truth          <- estimate_graph(P_cond_X_truth,        C_cond_X_truth,        V_cond_X_truth, p)    
    final_graph_estimates_coarse_truth     <- estimate_graph(P_cond_coarse_truth,   C_cond_coarse_truth,   V_cond_coarse_truth, p)
    final_graph_estimates_truth            <- estimate_graph(P_cond_truth,          C_cond_truth,          V_cond_truth, p)
    
    g_est_111 <- grid.arrange(visualize_pm_block_matrix_heatmap(adj_mat_i, 'Ground Truth'), 
                              textGrob("11. Final Estimates", gp = gpar(fontsize = 14)),
                              visualize_pm_block_matrix_heatmap(final_graph_estimates$w_mat, 'Estimate'),
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_g_ij$w_mat, 'Estimate G_ij'),
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_coarse_truth$w_mat, 'Coarse Truth X'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_truth$w_mat, 'Truth X'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_coarse_truth$w_mat, 'Coarse Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_truth$w_mat, 'Truth Theory'), 
                              layout_matrix = arr_mat_8)  
    
    
    
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
    w_mat_g_ij                <- hilbert_schmidt_norm_pm(P_cond_g_ij_full, p, m_est)
    
    

    g_est_112 <- grid.arrange(visualize_pm_block_matrix_heatmap(w_mat_ground_truth, 'Ground Truth'), 
                              visualize_pm_block_matrix_heatmap(w_mat_coarse_ground_truth, 'Coarse Ground Truth'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_truth$w_mat, 'Truth Theory'),   
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_coarse_truth$w_mat, 'Coarse Truth Theory'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_truth$w_mat, 'Truth X'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_X_coarse_truth$w_mat, 'Coarse Truth X'), 
                              visualize_pm_block_matrix_heatmap(final_graph_estimates$w_mat, 'Estimate'),
                              visualize_pm_block_matrix_heatmap(final_graph_estimates_g_ij$w_mat, 'Estimate G_ij'),
                              textGrob("11. Hilbert Schmidt\n Norm", gp = gpar(fontsize = 14)),                              
                              layout_matrix = arr_mat_10)     
    
    
    roc_ground_truth        <- roc_with_threshold(w_mat_ground_truth,        adj_mat_i, 'Ground Truth')
    roc_coarse_ground_truth <- roc_with_threshold(w_mat_coarse_ground_truth, adj_mat_i, 'Coarse Ground Truth')
    roc_truth               <- roc_with_threshold(w_mat_truth,               adj_mat_i, 'Truth Theory')   
    roc_coarse_truth        <- roc_with_threshold(w_mat_coarse_truth,        adj_mat_i, 'Coarse Truth Theory')
    roc_X_truth             <- roc_with_threshold(w_mat_X_truth,             adj_mat_i, 'Truth X')
    roc_X_coarse_truth      <- roc_with_threshold(w_mat_X_coarse_truth,      adj_mat_i, 'Coarse Truth X')
    roc_est                 <- roc_with_threshold(w_mat_est,                 adj_mat_i, 'Estimate')
    roc_g_ij                <- roc_with_threshold(w_mat_g_ij,                adj_mat_i, 'Estimate G_ij')
    
    g_est_113 <- grid.arrange(roc_ground_truth$plot,
                              roc_coarse_ground_truth$plot,
                              roc_truth$plot,
                              roc_coarse_truth$plot,
                              roc_X_truth$plot,
                              roc_X_coarse_truth$plot,
                              roc_est$plot,
                              roc_g_ij$plot,
                              textGrob("11. ROC Curve", gp = gpar(fontsize = 14)),
                              layout_matrix = arr_mat_10)
    
    
    
    
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
    estimated_graphs_part_1 <- list(g_01 = g_01, g_22 = g_est_22)
    
    return(list(estimated_graphs_part_1 = estimated_graphs_part_1,      
                estimated_graphs_part_2 = estimated_graphs_v2))      
    
  } else{
    estimated_graphs_part_1 <- list(g_01 = g_01, g_22 = g_est_22,
                                    g_31 = g_est_31,
                                    g_41 = g_est_41, g_42 = g_est_42,
                                    g_51 = g_est_51, g_52 = g_est_52, g_53 = g_est_53,
                                    g_71 = g_est_71, g_72 = g_est_72)
    
    return(list(estimated_graphs_part_1 = estimated_graphs_part_1,      
                estimated_graphs_part_2 = estimated_graphs_v2))        
    
  }
  
  
  
}
