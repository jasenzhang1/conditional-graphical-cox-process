

visualize_truths_from_est <- function(step_list, graph_ids, time_grid_est, time_grid, time_grid_both, p){
  
  
  # load step_list

  step_1 <- step_list$step_1 
  step_2 <- step_list$step_2 
  step_3 <- step_list$step_3 
  step_4 <- step_list$step_4 
  step_5 <- step_list$step_5 
  step_6 <- step_list$step_6 
  step_7 <- step_list$step_7 
  step_8 <- step_list$step_8 
  step_9 <- step_list$step_9 
  step_10 <- step_list$step_10 
  step_11 <- step_list$step_11 
  step_12 <- step_list$step_12 
  step_2b <- step_list$step_2b 
  
  # graph_results = output
  # graph_ids = vector of strings of graphs we want
  
  graphs = list()
  
  arr_mat <- matrix(1:6, nrow = 2, byrow = F)
  arr_mat_8 <- matrix(1:8, nrow = 2, byrow = F)
  arr_mat_10 <- matrix(1:10, nrow = 2, byrow = F)  
  
  m <- length(time_grid)
  m_est <- length(time_grid_est)
  
  if('11' %in% graph_ids){
    #   - X_k_est                   (p x m_est   x n)
    #   - X_k_truth                 (p x m_truth x n)
    #   - X_k_coarse_truth          (p x m_est   x n)
    #   - X_k_both_truth            (p x m_both  x n)
    #   - Lambda_k_truth            (p x m_truth x n)
    #   - Lambda_k_coarse_truth     (p x m_est   x n)
    graphs[['g_11']] <- grid.arrange(visualize_log_intensity(step_1[[1]][1:5,,1],   time_grid_est,  'Estimate'),
                                     visualize_log_intensity(step_1[[3]][1:5,,1],   time_grid_est,  'Coarser Truth'),
                                     visualize_log_intensity(step_1[[2]][1:5,,1],   time_grid,      'Finer Truth'),
                                     visualize_log_intensity(step_1[[4]][1:5,,1],   time_grid_both, 'Combined Truth'),
                                     textGrob("0. Log Intensity\n of first 5 processes\nof subject 1", gp = gpar(fontsize = 14)),
                                     layout_matrix = arr_mat) 
  }
  
  if('22' %in% graph_ids){
    #   - rho_i_truth               (p x m)
    #   - rho_i_coarse_truth        (p x m_est)
    #   - rho_i_X_truth             (p x m)
    #   - rho_i_X_coarse_truth      (p x m_est)
    #   - rho_i_est                 (p x m_est)
    #   - rho_list                  (list format, [[1]] = rho_i, [[2]] = rho_ij, rho_ij may be rho_ii only)  
    
    graphs[['g_22']] <- grid.arrange(visualize_log_intensity(step_2[[1]][1:5,],   time_grid,     'Truth Theory'),
                                     visualize_log_intensity(step_2[[2]][1:5,],   time_grid_est, 'Coarse Truth Theory'),
                                     visualize_log_intensity(step_2[[3]][1:5,],   time_grid,     'Truth X'),
                                     visualize_log_intensity(step_2[[4]][1:5,],   time_grid_est, 'Coarse Truth X'),
                                     textGrob("2. Rho_i\nEstimation for \n first 5 processes", gp = gpar(fontsize = 14)),
                                     visualize_log_intensity(step_2[[5]][1:5,],   time_grid_est, 'Estimate'),
                                     layout_matrix = arr_mat) 
  }
  
  if('24' %in% graph_ids){
    #   - rho_ii_truth               (list of m x m matrices)
    #   - rho_ii_coarse_truth        (list of m_est x m_est matrices)
    #   - rho_ii_X_truth             (list of m x m matrices)
    #   - rho_ii_X_coarse_truth      (list of m_est x m_est matrices)
    #   - rho_ii_est                 (list of m_est x m_est matrices)
    
    graphs[['g_24']] <- grid.arrange(visualize_matrix_heatmap(step_2b[[1]][['1_1']],   'Truth Theory',        40000, NULL, 200000),
                                     visualize_matrix_heatmap(step_2b[[2]][['1_1']],   'Coarse Truth Theory', 40000, NULL, 200000),
                                     visualize_matrix_heatmap(step_2b[[3]][['1_1']],   'Truth X',             40000, NULL, 200000),
                                     visualize_matrix_heatmap(step_2b[[4]][['1_1']],   'Coarse Truth X',      40000, NULL, 200000),
                                     textGrob("2. Rho ii\nEstimation \n for process 1", gp = gpar(fontsize = 14)),
                                     visualize_matrix_heatmap(step_2b[[5]][['1_1']],   'Estimate',            40000, NULL, 200000),
                                     layout_matrix = arr_mat)
  }
  
  if('25' %in% graph_ids){
    graphs[['g_25']] <- grid.arrange(visualize_matrix_heatmap(step_2b[[1]][['1_2']],   'Truth Theory',        40000, NULL, 100000),
                                     visualize_matrix_heatmap(step_2b[[2]][['1_2']],   'Coarse Truth Theory', 40000, NULL, 100000),
                                     visualize_matrix_heatmap(step_2b[[3]][['1_2']],   'Truth X',             40000, NULL, 100000),
                                     visualize_matrix_heatmap(step_2b[[4]][['1_2']],   'Coarse Truth X',      40000, NULL, 100000),
                                     textGrob("2. Rho ij\nEstimation \n for process 1 and 2", gp = gpar(fontsize = 14)),
                                     visualize_matrix_heatmap(step_2b[[5]][['1_2']],   'Estimate',            40000, NULL, 100000),
                                     layout_matrix = arr_mat)    
  }
  
  if('26' %in% graph_ids){
    graphs[['g_26']] <- grid.arrange(visualize_matrix_heatmap(step_2b[[1]][['3_6']],   'Truth Theory',        40000, NULL, 120000),
                                     visualize_matrix_heatmap(step_2b[[2]][['3_6']],   'Coarse Truth Theory', 40000, NULL, 120000),
                                     visualize_matrix_heatmap(step_2b[[3]][['3_6']],   'Truth X',             40000, NULL, 120000),
                                     visualize_matrix_heatmap(step_2b[[4]][['3_6']],   'Coarse Truth X',      40000, NULL, 120000),
                                     textGrob("2. Rho ij\nEstimation \n for process 3 and 6", gp = gpar(fontsize = 14)),
                                     visualize_matrix_heatmap(step_2b[[5]][['3_6']],   'Estimate',            40000, NULL, 120000),
                                     layout_matrix = arr_mat)     
  }
  
  if('31' %in% graph_ids){
    #   - g_ij_ground_truth                (list of m x m matrices for i_j entries)
    #   - g_ij_coarse_ground_truth         (list of m_est x m_est matrices for i_j entries)
    #   - g_ij_truth                       (list of m x m matrices for i_j entries)
    #   - g_ij_coarse_truth                (list of m_est x m_est matrices for i_j entries)
    #   - g_ij_X_truth                     (list of m x m matrices for i_j entries)
    #   - g_ij_X_coarse_truth              (list of m_est x m_est matrices for i_j entries)
    #   - g_ij_est                         (list of m_est x m_est matrices for i_j entries)      
    graphs[['g_31']] <- grid.arrange(visualize_matrix_heatmap(step_3[[1]][['1_1']], 'Ground Truth', -1, 0, 1),
                                     visualize_matrix_heatmap(step_3[[2]][['1_1']], 'Coarse Ground Truth', -1, 0, 1),
                                     visualize_matrix_heatmap(step_3[[3]][['1_1']], 'Truth Theory', -1, 0, 1),  
                                     visualize_matrix_heatmap(step_3[[4]][['1_1']], 'Coarse Truth Theory', -1, 0, 1),
                                     visualize_matrix_heatmap(step_3[[5]][['1_1']], 'Truth X', -1, 0, 1),    
                                     visualize_matrix_heatmap(step_3[[6]][['1_1']], 'Coarse Truth X', -1, 0, 1),
                                     textGrob("3. Covariance Function\nEstimation (G_ii)", gp = gpar(fontsize = 14)),
                                     visualize_matrix_heatmap(step_3[[7]][['1_1']], 'Estimate', -1, 0, 1),
                                     layout_matrix = arr_mat_8)     
  }
  

  # 4.1) do the eigenfunctions between coarse and regular truths differ? - yes - the g_ij is too different from the truth
  
  if('41' %in% graph_ids){
    #   - eigen_decomp_truth               (list of 3 things)
    #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
    #     - [[2]] eigenvectors (list of p matrices of m x d_i)
    #     - [[3]] n_dims       (list of p integers denoting d_i)
    #
    #   - eigen_decomp_coarse_truth
    #   - eigen_decomp_X_truth
    #   - eigen_decomp_X_coarse_truth
    #   - eigen_decomp_est    
    graphs[['g_41']] <- grid.arrange(visualize_log_intensity(t(step_4[[1]]$eigenfunctions[[1]]),    time_grid,      'Truth Theory'),
                                     visualize_log_intensity(t(step_4[[2]]$eigenfunctions[[1]]),    time_grid_est,  'Coarse Truth Theory'),
                                     visualize_log_intensity(t(step_4[[3]]$eigenfunctions[[1]]),    time_grid,      'Truth X'),
                                     visualize_log_intensity(t(step_4[[4]]$eigenfunctions[[1]]),    time_grid_est,  'Coarse Truth X'),
                                     textGrob("4. Eigenfunction\nBasis of\nProcess 1", gp = gpar(fontsize = 14)),
                                     visualize_log_intensity(t(step_4[[5]]$eigenfunctions[[1]]),    time_grid_est,  'Estimate'),
                                     layout_matrix = arr_mat)    
  }

  # visualize g_ij vs its eigendecomposition reconstruction
  
  if('42' %in% graph_ids){
    
    g_ii_truth          <- prep_eigendecomposition_ii(step_3[[3]], p)
    g_ii_coarse_truth   <- prep_eigendecomposition_ii(step_3[[4]], p)
    g_ii_X_truth        <- prep_eigendecomposition_ii(step_3[[5]], p)
    g_ii_X_coarse_truth <- prep_eigendecomposition_ii(step_3[[6]], p)
    g_ii_est            <- prep_eigendecomposition_ii(step_3[[7]], p) 
    
    # validate
    g_ii_truth_decomp           <- validate_eigendecomposition_ii(g_ii_truth,          step_4[[1]])
    g_ii_coarse_truth_decomp    <- validate_eigendecomposition_ii(g_ii_coarse_truth,   step_4[[2]]) 
    g_ii_X_truth_decomp         <- validate_eigendecomposition_ii(g_ii_X_truth,        step_4[[3]]) 
    g_ii_X_coarse_truth_decomp  <- validate_eigendecomposition_ii(g_ii_X_coarse_truth, step_4[[4]])  
    g_ii_est_decomp             <- validate_eigendecomposition_ii(g_ii_est,            step_4[[5]])
    
    
    
    graphs[['g_42']] <- grid.arrange(visualize_matrix_heatmap(g_ii_truth[,,1],               'Truth Theory (TT)', -1, 0, 1),  # ground truth 50 x 50 covariance
                                     visualize_matrix_heatmap(g_ii_truth_decomp[,,1],        'TT Reconstruct', -1,0, 1), # ground truth reconstructed
                                     visualize_matrix_heatmap(g_ii_coarse_truth[,,1],        'Coarse Truth Theory (CTT)', -1, 0, 1), # ground truth 19 x 19 covariance
                                     visualize_matrix_heatmap(g_ii_coarse_truth_decomp[,,1], 'CTT Reconstruct', -1, 0, 1), # ground truth 19 x 19 covariance
                                     visualize_matrix_heatmap(g_ii_est[,,1],                 'Estimate (E)', -1, 0, 1),    # estimate 19 x 19 covariance
                                     visualize_matrix_heatmap(g_ii_est_decomp[,,1],          'E Reconstruct', -1, 0, 1), # reconstructed estimate
                                     textGrob("4. Eigenfunction\nReconstruction", gp = gpar(fontsize = 14)),
                                     layout_matrix = arr_mat_8)
  }
  
  # checking eigenreconstruction
  if('43' %in% graph_ids){
    
    eigen_decomp_est <- step_4[[1]]
    eigen_decomp_est <- step_4[[1]]
    eigen_decomp_est <- step_4[[1]]
    eigen_decomp_est <- step_4[[1]]
    eigen_decomp_est <- step_4[[5]]
    
    # 4.2) are all eigenfunctions orthogonal? Orthonormal? 

    
    mat_truth          <- t(step_4[[1]]$eigenfunctions[[1]]) %*% step_4[[1]]$eigenfunctions[[1]]
    mat_coarse_truth   <- t(step_4[[2]]$eigenfunctions[[1]]) %*% step_4[[2]]$eigenfunctions[[1]]   
    mat_X_truth        <- t(step_4[[3]]$eigenfunctions[[1]]) %*% step_4[[3]]$eigenfunctions[[1]]
    mat_X_coarse_truth <- t(step_4[[4]]$eigenfunctions[[1]]) %*% step_4[[4]]$eigenfunctions[[1]]
    mat_est            <- t(step_4[[5]]$eigenfunctions[[1]]) %*% step_4[[5]]$eigenfunctions[[1]]
    
    graphs[['g_43']] <- grid.arrange(visualize_matrix_heatmap(mat_truth,   'Truth Theory', zmid = 0),
                                     visualize_matrix_heatmap(mat_truth,   'Coarse Truth Theory', zmid = 0),
                                     visualize_matrix_heatmap(mat_truth,   'Truth X', zmid = 0),
                                     visualize_matrix_heatmap(mat_truth,   'Coarse Truth X', zmid = 0),
                                     textGrob("4. Eigenfunction\nOrthogonality", gp = gpar(fontsize = 14)),
                                     visualize_matrix_heatmap(mat_truth,   'Estimate', zmid = 0),
                                     layout_matrix = arr_mat)
  }
  
  # reconstruction error
  if('44' %in% graph_ids){
    
    g_ii_truth          <- prep_eigendecomposition_ii(step_3[[3]], p)
    g_ii_coarse_truth   <- prep_eigendecomposition_ii(step_3[[4]], p)
    g_ii_X_truth        <- prep_eigendecomposition_ii(step_3[[5]], p)
    g_ii_X_coarse_truth <- prep_eigendecomposition_ii(step_3[[6]], p)
    g_ii_est            <- prep_eigendecomposition_ii(step_3[[7]], p) 
    
    # validate
    g_ii_truth_decomp           <- validate_eigendecomposition_ii(g_ii_truth,          step_4[[1]])
    g_ii_coarse_truth_decomp    <- validate_eigendecomposition_ii(g_ii_coarse_truth,   step_4[[2]]) 
    g_ii_X_truth_decomp         <- validate_eigendecomposition_ii(g_ii_X_truth,        step_4[[3]]) 
    g_ii_X_coarse_truth_decomp  <- validate_eigendecomposition_ii(g_ii_X_coarse_truth, step_4[[4]])  
    g_ii_est_decomp             <- validate_eigendecomposition_ii(g_ii_est,            step_4[[5]]) 
    
    
    graphs[['g_44']] <- grid.arrange(visualize_error_histogram(g_ii_truth, g_ii_truth_decomp, 'Truth Theory', 20),
                                     visualize_error_histogram(g_ii_coarse_truth, g_ii_coarse_truth_decomp, 'Coarse Truth Theory', 20),
                                     visualize_error_histogram(g_ii_X_truth, g_ii_X_truth_decomp, 'Truth X', 20),
                                     visualize_error_histogram(g_ii_X_coarse_truth, g_ii_X_coarse_truth_decomp, 'Coarse Truth X', 20),
                                     textGrob("4. Eigendecomposition\nReconstruction Error", gp = gpar(fontsize = 14)),
                                     visualize_error_histogram(g_ii_est, g_ii_est_decomp, 'Estimate', 20),
                                     layout_matrix = arr_mat)  
  

  }
  

  
  
  # 5.3) Graph X_k (true, true_reconstruct, estimate, estimate_reconstruct)  for the first subject
  
  if (any(c('51', '52', '53') %in% graph_ids)) {
    #   - X_k_est                   (p x m_est   x n)
    #   - X_k_truth                 (p x m_truth x n)
    #   - X_k_coarse_truth          (p x m_est   x n)
    #   - X_k_both_truth            (p x m_both  x n)
    #   - Lambda_k_truth            (p x m_truth x n)
    #   - Lambda_k_coarse_truth     (p x m_est   x n)    
    


    df1 <- validate_kl_full(step_1[[1]], step_4[[5]]$eigenfunctions, time_grid_est,  'Estimate',            ncores)
    df2 <- validate_kl_full(step_1[[3]], step_4[[4]]$eigenfunctions, time_grid_est,  'Coarse Truth X',      ncores)
    df3 <- validate_kl_full(step_1[[2]], step_4[[3]]$eigenfunctions, time_grid,      'Truth X',             ncores)
    df4 <- validate_kl_full(step_1[[3]], step_4[[2]]$eigenfunctions, time_grid_est,  'Coarse Truth Theory', ncores)
    df5 <- validate_kl_full(step_1[[2]], step_4[[1]]$eigenfunctions, time_grid,      'Truth Theory',        ncores)

    df_both <- rbind(df1, df5)
    df_both$cat <- factor(df_both$cat)

    df_all <- rbind(df1, df2, df3, df4, df5)
    df_all$cat <- factor(df_all$cat)

    df_no_reconstruct <- df_all %>% filter(cat %in% c('Estimate', 'Coarse Truth X', 'Truth X', 'Coarse Truth Theory', 'Truth Theory'))


    # estimate and truth theory + reconstruct
    graphs[['g_51']] <- ggplot() +
      geom_line(data = df_both, aes(x = Var2, y = value, color = cat)) +
      facet_wrap(~ Var1) +
      ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
      ylab('Log Intensity') +
      xlab('Time') +
      labs(color = "Truth or Estimate") +
      theme_bw()

    # estimate, truth x and truth theory + reconstruct
    graphs[['g_52']] <- ggplot() +
      geom_line(data = df_all, aes(x = Var2, y = value, color = cat)) +
      facet_wrap(~ Var1) +
      ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
      ylab('Log Intensity') +
      xlab('Time') +
      labs(color = "Truth or Estimate") +
      theme_bw()

    # estimate, truth x and truth theory + no reconstruct
    graphs[['g_53']] <- ggplot() +
      geom_line(data = df_no_reconstruct, aes(x = Var2, y = value, color = cat)) +
      facet_wrap(~ Var1) +
      ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
      ylab('Log Intensity') +
      xlab('Time') +
      labs(color = "Truth or Estimate") +
      theme_bw()
  }
  
  
  # for pair 1_2, plot estimate and truths - this one should look correlated
  
  if('81' %in% graph_ids){
    #   - V_cond_ground_truth_full          (pm x pm matrix)
    #   - V_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
    #   - V_cond_truth_full                 (pm x pm matrix)
    #   - V_cond_coarse_truth_full          (pm_est x pm_est matrix)
    #   - V_cond_X_truth_full               (pm x pm matrix)
    #   - V_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
    #   - V_cond_est_full                   (pm_est x pm_est matrix)    
    graphs[['g_81']] <- grid.arrange(visualize_matrix_heatmap(extract_block_structure_ij(step_8[[1]], m,     1, 2), 'Ground Truth'), 
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[2]], m_est, 1, 2), 'Coarse Ground Truth'), 
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[3]], m,     1, 2), 'Truth Theory'),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[4]], m_est, 1, 2), 'Coarse Truth Theory'),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[5]], m,     1, 2), 'Truth X'), 
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[6]], m_est, 1, 2), 'Coarse Truth X'), 
                                     textGrob("8. Conditional\nCovariance Operator\nProcess 1 with 2", gp = gpar(fontsize = 14)),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[7]], m_est, 1, 2), 'Estimate'), 
                                     layout_matrix = arr_mat_8)      
  }
  
  # visualize again for process 3 with 6 - they should not look correlated  
  if('82' %in% graph_ids){
    graphs[['g_82']] <- grid.arrange(visualize_matrix_heatmap(extract_block_structure_ij(step_8[[1]], m,     3, 6), 'Ground Truth'), 
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[2]], m_est, 3, 6), 'Coarse Ground Truth'), 
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[3]], m,     3, 6), 'Truth Theory'),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[4]], m_est, 3, 6), 'Coarse Truth Theory'),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[5]], m,     3, 6), 'Truth X'), 
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[6]], m_est, 3, 6), 'Coarse Truth X'), 
                                     textGrob("8. Conditional\nCovariance Operator\nProcess 3 with 6", gp = gpar(fontsize = 14)),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_8[[7]], m_est, 3, 6), 'Estimate'), 
                                     layout_matrix = arr_mat_8)      
  }
  

 
  
  # est vs truth vs ground truth
  if('83' %in% graph_ids){
    graphs[['g_83']] <- grid.arrange(visualize_pm_block_matrix_heatmap(step_8[[1]], 'Ground Truth'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[3]], 'Truth Theory'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[7]], 'Estimate'),
                                     nrow = 1)    
  }

  
  # all 5 + ground truth
  if('84' %in% graph_ids){
    graphs[['g_84']] <- grid.arrange(visualize_pm_block_matrix_heatmap(step_8[[1]], 'Ground Truth'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[2]], 'Coarse Ground Truth'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[3]], 'Truth Theory'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[4]], 'Coarse Truth Theory'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[5]], 'Truth X'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[6]], 'Coarse Truth X'), 
                                     textGrob("8. Conditional\nCovariance Operator\nAll Processes", gp = gpar(fontsize = 14)),
                                     visualize_pm_block_matrix_heatmap(step_8[[7]], 'Estimate'),
                                     layout_matrix = arr_mat_8)      
  }


  
  if('91' %in% graph_ids){
    #   - C_cond_ground_truth_full          (pm x pm matrix)
    #   - C_cond_corase_ground_truth_full   (pm_est x pm_est matrix)
    #   - C_cond_truth_full                 (pm x pm matrix)
    #   - C_cond_coarse_truth_full          (pm_est x pm_est matrix)
    #   - C_cond_X_truth_full               (pm x pm matrix)
    #   - C_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
    #   - C_cond_est_full                   (pm_est x pm_est matrix)      
    graphs[['g_91']] <- grid.arrange(visualize_matrix_heatmap(extract_block_structure_ij(step_9[[1]], m,     1, 2), 'Ground Truth', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[2]], m_est, 1, 2), 'Coarse Ground Truth', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[3]], m,     1, 2), 'Truth Theory', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[4]], m_est, 1, 2), 'Coarse Truth Theory', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[5]], m,     1, 2), 'Truth X', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[6]], m_est, 1, 2), 'Coarse Truth X', zmid = 0),
                                     textGrob("9. Conditional\nCorrelation Operator\nProcess 1 with 2", gp = gpar(fontsize = 14)),  
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[7]], m_est, 1, 2), 'Estimate', zmid = 0),
                                     layout_matrix = arr_mat_8)      
  }

  if('92' %in% graph_ids){
    graphs[['g_92']] <- grid.arrange(visualize_matrix_heatmap(extract_block_structure_ij(step_9[[1]], m,     3, 6), 'Ground Truth', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[2]], m_est, 3, 6), 'Coarse Ground Truth', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[3]], m,     3, 6), 'Truth Theory', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[4]], m_est, 3, 6), 'Coarse Truth Theory', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[5]], m,     3, 6), 'Truth X', zmid = 0),
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[6]], m_est, 3, 6), 'Coarse Truth X', zmid = 0),
                                     textGrob("9. Conditional\nCorrelation Operator\nProcess 3 with 6", gp = gpar(fontsize = 14)),  
                                     visualize_matrix_heatmap(extract_block_structure_ij(step_9[[7]], m_est, 3, 6), 'Estimate', zmid = 0),
                                     layout_matrix = arr_mat_8)     
  }
  
  if('93' %in% graph_ids){
    graphs[['g_93']] <- grid.arrange(visualize_pm_block_matrix_heatmap(step_9[[1]], 'Ground Truth'), 
                                     visualize_pm_block_matrix_heatmap(step_9[[3]], 'Truth Theory'), 
                                     visualize_pm_block_matrix_heatmap(step_9[[7]], 'Estimate'),
                                     nrow = 1)      
  }
  
  
  
  # all 5 + ground truth
  
  if('94' %in% graph_ids){
    graphs[['g_94']] <- grid.arrange(visualize_pm_block_matrix_heatmap(step_9[[1]], 'Ground Truth'), 
                                     visualize_pm_block_matrix_heatmap(step_9[[2]], 'Coarse Ground Truth'), 
                                     visualize_pm_block_matrix_heatmap(step_9[[3]], 'Truth Theory'), 
                                     visualize_pm_block_matrix_heatmap(step_9[[4]], 'Coarse Truth Theory'), 
                                     visualize_pm_block_matrix_heatmap(step_9[[5]], 'Truth X'), 
                                     visualize_pm_block_matrix_heatmap(step_9[[6]], 'Coarse Truth X'), 
                                     textGrob("9. Conditional\nCorrelation Operator\nAll Processes", gp = gpar(fontsize = 14)),
                                     visualize_pm_block_matrix_heatmap(step_9[[7]], 'Estimate'),
                                     layout_matrix = arr_mat_8)      
  }

  
  # sanity check - hilbert schmidt norms of the C_Cond
  if('95' %in% graph_ids){
    
    C_HS_ground_truth_full          <- hilbert_schmidt_norm_pm(step_9[[1]], p, m) 
    C_HS_coarse_ground_truth_full   <- hilbert_schmidt_norm_pm(step_9[[2]], p, m_est)
    C_HS_truth_full                 <- hilbert_schmidt_norm_pm(step_9[[3]], p, m)
    C_HS_coarse_truth_full          <- hilbert_schmidt_norm_pm(step_9[[4]], p, m_est)
    C_HS_X_truth_full               <- hilbert_schmidt_norm_pm(step_9[[5]], p, m)
    C_HS_X_coarse_truth_full        <- hilbert_schmidt_norm_pm(step_9[[6]], p, m_est)
    C_HS_est                        <- hilbert_schmidt_norm_pm(step_9[[7]], p, m_est)
    
    diag(C_HS_ground_truth_full) <- 0
    diag(C_HS_coarse_ground_truth_full) <- 0
    diag(C_HS_truth_full) <- 0
    diag(C_HS_coarse_truth_full) <- 0
    diag(C_HS_X_truth_full) <- 0
    diag(C_HS_X_coarse_truth_full) <- 0
    diag(C_HS_est) <- 0
    
    graphs[['g_95']] <- grid.arrange(visualize_matrix_heatmap(C_HS_ground_truth_full,        'Ground Truth', zmid = 0),
                                     visualize_matrix_heatmap(C_HS_coarse_ground_truth_full, 'Coarse Ground Truth', zmid = 0),
                                     visualize_matrix_heatmap(C_HS_truth_full,               'Truth Theory', zmid = 0),
                                     visualize_matrix_heatmap(C_HS_coarse_truth_full,        'Coarse Truth Theory', zmid = 0),
                                     visualize_matrix_heatmap(C_HS_X_truth_full,             'Truth X', zmid = 0),
                                     visualize_matrix_heatmap(C_HS_X_coarse_truth_full,      'Coarse Truth X', zmid = 0),
                                     textGrob("9. Conditional\nCorrelation Operator\nHS Norm", gp = gpar(fontsize = 14)),
                                     visualize_matrix_heatmap(C_HS_est,                      'Estimate', zmid = 0),
                                     layout_matrix = arr_mat_8)     
  }
 
  
  if('101' %in% graph_ids){
    #   - P_cond_ground_truth_full          (pm x pm matrix)
    #   - P_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
    #   - P_cond_truth_full                 (pm x pm matrix)
    #   - P_cond_coarse_truth_full          (pm_est x pm_est matrix)
    #   - P_cond_X_truth_full               (pm x pm matrix)
    #   - P_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
    #   - P_cond_est_full                   (pm_est x pm_est matrix)    
    graphs[['g_101']] <- grid.arrange(visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[1]], m,     1, 2), 'Truth Ground'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[2]], m_est, 1, 2), 'Coarse Truth Ground'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[3]], m,     1, 2), 'Truth Theory'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[4]], m_est, 1, 2), 'Coarse Truth Theory'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[5]], m,     1, 2), 'Truth X'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[6]], m_est, 1, 2), 'Coarse Truth X'),
                                      textGrob("10. Conditional\nPrecision Operator\n Process 1 and 2", gp = gpar(fontsize = 14)),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[7]], m_est, 1, 2), 'Estimate'),
                                      layout_matrix = arr_mat_8)      
  }

  if('102' %in% graph_ids){
    graphs[['g_102']] <- grid.arrange(visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[1]], m,     3, 6), 'Truth Ground'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[2]], m_est, 3, 6), 'Coarse Truth Ground'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[3]], m,     3, 6), 'Truth Theory'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[4]], m_est, 3, 6), 'Coarse Truth Theory'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[5]], m,     3, 6), 'Truth X'),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[6]], m_est, 3, 6), 'Coarse Truth X'),
                                      textGrob("10. Conditional\nPrecision Operator\n Process 3 and 6", gp = gpar(fontsize = 14)),
                                      visualize_pm_block_matrix_heatmap(extract_block_structure_ij(step_10[[7]], m_est, 3, 6), 'Estimate'),
                                      layout_matrix = arr_mat_8)      
  }
  

  # 10.2) visualize the entire pm x pm block
  
  if('103' %in% graph_ids){
    graphs[['g_103']] <- grid.arrange(visualize_pm_block_matrix_heatmap(step_10[[1]], 'Ground Truth'), 
                                      visualize_pm_block_matrix_heatmap(step_10[[3]], 'Truth Theory'), 
                                      visualize_pm_block_matrix_heatmap(step_10[[7]], 'Estimate'),
                                      nrow = 1)      
  }
      
  
  
  # all 5 + ground truth
  if('104' %in% graph_ids){
    graphs[['g_104']] <- grid.arrange(visualize_pm_block_matrix_heatmap(step_10[[1]], 'Ground Truth'), 
                                      visualize_pm_block_matrix_heatmap(step_10[[2]], 'Coarse Ground Truth'), 
                                      visualize_pm_block_matrix_heatmap(step_10[[3]], 'Truth Theory'), 
                                      visualize_pm_block_matrix_heatmap(step_10[[4]], 'Coarse Truth Theory'), 
                                      visualize_pm_block_matrix_heatmap(step_10[[5]], 'Truth X'), 
                                      visualize_pm_block_matrix_heatmap(step_10[[6]], 'Coarse Truth X'), 
                                      textGrob("10. Conditional\nPrecision Operator\n All Processes", gp = gpar(fontsize = 14)),
                                      visualize_pm_block_matrix_heatmap(step_10[[7]], 'Estimate'),
                                      layout_matrix = arr_mat_8)     
  }
 

  # visualize HS norms with adj_mat (0 or 1)
  # g_111 <- grid.arrange(visualize_pm_block_matrix_heatmap(adj_mat_i, 'Adjacency Truth'), 
  #                       visualize_pm_block_matrix_heatmap(adj_mat_i, 'Coarse Adjacency Truth'), 
  #                       visualize_pm_block_matrix_heatmap(w_mat_truth, 'Truth Theory'), 
  #                       visualize_pm_block_matrix_heatmap(w_mat_coarse_truth, 'Coarse Truth Theory'), 
  #                       visualize_pm_block_matrix_heatmap(w_mat_X_truth, 'Truth X'), 
  #                       visualize_pm_block_matrix_heatmap(w_mat_X_coarse_truth, 'Coarse Truth X'), 
  #                       textGrob("11. Final Estimates\nvs Adj Truth", gp = gpar(fontsize = 14)),
  #                       visualize_pm_block_matrix_heatmap(w_mat_est, 'Estimate'),
  #                       layout_matrix = arr_mat_8)  
  
  
  
  # g_111b <- grid.arrange(visualize_pm_block_matrix_heatmap(adj_mat_i, 'Coarse Ground Truth'),
  #                        visualize_pm_block_matrix_heatmap(adj_mat_i, 'Ground Truth'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_coarse_truth, 'Coarse Truth Theory'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_truth, 'Truth Theory'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_X_coarse_truth, 'Coarse Truth X'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_X_truth, 'Truth X'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_est, 'Estimate'),
  #                        textGrob("11b. Normalized\nFinal Estimates\nvs Adj Truth", gp = gpar(fontsize = 14)),
  #                        layout_matrix = arr_mat_8)      
  
  
  # visualize HS norms with w_mat ground truth
  
  if('112' %in% graph_ids){
    #   - w_mat_ground_truth        (p x p)   matrix of HS norms of the pm x pm ground truth
    #   - w_mat_coarse_ground_truth (p x p)   matrix of HS norms of the pm_est x pm_est ground truth 
    #   - w_mat_X_coarse_truth      (p x p)   matrix of HS norms of ...
    #   - w_mat_X_truth             (p x p)   matrix of HS norms of ...
    #   - w_mat_coarse_truth        (p x p)   matrix of HS norms of ...
    #   - w_mat_truth               (p x p)   matrix of HS norms of ...
    #   - w_mat_est                 (p x p)   matrix of HS norms of ...    
    graphs[['g_112']] <- grid.arrange(visualize_pm_block_matrix_heatmap(step_11[[1]], 'Ground Truth'), 
                                      visualize_pm_block_matrix_heatmap(step_11[[2]], 'Coarse Ground Truth'), 
                                      visualize_pm_block_matrix_heatmap(step_11[[3]], 'Truth Theory'),   
                                      visualize_pm_block_matrix_heatmap(step_11[[4]], 'Coarse Truth Theory'), 
                                      visualize_pm_block_matrix_heatmap(step_11[[5]], 'Truth X'), 
                                      visualize_pm_block_matrix_heatmap(step_11[[6]], 'Coarse Truth X'),
                                      textGrob("11. Hilbert Schmidt\n Norm", gp = gpar(fontsize = 14)),   
                                      visualize_pm_block_matrix_heatmap(step_11[[7]], 'Estimate'),
                                      layout_matrix = arr_mat_8)     
  }
  
  # g_112b <- grid.arrange(visualize_pm_block_matrix_heatmap(w_mat_normalized_coarse_ground_truth, 'Coarse Ground Truth'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_ground_truth, 'Ground Truth'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_coarse_truth, 'Coarse Truth Theory'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_truth, 'Truth Theory'),   
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_X_coarse_truth, 'Coarse Truth X'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_X_truth, 'Truth X'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_est, 'Estimate'),
  #                        textGrob("11b. Normalized\nHilbert Schmidt\n Norm", gp = gpar(fontsize = 14)),                              
  #                        layout_matrix = arr_mat_8)  
  
  if('113' %in% graph_ids){
    
    #   - roc_ground_truth           (list of roc outputs)
    #   - roc_coarse_ground_truth
    #   - roc_truth
    #   - roc_coarse_truth
    #   - roc_X_truth
    #   - roc_X_coarse_truth
    #   - roc_est
    
    graphs[['g_113']] <- grid.arrange(step_12[[1]]$plot,
                                      step_12[[2]]$plot,
                                      step_12[[3]]$plot,
                                      step_12[[4]]$plot,
                                      step_12[[5]]$plot,
                                      step_12[[6]]$plot,
                                      textGrob("12. ROC Curve", gp = gpar(fontsize = 14)),
                                      step_12[[7]]$plot,
                                      layout_matrix = arr_mat_8)     
  }
  
  return(graphs)
 
  
}
