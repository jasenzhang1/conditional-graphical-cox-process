source('functions/00c_block_matrix_arrange.R')
source('functions/28z_Visualization_helpers.R')
source('functions/28_Simulation_Visualization.R')
source('functions/00b_matrix_norms.R')

library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)


# functions to compartmentalize all the g_XX graphs

# any function that ends in "prep" returns a list of plots ready to be printed
# all other functions store the canvas in g and return(g)

result_01_prep <- function(step_0){
  
  # return a list of plots
  
  # Build plots dynamically
  plots <- lapply(seq_len(length(step_0)), function(i) {
    visualize_matrix_heatmap(step_0[[i]],
                             paste0('y_c = ', as.character(i)), -10, NULL, 10)
  })  
  
  return(plots)
}

result_01 <- function(step_0){
  

  plots <- result_01_prep(step_0)
  
  # Add caption grob
  caption <- textGrob("0. Ground Truth\n Precision wrt Queried\nContinuous Covariate", gp = gpar(fontsize = 14))
  
  # Combine plots + caption
  grobs <- c(plots, list(caption))
  
  g <- grid.arrange(grobs = grobs, nrow = 1)  
  return(g)
}

result_11_prep <- function(step_1, time_grid, time_grid_est){
  #   - X_k_est                   (p x m_est   x n)
  #   - X_k_truth                 (p x m_truth x n)
  #   - X_k_coarse_truth          (p x m_est   x n)
  #   - X_k_both_truth            (p x m_both  x n)
  #   - Lambda_k_truth            (p x m_truth x n)
  #   - Lambda_k_coarse_truth     (p x m_est   x n)
  g_list <- list(visualize_log_intensity(step_1[[1]][1:5,,1],   time_grid_est,  'Estimate'),
                 visualize_log_intensity(step_1[[3]][1:5,,1],   time_grid_est,  'Coarser Truth'),
                 visualize_log_intensity(step_1[[2]][1:5,,1],   time_grid,      'Finer Truth'),
                 visualize_log_intensity(step_1[[4]][1:5,,1],   time_grid_both, 'Combined Truth'))
  
  return(g_list)
}

result_22_prep <- function(step_2, time_grid, time_grid_est, full = T){
  
  if(! full){
    g_list <- list(visualize_log_intensity(step_2$rho_i_est[1:5,],     time_grid_est, 'Estimate'),
                   visualize_log_intensity(step_2$rho_i_truth[1:5,],   time_grid,     'Truth'))
    return(g_list)
  }
  
  #   - rho_i_truth               (p x m)
  #   - rho_i_coarse_truth        (p x m_est)
  #   - rho_i_X_truth             (p x m)
  #   - rho_i_X_coarse_truth      (p x m_est)
  #   - rho_i_est                 (p x m_est)
  #   - rho_list                  (list format, [[1]] = rho_i, [[2]] = rho_ij, rho_ij may be rho_ii only)  
  
  g_list <- list(visualize_log_intensity(step_2[[1]][1:5,],   time_grid,     'Truth Theory'),
                 visualize_log_intensity(step_2[[2]][1:5,],   time_grid_est, 'Coarse Truth Theory'),
                 visualize_log_intensity(step_2[[3]][1:5,],   time_grid,     'Truth X'),
                 visualize_log_intensity(step_2[[4]][1:5,],   time_grid_est, 'Coarse Truth X'),
                 visualize_log_intensity(step_2[[5]][1:5,],   time_grid_est, 'Estimate'))
  
  return(g_list)
}

result_20s_prep <- function(step_2b, i, j, full = T){
  #   - rho_ii_truth               (list of m x m matrices)
  #   - rho_ii_coarse_truth        (list of m_est x m_est matrices)
  #   - rho_ii_X_truth             (list of m x m matrices)
  #   - rho_ii_X_coarse_truth      (list of m_est x m_est matrices)
  #   - rho_ii_est                 (list of m_est x m_est matrices)
  
  key <- paste0(i, '_', j)
  
  if(! full){
    g_list <- list(visualize_matrix_heatmap(step_2b$rho_ii_truth[[key]],   'Truth Theory',        40000, NULL, 200000),
                   visualize_matrix_heatmap(step_2b$rho_ii_est[[key]],     'Estimate',            40000, NULL, 200000))
    return(g_list)
  }
  
  g_list <- list(visualize_matrix_heatmap(step_2b[[1]][[key]],   'Truth Theory',        40000, NULL, 200000),
                 visualize_matrix_heatmap(step_2b[[2]][[key]],   'Coarse Truth Theory', 40000, NULL, 200000),
                 visualize_matrix_heatmap(step_2b[[3]][[key]],   'Truth X',             40000, NULL, 200000),
                 visualize_matrix_heatmap(step_2b[[4]][[key]],   'Coarse Truth X',      40000, NULL, 200000),
                 visualize_matrix_heatmap(step_2b[[5]][[key]],   'Estimate',            40000, NULL, 200000))
  
  return(g_list)
}

result_29 <- function(step_2, y_c_id){
  
  # plot the distribution of the weights in estimating rho_i(t)
  
  df_29 <- data.frame(x = step_2$y_c_s[,1], y = step_2$weights)
  g <- ggplot(data = df_29, aes(x = x, y = y)) + geom_line() + geom_point() + 
    ylab('weight') + 
    xlab('y_c_k') + 
    ggtitle(paste0('y_c_query=', y_c_id))
  
  return(g)
  
}

result_30s_prep <- function(step_3, i, j, full = T){
  #   - g_ij_ground_truth                (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_ground_truth         (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_truth                       (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_truth                (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_X_truth                     (list of m x m matrices for i_j entries)
  #   - g_ij_X_coarse_truth              (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_est                         (list of m_est x m_est matrices for i_j entries)   
  
  key <- paste0(i, '_', j)
  
  if(! full){
    g_list <- list(visualize_matrix_heatmap(step_3$g_ij_truth[[key]], 'Truth_v1',        -1, 0, 1),
                   visualize_matrix_heatmap(step_3$g_ij_truth_v2[[key]], 'Truth_v2',     -1, 0, 1),
                   visualize_matrix_heatmap(step_3$g_ij_est[[key]], 'Estimate',          -1, 0, 1))
    return(g_list)
  }
  
  g_list <- list(visualize_matrix_heatmap(step_3[[1]][[key]], 'Ground Truth',        -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[2]][[key]], 'Coarse Ground Truth', -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[3]][[key]], 'Truth Theory',        -1, 0, 1),  
                 visualize_matrix_heatmap(step_3[[4]][[key]], 'Coarse Truth Theory', -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[5]][[key]], 'Truth X',             -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[6]][[key]], 'Coarse Truth X',      -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[7]][[key]], 'Estimate',            -1, 0, 1))
  
  return(g_list)
}

result_41_prep <- function(step_4, time_grid, time_grid_est){
  
  g_list <- list(visualize_log_intensity(t(step_4[[1]]$eigenfunctions[[1]]),    time_grid,      'Truth Theory'),
                 visualize_log_intensity(t(step_4[[2]]$eigenfunctions[[1]]),    time_grid_est,  'Coarse Truth Theory'),
                 visualize_log_intensity(t(step_4[[3]]$eigenfunctions[[1]]),    time_grid,      'Truth X'),
                 visualize_log_intensity(t(step_4[[4]]$eigenfunctions[[1]]),    time_grid_est,  'Coarse Truth X'),
                 visualize_log_intensity(t(step_4[[5]]$eigenfunctions[[1]]),    time_grid_est,  'Estimate'))  
  
  return(g_list)
  
}

result_42_prep <- function(step_3, step_4, p){
  
  source('functions/04_eigendecomposition.R')
  source('functions/13_estimation_validation.R')
  
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
  
  
  
  g_list <- list(visualize_matrix_heatmap(g_ii_truth[,,1],               'Truth Theory (TT)', -1, 0, 1),  # ground truth 50 x 50 covariance
                 visualize_matrix_heatmap(g_ii_truth_decomp[,,1],        'TT Reconstruct', -1,0, 1), # ground truth reconstructed
                 visualize_matrix_heatmap(g_ii_coarse_truth[,,1],        'Coarse Truth Theory (CTT)', -1, 0, 1), # ground truth 19 x 19 covariance
                 visualize_matrix_heatmap(g_ii_coarse_truth_decomp[,,1], 'CTT Reconstruct', -1, 0, 1), # ground truth 19 x 19 covariance
                 visualize_matrix_heatmap(g_ii_est[,,1],                 'Estimate (E)', -1, 0, 1),    # estimate 19 x 19 covariance
                 visualize_matrix_heatmap(g_ii_est_decomp[,,1],          'E Reconstruct', -1, 0, 1)) # reconstructed estimate
  
  return(g_list)

}

result_43_prep <- function(step_4){
  
  mat_truth          <- t(step_4[[1]]$eigenfunctions[[1]]) %*% step_4[[1]]$eigenfunctions[[1]]
  mat_coarse_truth   <- t(step_4[[2]]$eigenfunctions[[1]]) %*% step_4[[2]]$eigenfunctions[[1]]   
  mat_X_truth        <- t(step_4[[3]]$eigenfunctions[[1]]) %*% step_4[[3]]$eigenfunctions[[1]]
  mat_X_coarse_truth <- t(step_4[[4]]$eigenfunctions[[1]]) %*% step_4[[4]]$eigenfunctions[[1]]
  mat_est            <- t(step_4[[5]]$eigenfunctions[[1]]) %*% step_4[[5]]$eigenfunctions[[1]]
  
  g_list <- list(visualize_matrix_heatmap(mat_truth,   'Truth Theory',        zmid = 0),
                 visualize_matrix_heatmap(mat_truth,   'Coarse Truth Theory', zmid = 0),
                 visualize_matrix_heatmap(mat_truth,   'Truth X',             zmid = 0),
                 visualize_matrix_heatmap(mat_truth,   'Coarse Truth X',      zmid = 0),
                 visualize_matrix_heatmap(mat_truth,   'Estimate',            zmid = 0))
  
  return(g_list)
}

result_44_prep <- function(step_3, step_4, p){
  
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
  
  
  g_list <- list(visualize_error_histogram(g_ii_truth,          g_ii_truth_decomp,          'Truth Theory',        20),
                 visualize_error_histogram(g_ii_coarse_truth,   g_ii_coarse_truth_decomp,   'Coarse Truth Theory', 20),
                 visualize_error_histogram(g_ii_X_truth,        g_ii_X_truth_decomp,        'Truth X',             20),
                 visualize_error_histogram(g_ii_X_coarse_truth, g_ii_X_coarse_truth_decomp, 'Coarse Truth X',      20),
                 visualize_error_histogram(g_ii_est,            g_ii_est_decomp,            'Estimate',            20))
                         
  
  return(g_list)
}

result_80s_prep <- function(step_8, m, m_est, i, j){
  
  #   - V_cond_ground_truth_full          (pm x pm matrix)
  #   - V_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - V_cond_truth_full                 (pm x pm matrix)
  #   - V_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - V_cond_X_truth_full               (pm x pm matrix)
  #   - V_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - V_cond_est_full                   (pm_est x pm_est matrix)  
  
  g_list <- list(visualize_matrix_heatmap(extract_block_structure_ij(step_8[[1]], m,     1, 2), 'Ground Truth'), 
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[2]], m_est, 1, 2), 'Coarse Ground Truth'), 
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[3]], m,     1, 2), 'Truth Theory'),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[4]], m_est, 1, 2), 'Coarse Truth Theory'),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[5]], m,     1, 2), 'Truth X'), 
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[6]], m_est, 1, 2), 'Coarse Truth X'), 
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[7]], m_est, 1, 2), 'Estimate')) 
  
  return(g_list)
  
}

result_90s_prep_ij <- function(step_9, m, m_est, i, j, full = T){
  
  # list of heatmaps of C_Xi_Xj 
  # i = 1
  # j = 2
  
  if(! full){
    g_list <- list(visualize_matrix_heatmap(extract_block_structure_ij(step_9$C_cond_truth_full, m,     i, j), 'Truth',               zmid = 0),
                   visualize_matrix_heatmap(extract_block_structure_ij(step_9$C_cond_est_full,   m_est, i, j), 'Estimate',            zmid = 0))
    
    return(g_list)
  }
  
  g_list <- list(visualize_matrix_heatmap(extract_block_structure_ij(step_9[[1]], m,     i, j), 'Ground Truth',           zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[2]], m_est, i, j), 'Coarse Ground Truth',    zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[8]], m,     i, j), 'Ground Truth v2',        zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[9]], m_est, i, j), 'Coarse Ground Truth v2', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[3]], m,     i, j), 'Truth Theory',           zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[4]], m_est, i, j), 'Coarse Truth Theory',    zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[5]], m,     i, j), 'Truth X',                zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[6]], m_est, i, j), 'Coarse Truth X',         zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[7]], m_est, i, j), 'Estimate',               zmid = 0))
  
  return(g_list)
}

result_90s_prep_pm <- function(step_9, m, m_est, full = T){
  
  if(! full){
    g_list <- list(visualize_matrix_heatmap(step_9$C_cond_truth_full, 'Ground Truth',           zmid = 0),
                   visualize_matrix_heatmap(step_9$C_cond_est_full,   'Estimate',               zmid = 0))
    return(g_list)
  }
  
  g_list <- list(visualize_matrix_heatmap(step_9[[1]], 'Ground Truth',           zmid = 0),
                 visualize_matrix_heatmap(step_9[[2]], 'Coarse Ground Truth',    zmid = 0),
                 visualize_matrix_heatmap(step_9[[8]], 'Ground Truth v2',        zmid = 0),
                 visualize_matrix_heatmap(step_9[[9]], 'Coarse Ground Truth v2', zmid = 0),
                 visualize_matrix_heatmap(step_9[[3]], 'Truth Theory',           zmid = 0),
                 visualize_matrix_heatmap(step_9[[4]], 'Coarse Truth Theory',    zmid = 0),
                 visualize_matrix_heatmap(step_9[[5]], 'Truth X',                zmid = 0),
                 visualize_matrix_heatmap(step_9[[6]], 'Coarse Truth X',         zmid = 0),
                 visualize_matrix_heatmap(step_9[[7]], 'Estimate',               zmid = 0))
  
  
  
  return(g_list)
}

result_95_prep <- function(step_9, m, m_est, p, full = T){
  
  
  if(! full){
    C_HS_truth_full                 <- hilbert_schmidt_norm_pm(step_9$C_cond_truth_full, p, m) 
    C_HS_est                        <- hilbert_schmidt_norm_pm(step_9$C_cond_est_full,   p, m_est)
    
    diag(C_HS_truth_full) <- 0
    diag(C_HS_est) <- 0
    
    g_list <- list(visualize_matrix_heatmap(C_HS_truth_full,        'Truth',        zmid = 0),
                   visualize_matrix_heatmap(C_HS_est,               'Estimate',     zmid = 0)) 
    
    return(g_list)
  }
  
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
  
  g_list <- list(visualize_matrix_heatmap(C_HS_ground_truth_full,        'Ground Truth',        zmid = 0),
                 visualize_matrix_heatmap(C_HS_coarse_ground_truth_full, 'Coarse Ground Truth', zmid = 0),
                 visualize_matrix_heatmap(C_HS_truth_full,               'Truth Theory',        zmid = 0),
                 visualize_matrix_heatmap(C_HS_coarse_truth_full,        'Coarse Truth Theory', zmid = 0),
                 visualize_matrix_heatmap(C_HS_X_truth_full,             'Truth X',             zmid = 0),
                 visualize_matrix_heatmap(C_HS_X_coarse_truth_full,      'Coarse Truth X',      zmid = 0),
                 visualize_matrix_heatmap(C_HS_est,                      'Estimate',            zmid = 0))        
  
  return(g_list)
  
}

result_100s_prep_ij <- function(step_10, m, m_est, i, j){
  #   - P_cond_ground_truth_full          (pm x pm matrix)
  #   - P_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - P_cond_truth_full                 (pm x pm matrix)
  #   - P_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - P_cond_X_truth_full               (pm x pm matrix)
  #   - P_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - P_cond_est_full                   (pm_est x pm_est matrix)    
  g_list <- list(visualize_matrix_heatmap(extract_block_structure_ij(step_10[[1]], m,     i, j), 'Truth Ground',         zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[2]], m_est, i, j), 'Coarse Truth Ground',  zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[3]], m,     i, j), 'Truth Theory',         zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[4]], m_est, i, j), 'Coarse Truth Theory',  zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[5]], m,     i, j), 'Truth X',              zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[6]], m_est, i, j), 'Coarse Truth X',       zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[7]], m_est, i, j), 'Estimate',             zmid = 0))
  
  return(g_list)
}

result_100s_prep_pm <- function(step_10, m, m_est){
  
  g_list <- list(visualize_matrix_heatmap(step_10[[1]], 'Ground Truth',           zmid = 0),
                 visualize_matrix_heatmap(step_10[[2]], 'Coarse Ground Truth',    zmid = 0),
                 # visualize_matrix_heatmap(step_10[[8]], 'Ground Truth v2',        zmid = 0),
                 # visualize_matrix_heatmap(step_10[[9]], 'Coarse Ground Truth v2', zmid = 0),
                 visualize_matrix_heatmap(step_10[[3]], 'Truth Theory',           zmid = 0),
                 visualize_matrix_heatmap(step_10[[4]], 'Coarse Truth Theory',    zmid = 0),
                 visualize_matrix_heatmap(step_10[[5]], 'Truth X',                zmid = 0),
                 visualize_matrix_heatmap(step_10[[6]], 'Coarse Truth X',         zmid = 0),
                 visualize_matrix_heatmap(step_10[[7]], 'Estimate',               zmid = 0))
  
  
  
  return(g_list)
}

result_112_prep <- function(step_11, remove_diag){
  #   - w_mat_ground_truth        (p x p)   matrix of HS norms of the pm x pm ground truth
  #   - w_mat_coarse_ground_truth (p x p)   matrix of HS norms of the pm_est x pm_est ground truth 
  #   - w_mat_X_coarse_truth      (p x p)   matrix of HS norms of ...
  #   - w_mat_X_truth             (p x p)   matrix of HS norms of ...
  #   - w_mat_coarse_truth        (p x p)   matrix of HS norms of ...
  #   - w_mat_truth               (p x p)   matrix of HS norms of ...
  #   - w_mat_est                 (p x p)   matrix of HS norms of ...  
  
  if(remove_diag){
    step_11 <- lapply(step_11, function(x) {
      diag(x) <- 0
      x
    })
  }
  
  g_list <- list(visualize_matrix_heatmap(step_11[[1]], 'Ground Truth',        zmid = 0), 
                 visualize_matrix_heatmap(step_11[[2]], 'Coarse Ground Truth', zmid = 0), 
                 visualize_matrix_heatmap(step_11[[3]], 'Truth Theory',        zmid = 0),   
                 visualize_matrix_heatmap(step_11[[4]], 'Coarse Truth Theory', zmid = 0), 
                 visualize_matrix_heatmap(step_11[[5]], 'Truth X',             zmid = 0), 
                 visualize_matrix_heatmap(step_11[[6]], 'Coarse Truth X',      zmid = 0),
                 visualize_matrix_heatmap(step_11[[7]], 'Estimate',            zmid = 0))
  
  return(g_list)
  
}

result_113_prep <- function(step_12){
  
  #   - roc_ground_truth           (list of roc outputs)
  #   - roc_coarse_ground_truth
  #   - roc_truth
  #   - roc_coarse_truth
  #   - roc_X_truth
  #   - roc_X_coarse_truth
  #   - roc_est
  
  # ROC plot - done on qrsh 
  roc_graphs <- list()
  
  title_names <- c('Ground Truth', 'Coarse Ground Truth', 'Truth Theory', 'Coarse Truth Theory', 'X Truth', 'X Coarse Truth',
                   'Estimate')
  
  for(i in 1:7){
    
    roc_df          <- step_12[[i]]$roc_df
    ideal_spec      <- step_12[[i]]$specificity
    ideal_sens      <- step_12[[i]]$sensitivity
    ideal_threshold <- step_12[[i]]$threshold
    auc_value       <- step_12[[i]]$auc
    
    roc_graphs[[i]] <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
      geom_step(direction = "vh", color = "blue", size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
      labs(title = title_names[i], x = "False Positive Rate", y = "True Positive Rate") +
      annotate("point", x = 1 - ideal_spec, y = ideal_sens, color = "red", size = 3) +
      annotate("text", x = 1 - ideal_spec, y = ideal_sens, 
               label = paste0("Threshold=", round(ideal_threshold, 3)),
               hjust = -0.1, vjust = -0.5, color = "red") +
      
      annotate("text", x = 0.6, y = 0.2,            # position for AUC label
               label = paste0("AUC = ", round(auc_value, 3)),
               color = "darkgreen", size = 5) +    
      theme_minimal()         
  }
  
  return(roc_graphs)
  
}

result_arr_mat <- function(g_list, grob_caption, arr_mat_i){
  
  # reused function to arrange all graphs in an arr_mat_8 fashion
  
  nrows <- dim(arr_mat_i)[1]
  ncols <- dim(arr_mat_i)[2]
  
  if(nrows == 2 & ncols == 2){
    g <- grid.arrange(g_list[[1]], 
                      g_list[[2]], 
                      textGrob(grob_caption, gp = gpar(fontsize = 14)),
                      g_list[[3]], layout_matrix = arr_mat_i)
  }
  
  if(nrows == 2 & ncols == 3){
    
    if(length(g_list) == 4){
      g <- grid.arrange(g_list[[1]], 
                        g_list[[2]], 
                        g_list[[3]], 
                        g_list[[4]], 
                        textGrob(grob_caption, gp = gpar(fontsize = 14)),
                        layout_matrix = arr_mat_i)      
    } else{
      g <- grid.arrange(g_list[[1]], 
                        g_list[[2]], 
                        g_list[[3]], 
                        g_list[[4]], 
                        textGrob(grob_caption, gp = gpar(fontsize = 14)),
                        g_list[[5]], 
                        layout_matrix = arr_mat_i)      
    }
    

  }  
  
  if(nrows == 2 & ncols == 4){
    
    if(length(g_list) == 6){
      g <- grid.arrange(g_list[[1]], g_list[[2]], g_list[[3]],
                        g_list[[4]], g_list[[5]], g_list[[6]],
                        textGrob(grob_caption, gp = gpar(fontsize = 14)),
                        layout_matrix = arr_mat_i)       
    } else{

      # often times, it's ground truth (1, 2) truth theory (3, 4) X truth (5, 6) and estimate (7)
      
      g <- grid.arrange(g_list[[1]], g_list[[2]], g_list[[3]],
                        g_list[[4]], g_list[[5]], g_list[[6]],
                        textGrob(grob_caption, gp = gpar(fontsize = 14)),
                        g_list[[7]],
                        layout_matrix = arr_mat_i)   
    }
  }
  
  if(nrows == 2 & ncols == 5){
    # often times, it's ground truth (1, 2) ground truth v2 (3, 4) truth theory (5, 6) X truth (7, 8) and estimate (9)
    g <- grid.arrange(g_list[[1]], g_list[[2]], g_list[[3]], g_list[[4]],
                      g_list[[5]], g_list[[6]], g_list[[7]], g_list[[8]],
                      textGrob(grob_caption, gp = gpar(fontsize = 14)),
                      g_list[[9]], 
                      layout_matrix = arr_mat_i)    
  }
  

  
  return(g)
}



