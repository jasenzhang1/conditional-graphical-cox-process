source('functions/00c_block_matrix_arrange.R')
source('functions/28z_Visualization_helpers.R')
source('functions/28_Simulation_Visualization.R')

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

result_11 <- function(step_1, time_grid, time_grid_est, arr_mat){
  
  g_list <- result_11_prep(step_1, time_grid, time_grid_est)
  
  
  g <- grid.arrange(g_list[[1]],
                    g_list[[2]],
                    g_list[[3]],
                    g_list[[4]],
                    textGrob("0. Log Intensity\n of first 5 processes\nof subject 1", gp = gpar(fontsize = 14)),
                    layout_matrix = arr_mat)   
  
  return(g)
}

result_22_prep <- function(step_2, time_grid, time_grid_est){
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

result_22 <- function(step_2, time_grid, time_grid_est, arr_mat){
  
  # arrange them 
  
  g_list <- result_22_prep(step_2, time_grid, time_grid_est)
  
  g <- grid.arrange(g_list[[1]], g_list[[2]], g_list[[3]], g_list[[4]],
                    textGrob("2. Rho_i\nEstimation for \n first 5 processes", gp = gpar(fontsize = 14)),
                    g_list[[5]], 
                    layout_matrix = arr_mat)
  
  return(g)
  
}

result_20s_prep <- function(step_2b, i, j){
  #   - rho_ii_truth               (list of m x m matrices)
  #   - rho_ii_coarse_truth        (list of m_est x m_est matrices)
  #   - rho_ii_X_truth             (list of m x m matrices)
  #   - rho_ii_X_coarse_truth      (list of m_est x m_est matrices)
  #   - rho_ii_est                 (list of m_est x m_est matrices)
  
  key <- paste0(i, '_', j)
  
  g_list <- list(visualize_matrix_heatmap(step_2b[[1]][[key]],   'Truth Theory',        40000, NULL, 200000),
                 visualize_matrix_heatmap(step_2b[[2]][[key]],   'Coarse Truth Theory', 40000, NULL, 200000),
                 visualize_matrix_heatmap(step_2b[[3]][[key]],   'Truth X',             40000, NULL, 200000),
                 visualize_matrix_heatmap(step_2b[[4]][[key]],   'Coarse Truth X',      40000, NULL, 200000),
                 visualize_matrix_heatmap(step_2b[[5]][[key]],   'Estimate',            40000, NULL, 200000))
  
  return(g_list)
}

result_20s <- function(step_2b, i, j, arr_mat){
  
  g_list <- result_20s_prep(step_2b, i, j)
  
  g <- grid.arrange(g_list[[1]],
                    g_list[[2]],
                    g_list[[3]],
                    g_list[[4]],
                    textGrob("2. Rho ii\nEstimation \n for process 1", gp = gpar(fontsize = 14)),
                    g_list[[5]],
                    layout_matrix = arr_mat)
  
  return(g)
  
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

result_30s_prep <- function(step_3, i, j){
  #   - g_ij_ground_truth                (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_ground_truth         (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_truth                       (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_truth                (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_X_truth                     (list of m x m matrices for i_j entries)
  #   - g_ij_X_coarse_truth              (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_est                         (list of m_est x m_est matrices for i_j entries)   
  
  key <- paste0(i, '_', j)
  
  g_list <- list(visualize_matrix_heatmap(step_3[[1]][[key]], 'Ground Truth',        -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[2]][[key]], 'Coarse Ground Truth', -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[3]][[key]], 'Truth Theory',        -1, 0, 1),  
                 visualize_matrix_heatmap(step_3[[4]][[key]], 'Coarse Truth Theory', -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[5]][[key]], 'Truth X',             -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[6]][[key]], 'Coarse Truth X',      -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[7]][[key]], 'Estimate',            -1, 0, 1))
  
  return(g_list)
}

result_30s <- function(step_3, i, j, arr_mat_8){
  
  g_list <- result_30s_prep(step_3, i, j)
  
  g <- grid.arrange(g_list[[1]],
                    g_list[[2]],
                    g_list[[3]],
                    g_list[[4]],
                    g_list[[5]],
                    g_list[[6]],
                    textGrob("3. Covariance Function\nEstimation (G_ii)", gp = gpar(fontsize = 14)),
                    g_list[[7]],
                    layout_matrix = arr_mat_8)  
  
  return(g)
  
}

result_90s_prep <- function(step_9, m, m_est, i, j){
  
  # list of heatmaps of C_Xi_Xj 
  # i = 1
  # j = 2
  
  g_list <- list(visualize_matrix_heatmap(extract_block_structure_ij(step_9[[1]], m,     i, j), 'Ground Truth', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[2]], m_est, i, j), 'Coarse Ground Truth', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[8]], m,     i, j), 'Ground Truth v2', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[9]], m_est, i, j), 'Coarse Ground Truth v2', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[3]], m,     i, j), 'Truth Theory', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[4]], m_est, i, j), 'Coarse Truth Theory', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[5]], m,     i, j), 'Truth X', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[6]], m_est, i, j), 'Coarse Truth X', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[7]], m_est, i, j), 'Estimate', zmid = 0))
  
  return(g_list)
}


result_90s <- function(step_9, m, m_est, arr_mat_10){
  
  # arrange them 
  
  g_list <- result_90s_prep(step_9, m, m_est)
  
  g <- grid.arrange(g_list[[1]], g_list[[2]], g_list[[3]], g_list[[4]],
                    g_list[[5]], g_list[[6]], g_list[[7]], g_list[[8]],
                    textGrob("9. Conditional\nCorrelation Operator\nProcess 1 with 2", gp = gpar(fontsize = 14)),  
                    g_list[[9]], 
                    layout_matrix = arr_mat_10)
  
  return(g)
  
}
