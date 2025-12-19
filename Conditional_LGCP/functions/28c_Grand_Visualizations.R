source('functions/28_Simulation_Visualization.R')
source('functions/28b_Simulation_Visualization_2.R')
source('functions/28y_Visualization_Blocks.R')
source('functions/28z_Visualization_helpers.R')
source('functions/00c_block_matrix_arrange.R')


# these functions open an entire dataset and search through all query_id 


visualize_prec_mat_over_time <- function(folder_name, n){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: plot true pxp matrix versus estimated pxp matrix (w_mat)
  #
  # - this is to be done for a single simulation (adj_method, est_method, n)
  # 
  # input: 
  #
  # - folder_name   (string)
  # - n             (integer)
  #
  # output:
  #
  # - 2 x n graph of prec_mat (row 1) and w_mat (row 2)
  #
  # 
  # ----------------------------------------------------------------------------
  
  
  # load 
  graph_results_i <- get_file_name_and_load(folder_name, n) 
  
  
  # 2) get all the step_0 adjacencies
  
  step_0 <- graph_results_i$estimated_graphs_part_1$step_0
  step_11 <- lapply(graph_results_i$estimated_graphs_part_2, function(x) x$step_11$w_mat_est)
  n_graphs <- length(step_0)
  
  
  
  # 3) Build all n plots of ground truth
  plots_truth <- lapply(seq_len(n_graphs), function(i) {
    visualize_matrix_heatmap(step_0[[i]],
                             paste0('y_c = ', as.character(i)), -2, NULL, 2)
  })  
  
  plots_w_est <- lapply(seq_len(n_graphs), function(i) {
    visualize_matrix_heatmap(step_11[[i]],
                             paste0('y_c = ', as.character(i)), zmid = 0)
  }) 
  
  return(list(graph = grid.arrange(grobs = c(plots_truth, plots_w_est), nrow = 2, ncol = n_graphs,
                                   left = textGrob("Est \t Ground Truth", rot = 90, gp = gpar(fontsize = 16))),
              values = list(truth = step_0,
                            est = step_11))
  )
  
  
  
}


# helper 
rearrange_plots <- function(g_list){
  
  # g_list is a list of lists
  # - first layer is for each y_c
  # - second layer is for (truth theory, est etc)
  
  
  n_time <- length(g_list)
  n_graphs <- length(g_list[[1]])
  N <- n_time * n_graphs
  
  # unlist the plots, reorganize them
  plots <- unlist(g_list, recursive = FALSE)
  indices <- as.vector(sapply(1:n_graphs, function(r) r + (0:(n_time-1))*n_graphs))
  plots_colwise <- plots[indices]
  
  
  # plot them
  arranged_plot <- do.call(grid.arrange, c(plots_colwise, nrow = n_graphs, ncol = n_time))  
  
  return(arranged_plot)
}

# everything over all query_id

visualize_over_time <- function(graph_results_i, graph_ids, beta_truth, X_truth){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: for a single dataset, plot a specific intermediate value over time
  #
  # - this is to be done for a single simulation (adj_method, est_method, n)
  # 
  # input: 
  #
  # - graph_results_i   (list of all step_X's)  --> (list of all y_c_queries) --> (list of 'est' or 'truth' or 'coarse_truth' etc)
  # - graph_ids         (vector of strings)    which graphs do we want?
  # - i
  # - j
  # - beta_truth            (boolean)             do our results have beta_truth values? 
  # - X_truth               (boolean)             do our results have X_truth values? 
  #
  # 
  # output:
  #
  # - graphs   (list)  list of graphs
  # 
  # ---------------------------------------------------------------------------- 
  
  
  # 1) load
  

  # step_0 <- step_list$step_0 
  # step_1 <- step_list$step_1 
  # step_2 <- graph_results_i$step_2
  # step_3 <- graph_results_i$step_3
  # step_4 <- graph_results_i$step_4
  # step_5 <- graph_results_i$step_5
  # step_6 <- graph_results_i$step_6
  # step_7 <- graph_results_i$step_7
  # step_8 <- graph_results_i$step_8 
  # step_9 <- graph_results_i$step_9
  # step_10 <- graph_results_i$step_10
  # step_11 <- graph_results_i$step_11
  # step_12 <- graph_results_i$step_12 
  # step_2b <- graph_results_i$step_2b
  
  # graph_results_i also contains time_grid, time_grid_est, y_c_query
  list2env(graph_results_i, envir = environment())  # does all of the above in one go
  
  # prep
  arr_mat_6 <- matrix(1:6, nrow = 2, byrow = F)
  arr_mat_8 <- matrix(1:8, nrow = 2, byrow = F)
  arr_mat_10 <- matrix(1:10, nrow = 2, byrow = F)  
  
  m <- length(time_grid)
  m_est <- length(time_grid_est)
  
  graphs <- list()
  
  y_c_names <- graph_results_i$y_c_query %>%
    as.numeric() %>%
    round(3) %>%
    formatC(format = "f", digits = 3)
  
  # 2) figure out which steps have unique results for y_c_query and which steps are constant
  
  # log-intensity of first 5 processes of subject i
  if('11' %in% graph_ids){
    
    i <- 1
    if(X_truth){
      graphs[['g_11']] <- grid.arrange(visualize_log_intensity(step_1$X_k_est[1:5,,i],            time_grid_est,  'Estimate',       step_1c$mu_t_coarse_truth ),
                                       visualize_log_intensity(step_1$X_k_coarse_truth[1:5,,i],   time_grid_est,  'Coarser Truth',  step_1c$mu_t_coarse_truth ),
                                       visualize_log_intensity(step_1$X_k_truth[1:5,,i],          time_grid,      'Finer Truth',    step_1c$mu_t_truth),
                                       visualize_log_intensity(step_1$X_k_both_truth[1:5,,i],     time_grid_both, 'Combined Truth', step_1c$mu_t_both_truth),
                                       textGrob("0. Log Intensity\n of first 5 processes\nof subject 1", gp = gpar(fontsize = 14)),
                                       layout_matrix = arr_mat_6) 
    } else{
      graphs[['g_11']] <- grid.arrange(visualize_log_intensity(step_1$X_k_est[1:5,,i],            time_grid_est,  'Estimate',       step_1c$mu_t_coarse_truth ),
                                       textGrob("0. Log Intensity\n of first 5 processes\nof subject 1", gp = gpar(fontsize = 14)),
                                       layout_matrix = arr_mat_6) 
    }
  }
  
  # average log-intensity of first 5 processes
  if('12' %in% graph_ids){
    
    if(X_truth){
      # p x m x n --> mean --> p x m --> choose first 5 processes --> 5 x m
      graphs[['g_12']] <- grid.arrange(visualize_log_intensity(apply(step_1$X_k_est,            c(1, 2), mean)[1:5, ],   time_grid_est,  'Estimate',        step_1c$mu_t_coarse_truth),
                                       visualize_log_intensity(apply(step_1$X_k_coarse_truth,   c(1, 2), mean)[1:5, ],   time_grid_est,  'Coarser Truth',   step_1c$mu_t_coarse_truth),
                                       visualize_log_intensity(apply(step_1$X_k_truth,          c(1, 2), mean)[1:5, ],   time_grid,      'Finer Truth',     step_1c$mu_t_coarse_truth),
                                       visualize_log_intensity(apply(step_1$X_k_both_truth,     c(1, 2), mean)[1:5, ],   time_grid_both, 'Combined Truth',  step_1c$mu_t_truth),
                                       textGrob("0. Average Log Intensity\n of first 5 processes", gp = gpar(fontsize = 14)),
                                       layout_matrix = arr_mat_6) 
    } else{
      graphs[['g_12']] <- grid.arrange(visualize_log_intensity(apply(step_1$X_k_est,          c(1, 2), mean)[1:5, ],   time_grid_est,  'Estimate',     step_1c$mu_t_coarse_truth),
                                       textGrob("0. Average Log Intensity\n of first 5 processes", gp = gpar(fontsize = 14)),
                                       layout_matrix = arr_mat_6) 
    }
  }
  
  # rho_i(t) for processes 1 through 5
  if('22' %in% graph_ids){
    graphs[['g_22']] <- result_line_graph_prep(step_2, 'rho_i', time_grid_est, num_processes = 5)
  }
  
  # rho_ij(s,t) for process pair 1_1
  if('24' %in% graph_ids){ 
    graphs[['g_24']] <- result_heatmap_ij_prep(step_2b, 'rho_ii', time_grid_est, i = 1, j = 1)
  }  
  
  # rho_ij(s,t) for process pair 1_2
  if('25' %in% graph_ids){ 
    graphs[['g_25']] <- result_heatmap_ij_prep(step_2b, 'rho_ii', time_grid_est, i = 1, j = 2)
  }  
  
  # rho_ij(s,t) for process pair 3_5
  if('26' %in% graph_ids){ 
    graphs[['g_26']] <- result_heatmap_ij_prep(step_2b, 'rho_ii', time_grid_est, i = 3, j = 5)
  }   
  
  # weights, all y_c_query settings in a row
  if('29' %in% graph_ids){
    graph_list <- lapply(seq_along(weights), function(i) {
      result_29(weights[[i]], Y_continuous, y_c_names[i])
    })
    
    graphs[['g_29']] <- do.call(grid.arrange, c(graph_list, nrow = 1))
  }
  
  # g_ij(s,t) at 1_1
  if('31' %in% graph_ids){ 
    graphs[['g_31']] <- result_heatmap_ij_prep(step_3, 'g_ij', time_grid_est, i = 1, j = 1, zmid = 0)
  }  
  
  # g_ij(s,t) at 1_2
  if('32' %in% graph_ids){ 
    graphs[['g_32']] <- result_heatmap_ij_prep(step_3, 'g_ij', time_grid_est, i = 1, j = 2, zmid = 0)
  }  
  
  # g_ij(s,t) at 3_5
  if('33' %in% graph_ids){ 
    graphs[['g_33']] <- result_heatmap_ij_prep(step_3, 'g_ij', time_grid_est, i = 3, j = 5, zmid = 0)
  }    
  
  # eigenfunctions
  if('41' %in% graph_ids){ 
    step_4_v2 <- result_41_prep(step_4)  
    
    graphs[['g_41']] <- result_line_graph_prep(step_4_v2, 'eigen_decomp', time_grid_est)
  }  
  
  # reconstructing g_11 from eigenfunctions
  if('42' %in% graph_ids){ 
    g_list <- lapply(1:length(step_3), function(i){result_42_prep(step_3[[i]], step_4[[i]], p)})
    graphs[['g_42']] <- rearrange_plots(g_list)
  }  
  
  # orthonormality of eigenfunctions
  if('43' %in% graph_ids){ 
    g_list <- lapply(step_4, function(x){result_43_prep(x)})
    graphs[['g_43']] <- rearrange_plots(g_list)
  }
  
  # reconstruction error histogram 
  if('44' %in% graph_ids){ 
    g_list <- lapply(1:length(step_3), function(i){result_44_prep(step_3[[i]], step_4[[i]], p, X_truth = F)})
    graphs[['g_44']] <- rearrange_plots(g_list)
  } 
  # eigenvalue decay - check if G_ii / m makes the eigenvalues similar between coarse/fine settings
  if('45' %in% graph_ids){
    g_list <- lapply(step_4, function(x) result_45_prep(x))
    graphs[['g_45']] <- rearrange_plots(g_list) 
  }
  
  # KL_cov of (1, 2) block
  if('55' %in% graph_ids){
    g_list <- lapply(step_5, function(x) result_55_prep(x, i = 1, j = 2))
    graphs[['g_55']] <- rearrange_plots(g_list) 
  }
  
  # KL_cor of (1, 2) block
  if('56' %in% graph_ids){
    g_list <- lapply(step_5b, function(x) result_56_prep(x, i = 1, j = 2))
    graphs[['g_56']] <- rearrange_plots(g_list) 
  }
  
  # KL_prec of (1, 2) block
  if('57' %in% graph_ids){
    g_list <- lapply(step_5c, function(x) result_57_prep(x, i = 1, j = 2))
    graphs[['g_57']] <- rearrange_plots(g_list) 
  }
  
  # KL_cor entire matrix
  if('58' %in% graph_ids){
    g_list <- lapply(step_5b, function(x) result_58_prep(x, p))
    graphs[['g_58']] <- rearrange_plots(g_list) 
  }
  
  # KL_prec entire matrix
  if('59' %in% graph_ids){
    g_list <- lapply(step_5c, function(x) result_59_prep(x, p))
    graphs[['g_59']] <- rearrange_plots(g_list) 
  }
  
  # V_Xi_Xj
  if('81' %in% graph_ids & 'step_8' %in% names(graph_results_i)){
    est_graphs <- lapply(graph_results_i$step_8, function(x) extract_block_structure_ij(x$step_8$V_cond_est_full, m_est, i = 1, j = 2))
  }
  
  # C_Xi_Xj for block (1, 1)
  
  if('90' %in% graph_ids){
    graphs[['g_90']] <- result_heatmap_ij_prep(step_9, 'C_cond', time_grid_est, i = 1, j = 1, zmid = 0)
  }  
  
  # C_Xi_Xj for block (1, 2)
  if('91' %in% graph_ids){
    graphs[['g_91']] <- result_heatmap_ij_prep(step_9, 'C_cond', time_grid_est, i = 1, j = 2, zmid = 0)
  }
  
  # C_Xi_Xj for block (3, 5)
  if('92' %in% graph_ids){
    graphs[['g_92']] <- result_heatmap_ij_prep(step_9, 'C_cond', time_grid_est, i = 3, j = 5, zmid = 0)
  }  
  
  # C for all blocks (pm x pm)
  if('93' %in% graph_ids){
    graphs[['g_93']] <- result_heatmap_nonblock_prep(step_9, 'C_cond', time_grid_est, data_format = 'list', zmid = 0)
  }  
  
  # C_HS for all pxp blocks
  if('95' %in% graph_ids){
    graphs[['g_95']] <- result_heatmap_nonblock_prep(step_11b, 'C_HS', time_grid_est, data_format = 'regular', rm_diag = T, zmid = 0)
  }
  
  # distribution of C_HS
  if('96' %in% graph_ids){
    graphs[['g_96']] <- result_histogram_prep(step_11b, 'C_HS', data_format = 'regular', nbins = 40)
  }
  
  # P_Xi_Xj for block (1, 1)
  
  if('100' %in% graph_ids){
    graphs[['g_100']] <- result_heatmap_ij_prep(step_10, 'P_cond', time_grid_est, i = 1, j = 1, zmid = 0)
  }  
  
  # P_Xi_Xj for block (1, 2)
  if('101' %in% graph_ids){
    graphs[['g_101']] <- result_heatmap_ij_prep(step_10, 'P_cond', time_grid_est, i = 1, j = 2, zmid = 0)
  }
  
  # P_Xi_Xj for block (3, 5)
  if('102' %in% graph_ids){
    graphs[['g_102']] <- result_heatmap_ij_prep(step_10, 'P_cond', time_grid_est, i = 3, j = 5, zmid = 0)
  }  
  
  # P for all blocks (pm x pm)
  if('103' %in% graph_ids){
    graphs[['g_103']] <- result_heatmap_nonblock_prep(step_10, 'P_cond', time_grid_est, data_format = 'list', rm_diag = T, zmid = 0)
  } 
  
  # P_HS (w_mat) for all pxp blocks
  if('112' %in% graph_ids){
    graphs[['g_112']] <- result_heatmap_nonblock_prep(step_11, 'w_mat', time_grid_est, data_format = 'regular', rm_diag = T, zmid = 0)
  }  
  
  # P_HS (w_mat) for all pxp blocks - delete est and X_truth
  if('112' %in% graph_ids){
    step_11_v2 <- step_11
    step_11_v2 <- lapply(step_11_v2, function(x) {
      x[!names(x) %in% c("w_mat_est", "w_mat_X_truth")]
    })
    graphs[['g_112b']] <- result_heatmap_nonblock_prep(step_11_v2, 'w_mat', time_grid_est, data_format = 'regular', rm_diag = T, zmid = 0)
  } 
  
  # ROC curves
  if('113' %in% graph_ids){
    
    y_names <- names(graph_results_i$step_12)
    x_names <- names(graph_results_i$step_12[[1]])
    
    all_plots <- list()
    
    for (i in seq_along(graph_results_i$step_12)) {
      for (j in seq_along(graph_results_i$step_12[[i]])) {
        
        roc_entry       <- graph_results_i$step_12[[i]][[j]]
        roc_df          <- roc_entry$roc_df
        ideal_spec      <- roc_entry$specificity
        ideal_sens      <- roc_entry$sensitivity
        ideal_threshold <- roc_entry$threshold
        auc_value       <- roc_entry$auc
        
        p <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
          geom_step(direction = "vh", color = "blue", size = 1) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
          labs(
            title = paste('y=', y_names[i], ', ', x_names[j]),
            x = "False Positive Rate", y = "True Positive Rate"
          ) +
          annotate("point", x = 1 - ideal_spec, y = ideal_sens, color = "red", size = 3) +
          annotate("text", x = 1 - ideal_spec, y = ideal_sens,
                   label = paste0("thr=", round(ideal_threshold, 3)),
                   hjust = -0.1, vjust = -0.5, color = "red") +
          annotate("text", x = 0.6, y = 0.2,
                   label = paste0("AUC=", round(auc_value, 3)),
                   color = "darkgreen", size = 5) +
          theme_minimal(base_size = 10)
        
        all_plots[[length(all_plots) + 1]] <- p
      }
    }
    
    # arrange all 70 plots in a 10x7 grid
    graphs[['113']] <- wrap_plots(all_plots, ncol = length(y_names), nrow = length(x_names), byrow = FALSE) 
    
  }
  
  return(graphs)
  
  # Final adj_mat for all pxp blocks
  if('114' %in% graph_ids){
    graphs[['g_114']] <- result_heatmap_nonblock_prep(step_12b, 'adj_mat', time_grid_est, data_format = 'regular', rm_diag = T, zmid = 0)
  }   
  
}









