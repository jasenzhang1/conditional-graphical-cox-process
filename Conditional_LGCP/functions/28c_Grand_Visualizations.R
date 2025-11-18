source('functions/28_Simulation_Visualization.R')
source('functions/28b_Simulation_Visualization_2.R')
source('functions/28y_Visualization_Blocks.R')
source('functions/28z_Visualization_helpers.R')
source('functions/00c_block_matrix_arrange.R')

library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)

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

visualize_over_time <- function(graph_results_i, graph_ids, full = T){
  
  
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
  # - full
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
    round(2) %>%
    formatC(format = "f", digits = 2)
  
  # 2) figure out which steps have unique results for y_c_query and which steps are constant
  
  if('11' %in% graph_ids){
    
    if(! full){
      graphs[['g_11']] <- grid.arrange(visualize_log_intensity(step_1$X_k_est[1:5,,1],            time_grid_est,  'Estimate',       step_1$mu_t_coarse_truth ),
                                       visualize_log_intensity(step_1$X_k_coarse_truth[1:5,,1],   time_grid_est,  'Coarser Truth',  step_1$mu_t_coarse_truth ),
                                       visualize_log_intensity(step_1$X_k_truth[1:5,,1],          time_grid,      'Finer Truth',    step_1$mu_t_truth),
                                       visualize_log_intensity(step_1$X_k_both_truth[1:5,,1],     time_grid_both, 'Combined Truth', step_1$mu_t_both_truth),
                                       textGrob("0. Log Intensity\n of first 5 processes\nof subject 1", gp = gpar(fontsize = 14)),
                                       layout_matrix = arr_mat_6) 
    }
    
    step_1 <- graph_results_i$step_1
    
    mu_t <- graph_results_i$kernel_params_i[[1]]$base_mean
    mu_t_est <- graph_results_i$kernel_params_i[[1]]$base_mean_est
    mu_t_both <- graph_results_i$kernel_params_i[[1]]$base_mean_both
    
    graphs[['g_11']] <- grid.arrange(visualize_log_intensity(step_1[[1]][1:5,,1],   time_grid_est,  'Estimate', mu_t_est ),
                                     visualize_log_intensity(step_1[[3]][1:5,,1],   time_grid_est,  'Coarser Truth', mu_t_est ),
                                     visualize_log_intensity(step_1[[2]][1:5,,1],   time_grid,      'Finer Truth', mu_t     ),
                                     visualize_log_intensity(step_1[[4]][1:5,,1],   time_grid_both, 'Combined Truth', mu_t_both),
                                     textGrob("0. Log Intensity\n of first 5 processes\nof subject 1", gp = gpar(fontsize = 14)),
                                     layout_matrix = arr_mat_6) 
  }
  
  if('12' %in% graph_ids){
    
    if(! full){
      
      # p x m x n --> mean --> p x m --> choose first 5 processes --> 5 x m
      
      graphs[['g_12']] <- grid.arrange(visualize_log_intensity(apply(step_1$X_k_est,          c(1, 2), mean)[1:5, ],   time_grid_est,  'Estimate',     step_1$mu_t_est),
                                       visualize_log_intensity(apply(step_1$X_k_truth,        c(1, 2), mean)[1:5, ],   time_grid,      'Truth',        step_1$mu_t_truth),
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
    
    # key <- '1_1'
    # 
    # zmin <- min(sapply(step_2b, function(x) c(as.numeric(x$rho_ii_est[[key]]), as.numeric(x$rho_ii_truth[[key]]))), na.rm = TRUE)   # global min
    # zmax <- max(sapply(step_2b, function(x) c(as.numeric(x$rho_ii_est[[key]]), as.numeric(x$rho_ii_truth[[key]]))), na.rm = TRUE)   # global min
    # 
    # zmin = zmin - 0.01 * (zmax - zmin)  # buffer area
    # zmax = zmax + 0.01 * (zmax - zmin)    
    # 
    # g_list <- lapply(step_2b, function(x){result_20s_prep(x, i = 1, j = 1, zmin, zmax, full)})
    # graphs[['g_24']] <- rearrange_plots(g_list)
    
    # new way
    
    graphs[['g_24']] <- result_heatmap_ij_prep(step_2b, 'rho_ii', i = 1, j = 1)
  }  
  
  # rho_ij(s,t) for process pair 1_2
  if('25' %in% graph_ids){ 
    
    graphs[['g_25']] <- result_heatmap_ij_prep(step_2b, 'rho_ii', i = 1, j = 2)
  }  
  
  # rho_ij(s,t) for process pair 3_5
  if('26' %in% graph_ids){ 
    
    graphs[['g_26']] <- result_heatmap_ij_prep(step_2b, 'rho_ii', i = 3, j = 5)
  }   
  
  # weights, all y_c_query settings in a row
  if('29' %in% graph_ids){
    
    graph_list <- lapply(seq_along(step_2), function(i) {
      result_29(step_2[[i]], y_c_names[i])
    })
    
    graphs[['g_29']] <- do.call(grid.arrange, c(graph_list, nrow = 1))
  }
  
  # g_ij(s,t) at 1_1
  if('31' %in% graph_ids){ 
    
    # g_list <- lapply(step_3, function(x){result_30s_prep(x, i = 1, j = 1, full)})
    # graphs[['g_31']] <- rearrange_plots(g_list)
    
    graphs[['g_31']] <- result_heatmap_ij_prep(step_3, 'g_ij', i = 1, j = 1, zmid = 0)
  }  
  
  # g_ij(s,t) at 1_2
  if('32' %in% graph_ids){ 
    graphs[['g_32']] <- result_heatmap_ij_prep(step_3, 'g_ij', i = 1, j = 2, zmid = 0)
  }  
  
  # g_ij(s,t) at 3_5
  if('33' %in% graph_ids){ 
    graphs[['g_33']] <- result_heatmap_ij_prep(step_3, 'g_ij', i = 3, j = 5, zmid = 0)
  }    
  
  # eigenfunctions
  if('41' %in% graph_ids){ 
    # g_list <- lapply(step_4, function(x){result_41_prep(x, time_grid, time_grid_est, full)})
    # graphs[['g_41']] <- rearrange_plots(g_list)
    
    # only keep eigenfunctions from first process
    
    step_4_v2 <- lapply(step_4, function(x) {
      list(
        eigen_decomp_est   = x$eigen_decomp_est$eigenfunctions[[1]] %>% t(),
        eigen_decomp_truth = x$eigen_decomp_truth$eigenfunctions[[1]] %>% t()
      )
    })   
    
    graphs[['g_41']] <- result_line_graph_prep(step_4_v2, 'eigen_decomp', time_grid_est)
  }  
  
  # reconstructing g_ij from eigenfunctions
  if('42' %in% graph_ids){ 
    g_list <- lapply(1:length(step_3), function(x){result_42_prep(step_3[[i]], step_4[[i]], p, full)})
    graphs[['g_42']] <- rearrange_plots(g_list)
  }  
  
  # orthonormality of eigenfunctions
  if('43' %in% graph_ids){ 
    g_list <- lapply(step_4, function(x){result_43_prep(x, full)})
    graphs[['g_43']] <- rearrange_plots(g_list)
  }
  
  # reconstruction error histogram 
  if('44' %in% graph_ids){ 
    g_list <- lapply(1:length(step_3), function(x){result_44_prep(step_3[[i]], step_4[[i]], p, full)})
    graphs[['g_44']] <- rearrange_plots(g_list)
  } 
  # eigenvalue decay - check if G_ii / m makes the eigenvalues similar between coarse/fine settings
  if('45' %in% graph_ids){
    g_list <- lapply(step_4, function(x) result_45_prep(x))
    graphs[['g_45']] <- rearrange_plots(g_list) 
  }
  
  # KL_cov of (1, 2) block
  if('55' %in% graph_ids){
    g_list <- lapply(step_5, function(x) result_55_prep(x, i = 1, j = 2, full))
    graphs[['g_55']] <- rearrange_plots(g_list) 
  }
  
  # KL_cor of (1, 2) block
  if('56' %in% graph_ids){
    g_list <- lapply(step_5b, function(x) result_56_prep(x, i = 1, j = 2, full))
    graphs[['g_56']] <- rearrange_plots(g_list) 
  }
  
  if('81' %in% graph_ids & 'step_8' %in% names(graph_results_i)){
    est_graphs <- lapply(graph_results_i$step_8, function(x) extract_block_structure_ij(x$step_8$V_cond_est_full, m_est, i = 1, j = 2))
  }
  
  # C_Xi_Xj for block (1, 1)
  
  if('90' %in% graph_ids){
    graphs[['g_90']] <- result_heatmap_ij_prep(step_9, 'C_cond', time_grid_est, is_full, i = 1, j = 1, zmid = 0)
  }  
  # C_Xi_Xj for block (1, 2)
  
  if('91' %in% graph_ids){
    # g_list <- lapply(graph_results_i$step_9, function(x){ result_90s_prep_ij(x, m, m_est, i = 1, j = 2, full) })
    # 
    # graphs[['g_91']] <- rearrange_plots(g_list) 
    # 
    # graphs[['g_91']] <- result_heatmap_ij_prep(step_9, 'C_cond', T, i = 1, j = 1)
    graphs[['g_91']] <- result_heatmap_ij_prep(step_9, 'C_cond', time_grid_est, is_full, i = 1, j = 2, zmid = 0)
  }
  
  # C_Xi_Xj for block (3, 5)
  
  if('92' %in% graph_ids){
    # graph_list <- lapply(graph_results_i$step_9, function(x){ result_90s_prep_ij(x, m, m_est, i = 3, j = 5, full) })
    # 
    # n_y_c_query <- length(graph_list)
    # n_settings <- length(graph_list[[1]])
    # 
    # # Flatten the nested list: row-wise
    # flat_graphs <- unlist(graph_list, recursive = FALSE)
    # 
    # # Create column-major index mapping
    # # R's matrix() fills column-wise by default, so we transpose to reorder properly
    # idx <- as.vector(t(matrix(seq_along(flat_graphs), nrow = n_settings, ncol = n_y_c_query)))
    # 
    # # Reorder the flat list
    # flat_graphs_colwise <- flat_graphs[idx]
    # 
    # # Arrange in n_settings rows x n_y_c_query columns
    # graphs[['g_92']] <- do.call(grid.arrange, c(flat_graphs_colwise, nrow = n_settings, ncol = n_y_c_query))
    
    graphs[['g_92']] <- result_heatmap_ij_prep(step_9, 'C_cond', time_grid_est, is_full, i = 3, j = 5, zmid = 0)
    
  }  
  
  # C for all blocks (pm x pm)
  if('93' %in% graph_ids){
    g_list <- lapply(graph_results_i$step_9, function(x) result_90s_prep_pm(x, m, m_est, full))
    graphs[['g_93']] <- rearrange_plots(g_list) 
  }  
  
  # C_HS for all pxp blocks
  if('95' %in% graph_ids){
    g_list <- lapply(graph_results_i$step_11, function(x) result_95_prep(x, m, m_est, p, full))
    graphs[['g_95']] <- rearrange_plots(g_list) 
  }
  
  # P for all blocks (pm x pm)
  if('103' %in% graph_ids){
    g_list <- lapply(graph_results_i$step_10, function(x) result_100s_prep_pm(x, m, m_est, full))
    graphs[['g_103']] <- rearrange_plots(g_list) 
  } 
  
  # P_HS for all pxp blocks
  if('112' %in% graph_ids){
    g_list <- lapply(graph_results_i$step_11, function(x) result_112_prep(x, remove_diag = T, full))
    graphs[['g_112']] <- rearrange_plots(g_list) 
  }  
  
  if('113' %in% graph_ids){
    
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
            title = paste("ROC:", i, "-", j),
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
    graphs[['113']] <- wrap_plots(all_plots, ncol = 7, nrow = 10) 
    
  }
  
  return(graphs)
  
  
}









