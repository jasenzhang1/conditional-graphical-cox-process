source('functions/28_Simulation_Visualization.R')
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


# everything over all query_id

visualize_over_time <- function(folder_name, n, est_only, graph_id, m_est, m, i, j){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: for a single dataset, plot a specific intermediate value over time
  #
  # - this is to be done for a single simulation (adj_method, est_method, n)
  # 
  # input: 
  #
  # - folder_name   (string)
  # - n             (integer)
  # - est_only      (boolean)  if true, only display est over time with ground truth 
  # 
  # output:
  #
  #
  # 
  # ---------------------------------------------------------------------------- 
  
  
  # 1) load
  
  graph_results_i <- get_file_name_and_load(folder_name, n) 
  graph_results_i <- old_to_new_graph_results_i(graph_results_i)
  
  # prep
  arr_mat <- matrix(1:6, nrow = 2, byrow = F)
  arr_mat_8 <- matrix(1:8, nrow = 2, byrow = F)
  arr_mat_10 <- matrix(1:10, nrow = 2, byrow = F)  
  
  
  graphs <- list()
  
  y_c_names <- names(graph_results_i$kernel_params_i) %>% as.numeric() %>% round(2)
  
  # 2) figure out which steps have unique results for y_c_query and which steps are constant
  
  if('11' %in% graph_ids){
    #   - X_k_est                   (p x m_est   x n)
    #   - X_k_truth                 (p x m_truth x n)
    #   - X_k_coarse_truth          (p x m_est   x n)
    #   - X_k_both_truth            (p x m_both  x n)
    #   - Lambda_k_truth            (p x m_truth x n)
    #   - Lambda_k_coarse_truth     (p x m_est   x n)
    
    step_1 <- graph_results_i$step_1
    
    mu_t <- graph_results_i$kernel_params_i[[1]]$base_mean
    mu_t_est <- graph_results_i$kernel_params_i[[1]]$base_mean_est
    mu_t_both <- graph_results_i$kernel_params_i[[1]]$base_mean_both
    
    graphs[['g_11']] <- grid.arrange(visualize_log_intensity(step_1[[1]][1:5,,1],   time_grid_est,  'Estimate', mu_t_est ),
                                     visualize_log_intensity(step_1[[3]][1:5,,1],   time_grid_est,  'Coarser Truth', mu_t_est ),
                                     visualize_log_intensity(step_1[[2]][1:5,,1],   time_grid,      'Finer Truth', mu_t     ),
                                     visualize_log_intensity(step_1[[4]][1:5,,1],   time_grid_both, 'Combined Truth', mu_t_both),
                                     textGrob("0. Log Intensity\n of first 5 processes\nof subject 1", gp = gpar(fontsize = 14)),
                                     layout_matrix = arr_mat) 
  }  
  
  if('81' %in% graph_ids & 'step_8' %in% names(graph_results_i)){
    est_graphs <- lapply(graph_results_i$step_8, function(x) extract_block_structure_ij(x$step_8$V_cond_est_full, m_est, i, j))
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









