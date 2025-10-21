source('functions/28z_Visualization_helpers.R')
source('functions/28y_Visualization_Blocks.R')
source('functions/28_Simulation_Visualization.R')

library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(purrr)

# visualize results across all available sample sizes and y_c

Visualize_over_n_and_time <- function(folder_name, graph_ids, m_est, m, i, j){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: for all n and y_c values, plot a specific intermediate value
  #
  # 
  # input: 
  #
  # - folder_name   (string)   "simu_results/single_c2/CPGM"
  # - graph_ids     (vector)   vector of graph ID's that we want
  # - m_est
  # - m
  # - i
  # - j
  # 
  # output:
  #
  #
  # 
  # ----------------------------------------------------------------------------   
  
  # prep - load all datasets
  
  results_list <- load_all_results(folder_name)
  

  arr_mat <- matrix(1:6, nrow = 2, byrow = F)
  arr_mat_8 <- matrix(1:8, nrow = 2, byrow = F)
  arr_mat_10 <- matrix(1:10, nrow = 2, byrow = F)  
  graphs <- list()

  
  # C_Xi_Xj
  if('91' %in% graph_ids){
    
    # 1) get a list of all C_cond_est_fulls
    
    
    C_cond_list <- map(results_list, function(n_list) {
      step9 <- n_list$graph_results_i$step_9
      map(step9, ~ .x$C_cond_est_full)
    })    
    
    # 2) extract the i, j-th part of them, plot them all in a heatmap, then arrange them
    i <- 1
    j <- 2
    plots_nested <- imap(C_cond_list, function(sublist, n_name) {
      imap(sublist, function(matrix, q_index) {
        mat_block <- extract_block_structure_ij(matrix, m_est, i, j)
        visualize_matrix_heatmap(mat_block,
                                 paste0("n=", n_name, ", q=", q_index),
                                 zmid = 0)
      })
    })
    
    # Step 3: Flatten to a single list of plots
    plots_flat <- flatten(plots_nested)
    
    # Step 4: Arrange in a 4 x q grid
    n_sizes <- length(results_list)
    q <- length(results_list[[1]]$graph_results_i$step_9)
    
    graphs[['g_91']] <- grid.arrange(grobs = plots_flat, nrow = n_sizes, ncol = q)
  }
  
}