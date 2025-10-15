
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
  
  # 2) figure out which steps have unique results for y_c_query and which steps are constant
  
  step_const <- names(graph_results_i$estimated_graphs_part_1)
  step_vary <- names(graph_results_i$estimated_graphs_part_2[[1]])
  y_c_query <- names(graph_results_i$estimated_graphs_part_2)
  
  if('81' %in% graph_ids & 'step_8' %in% step_vary){
    est_graphs <- lapply(graph_results_i$estimated_graphs_part_2, function(x) extract_block_structure_ij(x$step_8$V_cond_est_full, m_est, i, j))
  }
  
  
}









