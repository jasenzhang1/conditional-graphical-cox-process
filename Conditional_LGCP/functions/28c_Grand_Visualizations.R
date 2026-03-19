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
rearrange_plots <- function(g_list, main_title = NULL){
  
  n_time <- length(g_list)
  n_graphs <- length(g_list[[1]])
  
  # 1. Unlist and reorganize plots
  plots <- unlist(g_list, recursive = FALSE)
  indices <- as.vector(sapply(1:n_graphs, function(r) r + (0:(n_time-1))*n_graphs))
  plots_colwise <- plots[indices]
  
  # 2. Setup the arguments for grid.arrange
  # We use grobs = plots_colwise to explicitly separate the plots from the layout settings
  args_list <- list(
    grobs = plots_colwise,
    nrow = n_graphs,
    ncol = n_time
  )
  
  # 3. Add the title only if it exists
  if(!is.null(main_title)) {
    args_list$top <- grid::textGrob(main_title, gp = grid::gpar(fontsize = 16, fontface = "bold"))
  }
  
  # 4. Call grid.arrange using the cleaned argument list
  arranged_plot <- do.call(gridExtra::grid.arrange, args_list)
  
  return(arranged_plot)
}

# helper 
label_model_type <- function(vec) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: rename items to whether they are the truth or how they are thresholded
  # 
  #       valid names are 'truth', 'none', 'local', 'hybrid', 'global' to describe the thresholding
  #
  #       example: "C_HS_truth"                  "C_HS_KL_est_eig1"            "C_HS_KL_GIC_local_est_eig1"  "C_HS_KL_GIC_hybrid_est_eig1"
  #       turn this into 'truth', 'none', 'local', and 'hybrid'
  #
  # inputs:
  #
  # - vec   (vector of strings)
  #
  # - results (vector of renamed strings)
  #
  # ----------------------------------------------------------------------------
  
  # Initialize with "none" as the default
  results <- rep("none", length(vec))
  
  # Update based on keyword matches
  # Note: The order here determines priority if a string has multiple keywords
  results[grepl("truth",  vec)] <- "truth"
  results[grepl("local",  vec)] <- "local"
  results[grepl("hybrid", vec)] <- "hybrid"
  results[grepl("global", vec)] <- "global"
  
  return(results)
}

# everything over all query_id

visualize_over_time <- function(graph_results_i, graph_ids, ground_truth, beta_truth, X_truth, eigen_setting){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: for a single dataset, plot a specific intermediate value over time
  #
  # - this is to be done for a single simulation (adj_method, est_method, n)
  # 
  # input: 
  #
  # - graph_results_i   (list of all step_X's)  --> (list of all y_c_queries) --> (list of 'est' or 'truth' or 'coarse_truth' etc)
  #
  #   - step_2, 2b, 3, 4, 5, 5b, 5c, 5d, 9, 9b, 10, 11, 11b, 12, 12b
  #   - step_1, 1b, 1c, step_0_events
  #   - true_graphs
  #   - y_c_query
  #   - p
  #   - Y_continuous
  #   - weights
  #   - W_y
  #   - time_grid, time_grid_est, time_grid_both
  #
  # - graph_ids             (vector of strings)   which graphs do we want?
  # - ground_truth          (boolean)             do our results have the overall truth?
  # - beta_truth            (boolean)             do our results have beta_truth values? 
  # - X_truth               (boolean)             do our results have X_truth values? 
  # - eigen_setting         (string)              only_joint, trig_and_joint, trig_simple
  # 
  # output:
  #
  # - graphs   (list)  list of graphs
  # 
  # ---------------------------------------------------------------------------- 
  
  
  # 1) load everything from graph_results_i
  
  # graph_results_i also contains time_grid, time_grid_est, y_c_query
  list2env(graph_results_i, envir = environment())  # does all of the above in one go
  
  
  # If time_grid doesn't exist in the current environment, set it to NA
  if (!exists("time_grid", inherits = FALSE)) {
    time_grid <- NA
  }
  
  # prep
  arr_mat_6 <- matrix(1:6, nrow = 2, byrow = F)
  arr_mat_8 <- matrix(1:8, nrow = 2, byrow = F)
  arr_mat_10 <- matrix(1:10, nrow = 2, byrow = F)  
  
  m_est <- length(time_grid_est)
  
  graphs <- list()
  
  y_c_names <- graph_results_i$y_c_query %>%
    as.numeric() %>%
    round(3) %>%
    formatC(format = "f", digits = 3)
  
  # 2) figure out which steps have unique results for y_c_query and which steps are constant
  
  # raw data for subject k
  if('01' %in% graph_ids){
    k <- 1
    graphs[['g_01']] <- visualize_points_step_0_events(step_0_events, k) #28b
  }
  
  
  # histograms of event times with intensity (Lambda) overlay
  if('02' %in% graph_ids){
    k <- 1
    if(is.na(time_grid)){
      graphs[['g_02']] <- visualize_points_on_intensity(step_0_events, step_1b, k, time_grid_est, F, time_grid) #28b
    } else{
      graphs[['g_02']] <- visualize_points_on_intensity(step_0_events, step_1b, k, time_grid_est, T, time_grid) #28b
    }
    
  }
  
  # intensities for each subject across all processes
  if('03' %in% graph_ids){
    g_03_est <- visualize_log_intensity_bold_mean(step_1b$Lambda_k_est, time_grid_est)
    g_03_X_truth <- visualize_log_intensity_bold_mean(step_1b$Lambda_k_truth, time_grid)
  }
  
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
      
      step_1_prep <- lapply(step_1, function(x) x[1:5, , i])
      step_1_prep <- list(y = step_1_prep)
      
      graphs[['g_11b']] <- result_line_graph_prep(step_1_prep, 'X_k', time_grid_est, grouping = 'process')
      
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
                                       visualize_log_intensity(apply(step_1$X_k_truth,          c(1, 2), mean)[1:5, ],   time_grid,      'Finer Truth',     step_1c$mu_t_truth),
                                       visualize_log_intensity(apply(step_1$X_k_both_truth,     c(1, 2), mean)[1:5, ],   time_grid_both, 'Combined Truth',  step_1c$mu_t_both_truth),
                                       textGrob("0. Average Log Intensity\n of first 5 processes", gp = gpar(fontsize = 14)),
                                       layout_matrix = arr_mat_6) 
    } else{
      graphs[['g_12']] <- grid.arrange(visualize_log_intensity(apply(step_1$X_k_est,          c(1, 2), mean)[1:5, ],   time_grid_est,  'Estimate',     step_1c$mu_t_coarse_truth),
                                       textGrob("0. Average Log Intensity\n of first 5 processes", gp = gpar(fontsize = 14)),
                                       layout_matrix = arr_mat_6) 
    }
  }
  
  # rho_i(t) for processes 1 through 5 - grouped by estimand (so all processes are together)
  if('22' %in% graph_ids){
    graphs[['g_22']] <- result_line_graph_prep(step_2, 'rho_i', time_grid_est, grouping = 'estimand', num_processes = 5)
  }
  
  # rho_i(t) for processes 1 through 5 - grouped by process (so all different estimands are together)
  if('23' %in% graph_ids){
    graphs[['g_23']] <- result_line_graph_prep(step_2, 'rho_i', time_grid_est, grouping = 'process', num_processes = 12)
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
    graphs[['g_29']] <- result_29(weights, Y_continuous, y_c_names)
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
    g_list <- lapply(1:length(step_3), function(i){result_44_prep(step_3[[i]], step_4[[i]], p, ground_truth, X_truth, eigen_setting)}) 
    names(g_list) <- y_c_names
    
    graphs[['g_44']] <- result_histogram_prep(g_list, 'error', data_format = 'regular', nbins = 20)
  } 
  # eigenvalue decay - check if G_ii / m makes the eigenvalues similar between coarse/fine settings
  if('45' %in% graph_ids){
    g_list <- lapply(step_4, function(x) result_45_prep(x))
    graphs[['g_45']] <- rearrange_plots(g_list) 
  }
  
  # orthogonality - check if phi_a %*% G_ii %*% phi_b = 0 for each process
  if('46' %in% graph_ids){
    g_list <- lapply(1:length(step_3), function(i){
      result_46_prep(step_3[[i]], step_4[[i]])
    })
    graphs[['g_45']] <- rearrange_plots(g_list) 
  }
  
  # assume we know beta_truths, do their empirical correlations match the actual correlations?
  if('50' %in% graph_ids){
    
    graphs[['g_50']] <- lapply(1:length(step_5d), function(i){
      visualize_beta_corr(step_5d[[i]]$KL_coeffs_truth, 
                          true_graphs[[i]]$cov_mat_truth,
                          true_graphs[[i]]$cor_mat_truth,
                          true_graphs[[i]]$prec_mat_truth,
                          graph_type = 'heatmap')
    })
    names(graphs[['g_50']]) <- y_c_names
  }
  
  # assume we know beta_truths, do their empirical correlations match the actual correlations? With a histogram
  if('51' %in% graph_ids){
    
    graphs[['g_51']] <- lapply(1:length(step_5d), function(i){
      visualize_beta_corr(step_5d[[i]]$KL_coeffs_truth, 
                          true_graphs[[i]]$cov_mat_truth,
                          true_graphs[[i]]$cor_mat_truth,
                          true_graphs[[i]]$prec_mat_truth,
                          graph_type = 'histogram')
    })
    
    names(graphs[['g_51']]) <- y_c_names
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
  
  
  
  # KL_cor entire matrix + KL_cov 
  if('58' %in% graph_ids){
    g_list <- lapply(step_5b, function(x) result_58_prep(x, p))   # correlation
    graphs[['g_58']] <- rearrange_plots(g_list, main_title = 'KL Correlation') 
    
    g_list <- lapply(step_5, function(x) result_58_cov_prep(x, p))  # covariance
    graphs[['g_58b']] <- rearrange_plots(g_list, main_title = 'KL Covariance') 
    
    
  }
  
  # KL_prec entire matrix
  if('59' %in% graph_ids){
    g_list <- lapply(step_5c, function(x) result_59_prep(x, p))   # precision
    graphs[['g_59']] <- rearrange_plots(g_list, main_title = 'KL Precision') 
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
    graphs[['g_93']] <- result_heatmap_nonblock_prep(step_9, 'C_cond', data_format = 'list', time_grid_est = time_grid_est, zmid = 0)
  }  
  
  # C_HS for all pxp blocks
  if ('95' %in% graph_ids) {
    
    # 1. Generate the base heatmaps using the raw step_11b data
    # Note: Keeping this outside the filter logic to ensure g_94/g_95 always exist
    graphs[['g_94']] <- result_heatmap_nonblock_prep(step_11b, 'C_HS', data_format = 'regular', 
                                                     time_grid_est = time_grid_est, rm_diag = F, 
                                                     graph_title = 'C_HS values (display diag)', zmid = 0)
    graphs[['g_95']] <- result_heatmap_nonblock_prep(step_11b, 'C_HS', data_format = 'regular', 
                                                     time_grid_est = time_grid_est, rm_diag = T, 
                                                     graph_title = 'C_HS values (hide diag)', zmid = 0)
    
    # 2. Define the regex pattern based on eigen_setting
    # trig_simple   -> eig1
    # trig_and_joint -> eig1, eig2, or eig3
    # only_joint    -> eig3
    suffix_pattern <- switch(eigen_setting,
                             "trig_simple"    = "est_eig1$",
                             "trig_and_joint" = "est_eig[123]$",
                             "only_joint"     = "est_eig3$",
                             "est_eig[123]$"  # Default fallback
    )
    
    # 3. Identify and filter indices
    step_11b_names <- names(step_11b[[1]])
    suffixes       <- sub("^C_HS_", "", step_11b_names)
    
    target_indices <- grep(suffix_pattern, suffixes)
    truth_index    <- which(suffixes == 'truth')
    kept_indices   <- sort(unique(c(truth_index, target_indices)))
    
    # 4. Map names and transform the list
    kept_names <- step_11b_names[kept_indices]
    new_names  <- label_model_type(kept_names)
    name_map   <- setNames(new_names, kept_names)
    
    step_11b_v2 <- lapply(step_11b, function(x) {
      x_filtered <- x[kept_indices]
      names(x_filtered) <- name_map[names(x_filtered)]
      return(x_filtered)
    })
    
    # 5. Generate the filtered heatmaps (g_94b / g_95b)
    title_suffix <- paste("for", eigen_setting)
    graphs[['g_94b']] <- result_heatmap_nonblock_prep(step_11b_v2, NULL, data_format = 'regular', 
                                                      time_grid_est = time_grid_est, rm_diag = F, 
                                                      graph_title = paste('C_HS values', title_suffix, '(display diag)'), zmid = 0)
    graphs[['g_95b']] <- result_heatmap_nonblock_prep(step_11b_v2, NULL, data_format = 'regular', 
                                                      time_grid_est = time_grid_est, rm_diag = T, 
                                                      graph_title = paste('C_HS values', title_suffix, '(hide diag)'), zmid = 0)
    
    if('none' %in% names(step_11b_v2[[1]])){
      step_11b_v3 <- lapply(step_11b_v2, function(sublist) {
        sublist$none <- NULL
        return(sublist)
      })
      
      title_suffix <- paste("for", eigen_setting)
      graphs[['g_94c']] <- result_heatmap_nonblock_prep(step_11b_v3, NULL, data_format = 'regular', 
                                                        time_grid_est = time_grid_est, rm_diag = F, 
                                                        graph_title = paste('C_HS values', title_suffix, '(display diag)'), zmid = 0)
      graphs[['g_95c']] <- result_heatmap_nonblock_prep(step_11b_v3, NULL, data_format = 'regular', 
                                                        time_grid_est = time_grid_est, rm_diag = T, 
                                                        graph_title = paste('C_HS values', title_suffix, '(hide diag)'), zmid = 0)
      
    }
    
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
    graphs[['g_103']] <- result_heatmap_nonblock_prep(step_10, 'P_cond', data_format = 'list', time_grid_est = time_grid_est, rm_diag = T, zmid = 0)
  } 
  
  
  
  # P_HS (w_mat) for all pxp blocks
  # 110 = don't remove diagonals
  # 111 = remove diagonals
  
  if ('111' %in% graph_ids) {
    
    # 1. Define the regex pattern for filtering based on eigen_setting
    suffix_pattern <- switch(eigen_setting,
                             "trig_simple"    = "est_eig1$",
                             "trig_and_joint" = "est_eig[123]$",
                             "only_joint"     = "est_eig3$",
                             "est_eig[123]$"  # Default fallback
    )
    
    # 2. Generate base heatmaps using raw step_11 data
    graphs[['g_110']] <- result_heatmap_nonblock_prep(step_11, 'w_mat', data_format = 'regular', 
                                                      time_grid_est = time_grid_est, rm_diag = F, 
                                                      graph_title = 'P_HS values (display diag)', zmid = 0)
    graphs[['g_111']] <- result_heatmap_nonblock_prep(step_11, 'w_mat', data_format = 'regular', 
                                                      time_grid_est = time_grid_est, rm_diag = T, 
                                                      graph_title = 'P_HS values (hide diag)', zmid = 0)
    
    # 3. Identify indices and map names
    step_11_names <- names(step_11[[1]])
    suffixes      <- sub("^w_mat_", "", step_11_names)
    
    target_indices <- grep(suffix_pattern, suffixes)
    truth_index    <- which(suffixes == 'truth')
    kept_indices   <- sort(unique(c(truth_index, target_indices)))
    
    kept_names <- step_11_names[kept_indices]
    new_names  <- label_model_type(kept_names)
    name_map   <- setNames(new_names, kept_names)
    
    # 4. Create the filtered and renamed list (step_11_v2)
    step_11_v2 <- lapply(step_11, function(x) {
      x_filtered <- x[kept_indices]
      names(x_filtered) <- name_map[names(x_filtered)]
      return(x_filtered)
    })
    
    # 5. Generate the "b" version graphs (and handle 'none' filtering if needed)
    # We check if 'none' is in the mapped names of the FIRST element
    final_v2_names <- names(step_11_v2[[1]])
    
    
    title_tag <- paste("for", eigen_setting)
    graphs[['g_110b']] <- result_heatmap_nonblock_prep(step_11_v2, NULL, data_format = 'regular', 
                                                       time_grid_est = time_grid_est, rm_diag = F, 
                                                       graph_title = paste('P_HS values', title_tag, '(display diag)'), zmid = 0)
    graphs[['g_111b']] <- result_heatmap_nonblock_prep(step_11_v2, NULL, data_format = 'regular', 
                                                       time_grid_est = time_grid_est, rm_diag = T, 
                                                       graph_title = paste('P_HS values', title_tag, '(hide diag)'), zmid = 0)
    
    if('none' %in% names(step_11_v2[[1]])){
      step_11_v3 <- lapply(step_11_v2, function(sublist) {
        sublist$none <- NULL
        return(sublist)
      })
      
      title_suffix <- paste("for", eigen_setting)
      graphs[['g_110c']] <- result_heatmap_nonblock_prep(step_11_v3, NULL, data_format = 'regular', 
                                                         time_grid_est = time_grid_est, rm_diag = F, 
                                                         graph_title = paste('P_HS values', title_tag, '(display diag)'), zmid = 0)
      graphs[['g_111c']] <- result_heatmap_nonblock_prep(step_11_v3, NULL, data_format = 'regular', 
                                                         time_grid_est = time_grid_est, rm_diag = T, 
                                                         graph_title = paste('P_HS values', title_tag, '(hide diag)'), zmid = 0)
      
    }
    
  }
  
  # P_HS (w_mat) for all pxp blocks - delete est and X_truth
  if('112' %in% graph_ids){
    
    # groups
    # 1) truth_names     -  truths
    # 2) GIC_names       -  KL_GIC est and X_truth
    # 3) GIC_est_names   -  KL_GIC est 
    # 4) est_names       -  est and X_truth
    # 5) KL_est_names    -  just est
    
    step_11_names <- names(step_11[[1]])
    truth_names   <- step_11_names[grepl("truth$", step_11_names)]
    step_11_names <- setdiff(step_11_names, truth_names)
    
    
    KL_GIC_names  <- step_11_names[grepl("^w_mat_KL_GIC", step_11_names)]  
    KL_GIC_est_names <- grep("_est_", KL_GIC_names, value = TRUE)
    est_names <- setdiff(step_11_names, KL_GIC_names)
    KL_est_names <- est_names[startsWith(est_names, "w_mat_KL_est")]
    
    
    # truths
    step_11_v2 <- step_11
    step_11_v2 <- lapply(step_11_v2, function(x) {
      x[names(x) %in% truth_names]
    })
    graphs[['g_112b']] <- result_heatmap_nonblock_prep(step_11_v2, 'w_mat', data_format = 'regular', time_grid_est = time_grid_est, rm_diag = T, graph_title = 'P_HS values of truths', zmid = 0)
    
    # estimates
    step_11_v3 <- step_11
    step_11_v3 <- lapply(step_11_v3, function(x) {
      x[names(x) %in% est_names]
    })
    graphs[['g_112c']] <- result_heatmap_nonblock_prep(step_11_v3, 'w_mat', data_format = 'regular', time_grid_est = time_grid_est, rm_diag = T, graph_title = 'P_HS values of non-thresholded results', zmid = 0)
    
    # KL_est
    step_11_v6 <- step_11
    step_11_v6 <- lapply(step_11_v6, function(x) {
      x[names(x) %in% KL_est_names]
    })
    graphs[['g_112d']] <- result_heatmap_nonblock_prep(step_11_v6, 'w_mat', data_format = 'regular', time_grid_est = time_grid_est, rm_diag = T, graph_title = 'P_HS values of non-thresholded estimates', zmid = 0)
    
    # KL_GIC
    step_11_v4 <- step_11
    step_11_v4 <- lapply(step_11_v4, function(x) {
      x[names(x) %in% KL_GIC_names]
    })
    graphs[['g_112e']] <- result_heatmap_nonblock_prep(step_11_v4, 'w_mat', data_format = 'regular', time_grid_est = time_grid_est, rm_diag = T, graph_title = 'P_HS values of thresholded results', zmid = 0)
    
    # KL_GIC_est
    step_11_v5 <- step_11
    step_11_v5 <- lapply(step_11_v5, function(x) {
      x[names(x) %in% KL_GIC_est_names]
    })
    graphs[['g_112f']] <- result_heatmap_nonblock_prep(step_11_v5, 'w_mat', data_format = 'regular', time_grid_est = time_grid_est, rm_diag = T, graph_title = 'P_HS values of thresholded estimates', zmid = 0)
    
    
    
  } 
  
  # ROC curves
  if('113' %in% graph_ids){
    
    
    # groups
    # 1) truth_names     -  truths
    # 2) GIC_names       -  KL_GIC est and X_truth
    # 3) GIC_est_names   -  KL_GIC est 
    # 4) est_names       -  est and X_truth
    # 5) KL_est_names    -  just est
    
    step_12_names <- names(step_12[[1]])
    truth_names   <- step_12_names[grepl("truth$", step_12_names)]
    step_12_names <- setdiff(step_12_names, truth_names)
    
    
    KL_GIC_names  <- step_12_names[grepl("^roc_KL_GIC", step_12_names)]  
    KL_GIC_est_names <- grep("_est_", KL_GIC_names, value = TRUE)
    est_names <- setdiff(step_12_names, KL_GIC_names)
    KL_est_names <- est_names[startsWith(est_names, "roc_KL_est")]
    
    
    graphs[['g_113']] <- result_ROC_prep(step_12, NULL, NULL)                 # everything
    graphs[['g_113c']] <- result_ROC_prep(step_12, est_names, NULL)           # all non-threshold results
    graphs[['g_113d']] <- result_ROC_prep(step_12, KL_est_names, NULL)        # all non-threshold estimates
    graphs[['g_113e']] <- result_ROC_prep(step_12, KL_GIC_names, NULL)        # all threshold results
    graphs[['g_113f']] <- result_ROC_prep(step_12, KL_GIC_est_names, NULL)    # all threshold estimates
    
  }
  
  # Final adj_mat for all pxp blocks
  if ('114' %in% graph_ids) {
    
    # 1. Map eigen_setting to the correct regex pattern
    suffix_pattern <- switch(eigen_setting,
                             "trig_simple"    = "est_eig1$",
                             "trig_and_joint" = "est_eig[123]$",
                             "only_joint"     = "est_eig3$",
                             "est_eig[123]$"  # Default fallback
    )
    
    # 2. Base heatmaps using raw step_12b data (adj_mat prefix)
    graphs[['g_114']] <- result_heatmap_nonblock_prep(step_12b, 'adj_mat', data_format = 'regular', 
                                                      time_grid_est = time_grid_est, rm_diag = T, 
                                                      graph_title = 'edge sets', zmid = 0)
    
    # 3. Identify indices to keep (including truth and none)
    step_12b_names <- names(step_12b[[1]])
    suffixes       <- sub("^adj_mat_", "", step_12b_names)
    
    target_indices <- grep(suffix_pattern, suffixes)
    truth_indices  <- which(grepl("truth$", suffixes))
    none_indices   <- which(suffixes == 'none')
    
    kept_indices   <- sort(unique(c(truth_indices, target_indices, none_indices)))
    
    # 4. Map names and transform the list for the 'b' versions
    kept_names <- step_12b_names[kept_indices]
    new_names  <- label_model_type(kept_names)
    name_map   <- setNames(new_names, kept_names)
    
    step_12b_v2 <- lapply(step_12b, function(x) {
      x_filtered <- x[kept_indices]
      names(x_filtered) <- name_map[names(x_filtered)]
      return(x_filtered)
    })
    
    # 5. Generate the sub-group heatmaps following the original logic but filtered
    # Use the mapped names in step_12b_v2 to define subsets
    v2_names <- names(step_12b_v2[[1]])
    
    # b) Truths
    truth_v2_names <- v2_names[grepl("truth", v2_names)]
    
    if(length(truth_v2_names) > 0){
      step_12b_v2_truths <- lapply(step_12b_v2, function(x) x[names(x) %in% truth_v2_names])
      graphs[['g_114b']] <- result_heatmap_nonblock_prep(step_12b_v2_truths, NULL, data_format = 'regular', 
                                                         time_grid_est = time_grid_est, rm_diag = T, 
                                                         graph_title = 'edge set of truths', zmid = 0)
    }
    
    
    # c) Estimates (non-GIC)
    # Based on your original code's logic for est_names
    est_v2_names <- v2_names[!grepl("truth", v2_names) & !grepl("GIC", v2_names)]
    step_12b_v2_est <- lapply(step_12b_v2, function(x) x[names(x) %in% est_v2_names])
    graphs[['g_114c']] <- result_heatmap_nonblock_prep(step_12b_v2_est, NULL, data_format = 'regular', 
                                                       time_grid_est = time_grid_est, rm_diag = T, 
                                                       graph_title = 'edge set of non-thresholded results', zmid = 0)
    
  }
  
  # Final adj_mat for local, hybrid, global with x's on blocks that are false
  if ('115' %in% graph_ids){
    
    # 1. Map eigen_setting to the correct regex pattern
    suffix_pattern <- switch(eigen_setting,
                             "trig_simple"    = "est_eig1$",
                             "trig_and_joint" = "est_eig[123]$",
                             "only_joint"     = "est_eig3$",
                             "est_eig[123]$"  # Default fallback
    )
    

    
    # 3. Identify indices to keep (including truth and none)
    step_12b_names <- names(step_12b[[1]])
    suffixes       <- sub("^adj_mat_", "", step_12b_names)
    
    target_indices <- grep(suffix_pattern, suffixes)
    truth_indices  <- which(grepl("truth$", suffixes))

    kept_indices   <- sort(unique(c(truth_indices, target_indices)))
    
    # 4. Map names and transform the list for the 'b' versions
    kept_names <- step_12b_names[kept_indices]
    new_names  <- label_model_type(kept_names)
    name_map   <- setNames(new_names, kept_names)
    
    step_12b_v2 <- lapply(step_12b, function(x) {
      x_filtered <- x[kept_indices]
      names(x_filtered) <- name_map[names(x_filtered)]
      return(x_filtered)
    })
    
    # 5. Generate the sub-group heatmaps following the original logic but filtered
    # Use the mapped names in step_12b_v2 to define subsets
    v2_names <- names(step_12b_v2[[1]])
    
    
    # c) Estimates (non-GIC)
    # Based on your original code's logic for est_names
    est_v2_names <- v2_names[!grepl("none", v2_names)]
    step_12b_v3 <- lapply(step_12b_v2, function(x) x[names(x) %in% est_v2_names])
    graphs[['g_115']] <- result_heatmap_mismatch_x(step_12b_v3, 'truth', 1, rm_diag = T, zmid = 0)
    

  }
  
  # accuracy/f1/sens/spec/ppv/npv
  if('120' %in% graph_ids){
    
    step_12_names <- names(step_12[[1]])
    truth_names   <- step_12_names[grepl("truth$", step_12_names)]
    step_12_names <- setdiff(step_12_names, truth_names)
    
    
    KL_GIC_names  <- step_12_names[grepl("^roc_KL_GIC", step_12_names)]  
    KL_GIC_est_names <- grep("_est_", KL_GIC_names, value = TRUE)
    est_names <- setdiff(step_12_names, KL_GIC_names)
    KL_est_names <- est_names[startsWith(est_names, "roc_KL_est")]
    
    
    graphs[['g_120']]  <- result_step_12_prep(step_12, step_12_names)
    graphs[['g_120c']] <- result_step_12_prep(step_12, est_names)
    graphs[['g_120d']] <- result_step_12_prep(step_12, KL_est_names)
    graphs[['g_120e']] <- result_step_12_prep(step_12, KL_GIC_names)
    graphs[['g_120f']] <- result_step_12_prep(step_12, KL_GIC_est_names)
    
    
  }
  
  # tau_c and tau_p values
  if('121' %in% graph_ids){
    
    step_11x_names <- names(step_11x[[1]])
    
    suffix_names <- gsub("^tau_c_", "", step_11x_names)
    est_suffix_names <- grep("_est_", suffix_names, value = TRUE)
    
    
    graphs[['g_121']]  <- result_121_prep(step_11x, step_11y, est_suffix_names)
  }
  
  
  return(graphs)
  
}



visualize_adj_grid <- function(sparse_data_list, all_weeks, absent_week_list, output) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: arrange edge set adjacency matrices in a 2d grid, 
  #       helper function for visualize_discrete_comparison
  #
  # 
  # inputs:
  #
  # - sparse_data_list   (list of lists)      each item is a list of edge coordinates without repeating (j, i) since we have (i, j)
  # - all_weeks          (vector)             all weeks in vector form
  # - absent_week_list   (list of vectors)    for each setting, which weeks are absent so we can gray them out 
  # - output             (string)             'adj', or 'P_HS', or 'C_HS'
  #
  #
  # ----------------------------------------------------------------------------
  
  plot_data_list <- list()
  row_names <- names(sparse_data_list)
  
  new_palette <- hcl.colors(3, palette = 'Blue-Red 2')
  c_low <- new_palette[1]
  c_mid <- new_palette[2]
  c_high <- new_palette[3]  
  
  # 1. Process the Edge Data
  for (r_name in row_names) {
    row_content <- sparse_data_list[[r_name]]
    
    for (c_idx in seq_along(row_content)) {
      df_coords <- row_content[[c_idx]]
      if (is.null(df_coords) || nrow(df_coords) == 0) next
      
      # Determine if we are handling weighted data
      is_weighted <- output %in% c("P_HS", "C_HS")
      
      # Mirroring logic
      # For 'adj', we just have coords. For weighted, we have Node_Row, Node_Col, and Value.
      original <- df_coords
      
      # Create mirrored copy: Swap Row and Col
      # For weighted data, this preserves the 'Value' column correctly
      mirrored <- original
      mirrored$Node_Row <- original$Node_Col
      mirrored$Node_Col <- original$Node_Row
      
      # Combine and remove duplicates (especially important for the diagonal entries)
      combined_df <- unique(rbind(original, mirrored))
      
      combined_df$Row_ID <- r_name
      combined_df$Col_ID <- as.numeric(names(row_content)[c_idx])
      
      plot_data_list[[length(plot_data_list) + 1]] <- combined_df
    }
  }
  
  plot_data <- do.call(rbind, plot_data_list)
  
  # 2. Create the Gray-out Data
  bg_gray_list <- list()
  for (r_name in names(absent_week_list)) {
    absent_weeks <- absent_week_list[[r_name]]
    if (length(absent_weeks) > 0) {
      bg_gray_list[[r_name]] <- data.frame(
        Row_ID = r_name,
        Col_ID = absent_weeks
      )
    }
  }
  bg_gray_data <- do.call(rbind, bg_gray_list)
  
  # 3. Build the Plot
  plot_data$Col_ID <- factor(plot_data$Col_ID, levels = all_weeks)
  if (!is.null(bg_gray_data)) {
    bg_gray_data$Col_ID <- factor(bg_gray_data$Col_ID, levels = all_weeks)
  }
  
  # Initialize the ggplot object
  g <- ggplot() + 
    # Layer 1: Gray out absent weeks
    geom_rect(data = bg_gray_data, 
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              fill = "gray80", alpha = 0.5)
  
  # Layer 2: Plot the actual edges (Binary vs Weighted)
  if (output == "adj") {
    g <- g + geom_tile(data = plot_data, aes(x = Node_Col, y = -Node_Row), fill = "red")
  } else {
    # Weighted plots: Use the Value column for fill
    # viridis or scale_fill_gradient2 (useful if values are log-probs or centered at 0)
    

    zmin <- min(plot_data$Value)
    zmax <- max(plot_data$Value)
    
    zmax <- zmax + 0.1 * (zmax - zmin)
    zmin <- zmin - 0.1 * (zmax - zmin)
    
    
    g <- g + geom_tile(data = plot_data, aes(x = Node_Col, y = -Node_Row, fill = Value)) +
      # scale_fill_viridis_c(option = "magma") + # Good for log densities
      scale_fill_gradient2(low = c_low, mid = c_mid, high = c_high,
                           midpoint = 0,
                           limits = c(zmin, zmax)) + 
      theme(legend.position = "right")      # We actually want the legend here
  }
  
  # 4. Final Formatting
  g <- g + 
    facet_grid(Row_ID ~ Col_ID, drop = FALSE) + 
    coord_fixed() +
    theme_minimal(base_size = 15) + 
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, color = "gray90"), 
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold", size = rel(2))
    )
  
  if(output == "adj") g <- g + theme(legend.position = "none")
  
  return(g)

}

visualize_discrete_comparison <- function(results_folder, ID, time_scale, discrete_levels, output) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: plot adjacency matrix results in a grid-like manner where rows represent different discrete levels
  #       and columns represent weeks. Weeks must be integers. 
  #
  # 
  # inputs:
  #
  # - results_folder   (string)
  # - ID               (string)  'Tau3'
  # - time_scale       (integer)   10 
  # - discrete_levels  (vector of strings)  'm0vr0', 'm1vr1' etc
  # - output           (string)  what to look at. For example 'adj', 'P_HS', 'C_HS'
  #
  #
  #
  # ----------------------------------------------------------------------------
  
  # Prepare to store data for the grid
  
  sparse_data_list <- list()
  absent_week_list <- list()
  
  # 1. Loop through each discrete level (Rows of the grid)
  for (d_level in discrete_levels) {
    
    sparse_data_list_2 <- list()
    
    # Construct file path and load
    file_path <- paste0(results_folder, '/', ID, '_', d_level, '_t', time_scale, '.RData')
    
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    # Load into a temporary environment to avoid overwriting loop variables
    tmp_env <- new.env()
    load(file_path, envir = tmp_env)
    res_i <- tmp_env$graph_results_i
    
    if(output == 'adj'){
      
      # Extract the estimated adjacency matrices (assuming they are in step_12b)
      # We filter for 'est_eig1' as requested
      step_data <- res_i$step_12b
      y_c_weeks <- as.numeric(res_i$y_c_query) # Convert "001.000" etc to integers
      
      absent_weeks <- setdiff(17:38, y_c_weeks)
      
      # 2. Loop through each week (Columns)
      for (idx in seq_along(y_c_weeks)) {
        current_week <- y_c_weeks[idx]
        
        # Target the specific adjacency matrix
        adj_mat <- step_data[[idx]][["adj_mat_KL_GIC_local_est_eig3"]]
        
        if (is.null(adj_mat)) next
        
        # 3. SPARSITY STEP: Get coordinates of 1s only
        coords <- which(adj_mat == 1, arr.ind = TRUE)
        
        if (nrow(coords) > 0) {
          # Keep only upper triangular (Column > Row)
          # This removes symmetry and the diagonal
          coords <- coords[coords[, 2] > coords[, 1], , drop = FALSE]
        }
        
        if (nrow(coords) > 0) {
          df_coords <- as.data.frame(coords)
          colnames(df_coords) <- c("Node_Row", "Node_Col")
          
          # Store for later binding
          sparse_data_list_2[[as.character(current_week)]] <- df_coords
        }
      }
    } else if(output == 'P_HS'){
      
      step_data <- res_i$step_11
      y_c_weeks <- as.numeric(res_i$y_c_query) # Convert "001.000" etc to integers
      
      absent_weeks <- setdiff(17:38, y_c_weeks)
      
      # 2. Loop through each week (Columns)
      for (idx in seq_along(y_c_weeks)) {
        current_week <- y_c_weeks[idx]
        
        # Target the specific matrix
        P_HS_log <- log(step_data[[idx]][["w_mat_KL_est_eig3"]])
        
        if (is.null(P_HS_log)) next
        
        # 1. Get matrix dimensions
        n_size <- nrow(P_HS_log)
        
        # 2. Generate all Upper Triangular indices (including diagonal)
        # row() and col() generate matrices of indices; we filter where col >= row
        upper_tri_idx <- which(col(P_HS_log) >= row(P_HS_log), arr.ind = TRUE)
        
        # 3. Create the data frame
        df_coords <- as.data.frame(upper_tri_idx)
        colnames(df_coords) <- c("Node_Row", "Node_Col")
        
        # 4. Extract the actual values
        # Since you mentioned it's real numbers, you likely want to store them!
        df_coords$Value <- P_HS_log[upper_tri_idx]
        
        # Store for later binding
        sparse_data_list_2[[as.character(current_week)]] <- df_coords
      }
      
    } else if(output == 'C_HS'){
      
      step_data <- res_i$step_11b
      y_c_weeks <- as.numeric(res_i$y_c_query) # Convert "001.000" etc to integers
      
      absent_weeks <- setdiff(17:38, y_c_weeks)
      
      # 2. Loop through each week (Columns)
      for (idx in seq_along(y_c_weeks)) {
        current_week <- y_c_weeks[idx]
        
        # Target the specific matrix
        C_HS_log <- log(step_data[[idx]][["C_HS_KL_est_eig3"]])
        
        if (is.null(C_HS_log)) next
        
        # 1. Get matrix dimensions
        n_size <- nrow(C_HS_log)
        
        # 2. Generate all Upper Triangular indices (including diagonal)
        # row() and col() generate matrices of indices; we filter where col >= row
        upper_tri_idx <- which(col(C_HS_log) >= row(C_HS_log), arr.ind = TRUE)
        
        # 3. Create the data frame
        df_coords <- as.data.frame(upper_tri_idx)
        colnames(df_coords) <- c("Node_Row", "Node_Col")
        
        # 4. Extract the actual values
        # Since you mentioned it's real numbers, you likely want to store them!
        df_coords$Value <- C_HS_log[upper_tri_idx]
        
        # Store for later binding
        sparse_data_list_2[[as.character(current_week)]] <- df_coords
      }
      
    } else{
      stop('error in 28C: visualize_discrete_comparison')
    }
    
    
  
    
    # Clean up large environment immediately
    rm(tmp_env)
    rm(res_i)
    rm(step_data)
    
    sparse_data_list[[d_level]] <- sparse_data_list_2
    
    absent_week_list[[d_level]] <- absent_weeks
    
  }
  

  
  # 4. Assembly and returning
  
  return(visualize_adj_grid(sparse_data_list, 17:38, absent_week_list, output))
  
}





