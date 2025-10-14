

# these functions open an entire dataset and search through all query_id 

# helper function

get_ns <- function(folder_name){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: within a folder, there are several datasets that end in a number then .RData
  #       get all possible n's 
  #
  # inputs:
  #
  # - folder_name     (string)    'simu_results/banded_trig2/OG'
  # 
  #
  # outputs:
  #
  # - ns              (vector)
  #
  # ----------------------------------------------------------------------------  
  
  # 1) find the file in the folder and load it 
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <- suppressWarnings({
    as.numeric(sub(".*_(.*)\\.RData$", "\\1", files))   # retrieve number before .RData and after recent underscore
  })
  
  ns <- ns[!is.na(ns)]
  
  return(ns)
}

# helper function
get_file_name_and_load <- function(folder_name, n){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: within a folder, there are several datasets that end in a number then .RData
  #       load that file
  #
  # inputs:
  #
  # - folder_name     (string)    'simu_results/banded_trig2/OG'
  # - n               (integer)   100
  # 
  #
  # outputs:
  #
  # - NULL, just loads that file
  #
  # ----------------------------------------------------------------------------
  
  # 1) find the file in the folder and load it 
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <- suppressWarnings({
    as.numeric(sub(".*_(.*)\\.RData$", "\\1", files))   # retrieve number before .RData and after recent underscore
  })
  
  
  if(n %in% ns){
    idx <- which(ns == n)
  } else{
    stop('ERROR on get_file_name: n not available')
  }  
  
  env <- new.env()
  
  load(files[idx], envir = env) 
  
  return(as.list(env)[[1]])
  
}

# helper function
get_y_c_query <- function(folder_name, n){
  
  graph_results_i <- get_file_name_and_load(folder_name, n)
  
  return(names(graph_results_i$estimated_graphs_part_2))
  
}

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
  
  return(grid.arrange(grobs = c(plots_truth, plots_w_est), nrow = 2, ncol = n_graphs,
                      left = textGrob("Est \t Ground Truth", rot = 90, gp = gpar(fontsize = 16))))
  
  
  
}
