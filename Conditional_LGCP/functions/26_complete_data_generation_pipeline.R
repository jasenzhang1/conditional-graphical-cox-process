simulate_conditional_cox_data_v4 <- function(
  n,                     # Sample size (n)
  p,                     # Number of processes (p)  
  T_max,                 # Time horizon (T)
  query_y_cs,            # matrix of query y_cs (num_query x q_c)
  adj_type,              # pxp precion matrix structure
  adj_params,            # associated parameters
  time_grid,             # Time discretization (m-dim vec)
  time_grid_est,         # Time discretization of the estimate (m_est-dim vec)
  base_kernel_params,
  ncores,
  seed = NULL,
  verbose = FALSE,        # do we return everything?
  parallel = TRUE
){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: simulate data
  #
  # - same as v3 but parallelized
  # 
  # - 9/10/2025 - cutting down lots of parameters
  #
  # ----------------------------------------------------------------------------
  
  if(is.null(seed)){
    seed_mc <- TRUE  # if seed is null, let mclapply take over
  } else{
    seed_mc <- FALSE # if we specify a seed, don't let mclapply take over
  }
  
  
  # time_grid
  
  
  time_grid_both <- sort(union(time_grid, time_grid_est))
  
  m <- length(time_grid)
  m_est <- length(time_grid_est)
  m2 <- length(time_grid_both)
  
  
  # 1) generate y_c

  Y_continuous <- generate_y_c_adj_type(n, adj_type, adj_params)
  
  
  # 2) for each subject, generate their parameters and event data
  
  print('at subject generation')
  
  subject_data <- simulate_subject_data(n, p, Y_continuous, adj_type, adj_params,
                                        time_grid, time_grid_est,
                                        base_kernel_params, ncores, parallel)

  # 3) for each query point, generate parameters
  print('at query point generation')
  
  apply_fun <- function(X, FUN, ...) {
    if (parallel) {
      mclapply(X, FUN, mc.cores = ncores, ...)
    } else {
      lapply(X, FUN, ...)
    }
  }
  query_data <- apply_fun(
    1:nrow(query_y_cs), 
    function(k){
      y_c_k <- query_y_cs[k, ]  
      
      # ground truth precision matrix
      prec_mat_truth <- generate_sparse_precision_matrix(y_c_k, p, adj_type, adj_params)
      
      # Generate precision operator P^{(y_c^k, y_d^k)} and adjacency matrix
      mats_k <- sample_conditional_precision_v3(
        time_grid, time_grid_est, 
        base_kernel_params,
        prec_mat_truth,
        y_c_k
      )
      
      mats_k
    }
  )
  
  # 5) get list of event_times for each subject (n) and process (p)
  
  event_times_list <- list()
  
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      key <- paste(i, j, sep = "_")
      event_times_list[[key]] <- subject_data[[i]]$event_times[[j]]
    }
  }
  
  # 5.1) get every single X_k vector (n x p x m)
  
  X_k_truth <- array(
    unlist(lapply(subject_data, function(x) x$X_functions)),
    dim = c(p, m, length(subject_data))
  )
  
  X_k_coarse_truth <- array(
    unlist(lapply(subject_data, function(x) x$X_functions_coarse)),
    dim = c(p, m_est, length(subject_data))
  )  
  
  X_k_both_truth <- array(
    unlist(lapply(subject_data, function(x) x$X_functions_full)),
    dim = c(p, m2, length(subject_data))
  )    
  
  # 7) Compute summary statistics
  
  total_events <- sum(sapply(event_times_list, length))
  avg_events_per_process <- total_events / (n * p)
  
  cat("Data generation completed.\n")
  cat("  Total events:", total_events, "\n")
  cat("  Total replicates:", n, "\n")
  cat("  Total processes:", p, "\n")
  cat("  Average events per replicate per process:", round(avg_events_per_process, 2), "\n")
  
  
  result <- list(
    
    # event times and subject data
    event_times = event_times_list,      # List: key = "k_i" -> event times
    X_k_truth = X_k_truth,               # all log intensities used to generate data (p x m x n)
    X_k_coarse_truth = X_k_coarse_truth, # all log intensities used to generate data coarsely (p x m_est x n)
    X_k_both_truth = X_k_both_truth,
    Y_continuous = Y_continuous,         # n x q_c matrix
    
    # Simulation parameters
    simulation_params = list(
      n = n, p = p, T_max = T_max, query_y_cs = query_y_cs,
      adj_type = adj_type,
      adj_params = adj_params,     
      time_grid = time_grid,               
      time_grid_est = time_grid_est,
      time_grid_both = time_grid_both,      
      seed = seed
    ),
    
    # ground truths - adj_mat, prec_mat for all queries
    true_graphs = query_data,             
    
    summary_stats = list(
      total_events = total_events,
      avg_events_per_process = avg_events_per_process
    )
  )
  
  if(!verbose){
    result$subject_data = subject_data         # Complete subject-level data
  }
  
  return(result)  
}

simulate_subject_data <- function(n, p, Y_continuous, adj_type, adj_params, 
                                  time_grid, time_grid_est,
                                  base_kernel_params, ncores, parallel){
  
  # Wrapper function to unify interface
  apply_fun <- function(X, FUN, ...) {
    if (parallel) {
      pbmclapply(X, FUN, mc.cores = ncores, ...)
    } else {
      lapply(X, FUN, ...)
    }
  }  
  
  subject_data <- apply_fun(1:n, function(k){
    
    if(k %% 50 == 0){
      print(paste0(k, ' out of ', n))
    }
    
    
    
    
    y_c_k <- Y_continuous[k, ]  
    
    
    # 4.2) ground truth precision matrix
    prec_mat_truth <- generate_sparse_precision_matrix(y_c_k, p, adj_type, adj_params)
    
    
    
    # 4.3) Generate precision operator P^{(y_c^k, y_d^k)} and adjacency matrix E_{y_c^k, y_d^k}
    mats_k <- sample_conditional_precision_v3(time_grid, time_grid_est,
                                              base_kernel_params,
                                              prec_mat_truth,
                                              y_c_k)
    
    
    
    # generate log-intensity functions
    
    
    
    X_k_list <- generate_log_intensity_functions_full_mat(
      kronecker(prec_mat_truth$cor_mat, mats_k$P_block_kronecker$base_cov_both),
      time_grid,
      time_grid_est,
      baseline_mean = mats_k$P_block_kronecker$base_mean_both,
      sample_mode = 'multi')
    
    # [[1]] = full
    # [[2]] = simulation step size (50)
    # [[3]] = estimation step size (19)
    
    
    
    # Generate point process events
    events_k <- generate_cox_process_events(
      X_k_list[[2]],
      time_grid,
      T_max,
      max_intensity = Inf
    )
    
    
    
    # Store complete subject information
    
    list(
      Y_continuous = y_c_k,
      X_functions_full = X_k_list[[1]],
      X_functions = X_k_list[[2]],
      X_functions_coarse = X_k_list[[3]],
      precision_and_graph = mats_k,
      event_times = events_k$event_times,
      event_counts = events_k$event_counts
    )
    
    
  }) # end of apply_fun
}



extract_event_times_df <- function(subject_list) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: helper function for convert_data_for_estimation
  # 
  # ----------------------------------------------------------------------------
  
  
  do.call(rbind, lapply(seq_along(subject_list), function(subject_id) {
    event_times <- subject_list[[subject_id]]$event_times
    
    # Handle if event_times is NULL or missing
    if (is.null(event_times)) return(NULL)
    
    do.call(rbind, lapply(seq_along(event_times), function(event_id) {
      values <- event_times[[event_id]]
      
      if (length(values) == 0) return(NULL)  # skip empty vectors
      
      data.frame(
        value = values,
        event_id = event_id,
        subject_id = subject_id
      )
    }))
  }))
}
