# source('21_graph_structure_generation.R')
# source('22_generate_random_variables.R')
# source('23_conditional_dependence_functions.R')
# source('24_log_intensity_generation.R')
# source('25_point_generation_process.R')


simulate_conditional_cox_data <- function(
  n = 100,                    # Sample size (n)
  p = 10,                     # Number of processes (p)  
  T_max = 10,                 # Time horizon (T)
  q_c = 2,                    # Continuous conditioning dimension (q_c)
  K = 4,                      # Discrete combinations (K)
  sparsity = 0.2,             # Graph sparsity (s)
  theta = 1.0,                # Signal strength (theta)
  graph_type = "random",      # Graph topology
  dependence_type = "linear", # Conditional dependence type
  time_grid_size = 50,        # Time discretization (m)
  seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  cat("Generating simulation data with parameters:\n")
  cat("  n =", n, ", p =", p, ", T =", T_max, "\n")
  cat("  q_c =", q_c, ", K =", K, ", sparsity =", sparsity, "\n")
  cat("  theta =", theta, ", dependence =", dependence_type, "\n")
  
  # Time grid: {t_1, ..., t_m}
  time_grid <- seq(0, T_max, length.out = time_grid_size)
  
  # Step 1: Code 22, Generate conditioning variables {Y_c^k, Y_d^k}_{k=1}^n
  conditioning_vars <- generate_conditioning_variables(n, q_c, K, seed)
  Y_continuous <- conditioning_vars$Y_continuous  # n x q_c
  Y_discrete <- conditioning_vars$Y_discrete      # n x q_d
  Y_discrete_keys <- paste0(Y_discrete[,1], '_', Y_discrete[,2])

  strata_keys <- unique(Y_discrete_keys)
  
  
  
  adj_mats <- list()
  n_edges <- list()
  precision_spec <- list()
  for(i in 1:length(strata_keys)){
    
    strata_key <- strata_keys[i]
  
    # Step 2: For each strata, generate graph structure E_{y_d}
    if (graph_type == "random") {
      adj_matrix <- generate_random_graph(p, sparsity, seed - i)
    } else if (graph_type == "scale_free") {
      adj_matrix <- generate_scale_free_graph(p, sparsity, seed - i)
    } else if (graph_type == "block") {
      adj_matrix <- generate_block_diagonal_graph(p, sparsity, seed = seed - i)
    } else {
      stop("Unknown graph_type")
    }
    
    adj_mats[[strata_key]] <- adj_matrix
    
    n_edges[[strata_key]] <- sum(adj_matrix) / 2
    cat("  Generated", graph_type, "graph with", sum(adj_matrix) / 2, "edges\n")
    
    # Step 3: Code 23, Generate precision operator specification = generate betas
    precision_spec[[strata_key]] <- generate_precision_operators(
      adj_mats[[strata_key]], theta, q_c, dependence_type, seed
    )    
    
  }  # exit discrete strata loop  
  

  
  # Step 4: Generate data for each subject k = 1, ..., n
  event_times_list <- list()  # For main algorithm format
  subject_data <- list()      # For analysis
  
  for (k in 1:n) {
    if (k %% 50 == 0) cat("  Generating subject", k, "/", n, "\n")
    
    # Extract conditioning values for subject k
    y_c_k <- Y_continuous[k, ]  # Y_c^k
    y_d_k <- Y_discrete[k, ]    # Y_d^k
    ydk_key <- paste0(y_d_k[1], '_', y_d_k[2])
    
    # Generate precision operator P^{(y_c^k, y_d^k)}
    # P_k = (p x p x m x m)
    P_k <- sample_conditional_precision(precision_spec[[ydk_key]], y_c_k, time_grid) # code 24
    
    # Generate log-intensity functions X^k = {X_i^k(t)}_{i=1}^p
    # X_k = (p x m matrix)
    X_k <- generate_log_intensity_functions(P_k, time_grid, seed = seed + k)
    
    # Generate point process events {N_i^k}_{i=1}^p
    # events_k = (p x m matrix)
    events_k <- generate_cox_process_events(X_k, time_grid, T_max, seed = seed + k + n)
    
    # Store in event_times format expected by main algorithm
    # key = subject_process
    for (i in 1:p) {
      key <- paste(k, i, sep = "_")
      event_times_list[[key]] <- events_k$event_times[[i]]
    }
    
    # Store complete subject information
    subject_data[[k]] <- list(
      Y_continuous = y_c_k,
      Y_discrete = y_d_k,
      X_functions = X_k,
      precision_operator = P_k,
      event_times = events_k$event_times,
      event_counts = events_k$event_counts
    )
  }
  
  # Compute summary statistics
  total_events <- sum(sapply(event_times_list, length))
  avg_events_per_process <- total_events / (n * p)
  
  cat("Data generation completed.\n")
  cat("  Total events:", total_events, "\n")
  cat("  Average events per process:", round(avg_events_per_process, 2), "\n")
  
  return(list(
    # Data in format expected by main algorithm
    event_times = event_times_list,      # List: key = "k_i" -> event times
    Y_continuous = Y_continuous,         # n x q_c matrix
    Y_discrete = Y_discrete,             # n x q_d matrix  
    time_grid = time_grid,               # m x 1 vector
    
    # Simulation parameters and ground truth
    simulation_params = list(
      n = n, p = p, T_max = T_max, q_c = q_c, K = K,
      sparsity = sparsity, theta = theta,
      graph_type = graph_type, dependence_type = dependence_type,
      time_grid_size = time_grid_size, seed = seed
    ),
    true_graph = adj_mats,             # p x p adjacency matrix for each strata
    precision_specification = precision_spec, # one for each strata
    subject_data = subject_data,         # Complete subject-level data
    summary_stats = list(
      total_events = total_events,
      avg_events_per_process = avg_events_per_process,
      n_edges = n_edges
    )
  ))
}

simulate_conditional_cox_data_og <- function(
  n = 100,                    # Sample size (n)
  p = 10,                     # Number of processes (p)  
  T_max = 10,                 # Time horizon (T)
  q_c = 2,                    # Continuous conditioning dimension (q_c)
  K = 4,                      # Discrete combinations (K)
  sparsity = 0.2,             # Graph sparsity (s)
  theta = 1.0,                # Signal strength (theta)
  graph_type = "random",      # Graph topology
  dependence_type = "linear", # Conditional dependence type
  time_grid_size = 50,        # Time discretization (m)
  seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  cat("Generating simulation data with parameters:\n")
  cat("  n =", n, ", p =", p, ", T =", T_max, "\n")
  cat("  q_c =", q_c, ", K =", K, ", sparsity =", sparsity, "\n")
  cat("  theta =", theta, ", dependence =", dependence_type, "\n")
  
  # Time grid: {t_1, ..., t_m}
  time_grid <- seq(0, T_max, length.out = time_grid_size)
  
  # Step 1: Generate graph structure E_{y_d}
  if (graph_type == "random") {
    adj_matrix <- generate_random_graph(p, sparsity, seed)
  } else if (graph_type == "scale_free") {
    adj_matrix <- generate_scale_free_graph(p, sparsity, seed)
  } else if (graph_type == "block") {
    adj_matrix <- generate_block_diagonal_graph(p, sparsity, seed = seed)
  } else {
    stop("Unknown graph_type")
  }
  
  n_edges <- sum(adj_matrix) / 2
  cat("  Generated", graph_type, "graph with", n_edges, "edges\n")
  
  # Step 2: Code 22, Generate conditioning variables {Y_c^k, Y_d^k}_{k=1}^n
  conditioning_vars <- generate_conditioning_variables(n, q_c, K, seed)
  Y_continuous <- conditioning_vars$Y_continuous  # n x q_c
  Y_discrete <- conditioning_vars$Y_discrete      # n x q_d
  
  # Step 3: Code 23, Generate precision operator specification = generate betas
  precision_spec <- generate_precision_operators(
    adj_matrix, theta, q_c, dependence_type, seed
  )
  
  # Step 4: Generate data for each subject k = 1, ..., n
  event_times_list <- list()  # For main algorithm format
  subject_data <- list()      # For analysis
  
  for (k in 1:n) {
    if (k %% 50 == 0) cat("  Generating subject", k, "/", n, "\n")
    
    # Extract conditioning values for subject k
    y_c_k <- Y_continuous[k, ]  # Y_c^k
    y_d_k <- Y_discrete[k, ]    # Y_d^k
    
    # Generate precision operator P^{(y_c^k, y_d^k)}
    # P_k = (p x p x m x m)
    P_k <- sample_conditional_precision(precision_spec, y_c_k, time_grid) # code 24
    
    # Generate log-intensity functions X^k = {X_i^k(t)}_{i=1}^p
    X_k <- generate_log_intensity_functions(P_k, time_grid, seed = seed + k)
    
    # Generate point process events {N_i^k}_{i=1}^p
    # events_k = (p x m matrix)
    events_k <- generate_cox_process_events(X_k, time_grid, T_max, seed = seed + k + n)
    
    # Store in event_times format expected by main algorithm
    # key = subject_process
    for (i in 1:p) {
      key <- paste(k, i, sep = "_")
      event_times_list[[key]] <- events_k$event_times[[i]]
    }
    
    # Store complete subject information
    subject_data[[k]] <- list(
      Y_continuous = y_c_k,
      Y_discrete = y_d_k,
      X_functions = X_k,
      precision_operator = P_k,
      event_times = events_k$event_times,
      event_counts = events_k$event_counts
    )
  }
  
  # Compute summary statistics
  total_events <- sum(sapply(event_times_list, length))
  avg_events_per_process <- total_events / (n * p)
  
  cat("Data generation completed.\n")
  cat("  Total events:", total_events, "\n")
  cat("  Average events per process:", round(avg_events_per_process, 2), "\n")
  
  return(list(
    # Data in format expected by main algorithm
    event_times = event_times_list,      # List: key = "k_i" -> event times
    Y_continuous = Y_continuous,         # n x q_c matrix
    Y_discrete = Y_discrete,             # n x q_d matrix  
    time_grid = time_grid,               # m x 1 vector
    
    # Simulation parameters and ground truth
    simulation_params = list(
      n = n, p = p, T_max = T_max, q_c = q_c, K = K,
      sparsity = sparsity, theta = theta,
      graph_type = graph_type, dependence_type = dependence_type,
      time_grid_size = time_grid_size, seed = seed
    ),
    true_graph = adj_matrix,             # p x p adjacency matrix
    precision_specification = precision_spec,
    subject_data = subject_data,         # Complete subject-level data
    summary_stats = list(
      total_events = total_events,
      avg_events_per_process = avg_events_per_process,
      n_edges = n_edges
    )
  ))
}

simulate_conditional_cox_data_v2 <- function(
  n = 100,                    # Sample size (n)
  p = 10,                     # Number of processes (p)  
  T_max = 10,                 # Time horizon (T)
  q_c = 1,                    # Continuous conditioning dimension (q_c)
  y_c_borders = list(1:9),    # Border values
  threshold = 0.4,            # Threshold
  theta = 1.0,                # Signal strength (theta)
  dependence_type = "constant", # Conditional dependence type
  time_grid_size = 50,        # Time discretization (m)
  seed = NULL  
){

  # generate all ground truth graphs
  
  graphs <- generate_precision_operators_and_matrix(p, theta, q_c, 
                                                    y_c_borders,
                                                    threshold,
                                                    dependence_type = "constant", seed = seed)
  
  simu_settings <- collect_beta_and_parameters(p, theta, q_c, 
                                               y_c_borders,
                                               threshold,
                                               dependence_type = "constant", seed = NULL)
  
  # 2) generate n samples whose continuous variables vary. Each of them belong in a bucket and are assigned a graph
  
  time_grid <- seq(0, T_max, length.out = time_grid_size)
  
  # Code 22, Generate conditioning variables {Y_c^k}_{k=1}^n
  week_start <- y_c_borders[[1]][1]
  week_end <- max(y_c_borders[[1]]) + 1
  Y_list <- generate_conditioning_variables_one_strata_weeks(n, week_start, week_end, seed)
  Y_continuous <- Y_list$Y_continuous
  
  if(q_c == 1){
    Y_continuous <- matrix(Y_continuous)
  }
  
  # Step 4: Generate data for each subject k = 1, ..., n
  event_times_list <- list()  # For main algorithm format
  subject_data <- list()      # For analysis
  
  for (k in 1:n) {
    if (k %% 50 == 0) cat("  Generating subject", k, "/", n, "\n")
    
    # Extract conditioning values for subject k
    y_c_k <- Y_continuous[k, ]  # Y_c^k
    
    # find its region and its adj_mat
    
    region <- find_yc_group(y_c_k, y_c_borders)   # less than or equal
    region_id <- paste(region, collapse = "_")
    graph_k <- find_E_yc_yd_graph(region, graphs$adjacency)
    
    
    
    # Generate precision operator P^{(y_c^k, y_d^k)}
    # P_k = (p x p x m x m)
    P_k <- sample_conditional_precision(graphs, y_c_k, time_grid, region_id) # code 24
    
    # Generate log-intensity functions X^k = {X_i^k(t)}_{i=1}^p
    # X_k = (p x m matrix)
    X_k <- generate_log_intensity_functions(P_k, time_grid, seed = seed + k)
    
    # Generate point process events {N_i^k}_{i=1}^p
    # events_k = (p x m matrix)
    events_k <- generate_cox_process_events(X_k, time_grid, T_max, seed = seed + k + n)
    
    # Store in event_times format expected by main algorithm
    # key = subject_process
    for (i in 1:p) {
      key <- paste(k, i, sep = "_")
      event_times_list[[key]] <- events_k$event_times[[i]]
    }
    
    # Store complete subject information
    subject_data[[k]] <- list(
      Y_continuous = y_c_k,
      # Y_discrete = y_d_k,
      X_functions = X_k,
      precision_operator = P_k,
      event_times = events_k$event_times,
      event_counts = events_k$event_counts
    )
  }
  
  # Compute summary statistics
  total_events <- sum(sapply(event_times_list, length))
  avg_events_per_process <- total_events / (n * p)
  
  cat("Data generation completed.\n")
  cat("  Total events:", total_events, "\n")
  cat("  Average events per process:", round(avg_events_per_process, 2), "\n")
  
  return(list(
    # Data in format expected by main algorithm
    event_times = event_times_list,      # List: key = "k_i" -> event times
    Y_continuous = Y_continuous,         # n x q_c matrix
    time_grid = time_grid,               # m x 1 vector
    
    # Simulation parameters and ground truth
    simulation_params = list(
      n = n, p = p, T_max = T_max, q_c = q_c, y_c_borders = y_c_borders,
      threshold = threshold, theta = theta,
      dependence_type = dependence_type,
      time_grid_size = time_grid_size, seed = seed
    ),
    true_graphs = graphs$adjacency,             # p x p adjacency matrix for each strata
    precision_specification = graphs, # one for each strata
    subject_data = subject_data,         # Complete subject-level data
    summary_stats = list(
      total_events = total_events,
      avg_events_per_process = avg_events_per_process
      # n_edges = 
    )
  ))  
}

simulate_conditional_cox_data_v3 <- function(
  n = 100,                    # Sample size (n)
  p = 10,                     # Number of processes (p)  
  T_max = 10,                 # Time horizon (T)
  q_c = 1,                    # Continuous conditioning dimension (q_c)
  y_c_borders = list(1:9),    # Border values
  sparsity = 0.2,             # Sparsity
  theta = 1.0,                # Signal strength (theta)
  dependence_type = "constant", # Conditional dependence type
  time_grid,                    # Time discretization (m)
  base_kernel_params,
  seed = NULL
){
  
  # 1) collect beta coefficients
  # 2) generate h_ij, and create adj_mat + prec_mat together with value of tau (threshold)
  
  
  m <- length(time_grid)
  
  true_graphs <- list()

  
  simu_settings <- collect_beta_and_parameters(p, time_grid, theta, q_c, 
                                               y_c_borders,
                                               sparsity,
                                               base_kernel_params,
                                               dependence_type = "constant", seed = seed)
  
  # 2) generate n samples whose continuous variables vary. Each of them belong in a bucket and are assigned a graph
  
  
  
  # Code 22, Generate conditioning variables {Y_c^k}_{k=1}^n
  week_start <- y_c_borders[[1]][1]
  week_end <- max(y_c_borders[[1]]) + 1
  Y_list <- generate_conditioning_variables_one_strata_weeks(n, week_start, week_end, seed)
  Y_continuous <- Y_list$Y_continuous
  
  if(q_c == 1){
    Y_continuous <- matrix(Y_continuous)
  }
  
  # Step 4: Generate data for each subject k = 1, ..., n
  event_times_list <- list()  # For main algorithm format
  subject_data <- list()      # For analysis
  
  for (k in 1:n) {
    if (k %% 10 == 0) cat("  Generating subject", k, "/", n, "\n")
    
    # 4.1) Extract conditioning values for subject k
    y_c_k <- Y_continuous[k, ]  # Y_c^k   
  
    region <- find_yc_group(y_c_k, y_c_borders)   # less than or equal
    region_id <- paste(region, collapse = "_")    
    
    # 4.2) ground truth precision matrix
    prec_mat_truth <- generate_sparse_precision_matrix(p = p, graph_type = 'banded')
    
    
    # 4.3) Generate precision operator P^{(y_c^k, y_d^k)} and adjacency matrix E_{y_c^k, y_d^k}
    # P_k = (p x p x m x m)
    mats_k <- sample_conditional_precision_v3(simu_settings,  
                                              prec_mat_truth,
                                              y_c_k) # code 24    
    
    
    if(! region_id %in% names(true_graphs)){
      true_graphs[[region_id]] <- mats_k
    }
    
    
    # Generate log-intensity functions X^k = {X_i^k(t)}_{i=1}^p
    # X_k = (p x m matrix)
    
    # X_k_simple <- generate_log_intensity_functions_full_mat(GP_cov, time_grid, baseline_mean = 5, seed = seed + k)
    X_k <-        generate_log_intensity_functions_full_mat(mats_k$P_block_kronecker$GP_simu_var, time_grid, 
                                                            baseline_mean = mats_k$P_block_kronecker$GP_simu_mean, 
                                                            sample_mode = 'multi', seed = seed + k)
    
    # visualize_log_intensity(X_k_simple, time_grid)
    # visualize_log_intensity(X_k, time_grid)

    
    # Generate point process events {N_i^k}_{i=1}^p
    # events_k = (p x m matrix)
    # max intensity may matter??
    events_k <- generate_cox_process_events(X_k, time_grid, T_max, max_intensity = Inf, seed = seed + k + n)
    
    # Store in event_times format expected by main algorithm
    # key = subject_process
    for (i in 1:p) {
      key <- paste(k, i, sep = "_")
      event_times_list[[key]] <- events_k$event_times[[i]]
    }
    
    # Store complete subject information
    subject_data[[k]] <- list(
      Y_continuous = y_c_k,
      X_functions = X_k,
      precision_and_graph = mats_k,
      event_times = events_k$event_times,
      event_counts = events_k$event_counts
    )
  }
  
  # Compute summary statistics
  
  threshold_p <- mats_k$threshold_p
  
  total_events <- sum(sapply(event_times_list, length))
  avg_events_per_process <- total_events / (n * p)
  
  cat("Data generation completed.\n")
  cat("  Total events:", total_events, "\n")
  cat("  Average events per process:", round(avg_events_per_process, 2), "\n")
  
  return(list(
    
    # event times and subject data
    event_times = event_times_list,      # List: key = "k_i" -> event times
    subject_data = subject_data,         # Complete subject-level data
    Y_continuous = Y_continuous,         # n x q_c matrix
    
    
    time_grid = time_grid,               # m x 1 vector
    
    # Simulation parameters
    simulation_params = list(
      n = n, p = p, T_max = T_max, q_c = q_c, y_c_borders = y_c_borders,
      threshold_p = threshold_p,
      theta = theta,
      dependence_type = dependence_type,
      time_grid_size = time_grid_size, 
      seed = seed
    ),
    
    # ground truths - adj_mat, prec_mat for all y_c levels
    true_graphs = true_graphs,             
    
    summary_stats = list(
      total_events = total_events,
      avg_events_per_process = avg_events_per_process
    )
  ))  
}

simulate_conditional_cox_data_v4 <- function(
  n = 100,                    # Sample size (n)
  p = 10,                     # Number of processes (p)  
  T_max = 10,                 # Time horizon (T)
  q_c = 1,                    # Continuous conditioning dimension (q_c)
  y_c_borders = list(1:9),    # Border values
  sparsity = 0.2,             # Sparsity
  theta = 1.0,                # Signal strength (theta)
  dependence_type = "constant", # Conditional dependence type
  adj_type = 'banded_c1',       # pxp precion matrix structure
  time_grid,                    # Time discretization (m)
  time_grid_est, 
  base_kernel_params,
  seed = NULL,
  ncores
){
  
  # same as v3 but parallelized
  
  m <- length(time_grid)
  m_est <- length(time_grid_est)
  time_grid_both <- sort(union(time_grid, time_grid_est))
  m2 <- length(time_grid_both)
  
  # 1) collect beta coefficients
  
  simu_settings <- collect_beta_and_parameters(p, time_grid, time_grid_est, theta, q_c, 
                                               y_c_borders,
                                               sparsity,
                                               base_kernel_params,
                                               dependence_type, 
                                               adj_type,
                                               seed = seed)
  
  
  # 2) generate n samples whose continuous variables vary. Each of them belong in a bucket and are assigned a graph
  
  
  
  # Code 22, Generate conditioning variables {Y_c^k}_{k=1}^n
  week_start <- y_c_borders[[1]][1]
  week_end <- max(y_c_borders[[1]]) + 1
  Y_list <- generate_conditioning_variables_one_strata_weeks(n, week_start, week_end, seed)
  Y_continuous <- Y_list$Y_continuous
  
  if(q_c == 1){
    Y_continuous <- matrix(Y_continuous)
  }
  
  # Step 4: Generate data for each subject k = 1, ..., n
  subject_data <- pbmclapply(1:n, function(k){
    
    
    # 4.1) Extract conditioning values for subject k
    y_c_k <- Y_continuous[k, ]  # Y_c^k   
    
    region <- find_yc_group(y_c_k, y_c_borders)   # less than or equal
    region_id <- paste(region, collapse = "_")    
    
    # 4.2) ground truth precision matrix
    prec_mat_truth <- generate_sparse_precision_matrix(p, simu_settings$adj_type, y_c_k)
    
    
    # 4.3) Generate precision operator P^{(y_c^k, y_d^k)} and adjacency matrix E_{y_c^k, y_d^k}
    # P_k = (p x p x m x m)
    mats_k <- sample_conditional_precision_v3(simu_settings,  
                                              prec_mat_truth,
                                              y_c_k) # code 24    
    
  
    
    # generate log-intensity functions
    X_k_list <- generate_log_intensity_functions_full_mat(
      mats_k$P_block_kronecker$GP_simu_var_both, 
      time_grid, 
      time_grid_est,
      baseline_mean = mats_k$P_block_kronecker$GP_simu_mean_both,
      sample_mode = 'multi', 
      seed = seed + k)
    
    X_k_full <- X_k_list[[1]]
    X_k <- X_k_list[[2]]
    

    # Generate point process events
    events_k <- generate_cox_process_events(
      X_k, 
      time_grid, 
      T_max, 
      max_intensity = Inf, 
      seed = seed + k + n
    )
    

    
    # Store complete subject information
    
    list(
      region_id = region_id,
      Y_continuous = y_c_k,
      X_functions_full = X_k_full,
      X_functions = X_k,
      X_functions_coarse = X_k_list[[3]],
      precision_and_graph = mats_k,
      event_times = events_k$event_times,
      event_counts = events_k$event_counts
    )
  }, mc.cores = ncores) # end of pbmclapply
  
  
  
  # 5) extract p*n length list of event_times for each subject and process
  
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
  
  # 6) get true graphs
  
  true_graphs <- list()
  for (i in seq_len(n)) { 
    if(! subject_data[[i]]$region_id %in% names(true_graphs)){
      true_graphs[[subject_data[[i]]$region_id]] <- subject_data[[i]]$precision_and_graph
    }
  }
  
  # 7) Compute summary statistics
  
  threshold_p <- subject_data[[1]]$precision_and_graph$threshold_p
  
  total_events <- sum(sapply(event_times_list, length))
  avg_events_per_process <- total_events / (n * p)
  
  cat("Data generation completed.\n")
  cat("  Total events:", total_events, "\n")
  cat("  Average events per process:", round(avg_events_per_process, 2), "\n")
  
  return(list(
    
    # event times and subject data
    event_times = event_times_list,      # List: key = "k_i" -> event times
    subject_data = subject_data,         # Complete subject-level data
    X_k_truth = X_k_truth,               # all log intensities used to generate data (p x m x n)
    X_k_coarse_truth = X_k_coarse_truth, # all log intensities used to generate data coarsely (p x m_est x n)
    X_k_both_truth = X_k_both_truth,
    Y_continuous = Y_continuous,         # n x q_c matrix
    time_grid = time_grid,               # m x 1 vector
    time_grid_est = time_grid_est,
    time_grid_both = time_grid_both,
    
    # Simulation parameters
    simulation_params = list(
      n = n, p = p, T_max = T_max, q_c = q_c, y_c_borders = y_c_borders,
      threshold_p = threshold_p,
      theta = theta,
      dependence_type = dependence_type,
      time_grid_size = time_grid_size, 
      seed = seed
    ),
    
    # ground truths - adj_mat, prec_mat for all y_c levels
    true_graphs = true_graphs,             
    
    summary_stats = list(
      total_events = total_events,
      avg_events_per_process = avg_events_per_process
    )
  ))  
}

extract_event_times_df <- function(subject_list) {
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

convert_data_for_estimation <- function(subject_list){
  
  # subject_list (list)
  # - Y_continuous
  # - X_functions
  # - precision_operator
  # - event_times
  # - event_counts
  #
  # 
  # Output:
  #
  # - df (data.frame with 'feature_id', 'time', and 'subject_num')
  #
  #   - feature_id
  #   - time
  #   - subject_num
  
  df <- extract_event_times_df(subject_list)
  colnames(df) <- c('time', 'feature_id', 'subject_num')

  return(df)  
  
}
