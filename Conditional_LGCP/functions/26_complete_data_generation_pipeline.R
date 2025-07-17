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
