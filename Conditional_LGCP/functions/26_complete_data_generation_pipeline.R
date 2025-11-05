source('functions/21b_generate_finite_basis_expansion.R')
source('functions/22_generate_random_variables.R')
source('functions/25_point_generation_process.R')

library(abind)

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
  # GOAL: simulate data according to
  #
  # - a specific adj_type with adj_params
  # - a specific base_kernel with base_kernel_params
  #
  # 
  # inputs:
  #
  # - n               (integer)          sample size
  # - p               (integer)          number of processes
  # - T_max           (number)           point processes lie in [0, T_max], usually 1
  # - query_y_cs      (n x q_c matrix)   continuous covariates
  # - adj_type        (string)           precision matrix adjacency type
  # - adj_params      (vector)           associated vector of parameters of the adjacency type
  # - time_grid       (m-dim vec)        discretized time points from 0 to T_max
  # - time_grid_est   (m_est-dim vec)    discretized time points from 0 to T_max for estimation
  # - ncores          (integer)          number of cores
  # - seed            (integer)          seed number for reproducibility 
  # - verbose         (boolean)          do we want to see extra outputs?
  # - parallel        (boolean)          do we want to enable parallel computing?
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
  
  # package results
  
  result <- package_simulation_results(event_times_list, n, p, T_max, query_y_cs,
                                       adj_type, adj_params, time_grid, time_grid_est, time_grid_both, seed,
                                       X_k_truth, X_k_coarse_truth, X_k_both_truth, Y_continuous)

  
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


package_simulation_results <- function(event_times_list, n, p, T_max, query_y_cs,
                                       adj_type, adj_params, time_grid, time_grid_est, time_grid_both, seed,
                                       X_k_truth, X_k_coarse_truth, X_k_both_truth, Y_continuous,
                                       query_data){
  
  # ----------------------------------------------------------------------------
  #
  # 
  # GOAL: after generating points and log-intensities, package them to be ready to be estimated
  #
  #       - helper function for both simulation functions
  # 
  # inputs:
  #
  # - event_times_list (np-dim list of event times)   'k_i' item names for subject k and process i
  # - n
  # - p
  # - T_max
  # - query_y_cs            (n_query x q_c matrix)
  # - adj_type
  # - adj_params
  # - time_grid
  # - time_grid_est
  # - time_grid_both
  # - seed
  # 
  # - X_k_truth             (p x m x n)        all log-intensities
  # - X_k_coarse_truth      (p x m_est x n)  
  # - X_k_both_truth        (p x m_both x n)
  # - Y_continuous          (n x q_c matrix) 
  #
  # outputs:
  #
  # - result (list of results)
  #
  # ----------------------------------------------------------------------------
  
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
  
  return(result)    
}

# finite basis approach, all helper functions in 21b

simulate_finite_basis_cox_data <- function(n, d, p, adj_type, adj_params, beta_0, 
                                           time_grid, time_grid_est, time_grid_both,
                                           T_max, y_c_query, seed){
  
  #
  #
  #
  # - beta_0        (scalar)                      value for mu(t) which is constant across time and for all processes
  #
  #
  #
  

  
  # 0) get y_c_k
  
  Y_c <- generate_y_c_adj_type(n, adj_type, adj_params, seed = NULL)
  
  # 1) get log intensities for all n subjects and p processes
  
  basis_list <- trig_basis(d)
  basis_mat <- trig_basis_realization(basis_list, time_grid)
    
  cov_mat_list <- lapply(1:n, function(i) {
    trig_basis_cov_mat(d, p, Y_c[i, ], adj_type, adj_params)
  })
  
  # 2) set up for generating X_i(t)
  
 
  mu_t <- rep(beta_0, m)    # mu(t) = beta_0 (constant)
  mu_t_both <- rep(beta_0, length(time_grid_both))
  mu_t_coarse <- rep(beta_0, length(time_grid_est))
  mean_vec <- rep(0, p*d)   # m(t)  = all 0's, where beta ~ N(mean_vec, cov_mat) = pd-dim vec
  
  # 3) for each subject, obtain realizations of log intensities and beta coefficients that got them there
  
  result_both <- trig_basis_log_intensity(cov_mat_list, basis_list, mu_t_both, time_grid_both)
  
  log_intensities_both <- result_both$log_intensities %>% simplify2array()           # (p x m_both x n)
  log_intensities_est <- log_intensities_both[, time_grid_both %in% time_grid_est ,] # (p x m_est x n)
  log_intensities <- log_intensities_both[, time_grid_both %in% time_grid ,]         # (p x m x n)
  
  beta_coeffs <- result_both$beta_coefficients %>% simplify2array() # (p x d x n)

  # 4) Generate point process events
  
  print('at subject generation')
  events <- lapply(1:n, function(i){generate_cox_process_events(log_intensities[, , i], time_grid, T_max, max_intensity = Inf)})
  
  
  
  
  # 5) get list of event_times for each subject (n) and process (p)
  
  event_times_list <- list()
  
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      key <- paste(i, j, sep = "_")
      event_times_list[[key]] <- events[[i]]$event_times[[j]]
    }
  }  
  

  
  # query points
  print('at query point generation')

  cov_mat_query <- lapply(1:nrow(y_c_query), function(i) { 
    trig_basis_cov_mat(d, p, y_c_query[i, ], adj_type, adj_params)
  })
  
  # Store complete subject information
  result <- package_simulation_results(event_times_list, n, p, T_max, y_c_query,
                                       adj_type, adj_params, time_grid, time_grid_est, time_grid_both, seed,
                                       log_intensities, log_intensities_est, log_intensities_both, Y_c,
                                       cov_mat_query)
  
  # ----------------------------------------------------------------------------
  # find ground truths
  # ----------------------------------------------------------------------------
  
  
  print('at rho truth recovery')
  rho_truths <- lapply(1:nrow(y_c_query), function(i) trig_basis_rho_truth(basis_list, mean_vec, time_grid, mu_t, y_c_query[i,], adj_type, adj_params))
  
  g_ij_truth_v2 <- lapply(1:length(cov_mat_query), function(i) trig_basis_cross_covariance_truth(basis_list, cov_mat_query[[i]], time_grid))
  
  print('at eigen truth recovery')
  G <- trig_basis_gram_matrix(basis_list, 0, T_max)
  eigen_truths <- lapply(1:length(cov_mat_query), function(i) trig_basis_eigendecomposition(G, cov_mat_query[[i]], basis_list, time_grid, beta_coeffs))
  
  # ----------------------------------------------------------------------------
  # merge truths - layer 1 = item - layer 2 = y_c_query
  # ----------------------------------------------------------------------------
  
  step_1 <- list(X_k_truth = log_intensities,
                 X_k_coarse_truth = log_intensities_est,
                 X_k_both_truth = log_intensities_both,
                 Lambda_k_truth = exp(log_intensities),
                 Lambda_k_coarse_truth = exp(log_intensities_est),
                 Lambda_k_both_truth = exp(log_intensities_both),
                 mu_t_truth = mu_t,
                 mu_t_both_truth = mu_t_both,
                 mu_t_coarse_truth = mu_t_coarse)
  
  step_2 <- lapply(rho_truths, function(x) {
    list(rho_i_truth = x$rho_i_truth)
  })
  
  step_2b <- lapply(rho_truths, function(x) {
    list(rho_ii_truth = x$rho_ij_truth)
  })
  
  step_3 <- lapply(1:length(rho_truths), function(x) {
    list(g_ij_truth = rho_truths[[x]]$rho_ij_truth,
         g_ij_truth_v2 = g_ij_truth_v2[[x]])
  })
  
  step_4 <- lapply(eigen_truths, function(x) {
    list(eigen_decomp_truth = x$eigen_decomp)
  })
  
  step_5 <- lapply(eigen_truths, function(x) {
    list(KL_cov_truth = x$KL_cov)
  })  
  
  step_9 <- lapply(eigen_truths, function(x) {
    list(C_cond_truth_full = x$corr_op)
  }) 
  


  
  all_truths <- list(step_1 = step_1,
                     step_2 = step_2,
                     step_2b = step_2b,
                     step_3 = step_3,
                     step_4 = step_4,
                     step_5 = step_5,
                     step_9 = step_9)
  
  # rename to be steps
  
  
  
  return(list(dataset = result,
              all_truths = all_truths)) 
}





