source('functions/21b_generate_finite_basis_expansion.R')
source('functions/22_generate_random_variables.R')
source('functions/25_point_generation_process.R')



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
                                       X_k_truth, X_k_coarse_truth, X_k_both_truth, Y_continuous, beta_coeffs = NULL){
  
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
    
    beta_coeffs = beta_coeffs,
    
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
                                           T_max, y_c_query, min_limit, max_limit, seed){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: simulate finite basis data
  #
  # 
  # inputs:
  #
  # - n                 (scalar)      number of replicates
  # - d                 (scalar)      number of eigenfunctions
  # - p                 (scalar)      number of processes
  # - adj_type          (string)      "block_banded_v2"
  # - adj_params        (vector)      parameters associated with the adj_type
  # - beta_0            (scalar)      value for mu(t) which is constant across time and for all processes
  # - time_grid         (m-dim vec)
  # - time_grid_est     (m_est-dim vec)
  # - time_grid_both    (m_both-dim vcec)
  # - T_max             (scalar)    
  # - y_c_query         (n_query x q_c dim vec)
  # - seed              (integer)
  #
  #
  # outputs:
  #
  # - list of 2 items:
  #
  #   - dataset
  #   - all_truths
  #
  #
  # ----------------------------------------------------------------------------
  
  
  
  # 0) get y_c_k
  
  Y_c <- generate_y_c_adj_type(n, adj_type, adj_params, seed = NULL)
  
  # 1) generate cov_mat (pd x pd) for all n subjects according to their y_c_k
  
  basis_list <- trig_basis(d)
  basis_mat <- trig_basis_realization(basis_list, time_grid)
  
  cov_mat_list <- lapply(1:n, function(i) {
    trig_basis_cov_mat(d, p, Y_c[i, ], adj_type, adj_params)
  })
  
  # 2) set up for generating X_i(t)
  
  m <- length(time_grid)
  m_est <- length(time_grid_est)
  m_both <- length(time_grid_both)
  
  mu_t        <- rep(beta_0, m)    # mu(t) = beta_0 (constant)
  mu_t_both   <- rep(beta_0, m_both)
  mu_t_coarse <- rep(beta_0, m_est)
  
  mean_vec <- rep(0, p*d)   # m(t)  = all 0's, where beta ~ N(mean_vec, cov_mat) = pd-dim vec
  
  # 3) for each subject, obtain realizations of log intensities and beta coefficients that got them there
  
  result_both <- trig_basis_log_intensity(cov_mat_list, basis_list, mu_t_both, time_grid_both)
  
  log_intensities_both <- result_both$log_intensities %>% simplify2array()           # (p x m_both x n)
  log_intensities_est <- log_intensities_both[, time_grid_both %in% time_grid_est ,] # (p x m_est x n)
  log_intensities <- log_intensities_both[, time_grid_both %in% time_grid ,]         # (p x m x n)
  
  beta_coeffs <- result_both$beta_coefficients %>% simplify2array() # (p x d x n)
  
  
  
  # 4) Generate point process events
  max_events = Inf
  min_events = -1
  while(min_events < min_limit | max_events > max_limit){
    
    print('at event time generation')
    events <- lapply(1:n, function(i){generate_cox_process_events(log_intensities[, , i], time_grid, T_max, max_intensity = Inf)})
    
    print('finished event time generation')
    max_events <- max(sapply(events, function(i) max(i$event_counts)))
    min_events <- min(sapply(events, function(i) min(i$event_counts)))
    
    print(paste0('most events on a process: ', max_events))
    print(paste0('least events on a process: ', min_events))
  }
  
  
  
  # 5) get list of event_times for each subject (n) and process (p)
  
  event_times_list <- list()
  
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      key <- paste(i, j, sep = "_")
      event_times_list[[key]] <- events[[i]]$event_times[[j]]
    }
  }  
  
  
  
  # 6) query points
  print('at query point generation')
  
  cov_mat_query <- lapply(1:nrow(y_c_query), function(i) { 
    trig_basis_cov_mat(d, p, y_c_query[i, ], adj_type, adj_params)
  })
  
  cor_mat_query <- lapply(cov_mat_query, function(x) assemble_blockwise_correlation(x, p, d))
  prec_mat_query <- lapply(cor_mat_query, function(x) sym(solve(x)))
  
  
  
  
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
  
  # g_ij_truth_v2 <- lapply(1:length(cov_mat_query), function(i) trig_basis_cross_covariance_truth(basis_list, cov_mat_query[[i]], time_grid))
  
  print('at eigen truth recovery')
  G <- trig_basis_gram_matrix(basis_list, 0, T_max)
  eigen_truths <- lapply(1:length(cov_mat_query), function(i) trig_basis_eigendecomposition(G, cov_mat_query[[i]], cor_mat_query[[i]], prec_mat_query[[i]], basis_list, time_grid))
  
  # ----------------------------------------------------------------------------
  # find eigentruths with realized beta values and ground truth eigenfunctions
  # ----------------------------------------------------------------------------
  
  KL_X_truth <- lapply(1:nrow(y_c_query), function(i) { 
    weights_i <- KDE_weights(Y_c, y_c_query[i, ])
    
    X <- aperm(beta_coeffs, c(2, 1, 3)) %>% 
      matrix(nrow = p*d, ncol = n) %>% # reshape to p*d × n
      t()  # (n x pd)
    
    # weighted mean (length p*d)
    mu <- colSums(weights_i * X)
    
    # centered data
    XC <- sweep(X, 2, mu)
    
    # weighted covariance: sum_i w_i (x_i - mu)(x_i - mu)^T
    cov_w <- t(XC * weights_i) %*% XC
    
    kappa <- 1 - sum(weights_i^2)
    cov_w_unbiased <- cov_w / kappa 
    
    # convert to correlation using blockwise function
    corr_w_unbiased <- assemble_blockwise_correlation(cov_w_unbiased, p, d)
    prec_w_unbiased <- sym(solve(corr_w_unbiased))
    
    list(KL_cov = cov_w_unbiased,
         KL_corr = corr_w_unbiased,
         KL_prec = prec_w_unbiased)
    
  })
  
  # 6c) compute C_cond_X_truth 
  
  
  step_9_10_11_X_truth <- lapply(1:length(eigen_truths), function(i){
    
    eigen_decomp_truth    <- eigen_truths[[i]]$eigen_decomp
    KL_cor_X_truth        <- KL_X_truth[[i]]$KL_corr %>% extract_block_structure_v2(p, d)
    KL_prec_X_truth       <- KL_X_truth[[i]]$KL_prec %>% extract_block_structure_v2(p, d)
    C_cond_list_X_truth   <- correlation_estimation_KL_cor(eigen_decomp_truth, KL_cor_X_truth) 
    P_cond_list_X_truth   <- correlation_estimation_KL_cor(eigen_decomp_truth, KL_prec_X_truth) 
    
    C_cond_X_truth_full          <- C_cond_list_X_truth$C_cond %>% assemble_block_matrix_v2(p, m)
    C_cond_X_truth_unnorm_full   <- C_cond_list_X_truth$C_cond_unnorm %>% assemble_block_matrix_v2(p, m)
    
    # step 10
    P_cond_X_truth_full          <- P_cond_list_X_truth$C_cond %>% assemble_block_matrix_v2(p, m)
    P_cond_X_truth_unnorm_full   <- P_cond_list_X_truth$C_cond_unnorm %>% assemble_block_matrix_v2(p, m)
    
    # step 11
    C_HS_X_truth <- hilbert_schmidt_norm_pm(KL_X_truth[[i]]$KL_corr, p, d)
    P_HS_X_truth <- hilbert_schmidt_norm_pm(KL_X_truth[[i]]$KL_prec, p, d)
    
    list(C_cond_X_truth_full = C_cond_X_truth_full,
         C_cond_X_truth_unnorm_full = C_cond_X_truth_unnorm_full,
         P_cond_X_truth_full = P_cond_X_truth_full,
         P_cond_X_truth_unnorm_full = P_cond_X_truth_unnorm_full,
         C_HS_X_truth = C_HS_X_truth,
         P_HS_X_truth = P_HS_X_truth)
  })  
  
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
  
  step_3 <- lapply(rho_truths, function(x) {
    list(g_ij_truth = x$g_ij_truth)
  })
  
  step_4 <- lapply(eigen_truths, function(x) {
    list(eigen_decomp_truth = x$eigen_decomp)
  })
  
  step_5 <- lapply(1:length(eigen_truths), function(i) {
    list(KL_coeffs_truth = beta_coeffs,                          # (p x d x n)
         KL_cov_truth = eigen_truths[[i]]$KL_cov,
         KL_cov_X_truth = KL_X_truth[[i]]$KL_cov %>% extract_block_structure_v2(p, d))    # (pc2 list of dxd matrices)
  })  
  
  step_5b <- lapply(1:length(eigen_truths), function(i) {
    list(KL_cor_truth    = eigen_truths[[i]]$KL_cor,               # (pc2 list of dxd matrices)
         KL_cor_X_truth  = KL_X_truth[[i]]$KL_corr %>% extract_block_structure_v2(p, d),
         KL_prec_truth   = eigen_truths[[i]]$KL_prec,
         KL_prec_X_truth = KL_X_truth[[i]]$KL_prec %>% extract_block_structure_v2(p, d))              # (pc2 list of dxd matrices)
    
  })  
  
  step_9 <- lapply(1:length(eigen_truths), function(i) {
    list(C_cond_truth_full          = eigen_truths[[i]]$C_cond_full,
         C_cond_truth_unnorm_full   = eigen_truths[[i]]$C_cond_full_unnorm,
         C_cond_X_truth_full        = step_9_10_11_X_truth[[i]]$C_cond_X_truth_full,
         C_cond_X_truth_unnorm_full = step_9_10_11_X_truth[[i]]$C_cond_X_truth_unnorm_full)          # (pm x pm matrix)
  }) 
  
  step_9b <- lapply(eigen_truths, function(x) {
    list(efunc_outer_truth = x$efunc_outer,
         efunc_outer_unnorm_truth = x$efunc_outer_unnorm)          # (pc2 list of mxm matrices)
  }) 
  
  step_10 <- lapply(1:length(eigen_truths), function(i) {
    list(P_cond_truth_full           = eigen_truths[[i]]$P_cond_full,
         P_cond_truth_unnorm_full    = eigen_truths[[i]]$P_cond_full_unnorm,
         P_cond_X_truth_full         = step_9_10_11_X_truth[[i]]$P_cond_X_truth_full,
         P_cond_X_truth_unnorm_full  = step_9_10_11_X_truth[[i]]$P_cond_X_truth_unnorm_full)          # (pm x pm matrix)
  }) 
  
  step_11 <- lapply(1:length(eigen_truths), function(i) {
    list(w_mat_truth         = eigen_truths[[i]]$P_HS,
         C_HS_truth          = eigen_truths[[i]]$C_HS,
         w_mat_truth_unnorm  = eigen_truths[[i]]$P_HS_unnorm,
         C_HS_truth_unnorm   = eigen_truths[[i]]$C_HS_unnorm,
         w_mat_X_truth       = step_9_10_11_X_truth[[i]]$P_HS_X_truth,
         C_HS_X_truth        = step_9_10_11_X_truth[[i]]$C_HS_X_truth)                       # (pxp matrix)
  }) 
  
  
  all_truths <- list(step_1 = step_1,
                     step_2 = step_2,
                     step_2b = step_2b,
                     step_3 = step_3,
                     step_4 = step_4,
                     step_5 = step_5,
                     step_5b = step_5b,
                     step_9 = step_9,
                     step_9b = step_9b,
                     step_10 = step_10,
                     step_11 = step_11)
  
  # rename to be steps
  
  
  
  return(list(dataset = result,
              all_truths = all_truths)) 
}


# simulate log-intensities
simulate_finite_basis_cox_data_part1 <- function(temp_file_dir, setting_info_list, group_idx, n_group){
  
  # generate log-intensities for each batch
  
  # 0) load
  
  list2env(setting_info_list, envir = environment())
  
  part_0_info_list <- paste0("part0_", adj_type, '_n_', n, '.rds')
  results <- readRDS(file.path(temp_file_dir, part_0_info_list))
  list2env(results, envir = environment())
  
  # 1) with y_c, n_group, group_idx, find the specific n's that we care about
  
  n_start <- n_group * (group_idx - 1) + 1
  n_end <- n_group * (group_idx) 
  
  n_batch <- n_start:n_end
  
  cov_mat_list <- lapply(n_batch, function(i) {
    trig_basis_cov_mat(d, p, Y_c[i, ], adj_type, adj_params)
  })
  
  # 2) set up for generating X_i(t)
  
  m <- length(time_grid)
  m_est <- length(time_grid_est)
  m_both <- length(time_grid_both)
  
  mu_t        <- rep(beta_0, m)    # mu(t) = beta_0 (constant)
  mu_t_both   <- rep(beta_0, m_both)
  mu_t_coarse <- rep(beta_0, m_est)
  
  mean_vec <- rep(0, p*d)   # m(t)  = all 0's, where beta ~ N(mean_vec, cov_mat) = pd-dim vec
  
  # 3) for each subject, obtain realizations of log intensities and beta coefficients that got them there
  
  result_both <- trig_basis_log_intensity(cov_mat_list, basis_list, mu_t_both, time_grid_both)
  
  log_intensities_both <- result_both$log_intensities %>% simplify2array()           # (p x m_both x n)
  log_intensities_est <- log_intensities_both[, time_grid_both %in% time_grid_est ,] # (p x m_est x n)
  log_intensities <- log_intensities_both[, time_grid_both %in% time_grid ,]         # (p x m x n)
  
  beta_coeffs <- result_both$beta_coefficients %>% simplify2array() # (p x d x n)
  
  
  # save 
  
  out_list <- list(
    cov_mat_list = cov_mat_list,
    log_intensities_both = log_intensities_both,
    log_intensities_est = log_intensities_est,
    log_intensities = log_intensities,
    beta_coeffs = beta_coeffs
  )
  
  part_1_info_list <- paste0("part1_", adj_type, '_n_', n, '_group', group_idx, '.rds')
  
  
  saveRDS(out_list, file = file.path(temp_file_dir, part_1_info_list))    
  
}

# simulate events
simulate_finite_basis_cox_data_part2 <- function(temp_file_dir, setting_info_list, group_idx, n_group, min_events, max_events){
  
  # generate events for n_group 
  
  list2env(setting_info_list, envir = environment())
  
  # load items from part 0
  # - T_max

  part_0_info_list <- paste0("part0_", adj_type, '_n_', n, '.rds')
  results <- readRDS(file.path(temp_file_dir, part_0_info_list))
  list2env(results, envir = environment())
  
  # load items from part 1
  # - cov_mat_list
  # - log_intensities_both
  # - log_intensities_est
  # - log_intensities
  # - beta_coeffs
  part_1_info_list <- paste0("part1_", adj_type, '_n_', n, '_group', group_idx, '.rds')
  results <- readRDS(file.path(temp_file_dir, part_1_info_list))
  list2env(results, envir = environment())
  
  
  
  # 4) Generate point process events
  # max_events_obs = Inf
  # min_events_obs = -1
  # while(min_events_obs < min_events | max_events_obs > max_events){
    
  events <- lapply(1:n_group, function(i){generate_cox_process_events(log_intensities[, , i], time_grid, T_max, max_intensity = Inf)})
  
  max_events_obs <- max(sapply(events, function(i) max(i$event_counts)))
  min_events_obs <- min(sapply(events, function(i) min(i$event_counts)))
  #   print(max_events_obs)
  #   print(min_events_obs)
  # }
  
  # save 
  
  part_2_info_list <- paste0('events_', adj_type, '_n_', n, '_group', group_idx, '.rds')
  saveRDS(events, file = file.path(temp_file_dir, part_2_info_list))  
  
}

# merge data
simulate_finite_basis_cox_data_part3 <- function(temp_file_dir, setting_info_list, group_nums){
  
  # merge all the events and store in results
  
  # 0) load 
  
  # 0a) load the setting values
  list2env(setting_info_list, envir = environment())
  
  # 0b) load part 0's stuff
  part_0_info_list <- paste0("part0_", adj_type, '_n_', n, '.rds')
  results <- readRDS(file.path(temp_file_dir, part_0_info_list))
  list2env(results, envir = environment())
  
  # 0b) load all of the log-intensities and group them
  part_1_info_lists <- paste0(temp_file_dir, "/part1_", adj_type, '_n_', n, '_group', 1:group_nums, '.rds')
  all_part_1_loaded <- lapply(part_1_info_lists, readRDS)
  
  cov_mat_list <- do.call(c, lapply(all_part_1_loaded, `[[`, "cov_mat_list"))
  log_intensities_both <- abind(lapply(all_part_1_loaded, `[[`, "log_intensities_both"), along = 3)
  log_intensities_est  <- abind(lapply(all_part_1_loaded, `[[`, "log_intensities_est"), along = 3)
  log_intensities      <- abind(lapply(all_part_1_loaded, `[[`, "log_intensities"), along = 3)
  beta_coeffs          <- abind(lapply(all_part_1_loaded, `[[`, "beta_coeffs"), along = 3)  
  

  # 0c) load all of the events and group them
  part_2_info_lists <- paste0(temp_file_dir, '/events_', adj_type, '_n_', n, '_group', 1:group_nums, '.rds')
  all_part_2_loaded <- lapply(part_2_info_lists, readRDS)
  events <- do.call(c, all_part_2_loaded)
  
  # ------------------------
  # resume regular function
  # ------------------------
  
  event_times_list <- list()
  
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      key <- paste(i, j, sep = "_")
      event_times_list[[key]] <- events[[i]]$event_times[[j]]
    }
  }  
  
  

  # Store complete subject information
  dataset <- package_simulation_results(event_times_list, n, p, T_max, y_c_query,
                                        adj_type, adj_params, time_grid, time_grid_est, time_grid_both, seed,
                                        log_intensities, log_intensities_est, log_intensities_both, Y_c, beta_coeffs)
                                       
  part_3_info_list <- paste0('dataset_', adj_type, '_n_', n, '.rds')
  saveRDS(dataset, file = file.path(temp_file_dir, part_3_info_list))    

}

# get truths
simulate_finite_basis_cox_data_part4 <- function(temp_file_dir, setting_info_list, cont_ind){
  
  # obtain truths for y_c_query_k
  
  # 0) load 
  list2env(setting_info_list, envir = environment())
  
  
  # 0a) load part0 

  part_0_info_list <- paste0("part0_", adj_type, '_n_', n, '.rds')
  results <- readRDS(file.path(temp_file_dir, part_0_info_list))
  list2env(results, envir = environment())
  
  # 0b) load things from `dataset`
  
  part_3_info_list <- paste0('dataset_', adj_type, '_n_', n, '.rds')
  dataset <- readRDS(file.path(temp_file_dir, part_3_info_list))   
  beta_coeffs <- dataset$beta_coeffs
  
  log_intensities <- dataset$X_k_truth
  log_intensities_est <- dataset$X_k_coarse_truth
  log_intensities_both <- dataset$X_k_both_truth
  

  rm(dataset)
  
  # 1) estimation
  
  cov_mat_query <- trig_basis_cov_mat(d, p, y_c_query[cont_ind, ], adj_type, adj_params)
  cor_mat_query <- assemble_blockwise_correlation(cov_mat_query, p, d)
  prec_mat_query <- sym(solve(cor_mat_query))
  
  
  rho_truths <- trig_basis_rho_truth(basis_list, mean_vec, time_grid, mu_t, y_c_query[cont_ind,], adj_type, adj_params)
  

  G <- trig_basis_gram_matrix(basis_list, 0, T_max)
  eigen_truths <- trig_basis_eigendecomposition(G, cov_mat_query, cor_mat_query, prec_mat_query, basis_list, time_grid)
  
  # ----------------------------------------------------------------------------
  # find eigentruths with realized beta values and ground truth eigenfunctions
  # ----------------------------------------------------------------------------
  

  weights_k <- KDE_weights(Y_c, y_c_query[cont_ind, ])
  
  X <- aperm(beta_coeffs, c(2, 1, 3)) %>% 
    matrix(nrow = p*d, ncol = n) %>% # reshape to p*d × n
    t()  # (n x pd)
  
  # weighted mean (length p*d)
  mu <- colSums(weights_k * X)
  
  # centered data
  XC <- sweep(X, 2, mu)
  
  # weighted covariance: sum_i w_i (x_i - mu)(x_i - mu)^T
  cov_w <- t(XC * weights_k) %*% XC
  
  kappa <- 1 - sum(weights_k^2)
  cov_w_unbiased <- cov_w / kappa 
  
  # convert to correlation using blockwise function
  corr_w_unbiased <- assemble_blockwise_correlation(cov_w_unbiased, p, d)
  prec_w_unbiased <- sym(solve(corr_w_unbiased))
  
  KL_X_truth <- list(KL_cov = cov_w_unbiased,
                     KL_corr = corr_w_unbiased,
                     KL_prec = prec_w_unbiased)
         
    

  
  # 6c) compute C_cond_X_truth 
  
  

    
  eigen_decomp_truth    <- eigen_truths$eigen_decomp
  KL_cor_X_truth        <- KL_X_truth$KL_corr %>% extract_block_structure_v2(p, d)
  KL_prec_X_truth       <- KL_X_truth$KL_prec %>% extract_block_structure_v2(p, d)
  C_cond_list_X_truth   <- correlation_estimation_KL_cor(eigen_decomp_truth, KL_cor_X_truth) 
  P_cond_list_X_truth   <- correlation_estimation_KL_cor(eigen_decomp_truth, KL_prec_X_truth) 
  
  C_cond_X_truth_full          <- C_cond_list_X_truth$C_cond %>% assemble_block_matrix_v2(p, m)
  C_cond_X_truth_unnorm_full   <- C_cond_list_X_truth$C_cond_unnorm %>% assemble_block_matrix_v2(p, m)
  
  # step 10
  P_cond_X_truth_full          <- P_cond_list_X_truth$C_cond %>% assemble_block_matrix_v2(p, m)
  P_cond_X_truth_unnorm_full   <- P_cond_list_X_truth$C_cond_unnorm %>% assemble_block_matrix_v2(p, m)
  
  # step 11
  C_HS_X_truth <- hilbert_schmidt_norm_pm(KL_X_truth$KL_corr, p, d)
  P_HS_X_truth <- hilbert_schmidt_norm_pm(KL_X_truth$KL_prec, p, d)
  
  step_9_10_11_X_truth <- list(C_cond_X_truth_full = C_cond_X_truth_full,
                               C_cond_X_truth_unnorm_full = C_cond_X_truth_unnorm_full,
                               P_cond_X_truth_full = P_cond_X_truth_full,
                               P_cond_X_truth_unnorm_full = P_cond_X_truth_unnorm_full,
                               C_HS_X_truth = C_HS_X_truth,
                               P_HS_X_truth = P_HS_X_truth)
         
 
  
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
  
  step_2 <- list(rho_i_truth = rho_truths$rho_i_truth)

  
  step_2b <- list(rho_ii_truth = rho_truths$rho_ij_truth)

  
  step_3 <- list(g_ij_truth = rho_truths$g_ij_truth)
  
  step_4 <- list(eigen_decomp_truth = eigen_truths$eigen_decomp)
  
  step_5 <- list(KL_coeffs_truth = beta_coeffs,                          # (p x d x n)
                 KL_cov_truth = eigen_truths$KL_cov,
                 KL_cov_X_truth = KL_X_truth$KL_cov %>% extract_block_structure_v2(p, d))    # (pc2 list of dxd matrices)

  
  step_5b <- list(KL_cor_truth    = eigen_truths$KL_cor,               # (pc2 list of dxd matrices)
                  KL_cor_X_truth  = KL_X_truth$KL_corr %>% extract_block_structure_v2(p, d),
                  KL_prec_truth   = eigen_truths$KL_prec,
                  KL_prec_X_truth = KL_X_truth$KL_prec %>% extract_block_structure_v2(p, d))              # (pc2 list of dxd matrices)
    
    

  
  step_9 <- list(C_cond_truth_full          = eigen_truths$C_cond_full,
                 C_cond_truth_unnorm_full   = eigen_truths$C_cond_full_unnorm,
                 C_cond_X_truth_full        = step_9_10_11_X_truth$C_cond_X_truth_full,
                 C_cond_X_truth_unnorm_full = step_9_10_11_X_truth$C_cond_X_truth_unnorm_full)          # (pm x pm matrix)
  
  step_9b <- list(efunc_outer_truth = eigen_truths$efunc_outer,
                  efunc_outer_unnorm_truth = eigen_truths$efunc_outer_unnorm)          # (pc2 list of mxm matrices)
  
  step_10 <- list(P_cond_truth_full           = eigen_truths$P_cond_full,
                  P_cond_truth_unnorm_full    = eigen_truths$P_cond_full_unnorm,
                  P_cond_X_truth_full         = step_9_10_11_X_truth$P_cond_X_truth_full,
                  P_cond_X_truth_unnorm_full  = step_9_10_11_X_truth$P_cond_X_truth_unnorm_full)          # (pm x pm matrix)
  
  step_11 <- list(w_mat_truth         = eigen_truths$P_HS,
                  C_HS_truth          = eigen_truths$C_HS,
                  w_mat_truth_unnorm  = eigen_truths$P_HS_unnorm,
                  C_HS_truth_unnorm   = eigen_truths$C_HS_unnorm,
                  w_mat_X_truth       = step_9_10_11_X_truth$P_HS_X_truth,
                  C_HS_X_truth        = step_9_10_11_X_truth$C_HS_X_truth)                       # (pxp matrix)
  
  
  all_truths <- list(step_1 = step_1,
                     step_2 = step_2,
                     step_2b = step_2b,
                     step_3 = step_3,
                     step_4 = step_4,
                     step_5 = step_5,
                     step_5b = step_5b,
                     step_9 = step_9,
                     step_9b = step_9b,
                     step_10 = step_10,
                     step_11 = step_11)
  
  # store 
  
  part_4_info_list <- paste0('truths_', adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  saveRDS(all_truths, file = file.path(temp_file_dir, part_4_info_list))  
}

# merge truths and data
simulate_finite_basis_cox_data_part5 <- function(temp_file_dir, setting_info_list, cont_inds, group_nums){
  
  # merge all truths together, and merge this with all the results
  
  # 0) load 
  
  # 0a) load the setting values
  list2env(setting_info_list, envir = environment())
  
  # 0b) load the dataset
  
  part_3_info_list <- paste0('dataset_', adj_type, '_n_', n, '.rds') 
  dataset <- readRDS(file.path(temp_file_dir, part_3_info_list)) 
  
  # 0c) load the truths
  part_4_info_lists <- paste0(temp_file_dir, '/truths_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '.rds')
  all_part_4_loaded <- lapply(part_4_info_lists, readRDS)
  
  # Get all element names (step_1, step_2, ..., step_11)
  step_names <- names(all_part_4_loaded[[1]])
  
  # Build nested structure
  all_truths <- lapply(step_names, function(step) {
    lapply(all_part_4_loaded, `[[`, step)
  })
  
  names(all_truths) <- step_names
  
  # ------------------------
  # REMOVE ALL FILES HERE 
  # ------------------------
  

  part0_file_name          <- paste0('part0_', adj_type, '_n_', n, '.rds')
  part1_file_name          <- paste0('part1_', adj_type, '_n_', n, '_group', 1:group_nums, '.rds')
  events_file_name         <- paste0('events_', adj_type, '_n_', n, '_group', 1:group_nums, '.rds')
  truths_file_name         <- paste0('truths_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '.rds')
  dataset_file_name        <- paste0('dataset', adj_type, '_n_', n, '.rds')
  
  file.remove(file.path(temp_file_dir, part0_file_name)) 
  file.remove(file.path(temp_file_dir, part1_file_name)) 
  file.remove(file.path(temp_file_dir, events_file_name)) 
  file.remove(file.path(temp_file_dir, truths_file_name)) 
  file.remove(file.path(temp_file_dir, dataset_file_name)) 
  
  # -------------
  # return so it can be saved
  # -------------
  
  return(list(dataset = dataset,
              all_truths = all_truths)) 
  
}
