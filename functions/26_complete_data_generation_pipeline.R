source('functions/21b_generate_finite_basis_expansion.R')
source('functions/22_generate_random_variables.R')
source('functions/25_point_generation_process.R')



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
  
  # distribution of event counts
  
  over_10000 <- sum(sapply(event_times_list, function(x) length(x) > 10000))
  over_5000 <- sum(sapply(event_times_list, function(x) length(x) > 5000))
  under_10 <- sum(sapply(event_times_list, function(x) length(x) < 10))
  
  event_summary <- summary(sapply(event_times_list, length))
  
  cat("Data generation completed.\n")
  cat("  Total events:", total_events, "\n")
  cat("  Total replicates:", n, "\n")
  cat("  Total processes:", p, "\n")
  cat("  Total times over 10000:", over_10000, "\n")
  cat("  Total times over 5000:", over_5000, "\n")
  cat("  Total times under 10:", under_10, "\n")
  cat("  Event Count Summary:", names(event_summary), "\n")
  cat("  Event Count Summary:", event_summary, "\n")
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



# simulate log-intensities AND events
simulate_finite_basis_cox_data_parts1_and_2 <- function(temp_file_dir, setting_info_list, group_idx, n_group, min_events, max_events){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: generate log-intensities for each batch
  # 
  #       stored as: temp_data/simu_data/parts1_and_2...
  #
  # inputs:
  #
  # - temp_file_dir         (string)   'temp_data/simu_data'
  # - setting_info_list     (list)     list of 'n' and 'adj_type'
  # - group_idx             (integer)  batch number
  # - n_group               (integer)  number of processes per group
  # - min_events            (integer)  minimum number of events
  # - max_events            (integer)  maximum number of events for a subject's process
  # 
  #
  # outputs:
  #
  # - 'parts1_and_2_block_banded_v2_n_1000_group1.rds'
  # 
  # - results   (list of the following)
  #
  #   - events                    (list of n_group lists --> list of `event_times` and event_counts)
  #     - event_times             (list of p vectors of timestamps)
  #     - event_counts            (p-dim vector of # of events)
  #   - cov_mat_list              (n_group-dim list of pd x pd matrices)
  #   - log_intensities_both      (p x m_both x n_group matrix) 
  #   - log_intensities_est       (p x m_est x n_group matrix)
  #   - log_intensities           (p x m x n_group matrix)
  #   - beta_coeffs               (p x d x n_group matrix)
  #   - counter                   (integer)  how many times did we draw?
  #
  #
  # ----------------------------------------------------------------------------
  
  
  # 0) load
  
  list2env(setting_info_list, envir = environment())
  
  part_0_info_list <- paste0("part0_", adj_type, '_n_', n, '_rep_', rep_i, '.rds')
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
  # 3b) then, generate point process events
  
  max_events_obs = Inf
  min_events_obs = -1
  counter = 0
  while(min_events_obs < min_events | max_events_obs > max_events){
    
    result_both <- trig_basis_log_intensity(cov_mat_list, basis_list, mu_t_both, time_grid_both, mean_vec)
    
    log_intensities_both <- result_both$log_intensities %>% simplify2array()           # (p x m_both x n)
    log_intensities_est <- log_intensities_both[, time_grid_both %in% time_grid_est ,] # (p x m_est x n)
    log_intensities <- log_intensities_both[, time_grid_both %in% time_grid ,]         # (p x m x n)
    
    beta_coeffs <- result_both$beta_coefficients %>% simplify2array() # (p x d x n)
  
    # 3b) Generate point process events

    events <- lapply(1:n_group, function(i){generate_cox_process_events(log_intensities[, , i], time_grid, T_max, max_intensity = Inf)})
    
    max_events_obs <- max(sapply(events, function(i) max(i$event_counts)))
    min_events_obs <- min(sapply(events, function(i) min(i$event_counts)))
    
    counter <- counter + 1
    
    print(paste0('group: ', group_idx, ', counter: ', counter, ', min: ', min_events_obs, ', max: ', max_events_obs))
    
  }
  
  # package
  
  results <- list(events = events,
                  cov_mat_list = cov_mat_list,
                  log_intensities_both = log_intensities_both,
                  log_intensities_est = log_intensities_est,
                  log_intensities = log_intensities,
                  beta_coeffs = beta_coeffs,
                  counter = counter
                  )
  
  # save 
  
  parts_1_and_2_info_list <- paste0('parts1_and_2_', adj_type, '_n_', n, '_group', group_idx, '_rep_', rep_i, '.rds')
  saveRDS(results, file = file.path(temp_file_dir, parts_1_and_2_info_list))  
  
}

# merge data
simulate_finite_basis_cox_data_part3 <- function(temp_file_dir, setting_info_list, group_nums){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: merge all the events and store in results
  #
  #       stored as: temp_data/simu_data/dataset...rds
  # 
  # inputs:
  #
  # - temp_file_dir
  # - setting_info_list
  # - group_nums            (integer)  how many groups
  #
  # 
  # outputs:
  #
  # - dataset    (list of the following)
  #
  #   - event_times (pxn length list) each entry is 'k_i', where k = 1, ..., n, and i = 1, ...  p
  #                                   each entry is a list of timestamps of events for subject k on process i
  #   - X_k_truth           (p x m x n matrix)
  #   - X_k_coarse_truth    (p x m_est x n matrix)
  #   - X_k_both_truth      (p x m_both x n matrix)
  #   - Y_continuous        (n x q_c matrix)
  #
  #   - simulation_params  (list of the following)
  #
  #     - n                 (scalar)      number of replicates
  #     - p                 (scalar)      number of processes
  #     - T_max             (scalar)  
  #     - query_y_cs        (n_query x q_c dim vec)
  #     - adj_type          (string)      "block_banded_v2"
  #     - adj_params        (vector)      parameters associated with the adj_type
  #     - time_grid         (m-dim vec)
  #     - time_grid_est     (m_est-dim vec)
  #     - time_grid_both    (m_both-dim vcec)
  #     - seed              (integer)
  #
  #   - beta_coeffs         (p x d x n matrix)
  #
  # ----------------------------------------------------------------------------
  
  # 0) load 
  
  # 0a) load the setting values
  list2env(setting_info_list, envir = environment())
  
  # 0b) load part 0's stuff
  part_0_info_list <- paste0("part0_", adj_type, '_n_', n, '_rep_', rep_i, '.rds')
  results <- readRDS(file.path(temp_file_dir, part_0_info_list))
  list2env(results, envir = environment())
  
  
  # 0bc) load parts 1 and 2 together
  part_1_and_2_info_lists <- paste0(temp_file_dir, '/parts1_and_2_', adj_type, '_n_', n, '_group', 1:group_nums, '_rep_', rep_i, '.rds')
  all_parts_1_and_2_loaded <- lapply(part_1_and_2_info_lists, readRDS)
  
  events               <- do.call(c, lapply(all_parts_1_and_2_loaded, `[[`, "events"))
  cov_mat_list         <- do.call(c, lapply(all_parts_1_and_2_loaded, `[[`, "cov_mat_list"))
  log_intensities_both <- abind(lapply(all_parts_1_and_2_loaded, `[[`, "log_intensities_both"), along = 3)
  log_intensities_est  <- abind(lapply(all_parts_1_and_2_loaded, `[[`, "log_intensities_est"), along = 3)
  log_intensities      <- abind(lapply(all_parts_1_and_2_loaded, `[[`, "log_intensities"), along = 3)
  beta_coeffs          <- abind(lapply(all_parts_1_and_2_loaded, `[[`, "beta_coeffs"), along = 3)  
  total_counter        <- sum(sapply(all_parts_1_and_2_loaded, function(x) x$counter))
  
  cat(paste0('Total Simulation Attempts / Total Groups: ', total_counter, ' / ', group_nums))
  cat('\n')
      

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
                                       
  part_3_info_list <- paste0('dataset_', adj_type, '_n_', n, '_rep_', rep_i, '.rds')
  saveRDS(dataset, file = file.path(temp_file_dir, part_3_info_list))  
  


}

# get truths
simulate_finite_basis_cox_data_part4 <- function(temp_file_dir, setting_info_list, cont_ind){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: obtain truths for y_c_query_k
  #
  #       stored as: temp_data/simu_data/truths...
  # 
  # inputs:
  #
  # - temp_file_dir
  # - setting_info_list
  #
  #   - beta_truth        (boolean)   do we use do estimation with realized beta coefficients and ground truth eigenfunctions?
  # 
  # - cont_ind            (integer)   query number
  #
  # 
  # outputs:
  #
  # - all_truths    (list of the following)
  #
  #   - step_2    (list of 'rho_i_truth' --> p x m matrix)
  #   - step_2b   (list of 'rho_ii_truth' --> list of mxm matrices, one for each i_j)
  #   - step_3    (list of 'g_ij_truth' --> list of mxm matrices, one for each i_j)
  #   - step_4    (list of 'eigen_decomp_truth' --> three smaller lists 'eigenvalues', 'eigenfunctions', 'n_dims')
  #
  #     - eigenvalues  (p-dim list --> d-dim vectors)
  #     - eigenvectors (p-dim list --> mxd matrices)
  #     - n_dims       (p-dim list --> scalars)
  # 
  #   - step_5    (list of 'KL_coeffs_truth', 'KL_cov_truth', and 'KL_cov_X_truth')
  #
  #     - KL_coeffs_truth  (p x d x n matrix)
  #     - KL_cov_truth     (list of dxd matrices, one for each i_j pair)
  #     - KL_cov_X_truth   (list of dxd matrices, one for each i_j pair)
  #
  #   - step_5b   (list of the following)
  #
  #     - KL_cor_truth      (list of dxd matrices, one for each i_j pair)
  #     - KL_cor_X_truth    (list of dxd matrices, one for each i_j pair)
  #     - KL_prec_truth     (list of dxd matrices, one for each i_j pair)
  #     - KL_prec_X_truth   (list of dxd matrices, one for each i_j pair)
  #
  #   - step_9    (list of the following)
  #
  #     - C_cond_truth_full          (pm x pm matrix)
  #     - C_cond_truth_unnorm_full   (pm x pm matrix)
  #     - C_cond_X_truth_full        (pm x pm matrix)
  #     - C_cond_X_truth_unnorm_full (pm x pm matrix)
  #
  #   - step_9b  (list of the following)
  #
  #     - efunc_outer_truth         (list of mxm matrices, one for each i_j pair)
  #     - efunc_outer_unnorm_truth  (list of mxm matrices, one for each i_j pair)
  # 
  #   - step_10  (list of the following)
  #
  #     - P_cond_truth_full          (pm x pm matrix)
  #     - P_cond_truth_unnorm_full   (pm x pm matrix)
  #     - P_cond_X_truth_full        (pm x pm matrix)
  #     - P_cond_X_truth_unnorm_full (pm x pm matrix)
  # 
  #   - step_11  (list of the following)
  #
  #     - w_mat_truth          (p x p matrix)
  #     - C_HS_truth           (p x p matrix)
  #     - w_mat_truth_unnorm   (p x p matrix)
  #     - C_HS_truth_unnorm    (p x p matrix)
  #     - w_mat_X_truth        (p x p matrix)
  #     - C_HS_X_truth         (p x p matrix)
  #
  # ----------------------------------------------------------------------------
  
  
  
  # 0) load 
  list2env(setting_info_list, envir = environment())
  
  
  # 0a) load part0 

  part_0_info_list <- paste0("part0_", adj_type, '_n_', n, '_rep_', rep_i, '.rds')
  results <- readRDS(file.path(temp_file_dir, part_0_info_list))
  list2env(results, envir = environment())
  
  # 0b) load things from `dataset`
  
  part_3_info_list <- paste0('dataset_', adj_type, '_n_', n, '_rep_', rep_i, '.rds')
  dataset <- readRDS(file.path(temp_file_dir, part_3_info_list))   
  beta_coeffs <- dataset$beta_coeffs
  
  log_intensities <- dataset$X_k_truth
  log_intensities_est <- dataset$X_k_coarse_truth
  log_intensities_both <- dataset$X_k_both_truth
  

  rm(dataset)
  
  # 1) queried covariance matrix truth
  # 
  # - prec_mat is the inverse of the cor_mat
  
  cov_mat_query <- trig_basis_cov_mat(d, p, y_c_query[cont_ind, ], adj_type, adj_params)  # pd x pd matrix of covariances
  cor_mat_query <- assemble_blockwise_correlation(cov_mat_query, p, d)                    # pd x pd matrix of correlations
  prec_mat_query <- sym(solve(cor_mat_query))                                             # pd x pd matrix of precisions
  
  adj_mat_truth <- trig_basis_adj_mat(d, p, y_c_query[cont_ind, ], adj_type, adj_params, thresh = 1e-3)

  
  # 2) estimation
  
  rho_truths <- trig_basis_rho_truth(basis_list, mean_vec, cov_mat_query, time_grid, mu_t)
  

  G <- trig_basis_gram_matrix(basis_list, 0, T_max)
  eigen_truths <- trig_basis_eigendecomposition(G, cov_mat_query, cor_mat_query, prec_mat_query, basis_list, time_grid)
  
  # ----------------------------------------------------------------------------
  # 3) find eigentruths with realized beta values and ground truth eigenfunctions
  #
  # denote this ax beta_truth
  # 
  # ----------------------------------------------------------------------------
  
  if(beta_truth){
    y_c_query_k <- y_c_query[cont_ind, ]
    eigen_decomp_truth       <- eigen_truths$eigen_decomp
    
    eigen_beta_truths <- trig_basis_eigendecomposition_beta_truths(beta_coeffs, eigen_decomp_truth, Y_c, y_c_query_k)
  }
         
    

  # ----------------------------------------------------------------------------
  # merge truths - layer 1 = item - layer 2 = y_c_query
  # ----------------------------------------------------------------------------
  
  # step_1 <- list(X_k_truth = log_intensities,
  #                X_k_coarse_truth = log_intensities_est,
  #                X_k_both_truth = log_intensities_both,
  #                Lambda_k_truth = exp(log_intensities),
  #                Lambda_k_coarse_truth = exp(log_intensities_est),
  #                Lambda_k_both_truth = exp(log_intensities_both),
  #                mu_t_truth = mu_t,
  #                mu_t_both_truth = mu_t_both,
  #                mu_t_coarse_truth = mu_t_coarse)
  
  step_2 <- list(rho_i_truth = rho_truths$rho_i_truth)
  
  step_2b <- list(rho_ii_truth = rho_truths$rho_ij_truth)

  step_3 <- list(g_ij_truth = rho_truths$g_ij_truth)
  
  step_4 <- list(eigen_decomp_truth = eigen_truths$eigen_decomp)
  
  step_5 <- list(KL_cov_truth = eigen_truths$KL_cov)           # (i_j list of dxd matrices)

  step_5b <- list(KL_cor_truth    = eigen_truths$KL_cor)       # (i_j list of dxd matrices)
  
  step_5c <- list(KL_prec_truth   = eigen_truths$KL_prec)      # (i_j list of dxd matrices)
  
  step_5d <- list(KL_coeffs_truth = beta_coeffs)               # (p x d x n)
    
  step_9 <- list(C_cond_truth                = eigen_truths$C_cond)
                 #C_cond_truth_unnorm         = eigen_truths$C_cond_unnorm)          # (pm x pm matrix)
  
  step_9b <- list(efunc_outer_truth          = eigen_truths$efunc_outer)
                  #efunc_outer_unnorm_truth   = eigen_truths$efunc_outer_unnorm)     # (pc2 list of mxm matrices)
  
  step_10 <- list(P_cond_truth               = eigen_truths$P_cond)
                  #P_cond_truth_unnorm        = eigen_truths$P_cond_unnorm)          # (pm x pm matrix)
  
  step_11 <- list(w_mat_truth         = eigen_truths$P_HS,                          # matrices of HS norms 
                  #w_mat_truth_unnorm  = eigen_truths$P_HS_unnorm,                   # (pxp matrix)
                  w_mat_KL_truth      = eigen_truths$P_HS_KL)                   
  
  step_11b <- list(C_HS_truth          = eigen_truths$C_HS,
                   #C_HS_truth_unnorm   = eigen_truths$C_HS_unnorm,
                   C_HS_KL_truth       = eigen_truths$C_HS_KL)
  
  

  w_mat_truth     <- eigen_truths$P_HS
  w_mat_KL_truth  <- eigen_truths$P_HS_KL
  
  step_12 <- list(roc_truth    = roc_for_raw_w_mat(w_mat_truth, adj_mat_truth),
                  roc_KL_truth = roc_for_raw_w_mat(w_mat_KL_truth, adj_mat_truth))
  
  step_12b <- list(adj_mat_truth = adj_mat_truth)
  
  # true_graphs
  
  true_graphs <- list(adj_mat_truth = adj_mat_truth,
                      cov_mat_truth = cov_mat_query,
                      cor_mat_truth = cor_mat_query,
                      prec_mat_truth = prec_mat_query)
  
  if(beta_truth){
    
    step_5[['KL_cov_beta_truth']]   <- eigen_beta_truths$KL_cov_beta_truth
    step_5b[['KL_cor_beta_truth']]  <- eigen_beta_truths$KL_cor_beta_truth
    step_5c[['KL_prec_beta_truth']] <- eigen_beta_truths$KL_prec_beta_truth
    
    step_9[['C_cond_beta_truth']]        <- eigen_beta_truths$C_cond_beta_truth
    #step_9[['C_cond_beta_truth_unnorm']] <- eigen_beta_truths$C_cond_beta_truth_unnorm

    step_10[['P_cond_beta_truth']]        <- eigen_beta_truths$P_cond_beta_truth
    #step_10[['P_cond_beta_truth_unnorm']] <- eigen_beta_truths$P_cond_beta_truth_unnorm
    

    step_11[['w_mat_beta_truth']] <- eigen_beta_truths$P_HS_beta_truth
    step_11b[['C_HS_beta_truth']] <- eigen_beta_truths$C_HS_beta_truth


    step_12[['roc_beta_truth']] <- roc_for_raw_w_mat(eigen_beta_truths$P_HS_beta_truth, adj_mat_truth)
  }
  
  
  all_truths <- list(step_2 = step_2,
                     step_2b = step_2b,
                     step_3 = step_3,
                     step_4 = step_4,
                     step_5 = step_5,
                     step_5b = step_5b,
                     step_5c = step_5c,
                     step_5d = step_5d,
                     step_9 = step_9,
                     step_9b = step_9b,
                     step_10 = step_10,
                     step_11 = step_11,
                     step_11b = step_11b,
                     step_12 = step_12,
                     step_12b = step_12b,
                     true_graphs = true_graphs
                     )
  
  # store 
  
  part_4_info_list <- paste0('truths_', adj_type, '_n_', n, '_nquery', cont_ind, '_rep_', rep_i, '.rds')
  saveRDS(all_truths, file = file.path(temp_file_dir, part_4_info_list))  
}

# merge truths and data
simulate_finite_basis_cox_data_part5 <- function(temp_file_dir, setting_info_list, cont_inds, group_nums){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: merge all truths together, and merge this with all the results
  #
  # inputs:
  # 
  # - temp_file_dir
  # - setting_info_list
  # - cont_inds            (integer)  how many n_querys
  # - group_nums           (integer)  how many groups
  #
  # 
  # ----------------------------------------------------------------------------
  
  # 0) load 
  
  # 0a) load the setting values
  list2env(setting_info_list, envir = environment())
  
  
  # 0b) load the truths
  part_4_info_lists <- paste0(temp_file_dir, '/truths_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_rep_', rep_i, '.rds')
  all_part_4_loaded <- lapply(part_4_info_lists, readRDS)
  
  # 0c) load part 0 to get query_y_cs
  part_0_info_list <- paste0("part0_", adj_type, '_n_', n, '_rep_', rep_i, '.rds')
  results <- readRDS(file.path(temp_file_dir, part_0_info_list))
  list2env(results, envir = environment())
  
  names(all_part_4_loaded) <- round(y_c_query[,1], 3)
  
  # Get all element names (step_1, step_2, ..., step_11)
  step_names <- names(all_part_4_loaded[[1]])
  
  # Build nested structure
  all_truths <- lapply(step_names, function(step) {
    lapply(all_part_4_loaded, `[[`, step)
  })
  
  names(all_truths) <- step_names
  
  
  # 1) retrieve and insert step_1 (which doesn't vary based on y_c)
  
  # 1a) load the dataset
  
  part_3_info_list <- paste0('dataset_', adj_type, '_n_', n, '_rep_', rep_i, '.rds') 
  dataset <- readRDS(file.path(temp_file_dir, part_3_info_list)) 
  
  log_intensities <- dataset$X_k_truth
  log_intensities_est <- dataset$X_k_coarse_truth
  log_intensities_both <- dataset$X_k_both_truth
  

  
  # 1c) aggregate results
  step_1 <- list(X_k_truth = log_intensities,
                 X_k_coarse_truth = log_intensities_est,
                 X_k_both_truth = log_intensities_both)
  
  step_1b <- list(Lambda_k_truth = exp(log_intensities),
                  Lambda_k_coarse_truth = exp(log_intensities_est),
                  Lambda_k_both_truth = exp(log_intensities_both))
  
  step_1c <- list(mu_t_truth = mu_t,
                  mu_t_both_truth = mu_t_both,
                  mu_t_coarse_truth = mu_t_coarse)
  
  all_truths[['step_1']] <- step_1
  all_truths[['step_1b']] <- step_1b
  all_truths[['step_1c']] <- step_1c
  
  # ------------------------
  # REMOVE ALL FILES HERE 
  # ------------------------
  

  part0_file_name           <- paste0('part0_', adj_type, '_n_', n, '_rep_', rep_i, '.rds')
  parts_1_and_2_file_name   <- paste0('parts1_and_2_', adj_type, '_n_', n, '_group', 1:group_nums, '_rep_', rep_i, '.rds')
  truths_file_name          <- paste0('truths_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_rep_', rep_i, '.rds')
  dataset_file_name         <- paste0('dataset_', adj_type, '_n_', n, '_rep_', rep_i, '.rds')
  
  # file.remove(file.path(temp_file_dir, part0_file_name))
  # file.remove(file.path(temp_file_dir, parts_1_and_2_file_name))
  # file.remove(file.path(temp_file_dir, truths_file_name))
  # file.remove(file.path(temp_file_dir, dataset_file_name))
  
  # -------------
  # return so it can be saved
  # -------------
  
  return(list(dataset = dataset,
              all_truths = all_truths)) 
  
}
