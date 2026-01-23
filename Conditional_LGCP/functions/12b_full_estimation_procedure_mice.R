# truncated from "with truths" for mice data


# DO NOT KEEP TRACK OF TRUTHS THROUGHOUT THE PROCESS

full_conditional_estimation_with_no_truth <- function(dataset, method, ncores, dir, mouse = F){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: regular estimation without keeping track of truths
  #
  # 
  #
  #
  # Supports:
  #   - "OG", "JASA"  → covariate loop after step_5 (v2-style)
  #   - "CPGM"        → covariate loop after step_1 (v3-style)
  #
  # Input:
  #
  # - dataset  (list of the following features)
  #
  #   - event_times         (n*p-dim list)      each item is a vector of timestamps for process i in subject k
  #   - X_k_truth           (p x m x n matrix)
  #   - X_k_coarse_truth    (p x m_est x n matrix)
  #   - X_k_both_truth      (p x m_both x n matrix)
  #   - Y_continuous        (n x q_c matrix)
  #   - simulation_params   (list of params such as n, p, query_y_cs, adj_type, adj_params, time_grids, seed)
  #   - true_graphs         (list of true settings for each query_y_cs)
  #   - summary_stats       (list of summary stats)
  #
  # 
  # - method   ('CPGM')  
  # - ncores   (integer)
  # - dir      (string)     folder name, such as "simu_results/block_banded_v2"
  # - mouse    (boolean)    are we using a mouse?
  # 
  #
  # ----------------------------------------------------------------------------
  
  if (!(method %in% c("CPGM"))) {
    stop("Error 12b: estimation method must be 'CPGM'")
  }
  
  # load  
  if(mouse){
    time_grid_est <- dataset$simulation_params$time_grid_est
    m_est <- length(time_grid_est) 
    n <- dim(dataset$Y_continuous)[1]
  } else{
    time_grid_est <- dataset$simulation_params$time_grid_est
    time_grid <- dataset$simulation_params$time_grid
    time_grid_both <- dataset$simulation_params$time_grid_both
    n <- dim(dataset$Y_continuous)[1]
    m_est <- length(time_grid_est) 
  }

  


  # ----------------------------------------------------------------------------
  # Step 0 - Check and preprocess data
  # ----------------------------------------------------------------------------
  
  step_0_events <- step_0_keep_events(dataset, k = 1, i_vec = 1:5)

  processed_data <- step_0_preprocess(dataset)
  data_df4     <- processed_data[[1]]
  y_c_strata   <- processed_data[[2]]
  query_y_cs   <- processed_data[[3]]
  patient_sel  <- processed_data[[4]]
  feature_sel  <- processed_data[[5]]
  rm(processed_data)
  
  p <- length(feature_sel)
  full <- F
  
  # ----------------------------------------------------------------------------
  # Step 1 - Log intensities
  # ----------------------------------------------------------------------------
  
  step_1 <- step_1_log_intensities(dataset, data_df4, time_grid_est, NA, NA, full)
  
  # start covariate loop 
  # ----------------------------------------------------------------------------
  
  print('at covariate loop')
  cont_inds <- 1:nrow(query_y_cs)
  estimated_graphs <- lapply(cont_inds, function(cont_ind) {
    
    print(paste0(cont_ind, ' out of ', length(cont_inds)))
    
    query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()
    
    if(method %in% c('OG', 'JASA')){
      rho_kernel <- F
      i_neq_j <- F
      query_y_c_step_2 <- 'hold'
      y_c_strata_step_2 <- 'hold'
    } else{
      rho_kernel <- T
      i_neq_j <- T
      query_y_c_step_2 <- query_y_c
      y_c_strata_step_2 <- y_c_strata     
    }
    
    # Step 2: rho_i estimation (per subject)
    
    
    step_2 <- step_2_rho_i(dataset, data_df4, kernel_params_i, rho_kernel,
                           patient_sel, feature_sel, time_grid, time_grid_est,
                           query_y_c_step_2, y_c_strata_step_2, ncores, full)           
    
    
    # Step 2b onward
    step_2b <- step_2_rho_ij(step_1, step_2, kernel_params_i, i_neq_j, full)
    step_3  <- step_3_g_ij(step_2, step_2b, kernel_params_i, i_neq_j, full)
    
    
    
    step_4 <- tryCatch({
      step_4_eigendecomp(step_3, p, time_grid, time_grid_est, full)
    }, error = function(e) {
      cat("Error in step_4, saving dataset...\n")
      save(dataset, file = file.path(dir, "dataset.RData"))
      stop(e)
    })
    
    # step 5 to 9 split:
    
    if(method %in% c('OG', 'JASA')){
      step_5 <- step_5_KL_expansion(step_1, step_4, kernel_params_i, time_grid, time_grid_est, ncores)
      step_8 <- steps_78(step_4, step_5, kernel_params_i, y_c_strata, query_y_c, method, ncores)
      step_9 <- step_9_C_cond_from_V_cond(step_8, kernel_params_i)
    } else{
      step_5 <- step_5_KL_covariance(step_3, step_4, full)
      step_5b <- step_5b_KL_correlation(step_5, p, full)
      step_9 <- step_9_C_cond_from_KL_cor(step_4, step_5b, kernel_params_i, full)
      step_9b <- step_9b_eigenfunction_outers(step_4, full)
    }
    
    # reunite at step 10 onwards
    
    block <- F
    MP <- F
    
    # step_10 <- step_10_P_cond(step_9, kernel_params_i, p, block, MP)
    step_10 <- tryCatch({
      step_10_P_cond(step_9, kernel_params_i, p, block, MP, full)
    }, error = function(e) {
      cat("Error occurred in step_10, saving dataset...\n")
      save(dataset, file = paste0(dir, "/dataset.RData"))
      cat("Dataset saved to dataset.RData\n")
      stop(e)  # Re-throw the error
    })
    
    
    step_11 <- step_11_HS_norms(step_9, step_10, adj_mat_i, p, full)

    
    list(step_2 = step_2, step_2b = step_2b, step_3 = step_3,
         step_4 = step_4, step_5 = step_5, step_5b = step_5b, step_9 = step_9, step_9b = step_9b,
         step_10 = step_10, step_11 = step_11)
  })
  
  # ----------------------------------------------------------------------------
  # Reorganize by steps instead of by y_c_query
  # ----------------------------------------------------------------------------
  
  steps <- unique(unlist(lapply(estimated_graphs, names)))
  reorganized <- setNames(lapply(steps, function(step) {
    sapply(estimated_graphs, `[[`, step, simplify = FALSE)
  }), steps)
  
  # Return step_0, step_1 (shared) + reorganized per-subject steps
  estimated_graphs_part_1 <- list(step_0_events = step_0_events,
                                  step_1 = step_1)
  
  all_results <- c(estimated_graphs_part_1, reorganized)
  all_results$y_c_query <- query_y_cs
  
  if(mouse){
    all_results$time_grid_est <- time_grid_est
  } else{
    all_results$time_grid <- time_grid
    all_results$time_grid_est <- time_grid_est
    all_results$time_grid_both <- time_grid_both    
  }

  all_results$p <- p
  
  return(all_results)
  
}

full_conditional_estimation_with_no_truth_part1 <- function(dataset, setting_info_list, ncores, temp_file_dir, mouse, X_truth){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: bundle all relevant parameters into a list called:
  #
  #       temp_data/simu/part1_block_banded_c0_n_100.rds
  #       temp_data/mice/part1_Tau1_m0vr0_t5.rds
  #
  # inputs:
  #
  # - dataset               (list of the following)
  #   - event_times           (list of vectors)  each vector is named k_i for subject k and process i
  #   - Y_continuous          (n x q_c matrix)   continuous covariates
  #   - simulation_params     (list of various parameters)
  # 
  # - setting_info_list    (list)
  # - ncores               (integer)
  # - temp_file_dir        (string)   'temp_data/simu'
  # - mouse                (boolean)  are we working with mice data?
  # - X_truth              (boolean)  do we want to get estimates starting from true log intensities?
  #
  # 
  # outputs:
  # 
  # - save a list of items that are relevant when we select n out of n_large subjects, including:
  #
  #   - gamma_c (varies with n)
  #   - 
  # 
  # ----------------------------------------------------------------------------

  list2env(setting_info_list, envir = environment())  

  
  
  if (!(method %in% c("CPGM"))) {
    stop("Error 12b: estimation method must be 'CPGM'")
  }
  
  # load  
  if(mouse){
    time_grid_est <- dataset$simulation_params$time_grid_est
    m_est <- length(time_grid_est) 
    n <- dataset$simulation_params$n
  } else{
    time_grid_est <- dataset$simulation_params$time_grid_est
    time_grid <- dataset$simulation_params$time_grid
    time_grid_both <- dataset$simulation_params$time_grid_both
    m_est <- length(time_grid_est) 
  }
  
  
  
  
  # ----------------------------------------------------------------------------
  # Step 0 - Check and preprocess data
  # ----------------------------------------------------------------------------
  
  # load `data_df4`, `y_c_strata_full`, `y_c_strata`, `query_y_cs`, `patient_sel`, `feature_sel`
  processed_data <- step_0_preprocess(dataset)
  list2env(processed_data, envir = environment()) 
  rm(processed_data)
  
  
  p <- length(feature_sel)
  
  step_0_events <- step_0_keep_events(dataset, k_vec = 1:n, i_vec = 1:p)
  


  # ----------------------------------------------------------------------------
  # Step 1 - Log intensities
  # ----------------------------------------------------------------------------
  
  step_1_bundle <- step_1_log_intensities(data_df4, time_grid_est) 
  step_1 <- step_1_bundle$step_1
  step_1b <- step_1_bundle$step_1b
  
  if(X_truth){
    step_1[['X_k_coarse_truth']] <- dataset$X_k_coarse_truth
    step_1[['X_k_truth']]        <- dataset$X_k_truth
    step_1[['X_k_both_truth']]   <- dataset$X_k_both_truth
    
    step_1b[['Lambda_k_coarse_truth']] <- exp(dataset$X_k_coarse_truth)
    step_1b[['Lambda_k_truth']]        <- exp(dataset$X_k_truth)
    step_1b[['Lambda_k_both_truth']]   <- exp(dataset$X_k_both_truth)
  }
  
  # 2) load gamma_c and i_j keys
  
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  
  keys <- expand.grid(i = 1:p, j = 1:p) %>%
    subset(i <= j) %>%
    with(paste0(i, "_", j))
  
  # print
  cat("n_bivariate_processes=", length(keys), "\n")
  cat("n_processes=", p, "\n")
  cat("n_queries=", nrow(query_y_cs), "\n")
  
  
  # Save results for Stage 2
  
  results <- list(
    dataset = dataset,
    step_0_events = step_0_events,
    step_1 = step_1_bundle$step_1,
    step_1b = step_1_bundle$step_1b,
    data_df4 = data_df4,
    y_c_strata_full = y_c_strata_full,
    y_c_strata = y_c_strata,
    query_y_cs = query_y_cs,
    patient_sel = patient_sel,
    feature_sel = feature_sel,
    p = p,
    gamma_c = gamma_c,
    keys = keys,
    mouse = mouse,
    ncores = ncores,
    method = method
  )
  
  if(mouse){
    results[['time_grid_est']] <- time_grid_est
    results[['m_est']] <- m_est
    results[['n']] <- n    
  } else{
    results[['time_grid_est']] <- time_grid_est
    results[['time_grid']] <- time_grid
    results[['time_grid_both']] <- time_grid_both
    results[['m_est']] <- m_est
    results[['n']] <- n      
    results[['d']] <- dim(dataset$beta_coeffs)[2]
    
    if(X_truth){
      results[['X_k_coarse_truth']] <- dataset$X_k_coarse_truth
      results[['X_k_truth']] <- dataset$X_k_truth
      results[['X_k_both_truth']] <- dataset$X_k_both_truth
    }    
    
  }
  
  if(mouse){
    datafile_name <- paste0('part1_', ID, '_', discrete_level, '_t', time_scale, '.rds')
  } else{
    datafile_name <- paste0('part1_', adj_type, '_n_', n, '.rds')
  }
  

  
  saveRDS(results, file = file.path(temp_file_dir, datafile_name))  
  
}

# NOT USED 
full_conditional_estimation_with_no_truth_part2 <- function(temp_file_dir, setting_info_list, cont_ind){

  # load 
  
  
  # ID <- setting_info_list[['ID']]
  # y_c_structure <- setting_info_list[['y_c_structure']]
  # time_scale <- setting_info_list[['time_scale']]
  # method <- setting_info_list[['method']]
  # discrete_level <- setting_info_list[['discrete_level']]
  
  list2env(setting_info_list, envir = environment())
  
  datafile_name <- paste0("part1_", ID, '_', discrete_level, '_t', time_scale, '.rds')
  datafile_error_name <- paste0("dataset_part2_", ID, '_', discrete_level, '_t', time_scale, '.RData')  # in case we need to quit and troubleshoot

  results <- readRDS(file.path(temp_file_dir, datafile_name))
  
  
  # load all variabels
  list2env(results, envir = .GlobalEnv)
  
  # start covariate loop 
  # ----------------------------------------------------------------------------
  
  

  cont_inds <- 1:nrow(query_y_cs)

  full <- F
  
  # Step 0: prep 
  
  print(paste0(cont_ind, ' out of ', length(cont_inds)))
  
  query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()
  
  if(method %in% c('OG', 'JASA')){
    rho_kernel <- F
    i_neq_j <- F
    query_y_c_step_2 <- 'hold'
    y_c_strata_step_2 <- 'hold'
  } else{
    rho_kernel <- T
    i_neq_j <- T
    query_y_c_step_2 <- query_y_c
    y_c_strata_step_2 <- y_c_strata     
  }
  
  # Step 2: rho_i estimation (per subject)
  
  
  step_2 <- step_2_rho_i(dataset, data_df4, kernel_params_i, rho_kernel,
                         patient_sel, feature_sel, time_grid, time_grid_est,
                         query_y_c_step_2, y_c_strata_step_2, ncores, full)           
  
  
  # Step 2b onward
  step_2b <- step_2_rho_ij(step_1, step_2, kernel_params_i, i_neq_j, full)
  step_3  <- step_3_g_ij(step_2, step_2b, kernel_params_i, i_neq_j, full)
  
  
  
  step_4 <- tryCatch({
    step_4_eigendecomp(step_3, p, time_grid, time_grid_est, full)
  }, error = function(e) {
    cat("Error in step_4, saving dataset...\n")
    save(dataset, file = file.path(temp_file_dir, datafile_error_name))
    stop(e)
  })
  
  # step 5 to 9 split:
  
  if(method %in% c('OG', 'JASA')){
    step_5 <- step_5_KL_expansion(step_1, step_4, kernel_params_i, time_grid, time_grid_est, ncores)
    step_8 <- steps_78(step_4, step_5, kernel_params_i, y_c_strata, query_y_c, method, ncores)
    step_9 <- step_9_C_cond_from_V_cond(step_8, kernel_params_i)
  } else{
    step_5 <- step_5_KL_covariance(step_3, step_4, full)
    step_5b <- step_5b_KL_correlation(step_5, p, full)
    step_9 <- step_9_C_cond_from_KL_cor(step_4, step_5b, kernel_params_i, full)
    step_9b <- step_9b_eigenfunction_outers(step_4, full)
  }
  
  # reunite at step 10 onwards
  
  block <- F
  MP <- F
  
  # step_10 <- step_10_P_cond(step_9, kernel_params_i, p, block, MP)
  step_10 <- tryCatch({
    step_10_P_cond(step_9, kernel_params_i, p, block, MP, full)
  }, error = function(e) {
    cat("Error occurred in step_10, saving dataset...\n")
    save(dataset, file = file.path(temp_file_dir, datafile_error_name))
    cat("Dataset saved to dataset.RData\n")
    stop(e)  # Re-throw the error
  })
  
  
  step_11 <- step_11_HS_norms(step_9, step_10, adj_mat_i, p, full)
  
  
  estimated_graphs <- list(step_2 = step_2, step_2b = step_2b, step_3 = step_3,
                           step_4 = step_4, step_5 = step_5, step_5b = step_5b, step_9 = step_9, step_9b = step_9b,
                           step_10 = step_10, step_11 = step_11)
  
  
  # save as part2_1, part2_2
  file_name <- paste0('part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  saveRDS(estimated_graphs, file = file.path(temp_file_dir, file_name))  
   
}


estimate_intensities_stratum_parallel_with_yc_part0 <- function(temp_file_dirs, setting_info_list, cont_ind, mouse){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: we are beginning the estimation method for a specific n_query.
  #       calculate weights and adj_mat (truth), and store it in:
  #
  #       temp_data/simu/part2_block_banded_c0_n_100_nquery1.rds
  #       temp_data/mice/part2_Tau1_m0vr0_t5_nquery1.rds
  # 
  #
  # inputs:
  #
  # - temp_file_dirs      (vector of strings) c('temp_data/simu', 'temp_data/simu_data')
  # - setting_info_list   (list)
  # - cont_ind            (integer) n_query id
  # - mouse               (boolean)  are we working with mice data?
  #
  #
  # outputs:
  #
  # - list of items from part1 as well as:
  #
  #   - weights
  #   - W_y
  #   - adj_mat_i
  # 
  # ----------------------------------------------------------------------------
  
  list2env(setting_info_list, envir = environment())
  
  

  
  if(mouse){
    step_1_info_list <- paste0('part1_', ID, '_', discrete_level, '_t', time_scale, '.rds')
  } else{
    step_1_info_list <- paste0('part1_', adj_type, '_n_', n, '.rds')
  }
  
  results <- readRDS(file.path(temp_file_dirs[1], step_1_info_list))
  list2env(results, envir = environment())
  
  # 1) get the query and get the weights
  
  y_c_query <- query_y_cs[cont_ind,]
  
  weights <- apply(y_c_strata, 1, function(row) {
    step_6_kernel(as.numeric(row), y_c_query, gamma_c) 
  })  
  
  W_y <- sum(weights)
  weights2 <- weights / W_y # normalize
  
  results[['weights']] <- weights2
  results[['W_y']] <- W_y
  
  # 2) get the ground truth adj_mat
  
  if(! mouse){
    truth_data_name <- paste0('truths_', adj_type, '_n_', n_large, '_nquery', cont_ind, '.rds')
    truths <- readRDS(file.path(temp_file_dirs[2], truth_data_name))

    results[['adj_mat_i']] <- truths$true_graphs$adj_mat_truth
  }
  
  
  # 3) save
  
  if(mouse){
    datafile_name <- paste0('part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    datafile_name <- paste0('part2_', adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }
  
  saveRDS(results, file = file.path(temp_file_dirs[1], datafile_name))
  

}

# rho_i
estimate_intensities_stratum_parallel_with_yc_part1 <- function(temp_file_dir, setting_info_list, cont_ind, i, mouse, X_truth) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: calculate rho_i for a single i and save it as
  #       
  #       temp_data/simu/step_2_rho_i_block_banded_v2_n_100_nqueryk_i.rds'
  #
  # 
  # inputs:
  # 
  #
  # - temp_file_dir         (string)
  # - setting_info_list     (list)
  # - cont_ind              (integer)   which n_query index
  # - i                     (integer)   process_id from 1 to p
  # - mouse                 (boolean)   mouse (T) or simulation (F)
  # - X_truth               (boolean)   do we want to do estimation from true log-intensities?
  #
  #
  # outputs:
  #
  # - rho_i_result          (list)
  #
  #   - rho_i_est           (m-dim vector)  rho_i_est
  #   - rho_i_X_truth       (m-dim vector)  only computed if we have X_truth = T
  #
  # ----------------------------------------------------------------------------
  
  list2env(setting_info_list, envir = environment())
  
  # 1) load
  
  if(mouse){
    step_2_info_list <- paste0('part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    step_2_info_list <- paste0('part2_', adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }

  
  # load `dataset$X_k_truth`, `weights`, `data_df4`
  results <- readRDS(file.path(temp_file_dir, step_2_info_list))
  list2env(results, envir = environment())
  
  n_time <- length(time_grid_est)  # 19
  
  
  # 2) filter process i and update subject ID's in case some don't have process i
  
  data_i <- data_df4[feature_id == feature_sel[i], ]
  kept_subjects <- unique(data_i$subject_num) %>% sort()  # in case any subjects do not have data for process i
  reweights = weights[kept_subjects] / sum(weights[kept_subjects])  # recalculate weights in case we discard subjects
  
  setDT(data_i)
  data_i[, subject_num := match(subject_num, sort(unique(subject_num)))]
  
  # 3) estimation
  if (nrow(data_i) == 0) { # no events, estimate is the zero intensity
    rho_i <- rep(0, n_time)
  } else {
    Gamma_i <- data_i[, estimate_density(time, time_grid_est), by = "subject_num"]
    rho_mat <- matrix(Gamma_i$rho_hat, nrow = n_time)
    
    
    rho_mat2 <- sweep(rho_mat, 2, reweights, `*`) # multiply each 19-dim vec by `weight` which was already filtered for n <= n_large and normalized
    
    rho_i <- apply(rho_mat2, 1, sum) # since these are normalized weights, just add them
  }
  
  # store rho_i
  
  rho_i_result <- list(rho_i_est = rho_i) 
  
  # 4) X_truth if requested
  
  if(X_truth){
    # take (12 x 30 x 100), index only the i-th process, then take sample mean across n
    rho_i_result[['rho_i_X_truth']] <- as.vector(exp(dataset$X_k_truth[i, , ]) %*% reweights)
  }
  
  
  # 5) save
  
  if(mouse){
    rho_i_file_name <- paste0('step_2_rho_i_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '_', i, '.rds')
  } else{
    rho_i_file_name <- paste0('step_2_rho_i_', adj_type, '_n_', n, '_nquery', cont_ind, '_', i, '.rds')
  }
  
  saveRDS(rho_i_result, file = file.path(temp_file_dir, rho_i_file_name))  
  
}

# rho_i without weights
estimate_intensities_stratum_parallel_with_yc_part1_v5 <- function(temp_file_dir, setting_info_list, i, mouse, X_truth) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: calculate rho_i for a single i and save it as
  #       
  #       temp_data/simu/step_2_v5_rho_i_block_banded_v2_n_100_nqueryk_i.rds'
  #
  # 
  # inputs:
  # 
  #
  # - temp_file_dir         (string)    'temp_data/simu'
  # - setting_info_list     (list)
  # - i                     (integer)   process_id from 1 to p
  # - mouse                 (boolean)   mouse (T) or simulation (F)
  # - X_truth               (boolean)   do we want to do estimation from true log-intensities?
  #
  #
  # outputs:
  #
  # - rho_i_result          (list)
  #
  #   - rho_i_est           (m x n matrix)  rho_i_est for each subject
  #   - rho_i_X_truth       (m x n matrix)  only computed if we have X_truth = T
  #
  # ----------------------------------------------------------------------------
  
  list2env(setting_info_list, envir = environment())
  
  # 1) load
  
  if(mouse){
    step_2_info_list <- paste0('part1_', ID, '_', discrete_level, '_t', time_scale, '.rds')
  } else{
    step_2_info_list <- paste0('part1_', adj_type, '_n_', n, '.rds')
  }
  
  
  # load `dataset$X_k_truth`, `weights`, `data_df4`
  results <- readRDS(file.path(temp_file_dir, step_2_info_list))
  list2env(results, envir = environment())
  
  n_time <- length(time_grid_est)  # 19
  
  
  # 2) filter process i and update subject ID's in case some don't have process i
  
  data_i <- data_df4[feature_id == feature_sel[i], ]
  kept_subjects <- unique(data_i$subject_num) %>% sort()  # in case any subjects do not have data for process i

  
  setDT(data_i)
  data_i[, subject_num := match(subject_num, sort(unique(subject_num)))]
  
  # 3) estimation
  if (nrow(data_i) == 0) { # no events, estimate is the zero intensity
    rho_mat <- matrix(0, nrow = n_time, ncol = n)
  } else {
    Gamma_i <- data_i[, estimate_density(time, time_grid_est), by = "subject_num"]
    rho_mat <- matrix(Gamma_i$rho_hat, nrow = n_time)  # m x n
  }
  
  # store rho_i
  
  rho_i_result <- list(rho_i_est = rho_mat) 
  
  # 4) X_truth if requested
  
  if(X_truth){
    # take (12 x 30 x 100), index only the i-th process, then take sample mean across n
    rho_i_result[['rho_i_X_truth']] <- exp(dataset$X_k_truth[i, , ])
  }
  
  
  # 5) save
  
  if(mouse){
    rho_i_file_name <- paste0('step_2_v5_rho_i_', ID, '_', discrete_level, '_t', time_scale, '_', i, '.rds')
  } else{
    rho_i_file_name <- paste0('step_2_v5_rho_i_', adj_type, '_n_', n, '_', i, '.rds')
  }
  
  saveRDS(rho_i_result, file = file.path(temp_file_dir, rho_i_file_name))  
  
}

# key_k
estimate_intensities_stratum_parallel_with_yc_part2 <- function(temp_file_dir, setting_info_list, cont_ind, k, mouse, X_truth) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: calculate rho_ij for a single i_j pair and save it as 
  #       
  #       temp_data/simu/step_2_rho_ij_block_banded_v2_n_100_nquery_cont_ind_k.rds
  #       temp_data/mice/step_2_rho_ij_Tau1_m0vr0_t5_nquery_cont_ind_k.rds
  # 
  # inputs:
  # 
  #
  # - temp_file_dir         (string)
  # - setting_info_list     (list)
  # - cont_ind              (integer)   which n_query index
  # - k                     (integer)   index number corresponding to a i_j pair
  # - mouse                 (boolean)   mouse (T) or simulation (F)
  # - X_truth               (boolean)   do we want X_truth?
  #
  #
  # outputs:
  #
  # - rho_ij_result       (list)
  # 
  #   - rho_ii_est        (m x m matrix)
  #   - rho_ii_X_truth    (m x m matrix)  only calcualted when X_truth = T
  #
  #
  # ----------------------------------------------------------------------------
  
  
  # 1) load
  
  
  list2env(setting_info_list, envir = environment())
  
  
  if(mouse){
    step_2_info_list <- paste0("part2_", ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    step_2_info_list <- paste0("part2_", adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }
  
  # load `dataset`
  results <- readRDS(file.path(temp_file_dir, step_2_info_list))
  list2env(results, envir = environment())
  
  
  key_ij <- keys[k]
  key_split <- strsplit(key_ij, "_")[[1]]
  
  i <- as.numeric(key_split[1])
  j <- as.numeric(key_split[2])
  
  # 2) prep
  
  t_seq <- time_grid_est
  
  p <- length(feature_sel) # 249
  n_time <- length(t_seq)  # 19
  
  
  
  
  # 2) filter process i and update subject ID's in case some don't have process i
  

  data_i <- data_df4[feature_id == feature_sel[i], ]
  data_j <- data_df4[feature_id == feature_sel[j], ]
  
  
  # fitting
  if (nrow(data_j) == 0 || nrow(data_i) == 0) { # if any entry has nothing, return the flat bivariate intensity
    rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
  } else {
    times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]  # dataframe where first col = subject, 2nd col = vector of observations 
    times_j <- data_j[, .(event_times_j = list(time)), by = subject_num]  # all of which are for process i
    times_ij <- merge(times_i, times_j, by = "subject_num", all = FALSE)   # now make it 3 columns: subject, process i, and process j
    
    kept_subjects <- unique(times_ij$subject_num) %>% sort()
    reweights = weights[kept_subjects] / sum(weights[kept_subjects])  # recalculate weights in case we discard subjects
    
    setDT(times_ij)
    times_ij[, subject_num := match(subject_num, sort(unique(subject_num)))]
    
    # 3) estimation
    
    if(nrow(times_ij) == 0){      # if they don't occur during the same replicates, return the flat bivariate intensity
      rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
    } else{
      Gamma_ij <- times_ij[, estimate_bivariate_density(
        event_times_i[[1]], event_times_j[[1]],
        t_seq, t_seq, 'i'), by = 'subject_num']
      
      bivariate_intensity <- matrix(Gamma_ij$V1, nrow = n_time^2)
      
      bivariate_intensity2 <- sweep(bivariate_intensity, 2, reweights, `*`) # multiply each 19-dim vec by `weight` which was already filtered for n <= n_large and normalized   

      rho_ij <- apply(bivariate_intensity2, 1, sum) # since these are normalized weights, just add them
      rho_ij_mat <- matrix(rho_ij, nrow = n_time)
    }
  }
  
  # store rho_ij_mat

  rho_ij_result <- list(rho_ii_est = rho_ij_mat)
  
  # 4) X_truth if requested
  if(X_truth){
    rho_ij_result[['rho_ii_X_truth']] <- estimate_rho_ij_from_Lambda(dataset$X_k_truth, i, j, reweights)
  }
  
  if(mouse){
    rho_ij_file_name <- paste0('step_2_rho_ij_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '_', k, '.rds')
  } else{
    rho_ij_file_name <- paste0("step_2_rho_ij_", adj_type, '_n_', n, '_nquery', cont_ind, '_', k, '.rds')
  }
  
  saveRDS(rho_ij_result, file = file.path(temp_file_dir, rho_ij_file_name))  
  
}

estimate_intensities_stratum_parallel_with_yc_part2_v5 <- function(temp_file_dir, setting_info_list, k, mouse, X_truth) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: calculate rho_ij for a single i_j pair and save it as 
  #       
  #       temp_data/simu/step_2_rho_ij_block_banded_v2_n_100_nquery_cont_ind_k.rds
  #       temp_data/mice/step_2_rho_ij_Tau1_m0vr0_t5_nquery_cont_ind_k.rds
  # 
  # inputs:
  # 
  #
  # - temp_file_dir         (string)
  # - setting_info_list     (list)
  # - cont_ind              (integer)   which n_query index
  # - k                     (integer)   index number corresponding to a i_j pair
  # - mouse                 (boolean)   mouse (T) or simulation (F)
  # - X_truth               (boolean)   do we want X_truth?
  #
  #
  # outputs:
  #
  # - rho_ij_result       (list)
  # 
  #   - rho_ii_est        (m^2 x n matrix)
  #   - rho_ii_X_truth    (m^2 x n matrix)  only calcualted when X_truth = T
  #
  #
  # ----------------------------------------------------------------------------
  
  
  # 1) load
  
  
  list2env(setting_info_list, envir = environment())
  
  
  if(mouse){
    step_2_info_list <- paste0("part1_", ID, '_', discrete_level, '_t', time_scale, '.rds')
  } else{
    step_2_info_list <- paste0("part1_", adj_type, '_n_', n, '.rds')
  }
  
  # load `dataset`
  results <- readRDS(file.path(temp_file_dir, step_2_info_list))
  list2env(results, envir = environment())
  
  
  key_ij <- keys[k]
  key_split <- strsplit(key_ij, "_")[[1]]
  
  i <- as.numeric(key_split[1])
  j <- as.numeric(key_split[2])
  
  # 2) prep
  
  t_seq <- time_grid_est
  
  p <- length(feature_sel) # 249
  n_time <- length(t_seq)  # 19
  
  
  
  
  # 2) filter process i and update subject ID's in case some don't have process i
  
  
  data_i <- data_df4[feature_id == feature_sel[i], ]
  data_j <- data_df4[feature_id == feature_sel[j], ]
  
  
  # fitting
  if (nrow(data_j) == 0 || nrow(data_i) == 0) { # if any entry has nothing, return the flat bivariate intensity
    rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
  } else {
    times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]  # dataframe where first col = subject, 2nd col = vector of observations 
    times_j <- data_j[, .(event_times_j = list(time)), by = subject_num]  # all of which are for process i
    times_ij <- merge(times_i, times_j, by = "subject_num", all = FALSE)   # now make it 3 columns: subject, process i, and process j
    
    kept_subjects <- unique(times_ij$subject_num) %>% sort()

    
    setDT(times_ij)
    times_ij[, subject_num := match(subject_num, sort(unique(subject_num)))]
    
    # 3) estimation
    
    if(nrow(times_ij) == 0){      # if they don't occur during the same replicates, return the flat bivariate intensity
      rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
    } else{
      Gamma_ij <- times_ij[, estimate_bivariate_density(
        event_times_i[[1]], event_times_j[[1]],
        t_seq, t_seq, 'i'), by = 'subject_num']
      
      bivariate_intensity <- matrix(Gamma_ij$V1, nrow = n_time^2) # m^2 x n
      

    }
  }
  
  # store rho_ij_mat
  
  rho_ij_result <- list(rho_ii_est = bivariate_intensity)
  
  # 4) X_truth if requested
  if(X_truth){
    rho_ij_result[['rho_ii_X_truth']] <- estimate_rho_ij_from_Lambda_v5(dataset$X_k_truth, i, j)
  }
  
  if(mouse){
    rho_ij_file_name <- paste0('step_2_v5_rho_ij_', ID, '_', discrete_level, '_t', time_scale, '_', k, '.rds')
  } else{
    rho_ij_file_name <- paste0("step_2_v5_rho_ij_", adj_type, '_n_', n, '_', k, '.rds')
  }
  
  saveRDS(rho_ij_result, file = file.path(temp_file_dir, rho_ij_file_name))  
  
}

estimate_intensities_stratum_parallel_with_yc_part3 <- function(temp_file_dir, setting_info_list, cont_ind, n_keys_univariate, n_keys_bivariate, mouse, X_truth) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: putting the rho_i, rho_ii, and rho_list results together
  #       used in script_step2_part3
  #
  #       output file:
  #
  #       step_2_rho_list_block_banded_c0_n_100_nquery1.rds
  #
  # 
  # inputs
  #
  # - temp_file_dir
  # - setting_info_list
  # - cont_ind              (integer)   n_query id
  # - n_keys_univariate     (integer)   p = 12
  # - n_keys_bivariate      (integer)   pc2 + p = 78
  # - mouse                 (boolean)   is it a mouse?
  #
  # 
  # outputs:
  #
  # - result  --> 'step_2_rho_list_block_banded_v2_n_1000_nquery1.rds'
  # 
  #   - step_2    (each item is a p x m matrix)  rho_i_est
  #   - step_2b   (each item is a i_j list of matrices)   rho_ii_est
  # 
  # ----------------------------------------------------------------------------
  
  
  
  # 1) retrieve data
  
  list2env(setting_info_list, envir = environment())
  
  if(mouse){
    # keys_univariate = vector of c(1, 2, 3, ..., p)
    # keys_bivariate = vector of c('1_1', '1_2', ..., 'p_p')
    rho_i_file_names  <- paste0(temp_file_dir, '/step_2_rho_i_',  ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '_', 1:n_keys_univariate, '.rds')
    rho_ij_file_names <- paste0(temp_file_dir, '/step_2_rho_ij_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '_', 1:n_keys_bivariate, '.rds') 
    part2_file_name   <- paste0('part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    rho_i_file_names  <- paste0(temp_file_dir, "/step_2_rho_i_",  adj_type, '_n_', n, '_nquery', cont_ind, '_', 1:n_keys_univariate, '.rds')
    rho_ij_file_names <- paste0(temp_file_dir, "/step_2_rho_ij_", adj_type, '_n_', n, '_nquery', cont_ind, '_', 1:n_keys_bivariate, '.rds')
    part2_file_name   <- paste0('part2_', adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }

  # load `keys` 
  results <- readRDS(file.path(temp_file_dir, part2_file_name))
  keys <- results$keys
  #list2env(results, envir = environment())
  
  # 2) Load all rho_i files into a list + reorganize them
  
  rho_i_list_raw <- lapply(rho_i_file_names, readRDS)  # list of p items --> `rho_i_est` `rho_i_X_truth` etc...
                                                       # we wish to create step_2 --> `rho_i_est` = pxm matrix, `rho_i_X_truth` = pxm matrix etc...
  

  vec_names <- names(rho_i_list_raw[[1]])
  step_2 <- lapply(vec_names, function(nm) {
    mat <- t(sapply(rho_i_list_raw, function(sub) sub[[nm]]))  # stack the vectors
    return(mat)   
  })
  
  names(step_2) <- vec_names
  
  
  
  # 3) Load all rho_ij files into a list + organize them
  
  
  rho_ij_list_raw <- lapply(rho_ij_file_names, readRDS) # list of p items --> `rho_ii_est`, `rho_ii_X_truth` etc...
                                                        # we wish to create step_2b --> `rho_ii_est` = list of mxm matrices, `rho_ii_X_truth` = list of mxm matrices
  
  names(rho_ij_list_raw) <- keys
  mat_names <- names(rho_ij_list_raw[[1]])
  
  step_2b <- list()
  
  for (nm in mat_names) {
    
    # for this result name, extract the mxm matrix from each of the p lists
    mats <- lapply(rho_ij_list_raw, function(x) x[[nm]])
    
    # store as a list of mxm matrices
    step_2b[[nm]] <- mats
  }
  

  

  
  
  # 4) store them

  
  result <- list(step_2 = step_2,
                 step_2b = step_2b)
  
  if(mouse){
    rho_list_name <- paste0('step_2_rho_list_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    rho_list_name <- paste0('step_2_rho_list_', adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }
  
  # delete files
  
  # file.remove(rho_i_file_names)
  # file.remove(rho_ij_file_names)
  
  saveRDS(result, file = file.path(temp_file_dir, rho_list_name))  
  
}

full_conditional_estimation_with_no_truth_part2b <- function(temp_file_dir, setting_info_list, cont_ind, mouse, X_truth){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimation of everything else after step 2 in series
  #       used in script_fit_mice_data_part2b
  #
  #       save file called:
  #
  #       temp_data/simu/part3_block_banded_c0_n_100_nquery1.rds
  # 
  # 
  # inputs
  #
  # - temp_file_dir
  # - setting_info_list
  # - cont_ind              (integer)   n_query id
  # - mouse                 (boolean)   is it a mouse?
  #
  # 
  # loading
  #
  # - part2             (list of various parameters)
  # - step_2_rho_list   (list of step_2 and step_2b)
  # 
  # outputs:
  #
  # - estimated_graphs  (list of steps 2 and later, each of these differs based on y_c_query)
  # 
  # estimated_graphs <- list(step_2 = step_2, step_2b = step_2b, step_3 = step_3,
  #                          step_4 = step_4, step_5 = step_5, step_5b = step_5b, step_9 = step_9, step_9b = step_9b,
  #                          step_10 = step_10, step_11 = step_11)
  # 
  # ----------------------------------------------------------------------------
  

  
  # 0) load 
  
  print(paste0('starting part 2b, cont_ind = ', cont_ind))
  
  list2env(setting_info_list, envir = environment())
  
  # load everything from step_1
  if(mouse){
    step_2_info_list <- paste0("part2_", ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
    rho_list_name <- paste0('step_2_rho_list_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
    datafile_error_name <- paste0("dataset_part2_", ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.RData')  # in case we need to quit and troubleshoot
  } else{
    step_2_info_list <- paste0("part2_", adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
    rho_list_name <- paste0("step_2_rho_list_", adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
    datafile_error_name <- paste0("dataset_part2_", adj_type, '_n_', n, '_nquery', cont_ind, '.RData')  # in case we need to quit and troubleshoot
  }
  
  # load `W_y`, `p`
  results <- readRDS(file.path(temp_file_dir, step_2_info_list))
  p <- results$p
  W_y <- results$W_y
  adj_mat_i <- results$adj_mat_i
  #list2env(results, envir = environment())
  steps_2_and_2b <- readRDS(file.path(temp_file_dir, rho_list_name))
  
  step_2 <- steps_2_and_2b$step_2
  step_2b <- steps_2_and_2b$step_2b
  
  # ------------------
  # Step 2b onward
  # ------------------
  
  i_neq_j <- T

  step_3  <- step_3_g_ij(step_2, step_2b, i_neq_j)
  
  
  
  step_4 <- tryCatch({
    #step_4_eigendecomp(step_3, p, same_basis, constant_d)
    step_4_eigendecomp_troubleshoot(step_3, p)
  }, error = function(e) {
    cat("Error in step_4, saving dataset...\n")
    save(dataset, file = file.path(temp_file_dir, datafile_error_name))
    stop(e)
  })
  
  # 3) steps 5: obtain KL_cor and KL_prec
  

  #step_5 <- step_5_KL_covariance(step_3, step_4)
  step_5 <- step_5_KL_covariance_eigencases(step_3, step_4)
  step_5b <- step_5b_KL_correlation(step_5, p)
  step_5c <- step_5c_KL_precision(step_5b, p)
  
  # estimate C_HS, w_mat from KL without thresh 
  
  step_11_KL_no_thresh  <- step_11_HS_norms_from_KL(step_5c, p)
  step_11b_KL_no_thresh <- step_11b_HS_norms_from_KL(step_5b, p)
  
  # estimate C_HS, w_mat from KL WITH thresh
  
  step_11_GIC_bundle <- steps_10_11_GIC_from_KL(step_5b, p, W_y)
  
  step_11_KL_yes_thresh  <- step_11_HS_norms_from_KL_GIC(step_11_GIC_bundle)
  step_11b_KL_yes_thresh <- step_11b_HS_norms_from_KL_GIC(step_11_GIC_bundle)
  
  # 4) step 9: Construct mxm object from KL correlation + get precision operator
  
  step_9 <- step_9_C_cond_from_KL_cor(step_4, step_5b)
  step_9b <- step_9b_eigenfunction_outers(step_4)

  
  block <- F
  MP <- F
  step_10 <- tryCatch({
    step_10_P_cond(step_9, p, block, MP)
  }, error = function(e) {
    cat("Error occurred in step_10, saving dataset...\n")
    save(dataset, file = file.path(temp_file_dir, datafile_error_name))
    cat("Dataset saved to dataset.RData\n")
    stop(e)  # Re-throw the error
  })
  

  # 5) estimate C_HS, w_mat from pxp without thresh
  
  step_11_mxm_no_thresh_bundle <- step_11_HS_norms(step_9, step_10, p)
  
  # 6) estiamte C_HS, w_mat from pxp WITH thresh

  # step_11_mxm_GIC_bundle <- steps_10_11_GIC(step_9, p, W_y)
  # 
  # step_11_KL_yes_thresh  <- step_11_HS_norms_from_KL_GIC(step_11_mxm_GIC_bundle)
  # step_11b_KL_yes_thresh <- step_11b_HS_norms_from_KL_GIC(step_11_mxm_GIC_bundle)

  
  # 7) collect all of 11 and 11b results
  
  step_11 <- c(step_11_KL_no_thresh,
               step_11_KL_yes_thresh,
               step_11_mxm_no_thresh_bundle$step_11)
  
  step_11b <- c(step_11b_KL_no_thresh,
                step_11b_KL_yes_thresh,
                step_11_mxm_no_thresh_bundle$step_11b)
  
  
  if(! mouse){
    step_12 <- step_12_ROC(step_11, adj_mat_i)
  }
  
  step_12b <- step_12b_adj_mat(step_11)
  
  estimated_graphs <- list(step_2 = step_2, step_2b = step_2b, step_3 = step_3,
                           step_4 = step_4, step_5 = step_5, step_5b = step_5b, step_5c = step_5c, 
                           step_9 = step_9, step_9b = step_9b, step_10 = step_10, 
                           step_11 = step_11, step_11b = step_11b, 
                           step_12 = step_12, step_12b = step_12b)
  
  

  if(mouse){
    file_name <- paste0('part3_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    estimated_graphs[['step_12']] <- step_12
    file_name <- paste0("part3_", adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }
  
  saveRDS(estimated_graphs, file = file.path(temp_file_dir, file_name))  
  
}

full_conditional_estimation_with_no_truth_part3 <- function(temp_file_dir, setting_info_list, cont_inds, mouse){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: merge all the estimates from each y_c_query
  #
  #       just return an object called 'estimated_graphs'
  #
  # inputs:
  # 
  # - temp_file_dir
  # - setting_info_list
  # - cont_inds           (scalar)  number of y_c_queries
  # - mouse               (boolean) are we working with mice data
  #
  # outputs:
  #
  # - estimated_graphs   (list of the following)
  #
  #   - step_0_events     (list of i_j entries --> each entry is a vector of timestamps)
  #   - step_1            X_k_suffix
  #   - step_1b           Lambda_k_suffix
  #   - step_1c           mu_t_suffix
  #   - step_2            (list for each y_c_query --> rho_i_est)
  #   - step_2b           (list for each y_c_query --> rho_ii_est)
  #   - step_3            (list for each y_c_query --> g_ij_est)
  #   - step_4            (list for each y_c_query --> eigen_decomp_est)
  #   - step_5            (list for each y_c_query --> KL_cov_est)
  #   - step_5b           (list for each y_c_query --> KL_cor_est)
  #   - step_5c           (list for each y_c_query --> KL_prec_est)
  #   - step_9            (list for each y_c_query --> C_cond_est_full, C_cond_est_unnorm_full)
  #   - step_9b
  #   - step_10           (list for each y_c_query --> P_cond_est_full, P_cond_est_unnorm_full)
  #   - step_11           (list for each y_c_query --> w_mat_est, w_mat_est_unnorm, C_HS_est, C_HS_est_unnorm)
  #   - step_12           (optional, not used for mice)
  #   - y_c_query
  #   - p
  #   - time_grid
  #   - time_grid_est
  #   - time_grid_both
  #
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # read from steps 1 and 2 and then remove everything
  # ----------------------------------------------------------------------------  
  
  # Vector of file names from part3
  
  list2env(setting_info_list, envir = environment())
  
  if(mouse){
    file_names <- paste0(temp_file_dir, '/part3_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
  } else{
    file_names <- paste0(temp_file_dir, '/part3_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '.rds')
  }
  
  
  
  # Step 1: Load all files into a list
  all_loaded <- lapply(file_names, readRDS)
  
  # Step 2: Get all step names (assumes all files have same names)
  step_names <- names(all_loaded[[1]])
  
  # Step 3: Reorganize by step
  estimated_graphs <- setNames(lapply(step_names, function(step) {
    lapply(all_loaded, `[[`, step)  # collect that step from all files
  }), step_names)
  

 
  
  # load step 1 and keep important items in `estimated_graphs`
  if(mouse){
    step_1_list_name <- paste0('part1_', ID, '_', discrete_level, '_t', time_scale, '.rds')
  } else{
    step_1_list_name <- paste0('part1_', adj_type, '_n_', n, '.rds')
  }
  
  results <- readRDS(file.path(temp_file_dir, step_1_list_name))
  list2env(results, envir = environment())
  
  # label the sublists with their respective query_y_cs
  
  estimated_graphs <- lapply(estimated_graphs, function(x) {
    names(x) <- round(query_y_cs[,1], 3)
    x
  })

  estimated_graphs[['step_0_events']] <- step_0_events
  estimated_graphs[['step_1']] <- step_1
  estimated_graphs[['step_1b']] <- step_1b
  estimated_graphs$y_c_query <- query_y_cs
  estimated_graphs$p <- p
  estimated_graphs$Y_continuous <- y_c_strata
  
  # load step 2 and keep weights in step_2
  
  if(mouse){
    step_2_list_names <- paste0(temp_file_dir, '/part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
  } else{
    step_2_list_names <- paste0(temp_file_dir, '/part2_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '.rds')
  }
  
  step_2_all_data <- lapply(step_2_list_names, readRDS)
  
  all_weights <- lapply(step_2_all_data, `[[`, "weights")
  names(all_weights) <- round(query_y_cs[,1], 3)
  all_W_y <- lapply(step_2_all_data, `[[`, "W_y")
  names(all_W_y) <- round(query_y_cs[,1], 3)
  
  estimated_graphs[['weights']] <- all_weights
  estimated_graphs[['W_y']] <- all_W_y
  
  if(mouse){
    estimated_graphs$time_grid_est <- time_grid_est
  } else{
    estimated_graphs$time_grid <- time_grid
    estimated_graphs$time_grid_est <- time_grid_est
    estimated_graphs$time_grid_both <- time_grid_both    
  }

  
  # ------------------------
  # REMOVE ALL FILES HERE 
  # ------------------------
  
  if(mouse){
    part1_file_name          <- paste0('part1_', ID, '_', discrete_level, '_t', time_scale, '.rds')
    part2_file_name          <- paste0('part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
    part2_rho_i_file_name    <- paste0('step_2_rho_i_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '_', 1:p, '.rds')
    part2_rho_ij_file_name   <- paste0('step_2_rho_ij_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '_', 1:length(keys), '.rds')
    part2_rho_list_file_name <- paste0('step_2_rho_list_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
    part3_file_name          <- paste0('part3_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
  } else{
    part1_file_name          <- paste0('part1_', adj_type, '_n_', n, '.rds')
    part2_file_name          <- paste0('part2_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '.rds')
    part2_rho_i_file_name    <- paste0('step_2_rho_i_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_', 1:p, '.rds')
    part2_rho_ij_file_name   <- paste0('step_2_rho_ij_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_', 1:length(keys), '.rds')
    part2_rho_list_file_name <- paste0('step_2_rho_list_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '.rds')
    part3_file_name          <- paste0('part3_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '.rds')
  }
  

  # file.remove(file.path(temp_file_dir, part1_file_name)) 
  # file.remove(file.path(temp_file_dir, part2_file_name)) 
  # file.remove(file.path(temp_file_dir, part2_rho_i_file_name)) 
  # file.remove(file.path(temp_file_dir, part2_rho_ij_file_name)) 
  # file.remove(file.path(temp_file_dir, part2_rho_list_file_name)) 
  # file.remove(file.path(temp_file_dir, part3_file_name)) 
  
  
  return(estimated_graphs)     
}

full_conditional_estimation_CPGM <- function(processed_data, time_grid_est, method, terse, ncores, dir){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: troubleshoot the full conditional estimation procedure with ground truths between steps
  #
  # - instead of previously, where we only used G_ii and discarded G_ij
  # - here we actively look to estimate G_ij and use it
  # - 9/3/2025 - use G_ij^k, for each subject
  #
  # - CPGM method
  #
  # Input:
  #
  # - processed_data
  # - time_grid_est (m-dim vec) time discretization
  # - method   ('CPGM')
  # - terse    (boolean) if true, return much less
  # - ncores
  # - dir      (string) folder name, such as "simu_results_banded_trig2_2"
  # 
  #
  # ----------------------------------------------------------------------------
  
  # check for errors 
  if (!(method %in% c("CPGM"))) {
    stop("Error: method must be CPGM")
  }    
  
  
  arr_mat <- matrix(1:6, nrow = 2, byrow = F)
  arr_mat_8 <- matrix(1:8, nrow = 2, byrow = F)
  arr_mat_10 <- matrix(1:10, nrow = 2, byrow = F)  
  
  # part 0 - check ground truth precisions and process data --------------------
  
  # process data
  processed_data <- step_0_preprocess(dataset)
  data_df4 = processed_data[[1]]
  y_c_strata = processed_data[[2]]
  query_y_cs = processed_data[[3]]
  patient_sel = processed_data[[4]]
  feature_sel = processed_data[[5]]
  rm(processed_data)

  # load parameters
  
  n <- dim(y_c_strata)[1]
  m_est <- length(time_grid_est) 
    
  
  i_neq_j <- T
  rho_kernel <- T
  
  # part 1 - log intensities ---------------------------------------------------
  
  step_1 <- step_1_log_intensities(dataset, data_df4, time_grid_est, time_grid, time_grid_both)
  
  
  # now, we regress on y_c  ----------------------------------------------------
  
  print('at covariate loop')
  cont_inds <- 1:nrow(query_y_cs)
  estimated_graphs_v2 <- lapply(cont_inds, function(cont_ind) {  
    
    
    print(paste0(cont_ind, ' out of ', nrow(query_y_cs)))
    
    # obtain this query's parameters
    kernel_params_i = dataset$true_graphs[[cont_ind]]$P_block_kronecker
    adj_mat_i <- dataset$true_graphs[[cont_ind]]$adj_mat
    query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()    
    
    # part 2 - intensities -----------------------------------------------------
    
    step_2 <- step_2_rho_i(dataset, data_df4, kernel_params_i, rho_kernel, patient_sel, feature_sel,
                           time_grid, time_grid_est, query_y_c, y_c_strata, ncores)
    
    
    # part 2.1 - bivariate intensities -----------------------------------------
    
    step_2b <- step_2_rho_ij(step_1, step_2, kernel_params_i, i_neq_j)
    
    
    
    # part 3 - GP covariance estimation ----------------------------------------
    
    step_3 <- step_3_g_ij(step_2, step_2b, kernel_params_i, i_neq_j)
    
    
    # part 4 - eigendecomposition of GP covariance -----------------------------
    
    # careful about some processes being empty
    
    # step_4 <- step_4_eigendecomp(step_3, p, time_grid, time_grid_est)
    result <- tryCatch({
      step_4 <- step_4_eigendecomp(step_3, p, time_grid, time_grid_est)
    }, error = function(e) {
      cat("Error occurred in step_4, saving dataset...\n")
      save(dataset, file = paste0(dir, "/dataset.RData"))
      cat("Dataset saved to dataset.RData\n")
      stop(e)  # Re-throw the error
    })    
    
    
    # part 5 - covariance of KL coefficients -----------------------------------
    
    step_5 <- step_5_KL_covariance(step_3, step_4)
    
    # part 9 - conditional_correlation from KL covariance ----------------------
    
    step_9 <- step_9_C_cond_from_KL_cov(step_4, step_5, kernel_params_i)
    
    # part 10 - precision operator ---------------------------------------------
    
    block <- F
    MP <- F
    
    # step_10 <- step_10_P_cond(step_9, kernel_params_i, p, block, MP)
    
    result <- tryCatch({
      step_10 <- step_10_P_cond(step_9, kernel_params_i, p, block, MP)
    }, error = function(e) {
      cat("Error occurred in step_10, saving dataset...\n")
      save(dataset, file = paste0(dir, "/dataset.RData"))
      cat("Dataset saved to dataset.RData\n")
      stop(e)  # Re-throw the error
    })       
    
    # part 11 - HS norms of precision matrix -----------------------------------
    
    step_11 <- step_11_HS_norms(step_10, adj_mat_i, p)
    
    # part 12 - ROC curve ------------------------------------------------------
    
    step_12 <- step_12_ROC(step_11, adj_mat_i)
    
    # part 13 - metrics
    
    metrics <- get_metrics(list(step_10$P_cond_est_full,   # P_cond_est_full
                                step_9$C_cond_est_full,    # C_cond_est_full
                                NULL,
                                step_10$P_cond_coarse_ground_truth_full,   # P_cond_coarse_ground_truth_full
                                step_9$C_cond_corase_ground_truth_full,    # C_cond_coarse_ground_truth_full
                                NULL,
                                step_11$w_mat_est,   # w_mat_est
                                adj_mat_i,
                                step_12$roc_est,   # roc_est
                                step_2$rho_i_est,    # rho_i_est
                                step_2$rho_i_coarse_truth,    # rho_i_coarse_truth
                                step_2b$rho_ii_est,   # rho_ii_est
                                step_2b$rho_ii_coarse_truth,   # rho_ii_coarse_truth
                                step_3$g_ij_est,    # g_ij_est
                                step_3$g_ij_coarse_truth))   # g_ij_coarse_truth     
    
    
    # 8) save graphs (recall that we have ground truths of the form a_b_c)
    
    list(step_2 = step_2,
         step_2b = step_2b,
         step_3 = step_3,
         step_4 = step_4,
         step_5 = step_5,
         step_9 = step_9,
         step_10 = step_10,
         step_11 = step_11,
         step_12 = step_12,
         metrics = metrics,
         kernel_params_i = kernel_params_i)
    
    
    
  })# , mc.cores = ncores) # done with all y_c levels
  
  
  
  names(estimated_graphs_v2) <- apply(query_y_cs, 1, function(x){paste(x, collapse = '_')})
  
  
  # what to output
  
  
  estimated_graphs_part_1 <- list(g_01 = g_01, step_1 = step_1)
  
  return(list(estimated_graphs_part_1 = estimated_graphs_part_1,      
              estimated_graphs_part_2 = estimated_graphs_v2))      
  
  
  
}

