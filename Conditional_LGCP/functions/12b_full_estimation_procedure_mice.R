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

full_conditional_estimation_with_no_truth_part1 <- function(dataset, setting_info_list, ncores, temp_file_dir, mouse){
  

  list2env(setting_info_list, envir = environment())  

  
  
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
  y_c_strata   <- processed_data[[2]]  # full 
  query_y_cs   <- processed_data[[3]]
  patient_sel  <- processed_data[[4]]
  feature_sel  <- processed_data[[5]]
  y_c_strata_sel <- processed_data$y_c_strata_sel
  
  rm(processed_data)
  
  p <- length(feature_sel)
  full <- F

  # ----------------------------------------------------------------------------
  # Step 1 - Log intensities
  # ----------------------------------------------------------------------------
  
  step_1 <- step_1_log_intensities(dataset, data_df4, time_grid_est, NA, NA, full) 
  
  # more things to store:
  # - weights2
  # - keys
  
  print('====TROUBLESHOOT 5====')
  print(paste0('Class of y_c_strata_sel should be matrix array: ', class(y_c_strata_sel)))
  print(paste0('Dimension of y_c_strata_sel should be n x 1: ', dim(y_c_strata_sel)))
  print('======================')
  
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata_sel)
  
  # get the keys 
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
    step_1 = step_1,
    data_df4 = data_df4,
    y_c_strata = y_c_strata,
    y_c_strata_sel = y_c_strata_sel,
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
  }
  
  if(mouse){
    datafile_name <- paste0('part1_', ID, '_', discrete_level, '_t', time_scale, '.rds')
  } else{
    datafile_name <- paste0('part1_', adj_type, '_n_', n, '.rds')
  }
  

  
  saveRDS(results, file = file.path(temp_file_dir, datafile_name))  
  
  print('======TROUBLESHOOT 4======')
  print('Currently in part 1')
  print(paste0('Did we save ', datafile_name, '? ', file.exists(file.path(temp_file_dir, datafile_name))))
  print('==========================')
  
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


estimate_intensities_stratum_parallel_with_yc_part0 <- function(temp_file_dir, setting_info_list, cont_ind, mouse){
  
  # for this strata, calculate weights 
  
  list2env(setting_info_list, envir = environment())
  
  

  
  if(mouse){
    step_1_info_list <- paste0('part1_', ID, '_', discrete_level, '_t', time_scale, '.rds')
  } else{
    step_1_info_list <- paste0('part1_', adj_type, '_n_', n, '.rds')
  }
  
  print('======TROUBLESHOOT 4b=====')
  print('Currently in step_2 part 0')
  print(paste0('Did we save ', step_1_info_list, '? ', file.exists(file.path(temp_file_dir, step_1_info_list))))
  print('==========================')
  
  results <- readRDS(file.path(temp_file_dir, step_1_info_list))
  list2env(results, envir = environment())
  
  # get the query and get the weights
  y_c_query <- query_y_cs[cont_ind,]
  
  weights <- apply(y_c_strata, 1, function(row) {
    step_6_kernel(as.numeric(row), y_c_query, gamma_c) 
  })    
  weights2 <- weights / sum(weights) # normalize
  
  results[['weights2']] <- weights2
  
  # save
  
  if(mouse){
    datafile_name <- paste0('part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    datafile_name <- paste0('part2_', adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }
  
  saveRDS(results, file = file.path(temp_file_dir, datafile_name))  
  

}

# rho_i
estimate_intensities_stratum_parallel_with_yc_part1 <- function(temp_file_dir, setting_info_list, cont_ind, i, mouse) {
  
  # we first calculate rho_i_list
  # cont_ind = which continuous covariate
  # i = process_id from 1 to p
  
  list2env(setting_info_list, envir = environment())
  
  
  if(mouse){
    step_2_info_list <- paste0('part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    step_2_info_list <- paste0('part2_', adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }

  
  results <- readRDS(file.path(temp_file_dir, step_2_info_list))
  list2env(results, envir = environment())
  
  
  n_time <- length(time_grid_est)  # 19
  
  # step 2: univariate case 
  
  data_i <- data_df4[feature_id == feature_sel[i], ]
  
  if (nrow(data_i) == 0) { # no events, estimate is the zero intensity
    rho_i <- rep(0, n_time)
  } else {
    Gamma_i <- data_i[, estimate_density(time, time_grid_est), by = "subject_num"]
    rho_mat <- matrix(Gamma_i$rho_hat, nrow = n_time)
    
    included_weights <- weights2[unique(Gamma_i$subject_num)]  # assume weights2 contains everyone so we have to filter here
    
    normalized_weights = included_weights / sum(included_weights)
    
    print('======TROUBLESHOOT 3=========')
    print(paste0('The current sum of weights is: ', sum(included_weights)))
    print(paste0('After normalizing, the sum of the weights is: ', sum(normalized_weights)))
    
    rho_mat2 <- sweep(rho_mat, 2, normalized_weights, `*`) # multiply each 19-dim vec by its normalized weight
    
    rho_i <- apply(rho_mat2, 1, sum) # since these are normalized weights, just add them
  }
  
  # store rho_i
  
  out_list <- list()
  out_list[[i]] <- rho_i
  
  
  if(mouse){
    rho_i_file_name <- paste0('step_2_rho_i_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '_', i, '.rds')
  } else{
    rho_i_file_name <- paste0('step_2_rho_i_', adj_type, '_n_', n, '_nquery', cont_ind, '_', i, '.rds')
  }
  
  saveRDS(out_list, file = file.path(temp_file_dir, rho_i_file_name))  
  
}
# key_k
estimate_intensities_stratum_parallel_with_yc_part2 <- function(temp_file_dir, setting_info_list, cont_ind, k, mouse) {
  
  
  # bivariate estimation
  
  # k = index of the `key` vector = c('1_1', '1_2', ...)
  
  # 1) retrieve data
  
  
  list2env(setting_info_list, envir = environment())
  
  
  if(mouse){
    step_2_info_list <- paste0("part2_", ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    step_2_info_list <- paste0("part2_", adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }
  
  
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
  
  
  
  
  # step 3: bivariate case - for just the k-th case
  
  
  
  
  data_i <- data_df4[feature_id == feature_sel[i], ]
  data_j <- data_df4[feature_id == feature_sel[j], ]
  
  
  # fitting
  if (nrow(data_j) == 0 || nrow(data_i) == 0) { # if any entry has nothing, return the flat bivariate intensity
    rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
  } else {
    times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]  # dataframe where first col = subject, 2nd col = vector of observations 
    times_j <- data_j[, .(event_times_j = list(time)), by = subject_num]  # all of which are for process i
    times_ij <- merge(times_i, times_j, by = "subject_num", all = FALSE)   # now make it 3 columns: subject, process i, and process j
    
    if(nrow(times_ij) == 0){      # if they don't occur during the same replicates, return the flat bivariate intensity
      rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
    } else{
      Gamma_ij <- times_ij[, estimate_bivariate_density(
        event_times_i[[1]], event_times_j[[1]],
        t_seq, t_seq, 'i'), by = 'subject_num']
      
      bivariate_intensity <- matrix(Gamma_ij$V1, nrow = n_time^2)
      
      included_weights <- weights2[unique(Gamma_ij$subject_num)]
      
      bivariate_intensity2 <- sweep(bivariate_intensity, 2, included_weights, `*`) # multiply each 19-dim vec by its normalized weight      
      
      rho_ij <- apply(bivariate_intensity2, 1, sum) # since these are normalized weights, just add them
      rho_ij_mat <- matrix(rho_ij, nrow = n_time)
    }
  }
  
  # store rho_ij_mat
  
  out_list <- list()
  out_list[[key_ij]] <- rho_ij_mat
  
  if(mouse){
    rho_ij_file_name <- paste0('step_2_rho_ij_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '_', k, '.rds')
  } else{
    rho_ij_file_name <- paste0("step_2_rho_ij_", adj_type, '_n_', n, '_nquery', cont_ind, '_', k, '.rds')
  }
  
  saveRDS(out_list, file = file.path(temp_file_dir, rho_ij_file_name))  
  
}

estimate_intensities_stratum_parallel_with_yc_part3 <- function(temp_file_dir, setting_info_list, cont_ind, n_keys_univariate, n_keys_bivariate, mouse) {
  # used in script_step2_part3
  
  # putting the results together
  
  
  
  # 1) retrieve data
  
  list2env(setting_info_list, envir = environment())
  
  if(mouse){
    # keys_univariate = vector of c(1, 2, 3, ..., p)
    # keys_bivariate = vector of c('1_1', '1_2', ..., 'p_p')
    rho_i_file_names  <- paste0(temp_file_dir, '/step_2_rho_i_',  ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '_', 1:n_keys_univariate, '.rds')
    rho_ij_file_names <- paste0(temp_file_dir, '/step_2_rho_ij_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '_', 1:n_keys_bivariate, '.rds')    
  } else{
    rho_i_file_names  <- paste0(temp_file_dir, "/step_2_rho_i_",  adj_type, '_n_', n, '_nquery', cont_ind, '_', 1:n_keys_univariate, '.rds')
    rho_ij_file_names <- paste0(temp_file_dir, "/step_2_rho_ij_", adj_type, '_n_', n, '_nquery', cont_ind, '_', 1:n_keys_bivariate, '.rds')
  }

  
  # Step 1: Load all files into a list - but there are some serious wrangling issues
  # - rho_i_list
  # - rho_ij_list

  rho_i_list_raw <- lapply(rho_i_file_names, readRDS) # it's a list of lists, only rho_i_list[[2]][[2]] is relevant
  
  rho_i_list <- mapply(function(x, idx) x[[idx]], 
                       rho_i_list_raw, 
                       seq_along(rho_i_list_raw),
                       SIMPLIFY = FALSE)
  names(rho_i_list) <- 1:length(rho_i_list)
  
  
  rho_ij_list_raw <- lapply(rho_ij_file_names, readRDS)
  rho_ij_list <- lapply(rho_ij_list_raw, `[[`, 1)
  names(rho_ij_list) <- sapply(rho_ij_list_raw, function(x) names(x)[1])
  
  # delete files
  
  file.remove(rho_i_file_names)
  file.remove(rho_ij_file_names)
  
  
  
  # store rho_list and rho_i_est
  
  rho_i_est <- do.call(rbind, rho_i_list)
  
  rho_list = list(rho_i_list, rho_ij_list)
  
  step_2 <- list(rho_i_est = rho_i_est,
                 rho_list = rho_list)
  
  if(mouse){
    rho_list_name <- paste0('step_2_rho_list_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    rho_list_name <- paste0("step_2_rho_list_", adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }
  
  saveRDS(step_2, file = file.path(temp_file_dir, rho_list_name))  
  
  print("Checking existence (after saving) of:")
  print(rho_list_name)
  print(file.exists(paste0(temp_file_dir, '/', rho_list_name)))
  
  
}

full_conditional_estimation_with_no_truth_part2b <- function(temp_file_dir, setting_info_list, cont_ind, mouse){
  # used in script_fit_mice_data_part2b
  
  # all the estimation after step_2
  
  # load 
  
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
  
  print('====TROUBLESHOOT 2====')
  print(paste0('currently in full_conditional_estimation_with_no_truth_part2b with cont_ind = ', cont_ind))
  print(paste0('named dataset: ', step_2_info_list))
  print(paste0('is the part2_ dataset present?: ', file.exists(file.path(temp_file_dir, step_2_info_list))))
  print(paste0('named dataset: ', rho_list_name))
  print(paste0('is the step_2_rho_list_ dataset present?: ', file.exists(file.path(temp_file_dir, rho_list_name))))
  print('======================')
  
  results <- readRDS(file.path(temp_file_dir, step_2_info_list))
  list2env(results, envir = environment())
  step_2 <- readRDS(file.path(temp_file_dir, rho_list_name))
  
  # ------------------
  # Step 2b onward
  # ------------------
  full <- F
  i_neq_j <- T
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
  
  

  if(mouse){
    file_name <- paste0('part3_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  } else{
    file_name <- paste0("part3_", adj_type, '_n_', n, '_nquery', cont_ind, '.rds')
  }
  
  saveRDS(estimated_graphs, file = file.path(temp_file_dir, file_name))  
  
}

full_conditional_estimation_with_no_truth_part3 <- function(temp_file_dir, setting_info_list, cont_inds, mouse){
  
  #
  # cont_inds = number
  # 
  
  # ----------------------------------------------------------------------------
  # read from steps 1 and 2 and then remove everything
  # ----------------------------------------------------------------------------  
  
  # Vector of file names
  
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
  

 
  
  # load step 1
  if(mouse){
    step_1_list_name <- paste0('part1_', ID, '_', discrete_level, '_t', time_scale, '.rds')
  } else{
    step_1_list_name <- paste0('part1_', adj_type, '_n_', n, '.rds')
  }
  
  results <- readRDS(file.path(temp_file_dir, step_1_list_name))
  # load all variables
  list2env(results, envir = .GlobalEnv)
  
  
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
  

  file.remove(file.path(temp_file_dir, part1_file_name)) 
  file.remove(file.path(temp_file_dir, part2_file_name)) 
  file.remove(file.path(temp_file_dir, part2_rho_i_file_name)) 
  file.remove(file.path(temp_file_dir, part2_rho_ij_file_name)) 
  file.remove(file.path(temp_file_dir, part2_rho_list_file_name)) 
  file.remove(file.path(temp_file_dir, part3_file_name)) 
  
  
  return(all_results)     
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

