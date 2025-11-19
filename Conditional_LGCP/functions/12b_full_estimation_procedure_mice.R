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

full_conditional_estimation_with_no_truth_part1 <- function(dataset, setting_info_list, method, ncores, temp_file_dir, mouse = F){
  
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
  
  # Save results for Stage 2
  
  results <- list(
    dataset = dataset,
    step_0_events = step_0_events,
    step_1 = step_1,
    data_df4 = data_df4,
    y_c_strata = y_c_strata,
    query_y_cs = query_y_cs,
    patient_sel = patient_sel,
    feature_sel = feature_sel,
    p = p,
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
  
  ID <- setting_info_list[['ID']]
  y_c_structure <- setting_info_list[['y_c_structure']]
  time_scale <- setting_info_list[['time_scale']]
  method <- setting_info_list[['method']]
  discrete_level <- setting_info_list[['discrete_level']]
  
  
  
  datafile_name <- paste0("part1_", ID, '_', discrete_level, '_t', time_scale, '.rds')
  
  saveRDS(results, file = file.path(temp_file_dir, datafile_name))  
  
}

full_conditional_estimation_with_no_truth_part2 <- function(temp_file_dir, setting_info_list, cont_ind){

  # load 
  
  print('I am in part 2 code')
  
  ID <- setting_info_list[['ID']]
  y_c_structure <- setting_info_list[['y_c_structure']]
  time_scale <- setting_info_list[['time_scale']]
  method <- setting_info_list[['method']]
  discrete_level <- setting_info_list[['discrete_level']]
  
  datafile_name <- paste0("part1_", ID, '_', discrete_level, '_t', time_scale, '.rds')
  print(paste0('the datafile is', datafile_name))
  
  full_path <- file.path(temp_file_dir, datafile_name)
  print(paste0('the full path for the datafile is', full_path))
  
  # Check if it exists
  if (file.exists(full_path)) {
    print("File exists, safe to read")
  } else {
    print("File does NOT exist in temp_file_dir!")
  }
  
  results <- readRDS(file.path(temp_file_dir, datafile_name))
  
  print(paste0('size of part 1 data: ', length(results)))
  print('names of part 1 data')
  print(names(results))
  
  # load all variabels
  list2env(results, envir = .GlobalEnv)
  
  # start covariate loop 
  # ----------------------------------------------------------------------------
  
  

  cont_inds <- 1:nrow(query_y_cs)

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
  
  
  estimated_graphs <- list(step_2 = step_2, step_2b = step_2b, step_3 = step_3,
                           step_4 = step_4, step_5 = step_5, step_5b = step_5b, step_9 = step_9, step_9b = step_9b,
                           step_10 = step_10, step_11 = step_11)
  
  
  # save as part2_1, part2_2
  file_name <- paste0('part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind, '.rds')
  saveRDS(estimated_graphs, file = file.path(temp_file_dir, file_name))  
   
}

full_conditional_estimation_with_no_truth_part3 <- function(temp_file_dir, setting_info_list, cont_inds){
  
  #
  # cont_inds = number
  # 
  
  # ----------------------------------------------------------------------------
  # read from steps 1 and 2 and then remove everything
  # ----------------------------------------------------------------------------  
  
  # Vector of file names
  
  ID <- setting_info_list[['ID']]
  y_c_structure <- setting_info_list[['y_c_structure']]
  time_scale <- setting_info_list[['time_scale']]
  method <- setting_info_list[['method']]
  discrete_level <- setting_info_list[['discrete_level']]
  
  file_names <- paste0(temp_file_dir, '/part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
  
  
  # Step 1: Load all files into a list
  all_loaded <- lapply(file_names, readRDS)
  
  # Step 2: Get all step names (assumes all files have same names)
  step_names <- names(all_loaded[[1]])
  
  # Step 3: Reorganize by step
  estimated_graphs <- setNames(lapply(step_names, function(step) {
    lapply(all_loaded, `[[`, step)  # collect that step from all files
  }), step_names)
  

 
  # load step 1
  results <- readRDS(file.path(temp_file_dir, "part1.rds"))
  # load all variables
  list2env(results, envir = .GlobalEnv)
  
  
  # remove them 
  file.remove(file_names)
  
  part1_file_name <- paste0("part1_", ID, '_', discrete_level, '_t', time_scale, '.rds')
  
  file.remove(file.path(temp_file_dir, part1_file_name)) 
  
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

