# 9/9/2025

# we're not using G_ij anymore, but rather calculating a different cross covariance operator

# we're also calculating G_ii with rho_i_hat and rho_ii_hat

# we're also modifying the V_YcXij cross-covariance step to include w(y_c) in it

full_conditional_estimation_with_truths_v2 <- function(dataset, method, terse, ncores, dir){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: troubleshoot the full conditional estimation procedure with ground truths between steps
  #
  # - instead of previously, where we only used G_ii and discarded G_ij
  # - here we actively look to estimate G_ij and use it
  # - 9/3/2025 - use G_ij^k, for each subject
  #
  # Input:
  #
  # - dataset
  # - method   ('OG', 'JASA', 'CPGM')
  # - terse    (boolean) if true, return much less
  # - ncores
  #
  # 
  #
  # ----------------------------------------------------------------------------

  # check for errors 
  if (!(method %in% c("OG", "JASA", "CPGM"))) {
    stop("Error: 12c_full_conditional_estimation_with_truths_v2 - method must be either 'OG' or 'JASA' or 'CPGM'")
  }    
  
  # 0) load parameters
  
  time_grid <- dataset$simulation_params$time_grid
  time_grid_est <- dataset$simulation_params$time_grid_est
  time_grid_both <- dataset$simulation_params$time_grid_both
  
  n <- dim(dataset$Y_continuous)[1]
  m_est <- length(time_grid_est)
  m     <- length(time_grid)
  p <- dim(dataset$true_graphs[[1]]$P_block_kronecker$prec_mat_truth$adj_mat)[1]

  
  # assuming this is the truth all the time
  kernel_params <- dataset$true_graphs[[1]]$P_block_kronecker
  
  
  # part 0 - check ground truth precisions and process data --------------------
  
  step_0 <- step_0_store(dataset)
  step_0_check(dataset)
  
  # process data
  processed_data <- step_0_preprocess(dataset)
  data_df4 = processed_data[[1]]
  y_c_strata = processed_data[[2]]
  query_y_cs = processed_data[[3]]
  patient_sel = processed_data[[4]]
  feature_sel = processed_data[[5]]
  rm(processed_data)
  
  
 
  # part 1 - log intensities ---------------------------------------------------
  
  step_1 <- step_1_log_intensities(dataset, data_df4, time_grid_est, time_grid, time_grid_both)
  
  # part 2 - intensities -------------------------------------------------------
  
  # warning here
  rho_kernel <- F
  
  # hold = nothing, just can't be NULL
  step_2 <- step_2_rho_i(dataset, data_df4, kernel_params, rho_kernel, patient_sel, feature_sel, time_grid, time_grid_est, 'hold', 'hold', ncores)

  # part 2.1 - bivariate intensities -------------------------------------------

  i_neq_j <- F
  step_2b <- step_2_rho_ij(step_1, step_2, kernel_params, i_neq_j)

  
  # part 3 - GP covariance estimation ------------------------------------------
  
  step_3 <- step_3_g_ij(step_2, step_2b, kernel_params, i_neq_j)

  
  # part 4 - eigendecomposition of GP covariance -------------------------------
  
  # careful about some processes being empty
  step_4 <- step_4_eigendecomp(step_3, p, time_grid, time_grid_est)
  
 
  
  # part 5 - KL expansion ------------------------------------------------------
  
  step_5 <- step_5_KL_expansion(step_1, step_4, kernel_params, time_grid, time_grid_est, ncores)
  
  
  # now, we regress on y_c  ----------------------------------------------------
  

  
  print('at covariate loop')
  
  cont_inds <- 1:nrow(query_y_cs)
  estimated_graphs_v2 <- lapply(cont_inds, function(cont_ind) {
    
    
    print(paste0(cont_ind, ' out of ', nrow(query_y_cs)))
    
    
    # obtain this query's parameters
    kernel_params_i = dataset$true_graphs[[cont_ind]]$P_block_kronecker
    adj_mat_i <- dataset$true_graphs[[cont_ind]]$adj_mat
    query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()
    
    step_8 <- steps_78(step_4, step_5, kernel_params_i, y_c_strata, query_y_c, method, ncores)
    

    # part 9 - correlation operator --------------------------------------------
    
    step_9 <- step_9_C_cond_from_V_cond(step_8, kernel_params_i)
    
    
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
                                step_8$V_cond_est_full,    # V_cond_est_full
                                step_10$P_cond_coarse_ground_truth_full,   # P_cond_coarse_ground_truth_full
                                step_9$C_cond_coarse_ground_truth_full,    # C_cond_coarse_ground_truth_full
                                step_8$V_cond_coarse_ground_truth_full,    # V_cond_coarse_ground_truth_full
                                step_11$w_mat_est,   # w_mat_est
                                adj_mat_i,
                                step_12$roc_est,   # roc_est
                                step_2$rho_i_est,    # rho_i_est
                                step_2$rho_i_coarse_truth,    # rho_i_coarse_truth
                                step_2b$rho_ii_est,   # rho_ii_est
                                step_2b$rho_ii_coarse_truth,   # rho_ii_coarse_truth
                                step_3$g_ij_est,    # g_ij_est
                                step_3$g_ij_coarse_truth,  # g_ij_coarse_truth
                                kernel_params_i$base_cov_est,
                                kernel_params_i$prec_mat_truth$cor_mat, 
                                kernel_params_i$base_precision_est,
                                kernel_params_i$prec_mat_truth$prec_mat)
                           )
    

    # 8) return all results
    

    list(step_8 = step_8,
         step_9 = step_9,
         step_10 = step_10,
         step_11 = step_11,
         step_12 = step_12,
         metrics = metrics)
   
    
    
    
  })# , mc.cores = ncores) # done with all y_c levels
  
  
  
  names(estimated_graphs_v2) <- apply(query_y_cs, 1, function(x){paste(x, collapse = '_')})
  
  
  # what to output
  
  estimated_graphs_part_1 <- list(step_0 = step_0,
                                  step_1 = step_1,
                                  step_2 = step_2,
                                  step_2b = step_2b,
                                  step_3 = step_3,
                                  step_4 = step_4,
                                  step_5 = step_5)
  
  true_K_base <- lapply(dataset$true_graphs, function(x) x$P_block_kronecker$base_cov)
  true_K_adj <- lapply(dataset$true_graphs, function(x) x$P_block_kronecker$prec_mat_truth$cor_mat)
  
  return(list(estimated_graphs_part_1 = estimated_graphs_part_1,      
              estimated_graphs_part_2 = estimated_graphs_v2,
              true_K_base = true_K_base,
              true_K_adj = true_K_adj))  
  

}

full_conditional_estimation_with_truths_v3 <- function(dataset, method, terse, ncores, dir){
  
  
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
  # - dataset
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
  
  # load parameters
  
  time_grid <- dataset$simulation_params$time_grid
  time_grid_est <- dataset$simulation_params$time_grid_est
  time_grid_both <- dataset$simulation_params$time_grid_both
  
  n <- dim(dataset$Y_continuous)[1]
  m_est <- length(time_grid_est)
  m     <- length(time_grid)
  p <- dim(dataset$true_graphs[[1]]$P_block_kronecker$prec_mat_truth$adj_mat)[1]
  
  # part 0 - check ground truth precisions and process data --------------------
  
  step_0 <- step_0_store(dataset)
  step_0_check(dataset)
  

  # process data
  processed_data <- step_0_preprocess(dataset)
  data_df4 = processed_data[[1]]
  y_c_strata = processed_data[[2]]
  query_y_cs = processed_data[[3]]
  patient_sel = processed_data[[4]]
  feature_sel = processed_data[[5]]
  rm(processed_data)
  

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
                                step_3$g_ij_coarse_truth,  # g_ij_coarse_truth
                                kernel_params_i$base_cov_est,
                                kernel_params_i$prec_mat_truth$cor_mat, 
                                kernel_params_i$base_precision_est,
                                kernel_params_i$prec_mat_truth$prec_mat))   # g_ij_coarse_truth     
    
    
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
         kernel_params = kernel_params_i)
    
    
    
  })# , mc.cores = ncores) # done with all y_c levels
  
  
  
  names(estimated_graphs_v2) <- apply(query_y_cs, 1, function(x){paste(x, collapse = '_')})
  
  # reverse the order
  # - estimated_graphs_v2 has steps first, then y_c's second
  
  steps <- unique(unlist(lapply(estimated_graphs_v2, names)))  # get all step names
  
  reorganized_estimated_graphs_v2 <- setNames(
    lapply(steps, function(step) {
      # for each step, collect that step from all numeric groups
      sapply(estimated_graphs_v2, `[[`, step, simplify = FALSE)
    }),
    steps
  )
  
  rm(estimated_graphs_v2)
    
    
  
  # what to output
  
  estimated_graphs_part_1 <- list(step_0 = step_0, step_1 = step_1)
  
  
  return(c(estimated_graphs_part_1, reorganized_estimated_graphs_v2))      
    

  
}

full_conditional_estimation_with_truths_v4 <- function(dataset, method, terse, ncores, dir) {
  
  # ----------------------------------------------------------------------------
  # Unified version of v2 and v3
  # 
  #
  # GOAL: full conditional estimation procedure with ground truths between steps
  #
  #
  # Supports:
  #   - "OG", "JASA"  → covariate loop after step_5 (v2-style)
  #   - "CPGM"        → covariate loop after step_1 (v3-style)
  #
  # Input:
  #
  # - dataset
  # - method   ('CPGM')
  # - terse    (boolean) if true, return much less
  # - ncores
  # - dir      (string) folder name, such as "simu_results_banded_trig2_2"
  # 
  #
  # ----------------------------------------------------------------------------
  
  # Validate method
  if (!(method %in% c("OG", "JASA", "CPGM"))) {
    stop("Error: method must be 'OG', 'JASA', or 'CPGM'")
  }
  
  # Load key parameters
  time_grid <- dataset$simulation_params$time_grid
  time_grid_est <- dataset$simulation_params$time_grid_est
  time_grid_both <- dataset$simulation_params$time_grid_both
  
  n <- dim(dataset$Y_continuous)[1]
  m_est <- length(time_grid_est)
  m <- length(time_grid)


  
  # ----------------------------------------------------------------------------
  # Step 0 - Check and preprocess data
  # ----------------------------------------------------------------------------
  
  step_0 <- step_0_store(dataset)
  step_0_check(dataset)
  
  processed_data <- step_0_preprocess(dataset)
  data_df4     <- processed_data[[1]]
  y_c_strata   <- processed_data[[2]]
  query_y_cs   <- processed_data[[3]]
  patient_sel  <- processed_data[[4]]
  feature_sel  <- processed_data[[5]]
  rm(processed_data)
  
  # ----------------------------------------------------------------------------
  # Step 1 - Log intensities
  # ----------------------------------------------------------------------------
  
  step_1 <- step_1_log_intensities(dataset, data_df4, time_grid_est, time_grid, time_grid_both)
  
  # start covariate loop 
  # ----------------------------------------------------------------------------
  
  print('at covariate loop')
  cont_inds <- 1:nrow(query_y_cs)
  estimated_graphs <- lapply(cont_inds, function(cont_ind) {
    
    kernel_params_i <- dataset$true_graphs[[cont_ind]]$P_block_kronecker
    p               <- dim(kernel_params_i$prec_mat_truth$adj_mat)[1]
    adj_mat_i       <- dataset$true_graphs[[cont_ind]]$adj_mat
    query_y_c       <- query_y_cs[cont_ind, ] %>% as.numeric()
    
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
                           query_y_c, y_c_strata, ncores)           

    
    # Step 2b onward
    step_2b <- step_2_rho_ij(step_1, step_2, kernel_params_i, i_neq_j)
    step_3  <- step_3_g_ij(step_2, step_2b, kernel_params_i, i_neq_j)
    
    step_4 <- tryCatch({
      step_4_eigendecomp(step_3, p, time_grid, time_grid_est)
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
      step_5 <- step_5_KL_covariance(step_3, step_4)
      step_9 <- step_9_C_cond_from_KL_cov(step_4, step_5, kernel_params_i)
    }
    
    # reunite at step 10 onwards
    
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
    
    
    step_11 <- step_11_HS_norms(step_10, adj_mat_i, p)
    step_12 <- step_12_ROC(step_11, adj_mat_i)
    
    metrics_list <- list(
      step_10$P_cond_est_full,
      step_9$C_cond_est_full,
      NULL,
      step_10$P_cond_coarse_ground_truth_full,
      step_9$C_cond_coarse_ground_truth_full,
      NULL,
      step_11$w_mat_est,
      adj_mat_i,
      step_12$roc_est,
      step_2$rho_i_est,
      step_2$rho_i_coarse_truth,
      step_2b$rho_ii_est,
      step_2b$rho_ii_coarse_truth,
      step_3$g_ij_est,
      step_3$g_ij_coarse_truth,
      kernel_params_i$base_cov_est,
      kernel_params_i$prec_mat_truth$cor_mat,
      kernel_params_i$base_precision_est,
      kernel_params_i$prec_mat_truth$prec_mat
    )   
    
    if(method %in% c('OG', 'JASA')){
      metrics_list[[3]] <- step_8$V_cond_est_full
      metrics_list[[6]] <- step_8$V_cond_coarse_ground_truth_full
    }  
    
    metrics <- get_metrics(metrics_list)
    
    
    
    list(step_2 = step_2, step_2b = step_2b, step_3 = step_3,
         step_4 = step_4, step_5 = step_5, step_9 = step_9,
         step_10 = step_10, step_11 = step_11, step_12 = step_12,
         metrics = metrics, kernel_params = kernel_params_i)
  })
  
  # ----------------------------------------------------------------------------
  # Optional: reorganize by steps
  # ----------------------------------------------------------------------------
  
  steps <- unique(unlist(lapply(estimated_graphs, names)))
  reorganized <- setNames(lapply(steps, function(step) {
    sapply(estimated_graphs, `[[`, step, simplify = FALSE)
  }), steps)
  
  # Return step_0, step_1 (shared) + reorganized per-subject steps
  estimated_graphs_part_1 <- list(step_0 = step_0, step_1 = step_1)
  return(c(estimated_graphs_part_1, reorganized))

}

