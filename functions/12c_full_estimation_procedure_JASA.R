# 9/9/2025

# we're not using G_ij anymore, but rather calculating a different cross covariance operator

# we're also calculating G_ii with rho_i_hat and rho_ii_hat

# we're also modifying the V_YcXij cross-covariance step to include w(y_c) in it



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
  step_0_events <- step_0_keep_events(dataset, k = 1, i_vec = 1:5)
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
                           query_y_c_step_2, y_c_strata_step_2, ncores)           

    
    # Step 2b onward
    step_2b <- step_2_rho_ij(step_1, step_2, kernel_params_i, i_neq_j)
    step_3  <- step_3_g_ij(step_2, step_2b, kernel_params_i, i_neq_j)
    
    
    norm_G <- T # try normalizing G_ii = G_ii * diag(delta_t)
    norm_vec <- F # normalize eta so that Delta * \eta^\top * \eta = 1

    step_4 <- tryCatch({
      step_4_eigendecomp(step_3, p, time_grid, time_grid_est, norm_G, norm_vec)
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
      step_5 <- step_5_KL_covariance(step_3, step_4, norm_G)
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
      kernel_params_i$prec_mat_truth$prec_mat,
      step_9$C_cond_coarse_ground_truth_full_v2
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
  estimated_graphs_part_1 <- list(step_0 = step_0, 
                                  step_0_events = step_0_events,
                                  step_1 = step_1)
  
  all_results <- c(estimated_graphs_part_1, reorganized)
  all_results$y_c_query <- query_y_cs
  return(all_results)

}

