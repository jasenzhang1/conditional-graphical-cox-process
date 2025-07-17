full_conditional_estimation <- function(discrete_strata, y_d_df, data_df){
  
  #
  #
  #
  #
  # - discrete_strata (K x q_d data.frame of discrete strata)
  # - y_d_df          (n x q_d data.frame of discrete covariates for all n subjects)
  # - data_df         (data.frame with 'feature_id', 'time', and 'subject_num')
  
  n <- dim(y_d_df)[1]

  for(y_ind in 1:nrow(discrete_strata)){ # 3) for each discrete variable level 
    
    # prepare times
    t6_total <- 0
    t7_total <- 0
    t8_total <- 0 
    
    t0 <- Sys.time()
    
    # 4) filtering by strata + included neurons to get data_df4
    
    yd_i <- discrete_strata[y_ind, ]
    
    subject_nums <- 1:n
    
    s_yd <- subject_nums[apply(y_d_df, 1, function(x) all(x == yd_i))]  # s_yd = indices
    data_df4 <- data_df %>% filter(subject_num %in% s_yd)

    
    
    # 5) counting neurons
    
    s_yd <- unique(data_df4$subject_num)
    y_c_strata_both <- y_c_df %>% filter(subject_num %in% s_yd) %>% dplyr::select(-subject_num) %>% as.matrix() # age and timestamp
    y_c_strata_age <- y_c_strata_both[, "age", drop = FALSE] # only age
    y_c_strata <- list(y_c_strata_both, y_c_strata_age)
    
    silent_neurons <- setdiff(1:neuron_count,
                              unique(data_df2$feature_id))
    
    missing_neurons <- setdiff(unique(data_df2$feature_id),
                               unique(data_df4$feature_id))
    
    inactive_neurons <- c(silent_neurons, missing_neurons)
    
    # size of dataset
    print(paste0('number of subjects: ', length(s_yd)))
    print(paste0('number of spikes: ', nrow(data_df4)))
    
    # print number of neurons
    print(paste('total neurons: ', unname(neuron_count)))                     # total neurons
    print(paste('active neurons: ', length(unique(data_df4$feature_id))))     # active neurons
    print(paste('discarded neurons: ', length(missing_neurons)))              # discarded neurons
    print(paste('silent neurons: ', length(silent_neurons)))                  # silent neurons
    print(paste('do they add up? ', length(unique(data_df4$feature_id)) + length(missing_neurons) + length(silent_neurons) == unname(neuron_count)))
    print('-------------------------------------')
    
    # 6) retroactively get parameters for subjects (ntrain, patient_sel) and neurons (p, feature_sel)
    
    ntrain <- length(unique(data_df4$subject_num))
    p <- length(unique(data_df4$feature_id))
    
    
    patient_sel = unique(data_df4$subject_num) %>% sort()
    feature_sel = unique(data_df4$feature_id) %>% sort()
    
    
    # 7) estimation 
    
    t1 <- Sys.time()
    
    rho_list <- estimate_intensities_stratum_parallel_v3(data_df4, patient_sel, feature_sel, Tseq, ncores)
    t2 <- Sys.time()
    
    
    time_elapsed <- round(as.numeric(difftime(t2, t1, units = 'mins')), 2)
    print(paste0('checkpoint 1: ', time_elapsed, ' mins'))
    
    # part 3
    g_ij_st <- estimate_covariance_functions_ii(rho_list)
    t3 <- Sys.time()
    
    # part 4
    eigen_decomp <- compute_eigendecomposition_ii(g_ij_st)
    t4 <- Sys.time()
    
    # part 5
    # redundant alphas are set to 0
    kl_coeffs <- estimate_kl_coefficients_parallel(data_df4, eigen_decomp$eigenfunctions, eigen_decomp$n_dims, 
                                                   patient_sel, feature_sel, Tseq, ncores)
    t5 <- Sys.time()
    
    time_elapsed <- round(as.numeric(difftime(t5, t2, units = 'mins')), 2)
    print(paste0('checkpoint 2: ', time_elapsed, ' mins'))    
    
    
    # here, we fit to y_c_strata and y_c_strata_age
    
    for(y_case in 1:length(y_c_strata)){
      
      t6a <- Sys.time()
      
      # part 6
      gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata[[y_case]])
      K_c <- construct_kernel_matrix_step_6(y_c_strata[[y_case]], gamma_c)
      
      
      t6 <- Sys.time()
      t6_total <- t6_total + round(as.numeric(difftime(t6, t6a, units = 'mins')), 2)
      
      # part 7
      
      # V_YcXij_og <- construct_cross_covariance_matrix(kl_coeffs)
      V_YcXij <- construct_cross_covariance_matrix_v2(kl_coeffs)
      # M_hat_og <- estimate_regression_operators(K_c, V_YcXij_og, gamma_c, p)
      M_hat <- estimate_regression_operators_v3(K_c, V_YcXij, gamma_c, p)
      
      t7 <- Sys.time()
      t7_total <- t7_total + round(as.numeric(difftime(t7, t6, units = 'mins')), 2)
      
      time_elapsed <- round(as.numeric(difftime(t7, t6a, units = 'mins')), 2)
      print(paste0('checkpoint 3: ', time_elapsed, ' mins'))  
      
      # part 8 - parallelize this
      
      cont_inds <- 1:nrow(query_y_cs[[y_case]])
      pbmclapply(cont_inds, function(cont_ind) {
        
        
        query_y_c <- query_y_cs[[y_case]][cont_ind, ] %>% as.numeric()
        
        # V_cond_og <- evaluate_regression_at_query(M_hat_og, y_c_strata, query_y_c, eigen_decomp$eigenfunctions, gamma_c, p)
        # V_cond <- evaluate_regression_at_query_diag_only(M_hat, y_c_strata, query_y_c, eigen_decomp$eigenfunctions, gamma_c, p)
        V_cond <- evaluate_regression_at_query_v2(M_hat, y_c_strata[[y_case]], query_y_c, eigen_decomp$eigenfunctions, gamma_c, p)
        
        
        # part 9
        gamma1 <- 0.01
        # C_cond_og <- estimate_conditional_correlation(V_cond_og, gamma1, p)
        C_cond <- estimate_conditional_correlation_v2(V_cond, gamma1, p)
        
        
        
        
        # part 10
        gamma2 <- 0.01
        # P_cond_og <- estimate_precision_operator(C_cond_og, gamma2, p)
        P_cond <- estimate_precision_operator_v2(C_cond, gamma2, p)
        
        
        
        # part 11
        threshold <- select_threshold_by_stability(P_cond, p)
        final_graph_estimates <- estimate_graph(P_cond, threshold, p, c(silent_neurons, missing_neurons))
        
        
        # 8) save graphs -------------------------------------------------------
        
        if(length(query_y_c) == 1){
          query_name <- paste0('age_', query_y_c[1], '_ts_NA')
        } else{
          query_name <- paste0('age_', query_y_c[1], '_ts_', query_y_c[2])
        }
        
        
        results_name <- paste0(discrete_strata_name, '_', query_name)
        
        
        
        save_file <- paste0(save_dir, ID2[i], '_', results_name, '.rda')
        save(final_graph_estimates, file=save_file)
        
        return(NULL)
        
      }, mc.cores = ncores) # done with all y_c levels
      
      t8 <- Sys.time()
      
      
      time_elapsed <- round(as.numeric(difftime(t8, t7, units = 'mins')), 2)
      print(paste0('checkpoint 4: ', time_elapsed, ' mins'))   
      t8_total <- t8_total + time_elapsed
      
    } # done with both cases (only age, both age and ts)
    
    t9 <- Sys.time()
    
    # 9) save runtimes ---------------------------------------------------------
    
    times_end   <- c(t1, t2, t3, t4, t5, t9)
    times_start <- c(t0, t1, t2, t3, t4, t0)    
    times_diff <- difftime(times_end, times_start, units = 'mins') %>% as.numeric() %>% round(2)
    
    times_diff <- c(times_diff[1:5], t6_total, t7_total, t8_total, times_diff[6])
    
    
    
    strata_stat_vec <- c(ID2[i],                              # mouse ID
                         yd_i$movement,                       # movement
                         yd_i$VR,                             # VR
                         length(s_yd),                        # num subjects
                         nrow(data_df4),                      # num total spikes
                         unname(neuron_count),                # num total neurons
                         length(unique(data_df4$feature_id)), # num active neurons
                         length(missing_neurons),             # discarded neurons
                         length(silent_neurons),              # silent neurons
                         times_diff)                          # all 8 runtimes + total
    
    
    
    strata_stats <- rbind(strata_stats, strata_stat_vec)
    
    time_elapsed <- round(as.numeric(difftime(t8, t1, units = 'mins')), 2)
    print(paste0('strata ', y_ind, ' total time: ', time_elapsed, ' mins'))   
    print('-------------------------------------')
  } # done with all discrete stratas
}