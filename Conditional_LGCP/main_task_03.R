t_start <- Sys.time()


# set seed, loading functions, data --------------------------------------------

set.seed(1)

source('functions/00_function_wrapper.R')

IDs <- c('346', '351', '366', '361', '362', '368')  # mouse ID
ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3') # our name 

num_neurons <- c(169, 250, 240, 249, 235, 294)

IDs <- c('346', '368')  # mouse ID
ID2 <- c('Tau1', 'WT3') # our name 

num_neurons <- c(169, 294)

task_name <- 'task_03'
data_folder_name <- 'data/with_ts/'
ncores <- parallel::detectCores() - 1
# ncores <- 6

# continuous variables to query ------------------------------------------------

age_query <- c(1/4, 2/4, 3/4)
timestamp_query <- c(1/4, 2/4, 3/4)
query_y_cs <- data.frame(rep(age_query, length(timestamp_query)),
                         rep(timestamp_query, each = length(age_query))) # age, timestamp
colnames(query_y_cs) <- c('age', 'timestamp')
  
included_neurons_vec <- c(1000)
Tseq = seq(0.05,0.95,length=19)

# starting ---------------------------------------------------------------------

print('Starting ----------------------------')
print(paste0('number of cores: ', ncores))

strata_stats <- data.frame()

for(i in 1:length(IDs)){ 
  
  # 1) # for every mouse, load data
  
  load(paste0(data_folder_name, ID2[i], '_data.rda'))
  neuron_count <- num_neurons[i]
  
  
  data_df <- LGCP_data[[1]] 
  y_d_df <- LGCP_data[[2]]
  y_c_df <- LGCP_data[[3]]
  
  
  
  # 1.1) obtain discrete strata from y_d_df
  
  cov_df <- y_d_df %>% dplyr::select(- subject_num)
  discrete_strata <- cov_df %>% unique()
  
  
  for(included_neurons in included_neurons_vec){ # 2) vary number of included neurons
    
    for(y_ind in 1:nrow(discrete_strata)){ # 3) for each discrete variable level 
      
      t0 <- Sys.time()
      
      # 4) filtering by strata + included neurons to get data_df4
      
      yd_i <- discrete_strata[y_ind, ]
      discrete_strata_name <- paste0('m', yd_i$movement, 'vr', yd_i$VR)
      
      s_yd <- y_d_df$subject_num[apply(cov_df, 1, function(x) all(x == yd_i))]  # s_yd = indices
      data_df2 <- data_df %>% filter(subject_num %in% s_yd)
      data_df4 <- data_df2 %>% filter(feature_id %in% 1:included_neurons)
    
      
      # 5) counting neurons
  
      s_yd <- unique(data_df4$subject_num)
      y_c_strata <- y_c_df %>% filter(subject_num %in% s_yd) %>% dplyr::select(-subject_num) %>% as.matrix()
      
      
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
      
      # checking
      # rho_list_old <- estimate_intensities_stratum_parallel_v2(data_df4, patient_sel, feature_sel, Tseq, ncores)
      # g_ij_st_old <- estimate_covariance_functions(rho_list_old$rho_hat, rho_list_old$rho_hat_pairs)
      # eigen_decomp_old <- compute_eigendecomposition(g_ij_st_old)
      # kl_coeffs_old <- estimate_kl_coefficients(data_df4, eigen_decomp_old$eigenfunctions, eigen_decomp_old$n_dims, 
      #                                           patient_sel, feature_sel, Tseq)      
      
      # part 6
      gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
      K_c <- construct_kernel_matrix_step_6(y_c_strata, gamma_c)
      t6 <- Sys.time()
      
      # part 7
  
      # V_YcXij_og <- construct_cross_covariance_matrix(kl_coeffs)
      V_YcXij <- construct_cross_covariance_matrix_v2(kl_coeffs)
      # M_hat_og <- estimate_regression_operators(K_c, V_YcXij_og, gamma_c, p)
      M_hat <- estimate_regression_operators_v3(K_c, V_YcXij, gamma_c, p)
      
      t7 <- Sys.time()
      
      
      
      
      time_elapsed <- round(as.numeric(difftime(t7, t2, units = 'mins')), 2)
      print(paste0('checkpoint 2: ', time_elapsed, ' mins'))  
      
      # part 8
      for(cont_ind in 1:nrow(query_y_cs)){
        
        query_y_c <- query_y_cs[cont_ind, ] %>% as.numeric()
  
        # V_cond_og <- evaluate_regression_at_query(M_hat_og, y_c_strata, query_y_c, eigen_decomp$eigenfunctions, gamma_c, p)
        # V_cond <- evaluate_regression_at_query_diag_only(M_hat, y_c_strata, query_y_c, eigen_decomp$eigenfunctions, gamma_c, p)
        V_cond <- evaluate_regression_at_query_v2(M_hat, y_c_strata, query_y_c, eigen_decomp$eigenfunctions, gamma_c, p)
  
        
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
        
        query_name <- paste0('age_', query_y_c[1], '_ts_', query_y_c[2])
        
        results_name <- paste0(discrete_strata_name, '_', query_name)
        
        save_dir <- 'results/'  
        if (!dir.exists(save_dir)) {
          dir.create(save_dir)
        }     
        
        save_dir <- paste0(save_dir, task_name, '/')  
        if (!dir.exists(save_dir)) {
          dir.create(save_dir)
        } 
        
        save_file <- paste0(save_dir, ID2[i], '_', results_name, '.rda')
        save(final_graph_estimates, file=save_file)
        
        print(paste0('done with ', results_name))
        
      } # done with all y_c levels
      
      t8 <- Sys.time()
    
      # 9) save runtimes ---------------------------------------------------------
      
      times_end   <- c(t1, t2, t3, t4, t5, t6, t7, t8, t8)
      times_start <- c(t0, t1, t2, t3, t4, t5, t6, t7, t0)
      times_diff <- difftime(times_end, times_start, units = 'mins') %>% as.numeric() %>% round(2)
      
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
    
    print(paste0('done with neuron inclusion: ', included_neurons))
    print('----------------------------------------------------------')
  } # done with included neuron strata
  
  print(paste0('done with mouse: ', ID2[i]))
  print('----------------------------------------------------------')  
  
} # done with mouse 

# 10) name the dataframe and save and end ---------------------------------------

colnames(strata_stats) <- c('Mouse', 'Movement', 'VR', 'Num_Rep', 'Num_Spikes', 
                            'Total_Neurons', 'Active_Neurons', 'Discarded_Neurons', 'Silent_Neurons', 
                            't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't_total')

save_dir <- paste0(save_dir, 'run_time/')
if (!dir.exists(save_dir)) {
  dir.create(save_dir)
} 

file_name <- paste0(save_dir, 'run_time_results.RData')

save(strata_stats, file = file_name)

t_end <- Sys.time()
t_total <- round(as.numeric(difftime(t_end, t_start, units = 'mins')), 2)
print(paste0('DONE, total time taken: ', t_total, ' mins'))






