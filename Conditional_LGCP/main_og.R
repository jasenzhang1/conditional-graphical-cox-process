# a fitting model with slow functions

t_start <- Sys.time()
ffs = 1
set.seed(ffs)


source('functions/00_function_wrapper.R')

load('data/WT1_data.rda')
ID2 <- 'WT1'



ncores <- parallel::detectCores() - 1
# ncores <- 6

print('Starting ----------------------------')
print(paste0('number of cores: ', ncores))

data_df <- LGCP_data[[1]] 
y_d_df <- LGCP_data[[2]]
y_c_df <- LGCP_data[[3]]

# what to query
query_y_cs <- c(17, 23, 28)
included_neurons_vec <- c(5, 20, 50, 100)

# extract the binary levels
cov_df <- y_d_df %>% dplyr::select(- subject_num)

discrete_strata <- cov_df %>% unique()

strata_stats <- data.frame()

for(included_neurons in included_neurons_vec){ # vary number of included neurons
  
  for(y_ind in 1:nrow(discrete_strata)){ # 1) for each discrete variable level 
    
    t0 <- Sys.time()
    
    yd_i <- discrete_strata[y_ind, ]
    discrete_strata_name <- paste0('m', yd_i$movement, 'vr', yd_i$VR)
    
    s_yd <- y_d_df$subject_num[apply(cov_df, 1, function(x) all(x == yd_i))]  # s_yd = indices
    data_df2 <- data_df %>% filter(subject_num %in% s_yd)
    
    y_c_strata <- y_c_df %>% filter(subject_num %in% s_yd) %>% dplyr::select(-subject_num) %>% as.matrix()
    
    
    # 3) counting neurons
    
    data_df4 <- data_df2 %>% filter(feature_id %in% 1:included_neurons)
    
    s_yd <- unique(data_df4$subject_num)
    y_c_strata <- y_c_df %>% filter(subject_num %in% s_yd) %>% dplyr::select(-subject_num) %>% as.matrix()
    
    neuron_count <- 249
    
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
    
    # retroactively get parameters
    
    ntrain <- length(unique(data_df4$subject_num))
    p <- length(unique(data_df4$feature_id))
    
    
    ### get corr matrix for full data
    patient_sel = unique(data_df4$subject_num) %>% sort()
    feature_sel = unique(data_df4$feature_id) %>% sort()
    
    #Tseq = seq(0.01,0.99,length=99)
    Tseq = seq(0.05,0.95,length=19)
    
    dmax = 4
    FVE_thre = 0.9 # originally 0.9  
    
    # estimation of rho
    
    t1 <- Sys.time()
    
    rho_list <- estimate_intensities_stratum_parallel(data_df4, patient_sel, feature_sel, Tseq, ncores)
    t2 <- Sys.time()
    
    
    time_elapsed <- round(as.numeric(difftime(t2, t1, units = 'mins')), 2)
    print(paste0('checkpoint 1: ', time_elapsed, ' mins'))
    
    # part 3
    g_ij_st <- estimate_covariance_functions(rho_list$rho_hat, rho_list$rho_hat_pairs)
    t3 <- Sys.time()
    
    # part 4
    eigen_decomp <- compute_eigendecomposition(g_ij_st)
    t4 <- Sys.time()
    
    # part 5
    # redundant alphas are set to 0
    kl_coeffs <- estimate_kl_coefficients(data_df4, eigen_decomp$eigenfunctions, eigen_decomp$n_dims, 
                                          patient_sel, feature_sel, Tseq)
    t5 <- Sys.time()
    
    # part 6
    gamma_c <- select_gamma_c_bandwidth(y_c_strata)
    K_c <- construct_kernel_matrix_step_6(y_c_strata, gamma_c)
    t6 <- Sys.time()
    
    # part 7


    V_YcXij <- construct_cross_covariance_matrix(kl_coeffs)
    M_hat <- estimate_regression_operators(K_c, V_YcXij, gamma_c, p)
    
    t7 <- Sys.time()
    
    
    
    
    time_elapsed <- round(as.numeric(difftime(t7, t2, units = 'mins')), 2)
    print(paste0('checkpoint 2: ', time_elapsed, ' mins'))  
    
    # part 8
    for(query_y_c in query_y_cs){
      

      V_cond <- evaluate_regression_at_query(M_hat, y_c_strata, query_y_c, eigen_decomp$eigenfunctions, gamma_c, p)

      
      # part 9
      gamma1 <- 0.01
      C_cond <- estimate_conditional_correlation(V_cond, gamma1, p)
      

      # part 10
      gamma2 <- 0.01
      P_cond <- estimate_precision_operator(C_cond, gamma2, p)
      
      # part 11
      threshold <- select_threshold_by_stability(P_cond, p)
      final_graph_estimates <- estimate_graph(P_cond, threshold, p)
      
      results_name <- paste0(discrete_strata_name, '_w', query_y_c)
      
      # save 
      
      save_dir <- 'results/'  
      if (!dir.exists(save_dir)) {
        dir.create(save_dir)
      }      
      # save_file <- paste0(save_dir, 'WT1_', results_name, '.rda')
      # save(final_graph_estimates, file=save_file)   
      
      print(paste0('done with ', results_name))
      
    }
    
    t8 <- Sys.time()
  
    
    times_end   <- c(t1, t2, t3, t4, t5, t6, t7, t8)
    times_start <- c(t0, t1, t2, t3, t4, t5, t6, t7)
    times_diff <- difftime(times_end, times_start, units = 'mins') %>% as.numeric() %>% round(2)
    
    strata_stat_vec <- c(ID2,            # mouse ID
                         yd_i$movement,   # movement
                         yd_i$VR,        # VR
                         length(s_yd),   # num subjects
                         nrow(data_df4), # num total spikes
                         unname(neuron_count),  # num total neurons
                         length(unique(data_df4$feature_id)), # num active neurons
                         length(missing_neurons),             # discarded neurons
                         length(silent_neurons),              # silent neurons
                         times_diff)                          # all 8 times
    
  
  
    strata_stats <- rbind(strata_stats, strata_stat_vec)
    
    time_elapsed <- round(as.numeric(difftime(t8, t1, units = 'mins')), 2)
    print(paste0('strata ', y_ind, ' total time: ', time_elapsed, ' mins'))   
    print('-------------------------------------')
  }
  
  print(paste0('done with neuron inclusion: ', included_neurons))
  print('----------------------------------------------------------')
} # done with included neuron strata


colnames(strata_stats) <- c('Mouse', 'Movement', 'VR', 'Num_Rep', 'Num_Spikes', 
                            'Total_Neurons', 'Active_Neurons', 'Discarded_Neurons', 'Silent_Neurons', 
                            't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8')

save(strata_stats, file = "results/runtime/run_time_results.RData")

t_end <- Sys.time()
t_total <- round(as.numeric(difftime(t_end, t_start, units = 'mins')), 2)
print(paste0('DONE, total time taken: ', t_total, ' mins'))






