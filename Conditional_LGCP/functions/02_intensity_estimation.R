library(future.apply)
library(dplyr)
library(data.table)
library(pbmcapply)


estimate_intensities_stratum_parallel_v4 <- function(data_all, patient_sel, feature_sel, 
                                                     t_seq, i_neq_j, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  #  
  # GOAL: estimate the first and second order intensities
  #
  # - using pbmclapply
  # - rho_i  (p x m matrix)
  # - rho_{ij} (m x m matrix) for all (p x p) pairs
  #
  # - in v3, we only care about i = j entries for bivariate estimates
  #
  # - 8/5/2025
  #   need to tune gamma for KDE 
  #  
  # - 8/27/2025
  #   in v4, we want to allow for i =/= j entries to be included
  #
  # Input: 
  #
  #
  # - data_all          (data.table)      'feature_id', 'time', 'subject_num' 
  # - patient_sel       (vector)          ID's of non-discarded subjects
  # - feature_sel       (vector)          ID's of non-discarded features
  # - t_seq             (m-dim vector) 
  # - i_neq_j           (boolean)         whether to include i =/= j bivariate entries
  #
  #
  #
  # Output: 
  #
  # - list (list of two lists, one for univariate, one for bivariate)
  #     - first list   (values are m-dim vector)           univariate intensity, items are called 'i'
  #     - second list  (values are m x m dim matrix)        bivariate intensity, items are called 'i_j'
  #
  # ---------------------------------------------------------------------------
  
  
  NN = length(patient_sel) # 89
  p <- length(feature_sel) # 249
  n_time <- length(t_seq)  # 19
  
  # Step 1: obtain bivariate keys
  if(i_neq_j){ # if we allow i_j , only store cases i <= j
    
    keys = c()
    key_df <- data.frame()
    
    for(i in 1:p){
      for(j in i:p){
        key <- paste0(i, '_', j)
        keys <- c(keys, key)
        key_df <- rbind(key_df, c(i, j))
      }
    }
  } else{ # if we only allow i_i
    keys <- paste0(1:p, '_', 1:p)
    key_df <- data.frame(1:p, 1:p)
  }
    

  
  
  # step 2: univariate case 
  rho_i_list <- pbmclapply(1:p, function(i) {
    data_i <- data_all[feature_id == feature_sel[i], ]
    
    if (nrow(data_i) == 0) {
      rho_i <- rep(0, n_time)
    } else {
      Gamma_i <- data_i[, estimate_density(time, t_seq), by = "subject_num"]
      rho_mat <- matrix(Gamma_i$rho_hat, nrow = n_time)
      rho_i <- apply(rho_mat, 1, sum) / NN
    }
    
    rho_i
  }, mc.cores = ncores)
  
  # step 3: bivariate case - loop over all keys 
  rho_ij_list <- pbmclapply(1:nrow(key_df), function(k){
    i <- key_df[k,1]
    j <- key_df[k,2]
    
    data_i <- data_all[feature_id == feature_sel[i], ]
    data_j <- data_all[feature_id == feature_sel[j], ]
    
    
    # fitting
    if (nrow(data_j) == 0 || nrow(data_i) == 0) {
      rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
    } else {
      times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]  # dataframe where first col = subject, 2nd col = vector of observations 
      times_j <- data_j[, .(event_times_j = list(time)), by = subject_num]  # all of which are for process i
      times_ij <- merge(times_i, times_j, by = "subject_num", all = TRUE)   # now make it 3 columns: subject, process i, and process j
      
      Gamma_ij <- times_ij[, estimate_bivariate_density(
        event_times_i[[1]], event_times_j[[1]],
        t_seq, t_seq, 'i'), by = 'subject_num']
      
      bivariate_intensity <- matrix(Gamma_ij$V1, nrow = n_time^2)
      rho_ij <- apply(bivariate_intensity, 1, sum) / NN
      rho_ij_mat <- matrix(rho_ij, nrow = n_time)
    }
    rho_ij_mat # return this
  }, mc.cores = ncores)
  
  names(rho_i_list) <- 1:p
  names(rho_ij_list) <- keys
  
  return(list(rho_i_list, rho_ij_list))  
  
}

estimate_intensities_stratum_parallel_with_yc <- function(data_all, y_c_query, y_c_strata, 
                                                          patient_sel, feature_sel, 
                                                          t_seq, i_neq_j, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  #  
  # GOAL: estimate the first and second order intensities
  #
  # - using pbmclapply
  # - rho_i  (p x m matrix)
  # - rho_{ij} (m x m matrix) for all (p x p) pairs
  #
  # - in v3, we only care about i = j entries for bivariate estimates
  #
  # - 8/5/2025
  #   need to tune gamma for KDE 
  #  
  # - 8/27/2025
  #   in v4, we want to allow for i =/= j entries to be included
  #
  # - 9/22/2025
  #   we add y_c and y_c_strata
  #
  #
  # Input: 
  #
  #
  # - data_all          (data.table)        'feature_id', 'time', 'subject_num' 
  # - y_c_query         (q_c-dim vector)    queried y_c vector
  # - y_c_strata        (n x q_c matrix)    matrix of y_c_k values for each subject
  # - patient_sel       (vector)            ID's of non-discarded subjects
  # - feature_sel       (vector)            ID's of non-discarded features
  # - t_seq             (m-dim vector)      time discretization
  # - i_neq_j           (boolean)           whether to include i =/= j bivariate entries
  # - ncores            (integer)           number of cores
  #
  #
  #
  # Output: 
  #
  # - list (list of two lists, one for univariate, one for bivariate)
  #     - first list   (values are m-dim vector)           univariate intensity, items are called 'i'
  #     - second list  (values are m x m dim matrix)        bivariate intensity, items are called 'i_j'
  #
  # ---------------------------------------------------------------------------
  
  
  NN = length(patient_sel) # 89
  p <- length(feature_sel) # 249
  n_time <- length(t_seq)  # 19
  
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  
  # Step 1: obtain bivariate keys
  if(i_neq_j){ # if we allow i_j , only store cases i <= j
    
    keys = c()
    key_df <- data.frame()
    
    for(i in 1:p){
      for(j in i:p){
        key <- paste0(i, '_', j)
        keys <- c(keys, key)
        key_df <- rbind(key_df, c(i, j))
      }
    }
  } else{ # if we only allow i_i
    keys <- paste0(1:p, '_', 1:p)
    key_df <- data.frame(1:p, 1:p)
  }
  
  
  # get weights in preparation 
  weights <- apply(y_c_strata, 1, function(row) {
    step_6_kernel(as.numeric(row), y_c_query, gamma_c) 
  })    
  weights2 <- weights / sum(weights) # normalize
  
  
  # step 2: univariate case 
  rho_i_list <- pbmclapply(1:p, function(i) {
    data_i <- data_all[feature_id == feature_sel[i], ]
    
    if (nrow(data_i) == 0) {
      rho_i <- rep(0, n_time)
    } else {
      Gamma_i <- data_i[, estimate_density(time, t_seq), by = "subject_num"]
      rho_mat <- matrix(Gamma_i$rho_hat, nrow = n_time)
      
      

      rho_mat2 <- sweep(rho_mat, 2, weights2, `*`) # multiply each 19-dim vec by its normalized weight
      
      rho_i <- apply(rho_mat2, 1, sum) # no longer divide by NN
    }
    
    rho_i
  }, mc.cores = ncores)
  
  # step 3: bivariate case - loop over all keys 
  rho_ij_list <- pbmclapply(1:nrow(key_df), function(k){
    i <- key_df[k,1]
    j <- key_df[k,2]
    
    data_i <- data_all[feature_id == feature_sel[i], ]
    data_j <- data_all[feature_id == feature_sel[j], ]
    
    
    # fitting
    if (nrow(data_j) == 0 || nrow(data_i) == 0) {
      rho_ij_mat <- matrix(0, nrow = n_time, ncol = n_time)
    } else {
      times_i <- data_i[, .(event_times_i = list(time)), by = subject_num]  # dataframe where first col = subject, 2nd col = vector of observations 
      times_j <- data_j[, .(event_times_j = list(time)), by = subject_num]  # all of which are for process i
      times_ij <- merge(times_i, times_j, by = "subject_num", all = TRUE)   # now make it 3 columns: subject, process i, and process j
      
      Gamma_ij <- times_ij[, estimate_bivariate_density(
        event_times_i[[1]], event_times_j[[1]],
        t_seq, t_seq, 'i'), by = 'subject_num']
      
      bivariate_intensity <- matrix(Gamma_ij$V1, nrow = n_time^2)
      
    
      bivariate_intensity2 <- sweep(bivariate_intensity, 2, weights2, `*`) # multiply each 19-dim vec by its normalized weight      
      
      rho_ij <- apply(bivariate_intensity2, 1, sum) # no longer divide by NN
      rho_ij_mat <- matrix(rho_ij, nrow = n_time)
    }
    rho_ij_mat # return this
  }, mc.cores = ncores)
  
  names(rho_i_list) <- 1:p
  names(rho_ij_list) <- keys
  
  return(list(rho_i_list, rho_ij_list))  
  
}

