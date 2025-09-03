# 8/28/2025
# - for the new estimation procedure, we need to get G_{i, j}(s,t)^k 
# - for each subject

subject_specific_log_intensity <- function(data_df4, Tseq_est){
  
  # ---------------------------------------------------------------------------- 
  #
  #
  # GOAL: calcualte subject-specific log intensities for all (n) subjects and (p) processes
  # 
  #
  # input:
  # 
  # - data_df4   ('time', 'feature_id', 'subject_num' dataframe)
  # - Tseq_est   (m-dim vector)   time discretizations
  # 
  # output:
  # 
  # - X_k_est  (p x m x n)  matrix of estimated log-intensities for each subject (n) and each process (p)
  #
  # 
  # ----------------------------------------------------------------------------
  
  m_est <- length(Tseq_est)
  n <- length(unique(data_df4$subject_num))
  p <- length(unique(data_df4$feature_id))
  
  result_list <- data_df4 %>%
    group_by(feature_id, subject_num) %>%
    summarise(
      vec = list(estimate_log_intensity_function(time, Tseq_est)),  
      .groups = "drop"
    )
  
  
  X_k_est <- array(
    data = unlist(result_list$vec),
    dim = c(m_est, n, p)
  ) %>% aperm(c(3, 1, 2))  
  
  return(X_k_est)
  
}

subject_specific_density <- function(data_df4, Tseq_est){
  
  # ---------------------------------------------------------------------------- 
  #
  #
  # GOAL: calculate subject-specific densities (Gamma_k(t)) for all (n) subjects and (p) processes
  #
  # 
  # data_df4 ('time', 'feature_id', 'subject_num')
  #
  # 
  # output:
  # 
  # - Gamma_k_est  (p x m x n)  matrix of estimated densities for each subject (n) and each process (p)
  #
  # 
  # ----------------------------------------------------------------------------
  
  m_est <- length(Tseq_est)
  n <- length(unique(data_df4$subject_num))
  p <- length(unique(data_df4$feature_id))
  
  result_list <- data_df4 %>%
    group_by(feature_id, subject_num) %>%
    summarise(
      vec = list(estimate_density(time, Tseq_est)$gamma_hat),  
      .groups = "drop"
    )
  
  # Gamma^i_k(t) = density for process i, subject k 
  Gamma_k_est <- array(
    data = unlist(result_list$vec),
    dim = c(m_est, n, p)
  ) %>% aperm(c(3, 1, 2))  
  
  return(Gamma_k_est)
  
}

subject_specific_intensity <- function(data_df4, Tseq_est){
  
  # ---------------------------------------------------------------------------- 
  #
  #
  # GOAL: calculate subject-specific intensities (rho_k(t)) for all (n) subjects and (p) processes
  #
  # 
  # data_df4 ('time', 'feature_id', 'subject_num')
  #
  # 
  # output:
  # 
  # - rho_k_est  (p x m x n)  matrix of estimated densities for each subject (n) and each process (p)
  #
  # 
  # ----------------------------------------------------------------------------
  
  m_est <- length(Tseq_est)
  n <- length(unique(data_df4$subject_num))
  p <- length(unique(data_df4$feature_id))
  
  result_list <- data_df4 %>%
    group_by(feature_id, subject_num) %>%
    summarise(
      vec = list(estimate_density(time, Tseq_est)$rho_hat),  
      .groups = "drop"
    )
  
  # rho_i^k(t) = intensity for process i, subject k 
  rho_k_est <- array(
    data = unlist(result_list$vec),
    dim = c(m_est, n, p)
  ) %>% aperm(c(3, 1, 2))  
  
  return(rho_k_est)
  
}

subject_specific_bivariate_intensity <- function(data_df4, Tseq_est, ncores){
  
  # ---------------------------------------------------------------------------- 
  #
  #
  # GOAL: calculate subject-specific intensities (rho_ij_k(t)) for all (n) subjects and pairs of processes
  #
  # 
  # data_df4 ('time', 'feature_id', 'subject_num')
  #
  # 
  # output:
  # 
  # - rho_k_est  (p x m x n)  matrix of estimated densities for each subject (n) and each process (p)
  #
  # 
  # ----------------------------------------------------------------------------
  
  
  m_est <- length(Tseq_est)
  n <- length(unique(data_df4$subject_num))
  p <- length(unique(data_df4$feature_id))  
  
  pairs <- which(outer(1:p, 1:p, function(i,j) i <= j), arr.ind = TRUE)
  
  rho_ij_k_est <- mclapply(1:nrow(pairs), function(k) {
    i <- pairs[k, 1]
    j <- pairs[k, 2]
    
    event_times_i <- data_df4$time[data_df4$subject_num == i]
    
    # 1) find rho_ij_k 
    data_ij <- data_df4 %>% filter(feature_id %in% c(i, j))
    event_times_i <- data_df4 %>% filter()
    
    paste0("i=", i, ", j=", j)
  }, mc.cores = ncores)
  
  estimate_bivariate_density <- function(event_times_i, event_times_j, 
                                         eval_grid_s, eval_grid_t, d_or_i)  
  

  
  result_list <- data_df4 %>%
    group_by(feature_id, subject_num) %>%
    summarise(
      vec = list(estimate_density(time, Tseq_est)$rho_hat),  
      .groups = "drop"
    )
  
  # rho_i^k(t) = intensity for process i, subject k 
  rho_k_est <- array(
    data = unlist(result_list$vec),
    dim = c(m_est, n, p)
  ) %>% aperm(c(3, 1, 2))  
  
  return(rho_k_est)
  
}

estimate_subject_specific_covariance <- function(rho_k_est, rho_ij_k_est) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate G_{i, j}(s,t)^k 
  #         - for each i <= j
  #         - for each subject k
  #
  #
  # inputs:
  #
  # - rho_k_est  (p x m x n matrix) matrix of intensities for p processes, m time discretizations, n subjects 
  # - rho_ij_k_est ()        bivariate intensities
  # 
  #
  # output:
  #
  # -G_ij_k_est (list of ) each item is 'i_j' where i <= j, then each additional item has 'k' for each subject, for which it's a m x m matrix
  #
  #
  #
  # ----------------------------------------------------------------------------
  
  
  
  return(G_subject_specific)
}

estimate_subject_specific_covariance_OLD <- function(event_times, stratum_subjects, 
                                                 time_grid, bandwidth, p) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate G_{i, j}(s,t)^k 
  #         - for each i <= j
  #         - for each subject k
  #
  #
  # inputs:
  #
  # - event_times (n-dim vector)
  # - 
  #
  # df_estimate <- convert_data_for_estimation(dataset$subject_data) %>% as.data.table()
  # subj_1 <- df_estimate %>% filter(subject_num == 1)
  # event_times <- subj_1$time
  #
  #
  # ----------------------------------------------------------------------------
  
  
  time_grid <- time_grid_est
  stratum_subjects <- 1:100
  
   
  # Estimate subject-specific covariance functions G_{ij}^{(y_d),k}
  
  n_stratum <- length(stratum_subjects)   # n_y_d
  n_time <- length(time_grid)            # m
  
  cat("Subject-specific covariance input dimensions:\n")
  cat("  n_stratum:", n_stratum, "\n")
  cat("  p (processes):", p, "\n")
  cat("  m (time grid):", n_time, "\n")
  
  # Output array: n_stratum x p x p x m x m
  G_subject_specific <- array(0, dim=c(n_stratum, p, p, n_time, n_time))
  
  cat("Initialized G_subject_specific:", paste(dim(G_subject_specific), collapse=" x "), "\n")
  cat("Memory usage (approx):", prod(dim(G_subject_specific)) * 8 / 1e6, "MB\n")
  
  for (k_idx in 1:n_stratum) {
    k <- stratum_subjects[k_idx]
    
    if (k_idx %% 10 == 1) {
      cat("  Processing subject", k_idx, "/", n_stratum, "(ID:", k, ")\n")
    }
    
    # Get subject-specific intensities
    rho_k <- matrix(0, nrow=p, ncol=n_time)                    # p x m
    rho_k_pairs <- array(0, dim=c(p, p, n_time, n_time))      # p x p x m x m
    
    # Compute subject-specific first-order intensities
    for (i in 1:p) {
      for (t_idx in 1:n_time) {
        t <- time_grid[t_idx]
        event_key <- paste(k, i, sep="_")
        event_times_ki <- event_times[[event_key]]
        if (is.null(event_times_ki)) event_times_ki <- c()
        xi_ki <- length(event_times_ki)  # scalar
        
        if (xi_ki > 0) {
          gamma_ki_t <- estimate_density(event_times_ki, c(t), bandwidth)[1]  # scalar
          rho_k[i, t_idx] <- xi_ki * gamma_ki_t  # scalar
        }
      }
    }
    
    # Compute subject-specific second-order intensities
    for (i in 1:p) {
      for (j in 1:p) {
        for (s_idx in 1:n_time) {
          for (t_idx in 1:n_time) {
            s <- time_grid[s_idx]
            t <- time_grid[t_idx]
            
            event_key_i <- paste(k, i, sep="_")
            event_key_j <- paste(k, j, sep="_")
            event_times_ki <- event_times[[event_key_i]]
            event_times_kj <- event_times[[event_key_j]]
            if (is.null(event_times_ki)) event_times_ki <- c()
            if (is.null(event_times_kj)) event_times_kj <- c()
            xi_ki <- length(event_times_ki)  # scalar
            xi_kj <- length(event_times_kj)  # scalar
            
            if (xi_ki > 0 && xi_kj > 0) {
              gamma_kij <- estimate_bivariate_density(
                event_times_ki, event_times_kj, c(s), c(t), bandwidth
              )  # 1 x 1 matrix
              
              if (i == j) {
                rho_k_pairs[i, j, s_idx, t_idx] <- xi_ki * (xi_ki - 1) * gamma_kij[1, 1]
              } else {
                rho_k_pairs[i, j, s_idx, t_idx] <- xi_ki * xi_kj * gamma_kij[1, 1]
              }
            }
          }
        }
      }
    }
    
    # Compute subject-specific log-covariances
    for (i in 1:p) {
      for (j in 1:p) {
        for (s_idx in 1:n_time) {
          for (t_idx in 1:n_time) {
            numerator <- rho_k_pairs[i, j, s_idx, t_idx]    # scalar
            denominator <- rho_k[i, s_idx] * rho_k[j, t_idx]  # scalar
            
            numerator <- max(numerator, 1e-10)    # scalar
            denominator <- max(denominator, 1e-10)  # scalar
            
            G_subject_specific[k_idx, i, j, s_idx, t_idx] <- log(numerator / denominator)  # scalar
          }
        }
      }
    }
  }
  
  cat("Subject-specific covariance computation completed.\n")
  cat("Output dimensions:", paste(dim(G_subject_specific), collapse=" x "), "\n")
  
  return(G_subject_specific)
}