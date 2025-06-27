library(abind)

estimate_log_intensity_function <- function(event_times, t_seq) {
  
  #
  # GOAL: Obtain X_ik(t) estimate
  #
  # Input: 
  #
  # - event_times    (vector of length xi_k_i)
  # - t_seq          (vector of length m)
  #
  # Output: 
  #
  # - X_hat (vector of length m)
  
  if (length(event_times) == 0) {
    return(rep(-10.0, length(t_seq)))  # Large negative value
  }
  
  # Estimate intensity: returns m x 1 vector
  lambda_hat <- estimate_density(event_times, t_seq)$rho_hat
  
  # Take logarithm: (m x 1) -> (m x 1)
  X_hat <- log(pmax(lambda_hat, 1e-10))
  
  return(X_hat)
}

estimate_kl_coefficients <- function(data_all, eigenfunctions, truncations, patient_sel, feature_sel, 
                                     t_seq) {

  # ----------------------------------------------------------------------------
  #
  # GOAL: Obtain KL coefficients for subject k, process i
  #
  #
  # Input: 
  # 
  # - data_all          (data frame for all processes, all subjects)  'feature_id', 'time', 'subject_num'
  # - eigenfunctions    (list of p matrices of dimension m x d_i) 
  # - truncations       (list of integers denoting how many eigenfunctions are needed)
  # - patient_sel       (vector of patient_ID's in this stratum)
  # - t_seq             (vector of length m)
  # - p                 (number of processes)
  #
  # Output: 
  #
  # - alpha_hat (n_stratum x p x d array)
  #
  # ----------------------------------------------------------------------------
  
  p <- length(feature_sel)
  
  dt <- if (length(t_seq) > 1) t_seq[2] - t_seq[1] else 1.0
  
  n_stratum <- length(patient_sel)
  max_components <- ncol(eigenfunctions[[1]])  # d
  
  max_components <- max(sapply(eigenfunctions, ncol) %>% unlist())
  
  # KL coefficients: n_stratum x p x d array
  alpha_hat <- array(0, dim=c(n_stratum, p, max_components))
  
  for (k_idx in 1:n_stratum) {
    k <- patient_sel[k_idx]
    
    for (i in 1:p) {

      event_times_ki <- data_all %>% filter(feature_id == i & subject_num == k) %>% dplyr::pull(time) # event times of subject k, process i
      if (is.null(event_times_ki)) event_times_ki <- c()
      
      # Estimate log-intensity function: returns m x 1 vector
      X_hat_ki <- estimate_log_intensity_function(event_times_ki, t_seq)
      
      # Compute KL coefficients for all components
      for (a in 1:truncations[[i]]) {
        eta_ia <- eigenfunctions[[i]][, a]  # m x 1 vector
        # Numerical integration: scalar result
        integrand <- X_hat_ki * eta_ia      # (m x 1) .* (m x 1) = (m x 1)
        alpha_hat[k_idx, i, a] <- sum(integrand) * dt  # scalar
      }
    }
  }
  
  return(alpha_hat)
}

estimate_kl_coefficients_parallel <- function(data_all, eigenfunctions, truncations, patient_sel, feature_sel, 
                                     t_seq, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Obtain KL coefficients for subject k, process i
  #
  # - use pbmclapply to parallelize for each subject + process
  #
  #
  # Input: 
  # 
  # - data_all          (data frame for all processes, all subjects)  'feature_id', 'time', 'subject_num'
  # - eigenfunctions    (list of p matrices of dimension m x d_i) 
  # - truncations       (list of integers denoting how many eigenfunctions are needed)
  # - patient_sel       (vector of patient_ID's in this stratum)
  # - t_seq             (vector of length m)
  # - p                 (number of processes)
  #
  # Output: 
  #
  # - alpha_hat (n_stratum x p x d array)
  #
  # ----------------------------------------------------------------------------
  
  p <- length(feature_sel)
  
  dt <- if (length(t_seq) > 1) t_seq[2] - t_seq[1] else 1.0
  
  n_stratum <- length(patient_sel)
  
  max_components <- max(sapply(eigenfunctions, ncol) %>% unlist())
  
  

  
  # list of n entries, each is a (p x d) matrix
  alpha_list <- pbmclapply(1:n_stratum, function(k_idx) {
  
    alpha_slice <- array(0, dim = c(p, max_components)) # p x d
    k <- patient_sel[k_idx]
    
    for (i in 1:p) {
      
      event_times_ki <- data_all %>% filter(feature_id == i & subject_num == k) %>% dplyr::pull(time) # event times of subject k, process i
      if (is.null(event_times_ki)) event_times_ki <- c()
      
      # Estimate log-intensity function: returns m x 1 vector
      X_hat_ki <- estimate_log_intensity_function(event_times_ki, t_seq)
      
      # Compute KL coefficients for all components

      eta_i <- eigenfunctions[[i]]  # m x d_i matrix
      # Numerical integration: scalar result
      integrand <- X_hat_ki * eta_i      # (m x 1 vector) .* (m x d_i matrix) = (m x d_i matrix)
      
      # pad some zeros if needed
      integrand_sum <- colSums(integrand)
      integrand_pad <- c(integrand_sum, rep(0, max_components - dim(eta_i)[2]))
      alpha_slice[i, ] <- integrand_pad * dt  # scalar

    }
    
    alpha_slice
  }, mc.cores = ncores)
  
  # KL coefficients: n_stratum x p x d array
  alpha_tensor <- abind(alpha_list, along = 0)
  
  return(alpha_tensor)
}
