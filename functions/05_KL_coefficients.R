

estimate_log_intensity_function <- function(event_times, t_seq) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Obtain X_ik(t) estimate
  #
  # Input: 
  #
  # - event_times    (vector of length xi_k_i)  subject k, process i 
  # - t_seq          (vector of length m)
  #
  # 
  # Output: 
  #
  # - X_hat           (vector of length m)    subject k, process i
  #
  #
  # ----------------------------------------------------------------------------
  
  if (length(event_times) == 0) {
    return(rep(-10.0, length(t_seq)))  # Large negative value
  }
  
  # Estimate intensity: returns m x 1 vector
  lambda_hat <- estimate_density(event_times, t_seq)$rho_hat
  
  # Take logarithm: (m x 1) -> (m x 1)
  X_hat <- log(pmax(lambda_hat, 1e-10))
  
  return(X_hat)
}



estimate_kl_coefficients_parallel_v2 <- function(X_k_est, eigenfunctions, 
                                                 t_seq, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Obtain KL coefficients for subject k, process i
  #
  # - use pbmclapply to parallelize for each subject + process
  #
  # - 8/11/2025 faster, assumes we already have log-intensities
  # - 8/12/2025 I DON'T NEED TO MULTIPLY BY DT HUH
  #
  #
  # Input: 
  # 
  # - X_k_est           (p x m x n matrix)
  # - eigenfunctions    (list of p matrices of dimension m x d_i) 
  # - t_seq             (vector of length m)
  # - ncores            
  #
  # Output: 
  #
  # - alpha_tensor (n_stratum x p x d array)
  #
  # ----------------------------------------------------------------------------
  
  if(length(t_seq) != dim(eigenfunctions[[1]])[1]){
    print('tseq and eigenfunction dimensions are incompatible')
    return(NULL)
  }
  
  p <- dim(X_k_est)[1]
  m <- dim(X_k_est)[2]
  n <- dim(X_k_est)[3]
  
  dt <- if (length(t_seq) > 1) t_seq[2] - t_seq[1] else 1.0
  
  max_components <- max(sapply(eigenfunctions, ncol) %>% unlist())
  
  
  
  # list of n entries, each is a (p x d) matrix
  alpha_list <- pbmclapply(1:n, function(k_idx) {
    
    alpha_slice <- array(0, dim = c(p, max_components)) # p x d
    
    for (i in 1:p) {
      
      # Estimate log-intensity function: returns m x 1 vector
      X_hat_ki <- X_k_est[i, , k_idx]
      
      # Compute KL coefficients for all components
      
      eta_i <- eigenfunctions[[i]]  # m x d_i matrix
      
      # Numerical integration: scalar result
      integrand <- X_hat_ki * eta_i      # (m x 1 vector) * (m x d_i matrix) = (m x d_i matrix)
      
      # pad some zeros if needed
      integrand_sum <- colSums(integrand) # d_i-dim vector
      integrand_pad <- c(integrand_sum, rep(0, max_components - length(integrand_sum))) # d-dim vector
      alpha_slice[i, ] <- integrand_pad   # NO DT 
      
    }
    
    alpha_slice
  }, mc.cores = ncores)
  
  # KL coefficients: n_stratum x p x d array
  alpha_tensor <- abind(alpha_list, along = 0)
  
  return(alpha_tensor)
}


