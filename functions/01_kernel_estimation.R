

get_gamma <- function(tseq=NULL, gamma_max = 100){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: custom bandwidth parameter for KDE depending on the data (taken from Sun Lee Li Cai)
  #
  #       - calculate the average distance between all pairs of events, then take
  #         the inverse fourth power
  #
  #
  # input:
  #
  # - tseq             (vector)               event times
  # - gamma_max        (number)               upper cap on gamma, also returned when there are 2 or fewer events
  #
  #
  # output:
  #
  # - gamma            (number)               KDE bandwidth parameter
  #
  # ----------------------------------------------------------------------------
  
  n = length(tseq)
  
  if(n>2){
    
    avg_dist <- sum(dist(tseq)) / (n * (n-1) / 2)
    
    gamma = 1/(avg_dist^4)
  }else{
    return(gamma_max)
  }
  
  if(is.na(gamma)){
    return(gamma_max)
  } else{
    return(min(gamma, gamma_max))
  }
  
}



get_gamma_adaptive <- function(event_times, gamma_max = 100){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: KDE bandwidth parameter from the Sheather-Jones bandwidth, gamma = 1 / (2 * bw^2)
  #
  #
  # input:
  #
  # - event_times      (vector)               event times
  # - gamma_max        (number)               upper cap on gamma
  #
  #
  # output:
  #
  # - gamma            (number)               KDE bandwidth parameter
  #
  # ----------------------------------------------------------------------------
  
  bw <- bw.SJ(event_times)  # or use bw.ucv(event_times)
  gamma <- 1 / (2 * bw^2)
  
  if(gamma > gamma_max){
    gamma <- gamma_max
  }
  
  
  return(gamma)
}

truncNorm_denom <- function(t=NULL, gamma=1 ,a=0, b=1){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: obtain 1 / w_h(t), the boundary correction factor for KDE
  #
  #
  # input:
  #
  # - t                (number)               discretized time
  # - gamma            (number)               bandwidth parameter for KDE
  # - a                (number)               start of time domain
  # - b                (number)               end of time domain
  #
  #
  # output:
  #
  # - denom            (number)               correction factor at t, integral of exp(-gamma (t - x)^2) over [a, b]
  #
  # ----------------------------------------------------------------------------
  
  integrate(function(x) {exp(-gamma*(t-x)^2)}, a, b)$value
}

gaussian_kernel <- function(t_1, t_2, gamma){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: return K(t_1, t_2) for an rbf kernel
  #
  #       K(t_1, t_2) = exp(-gamma ||t_1 - t_2||^2)
  #
  #
  # input:
  #
  # - t_1              (m-dim vector)         vector of times
  # - t_2              (n-dim vector)         vector of times
  # - gamma            (number)               bandwidth parameter
  #
  #
  # output:
  #
  # - kernel_evals     (m x n matrix)         matrix of kernel evaluations
  #
  # ----------------------------------------------------------------------------
  
  pair_diffs <- outer(t_1, t_2, '-')         # m x n matrix of differences
  
  kernel_evals = exp(-gamma * (pair_diffs^2)) # m x n matrix of kernel evals
  
  return(kernel_evals)
}

estimate_density <- function(t_event, t_seq){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: obtain \Gamma_i^k(t) density estimate and \rho_i^k(t) intensity estimate
  #
  #       Checked 6/20/2025
  #
  #       8/5/2025
  #       Need to tune gamma of KDE to be more adaptive to density. KDE is too smooth.
  #
  #       1/12/2026
  #       adapting for mice data, where a given process doesn't necessarily have data from all subjects
  #
  #
  # input:
  #
  # - t_event          (vector)               all timestamps for mark i, subject k
  # - t_seq            (m-dim vector)         evenly spaced out times between 0 and 1
  #
  #
  # output:
  #
  # - output           (list)                 intensity/density values for subject k (an m-dim zero vector if t_event is empty)
  #   - gamma_hat      (m-dim vector)         density estimates for each t_seq time for subject k
  #   - rho_hat        (m-dim vector)         intensity estimates for each t_seq time for subject k
  #   - denom          (m-dim vector)         boundary correction denominator
  #   - gamma          (number)               KDE gamma parameter, obtained adaptively
  #
  # ----------------------------------------------------------------------------
  
  if (length(t_event) == 0) { # moot case when the density is zero
    return(rep(0, length(t_seq))) 
  }  
  
  # 1) Calculate density estimate \Lambda_i^k 
  
  # gamma <- get_gamma(t_event)
  # gamma <- get_gamma_silverman(t_event)
  # gamma <- get_gamma_silverman_v2(t_event)
  gamma <- get_gamma_adaptive(t_event)    # currently the best gamma method, but need a limit of 100
  
  
  kernel_evals <- gaussian_kernel(t_seq, t_event, gamma) # 19 x 3805 matrix
  
  # w_h(t) denominator for each time in Tseq
  denom = sapply(t_seq, function(x){
    truncNorm_denom(x,gamma,0,1) 
  })
  

  
  gamma_hat = apply(kernel_evals,1,mean)/denom   # 19-dim vec / 19-dim vec
  names(gamma_hat) = paste("V",1:length(t_seq),sep="")
  
  rho_hat = apply(kernel_evals,1,sum)/denom   # 19-dim vec / 19-dim vec
  names(rho_hat) = paste("V",1:length(t_seq),sep="")  
  
  
  
  return(list(gamma_hat = gamma_hat, rho_hat = rho_hat, denom = denom, gamma = gamma)) 
  
}


estimate_bivariate_density <- function(event_times_i, event_times_j, 
                                       eval_grid_s, eval_grid_t, d_or_i) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate \Gamma_{i, j}^k(s,t) for given processes i and j at replicate k
  #
  #       checked 6/20 (very sure)
  #
  #
  # input:
  #
  # - event_times_i    (n_i-dim vector)       event times of process i
  # - event_times_j    (n_j-dim vector)       event times of process j
  # - eval_grid_s      (m_s-dim vector)       evaluation times s
  # - eval_grid_t      (m_t-dim vector)       evaluation times t
  # - d_or_i           (string)               'd' for density or 'i' for intensity
  #
  #
  # output:
  #
  # - density2         (m_s x m_t matrix)     bivariate density/intensity estimate
  #
  # ----------------------------------------------------------------------------
  
  # Estimate bivariate density function
  
  xi_i <- length(event_times_i)
  xi_j <- length(event_times_j)
  
  gamma_i <- get_gamma(event_times_i)
  gamma_j <- get_gamma(event_times_j)
  
  if (xi_i == 0 || xi_j == 0) { # moot case when the density is zero
    return(matrix(0, nrow=length(eval_grid_s), ncol=length(eval_grid_t))) # ms x mt matrix
  }
  
  if (identical(event_times_i, event_times_j)) { # second moot case when i = j, and there is only one event.
    if(xi_i == 1){
      return(matrix(0, nrow=length(eval_grid_s), ncol=length(eval_grid_t))) # ms x mt matrix
    }
  }
  
  density <- matrix(0, nrow=length(eval_grid_s), ncol=length(eval_grid_t))
  intensity <- matrix(0, nrow=length(eval_grid_s), ncol=length(eval_grid_t))
  
  # correction factors
  
  ws <- sapply(eval_grid_s, function(x){
    truncNorm_denom(x,gamma_i,0,1)
  })
  
  wt <- sapply(eval_grid_t, function(x){
    truncNorm_denom(x,gamma_j,0,1)
  })  
  
  ws_wt <- outer(ws, wt, '*')   # matrix of all m_s x m_t combinations

  # 2) iterate over all (s, t) combos
  
  for (i in 1:length(eval_grid_s)) {
    for (j in 1:length(eval_grid_t)) {  # for each point (s, t) in the grid
      
      s <- eval_grid_s[i]
      t <- eval_grid_t[j]
      
      ks_vec <- gaussian_kernel(s, event_times_i, gamma_i)    # [1 x n_i] mat
      kt_vec <- gaussian_kernel(t, event_times_j, gamma_j)    # [1 x n_j] mat
      
      #kernel_prod_mat <- outer(as.numeric(ks_vec), as.numeric(kt_vec), '*')  # n_i x n_j
      kernel_prod_mat <- t(ks_vec) %*% kt_vec # n_i x n_j
      
      if (identical(event_times_i, event_times_j)) {  # Same process (diagonal)
        kernel_sum <- sum(kernel_prod_mat) - sum(diag(kernel_prod_mat))
        kernel_mean <- kernel_sum / (xi_i * (xi_i - 1))
      } else{ # different process 
        kernel_sum <- sum(kernel_prod_mat)
        kernel_mean <- kernel_sum / (xi_i * xi_j)
      }
      
      
      density[i, j] <- kernel_mean
      intensity[i, j] <- kernel_sum
    }
  }
  
  density2 <- density/ws_wt 
  intensity2 <- intensity/ws_wt 
  
  if(d_or_i == 'd'){
    return(density2)
  } 
  if(d_or_i == 'i'){
    return(intensity2)
  }
  
  return(NULL)
}
