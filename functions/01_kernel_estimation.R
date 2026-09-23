

get_gamma <- function(tseq=NULL, gamma_max = 100){
  
  # we want a custom bandwidth parameter in KDE depending on the data
  # Taken from Sun Lee Li Cai
  #
  # calculate the average time between all pairs of events. 
  # then take the inverse square
  #
  # input:
  # 
  # - tseq   (vector of event times)
  #
  # output:
  # - gamma 
  
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

get_gamma_silverman <- function(event_times){
  
  n <- length(event_times)
  
  bw <- 1.06 * sd(event_times) * n^(-1/5)
  gamma <- 1 / (2 * bw^2)  
  
  return(gamma)
}

# https://stats.stackexchange.com/questions/6670/which-is-the-formula-from-silverman-to-calculate-the-bandwidth-in-a-kernel-densi
get_gamma_silverman_v2 <- function(event_times){
  
  m <- min(sd(event_times), IQR(event_times) / 1.349)
  
  
  n <- length(event_times)
  
  bw <- 0.9 * m * n^(-1/5)

  gamma <- 1 / (2 * bw^2)
  return(gamma)
}

get_gamma_adaptive <- function(event_times, gamma_max = 100){
  bw <- bw.SJ(event_times)  # or use bw.ucv(event_times)
  gamma <- 1 / (2 * bw^2)
  
  if(gamma > gamma_max){
    gamma <- gamma_max
  }
  
  
  return(gamma)
}

truncNorm_denom <- function(t=NULL, gamma=1 ,a=0, b=1){
  
  # Obtain 1 / w_h(t), the boundary correction factor for KDE
  #
  # input:
  # 
  # - t       (vector of times)   discretized times 
  # - gamma   (number)            bandwidth parameter for KDE
  # - a       (number)            start of time domain
  # - b       (number)            end of time domain
  # 
  # output:
  #
  # - 1 / w_h(t)  (vector of times)   each element corresponds to the correction factor at t_i
  # 
  
  integrate(function(x) {exp(-gamma*(t-x)^2)}, a, b)$value
}

gaussian_kernel <- function(t_1, t_2, gamma){
  
  # ------------------------------------------------------------------------
  # GOAL: return K(t_1, t_2) for an rbf kernel
  #
  #        K(t_1, t_2) = exp(-gamma ||t_1 - t_2||^2)
  #
  # 
  #
  # Input:
  #
  # - t_1             (m-dim vector)   vector of times      
  # - t_2             (n-dim vector)   vector of times
  #
  # output:
  #
  # - kernel_evals    (m x n matrix)   matrix of kernel evalulations
  # -----------------------------------------------
  
  pair_diffs <- outer(t_1, t_2, '-')         # m x n matrix of differences
  
  kernel_evals = exp(-gamma * (pair_diffs^2)) # m x n matrix of kernel evals
  
  return(kernel_evals)
}

estimate_density <- function(t_event, t_seq){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: 
  #
  # obtain \Gamma_i^k(t) density estimate and \rho_i^k(t) intensity estimate
  #
  # Checked 6/20/2025
  #
  # 8/5/2025
  # Need to tune gamma of KDE to be more adaptive to density. KDE is too smooth.
  #
  #
  # 1/12/2026
  # adapting for mice data, where a given process doesn't necessarily have data from all subjects
  #
  #
  # Input:
  #
  # t_event  (3805-dim vector)         all timestamps for mark i, subject k 
  # t_seq    (m-dim vector)           evenly spaced out times between  0 and 1
  #
  # 
  # Output:
  #
  # - list of intensity/density values for subject k 
  #
  #   - gamma_hat     (m-dim vector)           density estimates for each Tseq time for subject k 
  #   - rho_hat       (m-dim vector)           intensity estimates for each Tseq time for subject k
  #   - denom         (m-dim vector)           denominator 
  #   - gamma         (number)                 KDE gamma parameter, obtained adaptively
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

estimate_rho_y_c <- function(t_event, t_seq, y_c_query, y_c_k, gamma_c){
  
  # -----------------------------------------------
  # GOAL: 
  #
  # obtain \rho_i^{(y_c), k} 
  # 
  # - for subject k
  # - incorporating y_c_k and queried y_c
  #
  #
  # Input:
  #
  # - t_event     (vector)                 all timestamps for mark i, subject k 
  # - t_seq       (m-dim vector)           evenly spaced out times between  0 and 1
  # - y_c_query   (q_c-dim vector)         queried vector
  # - y_c_k       (q_c-dim vector)         subject k vector 
  # - gamma_c     (scalar)                 K_c bandwidth
  #
  # 
  # Output:
  #
  # - rho_hat       (m-dim vector)           intensity estimates for each Tseq time for subject k
  # 
  # -----------------------------------------------
  
  if (length(t_event) == 0) { # moot case when the density is zero
    return(rep(0, length(t_seq))) 
  }  
  
  # 1) Calculate density estimate \Lambda_i^k 
  
  # gamma <- get_gamma(t_event)
  # gamma <- get_gamma_silverman(t_event)
  # gamma <- get_gamma_silverman_v2(t_event)
  gamma <- get_gamma_adaptive(t_event)    # currently the best gamma method 
  
  
  kernel_evals <- gaussian_kernel(t_seq, t_event, gamma) # 19 x 3805 matrix
  
  # w_h(t) denominator for each time in Tseq
  denom = sapply(t_seq, function(x){
    truncNorm_denom(x,gamma,0,1) 
  })
  
  
  
  rho_hat <- apply(kernel_evals,1,sum)/denom   # 19-dim vec
  
  
  return(rho_hat) 
  
}

estimate_bivariate_density <- function(event_times_i, event_times_j, 
                                       eval_grid_s, eval_grid_t, d_or_i) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: estimate \Gamma_{i, j}^k(s,t) for given processes i and j at replicate k
  #
  # cheked 6/20 (very sure)
  #
  # input:
  #
  # - event_times_i (n_i vector of times)
  # - event_times_j (n_j vector of times)
  # - eval_grid_s   (m_s-dimensional vector of realized times)
  # - eval_grid_t   (m_t-dimensional vector of realized times)
  # - d_or_i        (string):  density or intensity
  # 
  # output:
  # 
  # - density2   (m_s x m_t matrix) : bivariate density/intensity estimate
  #
  # ------------------------------------------------------------------------
  
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
