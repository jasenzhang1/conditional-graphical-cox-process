
library(mvtnorm)

# all functions to generate y_c

# 1) generate continuous and discrete RV's
# 2) generate only continuous RV's because we assume 1 strata
# 3) generate continuous RV's that reflect week; assume 1 strata

generate_conditioning_variables <- function(n, q_c, K, seed = NULL) {
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: Generate mixed continuous-discrete conditioning variables
  #
  # - continuous variables sampled from multivariate normal
  # - discrete variables sampled from multinomial
  #
  # 
  # Input: 
  #
  # - n     (integer)  sample size
  # - q_c   (integer)  dimension of continuous covariates
  # - K     (integer)  number of discrete covariate strata IN TOTAL (ABC12 = 5, not 6)
  # - seed  (integer)  randomization seed
  # 
  # Output: 
  # 
  # - list with Y_continuous (n x q_c) and Y_discrete (n x q_d)
  #
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1) Continuous variables: multivariate normal with moderate correlation -----
  
  # Sigma_c = sigma_c^2 * I + rho_c * 1*1^T  
  sigma_c <- 1.0  # Marginal variance
  rho_c <- 0.3    # Off-diagonal correlation
  Sigma_c <- sigma_c^2 * diag(q_c) + rho_c * matrix(1, q_c, q_c)
  mu_c <- rep(0, q_c)  # Mean vector
  
  Y_continuous <- rmvnorm(n, mean = mu_c, sigma = Sigma_c)
  
  
  # 2) Discrete variables: assume q_d = 2 for simplicity -----------------------
  q_d <- 2
  
  # Split K combinations across q_d variables
  K_per_var <- c(ceiling(sqrt(K)), ceiling(K / ceiling(sqrt(K))))
  
  Y_discrete <- matrix(0, n, q_d)
  
  # Equal probability multinomial
  Y_discrete[, 1] <- sample(1:K_per_var[1], n, replace = TRUE)
  Y_discrete[, 2] <- sample(1:K_per_var[2], n, replace = TRUE)
  
  return(list(
    Y_continuous = Y_continuous,
    Y_discrete = Y_discrete,
    q_c = q_c,
    q_d = q_d,
    parameters = list(
      sigma_c = sigma_c,
      rho_c = rho_c,
      mu_c = mu_c,
      Sigma_c = Sigma_c,
      K_per_var = K_per_var
    )
  ))
}

generate_conditioning_variables_one_strata <- function(n, q_c, seed = NULL) {
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: Generate continuous conditioning variables (assuming one discrete strata)
  #
  # - continuous variables sampled from multivariate normal
  #
  # 
  # Input: 
  #
  # - n     (integer)  sample size
  # - q_c   (integer)  dimension of continuous covariates
  # - seed  (integer)  randomization seed
  # 
  # Output: 
  # 
  # - list with Y_continuous (n x q_c)
  #
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1) Continuous variables: multivariate normal with moderate correlation -----
  
  # Sigma_c = sigma_c^2 * I + rho_c * 1*1^T  
  sigma_c <- 1.0  # Marginal variance
  rho_c <- 0.3    # Off-diagonal correlation
  Sigma_c <- sigma_c^2 * diag(q_c) + rho_c * matrix(1, q_c, q_c)
  mu_c <- rep(0, q_c)  # Mean vector
  
  Y_continuous <- rmvnorm(n, mean = mu_c, sigma = Sigma_c)
  
  

  
  return(list(
    Y_continuous = Y_continuous,
    q_c = q_c,
    parameters = list(
      sigma_c = sigma_c,
      rho_c = rho_c,
      mu_c = mu_c,
      Sigma_c = Sigma_c
    )
  ))
}

generate_conditioning_variables_one_strata_weeks <- function(n, week_vec, seed = NULL) {
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: Generate variables that are similar to weeks (assuming one discrete strata)
  #
  # - continuous variables sampled from multimodal dist
  #
  # 
  # Input: 
  #
  # - n             (integer)  sample size
  # - week_vec      (integer)  vector of possible weeks
  # - seed          (integer)  randomization seed
  # 
  # Output: 
  # 
  # - list with Y_continuous (n x 1)
  #
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1) Continuous variables: multimodal dist 

  
  Y_continuous <-  sample(week_vec, n, replace = TRUE)
  
  
  return(list(
    Y_continuous = Y_continuous,
    q_c = 1
  ))
}

generate_y_c_adj_type <- function(n, adj_type, params, seed = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: depending on the adj_type, generate Y_continuous
  #
  # - for now, assume q_c = 1, so Y_continuous is just a scalar for each subject
  #
  # inputs:
  #
  # - n               (integer)          sample size
  # - adj_type        (string)           precision matrix adjacency type (e.g. single_c2)
  # - params          (vector)           associated vector of parameters of the adjacency type
  # - seed            (integer)          reproducibility seed 
  #
  # 
  # output:
  #
  # - Y_continuous (n x q_c matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # extract last character
  
  y_c_type <- substring(adj_type, nchar(adj_type))
  
  if(y_c_type == '0'){
    
    # params[1] = value
    # all covariates are the same
    
    Y_continuous <- matrix(rep(params[1], n), nrow = n, ncol = 1)
    
    return(Y_continuous)
    
  }   
  
  if(y_c_type == '1'){
    
    # params[1] = min
    # params[2] = max
    # y_c ~ uniform[min, max]
    
    # the premise is that regardless of the value of Y, the underlying graph is the same
    # randomly sample between min and max
    
    Y_continuous <- matrix(runif(n, params[1], params[2]), nrow = n, ncol = 1)
    
    return(Y_continuous)
    
  }
  
  if(y_c_type == '2'){
    
    # params[1] = min
    # params[2] = max
    # equal spacing between all y_c_k's
    
    Y_continuous <- matrix(seq(params[1], params[2], length.out = n), nrow = n, ncol = 1)
    
    return(Y_continuous)
    
  }  
  
  if(y_c_type == '3'){
    
    # params[1] = min
    # params[2] = max
    # params[3] = num_groups
    # mimic weeks. have y_c take on repeated values spread throughtout [min, max]
    m <- params[3]
    values <- seq(params[1], params[2], length.out = m)
    
    counts <- rep(floor(n / m), m)
    counts[1:(n %% m)] <- counts[1:(n %% m)] + 1
    
    # repeat each value accordingly
    Y_continuous <- matrix(rep(values, counts), nrow = n, ncol = 1)
    
  }
  

}  
    
