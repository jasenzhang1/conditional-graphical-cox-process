
library(mvtnorm)

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
