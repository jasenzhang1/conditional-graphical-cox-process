

conditional_dependence_function <- function(y_c, beta_ij, type = "linear") {
  
  
  # ----------------------------------------------------------------------------
  #
  # 
  # GOAL: Compute conditional dependence h_ij(y_c)
  #
  # - we offer linear, quadratic, and non-linear dependencies
  # 
  # 
  # Input: 
  #
  # - y_c        (q_c x 1 - dim vector)      vector of continuous covariates
  # - beta_ij    (q_c x 1 - dim vector)      vector of beta coefficients
  # - type       (string)                    either "linear", "quadratic", or "nonlinear"
  #
  # 
  #  Output: 
  #
  # - scalar dependence strength
  #
  #
  # ----------------------------------------------------------------------------
  
  switch(type,
         "linear" = {
           # h_ij(y_c) = 1 + beta_ij^T * y_c
           1 + sum(beta_ij * y_c)
         },
         "quadratic" = {
           # h_ij(y_c) = 1 + beta_ij^T * y_c + 0.1 * ||y_c||^2
           1 + sum(beta_ij * y_c) + 0.1 * sum(y_c^2)
         },
         "nonlinear" = {
           # h_ij(y_c) = 1 + sin(beta_ij^T * y_c) * exp(-||y_c||^2/2)
           linear_term <- sum(beta_ij * y_c)
           1 + sin(linear_term) * exp(-sum(y_c^2)/2)
         },
         1  # Default: no dependence
  )
}


generate_precision_operators <- function(adj_matrix, theta, q_c, 
                                         dependence_type = "linear", seed = NULL) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Generate true precision operators for all discrete strata
  #
  # - here, we just sample beta coefficients uniformly 
  # 
  # Input: 
  #
  # - adj_matrix  (p x p matrix)    matrix of 1's and 0's, with 0 on diags
  # - theta       (number)          signal strength
  # - q_c         (integer)         dimension of continuous covariates
  #
  # 
  # Output: 
  #
  # - specification for precision operators
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  p <- nrow(adj_matrix)
  
  # Generate regression coefficients beta_ij for continuous dependence
  # beta_ij ~ Uniform([-0.5, 0.5]^{q_c}) for each edge (i,j)
  beta_coeffs <- array(0, dim = c(p, p, q_c))
  for (i in 1:p) {
    for (j in 1:p) {
      if (adj_matrix[i, j] == 1) {
        beta_coeffs[i, j, ] <- runif(q_c, -0.5, 0.5)  # why are beta_ij and beta_ji different?
      }
    }
  }
  
  return(list(
    adjacency = adj_matrix,              # E_{y_d}
    beta_coefficients = beta_coeffs,     # {beta_ij}
    signal_strength = theta,             # theta
    dependence_type = dependence_type,   # Type of h_ij(y_c)
    p = p,
    q_c = q_c
  ))
}
