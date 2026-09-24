# estimate G_{ij}(s, t)



# data reformat function


estimate_covariance_functions_ii <- function(rho_i, rho_ii, regularization=1e-10) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Want to estimate G_{ij}(s, t)  across all p x p pairs of processes
  #
  # - but we only care about G_{i,i}(s,t)
  #
  # 
  # Input: 
  #
  #   - rho_i      (p x m matrix)                    univariate intensity
  #   - rho_ii     (list of p m x m matrices)        bivariate intensity
  #
  #
  # Output: 
  #
  # - G_hat           (list of m x m matrices for i_i process)
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(rho_i)[1]
  m <- dim(rho_i)[2]
  G_hat <- list()
  
  rho_list <- lapply(1:p, function(i) {
    list(
      rho_i  = rho_i[i, ],     
      rho_ii_mat = rho_ii[[i]]  
    )
  })   
  
  for (i in 1:p) {
    key <- paste0(i, '_', i)
    
    numerator <- rho_ii[[i]]
    
    denominator <- outer(rho_i[i,], rho_i[i,])
        
    # Avoid log(0) by adding regularization
    numerator <- pmax(numerator, regularization)
    denominator <- pmax(denominator, regularization)
    
    G_hat[[key]] <- log(numerator / denominator)

    
  }
  
  return(G_hat)
}

estimate_covariance_functions_ij <- function(rho_i, rho_ii_mat, regularization=1e-10) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Want to estimate G_{ij}(s, t)  across all p x p pairs of processes
  #
  # 
  # Input: 
  # 
  #   - rho_i      (p x m matrix)                      univariate intensity
  #   - rho_ii_mat (i_j list of m x m matrices)        bivariate intensity
  #
  #
  # Output: 
  #
  # - G_hat           (i_j list of m x m matrices) 
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(rho_i)[1]

  G_hat <- list()
  
  for (i in 1:p) {
    for (j in i:p){
      key <- paste0(i, '_', j)
    
      # rho_ij
      numerator <- rho_ii_mat[[key]]                             # m x m
      
      # outer(rho_i, rho_j)
      denominator <- outer(rho_i[i,], rho_i[j,])  # m x m
      
      # Avoid log(0) by adding regularization
      numerator <- pmax(numerator, regularization)
      denominator <- pmax(denominator, regularization)
      
      G_hat[[key]] <- log(numerator / denominator)
    }
    
  }
  
  return(G_hat)

}

# section 2: functions to get cross-informed G_ii

# - start with rho's to get G_ii's and G_ij's 


