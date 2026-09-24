


estimate_rho_ij_from_Lambda_v5 <- function(X_mat, i, j){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: with log-intensity values (X_mat) for i_j pair, calculate its weighted bivariate intensity (rho_ij)
  #
  #
  # input:
  #
  # - X_mat         (p x m x n) 
  # - i             (integer) 
  # - j             (integer)
  #
  # outputs:
  #
  # - rho_ij    (m^2 x n) matrix
  #
  # ----------------------------------------------------------------------------
  
  
  Lambda_i <- exp(X_mat[i,,])       # (m x n)
  Lambda_j <- exp(X_mat[j,,])       # (m x n)
  
  m <- dim(Lambda_i)[1]
  
  # 2. Perform column-wise Kronecker product
  # Each column k of the result will be as.vector(L_i[,k] %*% t(L_j[,k]))
  # This results in an (m^2 x n) matrix
  rho_unrolled <- Lambda_i[rep(1:m, times = m), ] * Lambda_j[rep(1:m, each = m), ]  # (m^2 x n)
  
  
  return(rho_unrolled)
  
}
