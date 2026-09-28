

# eigendecomposition



solve_sym <- function(A){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: invert a matrix and symmetrize the result to remove numerical asymmetry
  #
  #
  # input:
  #
  # - A                (square matrix)        matrix to invert
  #
  #
  # output:
  #
  # - A_inv            (square matrix)        0.5 * (A^{-1} + t(A^{-1}))
  #
  # ----------------------------------------------------------------------------
  
  A_hat = solve(A)
  
  return(0.5 * (A_hat + t(A_hat)))
}

sym <- function(A){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: symmetrize a matrix to remove numerical asymmetry
  #
  #
  # input:
  #
  # - A                (square matrix)        matrix to symmetrize
  #
  #
  # output:
  #
  # - A_sym            (square matrix)        0.5 * (A + t(A))
  #
  # ----------------------------------------------------------------------------
  
  return(0.5 * (A + t(A)))
}



vec <- function(X){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: vectorize a matrix by stacking its columns
  #
  #
  # input:
  #
  # - X                (m x n matrix)         matrix to vectorize
  #
  #
  # output:
  #
  # - vec_X            (mn-dim vector)        columns of X stacked on top of each other
  #
  # ----------------------------------------------------------------------------
  
  return(as.vector(X))
}



# -- Main solver -------------------------------------------------------------

