

# eigendecomposition



solve_sym <- function(A){
  A_hat = solve(A)
  
  return(0.5 * (A_hat + t(A_hat)))
}

sym <- function(A){
  return(0.5 * (A + t(A)))
}



vec <- function(X) as.vector(X)



# -- Main solver -------------------------------------------------------------

