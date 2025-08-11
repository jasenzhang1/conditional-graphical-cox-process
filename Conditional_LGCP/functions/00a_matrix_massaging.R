adjust_R <- function(Rmat=NULL){
  
  # -------------------------------------------------------------------------
  # 
  # goal: when we truncate Rmat because we don't need all those eigenfunctions,
  #       we need to adjust the shrunken matrix 
  #  
  #
  # input: 
  # - Rmat (80 x 80 matrix), that used to be (120 x 120) but we truncated
  #
  # -------------------------------------------------------------------------
  
  Rmat = cov2cor(Rmat)
  Rmat[Rmat >= 0.99] = 0.99
  Rmat[Rmat <= -0.99] = -0.99
  diag(Rmat) = 1
  
  min_eig_val = min(find_min_eigen(c(list(Rmat)))) # find min eigenvalue
  
  # remember that covariance matrices must have nonnegative eigenvalues
  if(min_eig_val < pinv_eps){
    min_eig_val = abs(min_eig_val) + pinv_eps*(1+abs(min_eig_val))/0.99
    Rmat = adjust_S(Rmat,min_eig_val)
  }
  Rmat
}

adjust_S <- function(S=NULL, val=NULL){
  
  # remember that adding c * identity increases all the eigenvalues by c
  diag(S) = diag(S)+val
  cov2cor(S)
}


find_min_eigen <- function(S=NULL){
  
  # ----------------------------------
  #
  # input is a list of matrices
  #
  # return a list of min eigenvalues for each matrix
  #
  # ---------------------------
  eig = sapply(S, function(s){
    min(eigen(s)$val)
  })
  eig
}

eigendecomp_sampling <- function(precision_matrix){
  
  # eigendecomposition, and retain only the components that come from positive eigenvalues
  
  eig <- eigen(precision_matrix, symmetric = TRUE)
  # Keep only positive eigenvalues
  positive <- eig$values > 1e-8
  values_inv_sqrt <- 1 / sqrt(eig$values[positive])
  basis <- eig$vectors[, positive]
  z <- rnorm(sum(positive))
  x <- basis %*% (values_inv_sqrt * z)
}


psd_jitter <- function(A, pinv_eps = 1e-6){
  
  #
  # MAIN MASSAGING TECHNIQUE DURING ESTIMATION
  #
  #
  # A is a symmetric, matrix that we want to make positive definite
  # pinv_eps = value to make our smallest eigenvalue
  # 
  
  eigvals <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(eigvals)

  
  if (lambda_min <= 0) {
    epsilon <- -lambda_min + pinv_eps
    A_fixed <- A + epsilon * diag(nrow(A))
  } else {
    A_fixed <- A  # already PD
  }  
  
  return(A_fixed)
}


