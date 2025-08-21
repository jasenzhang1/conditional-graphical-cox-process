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

solve_sym <- function(A){
  A_hat = solve(A)
  
  return(0.5 * (A_hat + t(A_hat)))
}

psd_jitter <- function(A, pinv_eps = 1e-6){
  
  #
  # MAIN MASSAGING TECHNIQUE DURING ESTIMATION
  #
  #
  # A is a symmetric, matrix that we want to make positive definite
  # pinv_eps = value to make our smallest eigenvalue
  # 
  
  A <- (t(A) + A) / 2 # guarantee symmetry
  
  eigvals <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(eigvals)

  
  if (lambda_min <= 0) {
    epsilon <- -lambda_min + pinv_eps
    A_fixed <- A + epsilon * diag(nrow(A))
    
    A_fixed <- (t(A_fixed) + A_fixed) / 2 # guarantee symmetry again
    
  } else {
    A_fixed <- A  # already PD
  }  
  
  return(A_fixed)
}

kronecker_decomp <- function(A, p, m) {
  
  # we want to find the "Least squares" kronecker product decomposition 
  #
  # input:
  # 
  # - A  (pm x pm matrix)
  # - p  (value)
  # - m  (value)
  
  stopifnot(nrow(A) == p*m, ncol(A) == p*m)
  
  # Rearrange
  R <- matrix(0, p^2, m^2)
  for (i in 1:p) {
    for (j in 1:p) {
      block <- A[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)]
      R[(i-1)*p + j, ] <- as.vector(block)
    }
  }
  
  # SVD
  svd_res <- svd(R)
  
  # Extract X and Y
  sigma1 <- svd_res$d[1]
  u1 <- svd_res$u[,1] * sqrt(sigma1)
  v1 <- svd_res$v[,1] * sqrt(sigma1)
  
  X_hat <- matrix(u1, p, p)
  Y_hat <- matrix(v1, m, m)
  
  list(p_mat = X_hat, m_mat = Y_hat)
}


kronecker_decomp_diag1 <- function(A, p, m, maxit = 200, tol = 1e-10, restarts = 3, seed = NULL) {
  
  # 
  #
  # GOAL: want to estimate the kronecker decomposition of a pm x pm matrix
  #       but we also want to condition the mxm matrix to have 1's on its diagonal
  #
  #
  # 
  # input:
  # - A    (pm x pm matrix)
  # - p    (number)
  # - m    (number)
  # - maxit 
  # - tol
  # - restarts
  # - seed 
  #
  # 
  
  stopifnot(is.matrix(A), nrow(A) == p*m, ncol(A) == p*m)
  if (!is.null(seed)) set.seed(seed)
  
  # Build rearranged matrix R \in R^{p^2 x m^2}
  R <- matrix(0, p^2, m^2)
  for (i in 1:p) for (j in 1:p) {
    block <- A[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)]
    R[(i-1)*p + j, ] <- as.vector(block)
  }
  
  diag_idx <- (1:m) + (0:(m-1))*m  # indices in vec(Y) for the diagonal (column-major)
  
  best <- NULL
  best_err <- Inf
  
  one_run <- function(y_init) {
    y <- y_init
    # enforce diag=1 at start
    y[diag_idx] <- 1
    
    prev_err <- Inf
    for (it in 1:maxit) {
      # x update
      denom_y <- sum(y*y)
      if (denom_y == 0) break
      x <- drop(R %*% y) / denom_y
      
      # y free update
      denom_x <- sum(x*x)
      if (denom_x == 0) break
      y_free <- drop(t(R) %*% x) / denom_x
      
      # project onto diag(Y)=1: simply overwrite those entries
      y <- y_free
      y[diag_idx] <- 1
      
      # objective value
      # ||R - x y^T||_F^2 = ||R||_F^2 - 2 x^T R y + ||x||^2 ||y||^2
      Ry <- drop(R %*% y)
      err <- sum(R*R) - 2 * sum(x * Ry) + sum(x*x) * sum(y*y)
      
      if (abs(prev_err - err) < tol * (1 + prev_err)) break
      prev_err <- err
    }
    
    list(x = x, y = y, err = prev_err, iters = it)
  }
  
  # Restarts: identity, randoms
  inits <- list(as.vector(diag(m)))
  inits <- c(inits, replicate(restarts, rnorm(m*m), simplify = FALSE))
  for (y0 in inits) {
    run <- one_run(y0)
    if (run$err < best_err) {
      best_err <- run$err
      best <- run
    }
  }
  
  X_hat <- matrix(best$x, p, p)
  Y_hat <- matrix(best$y, m, m)  # will have diag=1 by construction
  list(p_mat = X_hat, m_mat = Y_hat, frob_error = best$err, iterations = best$iters)
}

vec <- function(X) as.vector(X)
unvec <- function(v, nrow, ncol) matrix(v, nrow = nrow, ncol = ncol)

# Rearrangement operator R(M): (p^2 x m^2) matrix whose r=(i-1)p+j row is vec( M_ij )^T
rearrange_pm <- function(M, p, m) {
  stopifnot(is.matrix(M), nrow(M) == p*m, ncol(M) == p*m)
  R <- matrix(0, nrow = p*p, ncol = m*m)
  r <- 1
  for (i in 1:p) {
    rows <- ((i - 1) * m + 1):(i * m)
    for (j in 1:p) {
      cols <- ((j - 1) * m + 1):(j * m)
      Rij <- as.vector(M[rows, cols])
      R[r, ] <- Rij
      r <- r + 1
    }
  }
  R
}

# Symmetrize a square matrix
symmetrize <- function(A) 0.5 * (A + t(A))

# PSD projection via eigenvalue clipping
project_psd <- function(A, jitter = 0) {
  A <- symmetrize(A)
  ev <- eigen(A, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  Apsd <- ev$vectors %*% (lam * t(ev$vectors))
  if (jitter > 0) Apsd <- Apsd + diag(jitter, nrow(Apsd))
  Apsd
}

# -- Main solver -------------------------------------------------------------

# Minimize || M - B ⊗ C ||_F with C ≽ 0 (covariance-viable), using ALS.
# Scale is non-unique; we enforce ||C||_F = 1 each iteration.
kronecker_psd_factor <- function(M, p, m, maxit = 200, tol = 1e-8, 
                                 init = c("svd", "random"),
                                 enforce_trace = FALSE) {
  stopifnot(is.matrix(M), nrow(M) == p*m, ncol(M) == p*m)
  init <- match.arg(init)
  
  # Rearranged matrix
  Rmat <- rearrange_pm(M, p, m)
  
  # Initialize via rank-1 SVD of R(M)
  if (init == "svd") {
    sv <- svd(Rmat, nu = 1, nv = 1)
    b <- sv$u[, 1] * sqrt(sv$d[1])
    c <- sv$v[, 1] * sqrt(sv$d[1])
  } else {
    b <- rnorm(p*p)
    c <- rnorm(m*m)
  }
  
  # Normalize C factor (to fix scale)
  C <- project_psd(unvec(c, m, m))
  Cf <- norm(C, type = "F")
  if (Cf == 0) Cf <- 1
  C <- C / Cf
  c <- vec(C)
  
  obj_prev <- Inf
  history <- numeric(0)
  
  for (it in 1:maxit) {
    # --- Update B given C (closed-form least squares on rearranged system)
    denom_c <- sum(c * c)
    if (denom_c == 0) denom_c <- 1e-16
    b <- as.vector(Rmat %*% c) / denom_c
    B <- unvec(b, p, p)
    
    # Optional: symmetrize B (often desired but not strictly required)
    B <- symmetrize(B)
    b <- vec(B)
    
    # --- Update C given B
    denom_b <- sum(b * b)
    if (denom_b == 0) denom_b <- 1e-16
    c <- as.vector(t(Rmat) %*% b) / denom_b
    C <- unvec(c, m, m)
    
    # Enforce covariance-viability: symmetric PSD; optional trace normalization
    C <- project_psd(C)
    if (enforce_trace) {
      tr <- sum(diag(C))
      if (tr > 0) C <- C * (m / tr)  # average variance ~ 1
    } else {
      # Default: Frobenius normalization
      Cf <- norm(C, type = "F")
      if (Cf > 0) C <- C / Cf
    }
    c <- vec(C)
    
    # --- Objective value ||M - B ⊗ C||_F
    Mhat <- kronecker(B, C)
    res <- M - Mhat
    obj <- sum(res * res)
    history <- c(history, obj)
    
    if (abs(obj_prev - obj) <= tol * max(1, obj_prev)) break
    obj_prev <- obj
  }
  
  list(
    B = B,
    C = C,             # guaranteed symmetric PSD
    objective = obj,
    iters = length(history),
    history = history
  )
}
