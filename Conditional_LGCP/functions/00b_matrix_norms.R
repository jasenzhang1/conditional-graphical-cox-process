suppressPackageStartupMessages(library(akima))

hilbert_schmidt_norm <- function(A) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: Compute Hilbert-Schmidt norm, which is just frobenius norm 
  #
  # - ||A||_{HS} = sqrt(sum of squares of elements)
  # 
  # 
  # Input: 
  #
  # - A       (m x m matrix)
  #
  # 
  # Output: 
  #
  # - HS norm (scalar)
  #
  # ------------------------------------------------------------------------
  
  return(sqrt(sum(A^2)))  # sqrt(sum of all squared elements)
}

hilbert_schmidt_norm_rmse <- function(A){
  
  sqrt(sum(A^2)) / sqrt(length(A)) 
  
}

hilbert_schmidt_norm_normalize <- function(A, delta_t) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: Compute Hilbert-Schmidt norm, which is just frobenius norm 
  #
  # - ||A||_{HS} = sqrt(sum of squares of elements)
  # 
  # 
  # Input: 
  #
  # - A       (m x m matrix)
  #
  # 
  # Output: 
  #
  # - HS norm (scalar)
  #
  # ------------------------------------------------------------------------
  
  return(delta_t * sqrt(sum(A^2)))  # sqrt(sum of all squared elements)
}

hilbert_schmidt_norm_pm <- function(A, p, m) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: Compute Hilbert-Schmidt norm of all p^2 mxm block matrices
  #
  # - ||A||_{HS} = sqrt(sum of squares of elements)
  # 
  # 
  # Input: 
  #
  # - A       (pm x pm matrix)
  #
  # 
  # Output: 
  #
  # - norms (p x p matrix)
  #
  # ------------------------------------------------------------------------
  
  
  if (!is.matrix(A) || nrow(A) != p*m || ncol(A) != p*m) {
    stop("matrix input must be a (p*m) x (p*m) matrix")
  }
  
  norms <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      # Row indices for block (i, j)
      row_idx <- ((i - 1) * m + 1):(i * m)
      # Column indices for block (j, j)
      col_idx <- ((j - 1) * m + 1):(j * m)
      block <- A[row_idx, col_idx]
      norms[i, j] <- sqrt(sum(block^2))  # HS (Frobenius) norm
    }
  }
  
  return(norms)
  
}

hilbert_schmidt_norm_pm_rmse <- function(A, p, m) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: Compute Hilbert-Schmidt norm of all p^2 mxm block matrices, then take RMSE
  # 
  #       RMSE is dividing by sqrt(# of entries) which is m
  #
  # - ||A||_{HS} = sqrt(sum of squares of elements)
  # 
  # 
  # Input: 
  #
  # - A       (pm x pm matrix)
  #
  # 
  # Output: 
  #
  # - norms (p x p matrix)
  #
  # ------------------------------------------------------------------------
  
  
  return(hilbert_schmidt_norm_pm(A, p, m) / m)
}

hilbert_schmidt_norm_pm_normalize <- function(A, p, m) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Compute Hilbert-Schmidt norm of all p^2 mxm block matrices
  #
  # - ||A||_{HS} = sqrt(sum of squares of elements)
  # - multiply this by delta_t to normalize
  # 
  # 
  # Input: 
  #
  # - A       (pm x pm matrix)
  #
  # 
  # Output: 
  #
  # - norms (p x p matrix)
  #
  # ----------------------------------------------------------------------------
  
  
  if (!is.matrix(A) || nrow(A) != p*m || ncol(A) != p*m) {
    stop("matrix input must be a (p*m) x (p*m) matrix")
  }
  
  delta_t <- 1/m
  
  norms <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      # Row indices for block (i, j)
      row_idx <- ((i - 1) * m + 1):(i * m)
      # Column indices for block (j, j)
      col_idx <- ((j - 1) * m + 1):(j * m)
      block <- A[row_idx, col_idx]
      norms[i, j] <- sqrt(sum(block^2))  # HS (Frobenius) norm
    }
  }
  
  return(delta_t * norms)
  
}

# (i_j list --> p x p matrix)
hilbert_schmidt_norm_list_to_mat <- function(M_list, p){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Given a list of i_j matrices, compute the Hilbert-Schmidt norm of these blocks and arrange them in a pxp matrix
  #
  #
  # inputs:
  # 
  # - M_list  (list of i_j matrices, where i <= j)
  # - p       (integer)
  # 
  # 
  # Output: 
  #
  # - HS_mat (p x p matrix)
  #
  # ----------------------------------------------------------------------------
  
  # 1) convert the names to an n x 2 matrix
  nm <- names(M_list)
  parts <- strsplit(nm, "_")
  ij_mat <- do.call(rbind, parts)
  ij_mat <- apply(ij_mat, 2, as.integer)
  
  # 2) Compute HS norms
  HS_norms <- sapply(M_list, hilbert_schmidt_norm)
  
  # 3) Initialize and fill in matrix
  HS_mat <- matrix(0, p, p)
  
  for(k in 1:nrow(ij_mat)){
    i <- ij_mat[k,1]
    j <- ij_mat[k,2]
    

    HS_mat[i,j] <- HS_norms[k]
    HS_mat[j,i] <- HS_norms[k]
  }
  
  return(HS_mat)
}

hilbert_schmidt_norm_diff <- function(mat_A, mat_B, time_grid_A, time_grid_B, fine_n = 200){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Calculate the Hilbert Schmidt norm of two matrices with different discretizations
  #
  # - we partition their grids finely so that they are the same dimensions.
  #
  # inputs:
  #
  # - mat_A         (m1 x m1 matrix)   first matrix of values
  # - mat_B         (m2 x m2 matrix)   second matrix of values
  # - time_grid_A   (m1-dim vec)       time discretization wrt first matrix
  # - time_grid_B   (m2-dim vec)       time_discretization wrt second matrix
  # - fine_n        (scalar)           how many timepoints we want to approximate the integral with
  # 
  # outputs:
  #
  # - HS_norm      (scalar)          HS norm of their differences
  # 
  #
  # ----------------------------------------------------------------------------
  
  if (!requireNamespace("akima", quietly = TRUE)) {
    stop("Please install the 'akima' package: install.packages('akima')")
  }
  
  # Create coarse (x,y,z) data for A
  grid_A <- expand.grid(x = time_grid_A, y = time_grid_A)
  z_A <- as.vector(mat_A)
  
  # Create coarse (x,y,z) data for B
  grid_B <- expand.grid(x = time_grid_B, y = time_grid_B)
  z_B <- as.vector(mat_B)
  
  # Define fine grid (common resolution)
  time_fine <- seq(0, 1, length.out = fine_n)
  
  # Bilinear interpolation onto fine grid
  interp_A <- akima::interp(x = grid_A$x, y = grid_A$y, z = z_A,
                            xo = time_fine, yo = time_fine, linear = TRUE, extrap = TRUE)
  interp_B <- akima::interp(x = grid_B$x, y = grid_B$y, z = z_B,
                            xo = time_fine, yo = time_fine, linear = TRUE, extrap = TRUE)
  
  # Approximate integral via Riemann sum
  dx <- diff(time_fine)[1]
  diff_sq <- (interp_A$z - interp_B$z)^2
  HS_norm <- sqrt(sum(diff_sq, na.rm = TRUE) * dx^2)
  
  return(HS_norm)
  
}