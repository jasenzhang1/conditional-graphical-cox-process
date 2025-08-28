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