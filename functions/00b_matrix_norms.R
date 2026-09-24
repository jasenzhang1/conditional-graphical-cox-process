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
  
  # sum of squares --> divide by amount of elements --> then take sqrt
  
  sqrt(sum(A^2) / length(A)) 
  
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
