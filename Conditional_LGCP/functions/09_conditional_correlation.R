estimate_conditional_correlation_v4 <- function(V_cond_mat, p, identity = T, pinv_eps = 1e-6) {
  
  # ---------------------------------------------------------------------
  #
  # GOAL: estimate conditional correlation
  #
  # - v3: we adaptively estimate gamma depending on the negative eigenvalues
  # -     use psd_jitter in 00a_matrix_massaging
  #
  # - v4: 8/14/2025
  # -     we now let the i = j term be calculated like the rest 
  #
  # Input: 
  #
  # - V_cond_mat         (pm x pm matrix)
  # - p                  (scalar)
  # - identity           (boolean)  is i_i set to be the identity?
  # - pinv_eps           (number)   the psd_jitter constant
  #
  #
  # Output: 
  #
  # - C_conditional      (list of length p^2, each element m x m matrix)
  #
  #
  #------------------------------------------------------------------------
  
  m <- dim(V_cond_mat)[1] / p
  
  V_conditional <- extract_block_structure_v2(V_cond_mat, p, m)
  
  C_conditional <- list()
  
  for (i in 1:p) {
    for (j in i:p) {
      key <- paste(i, j, sep="_")
      

      key_ii <- paste(i, i, sep="_")
      key_jj <- paste(j, j, sep="_")
      
      V_ii <- psd_jitter(V_conditional[[key_ii]], pinv_eps = pinv_eps)  # m x m
      V_jj <- psd_jitter(V_conditional[[key_jj]], pinv_eps = pinv_eps)  # m x m  
      V_ij <- V_conditional[[key]]                 # m x m
      
      # Matrix operations: (m x m) %*% (m x m) %*% (m x m) = (m x m)
      V_ii_inv_sqrt <- matrix_inv_sqrt(V_ii)  # m x m
      V_jj_inv_sqrt <- matrix_inv_sqrt(V_jj)  # m x m
      
      C_conditional[[key]] <- V_ii_inv_sqrt %*% V_ij %*% V_jj_inv_sqrt

    }
  }
  
  return(C_conditional)
}



correlation_estimation_KL_cov <- function(eigendecomp, KL_cov, identity = T){
  
  # ----------------------------------------------------------------------------
  # 
  #
  # GOAL: construct the correlation operator from the KL covariance method
  #
  # input:
  #
  # - eigendecomp (list of 3 entries)
  #
  #   - eigenvalues
  #   - eigenfunctions
  #   - n_dims
  #
  # - KL_cov (list of i_j entries, each is d x d matrix) each value represents covariance between KL coeffs of components a and b in process i and j
  # - identity    (boolean)  are the i = j entries just the identity?
  #
  # output:
  #
  # - C_cond   (list of i_j entries)
  #
  #
  # ----------------------------------------------------------------------------
  
  
  p <- length(eigendecomp[[1]])
  d <- dim(KL_cov[[1]])[1]
  m <- dim(eigendecomp[[2]][[1]])[1]
  Delta <- 1/m
  
  C_cond <- list()
  
  for(i in 1:p){
    for(j in i:p){
      
      key <- paste0(i, '_', j)
      
      if(i == j & identity){
        C_ij <- diag(m)
      } else{
        
      
        
        
        eval_i <- eigendecomp[[1]][[i]]
        eval_j <- eigendecomp[[1]][[j]]
        
        evec_i <- eigendecomp[[2]][[i]]
        evec_j <- eigendecomp[[2]][[j]]
        
        d_i <- dim(evec_i)[2]
        d_j <- dim(evec_j)[2]
        
        cov_ij <- KL_cov[[key]]
        
        
        # cov / sqrt(var_i * var_j)
        coeffs <- cov_ij / sqrt(tcrossprod(eval_i, eval_j))
        
        if (any(is.nan(coeffs))) {
          stop("correlation construction process contains NaN values!")
        }      
        
        
        C_ij <- matrix(0, nrow = m, ncol = m)
        
        for(a in 1:d_i){
          for(b in 1:d_j){
            
            # REMEMBER, we could have eigenvectors defined as Delta * (v^\top v) = 1. 
            
            C_ij <- C_ij + coeffs[a,b] * tcrossprod(evec_i[, a], evec_j[, b])
            
          }
        }
      } # end of nonidentity case 
      
      C_cond[[key]] <- C_ij
      
    }
  }
  
  return(C_cond)
  
  
}

correlation_estimation_KL_cor <- function(eigendecomp, KL_cor, identity = T){
  
  # ----------------------------------------------------------------------------
  # 
  #
  # GOAL: construct the correlation operator from the KL covariance method, but we already computed KL_cor values
  #
  # input:
  #
  # - eigendecomp (list of 3 entries)
  #
  #   - eigenvalues
  #   - eigenfunctions
  #   - n_dims
  #
  # - KL_cor      (list of i_j entries, each is d x d matrix) each value represents correlation between KL coeffs of components a and b in process i and j
  # - identity    (boolean)  are the i = j entries just the identity?
  #
  # output:
  #
  # - C_cond   (list of i_j entries)
  #
  #
  # ----------------------------------------------------------------------------
  
  
  p <- length(eigendecomp$eigenvalues)
  d <- dim(KL_cov[[1]])[1]
  m <- dim(eigendecomp$eigenfunctions[[1]])[1]
  Delta <- 1/m
  
  C_cond <- list()
  
  for(i in 1:p){
    for(j in i:p){
      
      key <- paste0(i, '_', j)
      
      if(i == j & identity){
        C_ij <- diag(m)
      } else{
        
        
        evec_i <- eigendecomp$eigenfunctions[[i]]
        evec_j <- eigendecomp$eigenfunctions[[j]]
        
        d_i <- dim(evec_i)[2]
        d_j <- dim(evec_j)[2]
        
        cor_ij <- KL_cor[[key]]
        
        if (any(is.nan(cor_ij))) {
          stop("correlation construction process contains NaN values!")
        }      
        
        
        C_ij <- matrix(0, nrow = m, ncol = m)
        
        for(a in 1:d_i){
          for(b in 1:d_j){
            
            # REMEMBER, we could have eigenvectors defined as Delta * (v^\top v) = 1. 
            
            C_ij <- C_ij + cor_ij[a,b] * tcrossprod(evec_i[, a], evec_j[, b])
            
          }
        }
      } # end of nonidentity case 
      
      C_cond[[key]] <- C_ij
      
    }
  }
  
  return(C_cond)
  
  
}
