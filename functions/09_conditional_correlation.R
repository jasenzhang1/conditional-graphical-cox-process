
correlation_estimation_KL_cor <- function(eigendecomp, KL_cor, identity = F){
  
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
  #   - eigenfunctions   NOT NORMALIZED
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
  d <- dim(KL_cor[[1]])[1]
  m <- dim(eigendecomp$eigenfunctions[[1]])[1]
  delta <- 1/m
  
  C_cond <- list()
  C_cond_unnorm <- list()
  
  for(i in 1:p){
    for(j in i:p){
      
      key <- paste0(i, '_', j)
      
      if(i == j & identity){
        stop('Error 09: do not use the identity route')
        C_ij <- diag(m)
      } else{
        
        
        evec_i <- eigendecomp$eigenfunctions[[i]] * sqrt(delta)
        evec_j <- eigendecomp$eigenfunctions[[j]] * sqrt(delta)
        
        d_i <- dim(evec_i)[2]
        d_j <- dim(evec_j)[2]
        
        cor_ij <- KL_cor[[key]]
        
        if (any(is.nan(cor_ij))) {
          stop("correlation construction process contains NaN values!")
        }      
        
        
        # C_ij <- matrix(0, nrow = m, ncol = m)
        # 
        # for(a in 1:d_i){
        #   for(b in 1:d_j){
        #     
        #     # REMEMBER, EIGENFUNCTIONS ARE NOT NORMALIZED
        #     
        #     C_ij <- C_ij + cor_ij[a,b] * Delta^2 * tcrossprod(evec_i[, a], evec_j[, b])
        #     
        #   }
        # }
        
        C_ij <- evec_i %*% cor_ij %*% t(evec_j)
        C_ij_unnorm <- eigendecomp$eigenfunctions[[i]] %*% cor_ij %*% t(eigendecomp$eigenfunctions[[j]])
        
      } # end of nonidentity case 
      
      C_cond[[key]] <- C_ij
      C_cond_unnorm[[key]] <- C_ij_unnorm
      
    }
  }
  
  return(list(C_cond = C_cond,
              C_cond_unnorm = C_cond_unnorm))
  
  
}


