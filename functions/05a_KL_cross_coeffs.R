# 9/4/2025

# In our new method that uses G_ij, we need to get cross-KL coeffs for all subjects
#
#   - used to be alpha_i^k, for process i and subject k
#   - now it's c_ij^k, for processes i and j and subject k



# 9/22/2025

# we now estimate the covariance between KL coefficients as per CPGM paper

estimate_KL_covariance <- function(G_hat, eigenfunctions){
  
  # ----------------------------------------------------------------------------
  # 
  #
  # GOAL: estimate cov(alpha_i^a, alpha_j^b) between processes i and j and eigencomponents a and b
  #
  #
  # inputs:
  #
  # - G_hat            (list of i_j m x m matrices)              each i_j is a G_{i,j}(s,t) covariance matrix
  # - eigenfunctions   (list of p entries, m x d_i matrices)     each entry represents the eigenfunctions for the i-th process, we assume they are not normalized
  #
  #
  # output:
  #
  # - KL_cov           (list of i_j d_max x d_max matrices)     each entry represents the covariance of KL coefficients for process i and j
  #
  #
  # ----------------------------------------------------------------------------
  
  d <- max(sapply(eigenfunctions, ncol))
  p <- length(eigenfunctions)
  m <- dim(G_hat[[1]])[1]
  
  KL_cov <- list()
  
  for(i in 1:p){
    for(j in i:p){
      key <- paste0(i, '_', j)
      
      G_ij <- G_hat[[key]]
      
      eigen_i <- eigenfunctions[[i]]
      eigen_j <- eigenfunctions[[j]]
      
      W <- diag(1/m, m)

      cov_ij <- t(eigen_i) %*% W %*% G_ij %*% W %*% eigen_j   # normalize the eigenfunctions

      
      
      # in case we need to force a d_max x d_max result
      # cov_ij2 <- matrix(0, nrow = d, ncol = d)
      # 
      # cov_ij2[1:dim(eigen_i)[2], 1:dim(eigen_j)[2]] <- cov_ij
      
      KL_cov[[key]] <- cov_ij
    }
  }
  
  return(KL_cov)
  
}

estimate_KL_correlation <- function(KL_cov, p){
  
  # ----------------------------------------------------------------------------
  # 
  #
  # GOAL: estimate cor(alpha_i^a, alpha_j^b) between processes i and j and eigencomponents a and b
  #
  #       if i == j, correlation is the identity no matter what
  #
  # inputs:
  #
  # - KL_cov           (list of i_j d_max x d_max matrices)     each entry represents the covariance of KL coefficients for process i and j
  # - p                (integer)
  #
  # output:
  #
  # - KL_cor           (list of i_j d_max x d_max matrices)     each entry represents the correlation of KL coefficients for process i and j 
  #
  # ----------------------------------------------------------------------------
  
  # --- Compute correlations using the KL variances ---
  
  KL_cor <- list()
  
  for(i in 1:p){
    for(j in i:p){
      
      key <- paste0(i, '_', j)
      cov_ij <- KL_cov[[key]]
      
      if(i == j){ # if i == j, make the identity
        
        d <- dim(cov_ij)[1]
        KL_cor[[key]] <- diag(d)
        
      } else{
        
        var_i  <- diag(KL_cov[[paste0(i, '_', i)]])
        var_j  <- diag(KL_cov[[paste0(j, '_', j)]])
        
        # Outer product of std deviations
        denom <- sqrt(outer(var_i, var_j))
        cor_ij <- cov_ij / denom
        
        KL_cor[[key]] <- cor_ij
        
      }
    }
  }
  
  return(KL_cor)
  
}

estimate_KL_precision <- function(KL_cor, p){
  
  # ----------------------------------------------------------------------------
  # 
  #
  # GOAL: estimate KL_precision between processes i and j and eigencomponents a and b
  #
  #
  # inputs:
  #
  # - KL_cor           (list of i_j d_max x d_max matrices)     each entry represents the correlation of KL coefficients for process i and j
  # - p                (integer)
  #
  # output:
  #
  # - KL_prec           (list of i_j d_max x d_max matrices)     each entry represents the precision of KL coefficients for process i and j 
  #
  # ----------------------------------------------------------------------------
  
  # assemble --> inverse --> extract
  
  
  KL_cor_bundle <- assemble_block_matrix_irregular(KL_cor, p)
  KL_prec_full <- sym(ginv(KL_cor_bundle$block_matrix))
  KL_prec <- extract_block_matrix_irregular(KL_prec_full, KL_cor_bundle$row_borders, KL_cor_bundle$col_borders)

  
  return(KL_prec)
  
}
