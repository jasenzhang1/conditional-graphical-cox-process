# 9/4/2025

# In our new method that uses G_ij, we need to get cross-KL coeffs for all subjects
#
#   - used to be alpha_i^k, for process i and subject k
#   - now it's c_ij^k, for processes i and j and subject k


estimate_cross_kl_coefficients_parallel <- function(g_ij_k_est, df_pair, eigenfunctions, t_seq, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Obtain cross KL coefficients for subject k, process i and j
  #
  #
  #
  # Input: 
  # 
  # - g_ij_k_est        (list of all i_j pairs)  each entry is a (m x m x n) matrix
  # - df_pair           (dataframe)              dataframe that specifies i_j
  # - eigenfunctions    (list of p matrices of dimension m x d_i) 
  # - t_seq             (vector of length m)
  # - ncores            
  #
  # Output: 
  #
  # - alpha_tensor (n_stratum x p x d array)
  #
  # ----------------------------------------------------------------------------
  
  if(length(t_seq) != dim(eigenfunctions[[1]])[1]){
    print('tseq and eigenfunction dimensions are incompatible')
    return(NULL)
  }
  
  # Input dimensions
  n_stratum <- dim(g_ij_k_est[[1]])[3]  
  p <- length(unique(df_pair$feature_i))
  m <- dim(g_ij_k_est[[1]])[1] 
  dt <- diff(t_seq)[1]  # assume uniform grid

  c_hat <- pbmclapply(1:nrow(df_pair), function(l) { 
    i <- df_pair$feature_i[l]
    j <- df_pair$feature_j[l]
    key <- paste0(i, '_', j)
    
    d_i <- dim(eigenfunctions[[i]])[2]
    d_j <- dim(eigenfunctions[[j]])[2]
    
    g_ij <- g_ij_k_est[[l]]
    
    
    c_abk <- array(0, dim = c(d_i, d_j, n_stratum))
    for (k in 1:n_stratum) {
      g_ij_k <- g_ij[,,k]
      for (a in 1:d_i) {
        for (b in 1:d_j) {
          eta_ia <- eigenfunctions[[i]][,a]
          eta_jb <- eigenfunctions[[j]][,b]
          proj <- outer(eta_ia, eta_jb)
          integrand <- g_ij_k * proj  # m x m
          c_abk[a][b][k] <- sum(integrand) * dt * dt  # scalar
        }
      }
    }
    
    c_abk
    
  }, mc.cores = ncores)
  
  names(c_hat) <- paste0(df_pair$feature_i, '_', df_pair$feature_j)
  
  return(c_hat)
}


# 9/22/2025

# we now estimate the covariance between KL coefficients as per CPGM paper

estimate_KL_covariance <- function(G_hat, eigenfunctions, norm_G){
  
  # ----------------------------------------------------------------------------
  # 
  #
  # GOAL: estimate cov(alpha_i^a, alpha_j^b) between processes i and j and eigencomponents a and b
  #
  #
  # inputs:
  #
  # - G_hat            (list of i_j m x m matrices)              each i_j is a G_{i,j}(s,t) covariance matrix
  # - eigenfunctions   (list of p entries, m x d_i matrices)     each entry represents the eigenfunctions for the i-th process
  # - norm_G           (boolean)    do we apply G_ii <- G_ii / m to get constant eigenvalues?
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
      
      if(norm_G){
        cov_ij <- t(eigen_i) %*% (G_ij / m) %*% eigen_j 
      } else{
        cov_ij <- t(eigen_i) %*% G_ij %*% eigen_j 
      }
      
      
      # in case we need to force a d_max x d_max result
      # cov_ij2 <- matrix(0, nrow = d, ncol = d)
      # 
      # cov_ij2[1:dim(eigen_i)[2], 1:dim(eigen_j)[2]] <- cov_ij
      
      KL_cov[[key]] <- cov_ij
    }
  }
  
  return(KL_cov)
  
}
