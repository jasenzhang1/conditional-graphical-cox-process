# 9/4/2025
#
# - construct conditional covariance using information from G_ij 

evaluate_regression_at_query_g_ij <- function(M_hat, df_pair, y_c_strata, y_c_k, eigenfunctions, gamma_c) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: calculate M_hat regression operator from c_ij^ab^k
  #
  # 
  # input:
  #
  # - M_hat        (list of all i_j pairs) where each entry is (n x d_i x d_j matrix)
  #                                         - d_i = number of eigencomponents of process i
  #                                         - d_j = number of eigencomponents of process j
  #
  # - df_pair         (data.frame)             lookup table for i_j pair names
  # - y_c_strata      (n x q_c matrix)         matrix of continuous covariates
  # - y_c_k           (q_c dim vector)         covariate vector for subject k
  # - eigenfunctions  (list of p matrices)     where each matrix is (m x d_i)
  # - ncores          (number)
  #
  # output:
  #
  # - V_hat
  #
  # ----------------------------------------------------------------------------  
  
  # Input dimensions
  n_stratum <- dim(y_c_strata)[1]
  m <- dim(eigenfunctions[[1]])[1]

  
  # Compute kernel weights at query
  weights <- sapply(1:n_stratum, function(k) step_6_kernel(y_c_strata[k,], y_c_k, gamma_c))  # n_y_d x 1
  
  V_hat <- pbmclapply(1:nrow(df_pair), function(l) { 
  
    i <- df_pair$feature_i[l]
    j <- df_pair$feature_j[l]
    
    M_ij <- M_hat[[l]]
    
    d_i <- dim(M_ij)[2]
    d_j <- dim(M_ij)[3]
    V_ij <- matrix(0, m, m) 
    
    for (a in 1:d_i) {
      for (b in 1:d_j) {
        m_ab_query <- sum(weights * M_hat[[l]][, a, b])  # scalar = sum(nx1 * nx1)
        tensor_ab <- outer(eigenfunctions[[i]][, a], eigenfunctions[[j]][, b])  # m x m
        V_ij <- V_ij + m_ab_query * tensor_ab  #  m x m = scalar * mxm
      }
    }
    
    V_ij

  
  }, mc.cores = ncores)
  
  names(V_hat) <- paste0(df_pair$feature_i, '_', df_pair$feature_j)  
  
  return(V_hat)
}
