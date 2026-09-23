# 9/4/2025

# use the JASA approach to obtain regression operator

# use construct_kernel_matrix_step_6 from step 6

estimate_rkhs_regression_operators <- function(c_hat, df_pair, y_c_strata, gamma_c, ncores) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: calculate M_hat regression operator from c_ij^ab^k
  #
  # 
  # input:
  #
  # - c_hat        (list of all i_j pairs)  where each entry is (d_i x d_j x n matrix)
  #                                         - d_i = number of eigencomponents of process i
  #                                         - d_j = number of eigencomponents of process j
  #
  # - df_pair      (data.frame)             lookup table for i_j pair names
  # - y_c_strata   (n x q_c matrix)         matrix of continuous covariates
  # - gamma_c      (number)                 K_c bandwidth
  # - ncores       (number)
  #
  # output:
  #
  # - M_hat        (list of all i_j pairs) where each entry is (n x d_i x d_j matrix)
  #
  # ----------------------------------------------------------------------------
  
  # Input dimensions
  n_stratum <- dim(y_c_strata)[1]
  
  K <- construct_kernel_matrix_step_6(y_c_strata, gamma_c)  # n_y_d x n_y_d
  V_YcYc <- K / n_stratum
  
  
  V_inv <- solve(psd_jitter(V_YcYc))  # n_y_d x n_y_d

  
  M_hat <- pbmclapply(1:nrow(df_pair), function(l) { 
    i <- df_pair$feature_i[l]
    j <- df_pair$feature_j[l]
    
    c_ij <- c_hat[[l]]
    
    d_i <- dim(c_ij)[1]
    d_j <- dim(c_ij)[2]
    
    M_ij <- c_ij * 0 # d_i x d_j x n
    M_ij <- aperm(M_ij, perm = c(3, 1, 2))  # make it n x d_i x d_j
    
    
    for (a in 1:d_i) {
      for (b in 1:d_j) {
        response <- c_ij[a, b, ] / n_stratum  # n_y_d x 1  taking the column but dividing every value by n? CHECK PLEASE 9/4/2025
        M_ij[, a, b] <- V_inv %*% response    # n_y_d x 1
      }
    }
    
    M_ij

  }, mc.cores = ncores)
  
  names(M_hat) <- paste0(df_pair$feature_i, '_', df_pair$feature_j)
  
  return(M_hat)
}
