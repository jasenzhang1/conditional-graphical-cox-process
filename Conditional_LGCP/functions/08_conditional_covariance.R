# helper, in use
evaluate_kernel_weights_at_query <- function(Y_continuous_stratum, query_y_c, gamma_c) {
  
  # ---------------------------------------------------------------------------
  #
  #
  # GOAL: Compute kernel weights, w_k(y_c), for a query point
  #
  # - each row of Y_continuous_stratum is a vector of continuous covariates
  # - feed this vector into the kernel wrt query_y_c and gamma_c
  # 
  # Input:
  # 
  # - Y_continuous_stratum   (n_stratum x q_c matrix)
  # - query_y_c              (q_c dim vector)
  # - gamma_c                (scalar)
  #
  # 
  # Output:
  #
  # - kernel_weights (n_stratum dim vector)
  # 
  #
  # ---------------------------------------------------------------------------
  
  
  kernel_weights <- apply(Y_continuous_stratum, 1, step_6_kernel, y2 = query_y_c, gamma_c = gamma_c)
  
  return(kernel_weights)
}

# deprecated
evaluate_regression_at_query <- function(M_hat, Y_continuous_stratum, query_y_c, 
                                         eigenfunctions, gamma_c, p) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: evaluate regression at query
  #
  #
  # Input: 
  #
  # - M_hat                  (p^2 length list where each entry is 'i_j') 
  # - Y_continuous_stratum   (n_stratum x q_c matrix)
  # - query_y_c              (q_c dim vector)
  # - eigenfunctions         (list of p matrices, each of which is m x d_i)
  # - gamma_c                (number)
  # - p                      (integer)
  #
  #
  # Output: 
  #
  # - V_conditional (list of length p^2, each element is m x m matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  components <- sapply(eigenfunctions, ncol) %>% unlist()
  max_components <- max(components)            # d
  n_time <- nrow(eigenfunctions[[1]])          # m
  n_stratum <- nrow(Y_continuous_stratum)      # n 
  
  # Compute kernel weights: n_stratum x 1 vector
  kernel_weights <- evaluate_kernel_weights_at_query(Y_continuous_stratum, query_y_c, gamma_c)
  
  V_conditional <- list()
  
  for (i in 1:p) {
    for (j in 1:p) {
      key <- paste(i, j, sep="_")
      
      # Construct conditional covariance operator: m x m matrix
      V_cond_ij <- matrix(0, nrow=n_time, ncol=n_time)
      
      for (a in 1:components[i]) {
        for (b in 1:components[j]) {
          # Evaluate regression coefficient at query point: scalar
          M_ab_at_query <- sum(kernel_weights * M_hat[[key]][, a, b])
          
          eta_ia <- eigenfunctions[[i]][, a]  # m x 1 vector
          eta_jb <- eigenfunctions[[j]][, b]  # m x 1 vector
          # Tensor product: (m x 1) %*% (1 x m) = (m x m)
          tensor_prod <- outer(eta_ia, eta_jb)
          # Weighted sum: scalar * (m x m) + (m x m) = (m x m)
          V_cond_ij <- V_cond_ij + M_ab_at_query * tensor_prod
        }
      }
      
      V_conditional[[key]] <- V_cond_ij
    }
  }
  
  return(V_conditional)
}

# deprecated
evaluate_regression_at_query_diag_only <- function(M_hat, Y_continuous_stratum, query_y_c, eigenfunctions, gamma_c, p) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: evaluate regression at query
  #
  #   - but we only evaluate for j >= i
  #
  # Input: 
  #
  # - M_hat                  (p^2 length list where each entry is 'i_j') 
  # - Y_continuous_stratum   (n_stratum x q_c matrix)
  # - query_y_c              (q_c dim vector)
  # - eigenfunctions         (list of p matrices, each of which is m x d_i)
  # - gamma_c                (number)
  # - p                      (integer)
  #
  #
  # Output: 
  #
  # - V_conditional (list of length p^2, each element is m x m matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  components <- sapply(eigenfunctions, ncol) %>% unlist()
  max_components <- max(components)            # d
  n_time <- nrow(eigenfunctions[[1]])          # m
  n_stratum <- nrow(Y_continuous_stratum)      # n 
  
  # Compute kernel weights: n_stratum x 1 vector
  kernel_weights <- evaluate_kernel_weights_at_query(Y_continuous_stratum, query_y_c, gamma_c)
  
  V_conditional <- list()
  
  for (i in 1:p) {
    for (j in i:p) {
      key <- paste(i, j, sep="_")
      
      # Construct conditional covariance operator: m x m matrix
      V_cond_ij <- matrix(0, nrow=n_time, ncol=n_time)
      
      for (a in 1:components[i]) {
        for (b in 1:components[j]) {
          # Evaluate regression coefficient at query point: scalar
          M_ab_at_query <- sum(kernel_weights * M_hat[[key]][, a, b])
          
          eta_ia <- eigenfunctions[[i]][, a]  # m x 1 vector
          eta_jb <- eigenfunctions[[j]][, b]  # m x 1 vector
          # Tensor product: (m x 1) %*% (1 x m) = (m x m)
          tensor_prod <- outer(eta_ia, eta_jb)
          # Weighted sum: scalar * (m x m) + (m x m) = (m x m)
          V_cond_ij <- V_cond_ij + M_ab_at_query * tensor_prod
        }
      }
      
      V_conditional[[key]] <- V_cond_ij
    }
  }
  
  return(V_conditional)
}

# in use
evaluate_regression_at_query_v2 <- function(M_hat, Y_continuous_stratum, query_y_c, 
                                         eigenfunctions, gamma_c, p) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: evaluate regression at query
  #
  # - v2: j >= i only
  # - but also vectorize each (i, j) calculation wrt the eigenfunctions
  #
  # Input: 
  #
  # - M_hat                  (p^2 length list where each entry is 'i_j') 
  # - Y_continuous_stratum   (n_stratum x q_c matrix)
  # - query_y_c              (q_c dim vector)
  # - eigenfunctions         (list of p matrices, each of which is m x d_i)
  # - gamma_c                (number)
  # - p                      (integer)
  #
  #
  # Output: 
  #
  # - V_conditional (list of length p^2, each element is m x m matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  components <- sapply(eigenfunctions, ncol) %>% unlist()
  max_components <- max(components)            # d
  n_time <- nrow(eigenfunctions[[1]])          # m
  n_stratum <- nrow(Y_continuous_stratum)      # n 
  
  # Compute kernel weights: n_stratum x 1 vector
  kernel_weights <- evaluate_kernel_weights_at_query(Y_continuous_stratum, query_y_c, gamma_c)
  
  V_conditional <- list()
  
  for (i in 1:p) {
    for (j in i:p) {
      key <- paste(i, j, sep="_")
      
      # Construct conditional covariance operator: m x m matrix
      V_cond_ij <- matrix(0, nrow=n_time, ncol=n_time)
      
      # recall M_hat[[key]] = n x comps x comps
      # flatten it to be n x (comps^2)
      M_hat_ij_flat <- matrix(M_hat[[key]], nrow = n_stratum, ncol = max_components^2) 
      M_at_query <- apply(kernel_weights * M_hat_ij_flat, 2, sum) # (comps^2) length vector
      
      
      
      # Split columns into vectors
      A_cols <- asplit(eigenfunctions[[i]], 2)  # list of m vectors (length 19)
      B_cols <- asplit(eigenfunctions[[j]], 2)  # list of n vectors
      
      
      # append A_cols and B_cols with 0's so that they have a list of (comps) vectors
      
      while(length(A_cols) < max_components){
        A_cols[[length(A_cols) + 1]] <- rep(0, n_time)
      } 
      
      while(length(B_cols) < max_components){
        B_cols[[length(B_cols) + 1]] <- rep(0, n_time)
      }       
      
      # Compute all d × d outer products using tcrossprod
      # Result: list of d * d matrices
      
      # BE CAREFUL!!! B_cols needs to be repeated for the amount of components of i
      #               A_cols needs to be repeated for the amount of components of j
      outer_products <- Map(tcrossprod,
                            rep(A_cols, times = max_components),
                            rep(B_cols,  each = max_components))
      
      # Weighted sum of all outer products
      result <- Reduce(`+`, Map(`*`, M_at_query, outer_products))
      
      V_conditional[[key]] <- result
      
    }
  }
  
  return(V_conditional)
}

