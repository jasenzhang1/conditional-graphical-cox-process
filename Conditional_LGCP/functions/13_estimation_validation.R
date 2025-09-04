
# 1) want to validate eigendecomposition

validate_eigendecomposition_ii <- function(G_hat, eigen_decomp){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: validate the function `compute_eigendecomposition_ii`
  #
  #
  # Input: 
  # 
  # - G_hat             (m x m x p array) ground truth
  # - eigen_decomp      (list)            list of three lists, each of length p
  # 
  #   - eigenvalues     (list of d-dim vectors)  the top d eigenvalues
  #   - eigenfunctions  (list of mxd matrices)   the top d eigenfunctions (of length m)
  #   - n_dims          (list of integers)       d for each process. we stop at 90% var explained or dmax
  # 
  # 
  # Output: 
  #
  # - G_hat_approx (m x m x p array)  estimate
  #
  # ----------------------------------------------------------------------------  
  
  p <- dim(G_hat)[3]
  m <- dim(G_hat)[2]
  
  G_hat_approx <- array(0, dim = c(m, m, p))
  
  for (i in 1:p) {
    lambdas <- eigen_decomp$eigenvalues[[i]]       # vector length d_i
    etas <- eigen_decomp$eigenfunctions[[i]]       # m x d_i matrix
    d_i <- eigen_decomp$n_dims[[i]]
    
    # Reconstruct covariance matrix for i-th process
    
    if(d_i == 1){
      eta_vec <- matrix(etas[, 1:d_i], m, 1)
      lambda_mat <- matrix(lambdas)
      G_hat_approx[, , i] <- eta_vec %*% lambda_mat %*% t(eta_vec)
    } else{
      G_hat_approx[, , i] <- etas[, 1:d_i] %*% diag(lambdas[1:d_i]) %*% t(etas[, 1:d_i])
    }
    
    
  }
  
  return(G_hat_approx)  
  
}


validate_eigendecomposition_ii_visualization <- function(G_hat, G_hat_approx, i){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: visualize G_hat vs G_hat_approx during KL expansion estimation step
  #
  #
  # Input: 
  # 
  # - G_hat             (p x m x m array) ground truth
  # - G_hat_approx      (p x m x m array) estimate
  # 
  # 
  # 
  # Output: 
  #
  # - grid.arranged graphs of G_hat[i,,] and G_hat_approx[i,,]
  #
  # ----------------------------------------------------------------------------   
  
  # Convert matrix to long format for ggplot
  df1 <- reshape2::melt(G_hat[i,,])
  df2 <- reshape2::melt(G_hat_approx[i,,])
  
  # Heatmap 1
  p1 <- ggplot(df1, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    ggtitle("Heatmap 1") +
    coord_fixed()
  
  # Heatmap 2
  p2 <- ggplot(df2, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    ggtitle("Heatmap 2") +
    coord_fixed()
  
  # Arrange side by side
  return(grid.arrange(p1, p2, ncol = 2))  
  
}

# 2) want to validate KL expansion

validate_kl_coeffs <- function(eigenfunctions, kl_coeffs){
  
  
  
  # ----------------------------------------------------------------------------
  # 
  # 
  # GOAL: recreate the X_k matrix of log intensities with eigenfunctions nad KL coefficients
  # 
  # Inputs:
  # 
  # - eigenfunctions       (p-dim list)   each entry is a m x d matrix of eigenfunctions for process i.
  # - kl_coeffs            (n x p x d)    all KL coeffs for n subjets and p processes
  #
  #
  #
  # Output:
  #
  # - X_k_reconstruct      (p x m x n)    all m-dim vectors of log intensities for p processes and n subjects
  #
  # ----------------------------------------------------------------------------
  
  
  p <- length(eigenfunctions)
  n <- dim(kl_coeffs)[1]
  m <- nrow(eigenfunctions[[1]])  
  
  X_k_reconstruct <- array(NA_real_, dim = c(p, m, n))
  
  for (i in seq_len(p)) {
    phi_i <- eigenfunctions[[i]]  # m x d_i
    d_i <- dim(phi_i)[2]
    
    
    for (k in seq_len(n)) {
      xi_k <- kl_coeffs[k, i, ][1:d_i]    # length d_i
      X_k_reconstruct[i, , k] <- phi_i %*% matrix(xi_k, nrow = d_i, ncol = 1)
    }
  }
  
  return(X_k_reconstruct)

  
}


validate_kl_full <- function(X_k, eigenfunctions, Tseq, x_name, ncores){
  
  # -----------------------------------------------------
  # 
  # code to parallelize repetitiveness
  #
  # 1) find the sample mean and subtract it from X_k to get X_k_center
  # 2) calculate KL coefficients
  # 3) reconstruct X_k_center
  # 4) add the sample mean to get X_k_reconstruct
  #
  # 5) keep the first subject and reshape it into a df
  # 
  # --------------------------------------------------
  
  # 1)
  mean_mat <- apply(X_k, c(1, 2), mean)  # result is p x m matrix
  X_k_center <- sweep(X_k, c(1, 2), mean_mat, FUN = "-")
  
  # 2) 
  kl_coeffs <- estimate_kl_coefficients_parallel_v2(X_k_center, eigenfunctions, Tseq, ncores)
  
  # 3) 
  X_k_est_reconstruct <- validate_kl_coeffs(eigenfunctions, kl_coeffs)
  
  # 4)
  X_k_est_reconstruct <- sweep(X_k_est_reconstruct, c(1, 2), mean_mat, FUN = "+")   # add back mean
  
  # 5) 
  df1 <- reshape2::melt(X_k[,,1])
  df2 <- reshape2::melt(X_k_est_reconstruct[,,1])
  df1$Var2 <- Tseq[df1$Var2]
  df2$Var2 <- Tseq[df2$Var2]
  df1$cat <- x_name
  df2$cat <- paste0(x_name, ' Reconstruct')
  return(rbind(df1, df2))
  
}

# 3) want to validate V_ij conditional covariance operator

validate_V_ij <- function(alpha_true, eigenfunctions){
  
  # 
  # 
  # GOAL: construct the true (m x m) conditional covariance operator for all i_j process pairs
  #
  # Inputs:
  # 
  # - alpha_true: array(n, p, d)  # true KL coefficients from true eigendecomposition of GP_cov onto true X_k log-intensities



  d <- dim(alpha_true)[3]
  p <- dim(alpha_true)[2]
  n <- dim(alpha_true)[1]
  
  Cij_list <- list()
  
  for(i in 1:p){
    for(j in i:p){
      
      # build n_y x d matrix for processes i and j
      A <- alpha_true[, i, ]  # n x d
      B <- alpha_true[, j, ]  # n x d  
      
      # empirical average of outer products: d x d
      Cij <- crossprod(A, B) / n    # equals (t(A) %*% B) / n
      
      key <- paste0(i, '_', j)
      
      Cij_list[[key]] <- Cij
      
    }
  }
  

  # Cij_list =  all i_j pairs of d x d matrix: the truncated operator in coefficient basis

  V_grid <- Phi %*% Cij %*% t(Phi)   # m x m matrix representation
  
}


# Required: base R only (no extra packages)
# alpha: n x p x dmax array of KL coefficients (use NA or ignore extra cols if d_i < dmax)
# Phi_list: list of length p, Phi_list[[i]] is m x d_i matrix (eigenfunctions on grid)
# idx_y: integer indices (subset of replicates) OR NULL to use weighting via weights vector
# weights: numeric vector length n of nonnegative weights (will be normalized). If NULL and idx_y provided, uniform weights over idx_y.
# returns: list with Cij (d_i x d_j) and V_grid (m x m)

validate_V_ij_v2 <- function(alpha_true, Phi_list, i, j, y_c_strata, query_y_c, gamma_c){
  
  # 1) obtain kernel weights of query_y_c with all Y_mat values
  
  # gamma_c <- 1
  # query_y_c <- 3
  # Phi_list <- eigen_decomp_truth$eigenfunctions
  # alpha_true <- kl_coeffs_truth_v2 # n x p x d
  
  n <- dim(alpha_true)[1]
  p <- dim(alpha_true)[2]
  dmax <- dim(alpha_true)[3]
  
  Phi_i <- Phi_list[[i]]
  Phi_j <- Phi_list[[j]]
  d_i <- ncol(Phi_i)
  d_j <- ncol(Phi_j)
  
  
  w <- apply(y_c_strata, 1, step_6_kernel, y2 = query_y_c, gamma_c = gamma_c) # n-dim vector
  
  
  # Xi: n x d matrix of coeffs for process i (first d)
  # Xj: n x d matrix for process j
  
  Xi <- alpha_true[, i, 1:d_i, drop = TRUE] 
  Xj <- alpha_true[, j, 1:d_j, drop = TRUE] 
  
  m_i <- colSums(w * Xi) / sum(w)
  m_j <- colSums(w * Xj) / sum(w)
  
  Sigma_ij_y <- matrix(0, dmax, dmax)
  for (k in 1:n) {
    Sigma_ij_y <- Sigma_ij_y + w[k] * (Xi[k,]-m_i) %*% t(Xj[k,]-m_j) # outer product to get dxd matrix
  }
  Sigma_ij_y <- Sigma_ij_y / sum(w)
  

  
  
  C_ij_y <- Phi_i %*% Sigma_ij_y %*% t(Phi_j)  
  
  
  return(C_ij_y)
}

# Helper: compute ground-truth for all (i,j) with j>=i and return lists
compute_all_ground_truth <- function(alpha, Phi_list, idx_y = NULL, weights = NULL){
  p <- dim(alpha)[2]
  out_C <- list()
  out_V <- list()
  for(i in 1:p){
    for(j in i:p){
      key <- paste(i, j, sep = "_")
      res <- get_truncated_operator(alpha, Phi_list, i, j, idx_y = idx_y, weights = weights)
      out_C[[key]] <- res$Cij
      out_V[[key]] <- res$V_grid
    }
  }
  return(list(C_list = out_C, V_list = out_V))
}