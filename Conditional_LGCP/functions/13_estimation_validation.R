
# 1) want to validate KL expansion

validate_eigendecomposition_ii <- function(G_hat, eigen_decomp){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: validate the function `compute_eigendecomposition_ii`
  #
  #
  # Input: 
  # 
  # - G_hat             (p x m x m array) ground truth
  # - eigen_decomp      (list)            list of three lists, each of length p
  # 
  #   - eigenvalues     (list of d-dim vectors)  the top d eigenvalues
  #   - eigenfunctions  (list of mxd matrices)   the top d eigenfunctions (of length m)
  #   - n_dims          (list of integers)       d for each process. we stop at 90% var explained or dmax
  # 
  # 
  # Output: 
  #
  # - G_hat_approx (p x m x m array)  estimate
  #
  # ----------------------------------------------------------------------------  
  
  p <- dim(G_hat)[1]
  m <- dim(G_hat)[2]
  
  G_hat_approx <- array(0, dim = c(p, m, m))
  
  for (i in 1:p) {
    lambdas <- eigen_decomp$eigenvalues[[i]]       # vector length d_i
    etas <- eigen_decomp$eigenfunctions[[i]]       # m x d_i matrix
    d_i <- eigen_decomp$n_dims[[i]]
    
    # Reconstruct covariance matrix for i-th process
    G_hat_approx[i, , ] <- etas[, 1:d_i] %*% diag(lambdas[1:d_i]) %*% t(etas[, 1:d_i])
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
