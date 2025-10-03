compute_eigendecomposition <- function(G_hat, var_explained = 0.9) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: from G_{i,j}(s,t), estimate the eigendecomposition
  #
  # 
  # Input: 
  # 
  # - G_hat (m x m array) (p x p x m x m)
  # 
  # 
  # Output: 
  #
  # - eigenvalues    (list of p vectors)
  # - eigenfunctions (list of p matrices)
  # - n_dims         (list of p numbers denoting d_i)
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(G_hat)[1]
  n_time <- dim(G_hat)[3]  # m
  
  eigenvalues <- list()
  eigenfunctions <- list()
  n_dims <- list()
  
  for (i in 1:p) {
    # Extract marginal covariance matrix: m x m
    G_ii <- G_hat[i, i, , ]  
    
    # Ensure symmetry for numerical stability
    G_ii <- (G_ii + t(G_ii)) / 2
    
    # Compute eigendecomposition
    # eigen() returns: values (m x 1), vectors (m x m)
    eigen_result <- eigen(G_ii, symmetric=TRUE)
    lambdas <- eigen_result$values  # m x 1 vector
    lambdas2 <- lambdas
    lambdas2[lambdas2 < 0] <- 0
    
    cum_var <- cumsum(lambdas2) / sum(lambdas2)
    d_i <- which(cum_var > var_explained)[1]
    etas <- eigen_result$vectors    # m x m matrix
    

    
    # Keep only top d components

    eigenvalues[[i]] <- lambdas[1:d_i]        # d x 1 vector
    
    if(d_i == 1){
      eigenfunctions[[i]] <- matrix(etas[, 1:d_i], nrow = n_time)
    } else{
      eigenfunctions[[i]] <- etas[, 1:d_i]      # m x d matrix
    }
    
    n_dims[[i]] <- d_i
  }
  
  # keep track of how many components are needed to explain 90%
  
  return(list(eigenvalues = eigenvalues, eigenfunctions = eigenfunctions, n_dims = n_dims))
}

# prep

prep_eigendecomposition_ii <- function(G_hat, p){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: restructure G_ij so that it can be fed into `compute_eigendecomposition_ii`
  #
  # - only keep the i_i values from the list
  #
  #
  # input:
  #
  # - G_hat (list of i_j entries)
  # - p     (integer)
  #
  # 
  # output:
  #
  # - G_mat (m x m x p matrix)
  # 
  #
  # ----------------------------------------------------------------------------
  
  m <- dim(G_hat[[1]])[1]
  
  G_mat <- array(0, dim = c(m, m, p))
  
  for(i in 1:p){
    key = paste0(i, '_', i)
    
    G_mat[,,i] <- G_hat[[key]]
    
    
  }
  
  return(G_mat)
}

compute_eigendecomposition_ii <- function(G_hat, inner_1, var_explained = 0.999) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: from G_{i,i}(s,t), estimate the eigendecomposition
  #
  # - note that we only have G_{i, i} entries 
  #
  # Input: 
  # 
  # - G_hat (m x m array) (m x m x p)
  # - inner_1       (boolean)  does the inner product of eta^\top eta = 1? 
  #                            If false, Delta * eta^\top \eta = 1
  # - var_explained (percentage)
  # 
  # 
  # Output: 
  #
  # - eigenvalues    (list of p vectors)
  # - eigenfunctions (list of p matrices)
  # - n_dims         (list of p numbers denoting d_i)
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(G_hat)[3]
  m <- dim(G_hat)[1]
  Delta <- 1/m
  
  eigenvalues <- list()
  eigenfunctions <- list()
  n_dims <- list()
  
  for (i in 1:p) {
    # Extract marginal covariance matrix: m x m
    
    if(inner_1){
      G_ii <- G_hat[, , i]  # / m
    } else{
      G_ii <- G_hat[, , i]  # / m
    }
    
    
    # Ensure symmetry for numerical stability
    G_ii <- (G_ii + t(G_ii)) / 2
    
    # Compute eigendecomposition
    # eigen() returns: values (m x 1), vectors (m x m)
    eigen_result <- eigen(G_ii, symmetric=TRUE)
    lambdas <- eigen_result$values  # m x 1 vector
    lambdas2 <- lambdas
    lambdas2[lambdas2 < 0] <- 0
    
    cum_var <- cumsum(lambdas2) / sum(lambdas2)
    d_i <- which(cum_var > var_explained)[1]
    etas <- eigen_result$vectors    # m x m matrix
    
    
    
    # Keep only top d components
    
    eigenvalues[[i]] <- lambdas[1:d_i]        # d x 1 vector
    
    if(inner_1){
      if(d_i == 1){
        eigenfunctions[[i]] <- matrix(etas[, 1:d_i], nrow = m)
      } else{
        eigenfunctions[[i]] <- etas[, 1:d_i]      # m x d matrix
      }
    } else{
      if(d_i == 1){
        eigenfunctions[[i]] <- matrix(etas[, 1:d_i], nrow = m) # / sqrt(Delta)
      } else{
        eigenfunctions[[i]] <- etas[, 1:d_i] # / sqrt(Delta)     # m x d matrix
      }      
    }
    
    n_dims[[i]] <- d_i
  }
  
  # keep track of how many components are needed to explain 90%
  
  return(list(eigenvalues = eigenvalues, eigenfunctions = eigenfunctions, n_dims = n_dims))
}

compute_eigendecomposition_dmax <- function(G_hat, d_max = 10) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: from G_{i,j}(s,t), estimate the eigendecomposition
  #
  # 
  # Input: 
  # 
  # - G_hat (m x m array) (p x p x m x m)
  # 
  # 
  # Output: 
  #
  # - eigenvalues (list of p vectors)
  # - eigenfunctions (list of p matrices)
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(G_hat)[1]
  n_time <- dim(G_hat)[3]  # m
  
  eigenvalues <- list()
  eigenfunctions <- list()
  n_dims <- list()
  
  for (i in 1:p) {
    # Extract marginal covariance matrix: m x m
    G_ii <- G_hat[i, i, , ]  
    
    # Ensure symmetry for numerical stability
    G_ii <- (G_ii + t(G_ii)) / 2
    
    # Compute eigendecomposition
    # eigen() returns: values (m x 1), vectors (m x m)
    eigen_result <- eigen(G_ii, symmetric=TRUE)
    lambdas <- eigen_result$values  # m x 1 vector
    lambdas2 <- lambdas
    lambdas2[lambdas2 < 0] <- 0
    
    cum_var <- cumsum(lambdas2) / sum(lambdas2)
    etas <- eigen_result$vectors    # m x m matrix
    
    
    
    # Keep only d_max components
    
    eigenvalues[[i]] <- lambdas[1:d_max]        # d x 1 vector
    
    eigenfunctions[[i]] <- etas[, 1:d_max]      # m x d matrix
    
    
  }
  
  # keep track of how many components are needed to explain 90%
  
  return(list(eigenvalues = eigenvalues, eigenfunctions = eigenfunctions))
}

select_components_by_variance <- function(eigenvalues, variance_threshold=0.9) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Select number of components to explain given variance
  #
  #
  # input:
  #
  # - eigenvalues           (p-dim list of vectors of length d_max)   vector of sorted eigenvalues for the p-th 
  # - variance_threshold    (number)                                  percent of variance explained threshold         
  # 
  components_selected <- list()
  
  for (i in 1:length(eigenvalues)) {
    lambdas <- eigenvalues[[i]]
    # Only positive eigenvalues
    positive_lambdas <- pmax(lambdas, 0)
    total_var <- sum(positive_lambdas)
    cumsum_var <- cumsum(positive_lambdas)
    
    if (total_var > 0) {
      n_components <- which(cumsum_var >= variance_threshold * total_var)[1]
      if (is.na(n_components)) n_components <- length(lambdas)
    } else {
      n_components <- 1
    }
    
    components_selected[[i]] <- n_components
  }
  
  return(components_selected)
}
