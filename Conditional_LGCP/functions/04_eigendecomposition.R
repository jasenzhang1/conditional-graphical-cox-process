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
  # - G_hat        (list of i_j entries)
  # - p            (integer)
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

compute_eigendecomposition_ii <- function(G_hat, same_basis, constant_d, var_explained = 0.95) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: from G_{i,i}(s,t), estimate the eigendecomposition
  #
  # - note that we only have G_{i, i} entries 
  #
  # Input: 
  # 
  # - G_hat        (m x m x p matrix)
  # - same_basis   (boolean)            if true, use the trig basis in finite_basis procedure
  # - constant_d   (integer)            if not null, each process gets d eigencomponents                        
  # - var_explained (percentage)        used to calculate number of eigencomponents if constant_d = null
  # 
  # 
  # Output: 
  #
  # - eigenvalues    (list of p vectors)
  # - eigenfunctions (list of p matrices)
  # - n_dims         (list of p numbers denoting d_i)
  #
  # ----------------------------------------------------------------------------
  
  if(same_basis){
    source('functions/21b_generate_finite_basis_expansion.R')
  }
  
  p <- dim(G_hat)[3]
  m <- dim(G_hat)[1]
  
    
  Delta <- 1/m         # integration constant
  max_possible_d <- m  # Helper to ensure we don't exceed m components
  
  eigenvalues <- list()
  eigenfunctions <- list()
  eigenfunctions_regular <- list()
  n_dims <- list()
  

  
  for (i in 1:p) {
    
    G_ii <- G_hat[, , i]

    # Ensure symmetry
    G_ii <- (G_ii + t(G_ii)) / 2
    
    if(same_basis){
      # ---------------------------------------------------------
      # SCENARIOS 1 & 2: Fixed Trig Basis Projection
      # ---------------------------------------------------------
      
      # Determine how many trig functions to generate
      # If constant_d is NULL, we generate up to max and truncate later
      d_to_gen <- if(!is.null(constant_d)) constant_d else max_possible_d
      
      # Get basis functions from your helper
      basis_fns <- trig_basis(d_to_gen)
      
      # Convert function list to m x d matrix and normalize for the grid
      time_grid <- make_time_grid(m)
      Phi <- trig_basis_realization(basis_fns, time_grid) # m x d

      
      # We calculate eigenvalues via the quadratic form: lambda_a = phi_a' * G * phi_a * Delta^2
      # However, since G_hat usually incorporates one Delta in intensity estimation, 
      # we check the scaling. Standard discrete projection:
      lambdas_trig <- sapply(1:d_to_gen, function(a) {
        phi_a <- Phi[, a]
        as.numeric(t(phi_a) %*% G_ii %*% phi_a) * (Delta^2)
      })
      
      # Scenario 2: If constant_d is NULL, find d_i based on var_explained
      if (is.null(constant_d)) {
        lambdas_trig[lambdas_trig < 0] <- 0
        cum_var <- cumsum(lambdas_trig) / sum(lambdas_trig)
        d_i <- which(cum_var >= var_explained)[1]
        if(is.na(d_i)) d_i <- d_to_gen
      } else {
        d_i <- constant_d
      }
      
      eigenvalues[[i]] <- lambdas_trig[1:d_i]
      eigenfunctions[[i]] <- matrix(Phi[, 1:d_i], nrow = m)
      n_dims[[i]] <- d_i
    } else{
      # ---------------------------------------------------------
      # SCENARIOS 3 & 4: Empirical PCA (Standard Eigendecomposition)
      # ---------------------------------------------------------
      
      # Standardize G for eigen() to account for discretization
      G_ii_norm <- G_ii * Delta 
      
      eigen_result <- tryCatch({
        eigen(G_ii_norm, symmetric = TRUE)
      }, error = function(e) {
        stop(paste("Eigen error in process", i, ":", e$message))
      })
      
      etas <- eigen_result$vectors    # m x m
      lambdas <- eigen_result$values  # m x 1
      lambdas[lambdas < 0] <- 0
      
      # Determine truncation point d_i
      if (!is.null(constant_d)) {
        # Scenario 4 (Manual Truncation)
        d_i <- constant_d
      } else {
        # Scenario 3 (Variance Explained)
        cum_var <- cumsum(lambdas) / sum(lambdas)
        d_i <- which(cum_var >= var_explained)[1]
      }
      
      # Store results
      eigenvalues[[i]] <- lambdas[1:d_i]
      
      # Normalize eigenfunctions so that integral of phi^2 = 1
      # eigen() vectors have sum(v^2) = 1, so we divide by sqrt(Delta)
      phi_matrix <- matrix(etas[, 1:d_i], nrow = m) / sqrt(Delta)
      
      eigenfunctions[[i]] <- phi_matrix
      n_dims[[i]] <- d_i
    }
  }
  
  return(list(eigenvalues = eigenvalues, 
              eigenfunctions = eigenfunctions, 
              n_dims = n_dims))
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
