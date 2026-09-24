
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



compute_eigendecomposition_mfpca <- function(G_hat_list, p, same_basis, constant_d,
                                             var_explained = 0.95) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: from all G_{i,j}(s,t), estimate the MFPCA eigendecomposition
  #       via Happ & Greven (2018) two-step approach:
  #         Step A — univariate FPCA on each G_{ii} -> basis Phi_i (m x d_i)
  #         Step B — project all G_{ij} onto Phi_i x Phi_j -> matrix C (D x D)
  #         Step C — eigendecompose C -> joint eigenvalues & eigenfunctions
  #
  # Input:
  #   - G_hat_list   (named list of m x m matrices, keys "i_j", i <= j)
  #   - p            (integer) number of processes
  #   - same_basis   (boolean) if TRUE use trig basis; if FALSE use empirical PCA
  #   - constant_d   (integer or NULL) fixed d_i per process; if NULL use var_explained
  #   - var_explained (numeric) cumulative variance threshold when constant_d = NULL
  #
  # Output:  list with
  #   - eigenvalues     numeric(D)        joint eigenvalues, descending
  #   - eigenfunctions  list of p (m x D) multivariate eigenfunctions per process
  #   - univariate      list of p lists   per-process eigenvalues & eigenfunctions
  #   - n_dims          list of p integers d_i per process
  #   - C               (D x D) matrix    projected covariance matrix
  #
  # ----------------------------------------------------------------------------
  
  if (same_basis) source('functions/21b_generate_finite_basis_expansion.R')
  
  m     <- dim(G_hat_list[[1]])[1]
  Delta <- 1 / m
  
  # ----------------------------------------------------------------
  # Step A: Univariate FPCA on each diagonal block G_{ii}
  #         Reuses your existing compute_eigendecomposition_ii logic
  # ----------------------------------------------------------------
  
  # Build m x m x p array of diagonal blocks for the existing function
  G_diag_array <- array(0, dim = c(m, m, p))
  for (i in 1:p) {
    G_diag_array[,, i] <- G_hat_list[[paste0(i, '_', i)]]
  }
  
  uni <- compute_eigendecomposition_ii(
    G_hat        = G_diag_array,
    same_basis   = same_basis,
    constant_d   = constant_d,
    var_explained = var_explained
  )
  
  Phi    <- uni$eigenfunctions  # list of p: each (m x d_i)
  n_dims <- uni$n_dims          # list of p integers
  D      <- sum(unlist(n_dims))
  
  # Block offsets (1-indexed): block i lives in rows/cols (offsets[i]+1):offsets[i+1]
  offsets <- c(0L, cumsum(unlist(n_dims)))
  
  # ----------------------------------------------------------------
  # Step B: Build projected covariance matrix C (D x D)
  #   C[block_i, block_j] = Phi_i^T %*% G_{ij} %*% Phi_j * Delta^2
  # ----------------------------------------------------------------
  
  C <- matrix(0.0, nrow = D, ncol = D)
  
  for (i in 1:p) {
    ri <- (offsets[i] + 1L):offsets[i + 1L]
    for (j in i:p) {  # j starts at i, not 1
      rj  <- (offsets[j] + 1L):offsets[j + 1L]
      Gij <- G_hat_list[[paste0(i, '_', j)]]
      blk <- t(Phi[[i]]) %*% Gij %*% Phi[[j]] * Delta^2  # d_i x d_j
      C[ri, rj] <- blk
      C[rj, ri] <- t(blk)  # fill lower triangle via transpose
    }
  }
  
  C <- (C + t(C)) / 2  # enforce symmetry
  
  # ----------------------------------------------------------------
  # Step C: Eigendecompose C
  # ----------------------------------------------------------------
  
  eig_C  <- eigen(C, symmetric = TRUE)  # descending
  c_vals <- eig_C$values
  c_vecs <- eig_C$vectors               # D x D
  
  c_vals[c_vals < 0] <- 0
  
  # ----------------------------------------------------------------
  # Step D: Reconstruct multivariate eigenfunctions on the grid
  #   Psi_k^(i)(t) = Phi_i %*% c_vecs[block_i, k],  shape (m x D)
  # ----------------------------------------------------------------
  
  eigenfunctions <- vector("list", p)
  for (i in 1:p) {
    ri                <- (offsets[i] + 1L):offsets[i + 1L]
    eigenfunctions[[i]] <- Phi[[i]] %*% c_vecs[ri, , drop = FALSE]  # m x D
  }
  
  return(list(
    eigenvalues    = c_vals,
    eigenfunctions = eigenfunctions,
    univariate     = uni,
    n_dims         = n_dims,
    C              = C
  ))
}

