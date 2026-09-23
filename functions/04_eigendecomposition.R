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

compute_eigendecomposition_ii_v2 <- function(G_hat_list, p, m,
                                             var_explained_diag = 0.95,
                                             var_explained_offdiag = 0.90,
                                             d_max = 20,
                                             max_iter = 100) {
  
  # --------------------------------------------------------------------------
  #
  # GOAL: select d_i for each process i such that:
  #   1. d_i explains var_explained_diag of G_{ii}          (marginal criterion)
  #   2. all off-diagonal blocks (i,j) have VE_{ij} >= var_explained_offdiag
  #
  # Algorithm:
  #   - initialise d_i from marginal variance explained on G_{ii}
  #   - repeatedly find the worst (i,j) off-diagonal block
  #   - increase d_i or d_j by 1, whichever gives the larger total VE gain
  #   - stop when all off-diagonal VE >= threshold, d_max is hit, or max_iter reached
  #
  # Input:
  #   - G_hat_list              named list of m x m matrices, keys "i_j" (i <= j)
  #   - p                       number of processes
  #   - m                       grid size
  #   - var_explained_diag      marginal threshold for initialising d_i
  #   - var_explained_offdiag   off-diagonal HS threshold to satisfy
  #   - d_max                   hard cap on any d_i
  #   - max_iter                iteration cap
  #
  # Output: list with
  #   - eigenvalues     (list of p vectors)
  #   - eigenfunctions  (list of p matrices, each m x d_i)
  #   - n_dims          (list of p integers)
  #   - diagnostics     list with:
  #       - ve_offdiag  (p x p) matrix of final off-diagonal VE values
  #       - history     data.frame logging each iteration
  #
  # --------------------------------------------------------------------------
  
  Delta <- 1 / m
  
  # ----------------------------------------------------------------
  # Helper: get eigenfunctions of G_{ii} up to d components
  # ----------------------------------------------------------------
  get_phi <- function(i, d) {
    Gii  <- G_hat_list[[paste0(i, "_", i)]]
    Gii  <- (Gii + t(Gii)) / 2
    eig  <- eigen(Gii * Delta, symmetric = TRUE)
    vecs <- eig$vectors[, 1:d, drop = FALSE]
    norms <- sqrt(colSums(vecs^2) * Delta)
    sweep(vecs, 2, norms, "/")
  }
  
  # ----------------------------------------------------------------
  # Helper: compute VE_{ij} for given Phi_i, Phi_j
  # ----------------------------------------------------------------
  compute_ve <- function(i, j, Phi_i, Phi_j) {
    key <- paste0(min(i,j), "_", max(i,j))
    Gij <- G_hat_list[[key]]
    if (i > j) Gij <- t(Gij)
    denom <- norm(Gij * Delta, type = "F")^2
    if (denom < .Machine$double.eps) return(1.0)
    num <- norm(t(Phi_i) %*% Gij %*% Phi_j * Delta^2, type = "F")^2
    num / denom
  }
  
  # ----------------------------------------------------------------
  # Step 1: initialise d_i from marginal variance explained
  # ----------------------------------------------------------------
  d <- integer(p)
  
  for (i in 1:p) {
    Gii  <- G_hat_list[[paste0(i, "_", i)]]
    Gii  <- (Gii + t(Gii)) / 2
    eig  <- eigen(Gii * Delta, symmetric = TRUE)
    lams <- pmax(eig$values, 0)
    cumve <- cumsum(lams) / sum(lams)
    d[i] <- min(which(cumve >= var_explained_diag)[1], d_max)
  }
  
  cat("Initial d:", d, "\n")
  
  # ----------------------------------------------------------------
  # Step 2: cache current eigenfunctions
  # ----------------------------------------------------------------
  Phi <- lapply(1:p, function(i) get_phi(i, d[i]))
  
  # ----------------------------------------------------------------
  # Step 3: greedy loop
  # ----------------------------------------------------------------
  history <- data.frame(iter      = integer(),
                        i         = integer(),
                        j         = integer(),
                        ve        = numeric(),
                        increased = integer(),
                        d_new     = integer())
  
  for (iter in 1:max_iter) {
    
    # Compute VE for all off-diagonal blocks
    ve_mat <- matrix(1.0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in (1:p)[-i]) {
        ve_mat[i, j] <- compute_ve(i, j, Phi[[i]], Phi[[j]])
      }
    }
    
    # Find worst off-diagonal block
    offdiag_mask <- row(ve_mat) != col(ve_mat)
    min_ve <- min(ve_mat[offdiag_mask])
    worst  <- which(ve_mat == min_ve & offdiag_mask, arr.ind = TRUE)[1, ]
    wi <- worst[1]; wj <- worst[2]
    
    cat(sprintf("Iter %d | worst block (%d,%d) VE = %.4f | d = [%s]\n",
                iter, wi, wj, min_ve, paste(d, collapse = ",")))
    
    # Check convergence
    if (min_ve >= var_explained_offdiag) {
      cat("All off-diagonal blocks satisfy threshold. Done.\n")
      break
    }
    
    # Try increasing d_i vs d_j, pick whichever gives largest total VE gain
    best_gain <- -Inf
    best_proc <- NA
    
    for (cand in c(wi, wj)) {
      if (d[cand] >= d_max) next
      
      d_cand   <- d[cand] + 1L
      Phi_cand <- get_phi(cand, d_cand)
      
      # Sum VE gain across all blocks involving cand
      gain <- 0
      for (k in (1:p)[-cand]) {
        if (cand < k) {
          ve_old <- compute_ve(cand, k, Phi[[cand]], Phi[[k]])
          ve_new <- compute_ve(cand, k, Phi_cand,   Phi[[k]])
        } else {
          ve_old <- compute_ve(k, cand, Phi[[k]], Phi[[cand]])
          ve_new <- compute_ve(k, cand, Phi[[k]], Phi_cand)
        }
        gain <- gain + (ve_new - ve_old)
      }
      
      if (gain > best_gain) {
        best_gain <- gain
        best_proc <- cand
      }
    }
    
    # If both at d_max, cannot improve — exit
    if (is.na(best_proc)) {
      cat(sprintf("Block (%d,%d) cannot be improved — d_%d=%d and d_%d=%d both at d_max=%d.\n",
                  wi, wj, wi, d[wi], wj, d[wj], d_max))
      break
    }
    
    # Apply update
    d[best_proc]    <- d[best_proc] + 1L
    Phi[[best_proc]] <- get_phi(best_proc, d[best_proc])
    
    history <- rbind(history, data.frame(iter      = iter,
                                         i         = wi,
                                         j         = wj,
                                         ve        = min_ve,
                                         increased = best_proc,
                                         d_new     = d[best_proc]))
  }
  
  # ----------------------------------------------------------------
  # Final outputs
  # ----------------------------------------------------------------
  
  # Final VE matrix
  ve_final <- matrix(1.0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in (1:p)[-i]) {
      ve_final[i, j] <- compute_ve(i, j, Phi[[i]], Phi[[j]])
    }
  }
  
  # Eigenvalues per process
  evals <- lapply(1:p, function(i) {
    Gii  <- G_hat_list[[paste0(i, "_", i)]]
    Gii  <- (Gii + t(Gii)) / 2
    eig  <- eigen(Gii * Delta, symmetric = TRUE)
    pmax(eig$values[1:d[i]], 0)
  })
  
  return(list(
    eigenvalues    = evals,
    eigenfunctions = Phi,
    n_dims         = as.list(d),
    diagnostics    = list(
      ve_offdiag = ve_final,
      history    = history
    )
  ))
}
