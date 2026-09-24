# helpers
GIC_pseudo_logdet <- function(A, tol = 1e-8) {
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: pseudo logdet function in case our matrix is near singular
  #
  # 
  # inputs: 
  #
  # - A     (square matrix)
  # - tol   (value)           eigenvalue threshold when taking sum(log(lambda)) = logd
  #
  # ----------------------------------------------------------------------------
  
  # symmetrize to avoid numerical asymmetry
  A_sym <- (A + t(A))/2
  
  # eigenvalues
  ev <- eigen(A_sym, symmetric = TRUE, only.values = TRUE)$values
  
  # keep only sufficiently positive eigenvalues
  ev_pos <- ev[ev > tol]
  
  # handle fully singular case
  if (length(ev_pos) == 0) {
    return(-Inf)
  }
  
  return(sum(log(ev_pos)))
  
}

GIC_local_loss <- function(C_mat, Theta){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: calculate local loss component for GIC
  #
  #       L(Theta) = trace(C_mat %*% Theta) - log_det(Theta)
  #
  #
  # Input: 
  #
  # - C_mat        (pd x pd matrix)   correlation operator before any thresholding
  # - Theta        (pd x pd matrix)   precision operator after thresholding for tau_c and tau_p
  #
  # 
  # Output: 
  #
  # - local_loss   (scalar)  local loss evaluation
  #
  # ----------------------------------------------------------------------------
  
  Theta <- sym(Theta)
  
  # Check dimensions
  if (!all(dim(C_mat) == dim(Theta))) {
    stop("10b GIC_local_loss: C_mat and Theta must have the same dimensions.")
  }
  
  # Ensure Theta is square
  if (nrow(Theta) != ncol(Theta)) {
    stop("10b GIC_local_loss: Theta must be a square matrix.")
  }
  
  # Compute trace safely
  Xt <- C_mat %*% Theta
  if (anyNA(Xt) || any(!is.finite(Xt))) {
    stop("10b GIC_local_loss: Non-finite values encountered in C_mat %*% Theta.")
  }
  
  tr_val <- sum(diag(Xt))
  
  log_det_val <- GIC_pseudo_logdet(Theta)
  
  # Check for invalid log det
  if (!is.finite(log_det_val)) {
    return(Inf)
  }
  
  # Final check
  local_loss <- tr_val - log_det_val
  
  # if (!is.finite(local_loss)) {
  #   stop("10b GIC_local_loss: Local loss is not finite.")
  # }
  
  return(local_loss)
  
}


GIC_set_CXX_diag_identity <- function(C_cond){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: set C_{Xi Xi} blocks to be identity blocks
  #
  #
  # input:
  #
  # - C_cond    (list of i_j matrices)
  #
  # output:
  #
  # - C_cond2  (list of i_j matrices)
  # 
  # ----------------------------------------------------------------------------
  
  C_cond2 <- lapply(seq_along(C_cond), function(idx) {
    block <- C_cond[[idx]]
    nm    <- names(C_cond)[idx]
    
    # name must be "i_j"
    parts <- strsplit(nm, "_")[[1]]
    if (length(parts) != 2)
      stop("Each list name must be of form 'i_j'")
    
    i <- as.numeric(parts[1])
    j <- as.numeric(parts[2])
    
    # If this is a diagonal block, set to identity
    if (i == j) {
      nr <- nrow(block)
      nc <- ncol(block)
      
      # sanity check
      if (nr != nc)
        stop(sprintf("Block '%s' is not square (%d x %d)", nm, nr, nc))
      
      return(diag(1, nr))
    }
    
    # otherwise return unchanged
    return(block)
  })
  
  names(C_cond2) <- names(C_cond)
  return(C_cond2)
  
}



GIC_theta_to_adj <- function(M_list, p){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: Given a full theta block matrix, return a pxp adjacency matrix 
  #       along with an nx2 matrix detailing which vertices are connected
  #
  #
  # input:
  # 
  # - M_list  (list of i_j entries of d_i x d_j matrices) 
  # - p       (integer)
  #
  # 
  # output:
  # 
  # list of the following:
  #
  # - adj_list   (n x 2 matrix)   matrix of all (i, j) edges 
  # - adj_mat    (p x p matrix)   adjacency matrix where diag entries are 0
  #
  # ----------------------------------------------------------------------------
  
  # 0) zero out the diagonal matrices
  
  diag_indices <- paste0(1:p, '_', 1:p)
  
  # 3) zero out select matrices
  for (idx in diag_indices) {
    block <- M_list[[idx]]
    M_list[[idx]] <- matrix(0, nrow(block), ncol(block))
  }  
  
  # 1) Compute HS norms for each matrix on the list
  hs_vals <- sapply(M_list, hilbert_schmidt_norm)
  
  # 2) Get names of matrices with HS norm > 0
  names_nonzero <- names(M_list)[hs_vals > 0] 
  
  # 3) make the nx2 matrix of adjacencies
  if(length(names_nonzero) > 0){
    tmp <- strsplit(names_nonzero, "_")
    mat <- do.call(rbind, tmp)
    adj_list <- apply(mat, 2, as.numeric)
    
    # edge case where adj_list becomes a vector if dim = 1
    if(dim(mat)[1] == 1){
      adj_list <- matrix(adj_list, nrow = 1)
    }
    
    # 4) make the pxp matrix of adjacencies
    adj_mat <- matrix(0, p, p)
    for(k in 1:nrow(adj_list)){
      i <- adj_list[k, 1]
      j <- adj_list[k, 2]
      
      adj_mat[i, j] <- 1
      adj_mat[j, i] <- 1
    }
  } else{
    adj_list <- data.frame()
    adj_mat <- matrix(0, p, p)
  }
  
  
  # 5) return
  return(list(adj_mat = adj_mat,
              adj_list = adj_list))
}

GIC_evalulation <- function(C_cond, Theta_cond_thresh, W_y, num_edges){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: evalulate GIC for a specific value of tau_c, tau_p
  #
  #
  # Input: 
  #
  # - C_cond             (pd x pd matrix)
  # - Theta_cond_thresh  (pd x pd matrix)
  # - W_y                (scalar)           effective sample size
  # - num_edges          (scalar)
  #
  # 
  # Output: 
  #
  # - GIC                (scalar)
  #
  # ----------------------------------------------------------------------------
  

  GIC <- W_y * GIC_local_loss(C_cond, Theta_cond_thresh) + sqrt(W_y) * num_edges
  
  return(GIC)
  
}

GIC_edge_count <- function(M_list, p){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Given a list of matrices of the form i_j, find the edge count (Excluding diagonals)
  #
  #
  # inputs:
  #
  # - M_list (list of i_j entries that are d_i x d_j matrices)
  # - p      (integer, number of processes) 
  # 
  #
  # outputs:
  #
  # - num_edges (integer)
  #
  # 
  # ----------------------------------------------------------------------------
  
  hs_vals <- sapply(M_list, hilbert_schmidt_norm)
  
  # Parse names into i and j
  nm <- names(M_list)
  parts <- do.call(rbind, strsplit(nm, "_"))
  i_idx <- as.numeric(parts[,1])
  j_idx <- as.numeric(parts[,2])
  
  # Logical mask: off-diagonal & HS norm > 0
  num_edges <- sum((i_idx != j_idx) & (hs_vals > 0))
  
  return(num_edges)
  
}



GIC_get_percentile_info <- function(current_list, indices) {
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: Function to get threshold info based on percentiles
  #
  #
  # inputs:
  #
  # - current_list    (i_j list of any sized matrix, d_i x d_j or m x m)
  # - indices         (vector) vector if idx's of the off-diagonals of interest       
  #
  #
  # outputs:
  #
  # - hs_vals        (vector)    vector of thresholds that we chose
  # - index_path     (list)      list of off-diagonal values that we want to zero out for the i-th quantile
  #
  # ----------------------------------------------------------------------------
  
  # Calculate HS norms for off-diagonal blocks only
  hs_norms <- sapply(indices, function(idx) norm(current_list[[idx]], type = "F"))
  
  # Create 101 percentiles (0, 0.01, ..., 1.00)
  step_size = 1/100
  probs <- seq(0, 1, by = step_size)
  raw_quantiles <- quantile(hs_norms, probs = probs)
  
  # 3) CRITICAL: Keep only unique threshold values
  # This prevents "subscript out of bounds" when you have few blocks
  tau_levels <- unique(raw_quantiles)
  
  # If the smallest threshold is > 0, prepend 0 to allow a 'keep all' state
  if(min(tau_levels) > 0) {
    tau_levels <- c(0, tau_levels)
  }
  
  # 4) Map the path based on unique thresholds
  path <- lapply(tau_levels, function(t) indices[hs_norms <= t])
  
  return(list(hs_vals = as.numeric(tau_levels), index_path = path))
}


