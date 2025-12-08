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
  
  sum(log(ev_pos))
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
    stop("10b GIC_local_loss: log_det(Theta) is not finite (matrix may be singular or not positive definite).")
  }
  
  # Final check
  local_loss <- tr_val - log_det_val
  
  if (!is.finite(local_loss)) {
    stop("10b GIC_local_loss: Local loss is not finite.")
  }
  
  return(local_loss)
  
}

GIC_effective_sample_size <- function(Y_c_k, y_c_query){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: get effective sample size wrt y_c_query
  #
  # inputs:
  # 
  # - Y_c_k        (n x q_c dim matrix)
  # - y_c_query    (q_c-dim vector)
  #
  #
  # outputs:
  #
  # - eff_ss       (scalar)
  #
  #
  # ----------------------------------------------------------------------------
  
  # 1) get gamma 
  
  gamma_c <- select_gamma_c_bandwidth_v2(Y_c_k)
  
  weights <- apply(Y_c_k, 1, function(row) {
    step_6_kernel(as.numeric(row), y_c_query, gamma_c) 
  })    
  
  eff_ss <- sum(weights)
  
  return(eff_ss)
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

GIC_threshold_block_matrix <- function(M, thresh) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Threshold a block matrix by zeroing out i_j entries whose HS norms are < thresh, AND i =/= j
  #
  # 
  # inputs:
  #
  # - M        (list of i_j entries)
  # - thresh   (scalar)
  # 
  # 
  # outputs:
  #
  # - M2       (list of i_j entries) some of them are zeroed out 
  # 
  # ----------------------------------------------------------------------------
  
  M2 <- lapply(seq_along(M), function(idx) {
    block <- M[[idx]]
    nm    <- names(M)[idx]
    
    # Parse name "i_j"
    parts <- strsplit(nm, "_")[[1]]
    
    # Ensure it has two numeric parts
    if (length(parts) == 2) {
      i <- as.numeric(parts[1])
      j <- as.numeric(parts[2])
    } else {
      stop("List element names must be of form 'i_j'")
    }
    
    # Only threshold if i != j
    if (i != j) {
      hs <- hilbert_schmidt_norm(block)
      if (hs < thresh) {
        return(matrix(0, nrow(block), ncol(block)))
      }
    }
    
    # Otherwise return block unchanged
    return(block)
  })
  
  names(M2) <- names(M)   # preserve names
  return(M2)

}

GIC_get_thresholds <- function(M) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: get thresholds that go from zeroing nothing to zeroing all off-diagonals
  #
  #
  # inputs:
  #
  # - M     (list of i_j matrices)
  #
  # outputs:
  #
  # - list  (with the following entrys)
  # 
  #   - index_path  (list)      index of the i_j matrix that should be zeroed out
  #   - hs_vals     (vector)    threshold that zeros out blocks whose HS norms are <= this value
  # 
  # ----------------------------------------------------------------------------
  
  hs_vals   <- c()       # HS norms
  hs_index  <- c()       # corresponding list indices
  
  # 1) collect HS norms for off-diagonal blocks only
  for (idx in seq_along(M)) {
    nm <- names(M)[idx]
    block <- M[[idx]]
    
    parts <- strsplit(nm, "_")[[1]]
    if (length(parts) != 2)
      stop("Each name must be 'i_j'")
    
    i <- as.numeric(parts[1])
    j <- as.numeric(parts[2])
    
    if (i != j) {
      hs_vals  <- c(hs_vals, hilbert_schmidt_norm(block))
      hs_index <- c(hs_index, idx)
    }
  }
  
  # 2) sorting
  
  # sort by HS norm
  order_idx <- order(hs_vals)
  
  # sorted HS norms and their corresponding indices
  hs_sorted   <- hs_vals[order_idx]
  index_sorted <- hs_index[order_idx]
  
  # 3) prepend a dummy threshold (-1) for "no zeroing"
  hs_vals_aug      <- c(-1, hs_sorted)
  
  # 4) cumulative index path
  K <- length(index_sorted)
  index_path_aug <- vector("list", K + 1)
  
  index_path_aug[[1]] <- integer(0)  # no zeroing
  
  for (k in 1:K) {
    index_path_aug[[k + 1]] <- c(index_path_aug[[k]], index_sorted[k])
  }
  
  # 5) return
  list(
    hs_vals      = hs_vals_aug,      # augmented with -1
    index_path   = index_path_aug    # first = no zeroing index
  )
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


# final function used in 12z
GIC_algorithm <- function(C_cond, p, W_y){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: gridsearch across (tau_c, tau_p) for minimum GIC
  #
  #
  # Input: 
  #
  # - C_cond        (pd x pd matrix)
  # - p             (scalar)
  # - W_y           (scalar)           effective sample size
  #
  # 
  # Output: 
  #
  # - list of the following metrics:
  #
  #   - tau_c   (scalar)    threshold where if your HS norm is less than this, you are zeroed out 
  #   - tau_p   (scalar) 
  #
  # ----------------------------------------------------------------------------
  
  # 0) make diagonal entries identity + prep
  C_cond <- GIC_set_CXX_diag_identity(C_cond)
  best_k <- NA
  best_l <- NA
  lowest_GIC <- NA
  df_GIC <- data.frame(C_thresh_index = integer(0),
                       Theta_thresh_index = integer(0),
                       GIC = numeric(0))

  # 1) prepare for grid search: 
  #    get relevant threshold values for C_cond that zero out the off-diagonals
  
  threshold_list <- GIC_get_thresholds(C_cond)
  
  
  
  # 2) begin grid search on C_cond
  
  for(k in 1:length(threshold_list$index_path)){
    
    excluded_indices <- threshold_list$index_path[[k]]
    tau_c <- threshold_list$hs_vals[k]
    C_cond_thresh <- C_cond
    
    # 3) zero out select matrices
    for (idx in excluded_indices) {
      block <- C_cond_thresh[[idx]]
      C_cond_thresh[[idx]] <- matrix(0, nrow(block), ncol(block))
    }
    
    # 4) assemble this pd x pd matrix whose entries may not be all squares
    
    C_cond_full <- assemble_block_matrix_irregular(C_cond_thresh, p)
    
    # 5) MP inverse - uses MASS package
    Theta_full <- ginv(C_cond_full$block_matrix)
    
    # 6) make into a list again
    Theta_cond <- extract_block_matrix_irregular(Theta_full, C_cond_full$row_borders, C_cond_full$col_borders)
    
    # 7) thresholding for tau_p
    threshold_list_2 <- GIC_get_thresholds(Theta_cond)
    
    for(l in 1:length(threshold_list_2$index_path)){
      
      excluded_indies_2 <- threshold_list_2$index_path[[l]]
      tau_p <- threshold_list_2$hs_vals[l]
      Theta_cond_thresh <- Theta_cond
      
      # 8) zero out select matrices
      for (idx in excluded_indies_2) {
        block <- Theta_cond_thresh[[idx]]
        Theta_cond_thresh[[idx]] <- matrix(0, nrow(block), ncol(block))
      }
      
      # 9) evalulate GIC
      
      Theta_cond_full <- assemble_block_matrix_irregular(Theta_cond_thresh, p)
      num_edges <- GIC_edge_count(Theta_cond_thresh, p)
      
      GIC <- GIC_evalulation(C_cond_full$block_matrix, Theta_cond_full$block_matrix, W_y, num_edges)
      
      # 10) update best GIC and indices
      
      df_kl <- data.frame(C_thresh_index = k,
                          Theta_thresh_index = l,
                          GIC = GIC)
      
      df_GIC <- rbind(df_GIC, df_kl)
      
      if(is.na(lowest_GIC)){
        lowest_GIC <- GIC
        best_k <- k
        best_l <- l
      } else if(GIC < lowest_GIC){
        lowest_GIC <- GIC
        best_k <- k
        best_l <- l
      }
    }
  } # end grid search
  
  
  
  # ------------------------------------------------
  # repeat but use best settings
  # ------------------------------------------------
  
  # 2) obtain best adjacency outcome
  
  excluded_indices <- threshold_list$index_path[[best_k]]
  tau_c_final <- threshold_list$hs_vals[best_k]
  C_cond_thresh_final <- C_cond
  
  
  # 3) zero out select matrices
  for (idx in excluded_indices) {
    block <- C_cond_thresh_final[[idx]]
    C_cond_thresh_final[[idx]] <- matrix(0, nrow(block), ncol(block))
  }
  
  # 4) assemble this pd x pd matrix whose entries may not be all squares
  
  C_cond_full_final <- assemble_block_matrix_irregular(C_cond_thresh_final, p)
  
  # 5) MP inverse - uses MASS package
  Theta_full_final <- ginv(C_cond_full_final$block_matrix)
  
  # 6) make into a list again
  Theta_cond_final <- extract_block_matrix_irregular(Theta_full_final, C_cond_full_final$row_borders, C_cond_full_final$col_borders)
  
  # 7) thresholding for tau_p
  threshold_list_2 <- GIC_get_thresholds(Theta_cond_final)
  
  # use best_l
  excluded_indies_2 <- threshold_list_2$index_path[[best_l]]
  tau_p_final <- threshold_list_2$hs_vals[best_l]
  Theta_cond_thresh_final <- Theta_cond_final
  
  # 8) zero out select matrices
  for (idx in excluded_indies_2) {
    block <- Theta_cond_thresh_final[[idx]]
    Theta_cond_thresh_final[[idx]] <- matrix(0, nrow(block), ncol(block))
  }
  
  # 9) evalulate GIC + get number of edges
  
  Theta_cond_full_final <- assemble_block_matrix_irregular(Theta_cond_thresh_final, p)
  num_edges_final <- GIC_edge_count(Theta_cond_thresh_final, p)
  GIC_final <- GIC_evalulation(C_cond_full_final$block_matrix, Theta_cond_full_final$block_matrix, W_y, num_edges_final)
  
  # 10) oreganize and return
  
  adj_results <- GIC_theta_to_adj(Theta_cond_thresh_final, p)
  
  
  
  result <- list(# thresholds
                 tau_c = tau_c_final,                      
                 tau_p = tau_p_final,
                 tau_c_levels = threshold_list$hs_vals,
                 tau_p_levels = threshold_list_2$hs_vals,
                 # operators in list form
                 C_cond = C_cond_thresh_final,                    
                 Theta_cond = Theta_cond_thresh_final,
                 # operators in matrix form
                 C_cond_full = C_cond_full_final$block_matrix,
                 Theta_cond_full = Theta_cond_full_final$block_matrix,
                 # adjacencies
                 adj_mat = adj_results$adj_mat,
                 adj_list = adj_results$adj_list,
                 num_edges = num_edges_final,
                 # troubleshooting
                 df_GIC = df_GIC
                 )
  
  return(result)
  
}
