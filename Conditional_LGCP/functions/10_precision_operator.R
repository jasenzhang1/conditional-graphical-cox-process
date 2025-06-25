assemble_block_matrix <- function(operator_list, p, block_size) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: rearrange correlation list into a pm x pm matrix
  #
  #
  # Input: 
  #
  # - operator_list     (list of length p^2, each element block_size x block_size)
  # - p                 (scalar)
  # - block_size        (scalar, equals m)
  #
  # 
  # Output: 
  #
  # - block_matrix      (pm x pm matrix)
  #
  # ------------------------------------------------------------------------
  
  total_size <- p * block_size
  block_matrix <- matrix(0, nrow=total_size, ncol=total_size)
  
  for (i in 1:p) {
    for (j in 1:p) {
      # Calculate block indices
      row_start <- (i-1) * block_size + 1  # (i-1)*m + 1
      row_end <- i * block_size             # i*m
      col_start <- (j-1) * block_size + 1   # (j-1)*m + 1  
      col_end <- j * block_size             # j*m
      
      key <- paste(i, j, sep="_")
      # Insert m x m block into pm x pm matrix
      block_matrix[row_start:row_end, col_start:col_end] <- operator_list[[key]]
    }
  }
  
  return(block_matrix)
}

assemble_block_matrix_v2 <- function(operator_list, p, block_size) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: rearrange correlation list into a pm x pm matrix
  #
  # - v2: operator list only has j >= i entries
  # - we assume the i_j entry that is ommitted is the transpose of the existing one
  #
  #
  # Input: 
  #
  # - operator_list     (list of length p + pC2, each element is a matrix of block_size x block_size)
  # - p                 (scalar)
  # - block_size        (scalar, equals m)
  #
  # 
  # Output: 
  #
  # - block_matrix      (pm x pm matrix)
  #
  # ------------------------------------------------------------------------
  
  total_size <- p * block_size
  block_matrix <- matrix(0, nrow=total_size, ncol=total_size)
  
  for (i in 1:p) {
    for (j in i:p) {
      # Calculate block indices
      row_start <- (i-1) * block_size + 1  # (i-1)*m + 1
      row_end <- i * block_size             # i*m
      col_start <- (j-1) * block_size + 1   # (j-1)*m + 1  
      col_end <- j * block_size             # j*m
      
      key <- paste(min(i, j), max(i, j), sep="_")
      
      # Insert m x m block into pm x pm matrix

      block_matrix[row_start:row_end, col_start:col_end] <- operator_list[[key]]
      if(i != j){
        block_matrix[col_start:col_end, row_start:row_end] <- t(operator_list[[key]])
      }
    }
  }
  
  return(block_matrix)
}

extract_block_structure <- function(block_matrix, p, block_size) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: extract block sub-matrices from a block matrix
  #
  # - it's the reverse of assemble_block_matrix
  #
  # 
  # Input: 
  #
  # - block_matrix   (pm x pm matrix)
  # - p              (scalar)
  # - block_size     (scalar)
  #
  # 
  # Output: 
  #
  # - operator_list (list of length p^2, each element block_size x block_size)
  #
  # ------------------------------------------------------------------------
  
  operator_list <- list()
  
  for (i in 1:p) {
    for (j in 1:p) {
      row_start <- (i-1) * block_size + 1
      row_end <- i * block_size
      col_start <- (j-1) * block_size + 1
      col_end <- j * block_size
      
      key <- paste(i, j, sep="_")
      # Extract m x m block from pm x pm matrix
      operator_list[[key]] <- block_matrix[row_start:row_end, col_start:col_end]
    }
  }
  
  return(operator_list)
}

extract_block_structure_v2 <- function(block_matrix, p, block_size) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: extract block sub-matrices from a block matrix
  #
  # - it's the reverse of assemble_block_matrix
  #
  # 
  # Input: 
  #
  # - block_matrix   (pm x pm matrix)
  # - p              (scalar)
  # - block_size     (scalar)
  #
  # 
  # Output: 
  #
  # - operator_list (list of length p + pC2, each element block_size x block_size)
  #
  # ------------------------------------------------------------------------
  
  operator_list <- list()
  
  for (i in 1:p) {
    for (j in i:p) {
      row_start <- (i-1) * block_size + 1
      row_end <- i * block_size
      col_start <- (j-1) * block_size + 1
      col_end <- j * block_size
      
      key <- paste(i, j, sep="_")
      # Extract m x m block from pm x pm matrix
      operator_list[[key]] <- block_matrix[row_start:row_end, col_start:col_end]
    }
  }
  
  return(operator_list)
}

estimate_precision_operator <- function(C_conditional, gamma2, p) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: 
  #
  #
  # Input: 
  #
  # - C_conditional (list of length p^2, each element m x m)
  # - gamma2 (scalar)
  # - p (scalar)
  #
  # 
  # Output: 
  #
  # - P_conditional (list of length p^2, each element m x m)
  #
  # ------------------------------------------------------------------------
  
  n_time <- nrow(C_conditional[[1]])  # m
  
  # Assemble block correlation matrix: pm x pm
  C_block <- assemble_block_matrix(C_conditional, p, n_time)
  
  # Add regularization: (pm x pm) + (pm x pm) = (pm x pm)
  C_block_reg <- C_block + gamma2 * diag(p * n_time)
  
  # Compute precision operator: (pm x pm)^{-1} = (pm x pm)
  P_block <- solve(C_block_reg)
  
  # Extract block structure: pm x pm -> list of p^2 blocks (m x m each)
  P_conditional <- extract_block_structure(P_block, p, n_time)
  
  return(P_conditional)
}

estimate_precision_operator_v2 <- function(C_conditional, gamma2, p) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: 
  #
  #
  # Input: 
  #
  # - C_conditional (list of length p^2, each element m x m)
  # - gamma2 (scalar)
  # - p (scalar)
  #
  # 
  # Output: 
  #
  # - P_conditional (list of length p^2, each element m x m)
  #
  # ------------------------------------------------------------------------
  
  n_time <- nrow(C_conditional[[1]])  # m
  
  # Assemble block correlation matrix: pm x pm
  C_block <- assemble_block_matrix_v2(C_conditional, p, n_time)
  
  # Add regularization: (pm x pm) + (pm x pm) = (pm x pm)
  C_block_reg <- C_block + gamma2 * diag(p * n_time)
  
  # Compute precision operator: (pm x pm)^{-1} = (pm x pm)
  P_block <- solve(C_block_reg)
  
  # Extract block structure: pm x pm -> list of p^2 blocks (m x m each)
  P_conditional <- extract_block_structure_v2(P_block, p, n_time)
  
  return(P_conditional)
}
