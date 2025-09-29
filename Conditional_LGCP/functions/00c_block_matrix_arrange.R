# functions to help assemble and extract block matrices



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



extract_block_structure_v2 <- function(block_matrix, p, block_size) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: extract block sub-matrices from a block matrix
  #
  # - it's the reverse of assemble_block_matrix
  # - only want to keep the i_j entries for i <= j
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

extract_block_structure_ij <- function(block_matrix, block_size, i, j) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: extract block sub-matrices from a block matrix
  #
  # - we only want to extract a single i_j entry 
  #
  # 
  # Input: 
  #
  # - block_matrix   (pm x pm matrix)
  # - block_size     (scalar)  m
  # - i              (scalar)  row number
  # - j              (scalar)  column number
  #
  # 
  # Output: 
  #
  # -  the [i, j] block (m x m matrix)
  #
  # ------------------------------------------------------------------------
  
  
  row_start <- (i-1) * block_size + 1
  row_end <- i * block_size
  col_start <- (j-1) * block_size + 1
  col_end <- j * block_size
  
  
  # Extract m x m block from pm x pm matrix
  block_ij <- block_matrix[row_start:row_end, col_start:col_end]
  
  return(block_ij)
}
