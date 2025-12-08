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
  #                     (10/1/2025, we allow i_i entries only and assemble them accordingly)
  # - p                 (scalar)
  # - block_size        (scalar, equals m)
  #
  # 
  # Output: 
  #
  # - block_matrix      (pm x pm matrix)
  #
  # ------------------------------------------------------------------------
  
  check1 <- all(sapply(my_list, function(x) {
    is.matrix(x) && all(dim(x) == c(block_size, block_size))
  }))
  
  if(! check1){
    stop('STOP: assemble_block_matrix_v2: not a perfect square block matrix')
  }
  
  key_mat <- do.call(rbind, strsplit(names(operator_list), "_"))
  key_mat <- apply(key_mat, 2, as.numeric)
  
  total_size <- p * block_size
  block_matrix <- matrix(0, nrow=total_size, ncol=total_size)
  
  for(k in 1:nrow(key_mat)){
    
    i <- key_mat[k,1]
    j <- key_mat[k,2]
    
    # Calculate block indices
    row_start <- (i-1) * block_size + 1  # (i-1)*m + 1
    row_end <- i * block_size             # i*m
    col_start <- (j-1) * block_size + 1   # (j-1)*m + 1  
    col_end <- j * block_size             # j*m
    
    key <- paste(min(i, j), max(i, j), sep="_")
    
    if(key %in% names(operator_list)){
      # Insert m x m block into pm x pm matrix
      
      block_matrix[row_start:row_end, col_start:col_end] <- operator_list[[key]]
      if(i != j){
        block_matrix[col_start:col_end, row_start:row_end] <- t(operator_list[[key]])
      }
    }
  }

  
  return(block_matrix)
}

assemble_block_matrix_irregular <- function(operator_list, p) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Assemble a block matrix from a list of blocks operator_list
  # 
  # - Only contains i_j with i <= j. Off-diagonal blocks can be rectangular
  # - j_i entries are automatically set as transpose of i_j
  #
  #
  # input:
  # 
  # - operator_list  (list of entries named i_j)
  # - p              (integer)                      number of processes
  #
  #
  # output:
  #
  # - list of the following:
  #
  #   - block_matrix  (square pd x pd-ish matrix)
  #   - row_borders   (last index of the respective row block)
  #   - col_borders   (last index of the respective column block)
  # 
  # ----------------------------------------------------------------------------
  
  # parse names into i,j
  key_mat <- do.call(rbind, strsplit(names(operator_list), "_"))
  key_mat <- apply(key_mat, 2, as.numeric)
  
  # determine total number of rows and columns
  # by summing the rows of each diagonal block (i==j) and columns of each diagonal block
  row_sizes <- col_sizes <- numeric(p)
  for (k in 1:nrow(key_mat)) {
    i <- key_mat[k,1]
    j <- key_mat[k,2]
    block <- operator_list[[paste(i, j, sep="_")]]
    nr <- nrow(block)
    nc <- ncol(block)
    if (i == j) {
      row_sizes[i] <- nr
      col_sizes[j] <- nc
    }
  }
  
  total_rows <- sum(row_sizes)
  total_cols <- sum(col_sizes)
  
  block_matrix <- matrix(0, nrow = total_rows, ncol = total_cols)
  
  # compute row/column starting positions
  row_starts <- cumsum(c(0, row_sizes[-p])) + 1
  row_ends   <- cumsum(row_sizes)
  col_starts <- cumsum(c(0, col_sizes[-p])) + 1
  col_ends   <- cumsum(col_sizes)
  
  # fill in blocks
  for (k in 1:nrow(key_mat)) {
    i <- key_mat[k,1]
    j <- key_mat[k,2]
    block <- operator_list[[paste(i,j,sep="_")]]
    
    # indices in final matrix
    rs <- row_starts[i]; re <- row_ends[i]
    cs <- col_starts[j]; ce <- col_ends[j]
    
    # place block
    block_matrix[rs:re, cs:ce] <- block
    
    # for off-diagonal, fill j_i as transpose
    if (i != j) {
      block_matrix[cs:ce, rs:re] <- t(block)
    }
  }
  
  return(list(
    block_matrix = block_matrix,
    row_borders = row_ends,
    col_borders = col_ends
  ))
}

extract_block_matrix_irregular <- function(full_matrix, row_borders, col_borders) {
  
  
  # ----------------------------------------------------------------------------
  # GOAL:
  # Given a full matrix and row/column borders, extract blocks as a list.
  #
  # Inputs:
  # - full_matrix  (pm x qm matrix)
  # - row_borders: (vector of cumulative row ends, e.g. row_ends)
  # - col_borders: (vector of cumulative col ends, e.g. col_ends)
  #
  # Outputs:
  # - block_list: list of i_j blocks (names "i_j")
  # ----------------------------------------------------------------------------
  
  p <- length(row_borders)
  q <- length(col_borders)
  
  # compute row/col starts
  row_starts <- c(1, row_borders[-p] + 1)
  col_starts <- c(1, col_borders[-q] + 1)
  
  block_list <- list()
  
  for (i in 1:p) {
    for (j in i:q) {  # only upper triangular blocks (i <= j)
      rs <- row_starts[i]
      re <- row_borders[i]
      cs <- col_starts[j]
      ce <- col_borders[j]
      
      block <- full_matrix[rs:re, cs:ce, drop=FALSE]
      block_list[[paste(i,j,sep="_")]] <- block
    }
  }
  
  return(block_list)
}

extract_block_structure_v2 <- function(block_matrix, p, block_size) {
  
  # ----------------------------------------------------------------------------
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
  # - block_size     (scalar) m 
  #
  # 
  # Output: 
  #
  # - operator_list (list of length p + pC2, each element block_size x block_size)
  #
  # ----------------------------------------------------------------------------
  
  if(dim(block_matrix)[1] != dim(block_matrix)[2]){
    stop('STOP: extract_block_structure_v2: input is not a square matrix')
  }
  
  if(p * block_size != dim(block_matrix)[1]){
    stop('STOP: extract_block_structure_v2: parameter sizes do not match')
  }
  
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

assemble_blockwise_correlation <- function(Sigma, p, d){
  
  
  # ----------------------------------------------------------------------------
  # 
  # 
  # GOAL: for a pd x pd covariance matrix, return blocks that are correlation matrices
  #
  # - [Sigma]_{ij} is the (i,j)th block
  # - Let [R]_{ij} = diag(Sigma_{ii})^{-1/2} %*% Sigma_{ij} %*% diag(Sigma_{jj})^{-1/2}
  # - Then block them into pd x pd matrix R
  #
  #
  # inputs:
  #
  # - Sigma        (pd x pd matrix)  block covariance matrix
  # - p            (scalar)
  # - d            (scalar)
  #
  # output:
  #
  # - R            (pd x pd matrix)   block correlation matrix
  #
  # ----------------------------------------------------------------------------
  

  R <- matrix(0, nrow = p * d, ncol = p * d)
  
  # Precompute inverse sqrt of diagonal blocks
  D_inv_sqrt_list <- vector("list", p)
  for (i in 1:p) {
    idx_i <- ((i - 1) * d + 1):(i * d)
    Sigma_ii <- Sigma[idx_i, idx_i]
    D_inv_sqrt_list[[i]] <- diag(1 / sqrt(diag(Sigma_ii)))
  }
  
  # Fill in correlation blocks
  for (i in 1:p) {
    for (j in 1:p) {
      idx_i <- ((i - 1) * d + 1):(i * d)
      idx_j <- ((j - 1) * d + 1):(j * d)
      Sigma_ij <- Sigma[idx_i, idx_j]
      R[idx_i, idx_j] <- D_inv_sqrt_list[[i]] %*% Sigma_ij %*% D_inv_sqrt_list[[j]]
    }
  }
  
  return(R)
  
}
