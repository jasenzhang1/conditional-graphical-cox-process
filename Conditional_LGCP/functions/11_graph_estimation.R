hilbert_schmidt_norm <- function(A) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: Compute Hilbert-Schmidt norm, which is just frobenius norm 
  #
  # - ||A||_{HS} = sqrt(sum of squares of elements)
  # 
  # 
  # Input: 
  #
  # - A       (m x m matrix)
  #
  # 
  # Output: 
  #
  # - HS norm (scalar)
  #
  # ------------------------------------------------------------------------
  
  return(sqrt(sum(A^2)))  # sqrt(sum of all squared elements)
}

list_to_adj_mat <- function(edge_list, p){
  
  # I have a list of edge strengths
  # Create a weighted adjacency matrix 
  #
  
  
  adj_mat <- matrix(0, p, p)
  
  for(k in names(edge_list)){
    
    i_j <- as.integer(strsplit(k, "_")[[1]])
    
    i <- i_j[1]
    j <- i_j[2]
    
    adj_mat[i, j] <- edge_list[[k]]
    adj_mat[j, i] <- edge_list[[k]]
    
  }
  
  return(adj_mat)
  
}

estimate_graph <- function(P_conditional, C_conditional, V_conditional, threshold, p, discarded_neurons = NULL) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: get edge estiamtes for our graph
  # 
  # 
  # Input: 
  #
  # - P_conditional   (list of length p^2, each element m x m)
  # - threshold       (scalar)
  # - p               (scalar)
  # - discarded_neurons (list of integers)
  #
  # 
  # Output: 
  # - list with edges and edge_strengths
  #
  # ------------------------------------------------------------------------
  
  edges <- list()
  edge_strengths <- list()
  
  
  adj_mat <- matrix(0, p, p)
  w_mat <- matrix(0, p, p)
  
  # Only consider upper triangle: i < j (undirected graph)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key_ij <- paste(i, j, sep="_")
      P_ij <- P_conditional[[key_ij]]  # m x m matrix
      
      
      
      # Compute Hilbert-Schmidt norm: m x m -> scalar
      hs_norm <- hilbert_schmidt_norm(P_ij)
      edge_strengths[[key_ij]] <- hs_norm
      w_mat[i,j] <- hs_norm
      w_mat[j,i] <- hs_norm
      
      # Apply threshold: scalar > scalar -> boolean
      if (hs_norm > threshold) {
        edges <- append(edges, list(c(i, j)))
        
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
      }
    }
  }
  
  m <- dim(P_ij)[1]
  
  return(list(edges = edges,                       # list of [[i_j]] entries that satisfy the threshold
              edge_strengths = edge_strengths,     # list of [[i_j]] weights
              w_mat = w_mat,
              adj_mat = adj_mat,
              P_conditional = P_conditional,       # original precision operator
              C_conditional = C_conditional,
              V_conditional = V_conditional,
              P_mat = assemble_block_matrix_v2(P_conditional, p, m),
              C_mat = assemble_block_matrix_v2(C_conditional, p, m),
              V_mat = assemble_block_matrix_v2(V_conditional, p, m),
              discarded_neurons = discarded_neurons))
}

select_threshold_by_stability <- function(P_conditional, p) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: how to choose edge selection threshold?
  # 
  # 
  # Input: 
  #
  # - P_conditional    (p^2 length list)
  # - p                (scalar)
  #
  # Output: 
  #
  # - threshold        (scalar)
  #
  # ------------------------------------------------------------------------
  
  # Collect all edge strengths
  edge_strengths <- c()
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key <- paste(i, j, sep="_")
      P_ij <- P_conditional[[key]]
      edge_strengths <- c(edge_strengths, hilbert_schmidt_norm(P_ij))
    }
  }
  
  # Use median as heuristic threshold
  return(if (length(edge_strengths) > 0) median(edge_strengths) else 0.1)
}



