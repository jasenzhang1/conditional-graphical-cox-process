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

estimate_graph <- function(P_conditional, threshold, p) {
  
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
  #
  # 
  # Output: 
  # - list with edges and edge_strengths
  #
  # ------------------------------------------------------------------------
  
  edges <- list()
  edge_strengths <- list()
  
  # Only consider upper triangle: i < j (undirected graph)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key_ij <- paste(i, j, sep="_")
      P_ij <- P_conditional[[key_ij]]  # m x m matrix
      
      # Compute Hilbert-Schmidt norm: m x m -> scalar
      hs_norm <- hilbert_schmidt_norm(P_ij)
      edge_strengths[[key_ij]] <- hs_norm
      
      # Apply threshold: scalar > scalar -> boolean
      if (hs_norm > threshold) {
        edges <- append(edges, list(c(i, j)))
      }
    }
  }
  
  return(list(edges = edges, edge_strengths = edge_strengths))
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