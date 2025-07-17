library(Matrix)
library(mvtnorm)
library(igraph)
library(MASS)

# these functions create the adjacency matrix ground truth

generate_random_graph <- function(p, sparsity, seed = NULL) {
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: Generate Erdős-Rényi random graph
  #
  # - any edge has probability p of occurring independently
  # 
  # Input: 
  #
  # - p (number of nodes), sparsity (edge probability)
  #
  # 
  # Output: 
  #
  # - adjacency matrix (p x p) 
  #
  #
  # ----------------------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  adj_matrix <- matrix(0, p, p)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      if (runif(1) < sparsity) {
        adj_matrix[i, j] <- adj_matrix[j, i] <- 1
      }
    }
  }
  return(adj_matrix)
}

generate_scale_free_graph <- function(p, sparsity, seed = NULL) {
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: Generate scale-free graph using preferential attachment
  #
  # - rich get richer. 
  #
  # 
  # Input: 
  #
  # - p (number of nodes), sparsity (edge probability)
  #
  # 
  # Output: 
  #
  # - adjacency matrix (p x p) 
  #
  #
  # ----------------------------------------------------------------------------  
  
  if (!is.null(seed)) set.seed(seed)
  
  # Target number of edges
  target_edges <- round(sparsity * p * (p-1) / 2)
  
  # Generate scale-free network
  g <- sample_pa(p, power = 1, m = max(1, target_edges %/% p), directed = FALSE)
  adj_matrix <- as.matrix(as_adjacency_matrix(g))
  
  return(adj_matrix)
}

generate_block_diagonal_graph <- function(p, sparsity, n_blocks = 3, seed = NULL) {
  
  

  # ----------------------------------------------------------------------------
  # 
  # GOAL: Generate block-diagonal graph structure
  #
  # - within each block, generate erdos renyi random graph
  #
  # 
  # Input: 
  #
  # - p        (integer)                   number of nodes
  # - sparsity (number between 0 and 1)    edge probability
  # - n_blocks (integer)                   number of blocks. Should be a factor of p
  # - seed     (integer)                   seed number
  # 
  # Output: 
  #
  # - adjacency matrix (p x p) 
  #
  #
  # ----------------------------------------------------------------------------    
  
  if (!is.null(seed)) set.seed(seed)
  
  adj_matrix <- matrix(0, p, p)
  block_size <- p %/% n_blocks # if not a factor, n_blocks plays it safe 
  
  for (b in 1:n_blocks) {
    start_idx <- (b-1) * block_size + 1
    end_idx <- min(b * block_size, p)
    block_indices <- start_idx:end_idx
    
    # Generate edges within block
    for (i in block_indices) {
      for (j in block_indices) {
        if (i < j && runif(1) < sparsity) {
          adj_matrix[i, j] <- adj_matrix[j, i] <- 1
        }
      }
    }
  }
  
  return(adj_matrix)
}

