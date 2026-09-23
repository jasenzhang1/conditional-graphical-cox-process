# functions that:
# 
# - construct adjacency matrices that can then be calculated
# - do the calculation


get_adj_matrix <- function(final_graph_estimates, thresh_value, n_neurons_not_discarded){
  
  # ----------------------------------------------------------------------------
  #
  # Goal: from the estimated object, construct all kinds of matrices 
  #
  # - adjacency matrix
  # - adjacency matrix with 1's on the diagonal
  # - weighted adjacency matrix
  # - weighted adjacency matrix with 1's on the diagonal
  #
  #
  # input:
  #
  # - final_graph_estimates      (list)
  #   - edges
  #   - edge_strengths
  #   - discarded_neurons
  # - thresh_value               (integer)
  # - n_neurons_not_discarded    (integer) number of neurons considered in the list
  #
  # 
  # output:
  # 
  # - adjacency matrix 
  # - adjacacency matrix with 1's on diag
  # - weighted matrix (with 0's on diag)
  #
  # ----------------------------------------------------------------------------
  
  
  # thresholding
  
  weights <- unlist(final_graph_estimates$edge_strengths)
  kept_indices <- weights > thresh_value
  kept_edges <- final_graph_estimates$edge_strengths[kept_indices] # list of items whose names are i_j
  
  
  # make adjacency matrix (adj_mat, adj_mat_w_diag) and weighted adjacency matrix (weight_mat)
  
  v_i <- sub("_.*", "", names(kept_edges)) %>% as.numeric()
  v_j <- sub(".*_", "", names(kept_edges)) %>% as.numeric()
  e_ij <- unlist(kept_edges)
  
  df_adj <- data.frame(v_i, v_j, e_ij)
  
  adj_mat <- matrix(0, nrow = n_neurons_not_discarded, ncol = n_neurons_not_discarded)
  weight_mat <- matrix(0, nrow = n_neurons_not_discarded, ncol = n_neurons_not_discarded)
  
  if(nrow(df_adj) > 0){
    for(k in 1:nrow(df_adj)){
      i <- df_adj[k,1]
      j <- df_adj[k,2]
      ij <- df_adj[k,3]
      adj_mat[i,j] <- 1
      adj_mat[j,i] <- 1
      weight_mat[i,j] <- ij
      weight_mat[j,i] <- ij
    }
  }
  
  adj_mat_w_diag <- adj_mat + diag(rep(1, n_neurons_not_discarded))
  
  return(list(adj_mat=adj_mat, 
              adj_mat_w_diag=adj_mat_w_diag, 
              weight_mat=weight_mat))
}
