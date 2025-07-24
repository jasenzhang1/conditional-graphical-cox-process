
library(reshape2)
library(ggplot2)

visualize_precision_matrix <- function(precision_op){
  
  # precision_op = (p x p x m x m) matrix
  
  p <- dim(precision_op)[1]
  m <- dim(precision_op)[3]
  n <- p * m
  
  precision_mat <- assemble_block_matrix_24(precision_op)
  
  df <- reshape2::melt(precision_mat)
  colnames(df) <- c("row", "col", "value")
  
  heatmap_plot <- ggplot(df, aes(x = col, y = row, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +  # or use scale_fill_gradient2() for diverging color
    scale_y_reverse() +       # So origin is top-left like in a matrix
    coord_fixed() +           # Make cells square
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())
  
  # Compute block borders (at multiples of m)
  grid_lines <- seq(0, n, by = m)
  
  # Overlay lines
  for (g in grid_lines) {
    heatmap_plot <- heatmap_plot +
      geom_hline(yintercept = g + 0.5, color = "black", size = 0.3) +
      geom_vline(xintercept = g + 0.5, color = "black", size = 0.3)
  }  
  
  return(heatmap_plot)
  
}

HS_heatmap <- function(precision_op, delta_t){
  
  # precision_op = (p x p x m x m) matrix
  
  # recall that the precision matrix HS norm scales with delta_t^3
  HS_values <- apply(precision_op, c(1, 2), function(mat){sum(mat^2) * delta_t^3})
  
  HS_adj <- HS_values
  diag(HS_adj) <- 0
  HS_adj <- HS_adj > 0
  
  df <- reshape2::melt(HS_adj)
  colnames(df) <- c("row", "col", "value")
  df_nonzero <- df
  df_nonzero$value <- (df$value > 0)
  
  heatmap_plot <- ggplot(df, aes(x = col, y = row, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +  # or use scale_fill_gradient2() for diverging color
    scale_y_reverse() +       # So origin is top-left like in a matrix
    coord_fixed() +           # Make cells square
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank()) 
  
  heatmap_plot <- ggplot(df, aes(x = col, y = row, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue")) +
    scale_y_reverse() + 
    theme_minimal() +
    coord_fixed() +
    labs(title = "Logical Heatmap", fill = "Value")  
  
  
}


visualize_log_intensity <- function(X_k, time_grid){
  
  # X_k = (p x m matrix)
  # time_grid = (m-dim vec of timepoints)
  
  p <- dim(X_k)[1]
  m <- dim(X_k)[2]
  
  rownames(X_k) <- paste0(seq_len(p))
  
  # Convert to long format
  df <- reshape2::melt(X_k)
  colnames(df) <- c("Process", "TimeIndex", "Value")
  
  df$Time <- time_grid[df$TimeIndex]
  df$Process <- factor(df$Process)
  g <- ggplot(df, aes(x = Time, y = Value, color = Process)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Log Intensities", x = "Time", y = "Value") +
    theme(legend.position = "right")  
  
  return(g)
  
}


ground_truth_rho_ij <- function(my_list){
  
  # Assuming your list is called my_list
  # Each element has a matrix called X_functions of dimension p x m
  
  p <- nrow(my_list[[1]]$X_functions)
  m <- ncol(my_list[[1]]$X_functions)
  
  # Get all (i,j) index pairs for upper triangle including diagonal
  upper_pairs <- which(upper.tri(matrix(1, p, p), diag = TRUE), arr.ind = TRUE)
  
  # Initialize an empty list to store running sums of each outer product
  outer_sums <- vector("list", nrow(upper_pairs))
  names(outer_sums) <- paste0(upper_pairs[,1], "_", upper_pairs[,2])
  
  # Fill with zero matrices of dimension m x m
  outer_sums <- lapply(outer_sums, function(x) matrix(0, m, m))
  
  # Loop over list elements and accumulate outer products
  for (entry in my_list) {
    X <- exp(entry$X_functions)  # p x m
    
    for (k in seq_len(nrow(upper_pairs))) {
      i <- upper_pairs[k, 1]
      j <- upper_pairs[k, 2]
      
      v_i <- X[i, ]  # length m
      v_j <- X[j, ]
      
      outer_ij <- tcrossprod(v_i, v_j)  # m x m outer product
      outer_sums[[k]] <- outer_sums[[k]] + outer_ij
    }
  }
  
  # Average over n
  n <- length(my_list)
  outer_means <- lapply(outer_sums, function(mat) mat / n)
  
  return(outer_means)
}


block_matrix_HS <- function(pm_mat, p){
  
  #
  # GOAL: to check if HS norms of precision and covariace operators are the same 
  #       We have a pm x pm matrix, but want to take the HS norm of each m x m submatrix
  #
  #
  # input:
  #
  # - pm_matrix (pm x pm matrix)
  # - p         (number of processes)
  # - delta_t   (number to denote time spacing)
  #
  #
  
  m = dim(pm_mat)[1] / p
  HS_mat <- matrix(0, nrow = p, ncol = p)
  
  # Loop over block indices
  for (i in 1:p) {
    for (j in 1:p) {
      row_idx <- ((i - 1) * m + 1):(i * m)
      col_idx <- ((j - 1) * m + 1):(j * m)
      
      block <- pm_mat[row_idx, col_idx]
      HS_mat[i, j] <- sqrt(sum(block^2))   # Frobenius norm
    }
  }  
  
  return(HS_mat)
  
  
}
