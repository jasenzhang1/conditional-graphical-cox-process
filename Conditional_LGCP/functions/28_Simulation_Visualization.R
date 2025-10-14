
library(reshape2)
library(ggplot2)
#library(gganimate)
#library(magick)

visualize_pm_block_matrix_heatmap <- function(pm_block_matrix, g_title = NULL){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: visualize the pm block matrix in a heatmap 
  #
  # - force midpoint to be white = 0
  # 
  # input:
  #
  # - pm_block_matrix (pm x pm matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  df <- reshape2::melt(pm_block_matrix)
  colnames(df) <- c("Row", "Col", "Value")
  
  # Plot heatmap
  g <- ggplot(df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = "Matrix Heatmap", fill = "Value") + 
    ggtitle(g_title)
  
  return(g)
  
}

visualize_matrix_heatmap <- function(mat, g_title = NULL, zmin = NULL, zmid = NULL, zmax = NULL) {
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: plot a simple heatmap with optional arguments for max and min values
  #   - low = blue
  #   - mid = white
  #   - max = red
  #
  #
  # input:
  #
  # - mat       (p x p matrix)
  # - g_title   (string)
  # - zmin      (number) to denote the very smallest value
  # - zmid      (number) to denote the midpoint (white) value
  # - zmax      (number) to denote the highest value (red)
  #
  # output:
  #
  # - graph
  #
  # ----------------------------------------------------------------------------
  
  # Convert matrix to data frame for ggplot
  df <- reshape2::melt(mat)
  colnames(df) <- c("x", "y", "value")
  
  # Set defaults for color scale
  if (is.null(zmin)) zmin <- min(df$value, na.rm = TRUE)
  if (is.null(zmax)) zmax <- max(df$value, na.rm = TRUE)
  if (is.null(zmid)) zmid <- (zmin + zmax) / 2  
  
  if(zmin > zmid){
    zmin = zmid
  }
  
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue",     # negative values
      mid = "white",    # zero
      high = "red",     # positive values
      midpoint = zmid,     
      limits = c(zmin, zmax),
      oob = scales::squish
    ) +    
    # scale_fill_gradientn(
    #   colours = c("blue", "white", "red"),
    #   limits = c(zmin, zmax),
    #   oob = scales::squish
    # ) +
    coord_fixed() +
    theme_minimal() +
    scale_y_reverse() +  # So origin is at top-left like a matrix
    labs(x = NULL, y = NULL, fill = "Value") + 
    ggtitle(g_title)
  

}


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

visualize_error_histogram <- function(mat_est, mat_reconstruct, g_title, bin_count){
  
  # ----------------------------------------------------------------------------
  #
  # visualize the errors between mat_est and mat_reconstruct
  # 
  # mat_est = (matrix)
  # mat_reconstruct = (matrix)
  # g_title
  # n_bins = number of histogram bins
  #
  #
  # output:
  # 
  # histogram of elementwise differences
  #
  # ----------------------------------------------------------------------------
  
  errors <- as.numeric(mat_est - mat_reconstruct)
  
  df <- data.frame(errors = errors)
  
  g <- ggplot(df, aes(x = errors)) +
    geom_histogram(bins = bin_count, fill = "skyblue", color = "black") +
    labs(title = g_title) +
    theme_minimal()
  
  return(g)
}

visualize_log_intensity <- function(X_k, time_grid, g_title, legend_title = 'Process'){
  
  #
  # visualize the log intensities
  # 
  # X_k = (p x m matrix)
  # time_grid = (m-dim vec of timepoints)
  #
  #
  # output:
  # 
  # graph of all p log intensites at m timepoints
  
  p <- dim(X_k)[1]
  m <- dim(X_k)[2]
  
  rownames(X_k) <- paste0(seq_len(p))
  
  # Convert to long format
  df <- reshape2::melt(X_k)
  colnames(df) <- c("Process", "TimeIndex", "Value")
  
  df$Time <- time_grid[df$TimeIndex]
  df$Process <- factor(df$Process)
  g <- ggplot() +
    geom_line(data = df, aes(x = Time, y = Value, color = Process)) +
    theme_minimal() +
    labs(title = g_title, x = "Time", y = "Value", color = legend_title) +
    theme(legend.position = "right")  
  
  return(g)
  
}


ground_truth_rho_ij <- function(my_list){
  

  #
  # Assuming your list is called my_list
  # Each element has a matrix called X_functions of dimension p x m
  #   
  #
  # generate rho_ij estimate:
  #
  # 1) exponentiate log-intensity (X) to get intensity
  # 2) take cross product to get rho_ij
  #  
  #     mean( X %*% t(X) )
  #
  
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

# visualize how a precision matrix changes over time
# such as for banded_trig

visualize_precision_yc <- function(query_y_cs, p, adj_type, adj_params, ncores){
  
  #
  # GOAL: visualize how the partial correlation matrix changes over time with a gif
  #
  # input:
  #
  # - query_y_cs (n_query x q_c dim matrix)
  # - p 
  # - adj_type
  # - adj_params
  #
  #
  # output:
  #
  # - animated 
  #
  # 
  
  if(adj_type == 'banded_trig'){
    
    #
    # adj_params = [y_min, y_max, rho_max]
    #
    
    graphs <- pbmclapply(1:nrow(query_y_cs), function(k){
      
      y_c_k <- query_y_cs[k,]
      generate_sparse_precision_matrix(y_c_k, p, adj_type, adj_params)
    
    }, mc.cores = ncores)

    
    
  }
  
  # animate
  
  imgs <- list()
  
  for (i in seq_along(graphs)) {
    m <- graphs[[i]]$partial_cor_mat
    
    # Create a temporary image for each matrix heatmap
    tmpfile <- tempfile(fileext = ".png")
    png(tmpfile, width = 600, height = 600)
    
    visualize_matrix_heatmap(m, g_title = paste("Frame", i), -1, 1) %>% print()
    
    dev.off()
    
    imgs[[i]] <- image_read(tmpfile)
  }
  
  # parallel version - doesn't work
  #
  # imgs <- pbmclapply(1:nrow(query_y_cs), function(i){
  #   
  #   m <- graphs[[i]]$partial_cor_mat
  #   
  #   # Create a temporary image for each matrix heatmap
  #   tmpfile <- tempfile(fileext = ".png")
  #   png(tmpfile, width = 600, height = 600)
  #   
  #   visualize_matrix_heatmap(m, g_title = paste("Frame", i), -1, 1) %>% print()
  #   
  #   dev.off()
  #   
  #   image_read(tmpfile)    
  #   
  # }, mc.cores = ncores)   
  
  # Combine into an animated gif
  animation <- image_animate(image_join(imgs), fps = 100)
  image_write(animation, "heatmap_animation.gif")  
  
}