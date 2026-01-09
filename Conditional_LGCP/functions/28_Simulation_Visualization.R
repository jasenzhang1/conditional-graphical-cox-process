

#library(gganimate)
#library(magick)

source('functions/00a_matrix_massaging.R')
source('functions/24_log_intensity_generation.R')

visualize_matrix_heatmap <- function(mat, g_title = NULL, zmin = NULL, zmid = NULL, zmax = NULL, x_max_borders = NULL, y_max_borders = NULL, palette_ID = 'Blue-Red 2') {
  
  
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
  # - mat            (p x p matrix)
  # - g_title        (string)
  # - zmin           (number) to denote the very smallest value
  # - zmid           (number) to denote the midpoint (white) value
  # - zmax           (number) to denote the highest value (red)
  # - x_max_borders  (vector)  optional vector to denote block matrix borders
  # - y_max_borders  (vector)  optional vector to denote block matrix borders
  #
  # output:
  #
  # - graph
  #
  # ----------------------------------------------------------------------------
  
  # colors
  new_palette <- hcl.colors(3, palette = palette_ID)
  c_low <- new_palette[1]
  c_mid <- new_palette[2]
  c_high <- new_palette[3]
  
  # Convert matrix to data frame for ggplot
  df <- reshape2::melt(mat)
  colnames(df) <- c("y", "x", "value")
  
  # Set defaults for color scale
  if (is.null(zmin)) zmin <- min(df$value, na.rm = TRUE)
  if (is.null(zmax)) zmax <- max(df$value, na.rm = TRUE)
  if (is.null(zmid)) zmid <- (zmin + zmax) / 2  
  
  if(zmin > zmid){
    zmin = zmid
  }
  
  g <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      low = c_low, #"blue",     # negative values
      mid = c_mid, #"white",    # zero
      high = c_high, #"red",     # positive values
      midpoint = zmid,     
      limits = c(zmin, zmax),
      oob = scales::squish
    ) +    
    coord_fixed() +
    theme_minimal() +
    scale_y_reverse() +  # So origin is at top-left like a matrix
    labs(x = NULL, y = NULL, fill = "Value") + 
    ggtitle(g_title)
  
  # 3) if we have borders
  if (!is.null(x_max_borders)) {
    g <- g + geom_vline(xintercept = x_max_borders + 0.5, size = 0.5, alpha = 0.2)
  }
  
  if (!is.null(y_max_borders)) {
    g <- g + geom_hline(yintercept = y_max_borders + 0.5, size = 0.5, alpha = 0.2)
  }
  
  return(g)
  
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

visualize_histogram <- function(mat_est, mat_reconstruct = NULL, g_title = 'Title', bin_count = 20){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: 1) visualize the distribution of a matrix object
  #       2) if we supply two matrices, we wish to visualize the errors between mat_est and mat_reconstruct
  # 
  # 
  # inputs:
  # 
  # - mat_est           (matrix)
  # - mat_reconstruct   (matrix)
  # - g_title           (string)
  # - n_bins            (integer)  number of histogram bins
  #
  #
  # output:
  # 
  # histogram of elementwise differences
  #
  # ----------------------------------------------------------------------------
  
  if(is.null(mat_reconstruct)){
    plot_data = as.numeric(mat_est)
  } else{
    plot_data <- as.numeric(mat_est - mat_reconstruct)
  }
  
  df <- data.frame(plot_data = plot_data)
  
  g <- ggplot(df, aes(x = plot_data)) +
    geom_histogram(bins = bin_count, fill = "skyblue", color = "black") +
    labs(title = g_title) +
    theme_minimal()
  
  return(g)
}

visualize_log_intensity <- function(X_k, time_grid,  g_title = 'Title', mu_t = NULL, palette_ID = 'Dark 2', legend_title = 'Process', ymin = NULL, ymax = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # 
  # visualize the log intensities
  # 
  # inputs:
  # 
  # - X_k            (p x m matrix)
  # - time_grid      (m-dim vec of timepoints)
  # - g_title        (string)
  # - mu_t           (m-dim vec of the mean of the GP)
  # - legend_title   (string)
  #
  #
  # output:
  # 
  # graph of all p log intensites at m timepoints
  #
  # 
  # ----------------------------------------------------------------------------
  
  
  p <- dim(X_k)[1]
  m <- dim(X_k)[2]
  
  new_palette <- hcl.colors(p, palette = palette_ID)
  
  rownames(X_k) <- paste0(seq_len(p))
  
  # Convert to long format
  df <- reshape2::melt(X_k)
  colnames(df) <- c("Process", "TimeIndex", "Value")
  
  df$Time <- time_grid[df$TimeIndex]
  df$Process <- factor(df$Process)
  g <- ggplot() +
    geom_line(data = df, aes(x = Time, y = Value, color = Process)) +
    scale_color_manual(values = new_palette) + 
    theme_minimal() +
    labs(title = g_title, x = "Time", y = "Value", color = legend_title) +
    theme(legend.position = "right")  
  
  # add a dotted line to represent mean of the GP
  
  
  if(!is.null(mu_t)){
    mu_df <- data.frame(time_grid = time_grid, mu_t = mu_t)
    g <- g + geom_line(data = mu_df, aes(x = time_grid, y = mu_t), linetype = "dashed", size = 1, alpha = 0.7)
  }
  
  # ymin and ymax if offered
  if (!is.null(ymin) || !is.null(ymax)) {
    g <- g + coord_cartesian(ylim = c(ymin, ymax))
  }
  
  return(g)
  
}

visualize_intensity_with_points <- function(intensity, time_grid, event_times){
  
  # ----------------------------------------------------------------------------
  #
  # inputs:
  # 
  # - intensity    (m-dim vector)   y-values of the intensity
  # - time_grid    (m-dim vector)   x-values of the intensity
  # - event_times  (n-dim vector)   vector of timestamps of the events
  #
  #
  # outputs:
  #
  # - g            (ggplot object)  graph overlaying the intensity with step function of events, along with dots on the y = 0 line to represent realizations
  #
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # Prepare data frames
  # ----------------------------------------------------------------------------
  df_intensity <- data.frame(
    time = time_grid,
    intensity = intensity
  )
  
  df_step <- data.frame(
    time = sort(event_times),
    count = seq_along(event_times)
  )
  
  # Scale factor to map cumulative count to intensity range
  scale_factor <- max(intensity) / max(df_step$count)
  df_step$scaled_count <- df_step$count * scale_factor
  
  # ----------------------------------------------------------------------------
  # Plot
  # ----------------------------------------------------------------------------
  g <- ggplot() +
    # intensity curve on left y-axis
    geom_line(data = df_intensity, aes(x = time, y = intensity), color = "blue", size = 1) +
    
    # step function scaled to match left y-axis (will show right axis)
    geom_step(data = df_step, aes(x = time, y = scaled_count), color = "darkgreen", linetype = "dashed") +
    
    # event points at y=0
    geom_point(data = data.frame(time = event_times, y = rep(0, length(event_times))),
               aes(x = time, y = y), color = "red", size = 2) +
    
    # labels
    scale_y_continuous(
      name = "Intensity",
      sec.axis = sec_axis(~ . / scale_factor, name = "Cumulative Events")
    ) +
    
    labs(x = "Time", title = "Intensity Function with Event Points and Step Function") +
    theme_minimal()
  
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

# visualize how a precision matrix changes over time using a gif
# such as for banded_trig

visualize_precision_gif <- function(p, adj_type, adj_params, ncores){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: visualize how the partial correlation matrix changes over time with a gif
  #
  # - note that adj_params[1:2] denote the min and max y_c value
  # 
  #
  # input:
  #
  # - p             (integer)
  # - adj_type      (string)
  # - adj_params    (vector)
  #
  #
  # output:
  #
  # - animated gif
  #
  # ----------------------------------------------------------------------------
  
  
  # 1) generate precision and correlation (pxp) matrices
  
  y_min <- adj_params[1]
  y_max <- adj_params[2]
  n_times <- 100
  query_y_cs <- matrix(seq(y_min, y_max, length.out = n_times))
  

    
  graphs <- lapply(1:nrow(query_y_cs), function(k){
    
    y_c_k <- query_y_cs[k,]
    prec_mats <- generate_sparse_precision_matrix(y_c_k, p, adj_type, adj_params)
    
    list(adj_mat = prec_mats$adj_mat,
         prec_mat = prec_mats$prec_mat,
         cor_mat = prec_mats$cor_mat)
  
  })


  # 2) animate
  
  imgs <- list()
  
  for (i in seq_along(graphs)) {
    
    m1 <- graphs[[i]]$adj_mat
    m2 <- graphs[[i]]$prec_mat
    m3 <- graphs[[i]]$cor_mat
    
    y_c_print <- format(round(query_y_cs[i,], 2), nsmall = 2)
    
    # Create a temporary image for each matrix heatmap
    tmp1 <- tempfile(fileext = ".png")
    png(tmp1, width = 400, height = 400)
    gif_title1 <- paste0("Adj Matrix: y_c = ", y_c_print) 
    visualize_matrix_heatmap(m1, g_title = gif_title1, zmin = -1, zmid = 0, zmax = 1) %>% print()
    dev.off()
    
    # Middle plot
    
    tmp2 <- tempfile(fileext = ".png")
    png(tmp2, width = 400, height = 400)
    gif_title2 <- paste0("Prec Matrix: y_c = ", y_c_print) 
    prec_max <- max(sapply(graphs, function(x){max(as.numeric(x$prec_mat))}))
    prec_min <- min(sapply(graphs, function(x){min(as.numeric(x$prec_mat))}))
    visualize_matrix_heatmap(m2, g_title = gif_title2, zmin = prec_min, zmid = 0, zmax = prec_max) %>% print()
    dev.off()
    
    # Last plot
    
    tmp3 <- tempfile(fileext = ".png")
    png(tmp3, width = 400, height = 400)
    gif_title3 <- paste0("Cor Matrix: y_c = ", y_c_print) 
    visualize_matrix_heatmap(m3, g_title = gif_title3, zmin = -1, zmid = 0, zmax = 1) %>% print()
    dev.off()    
    
    # ---- read images and combine side by side ----
    img1 <- image_read(tmp1)
    img2 <- image_read(tmp2)
    img3 <- image_read(tmp3)
    
    combined <- image_append(c(img1, img2, img3))  # horizontal side-by-side
    imgs[[i]] <- combined
  }
  

  
  # Combine into an animated gif
  animation <- image_animate(image_join(imgs), fps = 100)
  image_write(animation, "heatmap_animation_3side_single_v2.gif")  
  
}

