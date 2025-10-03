
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
  # GOAL: visualize how the partial correlation matrix changes over time
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


# visualize ||V_cond - V_cond_est||_HS convergence

visualize_V_cond_convergence <- function(folder_name, mat_name, i, j){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: visualize V_cond_convergence at block matrix (i, j)
  #
  #       error =  || V_{X_i, X_j}^(y_c) - \hat{V}_{X_i, X_j}^(y_c) ||_HS
  #
  #
  # input:
  #
  # - folder_name (string)  'simu_results_banded_c1_3'
  # - mat_name    (string)  'V', 'C', or 'P'
  # - i, j        (scalar)  block matrix numbers
  #
  #
  # output:
  #
  # - graph of V_12 error as we change sample size 
  #
  #
  # ----------------------------------------------------------------------------
  

  files <- list.files(folder_name, full.names = TRUE)
  
  ns <-  as.numeric(sub(".*_(.*)\\.RData$", "\\1", files)) # retrieve numbers
  
  results_list <- lapply(files, function(f) {
    e <- new.env()          # create an isolated environment
    load(f, envir = e)      # load into that environment
    as.list(e)              # convert to a list (in case multiple objects)
  })
  
  names(results_list) <- ns
  
  
  # retrieve the 1_2 entry from metrics of V_cond 9written by chatgpt)
  
  results_df <- do.call(rbind, lapply(names(results_list), function(top_name) {
    
    top_entry <- results_list[[top_name]]$graph_results_i$estimated_graphs_part_2
    
    do.call(rbind, lapply(names(top_entry), function(low_name) {
      low_entry <- top_entry[[low_name]]
      
      # Dynamically select the matrix
      mat <- if (mat_name == "V") {
        low_entry$metrics$V_HS
      } else if (mat_name == "C") {
        low_entry$metrics$C_HS
      } else if (mat_name == "P") {
        low_entry$metrics$P_HS
      } else {
        stop("Unknown mat_name: must be 'V', 'C', or 'P'")
      }
      
      value <- mat[i, j]
      
      data.frame(
        top_level = top_name,
        low_level = low_name,
        value = value,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  colnames(results_df) <- c('n', 'y_c_query', 'V_12_error')
  
  results_df$n <- as.numeric(results_df$n)
  results_df$y_c_query <- as.numeric(results_df$y_c_query)
  
  if(mat_name == 'V'){
    title_name <- 'Covariance'
  } else if(mat_name == 'C'){
    title_name <- 'Correlation'
  } else if (mat_name == 'P'){
    title_name <- 'Precision'
  } else{
    stop("Unknown mat_name: must be 'V', 'C', or 'P'")
  }
  
  g <- ggplot() + geom_line(data = results_df, aes(x = y_c_query, y = V_12_error, group = n, color = n)) + 
    ylab(paste0(mat_name, '_' , i, '_', j, ' Error')) + 
    xlab('Y_c Query') + 
    ggtitle(paste0('Conditional ', title_name, ' Operator Convergence')) + 
    ylim(0, NA) + 
    theme_bw() 
  
  return(g)
  
    
}

# across all n's, plot Y_c_query vs AUC

visualize_AUC_across_n <- function(folder_name){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: visualize AUC metric across n and y_c_query
  #
  #
  #
  # input:
  #
  # - folder_name (string)  'simu_results_banded_c1_3'
  #
  #
  # output:
  #
  # - graph of AUC vs y_c_query (x-axis) and n (color)
  #
  #
  # ----------------------------------------------------------------------------
  
  
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <-  as.numeric(sub(".*_(.*)\\.RData$", "\\1", files)) # retrieve numbers
  
  results_list <- lapply(files, function(f) {
    e <- new.env()          # create an isolated environment
    load(f, envir = e)      # load into that environment
    as.list(e)              # convert to a list (in case multiple objects)
  })
  
  names(results_list) <- ns
  
  
  # retrieve the 1_2 entry from metrics of V_cond 9written by chatgpt)
  
  results_df <- do.call(rbind, lapply(names(results_list), function(top_name) {
    
    top_entry <- results_list[[top_name]]$graph_results_i$estimated_graphs_part_2
    
    do.call(rbind, lapply(names(top_entry), function(low_name) {
      low_entry <- top_entry[[low_name]]
      

      
      AUC <- low_entry$metrics$auc
      
      data.frame(
        top_level = top_name,
        low_level = low_name,
        value = AUC,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  colnames(results_df) <- c('n', 'y_c_query', 'AUC')
  
  results_df$n <- as.numeric(results_df$n)
  results_df$y_c_query <- as.numeric(results_df$y_c_query)
  
  title_name <- 'AUC'
  

  
  g <- ggplot() + geom_line(data = results_df, aes(x = y_c_query, y = AUC, group = n, color = n)) + 
    ylab('AUC') + 
    xlab('Y_c Query') + 
    ggtitle('AUC versus n and y_c_query') + 
    ylim(0, 1) + 
    theme_bw() 
  
  return(g)
  
  
}


visualize_metrics <- function(folder_name, metrics, i = NULL, j = NULL){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: visualize convergences of various metrics:
  #
  # - metrics
  #   - rho_i_dist (scalar)
  #   - rho_ij_dist (pxp matrix, each value is HS norm of m_est x m_est rho_ij)
  #   - g_ij_dist   (pxp matrix)
  #   - P_HS        (pxp matrix, each value is HS norm of difference of P_hat - P)
  #   - C_HS        (pxp matrix)
  #   - V_HS        (pxp matrix)
  #   - sens        (scalar)
  #   - spec        (scalar)
  #   - auc         (scalar)
  #   - accuracy    (scalar)
  #
  #
  #
  # input:
  #
  # - folder_name (string)  'simu_results_banded_c1_3'
  # - metric      (string)  
  #   - 'rho_i_dist'
  #   - 'rho_ij_dist'
  #   - 'g_ij_dist'
  #   - 'P_HS', 'C_HS', 'V_HS'
  # - i and j    (integers)  indices for matrix metrics
  #
  # output:
  #
  #
  #
  # ----------------------------------------------------------------------------
  
  
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <-  as.numeric(sub(".*_(.*)\\.RData$", "\\1", files)) # retrieve numbers
  
  results_list <- lapply(files, function(f) {
    e <- new.env()          # create an isolated environment
    load(f, envir = e)      # load into that environment
    as.list(e)              # convert to a list (in case multiple objects)
  })
  
  names(results_list) <- ns
  
  
  # retrieve the 1_2 entry from metrics of V_cond (written by chatgpt)
  
  # top_name = '100', '200' etc n
  # low_name = '0', '0.125', '0.25', etc y_c_query  
  results_df <- do.call(rbind, lapply(names(results_list), function(top_name) {
    
    top_entry <- results_list[[top_name]]$graph_results_i$estimated_graphs_part_2
    
    do.call(rbind, lapply(names(top_entry), function(low_name) {
      low_entry <- top_entry[[low_name]]
      
      values <- c()
      for(metric in metrics){
      
        metric_value <- low_entry[['metrics']][[metric]]
        
        if(metric %in% c('rho_ij_dist', 'g_ij_dist', 'P_HS', 'C_HS', 'V_HS')){
          metric_value = metric_value[i, j]
        }
        
        values <- c(values, metric_value)
        
      }
    
      
      c(top_name, low_name, values)
    }))
  }))
  
  # names 
  metric_names <- c()
  for(metric in metrics){
    
  
    if(metric %in% c('rho_ij_dist', 'g_ij_dist', 'P_HS', 'C_HS', 'V_HS')){
      metric_name <- paste0(metric, '_', i, '_', j)
    } else{
      metric_name <- metric  
    }
    
    metric_names <- c(metric_names, metric_name)
  }
  
  results_df <- data.frame(results_df) %>% mutate_all(as.numeric)
  colnames(results_df) <- c('n', 'y_c_query', metric_names)
  
  
  results_df$n <- as.factor(results_df$n)
  
  # graph
  graphs <- list()
  
  for(metric_name in metric_names){
    
    title_name <- paste0(metric_name, ' versus n and y_c_query')
    
    g <- ggplot() + geom_line(data = results_df, aes(x = y_c_query, y = .data[[metric_name]], group = n, color = n)) + 
      ylab(metric_name) + 
      xlab('Y_c Query') + 
      # ggtitle(title_name) + 
      theme_bw() 
    
    if(metric_name %in% c('auc')){
      g <- g + ylim(0, 1)
    } else{
      g <- g + ylim(0, NA)
    }
    graphs[[metric_name]] <- g
  }
  
  # grid arrange
  
  arranged_plots <- do.call(arrangeGrob, c(graphs, ncol = 3))
  
  # Display it
  return(grid.arrange(arranged_plots))
  
    
}
