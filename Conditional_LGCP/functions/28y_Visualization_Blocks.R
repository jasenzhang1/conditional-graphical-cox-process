source('functions/00c_block_matrix_arrange.R')
source('functions/28z_Visualization_helpers.R')
source('functions/28_Simulation_Visualization.R')
source('functions/00b_matrix_norms.R')


# functions to compartmentalize all the g_XX graphs

# any function that ends in "prep" returns a list of plots ready to be printed
# all other functions store the canvas in g and return(g)

result_01_prep <- function(step_0){
  
  # return a list of plots
  
  # Build plots dynamically
  plots <- lapply(seq_len(length(step_0)), function(i) {
    visualize_matrix_heatmap(step_0[[i]],
                             paste0('y_c = ', as.character(i)), -10, NULL, 10)
  })  
  
  return(plots)
}

result_01 <- function(step_0){
  

  plots <- result_01_prep(step_0)
  
  # Add caption grob
  caption <- textGrob("0. Ground Truth\n Precision wrt Queried\nContinuous Covariate", gp = gpar(fontsize = 14))
  
  # Combine plots + caption
  grobs <- c(plots, list(caption))
  
  g <- grid.arrange(grobs = grobs, nrow = 1)  
  return(g)
}

result_11_prep <- function(step_1, time_grid, time_grid_est){
  #   - X_k_est                   (p x m_est   x n)
  #   - X_k_truth                 (p x m_truth x n)
  #   - X_k_coarse_truth          (p x m_est   x n)
  #   - X_k_both_truth            (p x m_both  x n)
  #   - Lambda_k_truth            (p x m_truth x n)
  #   - Lambda_k_coarse_truth     (p x m_est   x n)
  g_list <- list(visualize_log_intensity(step_1[[1]][1:5,,1],   time_grid_est,  'Estimate'),
                 visualize_log_intensity(step_1[[3]][1:5,,1],   time_grid_est,  'Coarser Truth'),
                 visualize_log_intensity(step_1[[2]][1:5,,1],   time_grid,      'Finer Truth'),
                 visualize_log_intensity(step_1[[4]][1:5,,1],   time_grid_both, 'Combined Truth'))
  
  return(g_list)
}

result_22_prep <- function(step_2, time_grid, time_grid_est, ymin = NULL, ymax = NULL, full = T){
  
  if(! full){
    g_list <- list(visualize_log_intensity(step_2$rho_i_est[1:5,],     time_grid_est, 'Estimate', ymin = ymin, ymax = ymax),
                   visualize_log_intensity(step_2$rho_i_truth[1:5,],   time_grid,     'Truth', ymin = ymin, ymax = ymax))
    return(g_list)
  }
  
  #   - rho_i_truth               (p x m)
  #   - rho_i_coarse_truth        (p x m_est)
  #   - rho_i_X_truth             (p x m)
  #   - rho_i_X_coarse_truth      (p x m_est)
  #   - rho_i_est                 (p x m_est)
  #   - rho_list                  (list format, [[1]] = rho_i, [[2]] = rho_ij, rho_ij may be rho_ii only)  
  
  g_list <- list(visualize_log_intensity(step_2[[1]][1:5,],   time_grid,     'Truth Theory'),
                 visualize_log_intensity(step_2[[2]][1:5,],   time_grid_est, 'Coarse Truth Theory'),
                 visualize_log_intensity(step_2[[3]][1:5,],   time_grid,     'Truth X'),
                 visualize_log_intensity(step_2[[4]][1:5,],   time_grid_est, 'Coarse Truth X'),
                 visualize_log_intensity(step_2[[5]][1:5,],   time_grid_est, 'Estimate'))
  
  return(g_list)
}

result_20s_prep <- function(step_2b, i, j, zmin = NULL, zmax = NULL, full = T){
  #   - rho_ii_truth               (list of m x m matrices)
  #   - rho_ii_coarse_truth        (list of m_est x m_est matrices)
  #   - rho_ii_X_truth             (list of m x m matrices)
  #   - rho_ii_X_coarse_truth      (list of m_est x m_est matrices)
  #   - rho_ii_est                 (list of m_est x m_est matrices)
  
  key <- paste0(i, '_', j)
  
  if(! full){
    g_list <- list(visualize_matrix_heatmap(step_2b$rho_ii_truth[[key]],   'Truth Theory',        zmin, NULL, zmax),
                   visualize_matrix_heatmap(step_2b$rho_ii_est[[key]],     'Estimate',            zmin, NULL, zmax))
    return(g_list)
  }
  
  g_list <- list(visualize_matrix_heatmap(step_2b[[1]][[key]],   'Truth Theory',        zmin, NULL, zmax),
                 visualize_matrix_heatmap(step_2b[[2]][[key]],   'Coarse Truth Theory', zmin, NULL, zmax),
                 visualize_matrix_heatmap(step_2b[[3]][[key]],   'Truth X',             zmin, NULL, zmax),
                 visualize_matrix_heatmap(step_2b[[4]][[key]],   'Coarse Truth X',      zmin, NULL, zmax),
                 visualize_matrix_heatmap(step_2b[[5]][[key]],   'Estimate',            zmin, NULL, zmax))
  
  return(g_list)
}

result_heatmap_ij_prep <- function(my_list, entry_name, time_grid_est, i, j, palette_ID = 'Blue-Red 2', zmin = NULL, zmid = NULL, zmax = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: using facet_grid to plot all heatmaps together with one legend, one x-axis, one y-axis, etc...
  #
  #
  # inputs:
  #
  # - my_list     (list)     step_x; list of y_c_queries --> g_ij_est etc items
  # - entry_name  (string)   estimate prefix (e.g. rho_ii, g_ij)
  # - time_grid_est
  # - i 
  # - j
  #
  #
  # outputs:
  #
  # - ggplot object
  # 
  # ----------------------------------------------------------------------------
  

  key <- paste0(i, '_', j)
  
  # 1) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(my_list[[1]]), entry_name)
  
  queried_names <- names(my_list[[1]])
  y_c_names <- names(my_list)

  # 2) colors
  new_palette <- hcl.colors(3, palette = palette_ID)
  c_low <- new_palette[1]
  c_mid <- new_palette[2]
  c_high <- new_palette[3]  
  
  # Suppose your list is called `my_list` with length m
  m <- length(my_list)
  
  # Combine all matrices into one long dataframe
  df_all <- bind_rows(lapply(seq_len(m), function(k) {
    
    entry <- my_list[[k]]
    dfs <- list()
    
    for (t in seq_along(queried_names)) {
      
      this_name <- queried_names[t]   # e.g. rho_est, rho_X_truth, ...
      this_type <- suffix_names[t]    # e.g. est, X_truth, ...
      
      # check existence
      if (is.null(entry[[this_name]])) next
      
      # extract slice assuming it's in list form
      M <- entry[[this_name]][[key]]

      
      # melt into long format
      df_t <- reshape2::melt(M)
      
      df_t$Type <- this_type
      df_t$y_c <- y_c_names[k]
      
      dfs[[length(dfs) + 1]] <- df_t
    }
    
    bind_rows(dfs)
    
  }), .id = NULL)
  
  colnames(df_all)[1:3] <- c("Row", "Col", "Value")
  
  # Convert to factors for proper ordering
  df_all$Type <- factor(df_all$Type)
  df_all$y_c <- factor(df_all$y_c, levels = y_c_names)
  
  if(is.null(zmin)){
    zmin <- min(df_all$Value)
  }
  if(is.null(zmax)){
    zmax <- max(df_all$Value)
  }
  if(is.null(zmid)){
    zmid <- (zmin + zmax) / 2       
  }
  
  zmax <- zmax + 0.1 * (zmax - zmin)
  zmin <- zmin - 0.1 * (zmax - zmin)
  

  df_all$Col_val <- time_grid_est[df_all$Col]
  df_all$Row_val <- time_grid_est[df_all$Row]
  
  # Plot with facets
  ggplot(df_all, aes(x = Col_val, y = Row_val, fill = Value)) +
    geom_tile() +
    scale_y_reverse() + # matrix y-axis 
    scale_fill_gradient2(low = c_low, mid = c_mid, high = c_high,
                         midpoint = zmid,
                         limits = c(zmin, zmax)) +
    coord_fixed() +
    theme_minimal() +
    facet_grid(Type ~ y_c, scales = "fixed") +
    labs(x = "t", y = "s", fill = "f(s,t)") +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90),
      axis.text.y = element_text()
    )
}

result_heatmap_nonblock_prep <- function(my_list, entry_name, data_format, time_grid_est = NULL, rm_diag = F, palette_ID = 'Blue-Red 2', zmin = NULL, zmid = NULL, zmax = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: using facet_grid to plot all heatmaps together with one legend, one x-axis, one y-axis, etc...
  #
  #       instead of extracting the ij-th block, we look at the entire block
  #
  # inputs:
  #
  # - my_list          (list)     step_x; list of y_c_queries --> g_ij_est etc items
  # - entry_name       (string)   estimate prefix (e.g. rho_ii, g_ij)
  # - data_format      (string)   'full', 'regular', 'list'
  # - time_grid_est    (vector)   do we substitute indices with timestamps? If so, provide it
  # - rm_diag          (boolean)  do we remove the diag term?
  # - palette_ID       (string)
  # - zmin, zmid, zmax   (values)   do we manually decide on the bordering color values?
  #
  # outputs:
  # 
  # - ggplot object
  # 
  # ----------------------------------------------------------------------------
  
  # 1) obtain `queried_names` and `value_names`
  
  suffix_names <- step_00_grab_ID(names(my_list[[1]]), entry_name)
  
  queried_names <- names(my_list[[1]])
  
  sub_names <- names(my_list)
  
  # 2) colors
  new_palette <- hcl.colors(3, palette = palette_ID)
  c_low <- new_palette[1]
  c_mid <- new_palette[2]
  c_high <- new_palette[3]  
  
  # 3) combining entrys from `my_list`
  
  # Suppose your list is called `my_list` with length m
  m <- length(my_list)
  
  # Combine all matrices into one long dataframe
  df_all <- bind_rows(lapply(seq_len(m), function(k) {
    
    entry <- my_list[[k]]
    dfs <- list()
    
    for (t in seq_along(queried_names)) {
      
      this_name <- queried_names[t]     # e.g. rho_est, rho_X_truth, rho_truth
      this_type <- suffix_names[t]      # e.g. est, X_truth, truth
      
      # If an expected matrix isn't present, skip gracefully
      if (is.null(entry[[this_name]])) next
      
      # Pull matrix depending on the format
      if (data_format %in% c("full", "regular")) {
        M <- entry[[this_name]]
      } else {
        # infer dimension p from the i_j-block names
        i_j <- names(entry[[this_name]])
        j_j <- i_j[length(i_j)]
        p <- as.numeric(strsplit(j_j, "_")[[1]][1])
        M <- assemble_block_matrix_irregular(entry[[this_name]], p)$block_matrix
      }
      
      df_t <- reshape2::melt(M)
      df_t$Type <- this_type
      df_t$y_c <- sub_names[k]
      dfs[[length(dfs) + 1]] <- df_t
    }
    
    bind_rows(dfs)
    
  }), .id = NULL)
  
  colnames(df_all)[1:3] <- c("Row", "Col", "Value")
  
  # Convert to factors for proper ordering
  df_all$Type <- factor(df_all$Type)
  df_all$y_c <- factor(df_all$y_c)
  
  # now, df_all has 5 columns:
  #
  # - Row   (integer)
  # - Col   (integer)
  # - Value (number)
  # - Type  (factor) estimand such as est, truth, X_truth
  # - y_c   (factor) continuous covariate value
  
  
  # 4) Remove diagonal or diagonal blocks

  if (rm_diag) {
    
    if (data_format == "regular") {
      # Remove simple diagonal entries
      df_all <- df_all[df_all$Row != df_all$Col, ]
      
    } else if (data_format %in% c("list", "full")) {
      
      # Infer p: number of time grid points
      p <- length(time_grid_est)
      
      # Determine block index for each Row/Col
      df_all$Row_block <- ceiling(df_all$Row / p)
      df_all$Col_block <- ceiling(df_all$Col / p)
      
      # Remove diagonal blocks: block (i, i)
      df_all <- df_all[df_all$Row_block != df_all$Col_block, ]
      
      # Clean up helper columns
      df_all$Row_block <- NULL
      df_all$Col_block <- NULL
    }
  }
  
  # 5) min and max of legend 
  if(is.null(zmin)){
    zmin <- min(df_all$Value)
  }
  if(is.null(zmax)){
    zmax <- max(df_all$Value)
  }
  if(is.null(zmid)){
    zmid <- (zmin + zmax) / 2       
  }
  
  zmax <- zmax + 0.1 * (zmax - zmin)
  zmin <- zmin - 0.1 * (zmax - zmin)
  
  # 6) convert indices to timestamps if relevant and dimensions add up
  
  if(!is.null(time_grid_est) & length(unique(df_all$Col)) == length(time_grid_est)){
    df_all$Col_val <- time_grid_est[df_all$Col]
    df_all$Row_val <- time_grid_est[df_all$Row]    
  } else{
    df_all$Col_val <- df_all$Col
    df_all$Row_val <- df_all$Row    
  }

  
  # 7) Plot with facets
  g <- ggplot(df_all, aes(x = Col_val, y = Row_val, fill = Value)) +
    geom_tile() +
    scale_y_reverse() + # matrix y-axis 
    scale_fill_gradient2(low = c_low, mid = c_mid, high = c_high,
                         midpoint = zmid,
                         limits = c(zmin, zmax)) +
    coord_fixed() +
    theme_minimal() +
    facet_grid(Type ~ y_c, scales = "fixed") +
    labs(x = "t", y = "s", fill = "f(s,t)") +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90),
      axis.text.y = element_text()
    )
  
  return(g)
}

result_line_graph_prep <- function(my_list, entry_name, time_grid, grouping = 'estimand', palette_ID = 'Dark 2', num_processes = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: using facet_grid to plot all line graphs together with one legend, one x-axis, one y-axis, etc...
  #
  #
  # inputs:
  #
  # - my_list          (list)     step_x; list of y_c_queries --> g_ij_est etc items
  # - entry_name       (string)   estimate prefix (e.g. rho_i)
  # - time_grid        (vector)   vector of timepoints
  # - grouping         (string)   'estimand' or 'process' if we group by estimand, different processes will be placed together.
  #                                                       if we group by process, different estimands will be placed together.
  # - palette_ID       (string)   preset colors
  # - num_processes    (integer)  do we want to trim the number of processes?
  #
  #
  # ouptput:
  #
  # - ggplot object
  #
  # 
  # ----------------------------------------------------------------------------  
  
  # 0) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(my_list[[1]]), entry_name)

  queried_names <- names(my_list[[1]])
  
  
  # Suppose your list is called `my_list` with length m
  m <- length(my_list)
  y_c_names <- names(my_list)
  
  # 1) Combine all matrices into one long dataframe
  df_all <- bind_rows(lapply(seq_len(m), function(i) {
    
    entry <- my_list[[i]]
    dfs <- list()
    
    for (j in seq_along(queried_names)) {
      
      this_name <- queried_names[j]         # e.g. rho_i_truth, rho_i_est ...
      this_type <- suffix_names[j]          # e.g. truth, est ...
      
      # Skip if not found in list
      if (is.null(entry[[this_name]])) next
      
      A <- entry[[this_name]]
      
      # Optionally trim by number of processes
      if (!is.null(num_processes)) {
        A <- A[1:num_processes, ]
      }
      
      df_ij <- reshape2::melt(A)
      df_ij$Var2 <- time_grid[df_ij$Var2]  # we assume truth and est are all the same 
      
      df_ij$Type <- this_type
      df_ij$y_c <- y_c_names[i]
      
      dfs[[length(dfs) + 1]] <- df_ij
    }
    
    bind_rows(dfs)
  }), .id = NULL)
  
  # 2) renaming + convert to factors for proper ordering
  colnames(df_all)[1:3] <- c("process", "time", "Value")

  df_all$Type <- factor(df_all$Type)
  df_all$process <- as.factor(df_all$process)
  df_all$y_c     <- as.factor(df_all$y_c)
  
  # now, we have df_all with the following columns:
  #
  # - process  (factor)  1 through p
  # - time     (value)    which timestamp in time_grid_est 
  # - value    (number)   y-axis value
  # - type     (factor)   suffix such as truth, est, X_truth etc
  # - y_c      (factor)   y_c value 
  
  # 3) plot with facets depending on grouping
  
  if(grouping == 'estimand'){
    g <- ggplot(df_all, aes(x = time, y = Value, color = process, group = process)) +
      geom_line(alpha = 0.5, size = 1) +
      facet_grid(
        rows = vars(Type),
        cols = vars(y_c)
      ) +
      scale_color_brewer(palette = "Dark2") +
      theme_bw() +
      labs(
        x = "Time",
        y = "Intensity",
        color = "Process"
      ) +
      theme(
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        strip.placement = "outside"  # optional, puts labels outside the panel
      )
  } else if(grouping == 'process'){
    g <- ggplot(df_all, aes(x = time, y = Value, color = Type, group = Type)) +
      geom_line(alpha = 0.5, size = 1) +
      facet_grid(
        rows = vars(process),
        cols = vars(y_c)
      ) +
      scale_color_brewer(palette = "Dark2") +
      theme_bw() +
      labs(
        x = "Time",
        y = "Intensity",
        color = "Process"
      ) +
      theme(
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        strip.placement = "outside"  # optional, puts labels outside the panel
      )
  } else{
    g <- NULL
  }
  
  return(g)

}

result_histogram_prep <- function(my_list, entry_name, data_format, nbins = 20){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: using facet_grid to plot all histograms together with one legend, one x-axis, one y-axis, etc...
  #
  #
  # inputs:
  #
  # - my_list     (list)     step_x; list of y_c_queries --> g_ij_est etc items
  # - entry_name  (string)   estimate prefix (e.g. rho_i)
  # - data_format (string)   describes how the data is packaged (rho_i_est = matrix), options include 'full', 'regular', 'list', 'vector'
  # - nbins       (integer)  number of bins
  #
  #
  # outputs:
  #
  # - ggplot output
  # 
  # ----------------------------------------------------------------------------  
  
  # 1) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(my_list[[1]]), entry_name)
  
  queried_names <- names(my_list[[1]])
  y_c_names <- names(my_list)
  
  # 2) collect items
  
  # Suppose your list is called `my_list` with length m
  n_query <- length(my_list)
  
  
  # Combine all matrices into one long dataframe
  df_all <- bind_rows(lapply(seq_len(n_query), function(k) {
    
    entry <- my_list[[k]]
    dfs <- list()
    
    for (t in seq_along(queried_names)) {
      
      this_name <- queried_names[t]     # e.g. rho_est, rho_X_truth, rho_truth
      this_type <- suffix_names[t]      # e.g. est, X_truth, truth
      
      # If an expected matrix isn't present, skip gracefully
      if (is.null(entry[[this_name]])) next
      
      # Pull matrix depending on the format
      if (data_format %in% c("full", "regular")) {
        M <- entry[[this_name]]
      } else {
        # infer dimension p from the i_j-block names
        i_j <- names(entry[[this_name]])
        j_j <- i_j[length(i_j)]
        p <- as.numeric(strsplit(j_j, "_")[[1]][1])
        M <- assemble_block_matrix_irregular(entry[[this_name]], p)$block_matrix
      }
      
      df_t <- reshape2::melt(M)
      df_t$Type <- this_type
      df_t$y_c <- y_c_names[k]
      dfs[[length(dfs) + 1]] <- df_t
    }
    
    bind_rows(dfs)
    
  }), .id = NULL)
  
  
  # 3) padding and graphing
  
  df_all$Type <- factor(df_all$Type)   
  df_all$y_c <- factor(df_all$y_c)
  
  # Histogram plot
  ggplot(df_all, aes(x = value, fill = Type)) +
    geom_histogram(bins = nbins, fill = "skyblue", color = "black", position = "identity") +
    facet_grid(Type ~ y_c, scales = "free_y") +
    theme_minimal() +
    labs(
      x = "Value",
      y = "Count",
      fill = "Type"
    ) +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold")
    )
  
}

result_29 <- function(weights, y_c_values, y_c_id){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: plot the distribution of the weights in estimating rho_i(t)
  #
  # - weights    (n-dim vector)  weight vector
  # - y_c_values (n-dim vector)  time vector
  # - y_c_id     (value)         y_c_query
  #
  # ----------------------------------------------------------------------------

  
  df_29 <- data.frame(x = y_c_values, y = weights)
  g <- ggplot(data = df_29, aes(x = x, y = y)) + geom_line() + geom_point() + 
    ylab('weight') + 
    xlab('continuous covariate') + 
    ggtitle(paste0('y_c_query=', y_c_id))
  
  return(g)
  
}

result_30s_prep <- function(step_3, i, j, full = T){
  #   - g_ij_ground_truth                (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_ground_truth         (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_truth                       (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_truth                (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_X_truth                     (list of m x m matrices for i_j entries)
  #   - g_ij_X_coarse_truth              (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_est                         (list of m_est x m_est matrices for i_j entries)   
  
  key <- paste0(i, '_', j)
  
  if(! full){
    g_list <- list(visualize_matrix_heatmap(step_3$g_ij_truth[[key]], 'Truth_v1',        -1, 0, 1),
                   visualize_matrix_heatmap(step_3$g_ij_truth_v2[[key]], 'Truth_v2',     -1, 0, 1),
                   visualize_matrix_heatmap(step_3$g_ij_est[[key]], 'Estimate',          -1, 0, 1))
    return(g_list)
  }
  
  g_list <- list(visualize_matrix_heatmap(step_3[[1]][[key]], 'Ground Truth',        -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[2]][[key]], 'Coarse Ground Truth', -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[3]][[key]], 'Truth Theory',        -1, 0, 1),  
                 visualize_matrix_heatmap(step_3[[4]][[key]], 'Coarse Truth Theory', -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[5]][[key]], 'Truth X',             -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[6]][[key]], 'Coarse Truth X',      -1, 0, 1),
                 visualize_matrix_heatmap(step_3[[7]][[key]], 'Estimate',            -1, 0, 1))
  
  return(g_list)
}

result_41_prep <- function(step_4){
  
  # list --> n_query --> eigen_decomp_suffix
  
  full_names <- names(step_4[[1]])
  
  step_4_v2 <- lapply(step_4, function(x) {
    lapply(full_names, function(nm) {
      x[[nm]]$eigenfunctions[[1]] %>% t()
    }) %>% setNames(full_names)
  })
  
  
  return(step_4_v2)
  
}

# does eigenreconstruction give us our original g_ij?
result_42_prep <- function(step_3, step_4, p){
  
  # assume step_3 and step_4 are n_query specific lists
  
  source('functions/04_eigendecomposition.R')
  source('functions/13_estimation_validation.R')
  
  g_list <- list()
  
  # 0) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(step_3), 'g_ij')

  
  queried_names_3 <- names(step_3)
  queried_names_4 <- names(step_4)
  
  # 2) for each suffix name, do it 
  for(i in 1:length(suffix_names)){
    
    
    g_ii_suffix <- prep_eigendecomposition_ii(step_3[[queried_names_3[i]]], p)
    g_ii_suffix_decomp <- validate_eigendecomposition_ii(g_ii_suffix, step_4[[queried_names_4[i]]])
    
    name_i <- suffix_names[i]
    name_i_reconstruct <- paste0(name_i, ' reconstruct')
    
    g_list[[length(g_list) + 1]] <- visualize_matrix_heatmap(g_ii_suffix[,,1],               name_i,             -1, 0, 1)
    g_list[[length(g_list) + 1]] <- visualize_matrix_heatmap(g_ii_suffix_decomp[,,1],        name_i_reconstruct, -1, 0, 1)
                
  }
  
  
  return(g_list)

}

# orthogonality of eigenfunctions
result_43_prep <- function(step_4){
  
  g_list <- list()
  
  # 0) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(step_4), 'eigen_decomp')
  
  queried_names <- names(step_4)
  
  for(i in 1:length(suffix_names)){
    mat_ortho <- t(step_4[[queried_names[i]]]$eigenfunctions[[1]]) %*% step_4[[queried_names[i]]]$eigenfunctions[[1]]
    
    g_list[[length(g_list) + 1]] <- visualize_matrix_heatmap(mat_ortho, suffix_names[i], zmid = 0)
  }
  return(g_list)
}

# reconstruction error histogram
result_44_prep <- function(step_3, step_4, p, X_truth, eigen_troubleshoot){
  
  if(! eigen_troubleshoot){
    
  
    # prep
    g_ii_truth          <- prep_eigendecomposition_ii(step_3$g_ij_truth, p)
    g_ii_est            <- prep_eigendecomposition_ii(step_3$g_ij_est, p) 
    
    # validate
    g_ii_truth_decomp           <- validate_eigendecomposition_ii(g_ii_truth,          step_4$eigen_decomp_truth)
    g_ii_est_decomp             <- validate_eigendecomposition_ii(g_ii_est,            step_4$eigen_decomp_est) 
    
    error_list <- list(error_truth = as.numeric(g_ii_truth - g_ii_truth_decomp),
                       error_est   = as.numeric(g_ii_est - g_ii_est_decomp))
    
    if(X_truth){
      #prep and validate X_truth
      g_ii_X_truth        <- prep_eigendecomposition_ii(step_3$g_ij_X_truth, p) 
      g_ii_X_truth_decomp         <- validate_eigendecomposition_ii(g_ii_X_truth,        step_4$eigen_decomp_X_truth) 
      
      error_list[['error_X_truth']] <- as.numeric(g_ii_X_truth - g_ii_X_truth_decomp)
    } 
    
    return(error_list)
  } else{
    step_3_names <- names(step_3)
    step_4_names <- names(step_4)
    
    # do truth first and get it out of the way
    g_ii_truth         <- prep_eigendecomposition_ii(step_3$g_ij_truth, p)
    g_ii_truth_decomp  <- validate_eigendecomposition_ii(g_ii_truth,          step_4$eigen_decomp_truth)
    error_list         <- list(error_truth = as.numeric(g_ii_truth - g_ii_truth_decomp))
    
    step_4_names <- setdiff(step_4_names, 'eigen_decomp_truth')
    step_3_names <- setdiff(step_3_names, 'g_ij_truth')
    suffix_names <- step_00_grab_ID(step_4_names, prefix = 'eigen_decomp')
    
    
    step_3_names <- rep(step_3_names, each = 3) # for each eig1 eig2 eig3
    
    for(idx in 1:length(step_3_names)){
      g_ii_idx         <- prep_eigendecomposition_ii(step_3[[step_3_names[idx]]], p)
      g_ii_idx_decomp  <- validate_eigendecomposition_ii(g_ii_idx,          step_4[[step_4_names[idx]]])
      
      
      error_name <- paste0('error_', suffix_names[idx])
      error_list[[error_name]] <- as.numeric(g_ii_idx - g_ii_idx_decomp)
    }
    
    return(error_list)
    
  }
}

result_45_prep <- function(step_4){
  
  # keep eigenvalues of at most 5 processes
  step_4_evals <- lapply(step_4, function(x) x$eigenvalues[1:min(5, length(x$eigenvalues))])
  
  
  # reorganize dataframe of eigenvalues 
  df <- data.frame()
  
  for (outer_name in names(step_4_evals)) {
    inner_list <- step_4_evals[[outer_name]]
    
    for (j in seq_along(inner_list)) {
      vals <- inner_list[[j]]
      temp <- data.frame(
        description = outer_name,
        process = j,
        eigen_index = seq_along(vals),
        value = vals
      )
      df <- rbind(df, temp)
    }
  }
  
  # Plot: one plot per component, colored by process
  
  g_list <- list()
  
  for(i in 1:5){
    g_list[[i]] <- ggplot(data = df[df$process == i,], aes(x = eigen_index, y = value, color = description)) +
      geom_line() + 
      geom_point() + 
      labs(x = "Eigenvalue index", y = "Eigenvalue", color = "Quantity") +
      ggtitle(paste0('Process ', i)) + 
      theme_minimal() 
  }
  
  return(g_list)
  
  
}

result_46_prep <- function(step_3, step_4){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: For each process, check if Phi_a %*% G_ii %*% Phi_b = 0 for a \neq b
  # 
  #       We want to check if the covariance between different eigencomponents of the same process is indeed 0
  #
  # inputs:
  #
  # - step_3    (list) list of g_ij_suffix
  #   - (list of m x m matrices for i_j entries)
  # 
  # - step_4    (list) list of eigen_decomp_suffix with the following 3 entries
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  # 
  # ----------------------------------------------------------------------------
  
  # 1) arrange names due to eig_k
  step_3_sub <- step_3[names(step_3) != "g_ij_truth"]
  s3_names <- rep(names(step_3_sub), each = 3)
  s3_names <- c('g_ij_truth', s3_names)
  
  s4_names <- names(step_4)
  
  # 2) iterate 
  results_list <- list()
  off_diag_values <- c()
  
  # Iterate through the aligned scenarios
  for (idx in seq_along(s4_names)) {
    
    s3_key <- s3_names[idx]
    s4_key <- s4_names[idx]
    
    G_list <- step_3[[s3_key]]        # List of all i_j matrices (m x m)
    E_list <- step_4[[s4_key]][[2]]   # List of p matrices (m x d_i)
    p <- length(step_4[[s4_key]][[1]])
    
    # Store off-diagonal checks for each process i
    process_checks <- list()
    
    for (i in 1:p) {
      G_ii <- G_list[[paste0(i, "_", i)]]
      Phi  <- E_list[[i]] # m x d_i
      
      # The Quadratic Form: d_i x d_i matrix
      # Represents <Phi_a, G Phi_b>
      quad_form <- t(Phi) %*% G_ii %*% Phi
      
      # Extract off-diagonal elements
      if (ncol(quad_form) > 1) {
        off_diag_vals <- quad_form[row(quad_form) != col(quad_form)]
        max_off_diag  <- max(abs(off_diag_vals))
      } else {
        max_off_diag  <- 0 # Only one component, no off-diagonals to check
      }
      
      process_checks[[i]] <- list(
        matrix_check = quad_form,
        max_error = max_off_diag,
        is_near_zero = all.equal(max_off_diag, 0, tolerance = 1e-8)
      )
    }
    
    results_list[[s4_key]] <- process_checks
  }
  
  return(results_list)
  
}

# plot the (dxd) KL covariance values of the (i, j) block
result_55_prep <- function(step_5, i, j){
  
  key <- paste0(i, '_', j)
  
  g_list <- list()
  
  # 0) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(step_5), 'KL_cov')
  
  
  queried_names <- names(step_5)
  
  # 2) for each suffix name, do it 
  for(i in 1:length(suffix_names)){
    
    g_list[[length(g_list) + 1]] <- visualize_matrix_heatmap(step_5[[queried_names[i]]][[key]], g_title = suffix_names[i], zmid = 0)

    
  }
  
  return(g_list)
  
}

# plot the (dxd) KL correlation values of the (i, j) block
result_56_prep <- function(step_5b, i, j){
  
  
  key <- paste0(i, '_', j)
  
  g_list <- list()
  
  # 0) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(step_5b), 'KL_cor')
  
  
  queried_names <- names(step_5b)
  
  # 2) for each suffix name, do it 
  for(i in 1:length(suffix_names)){
    
    g_list[[length(g_list) + 1]] <- visualize_matrix_heatmap(step_5b[[queried_names[i]]][[key]], g_title = suffix_names[i], zmid = 0)
    
    
  }
  
  return(g_list)
  
}

# plot the (dxd) KL precision values of the (i, j) block
result_57_prep <- function(step_5c, i, j){
  
  
  key <- paste0(i, '_', j)
  
  g_list <- list()
  
  # 0) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(step_5c), 'KL_prec')
  
  
  queried_names <- names(step_5c)
  
  # 2) for each suffix name, do it 
  for(i in 1:length(suffix_names)){
    
    g_list[[length(g_list) + 1]] <- visualize_matrix_heatmap(step_5c[[queried_names[i]]][[key]], g_title = suffix_names[i], zmid = 0)
    
    
  }
  
  return(g_list)
  
}

# plot the assembled KL covariance values of all blocks
result_58_cov_prep <- function(step_5, p){
  
  g_list <- list()
  
  # 0) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(step_5), 'KL_cov')
  
  
  queried_names <- names(step_5)
  
  # 1) global max and min
  vals <- unlist(step_5, recursive = TRUE, use.names = FALSE)
  
  global_min <- min(vals, na.rm = TRUE)
  global_max <- max(vals, na.rm = TRUE)
  
  # 2) for each suffix name, do it 
  for(i in 1:length(suffix_names)){
    
    assembled_items <- assemble_block_matrix_irregular(step_5[[queried_names[i]]], p)
    
    assembled_cor <- assembled_items$block_matrix
    
    
    
    g_list[[length(g_list) + 1]] <- visualize_matrix_heatmap(assembled_items$block_matrix, 
                                                             g_title = suffix_names[i], 
                                                             x_max_borders = assembled_items$col_borders,
                                                             y_max_borders = assembled_items$row_borders,
                                                             zmid = 0)
    
    
  }
  
  return(g_list)
  
}

# plot the assembled KL correlation values of all blocks
result_58_prep <- function(step_5b, p){
  
  g_list <- list()
  
  # 0) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(step_5b), 'KL_cor')
  
  
  queried_names <- names(step_5b)
  
  # 1) global max and min
  vals <- unlist(step_5b, recursive = TRUE, use.names = FALSE)
  
  global_min <- min(vals, na.rm = TRUE)
  global_max <- max(vals, na.rm = TRUE)
  
  # 2) for each suffix name, do it 
  for(i in 1:length(suffix_names)){
    
    assembled_items <- assemble_block_matrix_irregular(step_5b[[queried_names[i]]], p)
    
    assembled_cor <- assembled_items$block_matrix
    
    
    
    g_list[[length(g_list) + 1]] <- visualize_matrix_heatmap(assembled_items$block_matrix, 
                                                             g_title = suffix_names[i], 
                                                             x_max_borders = assembled_items$col_borders,
                                                             y_max_borders = assembled_items$row_borders,
                                                             zmid = 0)
    
    
  }
  
  return(g_list)
  
}

# plot the assembled KL precision values of all blocks
result_59_prep <- function(step_5c, p){
  
  g_list <- list()
  
  # 0) borrow from 12z - get suffix names
  suffix_names <- step_00_grab_ID(names(step_5c), 'KL_prec')
  
  queried_names <- names(step_5c)
  
  # 2) for each suffix name, do it 
  for(i in 1:length(suffix_names)){
    
    assembled_items <- assemble_block_matrix_irregular(step_5c[[queried_names[i]]], p)
    
    assembled_cor <- assembled_items$block_matrix
    
    g_list[[length(g_list) + 1]] <- visualize_matrix_heatmap(assembled_items$block_matrix, 
                                                             g_title = suffix_names[i], 
                                                             x_max_borders = assembled_items$col_borders,
                                                             y_max_borders = assembled_items$row_borders,
                                                             zmid = 0)
  }
  
  return(g_list)
  
}

result_80s_prep <- function(step_8, m, m_est, i, j){
  
  #   - V_cond_ground_truth_full          (pm x pm matrix)
  #   - V_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - V_cond_truth_full                 (pm x pm matrix)
  #   - V_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - V_cond_X_truth_full               (pm x pm matrix)
  #   - V_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - V_cond_est_full                   (pm_est x pm_est matrix)  
  
  g_list <- list(visualize_matrix_heatmap(extract_block_structure_ij(step_8[[1]], m,     1, 2), 'Ground Truth'), 
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[2]], m_est, 1, 2), 'Coarse Ground Truth'), 
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[3]], m,     1, 2), 'Truth Theory'),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[4]], m_est, 1, 2), 'Coarse Truth Theory'),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[5]], m,     1, 2), 'Truth X'), 
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[6]], m_est, 1, 2), 'Coarse Truth X'), 
                 visualize_matrix_heatmap(extract_block_structure_ij(step_8[[7]], m_est, 1, 2), 'Estimate')) 
  
  return(g_list)
  
}

result_90s_prep_ij <- function(step_9, m, m_est, i, j, full = T){
  
  # list of heatmaps of C_Xi_Xj 
  # i = 1
  # j = 2
  
  if(! full){
    g_list <- list(visualize_matrix_heatmap(extract_block_structure_ij(step_9$C_cond_truth_full, m,     i, j), 'Truth',               zmid = 0),
                   visualize_matrix_heatmap(extract_block_structure_ij(step_9$C_cond_est_full,   m_est, i, j), 'Estimate',            zmid = 0))
    
    return(g_list)
  }
  
  g_list <- list(visualize_matrix_heatmap(extract_block_structure_ij(step_9[[1]], m,     i, j), 'Ground Truth',           zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[2]], m_est, i, j), 'Coarse Ground Truth',    zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[8]], m,     i, j), 'Ground Truth v2',        zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[9]], m_est, i, j), 'Coarse Ground Truth v2', zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[3]], m,     i, j), 'Truth Theory',           zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[4]], m_est, i, j), 'Coarse Truth Theory',    zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[5]], m,     i, j), 'Truth X',                zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[6]], m_est, i, j), 'Coarse Truth X',         zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_9[[7]], m_est, i, j), 'Estimate',               zmid = 0))
  
  return(g_list)
}

result_90s_prep_pm <- function(step_9, m, m_est, full = T){
  
  if(! full){
    g_list <- list(visualize_matrix_heatmap(step_9$C_cond_truth_full,        'Truth',            zmid = 0),
                   visualize_matrix_heatmap(step_9$C_cond_truth_unnorm_full, 'Truth Unnorm',     zmid = 0),
                   visualize_matrix_heatmap(step_9$C_cond_est_full,          'Estimate',         zmid = 0),
                   visualize_matrix_heatmap(step_9$C_cond_est_unnorm_full,   'Estimate Unnorm',  zmid = 0))
    return(g_list)
  }
  
  g_list <- list(visualize_matrix_heatmap(step_9[[1]], 'Ground Truth',           zmid = 0),
                 visualize_matrix_heatmap(step_9[[2]], 'Coarse Ground Truth',    zmid = 0),
                 visualize_matrix_heatmap(step_9[[8]], 'Ground Truth v2',        zmid = 0),
                 visualize_matrix_heatmap(step_9[[9]], 'Coarse Ground Truth v2', zmid = 0),
                 visualize_matrix_heatmap(step_9[[3]], 'Truth Theory',           zmid = 0),
                 visualize_matrix_heatmap(step_9[[4]], 'Coarse Truth Theory',    zmid = 0),
                 visualize_matrix_heatmap(step_9[[5]], 'Truth X',                zmid = 0),
                 visualize_matrix_heatmap(step_9[[6]], 'Coarse Truth X',         zmid = 0),
                 visualize_matrix_heatmap(step_9[[7]], 'Estimate',               zmid = 0))
  
  
  
  return(g_list)
}

result_95_prep <- function(step_11, m, m_est, p, full = T){
  
  
  if(! full){
    C_HS_truth                 <- step_11$C_HS_truth
    C_HS_est                   <- step_11$C_HS_est
    
    C_HS_truth_unnorm          <- step_11$C_HS_truth_unnorm
    C_HS_est_unnorm            <- step_11$C_HS_est_unnorm
    
    diag(C_HS_truth) <- 0
    diag(C_HS_est) <- 0
    diag(C_HS_truth_unnorm) <- 0
    diag(C_HS_est_unnorm) <- 0
    
    g_list <- list(visualize_matrix_heatmap(C_HS_truth,         'Truth',            zmid = 0),
                   visualize_matrix_heatmap(C_HS_truth_unnorm,  'Truth Unnorm',     zmid = 0),
                   visualize_matrix_heatmap(C_HS_est,           'Estimate',         zmid = 0),
                   visualize_matrix_heatmap(C_HS_est_unnorm,    'Estimate Unnorm',  zmid = 0)) 
    
    return(g_list)
  }
  
  C_HS_ground_truth_full          <- hilbert_schmidt_norm_pm(step_9[[1]], p, m) 
  C_HS_coarse_ground_truth_full   <- hilbert_schmidt_norm_pm(step_9[[2]], p, m_est)
  C_HS_truth_full                 <- hilbert_schmidt_norm_pm(step_9[[3]], p, m)
  C_HS_coarse_truth_full          <- hilbert_schmidt_norm_pm(step_9[[4]], p, m_est)
  C_HS_X_truth_full               <- hilbert_schmidt_norm_pm(step_9[[5]], p, m)
  C_HS_X_coarse_truth_full        <- hilbert_schmidt_norm_pm(step_9[[6]], p, m_est)
  C_HS_est                        <- hilbert_schmidt_norm_pm(step_9[[7]], p, m_est)
  
  diag(C_HS_ground_truth_full) <- 0
  diag(C_HS_coarse_ground_truth_full) <- 0
  diag(C_HS_truth_full) <- 0
  diag(C_HS_coarse_truth_full) <- 0
  diag(C_HS_X_truth_full) <- 0
  diag(C_HS_X_coarse_truth_full) <- 0
  diag(C_HS_est) <- 0
  
  g_list <- list(visualize_matrix_heatmap(C_HS_ground_truth_full,        'Ground Truth',        zmid = 0),
                 visualize_matrix_heatmap(C_HS_coarse_ground_truth_full, 'Coarse Ground Truth', zmid = 0),
                 visualize_matrix_heatmap(C_HS_truth_full,               'Truth Theory',        zmid = 0),
                 visualize_matrix_heatmap(C_HS_coarse_truth_full,        'Coarse Truth Theory', zmid = 0),
                 visualize_matrix_heatmap(C_HS_X_truth_full,             'Truth X',             zmid = 0),
                 visualize_matrix_heatmap(C_HS_X_coarse_truth_full,      'Coarse Truth X',      zmid = 0),
                 visualize_matrix_heatmap(C_HS_est,                      'Estimate',            zmid = 0))        
  
  return(g_list)
  
}



# only view the (i, j)-th block
result_100s_prep_ij <- function(step_10, m, m_est, i, j){
  #   - P_cond_ground_truth_full          (pm x pm matrix)
  #   - P_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - P_cond_truth_full                 (pm x pm matrix)
  #   - P_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - P_cond_X_truth_full               (pm x pm matrix)
  #   - P_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - P_cond_est_full                   (pm_est x pm_est matrix)    
  g_list <- list(visualize_matrix_heatmap(extract_block_structure_ij(step_10[[1]], m,     i, j), 'Truth Ground',         zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[2]], m_est, i, j), 'Coarse Truth Ground',  zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[3]], m,     i, j), 'Truth Theory',         zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[4]], m_est, i, j), 'Coarse Truth Theory',  zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[5]], m,     i, j), 'Truth X',              zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[6]], m_est, i, j), 'Coarse Truth X',       zmid = 0),
                 visualize_matrix_heatmap(extract_block_structure_ij(step_10[[7]], m_est, i, j), 'Estimate',             zmid = 0))
  
  return(g_list)
}

# view entire pm matrix
result_100s_prep_pm <- function(step_10, m, m_est, full = T){
  
  if(!full){
    g_list <- list(visualize_matrix_heatmap(step_10$P_cond_truth_full,        'Truth',            zmid = 0),
                   visualize_matrix_heatmap(step_10$P_cond_truth_unnorm_full, 'Truth Unnorm',     zmid = 0),
                   visualize_matrix_heatmap(step_10$P_cond_est_full,          'Estimate',         zmid = 0),
                   visualize_matrix_heatmap(step_10$P_cond_est_unnorm_full,   'Estimate Unnorm',  zmid = 0))
    
    
    
    return(g_list) 
  }
  
  g_list <- list(visualize_matrix_heatmap(step_10[[1]], 'Ground Truth',           zmid = 0),
                 visualize_matrix_heatmap(step_10[[2]], 'Coarse Ground Truth',    zmid = 0),
                 # visualize_matrix_heatmap(step_10[[8]], 'Ground Truth v2',        zmid = 0),
                 # visualize_matrix_heatmap(step_10[[9]], 'Coarse Ground Truth v2', zmid = 0),
                 visualize_matrix_heatmap(step_10[[3]], 'Truth Theory',           zmid = 0),
                 visualize_matrix_heatmap(step_10[[4]], 'Coarse Truth Theory',    zmid = 0),
                 visualize_matrix_heatmap(step_10[[5]], 'Truth X',                zmid = 0),
                 visualize_matrix_heatmap(step_10[[6]], 'Coarse Truth X',         zmid = 0),
                 visualize_matrix_heatmap(step_10[[7]], 'Estimate',               zmid = 0))
  
  
  
  return(g_list)
}

result_112_prep <- function(step_11, remove_diag, full = T){
  #   - w_mat_ground_truth        (p x p)   matrix of HS norms of the pm x pm ground truth
  #   - w_mat_coarse_ground_truth (p x p)   matrix of HS norms of the pm_est x pm_est ground truth 
  #   - w_mat_X_coarse_truth      (p x p)   matrix of HS norms of ...
  #   - w_mat_X_truth             (p x p)   matrix of HS norms of ...
  #   - w_mat_coarse_truth        (p x p)   matrix of HS norms of ...
  #   - w_mat_truth               (p x p)   matrix of HS norms of ...
  #   - w_mat_est                 (p x p)   matrix of HS norms of ...  
  

  
  if(remove_diag){
    step_11 <- lapply(step_11, function(x) {
      diag(x) <- 0
      x
    })
  }
  
  if(! full){
    g_list <- list(visualize_matrix_heatmap(step_11$w_mat_truth,        'Truth',            zmid = 0), 
                   visualize_matrix_heatmap(step_11$w_mat_truth_unnorm, 'Truth Unnorm',     zmid = 0),
                   visualize_matrix_heatmap(step_11$w_mat_est ,         'Estimate',         zmid = 0), 
                   visualize_matrix_heatmap(step_11$w_mat_est_unnorm,   'Estimate Unnorm',  zmid = 0))
    
    return(g_list) 
  }
  
  g_list <- list(visualize_matrix_heatmap(step_11[[1]], 'Ground Truth',        zmid = 0), 
                 visualize_matrix_heatmap(step_11[[2]], 'Coarse Ground Truth', zmid = 0), 
                 visualize_matrix_heatmap(step_11[[3]], 'Truth Theory',        zmid = 0),   
                 visualize_matrix_heatmap(step_11[[4]], 'Coarse Truth Theory', zmid = 0), 
                 visualize_matrix_heatmap(step_11[[5]], 'Truth X',             zmid = 0), 
                 visualize_matrix_heatmap(step_11[[6]], 'Coarse Truth X',      zmid = 0),
                 visualize_matrix_heatmap(step_11[[7]], 'Estimate',            zmid = 0))
  
  return(g_list)
  
}

result_113_prep <- function(step_12){
  
  #   - roc_ground_truth           (list of roc outputs)
  #   - roc_coarse_ground_truth
  #   - roc_truth
  #   - roc_coarse_truth
  #   - roc_X_truth
  #   - roc_X_coarse_truth
  #   - roc_est
  
  # ROC plot - done on qrsh 
  roc_graphs <- list()
  
  title_names <- c('Ground Truth', 'Coarse Ground Truth', 'Truth Theory', 'Coarse Truth Theory', 'X Truth', 'X Coarse Truth',
                   'Estimate')
  
  for(i in 1:7){
    
    roc_df          <- step_12[[i]]$roc_df
    ideal_spec      <- step_12[[i]]$specificity
    ideal_sens      <- step_12[[i]]$sensitivity
    ideal_threshold <- step_12[[i]]$threshold
    auc_value       <- step_12[[i]]$auc
    
    roc_graphs[[i]] <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
      geom_step(direction = "vh", color = "blue", size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
      labs(title = title_names[i], x = "False Positive Rate", y = "True Positive Rate") +
      annotate("point", x = 1 - ideal_spec, y = ideal_sens, color = "red", size = 3) +
      annotate("text", x = 1 - ideal_spec, y = ideal_sens, 
               label = paste0("Threshold=", round(ideal_threshold, 3)),
               hjust = -0.1, vjust = -0.5, color = "red") +
      
      annotate("text", x = 0.6, y = 0.2,            # position for AUC label
               label = paste0("AUC = ", round(auc_value, 3)),
               color = "darkgreen", size = 5) +    
      theme_minimal()         
  }
  
  return(roc_graphs)
  
}

result_arr_mat <- function(g_list, grob_caption, arr_mat_i){
  
  # reused function to arrange all graphs in an arr_mat_8 fashion
  
  nrows <- dim(arr_mat_i)[1]
  ncols <- dim(arr_mat_i)[2]
  
  if(nrows == 2 & ncols == 2){
    g <- grid.arrange(g_list[[1]], 
                      g_list[[2]], 
                      textGrob(grob_caption, gp = gpar(fontsize = 14)),
                      g_list[[3]], layout_matrix = arr_mat_i)
  }
  
  if(nrows == 2 & ncols == 3){
    
    if(length(g_list) == 4){
      g <- grid.arrange(g_list[[1]], 
                        g_list[[2]], 
                        g_list[[3]], 
                        g_list[[4]], 
                        textGrob(grob_caption, gp = gpar(fontsize = 14)),
                        layout_matrix = arr_mat_i)      
    } else{
      g <- grid.arrange(g_list[[1]], 
                        g_list[[2]], 
                        g_list[[3]], 
                        g_list[[4]], 
                        textGrob(grob_caption, gp = gpar(fontsize = 14)),
                        g_list[[5]], 
                        layout_matrix = arr_mat_i)      
    }
    

  }  
  
  if(nrows == 2 & ncols == 4){
    
    if(length(g_list) == 6){
      g <- grid.arrange(g_list[[1]], g_list[[2]], g_list[[3]],
                        g_list[[4]], g_list[[5]], g_list[[6]],
                        textGrob(grob_caption, gp = gpar(fontsize = 14)),
                        layout_matrix = arr_mat_i)       
    } else{

      # often times, it's ground truth (1, 2) truth theory (3, 4) X truth (5, 6) and estimate (7)
      
      g <- grid.arrange(g_list[[1]], g_list[[2]], g_list[[3]],
                        g_list[[4]], g_list[[5]], g_list[[6]],
                        textGrob(grob_caption, gp = gpar(fontsize = 14)),
                        g_list[[7]],
                        layout_matrix = arr_mat_i)   
    }
  }
  
  if(nrows == 2 & ncols == 5){
    # often times, it's ground truth (1, 2) ground truth v2 (3, 4) truth theory (5, 6) X truth (7, 8) and estimate (9)
    g <- grid.arrange(g_list[[1]], g_list[[2]], g_list[[3]], g_list[[4]],
                      g_list[[5]], g_list[[6]], g_list[[7]], g_list[[8]],
                      textGrob(grob_caption, gp = gpar(fontsize = 14)),
                      g_list[[9]], 
                      layout_matrix = arr_mat_i)    
  }
  

  
  return(g)
}



