type_1_graphs <- function(summary_df, graph_type, x_graph_title, x_var, y_var, color_var, group_var, color_palette=NULL, color_low='red', color_high='green'){
  
  # ----------------------------------------------------------------------------
  # 
  # goal: produce a graph of statistics over time
  #
  # - graph_type      (integer) 1: discrete color palette, 2: continuous color palette
  # - x_graph_title   (string) display name of the x-axis 
  # - x_var           (string) name of x-axis variable
  # - y_var           (string) name of y-axis variable 
  # - color_var       (string) name of the color variable, can admit NA values, for which they will be black
  # - group_var       (string) name of the grouping variable in facet_wrap 
  # 
  # - color_palette   (named vector) to denote color scheme for color_var
  #
  # ----------------------------------------------------------------------------
  

  
  if(graph_type == 1){
    g <- ggplot(summary_df, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[color_var]])) + 
      geom_point() + 
      geom_line() + 
      facet_wrap(vars(.data[[group_var]])) + 
      ylab(y_var) + 
      xlab(x_graph_title) + 
      theme_bw() + 
      scale_color_manual(values = color_palette) # discrete palette
  } 
  
  if(graph_type == 2){
    g <- ggplot(summary_df, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[color_var]], group = .data[[color_var]])) + 
      geom_point() + 
      geom_line() + 
      facet_wrap(vars(.data[[group_var]])) + 
      ylab(y_var) + 
      xlab(x_graph_title) + 
      theme_bw() + 
      scale_color_gradientn(colors = c(color_low, color_high), na.value = "black") # continuous palette
  }
    
  
  return(g)
}




plot_graph_stats_ts <- function(summary_df, graph_type, x_graph_title, x_var, color_var, group_var, mouse_var, final_results_dir, thresh_value){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: plot graph statistics over time
  #
  # - see how timestamp affects estimates
  # 
  # 
  # Input:
  # 
  # - thresh_value       (number)
  # - x_graph_title      (string)          name of x-axis to visualize 
  # - x_var              (string)          name of variable to plot along x axis
  # - color_var          (string)          name of variable to create different colors, can admit NA, and we make it black
  # - group_var          (string)          name of variable to facet_wrap() group by 
  # - mouse_var          (string)          name of variable to filter and iterate by
  # - final_results_dir  (string)   where to store final results
  # - thresh_value       (numeric)  value for which we chose to threshold by
  # - summary_df         (data.frame) each row captures information about a specific E(y_c, y_d) estimated graph 
  #
  #   - mouse_ID (string)
  #   - movement (0 = rest, 1 = run)
  #   - VR       (0 = off, 1 = on)
  #   - ew_num OR age_norm / ts_norm
  #
  #   - num_neurons
  #   - num_NA
  #   - num_candidates
  #   - num_islands
  #   - num_con_verts
  #   - num_edges
  #   - num_HE_edges
  #   - num_EE_edges
  #   - num_HH_edges
  #
  #   - avg_deg
  #   - avg_deg_HIP
  #   - avg_deg_EHC
  #   - avg_deg_normalized
  #   - avg_deg_HIP_normalized
  #   - avg_deg_EHC_normalized
  #   - num_comps
  #
  #   - avg_node_strength
  #   - avg_edge_strength
  #   - avg_dist
  #   - diameter
  #
  #   - num_zeros
  #   - fiedler_value
  #   - lambda_max
  #   - lambda_median
  #   - lambda_25
  #   - lambda_75
  #   - lambda_max_sym
  #   - lambda_median_sym
  #   - lambda_25_sym
  #   - lambda_75_sym
  #
  # Output:
  # 
  # - saved pdf of all the graphs
  # - save summary_df as well
  #
  # ----------------------------------------------------------------------------
  

  
  # 1) make everything but mouse_ID numeric, make mouse_ID factor
  # if "NA" is found, which is intended if ts_norm isn't used, we set it to NA
  
  summary_df[-which(names(summary_df) == "mouse_ID")] <- lapply(summary_df[-which(names(summary_df) == "mouse_ID")], function(x) as.numeric(na_if(x, "NA")))
  summary_df$mouse_ID <- factor(summary_df$mouse_ID)
  
  
  # 2) factors (movement_factor, VR_factor, m_vr_factor)
  
  summary_df <- get_mvr_interaction(summary_df)
  
  
  # 2.4) more statistics
  
  summary_df$pct_non_singletons <- 100 * summary_df$num_con_verts/summary_df$num_candidates
  summary_df$pct_inter_region_edge <- 100 * summary_df$num_HE_edges/summary_df$num_edges
  summary_df$pct_HH_edge <- 100 * summary_df$num_HH_edges/summary_df$num_edges
  summary_df$pct_EE_edge <- 100 * summary_df$num_EE_edges/summary_df$num_edges  
  summary_df$pct_candidate_neurons <- 100 * summary_df$num_candidates/summary_df$num_neurons
  
  ## 3) graphing the summary statistics over time ------------------------------
  
  mice_strains2 <- c(
    "Tau1" = "lightcoral",
    "Tau2" = "red",
    "Tau3" = "darkred",
    "WT1" = "lightgreen",
    "WT2" = "green",
    "WT3" = "darkgreen"
  )
  

  
  graphs <- list() # group by discrete strata
  graphs_no_ts <- list()
  

  value_names <- c('avg_deg', 'avg_deg_normalized',          # average_degree
                   'avg_deg_HIP', 'avg_deg_HIP_normalized',
                   'avg_deg_EHC', 'avg_deg_EHC_normalized',
                   'pct_non_singletons', 'num_comps',           # singletons and components
                   'pct_inter_region_edge', 'pct_HH_edge', 'pct_EE_edge', 'pct_candidate_neurons', 
                   'num_zeros', 'fiedler_value', 'lambda_max',   # laplacian
                   'lambda_median', 'lambda_25', 'lambda_75', 
                   'lambda_max_sym', 'lambda_median_sym', 'lambda_25_sym', 'lambda_75_sym',
                   'avg_node_strength', 'avg_edge_strength', 'avg_dist', 'diameter')  # fully connected weighted graph metrics

  
  # plot each outcome
  
  mouse_var_levels <- unique(summary_df[,mouse_var])  # vector
  for(k in 1:length(value_names)){
    for(k2 in mouse_var_levels){
      summary_df2 <- summary_df[dplyr::pull(summary_df, iter_var) == k2, ]      # all of one mouse
      
      
      graph_index <- paste0(value_names[k], '_', k2)
      
      graphs[[graph_index]] <-   type_1_graphs(summary_df2, 2, x_graph_title, x_var, value_names[k], color_var, group_var, color_palette=mice_strains2) + ggtitle(graph_index)
      
    }
    
    summary_df3 <- summary_df[is.na(dplyr::pull(summary_df, color_var)), ]  # keep all mice, filter for ts_norm = NA
    graphs_no_ts[[value_names[k]]] <- type_1_graphs(summary_df3, 1, x_graph_title, x_var, value_names[k], mouse_var, group_var, color_palette = mice_strains2)
     
  }

  
  # save data
  if (!dir.exists(final_results_dir)) {
    dir.create(final_results_dir)
  }
  
  # store summary statistics of each fit
  
  write.csv(summary_df, file = paste0(final_results_dir, 'graph_statistics_thresh_', thresh_value, '.csv'), row.names = FALSE)
  
  # graphs
  
  pdf(paste0(final_results_dir, 'graph_statistics_vary_ts_thresh_', thresh_value, '.pdf'), width = 8, height = 3)
  for(name_i in names(graphs)){
    print(graphs[[name_i]])
  }
  dev.off()  
  
  pdf(paste0(final_results_dir, 'graph_statistics_NA_ts_thresh_', thresh_value, '.pdf'), width = 8, height = 3)
  for(name_i in names(graphs_no_ts)){
    print(graphs_no_ts[[name_i]])
  }
  dev.off()    
}
