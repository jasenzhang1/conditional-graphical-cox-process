# I have a bunch of .rda results

# now I want to visualize them in a graph

library(ggplot2)
library(igraph) # visualize graphs

adjacency_heatmap <- function(adj_mat, ID, region_border=NULL){
  
  # adj_mat = matrix of 0's and 1's
  # ID = name of mouse for the title
  # region_border = neuron num cutoff between HIP and EHC
  
  n <- dim(adj_mat)[1]
  m2_melt <- reshape2::melt(adj_mat)

  
  border_id3 <- 0:n + 0.5
  
  g3 <- ggplot(m2_melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Column", y = "Row", fill = "Value") +
    scale_y_reverse() +  # reverse y-axis scale
    geom_vline(xintercept = border_id3, alpha = 0.3) +  # neuron grid
    geom_hline(yintercept = border_id3, alpha = 0.3) +  # neuron grid
    
    # titles and theme
    xlab('Neuron') + 
    ylab('Neuron') + 
    ggtitle(ID) +
    theme_minimal() +
    theme(legend.position = "none")
    
  
  if(!is.null(region_border)){
    g3 <- g3 + 
      geom_vline(xintercept = region_border, color = 'dodgerblue', alpha = 0.6) + 
      geom_hline(yintercept = region_border, color = 'dodgerblue', alpha = 0.6)     
  }
  
  g3
}

insert_na_symmetric <- function(mat, indices) {
  
  # reconstruct the adjacency matrix with all neurons, since we removed some originally
  # we insert NA rows and columns in ascending order of neuron ID number
  
  # mat = truncated adjacency matrix
  # indices = sorted array of neuron numbers that were removed during estimation
  
  # print('starting conditions')
  # print(dim(mat)[1])
  # print(indices)
  # print(dim(mat)[1] + length(indices))
  
  if(length(indices) == 0){
    return(mat)
  }
  
  for(index in indices){
    
    d_i <- dim(mat)[1]
    
    # print('index and d_i')
    # print(index)
    # print(d_i)    
    
    if(index <= d_i){ # case 1: missing neuron is inserted inside the matrix
    
      mat_1 <- rbind( mat[1:(index-1), ],   # first rows
                      NA,                   # pad NA
                      mat[index:d_i, ])     # last rows
      
      mat_2 <- cbind(mat_1[, 1:(index-1)],  # first columns
                     NA,                    # pad NA
                     mat_1[,index:d_i])     # last columns
      
      mat <- mat_2
      
    } else if(index == d_i + 1){ # case 2: missing neuron is the very next one, so we append
      
      mat_1 <- rbind(mat, NA)    # pad a row
      mat_2 <- cbind(mat_1, NA)  # pad a column
      mat <- mat_2
      
    } else{ # case 3: the index is way beyond d_i
      stop('ERROR IN insert_na_symmetric')
    }
    
  }
  
  return(mat)
}

get_network <- function(adj_mat, ID){
  
  # visualize the graph with nodes and edges
  # only keep vertices that have edges
  #
  # adj_mat = adjacency matrix
  # ID = 
  
  diag(adj_mat) <- 0
  verts <- which(apply(adj_mat,2,function(x){sum(x, na.rm = T)}) > 0)
  
  adj_mat2 <- adj_mat[verts, verts]
  colnames(adj_mat2) <- verts
  rownames(adj_mat2) <- verts
  
  
  g <- graph_from_adjacency_matrix(adj_mat2, mode = "undirected")  # or "directed"
  
  return(g)
  
  # plot(g, 
  #      vertex.label = V(g)$name,  # Optional: show vertex labels
  #      vertex.color = "lightblue",
  #      edge.color = "gray",
  #      edge.arrow.size = 0.5,
  #      main = ID)  
}

get_summary_statistics <- function(adj_mat_fill, adj_mat_og){
  
  # number statistics
  num_edges <- 0.5 * (sum(adj_mat_og) - dim(adj_mat_og)[1])     # number of edges            (121)
  num_neurons <- dim(adj_mat_fill)[1]                           # number of neurons in total (169)
  num_NA <- sum(is.na(diag(adj_mat_fill)))                      # number of NA neurons       (12)
  num_islands <- length(which(apply(adj_mat_og,2,sum) == 1))    # number of island neurons   (102)
  
  # simple statistics
  avg_deg <- round(2 * num_edges/num_neurons, 2)
  
  adj_mat2 <- adj_mat_og
  diag(adj_mat2) <- 0
  verts <- which(apply(adj_mat2,2,sum) > 0)
  num_con_verts <- length(verts)                                # number of connected neurons (55)
  
  
  
  # distribution of degrees
  adj_mat3 <- adj_mat2[verts, verts]
  
  degs <- table(apply(adj_mat3, 2, sum))
  
  # number of components
  g <- graph_from_adjacency_matrix(adj_mat3, mode = "undirected")
  num_comps <- components(g)$no  
  
  sum_stats <- list()
  sum_stats[['num_neurons']] <- num_neurons
  sum_stats[['num_NA']] <- num_NA
  sum_stats[['num_islands']] <- num_islands
  sum_stats[['num_con_verts']] <- num_con_verts
  sum_stats[['num_edges']] <- num_edges
  sum_stats[['num_comps']] <- num_comps
  sum_stats[['degs']] <- degs
  sum_stats[['verts']] <- verts
  sum_stats[['avg_deg']] <- avg_deg
  
  
  
  return(sum_stats)
  
}

save_summary_statistics <- function(figures_root, figure_name, settings_kept, sum_stats){
  
  txt_name <- paste(figure_name, '.txt', sep = '_summary')
  txt_name <- paste(figures_root, txt_name, sep = '')
  
  sink(txt_name)
  
  for(i in 1:length(sum_stats)){
    
    cat('-----------------------------------------\n')
    
    cat('Settings:\n')
    
    cat('Mouse: ')
    cat(settings_kept[[i]][['mouse']])
    cat('\n')
    cat('Week/epoch num: ')
    cat(settings_kept[[i]][['time']])
    cat('\n')
    cat('Movement: ')
    cat(settings_kept[[i]][['movement']])
    cat('\n') 
    cat('VR: ')
    cat(settings_kept[[i]][['VR']])
    cat('\n')
    
    cat('-----------------------------------------\n')
    
    cat('Total Neurons:\t')
    cat(sum_stats[[i]][['num_neurons']])
    cat('\n')
    cat('Discarded Neurons:\t')
    cat(sum_stats[[i]][['num_NA']])
    cat('\n')
    cat('Island Neurons:\t')
    cat(sum_stats[[i]][['num_islands']])
    cat('\n')
    cat('Connected Neurons:\t')
    cat(sum_stats[[i]][['num_con_verts']])
    cat('\n')   
    cat('Number of Edges:\t')
    cat(sum_stats[[i]][['num_edges']])
    cat('\n')  
    cat('Average Degree:\t')
    cat(sum_stats[[i]][['avg_deg']])
    cat('\n')      
    cat('Number of Components:\t')
    cat(sum_stats[[i]][['num_comps']])
    cat('\n\n')  
    cat('Degree Distribution: \n')
    cat('Degree:\t')
    cat(names(sum_stats[[i]][['degs']])) 
    cat('\n')
    cat('Count:\t')
    cat(unname(sum_stats[[i]][['degs']]))  
    
    
    cat('\n\n--------------\n\n')  
  }
  
  sink()  
}

unpack <- function(n, movement, VR, epoch_or_week, ew_num, min_edges, df_brain_region, first_root){
  
  
  
  # 1) obtain file directory
  if(epoch_or_week == 'week'){
    setting_ID <- paste('w', '_m', sep = as.character(ew_num))          # w1_m
  } else{
    setting_ID <- paste('e', '_m', sep = as.character(ew_num))          # e1_m
  }

  setting_ID <- paste(setting_ID, 'vr', sep = as.character(movement))    # e1_m2vr
  setting_ID <- paste(setting_ID, '_n', sep = as.character(VR))          # e1_m2vr0_n
  setting_ID <- paste(setting_ID, as.character(n), sep = '')             # e1_m2vr0_n400
  setting_ID <- paste(setting_ID, as.character(min_edges), sep = '_me')  # e1_m2vr0_n400_me1
  
  results_root <- paste(first_root, #"/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/"
                        '/',
                        sep = setting_ID)
  
  figures_root <- paste(first_root, #"/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/"
                        '/',
                        sep = 'figures')  
  
  # if the folder doesn't exist, create it 
  if (!dir.exists(figures_root)) {
    dir.create(figures_root)
  }     
  
  # 2) prepare to store results
  graphs <- list()
  g_vis <- list()
  edge_vec <- c()
  neuron_vec <- c()
  sum_stats <- list()
  
  file_names <- list.files(path = results_root, full.names = TRUE)
  
  print(file_names)
  
  m <- length(file_names)
  no_ext <- sub("\\.[^.]*$", "", file_names)
  rda_names  <- sub(".*/", "", no_ext)                     # rda names
  mice_names <-  sub(".*/([^/_]+)_.*", "\\1", file_names)  # get mice names, between last '/' and very next '_'
  
  for(i in 1:m){
    
    
    # 2) loading
    print(file_names[i])
    load(file_names[i])
    
    # 3) fill in NA's
    adj_mat_og <- graph_all[['GPP_BIC']]
    missing_neurons <- sort(graph_all[['missing_neurons']])
    adj_mat_fill <- insert_na_symmetric(adj_mat_og, missing_neurons)
    
    # 4) graph and store
    region_border <- max(df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == mice_names[i]]) + 0.5
    
    g_i <- adjacency_heatmap(adj_mat_fill, mice_names[i], region_border=region_border)
    
    graphs[[i]] <- g_i
    
    # 5) visualization of network
    g_vis[[i]] <- get_network(adj_mat_fill, rda_names[i])
    
    # 5) summary statistics
    
    sum_stats[[i]] <- get_summary_statistics(adj_mat_fill, adj_mat_og)
    
  }
  
  # save heatmap ---------------------------------------------------
  
  pdf_name <- paste(setting_ID, '.pdf', sep = '_heatmap')
  pdf_name <- paste(figures_root, pdf_name, sep = '')
  pdf(pdf_name, width = 10, height = 10)
  
  for(i in 1:m){
    print(graphs[[i]])
  }
  
  dev.off()
  
  # save graph visualization --------------------------------------
  
  pdf_name <- paste(setting_ID, '.pdf', sep = '_graph_vis')
  pdf_name <- paste(figures_root, pdf_name, sep = '')
  pdf(pdf_name, width = 10, height = 10)
  
  for(i in 1:m){
    
    HIP_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == mice_names[i]] 
    EHC_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Entorhinal_Cortex" & df_brain_region$ID2 == mice_names[i]] 
    
    plot_pts <- sum_stats[[i]][['verts']]
    
    plot_pts_color <- rep('gold', length(plot_pts))
    plot_pts_color[plot_pts %in% EHC_ID] <- 'dodgerblue'
    
    plot(g_vis[[i]], 
         vertex.size = 8,
         vertex.label = V(g_vis[[i]])$name,
         vertex.color = plot_pts_color,
         edge.color = "black",
         edge.arrow.size = 0.5,
         main = mice_names[i])
  }
  
  dev.off()
  
  # save summary -------------------------------------------------
  
  txt_name <- paste(setting_ID, '.txt', sep = '_summary')
  txt_name <- paste(figures_root, txt_name, sep = '')
  
  sink(txt_name)
  
  for(i in 1:m){
    
    cat('-----------------------------------------\n')
    
    cat('Mouse ID: ')
    cat(mice_names[i])
    cat('\n\n')
    
    cat('Total Neurons:\t')
    cat(sum_stats[[i]][['num_neurons']])
    cat('\n')
    cat('Discarded Neurons:\t')
    cat(sum_stats[[i]][['num_NA']])
    cat('\n')
    cat('Island Neurons:\t')
    cat(sum_stats[[i]][['num_islands']])
    cat('\n')
    cat('Connected Neurons:\t')
    cat(sum_stats[[i]][['num_con_verts']])
    cat('\n')   
    cat('Number of Edges:\t')
    cat(sum_stats[[i]][['num_edges']])
    cat('\n')      
    cat('Number of Components:\t')
    cat(sum_stats[[i]][['num_comps']])
    cat('\n\n')  
    cat('Degree Distribution: \n')
    cat('Degree:\t')
    cat(names(sum_stats[[i]][['degs']])) 
    cat('\n')
    cat('Count:\t')
    cat(unname(sum_stats[[i]][['degs']]))  
    
    
    cat('\n\n--------------\n\n')  
  }
  
  sink()
  
  
  return(NULL)
}

unpack_over_time <- function(n, movement, VR, epoch_or_week, ew_nums, mouse_ID, min_edges, df_brain_region, first_root){
  
  # 1) prepare to store results
  graphs <- list()
  g_vis <- list()
  edge_vec <- c()
  neuron_vec <- c()
  sum_stats <- list()  
  
  graph_index <- 0
  ew_kept <- c()      # times kept
  settings_kept <- list()
  
  for(ew_num in ew_nums){ # 2) FOR EACH TIMEPOINT

    # 2.1) obtain file directory
    if(epoch_or_week == 'week'){
      setting_ID <- paste('w', '_m', sep = as.character(ew_num))          # w1_m
    } else{
      setting_ID <- paste('e', '_m', sep = as.character(ew_num))          # e1_m
    }
    
    setting_ID <- paste(setting_ID, 'vr', sep = as.character(movement))    # e1_m2vr
    setting_ID <- paste(setting_ID, '_n', sep = as.character(VR))          # e1_m2vr0_n
    setting_ID <- paste(setting_ID, as.character(n), sep = '')             # e1_m2vr0_n400
    setting_ID <- paste(setting_ID, as.character(min_edges), sep = '_me')  # e1_m2vr0_n400_me1
    
    results_root <- paste(first_root, #"/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/"
                          '/',
                          sep = setting_ID)
    
    figures_root <- paste(first_root, #"/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/"
                          '/',
                          sep = 'figures')  
    
    # if the folder doesn't exist, create it 
    if (!dir.exists(figures_root)) {
      dir.create(figures_root)
    }     
    
    # 2.2) Obtain file names
    
    file_names <- list.files(path = results_root, full.names = TRUE)
    
    m <- length(file_names)
    no_ext <- sub("\\.[^.]*$", "", file_names)
    rda_names  <- sub(".*/", "", no_ext)                     # rda names
    mice_names <-  sub(".*/([^/_]+)_.*", "\\1", file_names)  # get mice names, between last '/' and very next '_'  
    
    if(mouse_ID %in% mice_names){ # 3) if this mouse is in this time, examine it
      
      ew_kept <- c(ew_kept, ew_num)
      mouse_index <- which(mice_names == mouse_ID)
      
      # 2) loading
      load(file_names[mouse_index])
      
      # 3) fill in NA's
      adj_mat_og <- graph_all[['GPP_BIC']]
      missing_neurons <- sort(graph_all[['missing_neurons']])
      adj_mat_fill <- insert_na_symmetric(adj_mat_og, missing_neurons)
      
      # 4) graph and store
      graph_index <- graph_index + 1
      region_border <- max(df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == mouse_ID]) + 0.5
      
      g_i <- adjacency_heatmap(adj_mat_fill, mouse_ID, region_border=region_border)
      
      graphs[[graph_index]] <- g_i
      
      # 5) visualization of network
      g_vis[[graph_index]] <- get_network(adj_mat_fill, rda_names[mouse_index])
      
      # 5) summary statistics
      sum_stats[[graph_index]] <- get_summary_statistics(adj_mat_fill, adj_mat_og)
      
      # 6) settings
      settings_kept[[graph_index]] <- list(mouse = mouse_ID,
                                           time = ew_num,
                                           movement = movement,
                                           VR = VR)
      
    }
        
  }

  
  # save heatmap ---------------------------------------------------
  
  figure_name <- paste(mouse_ID, movement, sep = '_m') #Tau1_m2
  figure_name <- paste(figure_name, VR, sep = 'vr') #Tau1_m2vr0
  
  pdf_name <- paste(figure_name, '.pdf', sep = '_heatmap')
  pdf_name <- paste(figures_root, pdf_name, sep = '')
  pdf(pdf_name, width = 10, height = 10)
  
  for(i in 1:graph_index){
    print(graphs[[i]])
  }
  
  dev.off()
  
  # save graph visualization --------------------------------------
  

  
  pdf_name <- paste(figure_name, '.pdf', sep = '_graph_vis')
  pdf_name <- paste(figures_root, pdf_name, sep = '')
  pdf(pdf_name, width = 20, height = 10)
  
  n <- length(g_vis) # start plotting
  plots_per_page <- 3
  
  for (i in seq(1, n, by=plots_per_page)) {
    # Start a new page
    par(mfrow = c(1, 3))
    
    # Plot up to 3 plots in this page
    for (j in i:min(i + plots_per_page - 1, n)) {
      
      HIP_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == mice_names[j]] 
      EHC_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Entorhinal_Cortex" & df_brain_region$ID2 == mice_names[j]] 
      
      plot_pts <- sum_stats[[j]][['verts']]
      
      plot_pts_color <- rep('gold', length(plot_pts))
      plot_pts_color[plot_pts %in% EHC_ID] <- 'dodgerblue'      
      
      
      plot(g_vis[[j]], 
           vertex.size = 8,
           vertex.label = V(g_vis[[j]])$name,
           vertex.color = plot_pts_color,
           edge.color = "black",
           edge.arrow.size = 0.5,
           main = ew_kept[j])
    }
    
    # Fill empty slots if the last page has fewer than 3 plots
    for (k in (j+1):(i + plots_per_page - 1)) {
      plot.new()  # blank plot
    }
  }
  
  dev.off()
  
  # save summary -------------------------------------------------
  
  save_summary_statistics(figures_root, figure_name, settings_kept, sum_stats)  
  
  
  return(NULL)
}