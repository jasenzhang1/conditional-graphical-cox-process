# I have a bunch of .rda results

# now I want to visualize them in a graph

library(ggplot2)
library(igraph) # visualize graphs

adjacency_heatmap <- function(adj_mat, ID, region_border=NULL){
  
  # adj_mat = matrix of 0's and 1's
  
  n <- dim(adj_mat)[1]
  m2_melt <- reshape2::melt(adj_mat)
  m2_melt$Var1 <- dim(adj_mat)[1] + 1 - m2_melt$Var1
  
  border_id3 <- 0:n + 0.5
  
  g3 <- ggplot(m2_melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Column", y = "Row", fill = "Value") +
    geom_vline(xintercept = border_id3, alpha = 0.3) + 
    geom_hline(yintercept = border_id3, alpha = 0.3) +    
    theme_minimal()  +
    ggtitle(ID)
  
  if(!is.null(region_border)){
    g3 <- g3 + 
      geom_vline(xintercept = region_border, color = 'dodgerblue', alpha = 0.6) + 
      geom_hline(yintercept = n + 1 - region_border, color = 'dodgerblue', alpha = 0.6)     
  }
  
  g3
}

insert_na_symmetric <- function(mat, indices) {
  
  # reconstruct the adjacency matrix with all neurons, since we removed some originally
  
  if(length(indices) == 0){
    return(mat)
  }
  
  for(index in indices){
    
    d_i <- dim(mat)[1]
    mat_1 <- rbind( mat[1:(index-1), ],   # first rows
                    NA,                   # pad NA
                    mat[index:d_i, ])     # last rows
    
    mat_2 <- cbind(mat_1[, 1:(index-1)],  # first columns
                   NA,                    # pad NA
                   mat_1[,index:d_i])     # last columns
    
    mat <- mat_2
    
  }
  
  
  return(mat)
}

get_network <- function(adj_mat, ID){
  
  # only keep vertices that have edges

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
  num_edges <- 0.5 * (sum(adj_mat_og) - dim(adj_mat_og)[1])     # number of edges (121)
  num_neurons <- dim(adj_mat_fill)[1]                           # number of neurons in total (169)
  num_NA <- sum(is.na(diag(adj_mat_fill)))                      # number of NA neurons       (12)
  num_islands <- length(which(apply(adj_mat_og,2,sum) == 1))    # number of island neurons   (102)
  
  
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
  

  
  return(sum_stats)
  
}


unpack <- function(n, movement, VR, epoch_num, min_edges, IDs, ID2, df_brain_region){
  

  setting_ID <- paste('e', '_m', sep = as.character(epoch_num))          # e1_m
  setting_ID <- paste(setting_ID, 'vr', sep = as.character(movement))    # e1_m2vr
  setting_ID <- paste(setting_ID, '_n', sep = as.character(VR))          # e1_m2vr0_n
  setting_ID <- paste(setting_ID, as.character(n), sep = '')             # e1_m2vr0_n400
  setting_ID <- paste(setting_ID, as.character(min_edges), sep = '_me')  # e1_m2vr0_n400_me1
  
  results_root <- paste("/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/result_simu/",
                        '/',
                        sep = setting_ID)
  
  graphs <- list()
  g_vis <- list()
  edge_vec <- c()
  neuron_vec <- c()
  sum_stats <- list()
  
  for(i in 1:length(IDs)){
  
    # 1) pasting
    file_ID <- paste(ID2[i], setting_ID, sep = '_')                   # Tau1_e1_m2vr0_n400_me1
    
    data_ID <- paste(results_root, '.rda', sep = file_ID)             # ../result_simu/e1_m2vr0_n400_me1/Tau1_e1_m2vr0_n400_me1.rda
    
    # 2) loading
    load(data_ID)
    
    # 3) fill in NA's
    adj_mat_og <- graph_all[['GPP_BIC']]
    missing_neurons <- graph_all[['missing_neurons']] %>% sort()
    adj_mat_fill <- insert_na_symmetric(adj_mat_og, missing_neurons)
    
    # 4) graph and store
    region_border <- max(df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == ID2[i]]) + 0.5
    
    g_i <- adjacency_heatmap(adj_mat_fill, ID2[i], region_border=region_border)
    
    graphs[[i]] <- g_i
    
    # 5) visualization of network
    g_vis[[i]] <- get_network(adj_mat_fill, file_ID)
    
    # 5) summary statistics
    
    sum_stats[[i]] <- get_summary_statistics(adj_mat_fill, adj_mat_og)
    
  }
  
  # save heatmap ---------------------------------------------------
  
  pdf_name <- paste(setting_ID, '.pdf', sep = '_heatmap')
  pdf_name <- paste(results_root, pdf_name, sep = '')
  pdf(pdf_name, width = 10, height = 10)
  
  for(i in 1:length(IDs)){
    print(graphs[[i]])
  }
  
  dev.off()
  
  # save graph visualization --------------------------------------
  
  pdf_name <- paste(setting_ID, '.pdf', sep = '_graph_vis')
  pdf_name <- paste(results_root, pdf_name, sep = '')
  pdf(pdf_name, width = 10, height = 10)
  
  for(i in 1:length(IDs)){
    
    HIP_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == ID2[i]] 
    EHC_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Entorhinal_Cortex" & df_brain_region$ID2 == ID2[i]] 
    
    plot_pts <- sum_stats[[i]][['verts']]
    
    plot_pts_color <- rep('gold', length(plot_pts))
    plot_pts_color[plot_pts %in% EHC_ID] <- 'dodgerblue'
    
    plot(g_vis[[i]], 
         vertex.size = 8,
         vertex.label = V(g_vis[[i]])$name,
         vertex.color = plot_pts_color,
         edge.color = "black",
         edge.arrow.size = 0.5,
         main = ID2[i])
  }
  
  dev.off()
  
  # save summary -------------------------------------------------
  
  txt_name <- paste(setting_ID, '.txt', sep = '_summary')
  txt_name <- paste(results_root, txt_name, sep = '')
  
  sink(txt_name)
  
  for(i in 1:length(IDs)){

    cat('-----------------------------------------\n')
    
    cat('Mouse ID: ')
    cat(ID2[i])
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


load('/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/spike_data/Brain_Region.RData')

n <- 400   
movement <- 2
VR <- 0
epoch_num <- 1
min_edges <- 1

IDs <- c('346', '351', '366', '361', '362', '368')  # mouse ID
ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3') # our name 

IDs <- c('346', '351', '366', '361')  # mouse ID
ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1') # our name 

unpack(n, movement, VR, epoch_num, min_edges, IDs, ID2, df_brain_region)
