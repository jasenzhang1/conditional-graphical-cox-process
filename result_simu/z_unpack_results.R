# I have a bunch of .rda results

# now I want to visualize them in a graph

library(ggplot2)

adjacency_heatmap <- function(adj_mat, ID){
  
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

unpack <- function(n, movement, VR, epoch_num, min_edges, IDs, ID2){
  

  setting_ID <- paste('e', '_m', sep = as.character(epoch_num))
  setting_ID <- paste(setting_ID, 'vr', sep = as.character(movement))
  setting_ID <- paste(setting_ID, '_n', sep = as.character(VR))
  setting_ID <- paste(setting_ID, as.character(n), sep = '')
  setting_ID <- paste(setting_ID, as.character(min_edges), sep = '_me')
  
  results_root <- paste("/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/result_simu/",
                        '/',
                        sep = setting_ID)
  
  graphs <- list()
  edge_vec <- c()
  neuron_vec <- c()
  for(i in 1:length(IDs)){
  
    # 1) pasting
    file_ID <- paste(ID2[i], setting_ID, sep = '_')
    
    data_ID <- paste(results_root, '.rda', sep = file_ID)   
    
    # 2) loading
    load(data_ID)
    
    # 3) fill in NA's
    adj_mat_og <- graph_all[['GPP_BIC']]
    missing_neurons <- graph_all[['missing_neurons']] %>% sort()
    adj_mat_fill <- insert_na_symmetric(adj_mat_og, missing_neurons)
    
    # 4) graph and save
    g_i <- adjacency_heatmap(adj_mat_fill, ID2[i])
    
    graphs[[i]] <- g_i
    
    # 5) summary statistics
    
    num_edges <- 0.5 * (sum(adj_mat_og) - dim(adj_mat_og)[1])
    num_neurons <- dim(adj_mat_og)[1]
    edge_vec <- c(edge_vec, num_edges)
    neuron_vec <- c(neuron_vec, num_neurons)
    
  }
  
  # graph ---------------------------------------------------
  
  pdf_name <- paste(setting_ID, '.pdf', sep = '')
  pdf(pdf_name, width = 10, height = 10)
  
  for(k in 1:length(graphs)){
    print(graphs[[k]])
  }
  
  dev.off()
  
  # summary -------------------------------------------------
  
  summary_df <- data.frame(ID2, neuron_vec, edge_vec)
  colnames(summary_df) <- c('Mouse', 'Num Neurons', 'Num Edges')
  
  csv_name <- paste(setting_ID, '.csv', sep = '')
  write.csv(summary_df, csv_name)
  return(NULL)
}

n <- 400   
movement <- 2
VR <- 0
epoch_num <- 1
min_edges <- 1

IDs <- c('346', '351', '366', '361', '362', '368')  # mouse ID
ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3') # our name 

IDs <- c('346')  # mouse ID
ID2 <- c('Tau1') # our name 

unpack(n, movement, VR, epoch_num, min_edges, IDs, ID2)
