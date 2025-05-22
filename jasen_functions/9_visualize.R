library(tidyverse)

adjacency_heatmap <- function(adj_mat){
  
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
    theme_minimal()  
  
  g3
}

gpp_heatmap <- function(m1, grpind, graph_true=NULL){
  
  # m1 = un-truncated matrix
  #
  #
  #
  
  m_size <- dim(m1)[1]
  
  m1_melt <- reshape2::melt(m1)
  
  m1_melt$Var1 <- m_size + 1 - m1_melt$Var1
  
  m1_melt$nonzero <- as.numeric(m1_melt$value != 0)
  
  
  
  # regular heatmap
  g2 <- ggplot(m1_melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Column", y = "Row", fill = "Value") +
    theme_minimal()
  

  
  # add square blocks
  
  border_id <- c(0, grpind[,2]) + 0.5
  
  g2 <- g2 + geom_vline(xintercept = border_id, alpha = 0.3) + 
    geom_hline(yintercept = border_id, alpha = 0.3)
  
  k <- 1
  graphs <- list()
  graphs[[k]] <- g2
  k <- k + 1  
  
  # nonzero heatmap
  g2_nonzero <- ggplot(m1_melt, aes(x = Var2, y = Var1, fill = nonzero)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Column", y = "Row", fill = "Value") +
    geom_vline(xintercept = border_id, alpha = 0.3) + 
    geom_hline(yintercept = border_id, alpha = 0.3) +
    theme_minimal()
  
  graphs[[k]] <- g2_nonzero
  k <- k + 1    
  
  # final adjacency heatmap
  m_final <- as.matrix((get_groupNorm(m1, grpind)!=0)+0)
  
  m2_melt <- reshape2::melt(m_final)
  m2_melt$Var1 <- dim(m_final)[1] + 1 - m2_melt$Var1
  
  border_id3 <- 0:40 + 0.5
  
  g3 <- ggplot(m2_melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Column", y = "Row", fill = "Value") +
    geom_vline(xintercept = border_id3, alpha = 0.3) + 
    geom_hline(yintercept = border_id3, alpha = 0.3) +    
    theme_minimal() 
  
  graphs[[k]] <- g3
  k <- k+1
  
  # true adjacency
  
  if(!is.null(graph_true)){
    m3_melt <- reshape2::melt(graph_true)
    m3_melt$Var1 <- dim(m_final)[1] + 1 - m3_melt$Var1
    
    g_true <- ggplot(m3_melt, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(x = "Column", y = "Row", fill = "Value") +
      geom_vline(xintercept = border_id3, alpha = 0.3) + 
      geom_hline(yintercept = border_id3, alpha = 0.3) +
      theme_minimal()  
    
    graphs[[k]] <- g_true
  }
  

  
  return(graphs)
}

# graphs_GMpp <- gpp_heatmap(res_all[["GPP_BIC"]], grpind, graph_true)
# graphs_COUNT <- gpp_heatmap(res_all[["count_BIC"]], grpind, graph_true)

get_event_counts <- function(data_all, patient_sel, feature_sel){
  N = length(patient_sel)
  p = length(feature_sel)
  
  # 1) ensure the dataset is proper + sorted
  data_all = data_all[is.element(feature_id, feature_sel),]
  data_all = data_all[is.element(subject_num, patient_sel),]
  
  setorder(data_all, subject_num, feature_id, time)
  
  # 2) create a table of counts
  data.cov = data_all[,.(count = .N),by=c("subject_num","feature_id")]
  
  # 2a) for each subject_num and feature_id, find its respective index in patient_sel or feature_sel
  #     but if patient_sel and feature_sel are just 1:n, 1:p, this is trivially just subject_num and feature_id again
  data.cov[, i:= match(subject_num, patient_sel) ]
  data.cov[, j:= match(feature_id , feature_sel)]  
  
  g <- ggplot(data.cov, aes(x = j, y = i, fill = count)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(x = "Feature", y = "Subject", fill = "Count") +
    theme_minimal()
}