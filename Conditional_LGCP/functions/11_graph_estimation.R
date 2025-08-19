hilbert_schmidt_norm <- function(A) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: Compute Hilbert-Schmidt norm, which is just frobenius norm 
  #
  # - ||A||_{HS} = sqrt(sum of squares of elements)
  # 
  # 
  # Input: 
  #
  # - A       (m x m matrix)
  #
  # 
  # Output: 
  #
  # - HS norm (scalar)
  #
  # ------------------------------------------------------------------------
  
  return(sqrt(sum(A^2)))  # sqrt(sum of all squared elements)
}

hilbert_schmidt_norm_pm <- function(A, p, m) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: Compute Hilbert-Schmidt norm of all p^2 mxm block matrices
  #
  # - ||A||_{HS} = sqrt(sum of squares of elements)
  # 
  # 
  # Input: 
  #
  # - A       (pm x pm matrix)
  #
  # 
  # Output: 
  #
  # - HS norm (scalar)
  #
  # ------------------------------------------------------------------------
  

  if (!is.matrix(A) || nrow(A) != p*m || ncol(A) != p*m) {
    stop("matrix input must be a (p*m) x (p*m) matrix")
  }
  
  norms <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      # Row indices for block (i, j)
      row_idx <- ((i - 1) * m + 1):(i * m)
      # Column indices for block (j, j)
      col_idx <- ((j - 1) * m + 1):(j * m)
      block <- A[row_idx, col_idx]
      norms[i, j] <- sqrt(sum(block^2))  # HS (Frobenius) norm
    }
  }
  
  return(norms)

}

list_to_adj_mat <- function(edge_list, p){
  
  # I have a list of edge strengths
  # Create a weighted adjacency matrix 
  #
  
  
  adj_mat <- matrix(0, p, p)
  
  for(k in names(edge_list)){
    
    i_j <- as.integer(strsplit(k, "_")[[1]])
    
    i <- i_j[1]
    j <- i_j[2]
    
    adj_mat[i, j] <- edge_list[[k]]
    adj_mat[j, i] <- edge_list[[k]]
    
  }
  
  return(adj_mat)
  
}

estimate_graph <- function(P_conditional, C_conditional, V_conditional, threshold, p, discarded_neurons = NULL) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: get edge estiamtes for our graph
  # 
  # 
  # Input: 
  #
  # - P_conditional   (list of length p^2, each element m x m)
  # - threshold       (scalar)
  # - p               (scalar)
  # - discarded_neurons (list of integers)
  #
  # 
  # Output: 
  # - list with edges and edge_strengths
  #
  # ------------------------------------------------------------------------
  
  edges <- list()
  edge_strengths <- list()
  
  
  adj_mat <- matrix(0, p, p)
  w_mat <- matrix(0, p, p)
  
  # Only consider upper triangle: i < j (undirected graph)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key_ij <- paste(i, j, sep="_")
      P_ij <- P_conditional[[key_ij]]  # m x m matrix
      
      
      
      # Compute Hilbert-Schmidt norm: m x m -> scalar
      hs_norm <- hilbert_schmidt_norm(P_ij)
      edge_strengths[[key_ij]] <- hs_norm
      w_mat[i,j] <- hs_norm
      w_mat[j,i] <- hs_norm
      
      # Apply threshold: scalar > scalar -> boolean
      if (hs_norm > threshold) {
        edges <- append(edges, list(c(i, j)))
        
        adj_mat[i,j] <- 1
        adj_mat[j,i] <- 1
      }
    }
  }
  
  m <- dim(P_ij)[1]
  
  return(list(edges = edges,                       # list of [[i_j]] entries that satisfy the threshold
              edge_strengths = edge_strengths,     # list of [[i_j]] weights
              w_mat = w_mat,
              adj_mat = adj_mat,
              P_conditional = P_conditional,       # original precision operator
              C_conditional = C_conditional,
              V_conditional = V_conditional,
              P_mat = assemble_block_matrix_v2(P_conditional, p, m),
              C_mat = assemble_block_matrix_v2(C_conditional, p, m),
              V_mat = assemble_block_matrix_v2(V_conditional, p, m),
              discarded_neurons = discarded_neurons))
}

select_threshold_by_stability <- function(P_conditional, p) {
  
  # ------------------------------------------------------------------------
  #
  # GOAL: how to choose edge selection threshold?
  # 
  # 
  # Input: 
  #
  # - P_conditional    (p^2 length list)
  # - p                (scalar)
  #
  # Output: 
  #
  # - threshold        (scalar)
  #
  # ------------------------------------------------------------------------
  
  # Collect all edge strengths
  edge_strengths <- c()
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      key <- paste(i, j, sep="_")
      P_ij <- P_conditional[[key]]
      edge_strengths <- c(edge_strengths, hilbert_schmidt_norm(P_ij))
    }
  }
  
  # Use median as heuristic threshold
  return(if (length(edge_strengths) > 0) median(edge_strengths) else 0.1)
}



library(pROC)

roc_with_threshold <- function(w_mat, adj_mat, g_title = NULL) {
  
  #
  #
  # GOAL: identify the ideal threshold for selecting an edge using ROC curves
  #
  #
  # input:
  # 
  # - w_mat (p x p matrix)
  # - adj_mat (p x p matrix)
  #
  #
  #
  #
  
  # 1) rearrange w_mat an adj_mat into vectors of the upper diagonal entries  
  
  
  scores <- w_mat[upper.tri(w_mat, diag = FALSE)]
  labels <- adj_mat[upper.tri(adj_mat, diag = FALSE)]
  
  # labels: vector of 0/1 or factor (true outcome)
  # scores: numeric vector of predicted probabilities or scores
  
  # Ensure labels are binary factors
  labels <- as.factor(labels)
  if (length(levels(labels)) != 2) stop("Labels must have two classes.")
  
  # Compute ROC
  roc_obj <- roc(labels, scores, quiet = TRUE)
  auc_value <- auc(roc_obj)   
  
  # Get thresholds, sensitivities, specificities
  coords_df <- coords(roc_obj, x = "all", ret = c("threshold", "sensitivity", "specificity"))
  
  # Compute Youden's J statistic
  youden <- coords_df$sensitivity + coords_df$specificity - 1
  best_idx <- which.max(youden)
  
  ideal_threshold <- coords_df$threshold[best_idx]
  ideal_sens <- coords_df$sensitivity[best_idx]
  ideal_spec <- coords_df$specificity[best_idx]
  
  # Create data frame for ggplot
  roc_df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )
  
  roc_df <- roc_df[order(roc_df$FPR, roc_df$TPR), ]
  
  # ROC plot
  p <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_step(direction = "vh", color = "blue", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
    annotate("point", x = 1 - ideal_spec, y = ideal_sens, color = "red", size = 3) +
    annotate("text", x = 1 - ideal_spec, y = ideal_sens, 
             label = paste0("Threshold=", round(ideal_threshold, 3)),
             hjust = -0.1, vjust = -0.5, color = "red") +
    
    annotate("text", x = 0.6, y = 0.2,            # position for AUC label
             label = paste0("AUC = ", round(auc_value, 3)),
             color = "darkgreen", size = 5) +    
    theme_minimal() + 
    ggtitle(g_title)
  
  list(
    threshold = ideal_threshold,
    sensitivity = ideal_sens,
    specificity = ideal_spec,
    auc = auc(roc_obj),
    plot = p
  )
}
