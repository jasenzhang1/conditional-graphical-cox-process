

evaluate_graph_recovery <- function(estimated_edges, true_adjacency) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Compute graph recovery metrics
  #
  #
  # Input: 
  #
  # - estimated_edges (list of edge pairs)
  # - true_adjacency  (p x p matrix)
  # 
  #
  # Output: 
  #
  # - list with sensitivity, specificity, precision, F1, MCC
  #
  #
  # ----------------------------------------------------------------------------
  
  p <- nrow(true_adjacency)
  
  # Convert estimated edges to adjacency matrix
  est_adj <- matrix(0, p, p)
  if (length(estimated_edges) > 0) {
    for (edge in estimated_edges) {
      i <- edge[1]
      j <- edge[2]
      est_adj[i, j] <- est_adj[j, i] <- 1
    }
  }
  
  # Extract upper triangle indices (undirected graph)
  true_edges <- which(upper.tri(true_adjacency) & true_adjacency == 1)
  est_edges <- which(upper.tri(est_adj) & est_adj == 1)
  all_possible_edges <- which(upper.tri(matrix(1, p, p)))
  
  # Confusion matrix elements
  TP <- length(intersect(true_edges, est_edges))    # True Positives
  FP <- length(setdiff(est_edges, true_edges))      # False Positives  
  TN <- length(setdiff(setdiff(all_possible_edges, true_edges), est_edges)) # True Negatives
  FN <- length(setdiff(true_edges, est_edges))      # False Negatives
  
  # Performance metrics
  sensitivity <- if (TP + FN > 0) TP / (TP + FN) else 0        # Recall
  specificity <- if (TN + FP > 0) TN / (TN + FP) else 0        # True Negative Rate
  precision <- if (TP + FP > 0) TP / (TP + FP) else 0         # Precision
  f1_score <- if (precision + sensitivity > 0) {              # F1-Score
    2 * precision * sensitivity / (precision + sensitivity)
  } else 0
  
  # Matthews Correlation Coefficient
  mcc_denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (mcc_denom > 0) (TP * TN - FP * FN) / mcc_denom else 0
  
  return(list(
    sensitivity = sensitivity,    # Sen = |E_hat ∩ E*| / |E*|
    specificity = specificity,    # Spec = |E_hat^c ∩ (E*)^c| / |(E*)^c|
    precision = precision,        # Prec = |E_hat ∩ E*| / |E_hat|
    f1_score = f1_score,         # F1 = 2 * Prec * Rec / (Prec + Rec)
    mcc = mcc,                   # Matthews Correlation Coefficient
    confusion_matrix = list(TP = TP, FP = FP, TN = TN, FN = FN)
  ))
}

compute_operator_error <- function(estimated_precision, true_precision) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: Compute operator norm and Frobenius norm errors
  # 
  #
  # Input: 
  #
  # - estimated_precision  (pm x pm matrix)
  # - true_precision       (pm x pm matrix)
  #
  #
  # Output: 
  #
  # - list with operator and Frobenius norm errors
  #
  # ---------------------------------------------------------------------------- 
  
  error_matrix <- estimated_precision - true_precision
  
  # Operator norm: largest singular value
  operator_norm_error <- norm(error_matrix, type = "2")
  
  # Frobenius norm: sqrt(sum of squared elements)
  frobenius_norm_error <- norm(error_matrix, type = "F")
  
  return(list(
    operator_norm = operator_norm_error,    # ||P_hat - P*||_op
    frobenius_norm = frobenius_norm_error   # ||P_hat - P*||_F
  ))
}