


roc_for_thresholded_w_mat <- function(w_mat, adj_mat) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: for a w_mat that's thresholded already, report ROC metrics
  #
  #       - AUC and roc_df are NA since the ranking info is lost after thresholding
  #
  #
  # input:
  #
  # - w_mat            (p x p matrix)         matrix of 0's and nonzero values
  # - adj_mat          (p x p matrix)         matrix of 0's and 1's with 0's on the diagonal
  #
  #
  # output:
  #
  # - output           (list)
  #   - threshold      (NA)
  #   - accuracy       (number)
  #   - f1_score       (number)
  #   - sensitivity    (number)
  #   - specificity    (number)
  #   - ppv            (number)
  #   - npv            (number)
  #   - auc            (NA)
  #   - roc_df         (NA)
  #   - counts         (list)                 TP, FP, TN, FN
  #
  # ----------------------------------------------------------------------------
  
  
  # 1) Extract upper triangle to avoid double-counting or diagonal bias
  scores <- w_mat[upper.tri(w_mat, diag = FALSE)]
  labels <- adj_mat[upper.tri(adj_mat, diag = FALSE)]
  
  # 2) Convert to binary predictions (Already zeroed out -> Non-zero is 1)
  predicted <- ifelse(scores != 0, 1, 0)
  true      <- as.numeric(as.character(labels))
  
  # 3) Tally Confusion Matrix components
  TP <- sum(predicted == 1 & true == 1)
  FP <- sum(predicted == 1 & true == 0)
  TN <- sum(predicted == 0 & true == 0)
  FN <- sum(predicted == 0 & true == 1)
  
  # 4) Calculate Metrics with NA handling 
  
  # Accuracy: (TP + TN) / Total
  accuracy <- (TP + TN) / length(true)
  
  # Sensitivity (Recall / TPR): TP / (TP + FN)
  # NA if truth has no positives (no edges)
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  
  # Specificity (TNR): TN / (TN + FP)
  # NA if truth has no negatives (fully connected)
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  
  # PPV (Positive Predictive Value / Precision): TP / (TP + FP)
  # NA if we predicted zero edges
  ppv <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  
  # NPV (Negative Predictive Value): TN / (TN + FN)
  # NA if we predicted everything is an edge
  npv <- if ((TN + FN) > 0) TN / (TN + FN) else NA
  
  # F1 Score: 2 * (Precision * Recall) / (Precision + Recall)
  f1 <- if (!is.na(ppv) && !is.na(sensitivity) && (ppv + sensitivity) > 0) {
    2 * (ppv * sensitivity) / (ppv + sensitivity)
  } else {
    NA
  }
  

  
  # Note: AUC is usually not sensible for a single fixed threshold 
  # but some define it as the area under the single-point 'curve'.
  # Here we return NA for AUC as the 'ranking' info is lost.
  
  return(list(
    threshold   = NA,
    accuracy    = accuracy,
    f1_score    = f1,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv         = ppv,
    npv         = npv,
    auc         = NA,
    roc_df      = NA,
    counts      = list(TP = TP, FP = FP, TN = TN, FN = FN)
  ))
}

roc_for_raw_w_mat <- function(w_mat, adj_mat) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: identify the ideal threshold for selecting an edge using ROC curves (Youden's J)
  #
  #
  # input:
  #
  # - w_mat            (p x p matrix)         edge weights
  # - adj_mat          (p x p matrix)         true adjacency matrix with 0's on the diagonal
  #
  #
  # output:
  #
  # - output           (list)                 same structure as roc_for_thresholded_w_mat (all NA if adj_mat has only one class)
  #   - threshold      (number)               ideal threshold
  #   - accuracy       (number)
  #   - f1_score       (number)
  #   - sensitivity    (number)
  #   - specificity    (number)
  #   - ppv            (number)
  #   - npv            (number)
  #   - auc            (number)
  #   - roc_df         (data.frame)           FPR and TPR columns for plotting
  #   - counts         (list)                 TP, FP, TN, FN
  #
  # ----------------------------------------------------------------------------
  
  # 1) rearrange w_mat an adj_mat into vectors of the upper diagonal entries  
  
  
  scores <- w_mat[upper.tri(w_mat, diag = FALSE)]
  labels <- adj_mat[upper.tri(adj_mat, diag = FALSE)]
  
  # labels: vector of 0/1 or factor (true outcome)
  # scores: numeric vector of predicted probabilities or scores
  
  # Ensure labels are binary factors
  labels_factor <- as.factor(labels)
  if (length(levels(labels_factor)) != 2){
    return(list(
      threshold = NA, sensitivity = NA, specificity = NA,
      ppv = NA, npv = NA, f1_score = NA,
      auc = NA, accuracy = NA, roc_df = NA, counts = NA
    ))
  }
  
  # 2) Compute ROC and AUC
  roc_obj   <- roc(labels_factor, scores, quiet = TRUE)
  auc_value <- as.numeric(auc(roc_obj))
  
  # 3) Identify ideal threshold using Youden's J
  coords_df <- coords(roc_obj, x = "all", ret = c("threshold", "sensitivity", "specificity"))
  youden    <- coords_df$sensitivity + coords_df$specificity - 1
  best_idx  <- which.max(youden)
  
  ideal_threshold <- coords_df$threshold[best_idx]
  ideal_sens      <- coords_df$sensitivity[best_idx]
  ideal_spec      <- coords_df$specificity[best_idx]
  
  # 4) Calculate binary metrics based on the ideal threshold
  predicted <- ifelse(scores >= ideal_threshold, 1, 0)
  true      <- as.numeric(as.character(labels_factor))
  
  TP <- sum(predicted == 1 & true == 1)
  FP <- sum(predicted == 1 & true == 0)
  TN <- sum(predicted == 0 & true == 0)
  FN <- sum(predicted == 0 & true == 1)
  
  # Accuracy
  accuracy <- mean(predicted == true)
  
  # PPV (Precision)
  ppv <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  
  # NPV
  npv <- if ((TN + FN) > 0) TN / (TN + FN) else NA
  
  # F1 Score
  f1 <- if (!is.na(ppv) && !is.na(ideal_sens) && (ppv + ideal_sens) > 0) {
    2 * (ppv * ideal_sens) / (ppv + ideal_sens)
  } else {
    NA
  }
  
  # 5) Create data frame for ggplot
  roc_df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )
  roc_df <- roc_df[order(roc_df$FPR, roc_df$TPR), ]
  
  # Return identical structure to roc_for_thresholded_w_mat
  return(list(
    threshold   = ideal_threshold,
    accuracy    = accuracy,
    f1_score    = f1,
    sensitivity = ideal_sens,
    specificity = ideal_spec,
    ppv         = ppv,
    npv         = npv,
    auc         = auc_value,
    roc_df      = roc_df,
    counts      = list(TP = TP, FP = FP, TN = TN, FN = FN)
  ))
}

