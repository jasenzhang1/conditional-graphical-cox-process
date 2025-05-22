

get_accuracy <- function(EstG=NULL, TrueG = NULL){
  
  # -------------------------
  #
  # - EstG (matrix): 40 x 40 binary adjacency matrix
  # - TrueG (matrix): 40 x 40 binary adjacency matrix
  #
  
  # 1) only compare the upper triangular (not including diagonals)
  est = EstG[upper.tri(EstG)]
  true = TrueG[upper.tri(TrueG)]
  accu = sum(est==true)/length(est)
  
  # 2) other statistics
  TP = sum(est==1 & true==1)
  FP = sum(est==1 & true==0)
  FN = sum(est==0 & true==1)
  TN = sum(est==0 & true==0) 
  TPR = TP / sum(true==1)
  TNR = TN / sum(true==0)
  FPR = FP / sum(true==0)
  FNR = FN / sum(true==1)
  PPV = TP / sum(est==1)
  NPV = TN / sum(est==0)
  F1 = 2*TP/(2*TP+FP+FN)
  res = c(accu,TPR,TNR,FPR,FNR, PPV, NPV,F1,TP,FP,FN,TN)
  names(res) = c('Accuracy','TPR','TNR', 'FPR', 'FNR','PPV','NPV','F1',
                 'TP','FP','FN','TN')
  res
}