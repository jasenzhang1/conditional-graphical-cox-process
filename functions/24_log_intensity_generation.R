


prec_mat_massager <- function(prec_mat, manual_thresh = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: if we manually generate a precision matrix, we need to:
  #
  # 1) invert it to be a covariance matrix
  # 2) normalize it to be a correlation matrix
  # 3) invert it again to be a standardized precision matrix
  # 4) make it a partial correlation matrix
  #
  #
  #
  # input: 
  #
  # - prec_mat (p x p matrix)  un-normalized precision matrix, where 0's mean no adjacency
  #
  # 
  # output:
  #
  # - list of precision matrix features
  #   - adj_mat         (p x p matrix)    adjacency matrix with 0's on the diagonal. 1 = dependent, 0 = independent
  #   - prec_mat        (p x p matrix)    precision matrix from the massaged correlation matrix
  #   - cor_mat         (p x p matrix)    correlation matrix 
  #   - partial_cor_mat (p x p matrix)    partial correlation matrix
  #   - threshold_p    (number)           smallest non-zero off-diagonal value. Anything less than this threshold will be assumed to be 0
  #
  # ----------------------------------------------------------------------------
  

  
  # 1) start with prec_mat, make it cov_mat, then normalize it to be cor_mat
  cov_mat <- solve_sym(prec_mat)
  cor_mat <- cov_mat %>% cov2cor() %>% sym()
  
  # 2) prec_mat2 = precision matrix from correlation matrix
  prec_mat2 <- solve_sym(cor_mat)
  
  
  # 3) partial correlation 
  D <- diag(1 / sqrt(diag(prec_mat))) 
  partial_cor_mat <- -D %*% prec_mat %*% D  
  diag(partial_cor_mat) <- 1  
  
  # 4) define the threshold as the min value that should not be nonzero
  threshold_p <- min(abs(prec_mat2[prec_mat != 0]))  
  
  # 5) adjacency_matrix takes on a value of 1 if original prec_mat is nonzero
  #    or, we define 1 with an indicator function
  
  if(!is.null(manual_thresh)){
    adj_mat <- matrix(0, p, p)
    adj_mat[abs(prec_mat) >= manual_thresh] <- 1
    diag(adj_mat) <- 0       
  } else{
    adj_mat <- matrix(0, p, p)
    adj_mat[prec_mat != 0] <- 1
    diag(adj_mat) <- 0      
  }
  

  
  return(list(adj_mat = adj_mat, 
              prec_mat_og = prec_mat,
              prec_mat = prec_mat2, 
              cor_mat = cor_mat, 
              cov_mat = cov_mat,
              partial_cor_mat = partial_cor_mat, 
              simu_mat = prec_mat2,
              threshold_p = threshold_p))
}



assemble_block_matrix_24 <- function(precision_operators){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: assemble block matrix
  #
  #
  # input:
  # 
  # - precision_operators (p x p x m x m)
  #
  # 
  # output:
  # 
  # - precision_matrix (pm x pm)
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(precision_operators)[1]
  m <- dim(precision_operators)[3]
  
  precision_matrix <- matrix(0, p*m, p*m)
  
  for(i in 1:p){
    i_start <- (i-1)*m + 1
    i_end <- i*m 
    for(j in i:p){
      j_start <- (j-1)*m + 1
      j_end <- j*m
      
      precision_matrix[i_start:i_end, j_start:j_end] <- precision_operators[i, j, ,]
      precision_matrix[j_start:j_end, i_start:i_end] <- t(precision_operators[i, j, ,])
    }
  }
  
  return(precision_matrix)
  
}


