
kernel_RBF <- function(s_yd, df_cont_cov, gamma){
  
  # s_yd    (vector of numbers)     replicate indices for a specific discrete variable strata
  # df_cov  (dataframe)             matrix storing vectors of continuous covariates, where the i-th row represents the i-th replicate
  # gamma   (number)                RBF hyperparameter: exp(- gamma * ||y_1 - y_2||^2)
  
  df_cov_2 <- df_cont_cov[s_yd, , drop = FALSE]
  n <- length(s_yd)
  
  kernel_mat <- matrix(0, n, n)
  
  # lower diag
  for(i in 2:n){
    for(j in 1:(i-1)){
      
      vec_i <- df_cov_2[i,] %>% as.numeric()
      vec_j <- df_cov_2[j,] %>% as.numeric()
      
      kernel_mat[i, j] <- exp(- gamma * sum( (vec_i - vec_j)^2 ) )
      
    }
  }
  
  # fill in upper diag
  kernel_mat <- kernel_mat + t(kernel_mat)
  
  # fill in diag
  diag(kernel_mat) <- 1
  
  return(kernel_mat)

}

## kernel density destimation --------------------------------------



