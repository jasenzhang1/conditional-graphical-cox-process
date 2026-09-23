
block_mult <- function(x, A_block_vals, n_s) {
  
  # 
  # performs Ax, where A is denoted by A_block_vals and n_s
  #
  # x (vector)
  # A_block_vals (matrix) 
  # n_i (sizes of each block)
  #
  # 
  
  m <- dim(A_block_vals)[1]
  stopifnot(length(x) == sum(n_s))
  result <- numeric(sum(n_s))
  
  row_start <- 1
  for (i in 1:m) {
    n_i <- n_s[i]
    col_start <- 1
    for (j in 1:m) {
      n_j <- n_s[j]
      
      a_ij <- A_block_vals[i, j]
      
      x_sub <- x[col_start:(col_start + n_j - 1)]
      block_contrib <- sum(x_sub) * a_ij  # since all entries are a_ij
      
      result[row_start:(row_start + n_i - 1)] <- 
        result[row_start:(row_start + n_i - 1)] + block_contrib
      
      col_start <- col_start + n_j
    }
    row_start <- row_start + n_i
  }
  
  return(result)
}

conjugate_gradient_solve <- function(A_block_vals, n_s, b){
  
  # 
  # Goal: solve Ax = b via conjugate gradient for a specific block matrix 
  #
  # 
  # input:
  # 
  # - A_block_vals (matrix)  each entry is a_ij, the value of the ijth block
  # - n_s          (vector)  a vector of numbers telling you the partition sizes of the matrix
  # - b            (vector)  a vector to solve for. Should be length of sum(n_s)
  #
  # 
  # output:
  # 
  # - x_sol        (vector) the solution to Ax = b
  # 
  
  
  
  A_func <- function(x) block_mult(x, A_block_vals, n_s)
  
  # Solve A x = b
  res <- pcg(A_func, b, tol = 1e-6, maxiter = 1000)
  
  x_sol <- res$x
  
  return(x_sol)
  
  A_op <- as(Matrix::Matrix(0, n_s, n_s, sparse = TRUE), "matrix")
  class(A_op) <- c("linop", class(A_op))
  attr(A_op, "mult") <- function(x) block_mult(x)
  
  # Use Krylov solver
  x_sol <- Matrix::Krylov(A_op, b, m = 100, method = "cg")  # m = max iterations  
  
}



A_block_vals <- matrix(runif(9), 3, 3)  # a_ij values
A_block_vals <- 0.5 * (A_block_vals + t(A_block_vals))

n_s <- sample(3:6, 3, replace = TRUE)
b <- runif(sum(n_s))

x_sol <- conjugate_gradient_solve(A_block_vals, n_s, b)