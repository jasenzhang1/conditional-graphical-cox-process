source('functions/22_generate_random_variables.R')
source('functions/00c_block_matrix_arrange.R')

# here are all the functions needed to generate X(t) from a finite basis

trig_basis <- function(d){
  
  # ----------------------------------------------------------------------------
  #
  # obtain the basis of trig functions up to d
  #
  # (1)  1
  # (2)  sin(pi t)
  # (3)  cos(pi t)
  # (4)  sin(2 pi t)
  # etc...
  # 
  #
  # inputs:
  #
  # - d  (integer)    number of components
  #
  #
  #
  # output:
  #
  # - basis_list   (list)   list of the first d basis functions
  #
  # ----------------------------------------------------------------------------
  
  basis_list <- list(function(t) rep(1, length(t)))
  if (d == 1) return(basis_list)
  
  # Use frequencies 1:(d-1) but actual argument is 2*pi*k*t
  trig_funcs <- lapply(1:(d-1), function(k) {
    function(t) sqrt(2) * sin(2 * k * pi * t)
  })
  
  basis_list <- c(basis_list, trig_funcs)
}


trig_basis_realization <- function(basis_list, time_grid){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: for all d basis functions, compute their realizations at timepoints in [tmin, tmax]
  #
  # inputs:
  #
  # - basis_list (list)             output from trig_basis
  # - time_grid  (m-dim vector)
  #
  # output:
  #
  # - basis_mat   (m x d matrix) matrix of all d discretized realizations of the eigenfunctions
  #
  #
  # ----------------------------------------------------------------------------
  

  
  # Evaluate each basis function at all timepoints
  B_mat <- sapply(basis_list, function(f) f(time_grid))  
  
  return(B_mat)
}


trig_basis_cov_mat <- function(d, p, y_c_k, adj_type, adj_params){
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: define the matrix to generate b's 
  #
  # inputs:
  #
  # - d                 (integer)
  # - p                 (integer)
  # - y_c_k             (q_c-dim vector)
  # - adj_type     (string)
  #
  #
  # outputs:
  #
  # cov_mat   (pd x pd matrix)
  #
  # ----------------------------------------------------------------------------
  
  if(! adj_type %in% c('block_banded_v2', 'block_banded_c2')){
    stop('Error 21b: adj_type not available')
  }
  if(! length(y_c_k) == 1){
    stop('Error 21b: y_c_k should be a scalar')
  }
  
  if(adj_type %in% c('block_banded_v2', 'block_banded_c2')){
    
    # calculating J_2_const
    if(adj_type == 'block_banded_v2'){
      y_c_min <- adj_params[1]
      y_c_max <- adj_params[2]
      c_min <- adj_params[3]
      c_max <- adj_params[4]
      
      J_2_const <- c_min + (c_max - c_min) * (y_c_k - y_c_min) / (y_c_max - y_c_min)
    } else{
      J_2_const <- adj_params[3]
    }
    
    J_2 <- J_2_const * (-1)^(1 + 1:d)
    off_block <- diag(J_2)
    
    on_block <- diag(d)
    
    full_pd_mat <- matrix(0, nrow = p*d, ncol = p*d)
    
    for(i in 1:p){
      for(j in 1:p){
        # Compute index ranges for block (i,j)
        row_idx <- ((i - 1) * d + 1):(i * d)
        col_idx <- ((j - 1) * d + 1):(j * d)
        

        if(i == j){
          full_pd_mat[row_idx, col_idx] <- on_block
        }
        
        if(abs(i-j) == 1){
          full_pd_mat[row_idx, col_idx] <- off_block
        }
      }
    }
    return(full_pd_mat)
  }
}

trig_basis_log_intensity <- function(cov_mat_list, basis_list, beta_0, betas, time_grid){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: draw n log-intensities from the trig-basis data generation process
  #
  #
  # input:
  #
  # - cov_mat_list  (n-dim list of pd x pd matrix)  covariance matrix for each subject
  # - basis_list    (list of d eigenfunctions)   
  # - beta_0        (number)                     constant added to the random function
  # - betas         (d-dim vector)               variance factor for each eigenfunction
  # - time_grid     (m-dim vector)               time discretization
  #
  #
  # output:
  # 
  # - log_intensities  (n-dim list)  each item is a (p x m) matrix of log-intensities for all p processes
  #
  # 
  # ----------------------------------------------------------------------------
  
  # dimensions
  n <- length(cov_mat_list)
  d <- length(basis_list)
  m <- length(time_grid)
  
  
  # precompute basis evaluated at time grid
  Phi <- sapply(basis_list, function(f) f(time_grid))  # m x d matrix
  
  
  # output container
  log_intensities <- vector("list", n)
  
  for (s in 1:n) {
    
    cov_mat <- cov_mat_list[[s]]
    pd <- nrow(cov_mat)
    p <- pd / d
    
    if (abs(p - round(p)) > 1e-10)
      stop("cov_mat dimensions are not compatible with number of basis functions.")
    
    # Cholesky for sampling from N(0, cov_mat)
    L <- chol(cov_mat)
    
    # draw b ~ N(0, cov_mat)
    z <- rnorm(pd)
    b <- as.vector(L %*% z)  # pd x 1
    
    # reshape b into p x d matrix of coefficients
    B_mat <- matrix(b, nrow = p, ncol = d, byrow = TRUE)
    
    # compute log-intensity for each process i
    log_int_mat <- matrix(0, nrow = p, ncol = m)
    for (i in 1:p) {
      # linear combination: 1 + sum_k beta_k * b_ik * phi_k(t)
      log_int_mat[i, ] <- beta_0 + colSums(matrix(betas * B_mat[i, ], nrow = d, ncol = m, byrow = FALSE) * t(Phi))
    }
    
    log_intensities[[s]] <- log_int_mat
  }
  
  return(log_intensities)
  
}

# truths

trig_basis_cross_covariance_truth <- function(basis_list, cov_mat, time_grid){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: calculate the ground truth cross-covariance (pm x pm) for any time discretization setting
  #
  # 
  # input:
  # 
  # - basis_list  (d-dim list)             output from trig_basis
  # - cov_mat     (pd x pd matrix)         covariance matrix
  # - time_grid   (m-dim vector)           time discretization
  #
  #
  # output:
  #
  # - C_ij_st      (pm x pm matrix)            discretized cross covariance between processes i and j
  #
  #
  # ----------------------------------------------------------------------------
  
  m <- length(time_grid)
  phi <- trig_basis_realization(basis_list, time_grid) # (m x d matrix)
  d <- dim(phi)[2]
  p <- dim(cov_mat)[1] / d
  
  blocks <- list()
  
  for(i in 1:p){
    for(j in i:p){
      key <- paste0(i, '_', j)
      blocks[[key]] <- phi %*% extract_block_structure_ij(cov_mat, d, i, j) %*% t(phi)
    }
  }
  
  C_mat <- assemble_block_matrix_v2(blocks, p, m)
  
  return(C_mat)
  
}

trig_basis_gram_matrix <- function(basis_list, t_min, t_max){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: obtain the (d x d) gram matrix of inner products of all basis functions
  #
  #
  # input: 
  #
  # - basis_list   (d-dim list)       each entry is a basis functions when we set up the log-intensities
  # - t_min        (number)           lower bound of integration
  # - t_max        (number)           upper bound of integration
  #
  #
  # output:
  #
  # - G         (d x d matrix)   Gram matrix, which contains the inner product between each pair of basis functions
  #
  # ----------------------------------------------------------------------------
  
  
  
  d <- length(basis_list)
  G <- matrix(0, nrow=d, ncol=d)
  
  # assuming t in [0,1] (change limits if needed)
  for(k in 1:d){
    for(l in k:d){  # symmetric
      f <- function(t) basis_list[[k]](t) * basis_list[[l]](t)
      G[k,l] <- integrate(f, lower=t_min, upper=t_max)$value
      G[l,k] <- G[k,l]
    }
  }
  
  # check if G is identity
  tol <- 1e-10
  
  diag_is_one <- max(abs(diag(G) - 1)) < tol

  off_diag_is_zero <- all(abs(G[!diag(nrow(G))  ]) < tol)
  
  if(diag_is_one & off_diag_is_zero){
    G <- diag(d)
  } 
  
  
  return(G)  # your d x d Gram matrix  
}

trig_basis_eigendecomposition <- function(G, cov_mat, basis_list, time_grid, betas){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: find ground truths in the trig_basis framework
  #
  #
  # inputs:
  #
  # - G             (d x d matrix)      Gram matrix
  # - cov_mat       (pd x pd matrix)    current covariance matrix for b_i's 
  # - basis_list
  # - time_grid
  # - betas         (d-dim vector)
  #
  #
  # outputs:
  #
  # - eigen_result (p-dim list)    each list has a list of eigenvalues and eigenvectors
  #
  # ----------------------------------------------------------------------------
  
  
  d <- dim(G)[1]
  p <- dim(cov_mat)[1] / d
  m <- length(time_grid)
  

  phi <- trig_basis_realization(basis_list, time_grid) # (m x d)
  
  eigen_result <- list()
  
  for(i in 1:p){ # for each process
    
    # 1) (dxd) eigendecomposition

    sigma_ii <- extract_block_structure_ij(cov_mat, d, i, i)
    
    
    
    eigen_mat <- sigma_ii %*% G
    if(all(eigen_mat == diag(d))){  # the case where our basis is orthonormal already
      eigen_result_i <- list(values = rep(1, d),
                             vectors = diag(d))
    } else{
      eigen_result_i <- eigen(sigma_ii %*% G)
    }
    
    
    # 2) retrieval of eigenfunctions - when sigma_ii = identity and G = identity, then eigenfunctions are just the basis functions

    eigenfunction_i <- phi %*% eigen_result_i$vectors   # (m x d) times (d x 2)
    
    # 3) KL coefficients - when basis functions are orthonormal, KL coeffs = betas
    
    betas_mat <- matrix(betas, nrow = length(betas), ncol = 1)
    KL_i <- eigen_result_i$vectors %*% betas
    
    # 4) Correlation between KL coefficients
    
    eigen_result_i$eigenfunction <- eigenfunction_i
    eigen_result_i$KL_coeffs <- KL_i
    
    eigen_result[[i]] <- eigen_result_i
    
  }
  
  # correlation of KL coefficients (KL_cov_result)
  
  KL_cov_result <- list()
  
  for(i in 1:p){
    for(j in i:p){
      key <- paste0(i, '_', j)
      evecs_i <- eigen_result[[i]]$vectors
      evecs_j <- eigen_result[[j]]$vectors
      cov_ij <- t(evecs_i) %*% extract_block_structure_ij(cov_mat, d, i, j) %*% evecs_j
      
      denom <- outer(eigen_result[[i]]$values, eigen_result[[j]]$values)
      
      cor_ij <- cov_ij/denom
      KL_cov_result[[key]] <- cor_ij
    }
  }
  
  # correlation operator (C_cond)
  
  C_cond <- list()
  
  for(i in 1:p){
    for(j in i:p){
      
      key <- paste0(i, '_', j)
      coeffs <- KL_cov_result[[key]]
      eigenfunction_i <- eigen_result[[i]]$eigenfunction
      eigenfunction_j <- eigen_result[[j]]$eigenfunction  
      
      C_ij <- matrix(0, nrow = m, ncol = m)

      
      for(a in 1:d){
        for(b in 1:d){
          
          
          C_ij <- C_ij + coeffs[a,b] * tcrossprod(eigenfunction_i[, a], eigenfunction_j[, b])
          
        }
      }
      C_cond[[key]] <- C_ij
    }
  }
  
  
  

  
  return(list(eigen_decomp = eigen_result,
              KL_cov = KL_cov_result,
              corr_op = C_cond))
}







