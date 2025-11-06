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
  # - adj_type          (string)
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
      beta_var <- adj_params[5] 
      
      J_2_const <- c_min + (c_max - c_min) * (y_c_k - y_c_min) / (y_c_max - y_c_min)
    } else{
      J_2_const <- adj_params[3]
    }
    
    J_2 <- J_2_const * (-1)^(1 + 1:d)
    off_block <- diag(J_2)
    
    on_block <- diag(d) * beta_var
    
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

trig_basis_log_intensity <- function(cov_mat_list, basis_list, mu_t, time_grid){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: draw n log-intensities from the trig-basis data generation process
  #
  #
  # input:
  #
  # - cov_mat_list  (n-dim list of pd x pd matrix)  covariance matrix for each subject
  # - basis_list    (list of d eigenfunctions)   
  # - mu_t          (m-dim vector)                baseline mean mu(t)
  # - time_grid     (m-dim vector)                time discretization
  #
  #
  # output:
  # 
  # - list of the following:
  #   - log_intensities      (n-dim list of p x m_both matrices) each item is for the k-th subject across all p processes
  #   - beta_coefficients    (n-dim list of p x d matrices)      each item is for the k-th subject across all p processes, realizations of beta
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
  beta_coefficients <- vector("list", n)   # realizations of beta ~ N(0, cov_mat)
  
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
      # linear combination: mu(t) + sum_k beta_k * b_ik * phi_k(t)
      #log_int_mat[i, ] <- mu_t + colSums(matrix(betas * B_mat[i, ], nrow = d, ncol = m, byrow = FALSE) * t(Phi))
      log_int_mat[i, ] <- mu_t + colSums(matrix(B_mat[i, ], nrow = d, ncol = m, byrow = FALSE) * t(Phi))
      
    }
    
    log_intensities[[s]] <- log_int_mat
    beta_coefficients[[s]] <- B_mat
  }
  
  return(list(log_intensities = log_intensities,
              beta_coefficients = beta_coefficients))
  
}



# truths



trig_basis_rho_truth <- function(basis_list, mean_vec, time_grid, mu_t, y_c_query, adj_type, adj_params){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: calculate the ground truth rho_i(t), rho_ij(s,t), and G_{ij}(s,t)
  #
  # - recall that X_i(t) = mu(t) + sum_{k=1}^d beta_{ik} * phi^k(t)
  # - recall that beta_i ~ N(m_i, [cov]_{ii})
  # - where phi^k(t) are the k orthonormal basis functions
  # - and beta_{ik} are normally distributed with a pd x pd matrix
  #
  # - E[exp(X_i(t))] = exp(mu(t) + \sum_{k=1}^d m_{ik} phi^k(t) +  0.5 \sum_{k=1}^d lambda_k [phi^k(t)]^2 )
  # 
  # input:
  # 
  # - basis_list       (d-dim list)                           output from trig_basis
  # - mean_vec         (pd-dim vector)                        mean vector for beta: beta ~ N(mean_vec, cov_mat)
  # - time_grid        (m-dim vector)                         time discretization
  # - mu_t             (m-dim vector)                         mu(t) when defining X_i(t)
  #
  #
  # output:
  #
  # - rho_i_truth       (p x m matrix)   marginal intensity for all p processes 
  # - rho_ij_Truth      (list of i_j entries)  each is mxm matrix
  # - g_ij_truth        (list of i_j entries)  each is mxm matrix
  #
  # ----------------------------------------------------------------------------
  
  # 1) generate the cov_mat (pd x pd matrix) for the following y_c_query
  

  d <- length(basis_list)
  p <- length(mean_vec) / d
  
  cov_mat <- trig_basis_cov_mat(d, p, y_c_query, adj_type, adj_params)  # (pd x pd)
  lambda_vec <- diag(cov_mat)

  # 2) get phi(t) and phi^2(t) ready
  
  phi <- trig_basis_realization(basis_list, time_grid)  # (m x d matrix) 
  phi_squared <- phi * phi # (m x d)
  phi_outer <- lapply(1:d, function(i) outer(phi[ ,i], phi[, i]))  # (d-dim list of mxm matrices)
  
  # 3) split mean_vec and lambda_vec into a list of p d-dim vectors
  #    lambda_vec = diag of Sigma
  #    cov_ii_list = list of the dxd diagonal blocks
  
  mean_list <- split(mean_vec, rep(1:p, each = d))       # (p-dim list of d-dim vectors)
  lambda_list <- split(lambda_vec, rep(1:p, each = d))
  cov_ii_list <- lapply(1:p, function(i) extract_block_structure_ij(cov_mat, d, i, i))  # p-dim list of (d x d) matrices
  
  # 4) obtain values for rho_i_truth
  
  term_2 <-  lapply(mean_list, function(v) as.numeric(t(v) %*% t(phi)))              # M_i(t)  (1 x d) * (d x)
  term_3b <- lapply(lambda_list, function(v) as.numeric(t(v) %*% t(phi_squared)))     # V_i(t)
  term_3 <- lapply(1:p, function(i) rowSums((phi %*% cov_ii_list[[i]]) * phi))        # each phi(t) is d x 1, so we have (1 x d) * (d x d) * (d x 1) for a scalar at each timepoint
  
  rho_i_truth <- lapply(1:p, function(i){exp(mu_t + term_2[[i]] + 0.5 * term_3[[i]])})  
  rho_i_truth_matrix <- do.call(rbind, rho_i_truth)
  
  # 5) obtain values for rho_ij_truth

  rho_ij_truth <- list()
  g_ij_truth   <- list()
  for(i in 1:p){
    for(j in i:p){
      key <- paste0(i, '_', j)
      
      # need the covariance between beta_k of processes i and j
      # this looks like taking the diagonal of the (i ,j)-th block
      
      # gamma_k_ij <- extract_block_structure_ij(cov_mat, d, i, j) %>% diag()
      sigma_ij <- extract_block_structure_ij(cov_mat, d, i, j)
      
      #g_ij - linear combination of phi_outer with the gamma's
      #g_ij_truth[[key]] <- Reduce("+", Map("*", gamma_k_ij, phi_outer))
      g_ij_truth[[key]] <- phi %*% sigma_ij %*% t(phi)          # (m x d) * (d x d) * (d x m)
      
      # rho_ij can use g_ij and rho_i results
      rho_ij_truth[[key]] <- tcrossprod(rho_i_truth[[i]], rho_i_truth[[j]]) * exp(g_ij_truth[[key]])

    }
  }
  
  

  
  return(list(rho_i_truth = rho_i_truth_matrix,
              rho_ij_truth = rho_ij_truth,
              g_ij_truth = g_ij_truth,
              mu_t = mu_t,
              term_2 = term_2,
              term_3 = term_3,
              term_3b = term_3b))
  

}



trig_basis_cross_covariance_truth <- function(basis_list, cov_mat, time_grid){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: calculate the ground truth cross-covariance (pm x pm) for any time discretization setting
  #
  # - Recall: Cov(X_i(s), X_j(t)) = \phi(s)^\top [Sigma]_{ij} \phi(t)
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
  # - G_ij_list      (list of mxm matrices) G_ij(s,t) for i_j where i <= j
  #
  #
  # ----------------------------------------------------------------------------
  
  m <- length(time_grid)
  phi <- trig_basis_realization(basis_list, time_grid) # (m x d matrix)
  d <- dim(phi)[2]
  p <- dim(cov_mat)[1] / d
  
  G_ij_list <- list()
  
  for(i in 1:p){
    for(j in i:p){
      key <- paste0(i, '_', j)
      G_ij_list[[key]] <- phi %*% extract_block_structure_ij(cov_mat, d, i, j) %*% t(phi)
    }
  }
  
  
  return(G_ij_list)
  
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
  # - time_grid     (m-dim vector)
  # - betas         (p x d x n)        realization of beta ~ N(0, cov_mat) for all n subjects 
  #
  # outputs:
  #
  # - eigen_result (list)  each list contains the following truths:
  #
  #   - eigen_decomp (list of 3 things)
  #     - eigenvalues         (p-dim list of d_i-dim vector of eigenvalues)
  #     - eigenfunctions      (p-dim list of m x d_i matrices of eigenfunctions)
  #     - n_dims              (p-dim list of d_i scalars)
  #
  #   - KL_cor      (list of i_j items) each item is a d x d correlation matrix
  #   - corr_op     (pm x pm matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  
  d <- dim(G)[1]
  p <- dim(cov_mat)[1] / d
  m <- length(time_grid)
  

  phi <- trig_basis_realization(basis_list, time_grid) # (m x d)
  
  eigen_result <- list()
  
  for(i in 1:p){ # for each process
    
    # 1) (dxd) eigendecomposition

    sigma_ii <- extract_block_structure_ij(cov_mat, d, i, i)  # (d x d)
    
    
    
    eigen_mat <- sigma_ii %*% G
    if(all(eigen_mat == diag(d))){  # the case where our basis is orthonormal already
      eigen_result_i <- list(values = rep(1, d),
                             vectors = diag(d))
    } else{
      eigen_result_i <- eigen(sigma_ii %*% G)
    }
    
    
    # 2) retrieval of eigenfunctions - when sigma_ii = identity and G = identity, then eigenfunctions are just the basis functions

    eigenfunction_i <- phi %*% eigen_result_i$vectors   # (m x d) times (d x 2)
    
    
    # 4) storing
    
    eigen_result_i$eigenfunctions <- eigenfunction_i
    
    eigen_result[[i]] <- eigen_result_i
    
  }
  
  # correlation of KL coefficients (KL_cor_result)
  
  # - for eigencomponents m and n
  # - for processes i and j
  # - cor(beta_i^m, beta_j^n) = [\Sigma_{ij}]_{mn} / sqrt([\Sigma_{ii}]_{mm} [\Sigma_{jj}]_{nn})
  
  KL_cov_result <- extract_block_structure_v2(cov_mat, p, d)
  
  KL_cor_result <- list()
  
  
  for(i in 1:p){
    for(j in i:p){
      
      key <- paste0(i, '_', j)
      
      # extract blocks
      rows_i <- ((i-1)*d + 1):(i*d)
      cols_j <- ((j-1)*d + 1):(j*d)
      
      Sigma_ij <- cov_mat[rows_i, cols_j]       # d x d
      Sigma_ii <- cov_mat[rows_i, rows_i]       # d x d
      Sigma_jj <- cov_mat[cols_j, cols_j]       # d x d
      
      # elementwise correlation
      R_ij <- Sigma_ij / sqrt(outer(diag(Sigma_ii), diag(Sigma_jj)))
      
      # store in list
      KL_cor_result[[key]] <- R_ij
    }
  }
  
  # correlation operator (C_cond)
  # - linear combination of KL covariances and tensor product of eigenfunctions
  
  
  C_cond <- list()
  
  for(i in 1:p){
    for(j in i:p){
      
      key <- paste0(i, '_', j)
      coeffs <- KL_cor_result[[key]]  #(d x d)
      eigenfunction_i <- eigen_result[[i]]$eigenfunction
      eigenfunction_j <- eigen_result[[j]]$eigenfunction  
      
      C_ij <- matrix(0, nrow = m, ncol = m)

      
      for(a in 1:d){
        for(b in 1:d){
          
          
          C_ij <- C_ij + coeffs[a,b] * tcrossprod(eigenfunction_i[, a], eigenfunction_j[, b]) # (m x m)
          
        }
      }
      C_cond[[key]] <- C_ij
    }
  }
  
  C_cond_full <- assemble_block_matrix_v2(C_cond, p, m) # (pm x pm)
  
  
  # reordering 
  
  # reorder for step_4
  
  eigen_result_v2 <- list(
    eigenvalues   = lapply(eigen_result, `[[`, "values"),
    eigenfunctions = lapply(eigen_result, `[[`, "eigenfunctions"),
    n_dims        = lapply(eigen_result, function(x) dim(x$eigenfunction)[2])
  )
  
  

  
  return(list(eigen_decomp = eigen_result_v2,
              KL_cor = KL_cor_result,
              KL_cov = KL_cov_result,
              corr_op = C_cond_full))
}







