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

# initialize the prec mat
trig_basis_prec_mat <- function(d, p, y_c_k, adj_type, adj_params){
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: define the precision matrix to generate beta's 
  #
  # inputs:
  #
  # - d                 (integer)
  # - p                 (integer)
  # - y_c_k             (q_c-dim vector)
  # - adj_type          (string)
  # - adj_params        (vector)
  #
  #
  # outputs:
  #
  # cov_mat   (pd x pd matrix)
  #
  # ----------------------------------------------------------------------------
  
  if(! adj_type %in% c('block_banded_c0', 
                       'block_banded_c2', 
                       'block_banded_v2',
                       'flexible_block_banded_c0', 
                       'flexible_block_banded_c2',
                       'flexible_block_banded_v2')){
    stop('Error 21b: adj_type not available')
  }
  if(! length(y_c_k) == 1){
    stop('Error 21b: y_c_k should be a scalar')
  }
  
  
  if(adj_type %in% c('block_banded_v2', 'block_banded_c2', 'block_banded_c0')){
    
    # Theta_{i,i}   = beta_var * I_d 
    # Theta_{i,i+1} = J_2_const * [1 0; 0 -1]
    
    # calculating J_2_const
    if(adj_type == 'block_banded_v2'){
      y_c_min <- adj_params[1]
      y_c_max <- adj_params[2]
      c_min <- adj_params[3]
      c_max <- adj_params[4]
      beta_var <- adj_params[5] 
      
      J_2_const <- c_min + (c_max - c_min) * (y_c_k - y_c_min) / (y_c_max - y_c_min)
    } else if(adj_type == 'block_banded_c2'){
      J_2_const <- adj_params[3]
      beta_var  <- adj_params[4] 
    } else{
      # block_banded_c0
      J_2_const <- adj_params[2]
      beta_var  <- adj_params[3] 
    }
    
    J_2 <- J_2_const * (-1)^(1 + 1:d)
    off_block <- diag(J_2)
    
    on_block <- diag(d) * beta_var
    
    theta_pd <- matrix(0, nrow = p*d, ncol = p*d)
    
    for(i in 1:p){
      for(j in 1:p){
        # Compute index ranges for block (i,j)
        row_idx <- ((i - 1) * d + 1):(i * d)
        col_idx <- ((j - 1) * d + 1):(j * d)
        
        
        if(i == j){
          theta_pd[row_idx, col_idx] <- on_block
        }
        
        if(abs(i-j) == 1){
          theta_pd[row_idx, col_idx] <- off_block
        }
      }
    }
    
    return(theta_pd)
  }
  
    

  if(adj_type %in% c('flexible_block_banded_v2', 'flexible_block_banded_c2', 'flexible_block_banded_c0')){
    
    # Theta_{i,i}   = [c1 0; 0 c2]
    # Theta_{i,i+1} = [c3 0; 0 c4]

    # 1) extract c1, c2, c3, c4 constants

    if(adj_type == 'flexible_block_banded_c0'){
      c1 <- adj_params[2]
      c2 <- adj_params[3]
      c3 <- adj_params[4] 
      c4 <- adj_params[5]
    } else if(adj_type == 'flexible_block_banded_c2'){
      c1 <- adj_params[3]
      c2 <- adj_params[4]
      c3 <- adj_params[5] 
      c4 <- adj_params[6]
    } else if(adj_type == 'flexible_block_banded_v2'){
      y_c_min <- adj_params[1]
      y_c_max <- adj_params[2]
      c1      <- adj_params[3]
      c2      <- adj_params[4]
      c3_min  <- adj_params[5]
      c3_max  <- adj_params[6]
      c4_min  <- adj_params[7]
      c4_max  <- adj_params[8]

      # interpolation
      c3 <- c3_min + (c3_max - c3_min) * (y_c_k - y_c_min) / (y_c_max - y_c_min)
      c4 <- c4_min + (c4_max - c4_min) * (y_c_k - y_c_min) / (y_c_max - y_c_min)
      
    } else{
      stop('Error 21b: adj_type not available')
    }
      
    # 2) assemble on and off-block matrics and the pd x pd matrix
    
    off_block <- diag(c(c3, c4))
    on_block  <- diag(c(c1, c2))
    
    theta_pd <- matrix(0, nrow = p*d, ncol = p*d)
    
    for(i in 1:p){
      for(j in 1:p){
        # Compute index ranges for block (i,j)
        row_idx <- ((i - 1) * d + 1):(i * d)
        col_idx <- ((j - 1) * d + 1):(j * d)
        
        
        if(i == j){
          theta_pd[row_idx, col_idx] <- on_block
        }
        
        if(abs(i-j) == 1){
          theta_pd[row_idx, col_idx] <- off_block
        }
      }
    }
    
    return(theta_pd)

  }
}

# invert prec_mat --> cov_mat
trig_basis_cov_mat <- function(d, p, y_c_k, adj_type, adj_params){
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: define the covariance matrix to generate beta's 
  #
  # inputs:
  #
  # - d                 (integer)
  # - p                 (integer)
  # - y_c_k             (q_c-dim vector)
  # - adj_type          (string)
  # - adj_params        (vector)
  #
  #
  # outputs:
  #
  # cov_mat   (pd x pd matrix)
  #
  # ----------------------------------------------------------------------------
  
  theta_pd <- trig_basis_prec_mat(d, p, y_c_k, adj_type, adj_params)
    
  cov_mat <- sym(solve(theta_pd))
    
  return(cov_mat)

}

# invert prec_mat --> cov_mat --> cor_mat
trig_basis_cor_mat <- function(d, p, y_c_k, adj_type, adj_params){
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: define the corrleation matrix that generated beta's 
  #
  # inputs:
  #
  # - d                 (integer)
  # - p                 (integer)
  # - y_c_k             (q_c-dim vector)
  # - adj_type          (string)
  # - adj_params        (vector)
  #
  #
  # outputs:
  #
  # cov_mat   (pd x pd matrix)
  #
  # ----------------------------------------------------------------------------
  
  cov_mat <- trig_basis_cov_mat(d, p, y_c_k, adj_type, adj_params)
  
  cor_mat <- assemble_blockwise_correlation(cov_mat, p, d)
  #cor_mat2 <- cov2cor(cov_mat)
  return(cor_mat)
  
}

# invert prec_mat --> cov_mat --> cor_mat --> invert to get prec_2
trig_basis_prec_mat_normalized <- function(d, p, y_c_k, adj_type, adj_params){
  
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: define the normalized precision matrix that generated beta's 
  #
  # inputs:
  #
  # - d                 (integer)
  # - p                 (integer)
  # - y_c_k             (q_c-dim vector)
  # - adj_type          (string)
  # - adj_params        (vector)
  #
  #
  # outputs:
  #
  # prec_mat_normalized   (pd x pd matrix)
  #
  # ----------------------------------------------------------------------------
  
  cor_mat <- trig_basis_cor_mat(d, p, y_c_k, adj_type, adj_params)
  
  prec_mat_normalized <- sym(solve(cor_mat))
  
  return(prec_mat_normalized)
  
}

# prec_mat --> HS norm --> adj_mat
trig_basis_adj_mat <- function(d, p, y_c_k, adj_type, adj_params, thresh = 1e-3){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: define the adjacency matrix behind a covariance matrix through its precision HS
  #
  # inputs:
  #
  # - d                 (integer)
  # - p                 (integer)
  # - y_c_k             (q_c-dim vector)
  # - adj_type          (string)
  # - adj_params        (vector)
  #
  #
  # outputs:
  #
  # adj_mat   (p x p matrix)  0's on the diagonal
  #
  # ----------------------------------------------------------------------------

  theta_pd <- trig_basis_prec_mat(d, p, y_c_k, adj_type, adj_params)
  
  HS_mat <- hilbert_schmidt_norm_pm(theta_pd, p, d)
  
  adj_mat <- matrix(0, p, p)
  adj_mat[HS_mat > thresh] <- 1
  diag(adj_mat) <- 0
  
  return(adj_mat)

}

trig_basis_log_intensity <- function(cov_mat_list, basis_list, mu_t, time_grid, mean_vec){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: draw n log-intensities from the trig-basis data generation process
  #
  #       beta_i ~ N(mean_vec, cov_mat[[i]])
  #
  # input:
  #
  # - cov_mat_list  (n-dim list of pd x pd matrix)  covariance matrix for each subject
  # - basis_list    (list of d eigenfunctions)   
  # - mu_t          (m-dim vector)                baseline mean mu(t)
  # - time_grid     (m-dim vector)                time discretization
  # - mean_vec      (pd-dim vector)               mean vector for beta generation
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
    b <- mean_vec + t(L) %*% z  # pd x 1
    
    # reshape b into p x d matrix of coefficients
    B_mat <- matrix(b, nrow = p, ncol = d, byrow = TRUE)
    
    # compute log-intensity for each process i
    log_int_mat <- matrix(0, nrow = p, ncol = m)
    for (i in 1:p) {
      # linear combination: mu(t) + sum_k beta_k * b_ik * phi_k(t)
      log_int_mat[i, ] <- mu_t + colSums(matrix(B_mat[i, ], nrow = d, ncol = m, byrow = FALSE) * t(Phi))
      
    }
    
    log_intensities[[s]] <- log_int_mat
    beta_coefficients[[s]] <- B_mat
  }
  
  return(list(log_intensities = log_intensities,
              beta_coefficients = beta_coefficients))
  
}



# truths

trig_basis_rho_truth <- function(basis_list, mean_vec, cov_mat, time_grid, mu_t){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: calculate the ground truth rho_i(t), rho_ij(s,t), and G_{ij}(s,t)
  #
  # - recall that X_i(t) = mu(t) + sum_{k=1}^d beta_{ik} * phi^k(t)
  # - recall that beta_i ~ N(m_i, [cov]_{ii})
  # - where phi^k(t) are the k orthonormal basis functions
  # - and beta_{ik} are normally distributed with a pd x pd matrix
  #
  # - E[exp(X_i(t))] = exp( mu(t) + \sum_{k=1}^d m_{ik} phi^k(t) +  0.5 \sum_{k=1}^d lambda_k [phi^k(t)]^2 )
  # 
  # input:
  # 
  # - basis_list       (d-dim list)                           output from trig_basis
  # - mean_vec         (pd-dim vector)                        mean vector for beta: beta ~ N(mean_vec, cov_mat)
  # - cov_mat          (pd x pm matrix)                       covariance matrix for beta: beta ~ N(mean_vec, cov_mat)
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
  
  lambda_vec <- diag(cov_mat) # eigenvalues

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
  #  
  #    term_2:   mu_beta * phi (usually 0)
  #    term_3b:  ????
  #    term_3:   phi %*% Sigma_ii %*% phi
  
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
  
  

  
  return(list(rho_i_truth = rho_i_truth_matrix,   # (p x m)
              rho_ij_truth = rho_ij_truth,        # (all pc2 lists of mxm matrices)
              g_ij_truth = g_ij_truth,            # (pc2 list of mxm matrices)
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

trig_basis_eigendecomposition <- function(G, cov_mat, cor_mat, prec_mat, basis_list, time_grid){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: find ground truths in the trig_basis framework
  #
  #
  # inputs:
  #
  # - G             (d x d matrix)      Gram matrix
  # - cov_mat       (pd x pd matrix)    current covariance matrix for betas 
  # - cor_mat       (pd x pd matrix)    current correlation matrix for betas 
  # - prec_mat      (pd x pd matrix)    current precision matrix for betas 
  # - basis_list
  # - time_grid     (m-dim vector)
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
  # - KL_cov        (pc2 list of dxd matrices)    cov(beta, beta)
  # - KL_cor        (pc2 list of dxd matrices)    cor(beta, beta)
  # - KL_prec       (pc2 list of dxd matrices)    prec(beta, beta)
  # - C_cond        (pc2 list of mxm matrices)    double sum of cor(beta, beta) * tensorprod(phi, phi)
  # - P_cond        (pc2 list of mxm matrices)    double sum of prec(beta, beta) * tensorprod(phi, phi)
  # - C_HS          (pxp matrix)                  HS norm of pm x pm matrix
  # - P_HS          (pxp matrix)                  HS norm of pm x pm matrix
  # - C_HS_unnorm   (pxp matrix)                  HS norm of pm x pm matrix
  # - P_HS_unnorm   (pxp matrix)                  HS norm of pm x pm matrix
  # - C_HS_KL       (pxp matrix)                  HS norm of pd x pd matrix of correlations
  # - P_HS_KL       (pxp matrix)                  HS norm of pd x pd matrix of precisions
  # - efunc_outer           (i_j list of mxm matrices)
  # - efunc_outer_unnorm    (i_j list of mxm matrices)
  #
  #
  # ----------------------------------------------------------------------------
  
  # 0) parameters
  d <- dim(G)[1]
  p <- dim(cov_mat)[1] / d
  m <- length(time_grid)
  delta <- 1/m

  # 1) obtain ground truth eigenfunction basis
  phi <- trig_basis_realization(basis_list, time_grid) # (m x d)
  
  eigen_result <- list()
  
  for(i in 1:p){ # for each process
    
    # 3) eigendecomposition of the (i, i) block, usually the identity matrix

    
    if(any(G != diag(d))){
      stop('Error: gram matrix says eigenfunctions are not orthogonal')
    }
    
    sigma_ii <- extract_block_structure_ij(cov_mat, d, i, i)  # (d x d)
    
    # 4) since G is identity, we manually create the eigendecomposition
    eigen_result_i <- list(values         = diag(sigma_ii),
                           eigenfunctions = phi)
    
    
    eigen_result[[i]] <- eigen_result_i
    
  }
  
  
  # 5) correlation operator (C_cond) in basis space (d-dim) and in regular space (m-dim)

  C_cond <- list()
  P_cond <- list()
  
  efunc_outer <- list()
  
  C_cond_unnorm <- list()
  P_cond_unnorm <- list()
  
  efunc_outer_unnorm <- list()
  
  KL_cov <- extract_block_structure_v2(cov_mat, p, d)   # (pd x pd) --> list of (d x d) matrices
  KL_cor  <- extract_block_structure_v2(cor_mat, p, d) 
  KL_prec <- extract_block_structure_v2(prec_mat, p, d) 
  
  for(i in 1:p){
    for(j in i:p){
      
      key <- paste0(i, '_', j)
      
      cor_ij <- KL_cor[[key]]  #(d x d)
      prec_ij <- KL_prec[[key]]  #(d x d)
      
      eigenfunction_i_norm <- eigen_result[[i]]$eigenfunctions * sqrt(delta)  # (m x d) normalized eigenfunction
      eigenfunction_j_norm <- eigen_result[[j]]$eigenfunctions * sqrt(delta)
      
      eigenfunction_i_unnorm <- eigen_result[[i]]$eigenfunctions # (m x d) normalized eigenfunction
      eigenfunction_j_unnorm <- eigen_result[[j]]$eigenfunctions 
      
      
      # take the outer products: sum_{a, b} phi_a * cor_ij(a, b) * phi_b
      
      C_ij <- eigenfunction_i_norm %*% cor_ij %*% t(eigenfunction_j_norm)  # (m x d) * (d x d) * (d x m)
      P_ij <- eigenfunction_i_norm %*% prec_ij %*% t(eigenfunction_j_norm) # (m x d) * (d x d) * (d x m)

      C_ij_unnorm <- eigenfunction_i_unnorm %*% cor_ij %*% t(eigenfunction_j_unnorm)  # (m x d) * (d x d) * (d x m)
      P_ij_unnorm <- eigenfunction_i_unnorm %*% prec_ij %*% t(eigenfunction_j_unnorm) # (m x d) * (d x d) * (d x m)      
            
      efunc_outer_ij_norm <- eigenfunction_i_norm %*% t(eigenfunction_j_norm)   # (m x d) * (d x m)
      efunc_outer_ij_unnorm <- eigenfunction_i_unnorm %*% t(eigenfunction_j_unnorm)   # (m x d) * (d x m)
      
      C_cond[[key]] <- C_ij                       # (m x m)
      P_cond[[key]] <- P_ij                       # (m x m)
      efunc_outer[[key]] <- efunc_outer_ij_norm   # (m x m)
      
      C_cond_unnorm[[key]] <- C_ij_unnorm
      P_cond_unnorm[[key]] <- P_ij_unnorm
      efunc_outer_unnorm[[key]] <- efunc_outer_ij_unnorm
    }
  }
  
  
  
  # 6) HS_truth from KL_cov (d x d)
  
  P_HS_KL <- hilbert_schmidt_norm_pm(prec_mat, p, d)
  C_HS_KL <- hilbert_schmidt_norm_pm(cor_mat, p, d)
  
  # 7) HS_truth from C_cond (m x m)
  
  C_cond_full <- assemble_block_matrix_v2(C_cond, p, m) 
  P_cond_full <- assemble_block_matrix_v2(P_cond, p, m)
  
  C_cond_full_unnorm <- assemble_block_matrix_v2(C_cond_unnorm, p, m) 
  P_cond_full_unnorm <- assemble_block_matrix_v2(P_cond_unnorm, p, m)
  
  C_HS <- hilbert_schmidt_norm_pm(C_cond_full, p, m)
  P_HS <- hilbert_schmidt_norm_pm(P_cond_full, p, m)
  C_HS_unnorm <- hilbert_schmidt_norm_pm(C_cond_full_unnorm, p, m)
  P_HS_unnorm <- hilbert_schmidt_norm_pm(P_cond_full_unnorm, p, m)
  

  
  # 8) reorder for step_4
  
  eigen_result_v2 <- list(
    eigenvalues   = lapply(eigen_result, `[[`, "values"),
    eigenfunctions = lapply(eigen_result, `[[`, "eigenfunctions"),
    n_dims        = lapply(eigen_result, function(x) dim(x$eigenfunctions)[2])
  )
  
  
  # 9) return
  
  return(list(eigen_decomp = eigen_result_v2,
              
              # (i_j list of dxd matrices)
              KL_cov = KL_cov,                           
              KL_cor = KL_cor,                           
              KL_prec = KL_prec,  
              
              # (i_j list of mxm matrices)  
              C_cond = C_cond,                          
              P_cond = P_cond,                           
              C_cond_unnorm = C_cond_unnorm,             
              P_cond_unnorm = P_cond_unnorm,             
              efunc_outer = efunc_outer,
              efunc_outer_unnorm = efunc_outer_unnorm,
              
              # (pxp matrix)
              C_HS = C_HS,                               
              P_HS = P_HS,
              C_HS_unnorm = C_HS_unnorm,
              P_HS_unnorm = P_HS_unnorm, 
              C_HS_KL = C_HS_KL,
              P_HS_KL = P_HS_KL
              ))               
}


trig_basis_eigendecomposition_beta_truths <- function(beta_coeffs, eigen_decomp_truth, Y_c, y_c_query_k){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: find ground truths from realized beta values
  #
  #
  # inputs:
  #
  #
  # - beta_coeffs        (p x d x n matrix of realized beta coefficients)
  # - eigen_decomp_truth (list of 3 things)
  #   - eigenvalues         (p-dim list of d_i-dim vector of eigenvalues)
  #   - eigenfunctions      (p-dim list of m x d_i matrices of eigenfunctions)
  #   - n_dims              (p-dim list of d_i scalars)
  #
  # - Y_c               (n x q_c matrix)
  # - y_c_query_k       (scalar)
  #
  # outputs:
  # 
  # - list of the following:
  #
  # - C_cond_beta_truth
  # - C_cond_beta_truth_unnorm
  # - P_cond_beta_truth
  # - P_cond_beta_truth_unnorm
  # - C_HS_beta_truth
  # - P_HS_beta_truth
  #
  #
  # ----------------------------------------------------------------------------
  
  # 1) get weights for this cont_ind
  
  p <- dim(beta_coeffs)[1]
  d <- dim(beta_coeffs)[2]
  n <- dim(beta_coeffs)[3]
  m <- dim(eigen_decomp_truth$eigenfunctions[[1]])[1]
  weights_k <- KDE_weights(Y_c, y_c_query_k)
  
  # 2) get unbiased weighted covariance matrix estimate (pd x pd) 
  
  X <- aperm(beta_coeffs, c(2, 1, 3)) %>% 
    matrix(nrow = p*d, ncol = n) %>% # reshape to p*d × n
    t()  # (n x pd)
  
  # weighted mean (length p*d)
  mu <- colSums(weights_k * X)
  
  # centered data
  XC <- sweep(X, 2, mu)
  
  # weighted covariance: sum_i w_i (x_i - mu)(x_i - mu)^T
  cov_w <- t(XC * weights_k) %*% XC
  
  kappa <- 1 - sum(weights_k^2)
  cov_w_unbiased <- cov_w / kappa 
  
  
  
  # 3) convert to correlation using blockwise function
  
  cor_w_unbiased <- assemble_blockwise_correlation(cov_w_unbiased, p, d)
  prec_w_unbiased <- sym(solve(cor_w_unbiased))
  
  KL_beta_truth <- list(KL_cov = cov_w_unbiased,
                        KL_cor = cor_w_unbiased,
                        KL_prec = prec_w_unbiased)
  
  # 3c) compute C_cond_X_truth 
  

  KL_cov_beta_truth           <- KL_beta_truth$KL_cov %>% extract_block_structure_v2(p, d)
  KL_cor_beta_truth           <- KL_beta_truth$KL_cor %>% extract_block_structure_v2(p, d)
  KL_prec_beta_truth          <- KL_beta_truth$KL_prec %>% extract_block_structure_v2(p, d)
  C_cond_list_beta_truth      <- correlation_estimation_KL_cor(eigen_decomp_truth, KL_cor_beta_truth) 
  P_cond_list_beta_truth      <- correlation_estimation_KL_cor(eigen_decomp_truth, KL_prec_beta_truth) 
  
  C_cond_beta_truth_full          <- C_cond_list_beta_truth$C_cond %>% assemble_block_matrix_v2(p, m)
  C_cond_beta_truth_unnorm_full   <- C_cond_list_beta_truth$C_cond_unnorm %>% assemble_block_matrix_v2(p, m)
  
  # step 10
  P_cond_beta_truth_full          <- P_cond_list_beta_truth$C_cond %>% assemble_block_matrix_v2(p, m)
  P_cond_beta_truth_unnorm_full   <- P_cond_list_beta_truth$C_cond_unnorm %>% assemble_block_matrix_v2(p, m)
  
  P_cond_beta_truth               <- extract_block_structure_v2(P_cond_beta_truth_full, p, m)
  P_cond_beta_truth_unnorm        <- extract_block_structure_v2(P_cond_beta_truth_unnorm_full, p, m)
  
  # step 11
  C_HS_beta_truth <- hilbert_schmidt_norm_pm(KL_beta_truth$KL_cor, p, d)
  P_HS_beta_truth <- hilbert_schmidt_norm_pm(KL_beta_truth$KL_prec, p, d)
  
  step_9_10_11_beta_truth <- list(KL_cov_beta_truth = KL_cov_beta_truth,
                                  KL_cor_beta_truth = KL_cor_beta_truth,
                                  KL_prec_beta_truth = KL_prec_beta_truth,
                                  C_cond_beta_truth = C_cond_list_beta_truth$C_cond,
                                  C_cond_beta_truth_unnorm = C_cond_list_beta_truth$C_cond_unnorm,
                                  P_cond_beta_truth = P_cond_beta_truth,
                                  P_cond_beta_truth_unnorm = P_cond_beta_truth_unnorm,
                                  C_HS_beta_truth = C_HS_beta_truth,
                                  P_HS_beta_truth = P_HS_beta_truth)
  
  return(step_9_10_11_beta_truth)
                               
}




