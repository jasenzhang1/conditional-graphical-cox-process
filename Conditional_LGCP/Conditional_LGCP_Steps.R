library(parallel)
library(kernlab)  # For Gaussian kernel
library(pracma)   # For trapz()

# Assumes 'get_rho_diag_pp_reproduce()' was already run
# Suppose the result is:
# out <- get_rho_diag_pp_reproduce(...)

# STEP 1: Stratify
stratify_by_discrete <- function(df, discrete_vars) {
  group_split(df, across(all_of(discrete_vars)))
}

# STEP 2–3: Compute rho_i and rho_ij and log-covariance
compute_rho_and_cov <- function(stratum, Tseq, feature_sel, res_list) {
  p <- length(feature_sel)
  n_stratum <- length(stratum)
  
  rho_list <- vector("list", p)
  rho_ij_list <- vector("list", p^2)
  G_list <- vector("list", p^2)
  
  for (i in 1:p) {
    xi_i <- sapply(stratum, function(k) res_list[[i]]$data_i_count$count[res_list[[i]]$data_i_count$subject_num == k])
    Gamma_i_list <- lapply(stratum, function(k) res_list[[i]]$intensity[, res_list[[i]]$res$subject_num == k])
    rho_list[[i]] <- Reduce(`+`, Map(`*`, xi_i, Gamma_i_list)) / n_stratum
  }
  
  for (i in 1:p) {
    for (j in 1:p) {
      idx <- (i - 1) * p + j
      xi_i <- sapply(stratum, function(k) res_list[[i]]$data_i_count$count[res_list[[i]]$data_i_count$subject_num == k])
      xi_j <- sapply(stratum, function(k) res_list[[j]]$data_i_count$count[res_list[[j]]$data_i_count$subject_num == k])
      Gamma_ij_list <- lapply(stratum, function(k) res_list[[i]]$diag_mat)  # for simplicity assume symmetry
      rho_ij_list[[idx]] <- Reduce(`+`, Map(function(a, b, mat) a * b * mat, xi_i, xi_j, Gamma_ij_list)) / n_stratum
      G_list[[idx]] <- log(rho_ij_list[[idx]] / (outer(rho_list[[i]], rho_list[[j]], "*") + 1e-8))
    }
  }
  
  list(rho = rho_list, rho_ij = rho_ij_list, G = G_list)
}

# STEP 4: Eigendecomposition
eigen_decompose <- function(G, dmax) {
  eig <- eigen(G, symmetric = TRUE)
  eig$values[eig$values < 0] <- 0
  list(values = eig$values[1:dmax], vectors = eig$vectors[, 1:dmax])
}

# STEP 5: KL coefficients
estimate_KL <- function(X_t, eigen_vecs, t_grid) {
  apply(eigen_vecs, 2, function(phi) trapz(t_grid, X_t * phi))
}

# STEP 6: Kernel matrix
rbf_kernel <- function(x, y, sigma = 1) {
  exp(-sum((x - y)^2) / (2 * sigma^2))
}

kernel_matrix <- function(Y_c, kernel_func = rbf_kernel) {
  n <- nrow(Y_c)
  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- kernel_func(Y_c[i, ], Y_c[j, ])
    }
  }
  K
}

# STEP 7: Regression Operator
estimate_operator <- function(K, V, gamma) {
  solve(K + gamma * diag(nrow(K)), V)
}

# STEP 8: Conditional Covariance
estimate_conditional_cov <- function(M, eta_i, eta_j, y_c) {
  d <- ncol(eta_i)
  cov_mat <- matrix(0, length(eta_i[,1])^2, 1)
  for (a in 1:d) {
    for (b in 1:d) {
      f_ab <- eta_i[, a] %o% eta_j[, b]
      m_ab <- M[, (a - 1) * d + b]
      cov_mat <- cov_mat + m_ab * as.vector(f_ab)
    }
  }
  matrix(cov_mat, nrow = length(eta_i[,1]))
}

# Step 9: Normalize covariance matrices to get correlation matrices
normalize_cov_to_corr <- function(V_list, gamma1) {
  p <- length(V_list)  # V_list is p x p list of d x d matrices
  C_list <- vector("list", p)
  for (i in 1:p) {
    C_list[[i]] <- vector("list", p)
    Vii <- V_list[[i]][[i]] + gamma1 * diag(nrow(V_list[[i]][[i]]))
    Vii_inv_sqrt <- pracma::sqrtm(solve(Vii))
    for (j in 1:p) {
      if (i == j) {
        C_list[[i]][[j]] <- diag(nrow(Vii))  # identity
      } else {
        Vjj <- V_list[[j]][[j]] + gamma1 * diag(nrow(V_list[[j]][[j]]))
        Vjj_inv_sqrt <- pracma::sqrtm(solve(Vjj))
        C_list[[i]][[j]] <- Vii_inv_sqrt %*% V_list[[i]][[j]] %*% Vjj_inv_sqrt
      }
    }
  }
  C_list
}

# Step 10: Compute precision operator from block correlation matrix
compute_precision_operator <- function(C_list, gamma2) {
  p <- length(C_list)
  d <- nrow(C_list[[1]][[1]])
  bigC <- matrix(0, p * d, p * d)
  
  for (i in 1:p) {
    for (j in 1:p) {
      row_idx <- ((i - 1) * d + 1):(i * d)
      col_idx <- ((j - 1) * d + 1):(j * d)
      bigC[row_idx, col_idx] <- C_list[[i]][[j]]
    }
  }
  bigP <- solve(bigC + gamma2 * diag(p * d))
  bigP
}


# Step 11: Extract graph edges from precision operator
extract_graph <- function(precision_matrix, p, d, threshold) {
  edges <- list()
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      row_idx <- ((i - 1) * d + 1):(i * d)
      col_idx <- ((j - 1) * d + 1):(j * d)
      block <- precision_matrix[row_idx, col_idx]
      hs_norm <- sqrt(sum(block^2))
      if (hs_norm > threshold) {
        edges <- append(edges, list(c(i, j)))
      }
    }
  }
  do.call(rbind, edges)
}

make_adjacency_matrix <- function(edges, p) {
  A <- matrix(0, p, p)
  for (edge in edges) {
    A[edge[1], edge[2]] <- 1
    A[edge[2], edge[1]] <- 1
  }
  A
}
