
step_0_check <- function(dataset){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: visualize the ground truth precision matrices
  #
  # 
  # input:
  #
  # - dataset (dataset generated from simulation)
  #
  # output: 
  #
  # - warnings and g_01
  #
  # ----------------------------------------------------------------------------
  
  # 0.1) check ground truth precision pxp matrices
  
  prec_ground_truths <- lapply(dataset$true_graphs, function(item) item$P_block_kronecker$prec_mat_truth$prec_mat)
  
  # check if any prec mat is too nonnegative!
  check_psd <- sapply(1:length(prec_ground_truths), function(i) {
    min(eigen(prec_ground_truths[[i]], symmetric = TRUE, only.values = TRUE)$values) > -1e-10
  })  
  
  if(! any(check_psd)){
    warning("code 01: some ground truth prec mats are not psd, proceeding anyway")
  } 
  
  
}

step_0_keep_events <- function(dataset, k, i_vec){
  
  
  # ----------------------------------------------------------------------------
  #
  # keep point process events for subject k and processes i
  #
  #
  # input:
  #
  # - dataset    (list)
  # - k          (integer)          ubject id
  # - i_vec      (q-dim vector)  vector of process ID's
  #
  # 
  # outputs:
  # 
  # - events     (list of q vectors)
  #
  # ----------------------------------------------------------------------------
  
  keys <- paste0(k, '_', i_vec)
  
  events <- dataset$event_times[keys]
  
  return(events)
  
}

step_0_store <- function(dataset){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: store values before we begin estimating:
  # 
  # - ground truth precision matrices
  # - P_block_kronecker
  # - 
  #
  # 
  # input:
  #
  # - dataset (dataset generated from simulation)
  #
  # output: 
  #
  # - warnings and g_01
  #
  # ----------------------------------------------------------------------------
  
  # 0.1) check ground truth precision pxp matrices
  
  prec_ground_truths <- lapply(dataset$true_graphs, function(item) item$P_block_kronecker$prec_mat_truth$prec_mat)
  
  
  
  
  return(prec_ground_truths)
}

step_0_preprocess <- function(dataset){

  # ----------------------------------------------------------------------------
  #
  # GOAL: create `data_df4` to pass into estimation functions
  # 
  #
  # input:
  #
  # - dataset (dataset created from simulation)
  #
  #
  #
  # output:
  # 
  # - processed_data (list)
  #
  #   - [[1]] data_df4
  #   - [[2]] y_c_strata        all n y_c_strata
  #   - [[3]] query_y_cs
  #   - [[4]] patient_sel
  #   - [[5]] feature_sel
  #   - [[6]] y_c_strata_sel    subset of stratas in case some subjects are discarded
  #
  # ----------------------------------------------------------------------------
  

  
  # 3.1) first, convert dataset into a format that can be used for estimation
  
  data_df4 <- convert_data_for_estimation_event_times(dataset$event_times) 

  y_c_strata <- dataset$Y_continuous
  y_c_strata_sel <- dataset$Y_continuous_k
  
  query_y_cs <- dataset$simulation_params$query_y_cs  
  

  
  # size of dataset
  print(paste0('number of subjects: ', dim(y_c_strata)[1]))
  print(paste0('number of spikes: ', nrow(data_df4)))
  
  
  # 3.3) in case there are some subjects or features with no events
  

  patient_sel = unique(data_df4$subject_num) %>% sort()
  feature_sel = unique(data_df4$feature_id) %>% sort()
  
  return(list(data_df4 = data_df4,
              y_c_strata = y_c_strata,
              query_y_cs = query_y_cs,
              patient_sel = patient_sel,
              feature_sel = feature_sel,
              y_c_strata_sel = y_c_strata_sel
              ))
  
  
  
}


step_1_log_intensities <- function(dataset, data_df4, time_grid_est, time_grid, time_grid_both, full = T){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: obtain log-intensities (X) and intensities (Lambda)
  #
  # 
  # inputs:
  # 
  # - dataset           (dataset generated from simulation)
  # - data_df4          (data.table with `time`, `feature_id`, `subject_num`)
  # - time_grid_est     (m_est-dim vec of discretized times)
  # - time_grid         (m-dim vec)
  # - time_grid_both    (union of the two)
  # - full              (boolean)      are we including truths?
  #
  #
  # outputs:
  # 
  # - list of:
  #   - X_k_est                   (p x m_est   x n)
  #   - X_k_truth                 (p x m_truth x n)
  #   - X_k_coarse_truth          (p x m_est   x n)
  #   - X_k_both_truth            (p x m_both  x n)
  #   - Lambda_k_truth            (p x m_truth x n)
  #   - Lambda_k_coarse_truth     (p x m_est   x n)
  #   - Lambda_k_est              (p x m_est   x n)
  #
  # ----------------------------------------------------------------------------
  
  if(! full){
    X_k_est <- subject_specific_log_intensity(data_df4, time_grid_est)
    return(list(X_k_est = X_k_est))
  } else{

    X_k_est <- subject_specific_log_intensity(data_df4, time_grid_est)
    X_k_truth <- dataset$X_k_truth
    X_k_coarse_truth <- dataset$X_k_coarse_truth
    X_k_both_truth <- dataset$X_k_both_truth
    
    Lambda_k_truth        <- exp(X_k_truth)
    Lambda_k_coarse_truth <- exp(X_k_coarse_truth)
    Lambda_k_est          <- exp(X_k_est)
  
    
    return(list(X_k_est = X_k_est,
                X_k_truth = X_k_truth,
                X_k_coarse_truth = X_k_coarse_truth,
                X_k_both_truth = X_k_both_truth,
                Lambda_k_truth = Lambda_k_truth,
                Lambda_k_coarse_truth = Lambda_k_coarse_truth,
                Lambda_k_est = Lambda_k_est))
  }
}

step_2_rho_i <- function(dataset, data_df4, kernel_params, rho_kernel, patient_sel, feature_sel, 
                         time_grid, time_grid_est, query_y_c, y_c_strata, ncores, full = T){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: obtain rho_i estimates
  #
  # - we allow for rho_i^(y_c) or rho_i without y_c via `rho_kenrel`
  #
  # inputs:
  #
  # - dataset
  # - kernel_params   (list of parameters for various kernels)
  #     only need base_kernel params (that don't change over time)
  # 
  # - rho_kernel      (binary) do we utilize y_c here?
  # - patient_sel     (id's of patients that are included)
  # - feature_sel     (id's of features that are included)
  # - time_grid       (vector of timepoints)
  # - time_grid_est   (vector of timepoints)
  # - query_y_c       (q-c dim vector)
  # - y_c_strata      (n x q-c matrix)
  # - ncores          (integer)
  # - full            (boolean)    are we including truths in our estimation?
  #
  # outputs:
  # 
  # - list of:
  #   - rho_i_truth               (p x m)
  #   - rho_i_coarse_truth        (p x m_est)
  #   - rho_i_X_truth             (p x m)
  #   - rho_i_X_coarse_truth      (p x m_est)
  #   - rho_i_est                 (p x m_est)
  #   - rho_list                  (list format, [[1]] = rho_i, [[2]] = rho_ij, rho_ij may be rho_ii only)
  #   - weights                   (n-dim vector)
  #   - y_c_s                     (n-dim vector)
  #
  # ----------------------------------------------------------------------------
  

  
  # prep
  
  if(full){
    p <- dim(kernel_params$prec_mat_truth$adj_mat)[1]
  
    # 1) rho_truth
    # rho_i_truth = exp(mu(t) + 0.5 * diag(GP_cov))
    
    
    
    rho_i_truth_value   <- exp(kernel_params$base_kernel_params$base_GP_mean + 0.5 * kernel_params$base_kernel_params$base_variance)
    rho_i_truth         <- replicate(p, exp(kernel_params$base_mean + 0.5 * diag(kernel_params$base_cov))) %>% t() 
    rho_i_coarse_truth  <- replicate(p, exp(kernel_params$base_mean_est + 0.5 * diag(kernel_params$base_cov_est))) %>% t()   
  }
  
  
  # 3) rho_i_est from data
  if(rho_kernel){
    
    gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
    
    # get weights in preparation 
    weights <- apply(y_c_strata, 1, function(row) {
      step_6_kernel(as.numeric(row), query_y_c, gamma_c) 
    })    
    weights2 <- weights / sum(weights) # normalize
    
    if(full){
      

      # 2) rho_X_truth - weighted average of X_truths
      mat_list <- lapply(dataset$subject_data, function(x) exp(x$X_functions))
      rho_i_X_truth <- Reduce(`+`, Map(function(m, wt) m * wt, mat_list, weights2))
      
      mat_list_coarse <- lapply(dataset$subject_data, function(x) exp(x$X_functions_coarse))
      rho_i_X_coarse_truth <- Reduce(`+`, Map(function(m, wt) m * wt, mat_list_coarse, weights2))
    }

    
    
    # i_neq_j is false because we want G_ij
    rho_list <- estimate_intensities_stratum_parallel_with_yc(data_df4, query_y_c, y_c_strata, 
                                                              patient_sel, feature_sel, 
                                                              time_grid_est, i_neq_j = T, ncores)
    
  } else{
    

    if(full){
      
      # 2) sample mean of the X_functions, beacuse there is no weight
      mat_list <- lapply(dataset$subject_data, function(x) exp(x$X_functions))
      weights2 <- rep(1/length(mat_list), length(mat_list))
      
      rho_i_X_truth <- Reduce(`+`, Map(function(m, wt) m * wt, mat_list, weights2))
      
      mat_list_coarse <- lapply(dataset$subject_data, function(x) exp(x$X_functions_coarse))
      rho_i_X_coarse_truth <- Reduce(`+`, Map(function(m, wt) m * wt, mat_list_coarse, weights2))
    }

    
    # i_neq_j is false because we don't care about G_ij
    rho_list <- estimate_intensities_stratum_parallel_v4(data_df4, patient_sel, feature_sel, time_grid_est, F, ncores)
  }
  
  rho_i_est <- do.call(rbind, lapply(rho_list[[1]], function(v) as.numeric(v))) # p x m matrix of estimated mean intensities     

  
  # choose what to output depending on if we have truths
  if(full){
    result_list <- list(rho_i_truth = rho_i_truth,
                        rho_i_coarse_truth = rho_i_coarse_truth,
                        rho_i_X_truth = rho_i_X_truth,
                        rho_i_X_coarse_truth = rho_i_X_coarse_truth,
                        rho_i_est = rho_i_est,
                        rho_list = rho_list,
                        weights = weights2)
  } else{
    result_list <- list(rho_i_est = rho_i_est,
                        rho_list = rho_list,
                        weights = weights2)
  }

  
  if(rho_kernel){
    result_list$y_c_s <- y_c_strata
  } 
  
  return(result_list)
}

step_2_rho_ij <- function(step_1, step_2, kernel_params, i_neq_j, full = T){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate rho_ij
  #
  #
  # inputs:
  #
  # - step_1
  #   - [[5]] Lambda_k_truth            (p x m_truth x n)
  #   - [[6]] Lambda_k_coarse_truth     (p x m_est   x n)
  # 
  # - step_2
  #   - rho_i_truth               (p x m)
  #   - rho_i_coarse_truth        (p x m_est)
  #   - rho_i_X_truth             (p x m)
  #   - rho_i_X_coarse_truth      (p x m_est)
  #   - rho_i_est                 (p x m_est)
  #   - rho_list                  (list format, [[1]] = rho_i, [[2]] = rho_ij)
  #
  # - kernel_params              (list of kernel params)
  # - i_neq_j                    (boolean) if true, calculate i =/= j
  # - full                       (boolean)    are we including truths in our estimation?
  #
  # outputs:
  # 
  # - list of:
  #   - rho_ii_truth               (list of m x m matrices)
  #   - rho_ii_coarse_truth        (list of m_est x m_est matrices)
  #   - rho_ii_X_truth             (list of m x m matrices)
  #   - rho_ii_X_coarse_truth      (list of m_est x m_est matrices)
  #   - rho_ii_est                 (list of m_est x m_est matrices)
  #
  # ----------------------------------------------------------------------------
  

  if(! full){
    rho_list  <- step_2$rho_list
    rho_ii_est <- rho_list[[2]]
    
    return(list(rho_ii_est = rho_ii_est))
  }
  
  # prep
  rho_i_truth          <- step_2[[1]]
  rho_i_coarse_truth   <- step_2[[2]]
  rho_i_X_truth        <- step_2[[3]]
  rho_i_X_coarse_truth <- step_2[[4]]
  rho_i_est            <- step_2[[5]]
  rho_list             <- step_2[[6]]
  
  Lambda_k_truth        <- step_1[[5]]
  Lambda_k_coarse_truth <- step_1[[6]]
  
  p <- dim(rho_i_truth)[1]
  m <- dim(rho_i_truth)[2]
  m_est <- dim(rho_i_coarse_truth)[2]
  n <- dim(Lambda_k_truth)[3]

  # est
  rho_ii_est <- rho_list[[2]]
  rho_ii_truth <- list()
  rho_ii_coarse_truth <- list()
  rho_ii_X_truth <- list()
  rho_ii_X_coarse_truth <- list()
  
  # create keys depending on if we want i =/= j
  if(i_neq_j){
    # i_j keys
    keys <- apply(which(upper.tri(matrix(1, p, p), diag = TRUE), arr.ind = TRUE), 1, 
                  function(x) paste0(x[1], "_", x[2]))
  } else{
    # i_i keys
    keys <- paste0(1:p, '_', 1:p)
  }
  
  
  for(key in keys){
    ij <- as.integer(strsplit(key, "_")[[1]])
    i <- ij[1]
    j <- ij[2]
    
    # theory
    
    GP_simu_var <- kronecker(kernel_params$prec_mat_truth$cor_mat, kernel_params$base_cov) # pm x pm matrix
    GP_simu_var_est <- kronecker(kernel_params$prec_mat_truth$cor_mat, kernel_params$base_cov_est)
    
    rho_ii_truth[[key]]        <- tcrossprod(rho_i_truth[i,],        rho_i_truth[j,])        * exp(extract_block_structure_ij(GP_simu_var, m, i, j))
    rho_ii_coarse_truth[[key]] <- tcrossprod(rho_i_coarse_truth[i,], rho_i_coarse_truth[j,]) * exp(extract_block_structure_ij(GP_simu_var_est, m_est, i, j))
    
    # X_k truths
    mats <- array(0, dim = c(m, m))
    coarse_mats <- array(0, dim = c(m_est, m_est))  
    
    for (k in 1:n) {
      vi <- Lambda_k_truth[i, , k]
      vj <- Lambda_k_truth[j, , k]
      mats <- mats + tcrossprod(vi, vj)  # faster than outer(v, v)
      
      vi <- Lambda_k_coarse_truth[i, , k]
      vj <- Lambda_k_coarse_truth[j, , k]
      coarse_mats <- coarse_mats + tcrossprod(vi, vj)  # faster than outer(v, v) 
      
      # Average across replicates
      rho_ii_X_truth[[key]] <- mats / n
      rho_ii_X_coarse_truth[[key]] <- coarse_mats / n      
    }
    
  }
  

  return(list(rho_ii_truth = rho_ii_truth,
              rho_ii_coarse_truth = rho_ii_coarse_truth,
              rho_ii_X_truth = rho_ii_X_truth,
              rho_ii_X_coarse_truth = rho_ii_X_coarse_truth,
              rho_ii_est = rho_ii_est))
}

step_3_g_ij <- function(step_2, step_2b, kernel_params, i_neq_j, full = T){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: covariance function (g_ij) estimation
  #
  # 
  # inputs:
  #
  # - step_2:
  #   - rho_i_truth               (p x m)
  #   - rho_i_coarse_truth        (p x m_est)
  #   - rho_i_X_truth             (p x m)
  #   - rho_i_X_coarse_truth      (p x m_est)
  #   - rho_i_est                 (p x m_est)
  #   - rho_list                  (list format, [[1]] = rho_i, [[2]] = rho_ij)
  #
  # - step_2b: 
  #   - rho_ii_truth               (list of m x m matrices, could have only i_i or i_j as well)
  #   - rho_ii_coarse_truth        (list of m_est x m_est matrices)
  #   - rho_ii_X_truth             (list of m x m matrices)
  #   - rho_ii_X_coarse_truth      (list of m_est x m_est matrices)
  #   - rho_ii_est                 (list of m_est x m_est matrices)
  #
  # - kernel_params              (list of kernel params)
  # - i_neq_j                    (boolean) do we include i =/= j terms?
  # - full                       (boolean)    are we including truths in our estimation?
  # 
  # outputs:
  #
  # - list of:
  #   - g_ij_ground_truth                (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_ground_truth         (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_truth                       (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_truth                (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_X_truth                     (list of m x m matrices for i_j entries)
  #   - g_ij_X_coarse_truth              (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_est                         (list of m_est x m_est matrices for i_j entries)  
  #
  # ----------------------------------------------------------------------------
  
  if(! full){ # just estimation
    rho_i_est <- step_2$rho_i_est
    rho_ii_est <- step_2b$rho_ii_est
    
    if(i_neq_j){
      g_ij_est <- estimate_covariance_functions_ij(rho_i_est, rho_ii_est)
    } else{
      g_ij_est <- estimate_covariance_functions_ii(rho_i_est, rho_ii_est)
    }
    
    return(list(g_ij_est = g_ij_est))
  }
  
  
  # prep
  rho_i_truth          <- step_2[[1]]
  rho_i_coarse_truth   <- step_2[[2]]
  rho_i_X_truth        <- step_2[[3]]
  rho_i_X_coarse_truth <- step_2[[4]]
  rho_i_est            <- step_2[[5]]
  
  rho_ii_truth          <- step_2b[[1]]
  rho_ii_coarse_truth   <- step_2b[[2]]
  rho_ii_X_truth        <- step_2b[[3]]
  rho_ii_X_coarse_truth <- step_2b[[4]]
  rho_ii_est            <- step_2b[[5]]
  
  p   <- dim(rho_i_truth)[1]
  m   <- dim(rho_i_truth)[2]
  m_est <- dim(rho_i_est)[2]
  
  # est
  
  if(i_neq_j){
    
    g_ij_truth          <- estimate_covariance_functions_ij(rho_i_truth,          rho_ii_truth)
    g_ij_coarse_truth   <- estimate_covariance_functions_ij(rho_i_coarse_truth,   rho_ii_coarse_truth)
    g_ij_X_truth        <- estimate_covariance_functions_ij(rho_i_X_truth,        rho_ii_X_truth)
    g_ij_X_coarse_truth <- estimate_covariance_functions_ij(rho_i_X_coarse_truth, rho_ii_X_coarse_truth)
    g_ij_est            <- estimate_covariance_functions_ij(rho_i_est,            rho_ii_est)
  } else{
    
    
    g_ij_truth          <- estimate_covariance_functions_ii(rho_i_truth,          rho_ii_truth)
    g_ij_coarse_truth   <- estimate_covariance_functions_ii(rho_i_coarse_truth,   rho_ii_coarse_truth)
    g_ij_X_truth        <- estimate_covariance_functions_ii(rho_i_X_truth,        rho_ii_X_truth)
    g_ij_X_coarse_truth <- estimate_covariance_functions_ii(rho_i_X_coarse_truth, rho_ii_X_coarse_truth)
    g_ij_est            <- estimate_covariance_functions_ii(rho_i_est,            rho_ii_est)
    
  }
  
  
  g_ij_ground_truth_pm <- kronecker(kernel_params$prec_mat_truth$cor_mat, kernel_params$base_cov)
  g_ij_coarse_ground_truth_pm <- kronecker(kernel_params$prec_mat_truth$cor_mat, kernel_params$base_cov_est)  
  
  g_ij_ground_truth        <- extract_block_structure_v2(g_ij_ground_truth_pm, p, m)
  g_ij_coarse_ground_truth <- extract_block_structure_v2(g_ij_coarse_ground_truth_pm, p, m_est)
  
  
  return(list(g_ij_ground_truth = g_ij_ground_truth,
              g_ij_coarse_ground_truth = g_ij_coarse_ground_truth,
              g_ij_truth = g_ij_truth,
              g_ij_coarse_truth = g_ij_coarse_truth,
              g_ij_X_truth = g_ij_X_truth,
              g_ij_X_coarse_truth= g_ij_X_coarse_truth,
              g_ij_est = g_ij_est))
}

step_4_eigendecomp <- function(step_3, p, time_grid, time_grid_est, full = T){
  

  # ----------------------------------------------------------------------------
  # 
  # GOAL: eigendecomposition of G_ii
  #
  # inputs:
  #
  # - step_3
  #   - g_ij_ground_truth                (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_ground_truth         (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_truth                       (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_truth                (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_X_truth                     (list of m x m matrices for i_j entries)
  #   - g_ij_X_coarse_truth              (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_est                         (list of m_est x m_est matrices for i_j entries)   
  #
  # - p                   (scalar)
  # - time_grid           (m-dim vector)
  # - time_grid_est       (m_est-dim vector)
  # - full                (boolean)    are we including truths in our estimation?
  #
  #
  # outputs:
  #
  # - list of:
  #   - eigen_decomp_truth             (list of 4 things)
  #     - [[1]] eigenvalues            (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors           (list of p matrices of m x d_i)
  #     - [[3]] n_dims                 (list of p integers denoting d_i)
  #
  #   - eigen_decomp_coarse_truth
  #   - eigen_decomp_X_truth
  #   - eigen_decomp_X_coarse_truth
  #   - eigen_decomp_est
  #
  # ----------------------------------------------------------------------------
  
  if(! full){
    g_ij_est <- step_3$g_ij_est
    g_ii_est <- prep_eigendecomposition_ii(g_ij_est, p)
    
    # perform eigendecomposition
    eigen_decomp_est <- compute_eigendecomposition_ii(g_ii_est)
    
    return(list(eigen_decomp_est = eigen_decomp_est))
  }

  # prep
  g_ij_truth          <- step_3[[3]]
  g_ij_coarse_truth   <- step_3[[4]]
  g_ij_X_truth        <- step_3[[5]]
  g_ij_X_coarse_truth <- step_3[[6]]
  g_ij_est            <- step_3[[7]]
  

  m     <- dim(g_ij_truth[[1]])[1]
  m_est <- dim(g_ij_truth[[2]])[1]

    
  
  # only keep g_ii

  g_ii_truth          <- prep_eigendecomposition_ii(g_ij_truth, p)
  g_ii_coarse_truth   <- prep_eigendecomposition_ii(g_ij_coarse_truth, p)
  g_ii_X_truth        <- prep_eigendecomposition_ii(g_ij_X_truth, p)
  g_ii_X_coarse_truth <- prep_eigendecomposition_ii(g_ij_X_coarse_truth, p)
  g_ii_est            <- prep_eigendecomposition_ii(g_ij_est, p)

  

  
  # perform eigendecomposition
  eigen_decomp_truth          <- compute_eigendecomposition_ii(g_ii_truth)
  eigen_decomp_coarse_truth   <- compute_eigendecomposition_ii(g_ii_coarse_truth)
  eigen_decomp_X_truth        <- compute_eigendecomposition_ii(g_ii_X_truth)
  eigen_decomp_X_coarse_truth <- compute_eigendecomposition_ii(g_ii_X_coarse_truth)
  eigen_decomp_est            <- compute_eigendecomposition_ii(g_ii_est)
  
  
  
  
  return(list(eigen_decomp_truth = eigen_decomp_truth,
              eigen_decomp_coarse_truth = eigen_decomp_coarse_truth,
              eigen_decomp_X_truth = eigen_decomp_X_truth,
              eigen_decomp_X_coarse_truth = eigen_decomp_X_coarse_truth,
              eigen_decomp_est = eigen_decomp_est))
  
  
  
}

step_5_KL_expansion <- function(step_1, step_4, kernel_params, time_grid, time_grid_est, ncores){

  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate KL coefficients
  #
  # inputs:
  #
  # - step_1
  #   - X_k_est                   (p x m_est   x n)
  #   - X_k_truth                 (p x m_truth x n)
  #   - X_k_coarse_truth          (p x m_est   x n)
  #   - X_k_both_truth            (p x m_both  x n)
  #   - Lambda_k_truth            (p x m_truth x n)
  #   - Lambda_k_coarse_truth     (p x m_est   x n)
  #
  # - step_4
  #   - eigen_decomp_truth               (list of 3 things)
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  #
  #   - eigen_decomp_coarse_truth
  #   - eigen_decomp_X_truth
  #   - eigen_decomp_X_coarse_truth
  #   - eigen_decomp_est
  #
  # - kernel_params    (list)
  # - time_grid        (m-dim vec)
  # - time_grid_est    (m_est-dim vec)
  # - ncores           (integer)
  #
  #
  # outputs:
  # - list of:
  #   - kl_coeffs_truth             (n x p x d_max matrix)
  #   - kl_coeffs_coarse_truth      (n x p x d_max matrix)
  #   - kl_coeffs_X_truth           (n x p x d_max matrix)
  #   - kl_coeffs_X_coarse_truth    (n x p x d_max matrix)
  #   - kl_coeffs_est               (n x p x d_max matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  # prep
  
  X_k_est          <- step_1[[1]]
  X_k_truth        <- step_1[[2]]
  X_k_coarse_truth <- step_1[[3]]
  
  
  eigen_decomp_truth          <- step_4[[1]]
  eigen_decomp_coarse_truth   <- step_4[[2]]
  eigen_decomp_X_truth        <- step_4[[3]]
  eigen_decomp_X_coarse_truth <- step_4[[4]]
  eigen_decomp_est            <- step_4[[5]]
  
  # 5.1) center X_k_est and X_k_truth
  

  # - X_k_est_center                (p x m_est x n)
  # - X_k_truth_center              (p x m x n)       X_k_truth - acutal mean
  # - X_k_truth_center_est          (p x m x n)       X_k_truth - empirical mean
  # - X_k_coarse_truth_center       (p x m_est x n)   X_k_coarse_truth - actual mean
  # - X_k_coarse_truth_center_est   (p x m_est x n)   X_k_coarse_truth - empirical mean
  
  mean_mat_est <- apply(X_k_est, c(1, 2), mean)  # result is p x m matrix
  X_k_est_center <- sweep(X_k_est, c(1, 2), mean_mat_est, FUN = "-")
  
  
  mean_mat_truth <- apply(X_k_truth, c(1, 2), mean)  # result is p x m matrix
  X_k_truth_center_est <- sweep(X_k_truth, c(1, 2), mean_mat_truth, FUN = "-") 
  
  X_k_truth_center <- sweep(X_k_truth, 2, kernel_params$base_mean, "-")
  #X_k_truth_center <- X_k_truth - kernel_params$base_mean
  
  
  mean_mat_coarse_truth <- apply(X_k_coarse_truth, c(1, 2), mean)  # result is p x m matrix
  X_k_coarse_truth_center_est <- sweep(X_k_coarse_truth, c(1, 2), mean_mat_coarse_truth, FUN = "-")   
  #X_k_coarse_truth_center <- X_k_coarse_truth - kernel_params$base_GP_mean
  X_k_coarse_truth_center <- sweep(X_k_coarse_truth, 2, kernel_params$base_mean_est, "-")
  
  
  # 5.2) get KL coeffs  
  # 
  # - kl_coeffs_v2                     (n x p x d)   coeffs from X_k_est
  # - kl_coeffs_X_coarse_truth_v2      (n x p x d)   coeffs from X_k_coarse_truth
  # - kl_coeffs_X_truth_v2             (n x p x d)   coeffs from X_k_truth
  # - kl_coeffs_coarse_truth_v2        (n x p x d)   coeffs from X_k_coarse_truth but eigendecomp from theory
  # - kl_coeffs_truth_v2               (n x p x d)   coeffs from X_k_truth        but eigendecomp from theory
  
  kl_coeffs_truth           <- estimate_kl_coefficients_parallel_v2(X_k_truth_center_est,         eigen_decomp_truth$eigenfunctions,          time_grid,      ncores)
  kl_coeffs_coarse_truth    <- estimate_kl_coefficients_parallel_v2(X_k_coarse_truth_center_est,  eigen_decomp_coarse_truth$eigenfunctions,   time_grid_est,  ncores)
  kl_coeffs_X_truth         <- estimate_kl_coefficients_parallel_v2(X_k_truth_center_est,         eigen_decomp_X_truth$eigenfunctions,        time_grid,      ncores)  
  kl_coeffs_X_coarse_truth  <- estimate_kl_coefficients_parallel_v2(X_k_coarse_truth_center_est,  eigen_decomp_X_coarse_truth$eigenfunctions, time_grid_est,  ncores)
  kl_coeffs_est             <- estimate_kl_coefficients_parallel_v2(X_k_est_center,               eigen_decomp_est$eigenfunctions,            time_grid_est,  ncores)

  
  
  return(list(
    kl_coeffs_truth = kl_coeffs_truth,
    kl_coeffs_coarse_truth = kl_coeffs_coarse_truth,
    kl_coeffs_X_truth = kl_coeffs_X_truth,
    kl_coeffs_X_coarse_truth = kl_coeffs_X_coarse_truth,
    kl_coeffs_est = kl_coeffs_est))
  
}

step_5_KL_covariance <- function(step_3, step_4, full = T){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: estimate covariances of the KL coefficients for CPGM method
  #
  # inputs:
  #
  # - step_3
  #   - g_ij_ground_truth                (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_ground_truth         (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_truth                       (list of m x m matrices for i_j entries)
  #   - g_ij_coarse_truth                (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_X_truth                     (list of m x m matrices for i_j entries)
  #   - g_ij_X_coarse_truth              (list of m_est x m_est matrices for i_j entries)
  #   - g_ij_est                         (list of m_est x m_est matrices for i_j entries)  
  # 
  # - step_4
  #   - eigen_decomp_truth               (list of 3 things)
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  #
  #   - eigen_decomp_coarse_truth
  #   - eigen_decomp_X_truth
  #   - eigen_decomp_X_coarse_truth
  #   - eigen_decomp_est
  #
  # - full         (boolean)    are we including truths in our estimation?
  #
  #
  # outputs:
  #
  # - list of:
  #   - KL_cov_truth           (list of d_i x d_j matrices for i_j entries)
  #   - KL_cov_coarse_truth    (list of d_i x d_j matrices for i_j entries)
  #   - KL_cov_X_truth         (list of d_i x d_j matrices for i_j entries)
  #   - KL_cov_X_coarse_truth  (list of d_i x d_j matrices for i_j entries)
  #   - KL_cov_est             (list of d_i x d_j matrices for i_j entries)
  #
  # ----------------------------------------------------------------------------
  
  if(! full){
    KL_cov_est <- estimate_KL_covariance(step_3$g_ij_est, step_4$eigen_decomp_est$eigenfunctions)
    
    return(list(KL_cov_est = KL_cov_est))
  }
  
  
  KL_cov_truth          <- estimate_KL_covariance(step_3[[3]], step_4[[1]]$eigenfunctions)
  KL_cov_coarse_truth   <- estimate_KL_covariance(step_3[[4]], step_4[[2]]$eigenfunctions)
  KL_cov_X_truth        <- estimate_KL_covariance(step_3[[5]], step_4[[3]]$eigenfunctions)
  KL_cov_X_coarse_truth <- estimate_KL_covariance(step_3[[6]], step_4[[4]]$eigenfunctions)
  KL_cov_est            <- estimate_KL_covariance(step_3[[7]], step_4[[5]]$eigenfunctions)
  
  return(list(KL_cov_truth = KL_cov_truth,
              KL_cov_coarse_truth = KL_cov_coarse_truth,
              KL_cov_X_truth = KL_cov_X_truth,
              KL_cov_X_coarse_truth = KL_cov_X_coarse_truth,
              KL_cov_est = KL_cov_est))
}

step_5b_KL_correlation <- function(step_5, p, full = T){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: estimate correlations of the KL coefficients for CPGM method
  #
  # inputs:
  #
  # - step_5
  #   - KL_cov_truth           (list of d_i x d_j matrices for i_j entries)
  #   - KL_cov_coarse_truth    (list of d_i x d_j matrices for i_j entries)
  #   - KL_cov_X_truth         (list of d_i x d_j matrices for i_j entries)
  #   - KL_cov_X_coarse_truth  (list of d_i x d_j matrices for i_j entries)
  #   - KL_cov_est             (list of d_i x d_j matrices for i_j entries)
  # 
  # - p            (integer)
  # - full         (boolean)    are we including truths in our estimation?
  #
  #
  # outputs:
  #
  # - list of:
  #   - KL_cor_truth           (list of d_i x d_j matrices for i_j entries)
  #   - KL_cor_coarse_truth    (list of d_i x d_j matrices for i_j entries)
  #   - KL_cor_X_truth         (list of d_i x d_j matrices for i_j entries)
  #   - KL_cor_X_coarse_truth  (list of d_i x d_j matrices for i_j entries)
  #   - KL_cor_est             (list of d_i x d_j matrices for i_j entries)
  #
  # ----------------------------------------------------------------------------
  
  if(! full){
    KL_cor_est <- estimate_KL_correlation(step_5$KL_cov_est, p)
    
    return(list(KL_cor_est = KL_cor_est))
  }
  
  
  KL_cor_truth          <- estimate_KL_correlation(step_5[[1]], p)
  KL_cor_coarse_truth   <- estimate_KL_correlation(step_5[[2]], p)
  KL_cor_X_truth        <- estimate_KL_correlation(step_5[[3]], p)
  KL_cor_X_coarse_truth <- estimate_KL_correlation(step_5[[4]], p)
  KL_cor_est            <- estimate_KL_correlation(step_5[[5]], p)
  
  return(list(KL_cor_truth = KL_cor_truth,
              KL_cor_coarse_truth = KL_cor_coarse_truth,
              KL_cor_X_truth = KL_cor_X_truth,
              KL_cor_X_coarse_truth = KL_cor_X_coarse_truth,
              KL_cor_est = KL_cor_est))
}

steps_78_OG <- function(kl_coeffs, y_c_strata, query_y_c, eigenfunctions, ncores){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: perform steps 7 and 8 for the OG method
  #
  # - v3, v3, v2
  # 
  # 
  
  
  p <- dim(kl_coeffs)[2]
  
  # step 6
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  K_c <- construct_kernel_matrix_step_6(y_c_strata, query_y_c, gamma_c)
  
  # steps 7 and 8
  V_YcXij <- construct_cross_covariance_matrix_v3(kl_coeffs, ncores)
  M_hat <- estimate_regression_operators_v3(K_c, V_YcXij, p)
  V_cond <- evaluate_regression_at_query_v2(M_hat, y_c_strata, query_y_c, eigenfunctions, gamma_c, p)  
  
  return(list(V_cond = V_cond,
              K_c = K_c))
  
}

steps_78_JASA <- function(kl_coeffs, y_c_strata, query_y_c, eigenfunctions, ncores) {
  
  # ----------------------------------------------------------------------------
  #
  #
  # 9/9/2025
  #
  # GOAL: obtain conditional covariance immediately, with the JASA routine
  #
  # v4, v3, v3
  #
  # input:
  # 
  # - kl_coeffs             (n_stratum x p x d matrix) tensor of KL coefficients
  # - y_c_strata            (n_stratum x q_c matrix)   matrix of continuous covariates
  # - query_y_c             (q_c dim vector)           vector of query covariates
  # - eigenfunctions        (list of p matrices, each of which is m x d_i)
  # - ncores                (integer)
  #
  #
  # output:
  #
  # - V_cond                (list of length p^2, each element is m x m matrix)
  #
  #
  # ----------------------------------------------------------------------------
  
  p <- dim(kl_coeffs)[2]
  
  # step 6
  gamma_c <- select_gamma_c_bandwidth_v2(y_c_strata)
  K_c <- construct_kernel_matrix_step_6(y_c_strata, query_y_c, gamma_c)
  
  V_YcXij <- construct_cross_covariance_matrix_v4(kl_coeffs, y_c_strata, query_y_c, gamma_c, ncores)
  M_hat <- estimate_regression_operators_v3(K_c, V_YcXij, p)
  V_cond <- evaluate_regression_at_query_v3(M_hat, y_c_strata, query_y_c, eigenfunctions, gamma_c, p)  
  
  return(list(V_cond = V_cond,
              K_c = K_c))
}

steps_78 <- function(step_4, step_5, kernel_params, y_c_strata, query_y_c, method, ncores){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: perform steps 7 and 8 (cross covariance and covariance) for OG and JASA methods
  #
  # input:
  # 
  # - step_4
  #   - eigen_decomp_truth               (list of 3 things)
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  #
  #   - eigen_decomp_coarse_truth
  #   - eigen_decomp_X_truth
  #   - eigen_decomp_X_coarse_truth
  #   - eigen_decomp_est
  #
  # - step_5
  #   - kl_coeffs_truth             (n x p x d_max matrix)
  #   - kl_coeffs_coarse_truth      (n x p x d_max matrix)
  #   - kl_coeffs_X_truth           (n x p x d_max matrix)
  #   - kl_coeffs_X_coarse_truth    (n x p x d_max matrix)
  #   - kl_coeffs_est               (n x p x d_max matrix)
  #
  # - kernel_params  (list)
  # - y_c_strata     (n x q_c matrix)  matrix of covariates
  # - query_y_c      (q_c dim vector)  queried vector
  # - method         (string)          OG or JASA
  # - ncores         (integer)         number of cores
  #
  # 
  # outputs
  # 
  # - list of: 
  #   - V_cond_ground_truth_full          (pm x pm matrix)
  #   - V_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - V_cond_truth_full                 (pm x pm matrix)
  #   - V_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - V_cond_X_truth_full               (pm x pm matrix)
  #   - V_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - V_cond_est_full                   (pm_est x pm_est matrix)
  # 
  # ----------------------------------------------------------------------------
  
  if (!method %in% c('OG', 'JASA')) {
    stop("method must be OG or JASA, stop at steps_78")
  }
  
  # prep:
  
  kl_coeffs_truth          <- step_5[[1]]
  kl_coeffs_coarse_truth   <- step_5[[2]]
  kl_coeffs_X_truth        <- step_5[[3]]
  kl_coeffs_X_coarse_truth <- step_5[[4]]
  kl_coeffs_est            <- step_5[[5]]
  
  eigen_decomp_truth          <- step_4[[1]]
  eigen_decomp_coarse_truth   <- step_4[[2]]
  eigen_decomp_X_truth        <- step_4[[3]]
  eigen_decomp_X_coarse_truth <- step_4[[4]]
  eigen_decomp_est            <- step_4[[5]]
  
  p     <- length(eigen_decomp_truth[[2]])
  m     <- dim(eigen_decomp_truth[[2]][[1]])[1]
  m_est <- dim(eigen_decomp_coarse_truth[[2]][[1]])[1]
  
  # estimate

  V_cond_ground_truth_full        <- kronecker(kernel_params$prec_mat_truth$cor_mat, kernel_params$base_cov)
  V_cond_coarse_ground_truth_full <- kronecker(kernel_params$prec_mat_truth$cor_mat, kernel_params$base_cov_est)
  V_cond_ground_truth             <- extract_block_structure_v2(V_cond_ground_truth_full, p, m)
  V_cond_coarse_ground_truth      <- extract_block_structure_v2(V_cond_coarse_ground_truth_full, p, m_est)
  

  if(method == 'OG'){
    V_cond_truth          <- steps_78_OG(kl_coeffs_truth,           y_c_strata, query_y_c, eigen_decomp_truth$eigenfunctions,          ncores)[[1]]
    V_cond_coarse_truth   <- steps_78_OG(kl_coeffs_coarse_truth,    y_c_strata, query_y_c, eigen_decomp_coarse_truth$eigenfunctions,   ncores)[[1]]
    V_cond_X_truth        <- steps_78_OG(kl_coeffs_X_truth,         y_c_strata, query_y_c, eigen_decomp_X_truth$eigenfunctions,        ncores)[[1]]
    V_cond_X_coarse_truth <- steps_78_OG(kl_coeffs_X_coarse_truth,  y_c_strata, query_y_c, eigen_decomp_X_coarse_truth$eigenfunctions, ncores)[[1]]
    V_cond_est            <- steps_78_OG(kl_coeffs_est,             y_c_strata, query_y_c, eigen_decomp_est$eigenfunctions,            ncores)[[1]]
    K_c                   <- steps_78_OG(kl_coeffs_est,             y_c_strata, query_y_c, eigen_decomp_est$eigenfunctions,            ncores)[[2]]
  } else if(method == 'JASA'){
    V_cond_truth          <- steps_78_JASA(kl_coeffs_truth,           y_c_strata, query_y_c, eigen_decomp_truth$eigenfunctions,          ncores)[[1]]
    V_cond_coarse_truth   <- steps_78_JASA(kl_coeffs_coarse_truth,    y_c_strata, query_y_c, eigen_decomp_coarse_truth$eigenfunctions,   ncores)[[1]]
    V_cond_X_truth        <- steps_78_JASA(kl_coeffs_X_truth,         y_c_strata, query_y_c, eigen_decomp_X_truth$eigenfunctions,        ncores)[[1]]
    V_cond_X_coarse_truth <- steps_78_JASA(kl_coeffs_X_coarse_truth,  y_c_strata, query_y_c, eigen_decomp_X_coarse_truth$eigenfunctions, ncores)[[1]]
    V_cond_est            <- steps_78_JASA(kl_coeffs_est,             y_c_strata, query_y_c, eigen_decomp_est$eigenfunctions,            ncores)[[1]]
    K_c                   <- steps_78_JASA(kl_coeffs_est,             y_c_strata, query_y_c, eigen_decomp_est$eigenfunctions,            ncores)[[2]]
  } else{
    stop("method must be OG or JASA, stop at steps_78")
  }
  

  # construct entire pm x pm matrix
  
  V_cond_est_full            <- assemble_block_matrix_v2(V_cond_est,            p, m_est)
  V_cond_X_coarse_truth_full <- assemble_block_matrix_v2(V_cond_X_coarse_truth, p, m_est)
  V_cond_X_truth_full        <- assemble_block_matrix_v2(V_cond_X_truth,        p, m)    
  V_cond_coarse_truth_full   <- assemble_block_matrix_v2(V_cond_coarse_truth,   p, m_est)
  V_cond_truth_full          <- assemble_block_matrix_v2(V_cond_truth,          p, m)
  

  return(list(
    V_cond_ground_truth_full = V_cond_ground_truth_full,
    V_cond_coarse_ground_truth_full = V_cond_coarse_ground_truth_full,
    V_cond_truth_full = V_cond_truth_full,
    V_cond_coarse_truth_full = V_cond_coarse_truth_full,
    V_cond_X_truth_full = V_cond_X_truth_full,
    V_cond_X_coarse_truth_full = V_cond_X_coarse_truth_full,
    V_cond_est_full = V_cond_est_full,
    K_c = K_c
  ))
  
}

step_9_C_cond_from_V_cond <- function(step_8, kernel_params_i){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate C_cond from V_cond
  #
  #
  # input:
  #
  # - step_8:
  #   - V_cond_ground_truth_full          (pm x pm matrix)
  #   - V_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - V_cond_truth_full                 (pm x pm matrix)
  #   - V_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - V_cond_X_truth_full               (pm x pm matrix)
  #   - V_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - V_cond_est_full                   (pm_est x pm_est matrix)
  #
  # - kernel_params_i    (list)
  #
  #
  # outputs:
  # 
  # - list of:
  #   - C_cond_ground_truth_full          (pm x pm matrix)
  #   - C_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - C_cond_truth_full                 (pm x pm matrix)
  #   - C_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - C_cond_X_truth_full               (pm x pm matrix)
  #   - C_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - C_cond_est_full                   (pm_est x pm_est matrix)  
  #
  # ----------------------------------------------------------------------------
  
  # prep
  
  V_cond_ground_truth_full        <- step_8[[1]]
  V_cond_coarse_ground_truth_full <- step_8[[2]]
  V_cond_truth_full               <- step_8[[3]]
  V_cond_coarse_truth_full        <- step_8[[4]]
  V_cond_X_truth_full             <- step_8[[5]]
  V_cond_X_coarse_truth_full      <- step_8[[6]]
  V_cond_est_full                 <- step_8[[7]]
  
  p <- dim(kernel_params_i$prec_mat_truth$adj_mat)[1]
  m <- dim(V_cond_ground_truth_full)[1] / p
  m_est <- dim(V_cond_coarse_ground_truth_full)[1] / p
  
  # ground truth is cor_mat \otimes I_m
  C_cond_ground_truth_full          <- kronecker(kernel_params_i$prec_mat_truth$cor_mat, diag(m))
  C_cond_coarse_ground_truth_full   <- kronecker(kernel_params_i$prec_mat_truth$cor_mat, diag(m_est))
  C_cond_ground_truth               <- extract_block_structure_v2(C_cond_ground_truth_full, p, m)
  C_cond_coarse_ground_truth        <- extract_block_structure_v2(C_cond_coarse_ground_truth_full, p, m_est)
  
  # v3: C_ii = indentity
  # v4: C_ii calculated like C_ij
  C_cond_truth          <- estimate_conditional_correlation_v4(V_cond_truth_full, p)
  C_cond_coarse_truth   <- estimate_conditional_correlation_v4(V_cond_coarse_truth_full, p)
  C_cond_X_truth        <- estimate_conditional_correlation_v4(V_cond_X_truth_full, p) 
  C_cond_X_coarse_truth <- estimate_conditional_correlation_v4(V_cond_X_coarse_truth_full, p)
  C_cond_est            <- estimate_conditional_correlation_v4(V_cond_est_full, p)
  

  
  # construct entire pm x pm block 
  
  C_cond_est_full            <- assemble_block_matrix_v2(C_cond_est,            p, m_est)
  C_cond_X_coarse_truth_full <- assemble_block_matrix_v2(C_cond_X_coarse_truth, p, m_est)
  C_cond_X_truth_full        <- assemble_block_matrix_v2(C_cond_X_truth,        p, m)    
  C_cond_coarse_truth_full   <- assemble_block_matrix_v2(C_cond_coarse_truth,   p, m_est)
  C_cond_truth_full          <- assemble_block_matrix_v2(C_cond_truth,          p, m)
  

  return(list(
    C_cond_ground_truth_full = C_cond_ground_truth_full,
    C_cond_coarse_ground_truth_full = C_cond_coarse_ground_truth_full,
    C_cond_truth_full = C_cond_truth_full,
    C_cond_coarse_truth_full = C_cond_coarse_truth_full,
    C_cond_X_truth_full = C_cond_X_truth_full,
    C_cond_X_coarse_truth_full = C_cond_X_coarse_truth_full,
    C_cond_est_full = C_cond_est_full))
  
  
}

step_9_C_cond_from_KL_cor <- function(step_4, step_5b, kernel_params, full = T){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: construct correlation matrix from CPGM method USING cor() instead of cov
  #
  #
  # inputs:
  #
  # - step_4
  #   - eigen_decomp_truth               (list of 3 things)
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)          THEY ARE NOT NORMALIZED
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  #
  #   - eigen_decomp_coarse_truth
  #   - eigen_decomp_X_truth
  #   - eigen_decomp_X_coarse_truth
  #   - eigen_decomp_est
  #
  # - step_5b
  #   - KL_cor_truth           (list of d_i x d_j matrices for i_j entries)
  #   - KL_cor_coarse_truth    (list of d_i x d_j matrices for i_j entries)
  #   - KL_cor_X_truth         (list of d_i x d_j matrices for i_j entries)
  #   - KL_cor_X_coarse_truth  (list of d_i x d_j matrices for i_j entries)
  #   - KL_cor_est             (list of d_i x d_j matrices for i_j entries)
  #
  # - kernel_params
  # - full              (boolean)    are we including truths in our estimation?
  # 
  # outputs:
  # 
  # - list of:
  #   - C_cond_ground_truth_full           (pm x pm matrix)
  #   - C_cond_coarse_ground_truth_full    (pm_est x pm_est matrix)
  #   - C_cond_truth_full                  (pm x pm matrix)
  #   - C_cond_coarse_truth_full           (pm_est x pm_est matrix)
  #   - C_cond_X_truth_full                (pm x pm matrix)
  #   - C_cond_X_coarse_truth_full         (pm_est x pm_est matrix)
  #   - C_cond_est_full                    (pm_est x pm_est matrix) 
  #   - C_cond_ground_truth_full_v2        (pm x pm matrix)
  #   - C_cond_coarse_ground_truth_full_v2 (pm_est x pm_est matrix)
  #
  # ----------------------------------------------------------------------------
  
  if(! full){
    eigen_decomp_est    <- step_4$eigen_decomp_est
    KL_cor_est          <- step_5b$KL_cor_est
    C_cond_list         <- correlation_estimation_KL_cor(eigen_decomp_est, KL_cor_est) 
    C_cond_est          <- C_cond_list$C_cond
    C_cond_est_unnorm   <- C_cond_list$C_cond_unnorm
    
    p <- length(eigen_decomp_est[[1]])
    m_est <- dim(C_cond_est[[1]])[1]
    
    C_cond_est_full            <- assemble_block_matrix_v2(C_cond_est, p, m_est)
    C_cond_est_unnorm_full     <- assemble_block_matrix_v2(C_cond_est_unnorm, p, m_est)
    
    return(list(C_cond_est_full = C_cond_est_full,
                C_cond_est_unnorm_full = C_cond_est_unnorm_full))
  }
  
  # 0) prep
  
  
  eigen_decomp_truth          <- step_4[[1]]
  eigen_decomp_coarse_truth   <- step_4[[2]]
  eigen_decomp_X_truth        <- step_4[[3]]
  eigen_decomp_X_coarse_truth <- step_4[[4]]
  eigen_decomp_est            <- step_4[[5]]
  
  KL_cor_truth          <- step_5b[[1]]
  KL_cor_coarse_truth   <- step_5b[[2]]
  KL_cor_X_truth        <- step_5b[[3]]
  KL_cor_X_coarse_truth <- step_5b[[4]]
  KL_cor_est            <- step_5b[[5]]
  
  p     <- length(eigen_decomp_truth[[1]])
  m     <- dim(eigen_decomp_truth[[2]][[1]])[1]
  m_est <- dim(eigen_decomp_coarse_truth[[2]][[1]])[1]
  
  # 1) estimation
  C_cond_truth          <- correlation_estimation_KL_cov(eigen_decomp_truth,          KL_cor_truth)
  C_cond_coarse_truth   <- correlation_estimation_KL_cov(eigen_decomp_coarse_truth,   KL_cor_coarse_truth)
  C_cond_X_truth        <- correlation_estimation_KL_cov(eigen_decomp_X_truth,        KL_cor_X_truth)
  C_cond_X_coarse_truth <- correlation_estimation_KL_cov(eigen_decomp_X_coarse_truth, KL_cor_X_coarse_truth)
  C_cond_est            <- correlation_estimation_KL_cov(eigen_decomp_est,            KL_cor_est)
  
  
  
  # 2) ground truth is cor_mat \otimes I_m
  C_cond_ground_truth_full          <- kronecker(kernel_params$prec_mat_truth$cor_mat, diag(m))
  C_cond_coarse_ground_truth_full   <- kronecker(kernel_params$prec_mat_truth$cor_mat, diag(m_est))
  
  # 3) or is ground truth cor_mat \otimes K_base?
  C_cond_ground_truth_full_v2          <- kronecker(kernel_params$prec_mat_truth$cor_mat, kernel_params$base_cov)
  C_cond_coarse_ground_truth_full_v2   <- kronecker(kernel_params$prec_mat_truth$cor_mat, kernel_params$base_cov_est)
  
  # 4) assemble and visualize the entire pm x pm block 
  
  C_cond_truth_full                 <- assemble_block_matrix_v2(C_cond_truth,               p, m)
  C_cond_coarse_truth_full          <- assemble_block_matrix_v2(C_cond_coarse_truth,        p, m_est)
  C_cond_X_truth_full               <- assemble_block_matrix_v2(C_cond_X_truth,             p, m)  
  C_cond_X_coarse_truth_full        <- assemble_block_matrix_v2(C_cond_X_coarse_truth,      p, m_est)
  C_cond_est_full                   <- assemble_block_matrix_v2(C_cond_est,                 p, m_est)
  
  
  
  
  return(list(C_cond_ground_truth_full = C_cond_ground_truth_full,
              C_cond_coarse_ground_truth_full = C_cond_coarse_ground_truth_full,
              C_cond_truth_full = C_cond_truth_full,
              C_cond_coarse_truth_full = C_cond_coarse_truth_full,
              C_cond_X_truth_full = C_cond_X_truth_full,
              C_cond_X_coarse_truth_full = C_cond_X_coarse_truth_full,
              C_cond_est_full = C_cond_est_full,
              C_cond_ground_truth_full_v2 = C_cond_ground_truth_full_v2,
              C_cond_coarse_ground_truth_full_v2 = C_cond_coarse_ground_truth_full_v2))
  
}

step_9b_eigenfunction_outers <- function(step_4, full = T){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: calculate outer products of normalized eigenfunctions
  #
  #
  # inputs:
  #
  # - step_4
  #   - eigen_decomp_truth               (list of 3 things)
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)          THEY ARE NOT NORMALIZED
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  #
  #   - eigen_decomp_coarse_truth
  #   - eigen_decomp_X_truth
  #   - eigen_decomp_X_coarse_truth
  #   - eigen_decomp_est
  #
  #
  # - full              (boolean)    are we including truths in our estimation?
  # 
  # outputs:
  # 
  # - list of:
  #
  #   - efunc_outer_truth                  (pc2 list of mxm matrices)
  #
  # ----------------------------------------------------------------------------
  
  if(! full){
    eigen_decomp_est <- step_4$eigen_decomp_est
    
    
    efunc_outer_list <- correlation_eigenfunction_outer(eigen_decomp_est)
    

    
    return(list(efunc_outer_truth = efunc_outer_list$efunc_outer_list,
                efunc_outer_truth_unnorm = efunc_outer_list$efunc_outer_list_unnorm))
  }
  
}

step_10_P_cond <- function(step_9, kernel_params, p, block, MP, full = T){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate precision operator once we have C_cond estiamtes
  #
  # 
  # inputs:
  #
  # - step_9
  #   - C_cond_ground_truth_full          (pm x pm matrix)
  #   - C_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - C_cond_truth_full                 (pm x pm matrix)
  #   - C_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - C_cond_X_truth_full               (pm x pm matrix)
  #   - C_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - C_cond_est_full                   (pm_est x pm_est matrix)  
  #
  # - kernel_params (list)     list of simulation parameters
  # - p             (integer)  number of processes
  # - block         (boolean)  do we take the inverse of each block? If false, take the inverse of the entire pm x pm matrix
  # - MP            (boolean)  do we use moore-penrose inverse? If false, to regular inverse.
  # - full          (boolean)  are we including truths in our estimation?
  #
  # outputs:
  #
  # - list of:
  #   - P_cond_ground_truth_full          (pm x pm matrix)
  #   - P_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - P_cond_truth_full                 (pm x pm matrix)
  #   - P_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - P_cond_X_truth_full               (pm x pm matrix)
  #   - P_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - P_cond_est_full                   (pm_est x pm_est matrix)
  #
  # ----------------------------------------------------------------------------
  
  if(! full){
    C_cond_est_full          <- step_9$C_cond_est_full
    C_cond_est_unnorm_full   <- step_9$C_cond_est_unnorm_full
    
    
    P_cond_est_full           <- estimate_precision_operator_v3(C_cond_est_full, p, block, MP)
    P_cond_est_unnorm_full    <- estimate_precision_operator_v3(C_cond_est_unnorm_full, p, block, MP)
    
    
    return(list(P_cond_est_full = P_cond_est_full,
                P_cond_est_unnorm_full = P_cond_est_unnorm_full))
  }
  
  # 1) prep
  C_cond_truth_full          <- step_9[[3]]
  C_cond_coarse_truth_full   <- step_9[[4]]
  C_cond_X_truth_full        <- step_9[[5]]
  C_cond_X_coarse_truth_full <- step_9[[6]]
  C_cond_est_full            <- step_9[[7]]
   
  m     <- dim(C_cond_truth_full)[1] / p
  m_est <- dim(C_cond_coarse_truth_full)[1] / p  
  
  
  # 2) estimate
  P_cond_truth_full          <- estimate_precision_operator_v3(C_cond_truth_full,          p, block, MP)
  P_cond_coarse_truth_full   <- estimate_precision_operator_v3(C_cond_coarse_truth_full,   p, block, MP)
  P_cond_X_truth_full        <- estimate_precision_operator_v3(C_cond_X_truth_full,        p, block, MP)  
  P_cond_X_coarse_truth_full <- estimate_precision_operator_v3(C_cond_X_coarse_truth_full, p, block, MP)
  P_cond_est_full            <- estimate_precision_operator_v3(C_cond_est_full,            p, block, MP)
  
  P_cond_truth          <- extract_block_structure_v2(P_cond_truth_full,          p, m)
  P_cond_coarse_truth   <- extract_block_structure_v2(P_cond_coarse_truth_full,   p, m_est)
  P_cond_X_truth        <- extract_block_structure_v2(P_cond_X_truth_full,        p, m)  
  P_cond_X_coarse_truth <- extract_block_structure_v2(P_cond_X_coarse_truth_full, p, m_est)
  P_cond_est            <- extract_block_structure_v2(P_cond_est_full,            p, m_est)
  
  # truth = prec_mat \otimes I_m
  P_cond_ground_truth_full           <- kronecker(kernel_params$prec_mat_truth$prec_mat, diag(m))
  P_cond_coarse_ground_truth_full    <- kronecker(kernel_params$prec_mat_truth$prec_mat, diag(m_est))
  P_cond_ground_truth                <- extract_block_structure_v2(P_cond_ground_truth_full, p, m)
  P_cond_coarse_ground_truth         <- extract_block_structure_v2(P_cond_coarse_ground_truth_full, p, m_est)
  
  
  return(list(
    P_cond_ground_truth_full = P_cond_ground_truth_full,
    P_cond_coarse_ground_truth_full = P_cond_coarse_ground_truth_full,
    P_cond_truth_full = P_cond_truth_full,
    P_cond_coarse_truth_full = P_cond_coarse_truth_full,
    P_cond_X_truth_full = P_cond_X_truth_full,
    P_cond_X_coarse_truth_full = P_cond_X_coarse_truth_full,
    P_cond_est_full = P_cond_est_full
  ))
  
}

step_11_HS_norms <- function(step_9, step_10, adj_mat_i, p, full = T){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: visualize HS norms of the P_cond estimates
  #
  # input:
  #
  # - step_9
  #   - C_cond_ground_truth_full          (pm x pm matrix)
  #   - C_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - C_cond_truth_full                 (pm x pm matrix)
  #   - C_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - C_cond_X_truth_full               (pm x pm matrix)
  #   - C_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - C_cond_est_full                   (pm_est x pm_est matrix)  
  #
  # - step_10
  #   - P_cond_ground_truth_full          (pm x pm matrix)
  #   - P_cond_coarse_ground_truth_full   (pm_est x pm_est matrix)
  #   - P_cond_truth_full                 (pm x pm matrix)
  #   - P_cond_coarse_truth_full          (pm_est x pm_est matrix)
  #   - P_cond_X_truth_full               (pm x pm matrix)
  #   - P_cond_X_coarse_truth_full        (pm_est x pm_est matrix)
  #   - P_cond_est_full                   (pm_est x pm_est matrix)
  #
  # - adj_mat_i    (p x p matrix of 0's and 1's) ground truth adjacencies for the i-th query
  # - p            (integer)
  # - full         (boolean)  are we including truths in our estimation?
  #
  # outputs:
  #
  # - list of:
  #   - w_mat_ground_truth        (p x p)   matrix of HS norms of the pm x pm ground truth
  #   - w_mat_coarse_ground_truth (p x p)   matrix of HS norms of the pm_est x pm_est ground truth 
  #   - w_mat_X_coarse_truth      (p x p)   matrix of HS norms of ...
  #   - w_mat_X_truth             (p x p)   matrix of HS norms of ...
  #   - w_mat_coarse_truth        (p x p)   matrix of HS norms of ...
  #   - w_mat_truth               (p x p)   matrix of HS norms of ...
  #   - w_mat_est                 (p x p)   matrix of HS norms of ...
  #
  # ----------------------------------------------------------------------------
  
  if(! full){
    
    C_cond_est_full <- step_9$C_cond_est_full
    C_cond_est_unnorm_full <- step_9$C_cond_est_unnorm_full
    P_cond_est_full <- step_10$P_cond_est_full
    P_cond_est_unnorm_full <- step_10$P_cond_est_unnorm_full
    
    m_est <- dim(P_cond_est_full)[1] / p
    w_mat_est <- hilbert_schmidt_norm_pm(P_cond_est_full, p, m_est)
    w_mat_est_unnorm <- hilbert_schmidt_norm_pm(P_cond_est_unnorm_full, p, m_est)
    C_HS_est <- hilbert_schmidt_norm_pm(C_cond_est_full, p, m_est)
    C_HS_est_unnorm <- hilbert_schmidt_norm_pm(C_cond_est_unnorm_full, p, m_est)
    
    return(
      list(w_mat_est = w_mat_est,
           w_mat_est_unnorm = w_mat_est_unnorm,
           C_HS_est = C_HS_est,
           C_HS_est_unnorm = C_HS_est_unnorm)
    )
  }
  
  # prep
  P_cond_ground_truth_full        <- step_10[[1]]
  P_cond_coarse_ground_truth_full <- step_10[[2]]
  P_cond_truth_full               <- step_10[[3]]
  P_cond_coarse_truth_full        <- step_10[[4]]
  P_cond_X_truth_full             <- step_10[[5]]
  P_cond_X_coarse_truth_full      <- step_10[[6]]
  P_cond_est_full                 <- step_10[[7]]
  
  m <- dim(P_cond_truth_full)[1] / p
  m_est <- dim(P_cond_coarse_truth_full)[1] / p
  
  # est
  w_mat_est            <- hilbert_schmidt_norm_pm(P_cond_est_full,            p, m_est)
  w_mat_X_coarse_truth <- hilbert_schmidt_norm_pm(P_cond_X_coarse_truth_full, p, m_est)
  w_mat_X_truth        <- hilbert_schmidt_norm_pm(P_cond_X_truth_full,        p, m)
  w_mat_coarse_truth   <- hilbert_schmidt_norm_pm(P_cond_coarse_truth_full,   p, m_est)
  w_mat_truth          <- hilbert_schmidt_norm_pm(P_cond_truth_full,          p, m)
  
  w_mat_ground_truth        <- hilbert_schmidt_norm_pm(P_cond_ground_truth_full, p, m)
  w_mat_coarse_ground_truth <- hilbert_schmidt_norm_pm(P_cond_coarse_ground_truth_full, p, m_est)

  
  # hilbert schmidt normalized HS norms 
  
  # w_mat_normalized_X_coarse_truth      <- hilbert_schmidt_norm_pm_normalize(P_cond_X_coarse_truth_full, p, m_est)
  # w_mat_normalized_X_truth             <- hilbert_schmidt_norm_pm_normalize(P_cond_X_truth_full, p, m)    
  # w_mat_normalized_coarse_truth        <- hilbert_schmidt_norm_pm_normalize(P_cond_coarse_truth_full, p, m_est)
  # w_mat_normalized_truth               <- hilbert_schmidt_norm_pm_normalize(P_cond_truth_full, p, m)
  # w_mat_normalized_est                 <- hilbert_schmidt_norm_pm_normalize(P_cond_est_full, p, m_est)  
  
  # diag(w_mat_normalized_X_coarse_truth) <- 0
  # diag(w_mat_normalized_X_truth) <- 0
  # diag(w_mat_normalized_coarse_truth) <- 0
  # diag(w_mat_normalized_truth) <- 0
  # diag(w_mat_normalized_est) <- 0
  
  # w_mat_normalized_ground_truth        <- hilbert_schmidt_norm_pm_normalize(P_cond_ground_truth_full, p, m)
  # w_mat_normalized_coarse_ground_truth <- hilbert_schmidt_norm_pm_normalize(P_cond_coarse_ground_truth_full, p, m_est)
  # diag(w_mat_normalized_ground_truth) <- 0
  # diag(w_mat_normalized_coarse_ground_truth) <- 0    
  

  
  return(list(
    w_mat_ground_truth = w_mat_ground_truth,
    w_mat_coarse_ground_truth = w_mat_coarse_ground_truth,
    w_mat_truth = w_mat_truth,
    w_mat_coarse_truth = w_mat_coarse_truth,
    w_mat_X_truth = w_mat_X_truth,
    w_mat_X_coarse_truth = w_mat_X_coarse_truth,
    w_mat_est = w_mat_est
  ))
  
}

step_12_ROC <- function(step_11, adj_mat_i, full = T){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estiamte ROC curves
  #
  # inputs:
  #
  # - step_11
  #   - w_mat_ground_truth        (p x p)   matrix of HS norms of the pm x pm ground truth
  #   - w_mat_coarse_ground_truth (p x p)   matrix of HS norms of the pm_est x pm_est ground truth 
  #   - w_mat_X_coarse_truth      (p x p)   matrix of HS norms of ...
  #   - w_mat_X_truth             (p x p)   matrix of HS norms of ...
  #   - w_mat_coarse_truth        (p x p)   matrix of HS norms of ...
  #   - w_mat_truth               (p x p)   matrix of HS norms of ...
  #   - w_mat_est                 (p x p)   matrix of HS norms of ...
  #
  #
  # - adj_mat_i    (p x p matrix of 0's and 1's)   denoting ground truth adjacencies with 0's on the diagonal
  #
  #
  # outputs:
  #
  # - list of:
  #   - roc_ground_truth           (list of roc outputs)
  #   - roc_coarse_ground_truth
  #   - roc_truth
  #   - roc_coarse_truth
  #   - roc_X_truth
  #   - roc_X_coarse_truth
  #   - roc_est
  #
  # ----------------------------------------------------------------------------
  
  # prep
  w_mat_ground_truth        <- step_11[[1]]
  w_mat_coarse_ground_truth <- step_11[[2]]
  w_mat_truth               <- step_11[[3]]
  w_mat_coarse_truth        <- step_11[[4]]
  w_mat_X_truth             <- step_11[[5]]
  w_mat_X_coarse_truth      <- step_11[[6]]
  w_mat_est                 <- step_11[[7]]
  
  # estimate 
  roc_ground_truth        <- roc_with_threshold(w_mat_ground_truth,        adj_mat_i, 'Ground Truth')
  roc_coarse_ground_truth <- roc_with_threshold(w_mat_coarse_ground_truth, adj_mat_i, 'Coarse Ground Truth')
  roc_truth               <- roc_with_threshold(w_mat_truth,               adj_mat_i, 'Truth Theory')   
  roc_coarse_truth        <- roc_with_threshold(w_mat_coarse_truth,        adj_mat_i, 'Coarse Truth Theory')
  roc_X_truth             <- roc_with_threshold(w_mat_X_truth,             adj_mat_i, 'Truth X')
  roc_X_coarse_truth      <- roc_with_threshold(w_mat_X_coarse_truth,      adj_mat_i, 'Coarse Truth X')
  roc_est                 <- roc_with_threshold(w_mat_est,                 adj_mat_i, 'Estimate')
  

  
  
  return(list(roc_ground_truth = roc_ground_truth,
              roc_coarse_ground_truth = roc_coarse_ground_truth,
              roc_truth = roc_truth,
              roc_coarse_truth = roc_coarse_truth,
              roc_X_truth = roc_X_truth,
              roc_X_coarse_truth = roc_X_coarse_truth,
              roc_est = roc_est))
  
}
