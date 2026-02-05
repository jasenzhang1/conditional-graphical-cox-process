step_00_grab_ID <- function(vec, prefix, suffix = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: grab names from a vector of names that share the same prefix and perhaps a suffix
  # 
  #       E.G. c('prefix_name1_suffix', 'prefix_name2_suffix', 'diff_prefix_name3')
  #
  #       output: c('name1', 'name2')
  #
  # inputs:
  #
  # - vec     (vector of untrimmed names)
  # - prefix  (string)                      prefix name, not including the underscore between prefix and name
  # - suffix  (string)                      suffix name, not including the underscore between name and suffix
  #
  # ouptuts:
  #
  # - out     (vector of names)
  #
  # ----------------------------------------------------------------------------
  
  # Pattern to match entries that start with the prefix
  prefix_pattern <- paste0("^", prefix, "_")
  
  # Keep only matching entries
  keep_idx <- grepl(prefix_pattern, vec)
  vec_filtered <- vec[keep_idx]
  
  # Remove prefix
  out <- sub(prefix_pattern, "", vec_filtered)
  
  # If suffix provided, remove it too
  if (!is.null(suffix)) {
    suffix_pattern <- paste0("_", suffix, "$")
    out <- sub(suffix_pattern, "", out)
  }
  
  return(out)
}

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

step_0_keep_events <- function(dataset, k_vec, i_vec){
  
  
  # ----------------------------------------------------------------------------
  #
  # keep point process events for subject k and processes i
  #
  #
  # input:
  #
  # - dataset    (list)
  # - k          (vector)            vector of subject id's
  # - i_vec      (vector)            vector of process ID's
  #
  # 
  # outputs:
  # 
  # - events     (list of q vectors)
  #
  # ----------------------------------------------------------------------------
  
  keys <- paste0(
    rep(k_vec, each = length(i_vec)),
    "_",
    rep(i_vec, times = length(k_vec))
  )
  
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
  #   - [[3]] y_c_strata_sel    subset of stratas in case some subjects are discarded
  #   - [[4]] query_y_cs
  #   - [[5]] patient_sel
  #   - [[6]] feature_sel
  #   
  #
  # ----------------------------------------------------------------------------
  

  
  # 3.1) first, convert dataset into a format that can be used for estimation
  
  data_df4 <- convert_data_for_estimation_event_times(dataset$event_times) 

  y_c_strata_full <- dataset$Y_continuous
  y_c_strata <- dataset$Y_continuous_k
  
  query_y_cs <- dataset$simulation_params$query_y_cs  
  

  
  # size of dataset
  print(paste0('number of subjects: ', dim(y_c_strata)[1]))
  print(paste0('number of spikes: ', nrow(data_df4)))
  
  
  # 3.3) in case there are some subjects or features with no events
  

  patient_sel = unique(data_df4$subject_num) %>% sort()
  feature_sel = unique(data_df4$feature_id) %>% sort()
  
  return(list(data_df4 = data_df4,
              y_c_strata_full = y_c_strata_full,
              y_c_strata = y_c_strata,
              query_y_cs = query_y_cs,
              patient_sel = patient_sel,
              feature_sel = feature_sel
              ))
  
  
  
}


step_1_log_intensities <- function(data_df4, time_grid_est){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: obtain log-intensities (X) and intensities (Lambda)
  #
  # 
  # inputs:
  # 
  # - data_df4          (data.table with `time`, `feature_id`, `subject_num`)
  # - time_grid_est     (m_est-dim vec of discretized times)
  #
  #
  # outputs:
  # 
  # - list of:
  #   - step_1:
  #     - X_k_est                 (p x m x n)
  #   - step_1b:
  #     - Lambda_k_est            (p x m x n)
  #
  # ----------------------------------------------------------------------------
  
  X_k_est <- subject_specific_log_intensity(data_df4, time_grid_est)
  Lambda_k_est <- exp(X_k_est)
  
  step_1 <- list(X_k_est = X_k_est)
  step_1b <- list(Lambda_k_est = Lambda_k_est)
  
  
  return(list(step_1 = step_1,
              step_1b = step_1b))
  
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
  # 
  # - if(full): list of:
  #   - rho_i_truth               (p x m)
  #   - rho_i_coarse_truth        (p x m_est)
  #   - rho_i_X_truth             (p x m)
  #   - rho_i_X_coarse_truth      (p x m_est)
  #   - rho_i_est                 (p x m_est)
  #   - rho_list                  (list format, [[1]] = rho_i, [[2]] = rho_ij, rho_ij may be rho_ii only)
  #   - weights                   (n-dim vector)
  #   - y_c_s                     (n-dim vector)
  #
  # - otherwise, list of:
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

step_3_g_ij <- function(step_2, step_2b, i_neq_j){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: covariance function (g_ij) estimation
  #
  # 
  # inputs:
  #
  # - step_2:
  #   - rho_i_suffix              (p x m)
  #   - suffixes = truth, X_truth, est, etc
  #
  #
  # - step_2b: 
  #   - rho_ii_suffix              (list of m x m matrices, could have only i_i or i_j as well)
  #
  # - i_neq_j                    (boolean) do we include i =/= j terms?
  # 
  # outputs:
  #
  # - list of:
  #   - g_ij_suffix                      (list of m x m matrices for i_j entries)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_2), 'rho_i')
  
  input_names_step_2 <- names(step_2)
  input_names_step_2b <- names(step_2b)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('g_ij_', core_names[i])
    
    if(i_neq_j){
      result[[name_i]] <- estimate_covariance_functions_ij(step_2[[input_names_step_2[i]]], step_2b[[input_names_step_2b[i]]])
    } else{
      result[[name_i]] <- estimate_covariance_functions_ii(step_2[[input_names_step_2[i]]], step_2b[[input_names_step_2b[i]]])
    }
  }
  
  return(result)
  
}

step_4_eigendecomp <- function(step_3, p, same_basis, constant_d){
  

  # ----------------------------------------------------------------------------
  # 
  # GOAL: eigendecomposition of G_ii
  #
  # inputs:
  #
  # - step_3
  #   - g_ij_suffix                (list of m x m matrices for i_j entries)
  #
  # - p                            (scalar)
  # - same_basis                   (boolean)  does each process use the same basis?
  # - constant_d                   (integer)  does each process use a constant amount of components, if not null
  #
  #
  # outputs:
  #
  # - list of:
  #   - eigen_decomp_suffix             (list of 3 things)
  #     - [[1]] eigenvalues            (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors           (list of p matrices of m x d_i)
  #     - [[3]] n_dims                 (list of p integers denoting d_i)
  #
  #
  # ----------------------------------------------------------------------------
  
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_3), 'g_ij')
  
  input_names <- names(step_3)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('eigen_decomp_', core_names[i])
    
    temp_var <- prep_eigendecomposition_ii(step_3[[input_names[i]]], p)  # prep
    result[[name_i]] <- compute_eigendecomposition_ii(temp_var, same_basis, constant_d)       # then compute

  }
  
  return(result)
}

step_4_eigendecomp_troubleshoot <- function(step_3, p){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: copy of 'step_4_eigendecomp', but we calculate all 3 settings:
  #
  # 1) same_basis = T, constant_d = 2
  # 2) same_basis = T, constant_d = NULL
  # 3) same_basis = F, constant_d = NULL
  #
  # inputs:
  #
  # - step_3
  #   - g_ij_suffix                (list of m x m matrices for i_j entries)
  #
  # - p                            (scalar)
  #
  #
  # outputs:
  #
  # - list of:
  #   - eigen_decomp_suffix             (list of 3 things)
  #     - [[1]] eigenvalues            (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors           (list of p matrices of m x d_i)
  #     - [[3]] n_dims                 (list of p integers denoting d_i)
  #
  #
  # ----------------------------------------------------------------------------
  
  basis_settings <- c(T, T, F)
  constant_d_settings <- c(2, NA, NA)

  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_3), 'g_ij')
  
  input_names <- names(step_3)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    for(j in 1:length(basis_settings)){
      same_basis <- basis_settings[j]
      constant_d <- constant_d_settings[j]
      
      name_i <- paste0('eigen_decomp_', core_names[i], '_eig', j)
      
      temp_var <- prep_eigendecomposition_ii(step_3[[input_names[i]]], p)  # prep
      
      if(is.na(constant_d)){
        result[[name_i]] <- compute_eigendecomposition_ii(temp_var, same_basis, constant_d = NULL)   
      } else{
        result[[name_i]] <- compute_eigendecomposition_ii(temp_var, same_basis, constant_d)     
      }
    }
  }
  
  return(result)
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

step_5_KL_covariance <- function(step_3, step_4){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: estimate covariances of the KL coefficients for CPGM method
  #
  # inputs:
  #
  # - step_3
  #   - g_ij_suffix                (list of m x m matrices for i_j entries)
  # 
  # - step_4
  #   - eigen_decomp_suffix               (list of 3 things)
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  #
  #
  #
  # outputs:
  #
  # - list of:
  #   - KL_cov_suffix           (list of d_i x d_j matrices for i_j entries)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_3), 'g_ij')
  
  input_names_step_3 <- names(step_3)
  input_names_step_4 <- names(step_4)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    step_3_i <- input_names_step_3[i]
    step_4_i <- input_names_step_4[i]
    
    name_i <- paste0('KL_cov_', core_names[i])
    
    result[[name_i]] <- estimate_KL_covariance(step_3[[step_3_i]], step_4[[step_4_i]]$eigenfunctions)     
    
  }
  
  return(result)
  
}

step_5_KL_covariance_eigencases <- function(step_3, step_4){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: step_5_KL_covariance troubleshooting
  #
  # inputs:
  #
  # - step_3
  #   - g_ij_suffix                (list of m x m matrices for i_j entries)
  # 
  # - step_4
  #   - eigen_decomp_suffix_eigk               (list of 3 things)
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  #
  # NOTE: eigk is either eig1, eig2, or eig3
  #
  #
  #
  # outputs:
  #
  # - list of:
  #   - KL_cov_suffix           (list of d_i x d_j matrices for i_j entries)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names_step_4 <- step_00_grab_ID(names(step_4), 'eigen_decomp')
  
  input_names_step_3 <- rep(names(step_3), each = 3)
  input_names_step_4 <- names(step_4)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names_step_4)){
    
    step_3_i <- input_names_step_3[i]
    step_4_i <- input_names_step_4[i]
    
    name_i <- paste0('KL_cov_', core_names_step_4[i])
    
    result[[name_i]] <- estimate_KL_covariance(step_3[[step_3_i]], step_4[[step_4_i]]$eigenfunctions)     
    
  }
  
  return(result)
  
}

step_5b_KL_correlation <- function(step_5, p){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: estimate correlations of the KL coefficients for CPGM method
  #
  # inputs:
  #
  # - step_5
  #   - KL_cov_suffix          (list of d_i x d_j matrices for i_j entries)
  # 
  # - p            (integer)
  #
  #
  # outputs:
  #
  # - list of:
  #   - KL_cor_suffix          (list of d_i x d_j matrices for i_j entries)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_5), 'KL_cov')
  
  input_names <- names(step_5)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('KL_cor_', core_names[i])
    
    result[[name_i]] <- estimate_KL_correlation(step_5[[input_names[i]]], p)     
    
  }
  
  return(result)
}

step_5c_KL_precision <- function(step_5b, p){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: estimate precisions of the KL coefficients for CPGM method
  #
  # inputs:
  #
  # - step_5b
  #   - KL_cor_suffix          (list of d_i x d_j matrices for i_j entries)
  # 
  # - p            (integer)
  #
  #
  # outputs:
  #
  # - list of:
  #   - KL_prec_suffix          (list of d_i x d_j matrices for i_j entries)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_5b), 'KL_cor')
  
  input_names <- names(step_5b)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('KL_prec_', core_names[i])
    
    result[[name_i]] <- estimate_KL_precision(step_5b[[input_names[i]]], p)     
    
  }
  
  return(result)
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

step_9_C_cond_from_KL_cor <- function(step_4, step_5b){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: construct correlation matrix from CPGM method USING cor() instead of cov
  #
  #
  # inputs:
  #
  # - step_4
  #   - eigen_decomp_suffix               (list of 3 things)
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)          THEY ARE NOT NORMALIZED
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  #
  #
  # - step_5b
  #   - KL_cor_suffix         (list of d_i x d_j matrices for i_j entries)
  #
  # 
  # outputs:
  # 
  # - list of:
  #   - C_cond_suffix           (list of i_j m x m matrices)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_4), 'eigen_decomp')
  
  input_names_step_4 <- names(step_4)
  input_names_step_5b <- names(step_5b)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    step_4_i <- input_names_step_4[i]
    step_5b_i <- input_names_step_5b[i]
    
    name_i <- paste0('C_cond_', core_names[i])
    #name_i_unnorm <- paste0(name_i, '_unnorm')
    
    temp_list <- correlation_estimation_KL_cor(step_4[[step_4_i]],
                                               step_5b[[step_5b_i]])
    
    result[[name_i]] <- temp_list$C_cond    
    #result[[name_i_unnorm]] <- temp_list$C_cond_unnorm  
  }
  
  return(result)
  
}

step_9b_eigenfunction_outers <- function(step_4){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: calculate outer products of normalized eigenfunctions
  #
  #
  # inputs:
  #
  # - step_4
  #   - eigen_decomp_suffix               (list of 3 things)
  #     - [[1]] eigenvalues  (list of p vectors of eigenvalues)
  #     - [[2]] eigenvectors (list of p matrices of m x d_i)          THEY ARE NOT NORMALIZED
  #     - [[3]] n_dims       (list of p integers denoting d_i)
  #
  # 
  # outputs:
  # 
  # - list of:
  #
  #   - efunc_outer_suffix                  (pc2 list of mxm matrices)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_4), 'eigen_decomp')
  
  input_names <- names(step_4)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    temp_result <- correlation_eigenfunction_outer(step_4[[input_names[i]]])
    names(temp_result) <- paste0(names(temp_result), '_', core_names[i])
    
    result <- c(result, temp_result)
  }
  
  return(result)
}

step_10_P_cond <- function(step_9, p, block, MP){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate precision operator once we have C_cond estiamtes
  #
  # 
  # inputs:
  #
  # - step_9
  #   - C_cond_suffix          (i_j list of mxm matrices)
  #
  # - kernel_params (list)     list of simulation parameters
  # - p             (integer)  number of processes
  # - block         (boolean)  do we take the inverse of each block? If false, take the inverse of the entire pm x pm matrix
  # - MP            (boolean)  do we use moore-penrose inverse? If false, to regular inverse.
  #
  # outputs:
  #
  # - list of:
  #   - P_cond_suffix                     (i_j list of mxm matrices)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_9), 'C_cond')
  
  input_names <- names(step_9)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    # make full --> invert --> make into list again
    
    name_i <- paste0('P_cond_', core_names[i])
    
    m <- dim(step_9[[input_names[i]]][[1]])[1]
    
    C_cond_full <- assemble_block_matrix_v2(step_9[[input_names[i]]], p, m)
    P_cond_full <- estimate_precision_operator_v3(C_cond_full, p, block, MP)
    result[[name_i]] <- extract_block_structure_v2(P_cond_full, p, m)
    
  }
  
  return(result)
}

step_11_HS_norms_from_KL <- function(step_5c, p){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: get w_mat from KL_prec
  #
  #
  # inputs:
  #
  # - step_5c (list)
  #   - KL_prec_suffix
  #
  # - p      (integer)
  #
  # outputs:
  #
  # - step_11 (list)
  #   - w_mat_KL_suffix
  # 
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_5c), 'KL_prec')
  input_names <- names(step_5c)

  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('w_mat_KL_', core_names[i])
    
    result[[name_i]]  <- hilbert_schmidt_norm_list_to_mat(step_5c[[input_names[i]]], p)   

  }
  
  return(result)
  
}

step_11b_HS_norms_from_KL <- function(step_5b, p){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: get C_HS from KL_cor
  #
  #
  # inputs:
  #
  # - step_5b (list)
  #   - KL_cor_suffix
  #
  # - p      (integer)
  #
  # outputs:
  #
  # - step_11 (list)
  #   - C_HS_KL_suffix
  # 
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_5b), 'KL_cor')
  input_names <- names(step_5b)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('C_HS_KL_', core_names[i])
    
    result[[name_i]]  <- hilbert_schmidt_norm_list_to_mat(step_5b[[input_names[i]]], p)   
    
  }
  
  return(result)
  
}

step_11_HS_norms <- function(step_9, step_10, p){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: visualize HS norms of the P_cond estimates
  #
  # input:
  #
  # - step_9
  #   - C_cond_suffix          (i_j list of mxm matrices)
  #
  # - step_10
  #   - P_cond_suffix          (i_j list of mxm matrices)
  #
  # - adj_mat_i    (p x p matrix of 0's and 1's) ground truth adjacencies for the i-th query
  # - p            (integer)
  # - full         (boolean)  are we including truths in our estimation?
  #
  # outputs:
  #
  # - list of the following:
  #   - list of 
  #     - w_mat_suffix        (p x p)   matrix of HS norms of the pm x pm matrices
  #
  #   - list of
  #     - C_HS_suffix
  #
  # ----------------------------------------------------------------------------
  
  result1 <- list()
  result2 <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_9), 'C_cond')
  
  input_names_step_9 <- names(step_9)
  input_names_step_10 <- names(step_10)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('w_mat_', core_names[i])
    name2_i <- paste0('C_HS_', core_names[i])
    
    m <- dim(step_9[[input_names_step_9[i]]][[1]])[1]
    
    C_cond_full <- step_9[[input_names_step_9[i]]]    %>% assemble_block_matrix_v2(p, m)
    P_cond_full <- step_10[[input_names_step_10[i]]]  %>% assemble_block_matrix_v2(p, m)
    
    
    result1[[name_i]]  <- hilbert_schmidt_norm_pm(P_cond_full, p, m)    
    result2[[name2_i]] <- hilbert_schmidt_norm_pm(C_cond_full, p, m)
  }
  
  return(list(step_11 = result1,
              step_11b = result2))
  
}

steps_10_11_GIC_from_KL <- function(step_5b, p, W_y){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate the adjacency structure with double thresholding
  #
  # 
  # inputs:
  #
  # - step_5b
  #   - KL_cor_suffix                        (list of i_j matrix)  
  #
  # - p             (integer)  number of processes
  # - W_y           (scalar)   weighted sample size of this y_c_query
  #
  # outputs:
  #
  # - list of:
  #   - GIC_KL_suffix                           (list of items)
  #
  # ----------------------------------------------------------------------------
  
  # 1) grab names
  result <- list()
  core_names <- step_00_grab_ID(names(step_5b), 'KL_cor')
  input_names <- names(step_5b)
  
  # for each name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('GIC_KL_', core_names[i])
    
    result[[name_i]] <- GIC_algorithm(step_5b[[input_names[i]]], p, W_y)
  }
  
  return(result)
  
  
}

steps_10_11_GIC <- function(step_9, p, W_y){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate the adjacency structure with double thresholding
  #
  # 
  # inputs:
  #
  # - step_9
  #   - C_cond_est                        (list of i_j matrix)  
  #
  # - p             (integer)  number of processes
  # - W_y           (scalar)   weighted sample size of this y_c_query
  #
  # outputs:
  #
  # - list of:
  #   - GIC_est                           (list of items)
  #
  # ----------------------------------------------------------------------------
  
  # 1) grab names
  result <- list()
  IDs <- step_00_grab_ID(names(step_9), 'C_cond', 'full')
  names_0 <- names(step_9)
  
  # for each name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(IDs)){
    
    name_i <- paste0('GIC_', IDs)
    
    result[[name_i]] <- GIC_algorithm(step_9[[names_0[i]]], p, W_y)
  }
    
  return(result)

  
}

step_11_HS_norms_from_KL_GIC <- function(step_11_GIC_bundle){
  
  
  # ----------------------------------------------------------------------------  
  #
  # GOAL: obtain HS_norms from the GIC bundle
  #
  # 
  # inputs:
  #
  # - step_11_GIC_bundle
  #   - GIC_KL_suffix                        (list of items)  
  #
  #
  # outputs:
  #
  # - list of:
  #   - w_mat_KL_GIC_suffix                  (list of pxp HS matrices)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_11_GIC_bundle), 'GIC_KL')
  
  input_names <- names(step_11_GIC_bundle)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    # obtain w_mat from the entry
    
    name_i <- paste0('w_mat_KL_GIC_', core_names[i])
    
    result[[name_i]] <- step_11_GIC_bundle[[input_names[i]]][['w_mat']]
    
  }
  
  return(result)
  
}

step_11b_HS_norms_from_KL_GIC <- function(step_11_GIC_bundle){
  
  
  # ----------------------------------------------------------------------------  
  #
  # GOAL: obtain HS_norms for correlation operator from the GIC bundle
  #
  # 
  # inputs:
  #
  # - step_11_GIC_bundle
  #   - GIC_KL_suffix                        (list of items)  
  #
  #
  # outputs:
  #
  # - list of:
  #   - C_HS_KL_GIC_suffix                  (list of pxp HS matrices)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_11_GIC_bundle), 'GIC_KL')
  
  input_names <- names(step_11_GIC_bundle)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    # obtain w_mat from the entry
    
    name_i <- paste0('C_HS_KL_GIC_', core_names[i])
    
    result[[name_i]] <- step_11_GIC_bundle[[input_names[i]]][['C_HS']]
    
  }
  
  return(result)
  
}

step_11x_HS_norms_from_KL_GIC <- function(step_11_GIC_bundle){
  
  
  # ----------------------------------------------------------------------------  
  #
  # GOAL: obtain tau_c and tau_p from the GIC bundle
  #
  # 
  # inputs:
  #
  # - step_11_GIC_bundle
  #   - GIC_KL_suffix                        (list of items)  
  #
  #
  # outputs:
  #
  # - list of:
  #   - step_11x
  #     - tau_c_local_est
  #     - tau_p_local_est
  #     - tau_c_local_X_truth etc...
  #
  # ----------------------------------------------------------------------------
  
  step_11x <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_11_GIC_bundle), 'GIC_KL')
  
  input_names <- names(step_11_GIC_bundle)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    # obtain tau_c and tau_p from the entry
    
    name_i_tau_c <- paste0('tau_c_local_', core_names[i])
    name_i_tau_p <- paste0('tau_p_local_', core_names[i])
    
    step_11x[[name_i_tau_c]] <- step_11_GIC_bundle[[input_names[i]]][['tau_c']]
    step_11x[[name_i_tau_p]] <- step_11_GIC_bundle[[input_names[i]]][['tau_p']]
    
  }
  
  return(step_11x)
  
  
}

step_12_ROC <- function(step_11, adj_mat_i){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estiamte ROC curves
  #
  # inputs:
  #
  # - step_11
  #   - w_mat_suffix        (p x p)   matrix of HS norms of the pm x pm ground truth
  #
  #
  # - adj_mat_i    (p x p matrix of 0's and 1's)   denoting ground truth adjacencies with 0's on the diagonal
  #
  #
  # outputs:
  #
  # - list of:
  #   - roc_suffix           (list of roc outputs)
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_11), 'w_mat')
  
  input_names <- paste0('w_mat_', core_names)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('roc_', core_names[i])
    
    result[[name_i]]  <- roc_with_threshold(step_11[[input_names[i]]], adj_mat_i, core_names[i])

  }
  
  return(result)
  
  
}

step_12b_adj_mat <- function(step_11){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: get pxp adjacency matrices from w_mat. 
  #
  #       label an edge if its HS_norm is > 0
  # 
  #
  # inputs:
  #
  # step_11 (list)
  #   - w_mat_suffix
  #
  # 
  # outputs:
  #
  # - step_12b (list)
  #   - adj_mat_suffix
  #
  #
  # ----------------------------------------------------------------------------
  
  result <- list()
  tol <- 1e-6
  p <- dim(step_11[[1]])[1]
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_11), 'w_mat')
  
  input_names <- paste0('w_mat_', core_names)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    name_i <- paste0('adj_mat_', core_names[i])
    

    adj_mat <- step_11[[input_names[i]]] > tol
    diag(adj_mat) <- 0
    
    result[[name_i]]  <- adj_mat
    
  }
  
  return(result)
  
}
