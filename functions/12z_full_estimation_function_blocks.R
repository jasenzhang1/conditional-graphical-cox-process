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
    
    name_i <- paste0('eigen_decomp_', core_names[i], '_eig3')
    
    temp_var <- prep_eigendecomposition_ii(step_3[[input_names[i]]], p)  # prep
    result[[name_i]] <- compute_eigendecomposition_ii(temp_var, same_basis, constant_d)       # then compute

  }
  
  return(result)
}

step_4_eigendecomp_mfpca <- function(step_3, p, same_basis, constant_d){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: copy of 'step_4_eigendecomp', but we assume mfpca setting
  #
  #
  # inputs:
  #
  # - step_3
  #   - g_ij_suffix                (list of m x m matrices for i_j entries)
  #
  # - p                            (scalar) number of processes
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
    
    name_i <- paste0('eigen_decomp_', core_names[i], '_eig4')
  
    result[[name_i]] <- compute_eigendecomposition_mfpca(step_3[[input_names[i]]], p, same_basis, constant_d)       # then compute
  }
  
  return(result)
}

step_4_eigendecomp_troubleshoot <- function(step_3, p, eigen_setting){
  
  
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
  # - eigen_setting                (string)  'trig_and_joint',  'trig_simple', 'mfpca
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
  
  if(eigen_setting == 'trig_simple'){
    basis_settings <- c(T)
    constant_d_settings <- c(2)
  } else{
    basis_settings <- c(T, T, F)
    constant_d_settings <- c(2, NA, NA) 
  }
  


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
    
    name_i <- paste0('w_mat_KL_GIC_local_', core_names[i])
    
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
    
    name_i <- paste0('C_HS_KL_GIC_local_', core_names[i])
    
    result[[name_i]] <- step_11_GIC_bundle[[input_names[i]]][['C_HS']]
    
  }
  
  return(result)
  
}

step_11xy_HS_norms_from_KL_GIC <- function(step_11_GIC_bundle){
  
  
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
  #     - tau_c_local_X_truth etc...
  #
  #   - step_11y
  #     - tau_p_local_est
  #     - tau_p_local_X_truth etc...
  #
  # ----------------------------------------------------------------------------
  
  step_11x <- list()
  step_11y <- list()
  
  # 1) grab names
  
  core_names <- step_00_grab_ID(names(step_11_GIC_bundle), 'GIC_KL')
  
  input_names <- names(step_11_GIC_bundle)
  
  # 2) for each core name `est`, `X_truth` etc... get the resulting name, apply the function on it, and store it
  for(i in 1:length(core_names)){
    
    # obtain tau_c and tau_p from the entry
    
    name_i_tau_c <- paste0('tau_c_local_', core_names[i])
    name_i_tau_p <- paste0('tau_p_local_', core_names[i])
    
    step_11x[[name_i_tau_c]] <- step_11_GIC_bundle[[input_names[i]]][['tau_c']]
    step_11y[[name_i_tau_p]] <- step_11_GIC_bundle[[input_names[i]]][['tau_p']]
    
  }
  
  return(list(step_11x = step_11x, 
              step_11y = step_11y))
  
  
}

step_12_ROC <- function(step_11, adj_mat_i){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: estimate ROC curves
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
    
    if(grepl("GIC", core_names[i])){
      result[[name_i]]  <- roc_for_thresholded_w_mat(step_11[[input_names[i]]], adj_mat_i)   # one procedure for w_mat that's zeroed out       (don't do ROC)
    } else{
      result[[name_i]]  <- roc_for_raw_w_mat(step_11[[input_names[i]]], adj_mat_i)          # one procedure for w_mat that's not zeroed out   (do ROC)
    }
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
