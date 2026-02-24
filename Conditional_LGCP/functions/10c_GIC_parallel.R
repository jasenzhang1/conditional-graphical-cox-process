# setup_gic.R
GIC_step1_precompute <- function(C_cond, p, W_y, save_path, id_suffix) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: obtain quantiles for tau_c to iterate through
  #
  #
  # Input: 
  #
  # - C_cond        (i_j list of any sized matrix, d_i x d_j or m x m)
  # - p             (scalar)
  # - W_y           (scalar)           effective sample size
  # - save_path     (string)          'temp_data/simu/GIC_local'
  # - id_suffix     (string)          'est' or 'X_truth'
  #
  # 
  # Output: 
  #
  # - list of the following metrics:
  #
  #   - tau_c   (scalar)    threshold where if your HS norm is less than this, you are zeroed out 
  #   - tau_p   (scalar) 
  #
  # ----------------------------------------------------------------------------
  
  # 1) Set diagonal identity and find off-diagonal blocks
  C_cond <- GIC_set_CXX_diag_identity(C_cond)
  block_names <- names(C_cond)
  off_diag_indices <- which(sapply(strsplit(block_names, "_"), function(x) x[1] != x[2]))
  
  # 2) Get threshold list for the outer loop (tau_c)
  threshold_list_c <- GIC_get_percentile_info(C_cond, off_diag_indices)
  
  # 3) naming
  
  file_name <- paste0('GIC_local_initial_data_', id_suffix, '.RData')
  
  save(C_cond, p, W_y, off_diag_indices, threshold_list_c, 
       file = paste0(save_path, "/", file_name))
  
  # 4) Return the number of tasks for Bash to know the loop range
  return(length(threshold_list_c$hs_vals))
}

GIC_step2_iterate_tau_c <- function(temp_file_dir, id_suffix, k){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: parallelize by working with a specific tau_c[k] value
  #
  # inputs:
  #
  # - temp_file_dir   (string)   'temp_data/simu'
  # - id_suffix       (string)   'est' or 'X_truth'
  # - k               (integer)   index of tau_c quantiles
  #
  # ----------------------------------------------------------------------------
  
  # 1. Load the ID-specific initial data created in Part 1
  # This contains C_cond, p, W_y, and threshold_list_c for this specific ID
  load(paste0(temp_file_dir, "/GIC_local_initial_data_", id_suffix, ".RData"))
  
  # 2. Run the Step 2 logic
  # Note: Ensure GIC_step2_iterate_tau_c is updated to accept/save with id_suffix
  # Or simply handle the save logic here:
  
  tau_c <- threshold_list_c$hs_vals[k]
  excluded_indices <- threshold_list_c$index_path[[k]]
  
  C_cond_thresh <- C_cond
  for (idx in excluded_indices) { C_cond_thresh[[idx]][] <- 0 }
  
  C_cond_full <- assemble_block_matrix_irregular(C_cond_thresh, p)
  Theta_full  <- ginv(C_cond_full$block_matrix)
  Theta_cond  <- extract_block_matrix_irregular(Theta_full, C_cond_full$row_borders, C_cond_full$col_borders)
  
  threshold_list_p <- GIC_get_percentile_info(Theta_cond, off_diag_indices)
  
  # Save intermediate data for the tau_p workers, specific to this ID and K
  save(C_cond_full, Theta_cond, tau_c, threshold_list_p, 
       file = paste0(temp_file_dir, "/GIC_local_tauc_", id_suffix, "_k", k, ".RData"))
  
  # Return num_l for Bash to capture
  cat(length(threshold_list_p$hs_vals))
  
}

GIC_step3_iterate_tau_p <- function(temp_file_dir, id_suffix, k, l){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: parallelize by working with a specific tau_p[l] value
  #
  # inputs:
  #
  # - temp_file_dir   (string)    'temp_data/simu'
  # - id_suffix       (string)    'est' or 'X_truth'
  # - k               (integer)   index of tau_c quantiles
  # - l               (integer)   index of tau_p quantiles
  #
  # ----------------------------------------------------------------------------
  
  # 1. Load both the ID-global data and the K-specific intermediate data
  load(paste0(temp_file_dir, "/GIC_local_initial_data_", id_suffix, ".RData"))
  load(paste0(temp_file_dir, "/GIC_local_tauc_", id_suffix, "_k", k, ".RData"))
  
  # 2. Apply tau_p thresholding
  tau_p <- threshold_list_p$hs_vals[l]
  excluded_indices_p <- threshold_list_p$index_path[[l]]
  
  Theta_cond_thresh <- Theta_cond
  for (idx in excluded_indices_p) { Theta_cond_thresh[[idx]][] <- 0 }
  
  # 3. Evaluate GIC
  Theta_cond_full_final <- assemble_block_matrix_irregular(Theta_cond_thresh, p)
  num_edges <- GIC_edge_count(Theta_cond_thresh, p)
  
  GIC_val <- GIC_evalulation(C_cond_full$block_matrix, Theta_cond_full_final$block_matrix, W_y, num_edges)
  
  # 4. Store result
  result_kl <- list(
    k = k, l = l,
    tau_c = tau_c, 
    tau_p = tau_p, 
    GIC = GIC_val,
    Theta_cond_thresh = Theta_cond_thresh 
  )
  
  save(result_kl, file = paste0(temp_file_dir, "/GIC_local_result_", id_suffix, "_k", k, "_l", l, ".RData"))
  
}

GIC_step2and3_serial_tau_c <- function(temp_file_dir, id_suffix, k) {
  
  
  # temp_file_dir = folder name
  # id_suffix = suffix name
  # k = tau_c index
  
  # --- Step 1: Logic from your original GIC_step2 ---
  # Load the ID-specific initial data (contains C_cond, p, W_y, threshold_list_c)
  load(paste0(temp_file_dir, "/GIC_local_initial_data_", id_suffix, ".RData"))
  
  tau_c <- threshold_list_c$hs_vals[k]
  excluded_indices <- threshold_list_c$index_path[[k]]
  
  C_cond_thresh <- C_cond
  for (idx in excluded_indices) { C_cond_thresh[[idx]][] <- 0 }
  
  C_cond_full <- assemble_block_matrix_irregular(C_cond_thresh, p)
  Theta_full  <- ginv(C_cond_full$block_matrix)
  Theta_cond  <- extract_block_matrix_irregular(Theta_full, C_cond_full$row_borders, C_cond_full$col_borders)
  
  # Get the tau_p candidates for this specific tau_c
  threshold_list_p <- GIC_get_percentile_info(Theta_cond, off_diag_indices)
  num_l <- length(threshold_list_p$hs_vals)
  
  # --- Step 2: Serial iteration over tau_p (l) ---
  best_GIC <- Inf
  best_result <- NULL
  
  for (l in 1:num_l) {
    tau_p <- threshold_list_p$hs_vals[l]
    excluded_indices_p <- threshold_list_p$index_path[[l]]
    
    Theta_cond_thresh <- Theta_cond
    for (idx in excluded_indices_p) { Theta_cond_thresh[[idx]][] <- 0 }
    
    # Evaluate GIC
    Theta_cond_full_final <- assemble_block_matrix_irregular(Theta_cond_thresh, p)
    num_edges <- GIC_edge_count(Theta_cond_thresh, p)
    
    # Calculate GIC
    current_GIC <- GIC_evalulation(C_cond_full$block_matrix, 
                                   Theta_cond_full_final$block_matrix, 
                                   W_y, 
                                   num_edges)
    
    # Keep only the best l for this k
    if (current_GIC < best_GIC) {
      best_GIC <- current_GIC
      best_result <- list(
        k = k, 
        l = l,
        tau_c = tau_c, 
        tau_p = tau_p, 
        GIC = current_GIC,
        Theta_cond_thresh = Theta_cond_thresh 
      )
    }
  }
  
  # --- Step 3: Save only the "Winner" for this k ---
  if (!is.null(best_result)) {
    save(best_result, file = paste0(temp_file_dir, "/GIC_local_best_k_", id_suffix, "_k", k, ".RData"))
  }
}

GIC_step4_finalize <- function(temp_file_dir) {
  # 1) Load the task map to know which IDs were processed
  task_map_path <- paste0(temp_file_dir, '/task_map.csv')
  
  if (!file.exists(task_map_path)) stop("Task map not found at: ", task_map_path)
  task_map <- read.csv(task_map_path)
  
  # 2) Initialize the final result list
  final_gic_results <- list()
  
  # 3) Loop through each ID (e.g., 1, 2... representing 'est', 'X_truth')
  for (i in 1:ncol(task_map)) {
    # In your Part 2/3 script, you used task_csv[1, id_suffix]
    # We follow that logic here to identify the suffix name
    id_name <- colnames(task_map)[i]
    name_entry <- task_map[1, i]
    name_i <- paste0('GIC_local_', name_entry)
    
    # Load ID-specific global data (using index i as the id_suffix)
    initial_data_path <- paste0(temp_file_dir, "/GIC_local_initial_data_", i, ".RData")
    if(!file.exists(initial_data_path)) next
    load(initial_data_path)
    
    # Identify result files - Updated pattern to match the "best_k" files
    result_files <- list.files(path = temp_file_dir, 
                               pattern = paste0("GIC_local_best_k_", i, "_k\\d+\\.RData"), 
                               full.names = TRUE)
    
    if (length(result_files) == 0) {
      cat("No result files found for suffix:", i, "\n")
      next
    }
    
    # Find best GIC among all k-winners
    lowest_GIC <- Inf
    best_data <- NULL
    df_GIC_id <- data.frame()
    
    for (f in result_files) {
      # Load 'best_result' (saved in script_GIC_local_part2and3_serial.R)
      load(f) 
      
      df_GIC_id <- rbind(df_GIC_id, data.frame(k = best_result$k, 
                                               l = best_result$l, 
                                               tau_c = best_result$tau_c, 
                                               tau_p = best_result$tau_p, 
                                               GIC = best_result$GIC))
      
      if (!is.na(best_result$GIC) && best_result$GIC < lowest_GIC) {
        lowest_GIC <- best_result$GIC
        best_data <- best_result
      }
    }
    
    if (is.null(best_data)) next
    
    # 4) Reconstruction
    # We use the threshold_list_c loaded from the initial_data_path
    C_cond_thresh_final <- C_cond
    excluded_indices_c <- threshold_list_c$index_path[[best_data$k]]
    for (idx in excluded_indices_c) {
      block <- C_cond_thresh_final[[idx]]
      C_cond_thresh_final[[idx]] <- matrix(0, nrow(block), ncol(block))
    }
    
    adj_results <- GIC_theta_to_adj(best_data$Theta_cond_thresh, p)
    num_edges_final <- GIC_edge_count(best_data$Theta_cond_thresh, p)
    
    # 5) Store
    final_gic_results[[name_i]] <- list(
      tau_c = best_data$tau_c, 
      tau_p = best_data$tau_p,
      tau_c_levels = threshold_list_c$hs_vals,
      C_cond = C_cond_thresh_final, 
      Theta_cond = best_data$Theta_cond_thresh,
      C_HS = hilbert_schmidt_norm_list_to_mat(C_cond_thresh_final, p),
      w_mat = hilbert_schmidt_norm_list_to_mat(best_data$Theta_cond_thresh, p),
      adj_mat = adj_results$adj_mat, 
      adj_list = adj_results$adj_list,
      num_edges = num_edges_final, 
      df_GIC = df_GIC_id
    )
  }
  
  # 6) Save combined results
  final_save_name <- "GIC_final_combined.RData"
  save(final_gic_results, file = paste0(temp_file_dir, "/", final_save_name))
}


# joint

GIC_joint_part1_setup <- function(C_cond_list, p, W_y_list, folder) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: gridsearch across (tau_c, tau_p) for minimum summed GIC for all n graphs
  #
  #
  # Input: 
  #
  # - C_cond_list   (list of i_j list of any sized matrix, d_i x d_j or m x m)   each item is for a y_c_query
  # - p             (scalar)
  # - W_y           (list of scalars) effective sample size
  # - folder        (string)   'temp_data/simu/GIC_global...'
  #
  # 
  # Output: 
  #
  # - list of the following metrics:
  #
  #   - tau_c   (scalar)    threshold where if your HS norm is less than this, you are zeroed out 
  #   - tau_p   (scalar) 
  #
  # ----------------------------------------------------------------------------
  
  # 0) Setup indices and norms
  n_datasets <- length(C_cond_list)
  block_names <- names(C_cond_list[[1]])
  off_diag_indices <- which(sapply(strsplit(block_names, "_"), function(x) x[1] != x[2]))
  
  C_cond_list <- lapply(C_cond_list, GIC_set_CXX_diag_identity)
  C_norms_list <- lapply(C_cond_list, function(dataset) {
    lapply(dataset, function(block) hilbert_schmidt_norm(block))
  })
  
  # 1) Tau_c Candidates
  all_c_norms <- unlist(lapply(C_cond_list, function(dataset) {
    sapply(off_diag_indices, function(idx) norm(dataset[[idx]], "F"))
  }))
  tau_c_levels <- unique(quantile(all_c_norms, probs = seq(0, 1, by = 0.01)))
  if(min(tau_c_levels) > 0) tau_c_levels <- c(0, tau_c_levels)
  
  # Save initial data for workers
  save(C_cond_list, C_norms_list, tau_c_levels, off_diag_indices, p, W_y_list, 
       file = paste0(folder, "/GIC_joint_initial.RData"))
  
  return(length(tau_c_levels))
}

GIC_joint_part2_tau_c <- function(k, folder) {
  

  # - k             (integer)  tau_c index
  # - folder        (string)   'temp_data/simu/GIC_global...'
  
  
  load(paste0(folder, "/GIC_joint_initial.RData"))
  tau_c <- tau_c_levels[k]
  n_datasets <- length(C_cond_list)
  
  current_Theta_cond_list <- list()
  current_C_full_matrices <- list()
  
  # Apply tau_c and Invert (The expensive part)
  for (i in 1:n_datasets) {
    C_thresh <- C_cond_list[[i]]
    dataset_norms <- C_norms_list[[i]]
    for (idx in off_diag_indices) {
      if (dataset_norms[[idx]] <= tau_c) C_thresh[[idx]][] <- 0
    }
    res <- assemble_block_matrix_irregular(C_thresh, p)
    current_C_full_matrices[[i]] <- res$block_matrix
    Theta_full <- ginv(res$block_matrix)
    current_Theta_cond_list[[i]] <- extract_block_matrix_irregular(Theta_full, res$row_borders, res$col_borders)
  }
  
  # Determine tau_p candidates for this specific k
  all_p_norms <- unlist(lapply(current_Theta_cond_list, function(Th) {
    sapply(off_diag_indices, function(idx) norm(Th[[idx]], "F"))
  }))
  tau_p_levels <- unique(quantile(all_p_norms, probs = seq(0, 1, by = 0.01)))
  if(length(tau_p_levels) > 0 && min(tau_p_levels) > 0) tau_p_levels <- c(0, tau_p_levels)
  
  # Save intermediate results for Part 3 workers
  save(current_Theta_cond_list, current_C_full_matrices, tau_p_levels, tau_c,
       file = paste0(folder, "/GIC_joint_k_", k, ".RData"))
  
  return(length(tau_p_levels))
}

GIC_joint_part3_evaluate <- function(k, l, folder) {
  
  
  # - k             (integer)  tau_c index
  # - l             (integer)  tau_p index
  # - folder        (string)   'temp_data/simu/GIC_global...'
  
  load(paste0(folder, "/GIC_joint_initial.RData")) # for W_y_list and p
  load(paste0(folder, "/GIC_joint_k_", k, ".RData")) # for Matrices and tau_p_levels
  
  tau_p <- tau_p_levels[l]
  total_GIC_at_pair <- 0
  n_datasets <- length(current_Theta_cond_list)
  
  for (i in 1:n_datasets) {
    Th_test <- current_Theta_cond_list[[i]]
    for (idx in off_diag_indices) {
      if (norm(Th_test[[idx]], "F") <= tau_p) Th_test[[idx]][] <- 0
    }
    
    TH_assembled <- assemble_block_matrix_irregular(Th_test, p)$block_matrix
    n_edges <- GIC_edge_count(Th_test, p)
    
    total_GIC_at_pair <- total_GIC_at_pair + 
      GIC_evalulation(current_C_full_matrices[[i]], TH_assembled, W_y_list[[i]], n_edges)
  }
  
  # Output small result file
  res <- data.frame(k = k, l = l, tau_c = tau_c, tau_p = tau_p, gic = total_GIC_at_pair)
  saveRDS(res, file = paste0(folder, "/GIC_joint_res_k", k, "_l", l, ".rds"))
}

GIC_joint_part4_finalize <- function(folder) {
  
  # - folder        (string)   'temp_data/simu/GIC_global...'
  
  # 1. Load the initial data to get original matrices and parameters
  initial_data_path <- paste0(folder, "/GIC_joint_initial.RData")
  if (!file.exists(initial_data_path)) stop("Initial data file not found.")
  load(initial_data_path) 
  # This loads: C_cond_list, C_norms_list, tau_c_levels, off_diag_indices, p, W_y_list
  
  # 2. Gather all individual GIC results from Part 3
  result_files <- list.files(path = folder, 
                             pattern = "GIC_joint_res_k\\d+_l\\d+\\.rds", 
                             full.names = TRUE)
  
  if (length(result_files) == 0) stop("No result files found in folder.")
  
  # Combine all small data frames into one master table
  all_results <- do.call(rbind, lapply(result_files, readRDS))
  
  # 3. Identify the optimal pair
  # We remove NAs just in case some inversions failed
  valid_results <- all_results[!is.na(all_results$gic), ]
  if (nrow(valid_results) == 0) stop("No valid GIC results found.")
  
  winner <- valid_results[which.min(valid_results$gic), ]
  best_tau_c <- winner$tau_c
  best_tau_p <- winner$tau_p
  lowest_gic <- winner$gic
  
  # 4. Final Reconstruction Pass
  n_datasets <- length(C_cond_list)
  final_Theta_list <- list()
  final_C_list     <- list()
  
  for (i in 1:n_datasets) {
    # 4a) Apply the optimal Global tau_c
    C_final <- C_cond_list[[i]]
    dataset_norms <- C_norms_list[[i]]
    for (idx in off_diag_indices) {
      if (dataset_norms[[idx]] <= best_tau_c) {
        C_final[[idx]][] <- 0
      }
    }
    
    # 4b) Perform the final inversion
    res_final      <- assemble_block_matrix_irregular(C_final, p)
    Theta_full_raw <- ginv(res_final$block_matrix)
    Th_cond_final  <- extract_block_matrix_irregular(Theta_full_raw, 
                                                     res_final$row_borders, 
                                                     res_final$col_borders)
    
    # 4c) Apply the optimal Global tau_p
    for (idx in off_diag_indices) {
      if (norm(Th_cond_final[[idx]], "F") <= best_tau_p) {
        Th_cond_final[[idx]][] <- 0
      }
    }
    
    final_C_list[[i]]     <- C_final
    final_Theta_list[[i]] <- Th_cond_final
  }
  
  # 5. Compile Final Output Object
  output <- list(
    joint_tau_c   = best_tau_c,
    joint_tau_p   = best_tau_p,
    total_min_GIC = lowest_gic,
    Theta_list    = final_Theta_list,
    Cond_list     = final_C_list,
    w_mat         = lapply(final_Theta_list, function(m) hilbert_schmidt_norm_list_to_mat(m, p)),
    C_HS          = lapply(final_C_list, function(m) hilbert_schmidt_norm_list_to_mat(m, p))
  )
  
  # 6. Cleanup (Optional: deletes thousands of tiny files to keep the disk clean)
  # file.remove(result_files)
  # file.remove(list.files(folder, pattern = "GIC_joint_k_.*\\.RData", full.names = TRUE))
  
  return(output)
}


