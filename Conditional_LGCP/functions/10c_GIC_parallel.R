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

GIC_step2and3_serial_tau_c <- function(temp_file_dir, id_suffix, k, min_connect_pct) {
  
  # temp_file_dir = folder name
  # id_suffix = suffix name
  # k = tau_c index
  # min_connect = min edges (0 for nothing, 1 for something)
  
  # --- Step 1: Logic from your original GIC_step2 ---
  # Load the ID-specific initial data (contains C_cond, p, W_y, threshold_list_c)
  load(paste0(temp_file_dir, "/GIC_local_initial_data_", id_suffix, ".RData"))
  
  # convert pct to minimum edge count
  max_edges <- p * (p - 1) / 2
  min_connect <- ceiling(min_connect_pct * max_edges)
  
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
    passes_connect <- (min_connect == 0) || (num_edges >= min_connect)
    
    if (current_GIC < best_GIC && passes_connect) {
      best_GIC <- current_GIC
      best_result <- list(
        k = k, 
        l = l,
        tau_c = tau_c, 
        tau_p = tau_p,
        GIC = current_GIC,
        num_edges = num_edges,
        Theta_cond_thresh = Theta_cond_thresh
      )
    }
  }
  
  if (is.null(best_result)) {
    cat(sprintf("[GIC WARNING] id_suffix=%s, k=%d: no valid l found across %d tau_p values with min_connect=%d\n",
                id_suffix, k, num_l, min_connect))
    cat(sprintf("[GIC WARNING] edge counts across l: %s\n",
                paste(sapply(1:num_l, function(l) {
                  Theta_cond_thresh <- Theta_cond
                  for (idx in threshold_list_p$index_path[[l]]) { Theta_cond_thresh[[idx]][] <- 0 }
                  GIC_edge_count(Theta_cond_thresh, p)
                }), collapse = ", ")))
  } else {
    cat(sprintf("[GIC] id_suffix=%s, k=%s: best l=%s, num_edges=%s, GIC=%.4f\n",
                as.character(id_suffix), as.character(k), 
                as.character(best_result$l), as.character(best_result$num_edges), 
                best_result$GIC))
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
  for (i in 1:nrow(task_map)) {
    
    name_entry <- task_map[i, 1]
    name_i <- paste0('GIC_KL_', name_entry)
    
    cat(sprintf("[GIC step4] Processing suffix: %s\n", name_entry))
    
    initial_data_path <- paste0(temp_file_dir, "/GIC_local_initial_data_", name_entry, ".RData")
    if (!file.exists(initial_data_path)) {
      cat(sprintf("[GIC step4] WARNING: initial data file not found for suffix=%s, skipping\n", name_entry))
      next
    }
    load(initial_data_path)
    
    # Identify result files - Updated pattern to match the "best_k" files
    result_files <- list.files(path = temp_file_dir,
                               pattern = paste0("GIC_local_best_k_", name_entry, "_k\\d+\\.RData"),
                               full.names = TRUE)
    
    cat(sprintf("[GIC step4] Found %d result files for suffix=%s\n", length(result_files), name_entry))
    
    if (length(result_files) == 0) {
      cat(sprintf("[GIC step4] WARNING: no best_k files found for suffix=%s - likely all k's had too few edges\n", name_entry))
      next
    }
    
    # Find best GIC among all k-winners
    lowest_GIC <- Inf
    best_data <- NULL
    df_GIC_id <- data.frame()
    
    for (f in result_files) {
      # Load 'best_result' (saved in script_GIC_local_part2and3_serial.R)
      load(f)
      
      cat(sprintf("[GIC step4] Loaded %s: k=%s, l=%s, num_edges=%s, GIC=%.4f\n",
                  basename(f), best_result$k, best_result$l, best_result$num_edges, best_result$GIC))
      
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
    
    if (is.null(best_data)) {
      cat(sprintf("[GIC step4] WARNING: best_data is NULL for suffix=%s after loading all files\n", name_entry))
      next
    }
    
    cat(sprintf("[GIC step4] Best for suffix=%s: k=%s, l=%s, num_edges=%s, GIC=%.4f\n",
                name_entry, best_data$k, best_data$l, best_data$num_edges, best_data$GIC))
    
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
  cat(sprintf("[GIC step4] final_gic_results has %d entries: %s\n",
              length(final_gic_results),
              paste(names(final_gic_results), collapse = ", ")))
  
  if (length(final_gic_results) == 0) {
    cat("[GIC step4] WARNING: final_gic_results is empty - GIC_final_combined.RData will be empty\n")
  }
  
  final_save_name <- "GIC_final_combined.RData"
  save(final_gic_results, file = paste0(temp_file_dir, "/", final_save_name))
}

# joint

GIC_joint_part1_setup <- function(C_cond_list, p, W_y_list, folder, id_suffix) {
  
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
  # - id_suffix     (string) (e.g., "est_eig1") to distinguish different estimation runs
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
  
  # Identify off-diagonal blocks (e.g., "1_2", "2_3") vs diagonal blocks ("1_1")
  off_diag_indices <- which(sapply(strsplit(block_names, "_"), function(x) x[1] != x[2]))
  
  # Ensure diagonal blocks are identity for GIC stability
  C_cond_list <- lapply(C_cond_list, GIC_set_CXX_diag_identity)
  
  # Pre-calculate HS norms for every block in every dataset
  C_norms_list <- lapply(C_cond_list, function(dataset) {
    lapply(dataset, function(block) hilbert_schmidt_norm(block))
  })
  
  # 1) Tau_c Candidate Generation (Aggregated across all datasets/queries)
  all_c_norms <- unlist(lapply(C_cond_list, function(dataset) {
    sapply(off_diag_indices, function(idx) norm(dataset[[idx]], "F"))
  }))
  
  # Create a grid based on percentiles of all observed off-diagonal norms
  tau_c_levels <- unique(quantile(all_c_norms, probs = seq(0, 1, by = 0.01)))
  
  # Ensure 0 is included as a candidate (no thresholding)
  if(length(tau_c_levels) == 0 || min(tau_c_levels) > 0) {
    tau_c_levels <- c(0, tau_c_levels)
  }
  
  # 2) Save Setup Data
  # Filename includes 'id' so that part 2/3 workers know which estimation type they are solving
  save_path <- paste0(folder, "/GIC_joint_initial_", id_suffix, ".RData")
  
  save(C_cond_list, 
       C_norms_list, 
       tau_c_levels, 
       off_diag_indices, 
       p, 
       W_y_list, 
       id_suffix,
       file = save_path)
  
  # Return the number of tau_c candidates for the Bash loop (max_k)
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

GIC_joint_part2and3_serialized <- function(k, suffix_name, folder) {
  
  # ----------------------------------------------------------------------------
  # GOAL: For a fixed tau_c (k) and estimation type (suffix_name), evaluate ALL tau_p levels.
  # This avoids reloading/re-inverting the large C matrices multiple times.
  # ----------------------------------------------------------------------------
  
  
  # 1. Load Initial Setup Data for this specific ID
  initial_file <- paste0(folder, "/GIC_joint_initial_", suffix_name, ".RData")
  if (!file.exists(initial_file)) stop("Initial data not found for: ", suffix_name)
  load(initial_file) 
  # Loads: C_cond_list, C_norms_list, tau_c_levels, off_diag_indices, p, W_y_list
  
  tau_c <- tau_c_levels[k]
  n_datasets <- length(C_cond_list)
  
  # 2. Part 2 Logic: Apply tau_c and Invert (The expensive step)
  # We do this once per k index.
  current_Theta_cond_list <- list()
  current_C_full_matrices <- list()
  
  for (i in 1:n_datasets) {
    C_thresh <- C_cond_list[[i]]
    dataset_norms <- C_norms_list[[i]]
    
    # Threshold C blocks
    for (idx in off_diag_indices) {
      if (dataset_norms[[idx]] <= tau_c) C_thresh[[idx]][] <- 0
    }
    
    # Assemble and Invert
    res <- assemble_block_matrix_irregular(C_thresh, p)
    current_C_full_matrices[[i]] <- res$block_matrix
    
    # Generate Theta (the precision matrix estimate)
    Theta_full <- ginv(res$block_matrix)
    current_Theta_cond_list[[i]] <- extract_block_matrix_irregular(
      Theta_full, res$row_borders, res$col_borders
    )
  }
  
  # 3. Determine tau_p candidates for this specific inversion state
  all_p_norms <- unlist(lapply(current_Theta_cond_list, function(Th) {
    sapply(off_diag_indices, function(idx) norm(Th[[idx]], "F"))
  }))
  tau_p_levels <- unique(quantile(all_p_norms, probs = seq(0, 1, by = 0.01)))
  if(length(tau_p_levels) > 0 && min(tau_p_levels) > 0) tau_p_levels <- c(0, tau_p_levels)
  
  # 4. Part 3 Logic: Loop through all tau_p (l) in memory
  k_results <- lapply(seq_along(tau_p_levels), function(l) {
    tau_p <- tau_p_levels[l]
    total_GIC_at_pair <- 0
    
    for (i in 1:n_datasets) {
      Th_test <- current_Theta_cond_list[[i]]
      
      # Apply tau_p thresholding to Theta blocks
      for (idx in off_diag_indices) {
        if (norm(Th_test[[idx]], "F") <= tau_p) Th_test[[idx]][] <- 0
      }
      
      # Re-assemble for GIC evaluation
      TH_assembled <- assemble_block_matrix_irregular(Th_test, p)$block_matrix
      n_edges <- GIC_edge_count(Th_test, p)
      
      # Calculate Joint GIC (summed across all queries/datasets)
      total_GIC_at_pair <- total_GIC_at_pair + 
        GIC_evalulation(current_C_full_matrices[[i]], TH_assembled, W_y_list[[i]], n_edges)
    }
    
    return(data.frame(suffix_name = suffix_name, k = k, l = l, tau_c = tau_c, tau_p = tau_p, gic = total_GIC_at_pair))
  })
  
  # 5. Save consolidated result file for this k and id
  combined_res <- do.call(rbind, k_results)
  save_path <- paste0(folder, "/GIC_joint_res_", suffix_name, "_k", k, ".rds")
  saveRDS(combined_res, file = save_path)
  
}

GIC_joint_part4_finalize <- function(GIC_folder, temp_file_dir, setting_info_list, cont_inds, mouse) {
  
  
  # 1) Setup original file paths for infusion
  list2env(setting_info_list, envir = environment())
  if(mouse){
    step_3_paths <- paste0(temp_file_dir, '/part3_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
  } else {
    step_3_paths <- paste0(temp_file_dir, '/part3_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_rep_', rep_i, '.rds')
  }
  original_results <- lapply(step_3_paths, readRDS)
  
  # 2) Load task map
  task_map <- read.csv(paste0(GIC_folder, '/task_map.csv'))
  
  for (i in 1:nrow(task_map)) {
    current_id <- as.character(task_map$id[i])
    suffix <- sub("^KL_cor_", "", current_id) 
    
    initial_data_path <- paste0(GIC_folder, "/GIC_joint_initial_", current_id, ".RData")
    if(!file.exists(initial_data_path)) next
    load(initial_data_path) 
    
    result_files <- list.files(path = GIC_folder, pattern = paste0("GIC_joint_res_", current_id, "_k\\d+\\.rds"), full.names = TRUE)
    if (length(result_files) == 0) next
    
    all_res_df <- do.call(rbind, lapply(result_files, readRDS))
    winner     <- all_res_df[which.min(all_res_df$gic), ]
    
    # 3) Reconstruction & Infusion Loop
    for (d in 1:cont_inds) {
      C_d <- C_cond_list[[d]]
      for (idx in off_diag_indices) {
        if (C_norms_list[[d]][[idx]] <= winner$tau_c) C_d[[idx]][] <- 0
      }
      res_inv    <- assemble_block_matrix_irregular(C_d, p)
      Theta_raw  <- ginv(res_inv$block_matrix)
      Theta_cond <- extract_block_matrix_irregular(Theta_raw, res_inv$row_borders, res_inv$col_borders)
      for (idx in off_diag_indices) {
        if (norm(Theta_cond[[idx]], "F") <= winner$tau_p) Theta_cond[[idx]][] <- 0
      }
      
      # Inject into original result structure
      original_results[[d]]$step_11[[paste0('w_mat_KL_GIC_global_', suffix)]]  <- hilbert_schmidt_norm_list_to_mat(Theta_cond, p)
      original_results[[d]]$step_11b[[paste0('C_HS_KL_GIC_global_', suffix)]] <- hilbert_schmidt_norm_list_to_mat(C_d, p)
      original_results[[d]]$step_11x[[paste0('tau_c_global_', suffix)]]       <- winner$tau_c
      original_results[[d]]$step_11y[[paste0('tau_p_global_', suffix)]]       <- winner$tau_p
    }
  }
  
  # 4) Save back to RDS
  for (d in 1:cont_inds) saveRDS(original_results[[d]], file = step_3_paths[d])
}

# hybrid can join in with joint

GIC_hybrid_part2and3_serialized <- function(k, suffix_name, folder) {
  
  
  # 1. Load Setup Data
  load(paste0(folder, "/GIC_joint_initial_", suffix_name, ".RData"))
  tau_c <- tau_c_levels[k]
  n_datasets <- length(C_cond_list)
  
  # 2. Apply Global tau_c and Invert
  current_Theta_list <- list()
  current_C_full_list <- list()
  res_list <- list()
  
  for (i in 1:n_datasets) {
    C_thresh <- C_cond_list[[i]]
    for (idx in off_diag_indices) {
      if (C_norms_list[[i]][[idx]] <= tau_c) C_thresh[[idx]][] <- 0
    }
    res <- assemble_block_matrix_irregular(C_thresh, p)
    
    res_list[[i]] <- res
    current_C_full_list[[i]] <- res$block_matrix
    current_Theta_list[[i]]  <- ginv(res$block_matrix)
  }
  
  # 3. Local Search for best tau_p per dataset
  hybrid_total_GIC <- 0
  local_winners <- list()
  
  for (i in 1:n_datasets) {
    # Generate tau_p candidates for THIS specific matrix
    Th_cond <- extract_block_matrix_irregular(current_Theta_list[[i]], res_list[[i]]$row_borders, res_list[[i]]$col_borders)
    p_norms <- sapply(off_diag_indices, function(idx) norm(Th_cond[[idx]], "F"))
    tau_p_levels <- unique(quantile(p_norms, probs = seq(0, 1, by = 0.01)))
    
    best_local_gic <- Inf
    best_local_tau_p <- 0
    
    for (tp in tau_p_levels) {
      # Apply tp locally
      Th_test <- Th_cond
      for (idx in off_diag_indices) {
        if (norm(Th_test[[idx]], "F") <= tp) Th_test[[idx]][] <- 0
      }
      
      # Evaluate GIC for this dataset only
      TH_assembled <- assemble_block_matrix_irregular(Th_test, p)$block_matrix
      n_edges <- GIC_edge_count(Th_test, p)
      current_gic <- GIC_evalulation(current_C_full_list[[i]], TH_assembled, W_y_list[[i]], n_edges)
      
      if (current_gic < best_local_gic) {
        best_local_gic <- current_gic
        best_local_tau_p <- tp
      }
    }
    
    hybrid_total_GIC <- hybrid_total_GIC + best_local_gic
    local_winners[[i]] <- list(tau_p = best_local_tau_p, gic = best_local_gic)
  }
  
  # 4. Save result for this global tau_c
  res <- list(k = k, tau_c = tau_c, total_gic = hybrid_total_GIC, locals = local_winners)
  saveRDS(res, file = paste0(folder, "/GIC_hybrid_res_", suffix_name, "_k", k, ".rds"))
  
}

GIC_hybrid_part4_finalize <- function(GIC_folder, temp_file_dir, setting_info_list, cont_inds, mouse) {
  # 1) Setup original file paths for infusion
  list2env(setting_info_list, envir = environment())
  if(mouse){
    step_3_paths <- paste0(temp_file_dir, '/part3_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
  } else {
    step_3_paths <- paste0(temp_file_dir, '/part3_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_rep_', rep_i, '.rds')
  }
  original_results <- lapply(step_3_paths, readRDS)
  
  # 2) Load task map
  task_map <- read.csv(paste0(GIC_folder, '/task_map.csv'))
  
  for (i in 1:nrow(task_map)) {
    current_id <- as.character(task_map$id[i])
    suffix <- sub("^KL_cor_", "", current_id) 
    
    initial_data_path <- paste0(GIC_folder, "/GIC_joint_initial_", current_id, ".RData")
    if(!file.exists(initial_data_path)) next
    load(initial_data_path) 
    
    result_files <- list.files(path = GIC_folder, pattern = paste0("GIC_hybrid_res_", current_id, "_k\\d+\\.rds"), full.names = TRUE)
    if (length(result_files) == 0) next
    
    best_overall_gic <- Inf
    winner_file <- NULL
    for (f in result_files) {
      tmp <- readRDS(f)
      if (!is.na(tmp$total_gic) && tmp$total_gic < best_overall_gic) {
        best_overall_gic <- tmp$total_gic
        winner_file <- tmp
      }
    }
    if (is.null(winner_file)) next
    
    # 3) Reconstruction & Infusion Loop
    for (d in 1:cont_inds) {
      C_d <- C_cond_list[[d]]
      for (idx in off_diag_indices) {
        if (C_norms_list[[d]][[idx]] <= winner_file$tau_c) C_d[[idx]][] <- 0
      }
      res_inv    <- assemble_block_matrix_irregular(C_d, p)
      Theta_raw  <- ginv(res_inv$block_matrix)
      Theta_cond <- extract_block_matrix_irregular(Theta_raw, res_inv$row_borders, res_inv$col_borders)
      
      # Use local tau_p for this specific dataset
      local_tp_d <- winner_file$locals[[d]]$tau_p
      for (idx in off_diag_indices) {
        if (norm(Theta_cond[[idx]], "F") <= local_tp_d) Theta_cond[[idx]][] <- 0
      }
      
      # Inject into original result structure
      original_results[[d]]$step_11[[paste0('w_mat_KL_GIC_hybrid_', suffix)]]  <- hilbert_schmidt_norm_list_to_mat(Theta_cond, p)
      original_results[[d]]$step_11b[[paste0('C_HS_KL_GIC_hybrid_', suffix)]] <- hilbert_schmidt_norm_list_to_mat(C_d, p)
      original_results[[d]]$step_11x[[paste0('tau_c_hybrid_', suffix)]]       <- winner_file$tau_c
      original_results[[d]]$step_11y[[paste0('tau_p_hybrid_', suffix)]]       <- local_tp_d
    }
  }
  
  # 4) Save back to RDS
  for (d in 1:cont_inds) saveRDS(original_results[[d]], file = step_3_paths[d])
}

