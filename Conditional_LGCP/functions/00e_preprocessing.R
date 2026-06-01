make_time_grid <- function(m){
  
  # make m bins from 0 to 1, and create each timepoint to be in the middle of each range.
  # 
  # ex: m = 10 --> (0.05, 0.15, ..., 0.95)
  
  return((2*(1:m) - 1) / (2 * m))
}

make_y_c_grid <- function(m){
  
  # make m bins from 0 to 1 including the borders
  # 
  # ex: m = 10 --> (0, 1/9, 2/9, ... , 9/9)
  
  return(seq(0, 1, length.out = m))
}

convert_data_for_estimation <- function(subject_list, Tmax){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: convert data that was generated in simulation to a format ready for estimation. 
  #
  #   - it uses data from `dataset$subject_data`
  #
  # 
  # subject_list (list)
  # - Y_continuous
  # - X_functions
  # - precision_operator
  # - event_times
  # - event_counts
  #
  # 
  # Output:
  #
  # - df (data.frame with 'feature_id', 'time', and 'subject_num')
  #
  #   - feature_id
  #   - time
  #   - subject_num
  #
  # ----------------------------------------------------------------------------
  
  df <- extract_event_times_df(subject_list)
  colnames(df) <- c('time', 'feature_id', 'subject_num')
  
  df$time <- df$time / Tmax
  
  return(as.data.table(df))  
  
}

extract_event_times_df <- function(subject_list) {
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: helper function for convert_data_for_estimation
  # 
  #
  # input:
  #
  # - subject_list   (list of the following)
  #
  #   - Y_continuous
  #   - X_functions
  #   - precision_operator
  #   - event_times
  #   - event_counts
  # 
  # 
  # ----------------------------------------------------------------------------
  
  
  do.call(rbind, lapply(seq_along(subject_list), function(subject_id) {
    event_times <- subject_list[[subject_id]]$event_times
    
    # Handle if event_times is NULL or missing
    if (is.null(event_times)) return(NULL)
    
    do.call(rbind, lapply(seq_along(event_times), function(event_id) {
      values <- event_times[[event_id]]
      
      if (length(values) == 0) return(NULL)  # skip empty vectors
      
      data.frame(
        value = values,
        event_id = event_id,
        subject_id = subject_id
      )
    }))
  }))
}

convert_data_for_estimation_event_times <- function(event_times){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: convert data that was generated in simulation to a format ready for estimation. 
  #
  #   - it uses data from `dataset$event_times`
  #
  # 
  # input:
  #
  # - event_times (n*p-dim list)  each item is named 'k_i' is a vector of timestamps for the i-th process and k-th subject
  # 
  #
  # 
  # Output:
  #
  # - df (data.frame with 'feature_id', 'time', and 'subject_num')
  #
  #   - feature_id
  #   - time
  #   - subject_num
  #
  # ----------------------------------------------------------------------------
  

    
  df <- do.call(rbind, lapply(names(event_times), function(name) {
    # Parse the "k_i" name into subject and feature IDs
    parts <- strsplit(name, "_")[[1]]
    subject_num <- as.integer(parts[1])
    feature_id  <- as.integer(parts[2])
    
    times <- event_times[[name]]
    if (length(times) == 0) return(NULL)  # skip empties
    
    data.frame(
      feature_id = feature_id,
      time = times,
      subject_num = subject_num
    )
  }))
  
  return(as.data.table(df))
  
}

convert_data_adj_check <- function(process_ids, subject_ids){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: in preparation for bivariate estimation, there must be at least one subject_ID connecting each pair of process IDs
  # 
  #       so we create a maximal clique by deleting processes that do not connect with everyone
  # 
  #
  # inputs:
  #
  # - process_ids   (vector)   vector of process ID's 
  # - subject_ids   (vector)   vector of subject ID's 
  #
  # outputs:
  #
  # - nodes_to_keep  (vector)    which vertices to keep of process_ids
  #
  # ----------------------------------------------------------------------------
  
  # part 1) obtain an adjacency matrix
  
  # 1. Create the incidence matrix (Binary: Process vs Subject)
  incidence_matrix <- table(process_ids, subject_ids)
  incidence_matrix[incidence_matrix > 1] <- 1  # Ensure it is binary
  
  # 2. Matrix Multiplication (P x S) * (S x P) = (P x P)
  adj_matrix <- incidence_matrix %*% t(incidence_matrix)
  
  # 3. Final touch: Binary adjacency (1 if shared, 0 otherwise)
  adj_matrix[adj_matrix > 0] <- 1
  

  diag(adj_matrix) <- 0
  

  # part 2) choose which vertices to delete, if any 
  
  # 1. Create the graph from your adjacency matrix
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  
  # 2. Find the Maximum Clique (the largest fully connected subset)
  max_clique_list <- largest_cliques(g)
  
  # 3. Get the names of the processes to KEEP
  nodes_to_keep <- names(V(g)[max_clique_list[[1]]])
  
  # 4. Identify which to DELETE
  all_nodes <- V(g)$name
  nodes_to_delete <- setdiff(all_nodes, nodes_to_keep)
  
  print(paste("Keep:", paste(nodes_to_keep, collapse=", ")))
  print(paste("Delete:", paste(nodes_to_delete, collapse=", ")))
  
  return(sort(as.numeric(nodes_to_keep)))
  
}

convert_data_for_storage <- function(LGCP_data, df_brain_region, ID, y_c_structure, movement_num, vr_num, region, time_scale,
                                     time_grid_est, min_events, n_weeks, max_processes = Inf, seed = NULL){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: Convert data from the mice pipeline and wrap it in a format ready
  #       for estimation.
  #
  #   - Re-numbers replicates if they are discarded due to movement and VR
  #     filtering.
  #   - Runs iterative pruning to ensure a valid, fully-connected neuron set
  #     (maximum clique) for the requested strata.
  #
  # input:
  #
  # - LGCP_data   (list of 3 items)
  #
  #   - [[1]] (data.frame with 'feature_id', 'time', and 'subject_num')
  #     - feature_id
  #     - time
  #     - subject_num
  #
  #   - [[2]] (nx3 data.frame with 'movement', 'VR', and 'subject_num')
  #   - [[3]] (nx3 data.frame with 'subject_num', 'age', and 'timestamp')
  #
  # - df_brain_region  (dataframe of)
  #
  #   - Neuron_Num      (integer)  i = 1, ..., p
  #   - Electrode_Num   (integer)  1 through 64
  #   - Brain_Region    (factor)   'Hippocampus' or 'Entorhinal_Cortex'
  #   - Mouse           (string)   '346' mouse ID in string form
  #   - Strain          (string)   'Tau' or 'WT'
  #   - ID2             (factor)   'Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3'
  #
  # - ID                (string)      mouse name like "Tau1"
  # - y_c_structure     (string)      "week_only" or "time_and_week"
  # - movement_num      (0 or 1)      movement filter for the target stratum
  # - vr_num            (0 or 1)      VR filter for the target stratum
  # - region            (string)      brain region selector:
  #                                     'HIP'               — hippocampus only
  #                                     'EHC'               — entorhinal cortex only
  #                                     'BOTH_100'          — top 50 per region by total spikes
  #                                     'BOTH_150'          — top 75 per region by total spikes
  #                                     'BOTH_100_NORMALIZED' — top 50 per region by spikes
  #                                       within the VR-on union (m0vr1 + m1vr1), subject to
  #                                       the constraint that every neuron pair has at least one
  #                                       co-active replicate in BOTH m0vr1 AND m1vr1 strata.
  #                                       A greedy maximum-clique approach is used jointly on
  #                                       both strata before trimming to the top 50 per region.
  # - time_scale        (integer)     seconds per replicate
  # - time_grid_est
  # - min_events        (integer)     minimum spikes for a (replicate, process) to be included
  # - n_weeks           (integer)     how many weeks to query
  # - max_processes     (integer)     cap on number of neurons (applied after region selection)
  # - seed              (integer)
  #
  # output:
  #
  # - output_list (list)
  #
  #   - event_times        (n*p-dim list)  each item named 'k_i' is a vector of
  #                                        timestamps for the i-th process and k-th subject
  #   - Y_continuous       (matrix)        continuous covariate matrix
  #   - simulation_params  (list)          run metadata and query grid
  #   - recovery_params    (list)          original neuron IDs, regions, and subjects kept
  #
  # ----------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Helper: given a set of subject_nums and feature_ids from LGCP_data[[1]],
  # run the iterative pruning loop (min_events + maximum clique) and return
  # the surviving data.table. Used both by the standard path and by the
  # BOTH_100_NORMALIZED joint-strata path.
  # --------------------------------------------------------------------------
  run_pruning <- function(dt_input) {
    
    dt        <- copy(dt_input)
    converged <- FALSE
    
    while (!converged) {
      
      n_start <- nrow(dt)
      
      # A) Minimum spikes per (process, replicate)
      dt <- dt[, n_spikes := .N, by = .(feature_id, subject_num)][n_spikes >= min_events]
      dt[, n_spikes := NULL]
      
      # B) Pairwise connectivity: keep the largest clique of processes such
      #    that every pair shares at least one common replicate. This ensures
      #    non-degenerate bivariate intensity estimates for all neuron pairs.
      if (nrow(dt) > 0) {
        incidence            <- table(dt$feature_id, dt$subject_num)
        incidence[incidence > 1] <- 1
        adj_matrix           <- incidence %*% t(incidence)
        diag(adj_matrix)     <- 0
        adj_matrix[adj_matrix > 0] <- 1
        
        g       <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
        cliques <- largest_cliques(g)
        
        if (length(cliques) > 0) {
          # Among all largest cliques, keep the one that retains the most rows
          clique_features <- lapply(cliques, function(cl) as.numeric(V(g)$name[cl]))
          best_idx        <- which.max(sapply(clique_features, function(feats) nrow(dt[feature_id %in% feats])))
          dt              <- dt[feature_id %in% clique_features[[best_idx]]]
        } else {
          dt <- dt[0]
        }
      }
      
      print(paste0('(4/6) After pruning, rows: ',      nrow(dt)))
      print(paste0('(5/6) After pruning, subjects: ',  length(unique(dt$subject_num))))
      print(paste0('(6/6) After pruning, processes: ', length(unique(dt$feature_id))))
      
      if (nrow(dt) == n_start) converged <- TRUE
      if (nrow(dt) == 0)       break
    }
    
    dt
  }
  
  # --- 1) Initial Extraction and Filtering ---
  p_og <- max(LGCP_data[[1]]$feature_id)
  
  if (region == 'HIP') {
    
    relevant_neurons <- df_brain_region$Neuron_Num[
      df_brain_region$ID2 == ID & df_brain_region$Brain_Region == 'Hippocampus']
    
  } else if (region == 'EHC') {
    
    relevant_neurons <- df_brain_region$Neuron_Num[
      df_brain_region$ID2 == ID & df_brain_region$Brain_Region == 'Entorhinal_Cortex']
    
  } else if (region == 'BOTH_100') {
    
    # Top 50 neurons per region by total spike count across all replicates.
    # Does not guarantee exactly 100 neurons — pruning may reduce this further.
    relevant_neurons <- do.call(c, lapply(c("Hippocampus", "Entorhinal_Cortex"), function(br) {
      ids <- df_brain_region$Neuron_Num[df_brain_region$ID2 == ID & df_brain_region$Brain_Region == br]
      tbl <- sort(table(LGCP_data[[1]]$feature_id[LGCP_data[[1]]$feature_id %in% ids]), decreasing = TRUE)
      as.integer(names(head(tbl, 50)))
    })) %>% sort()
    
  } else if (region == 'BOTH_150') {
    
    # Top 75 neurons per region by total spike count across all replicates.
    # Does not guarantee exactly 150 neurons — pruning may reduce this further.
    relevant_neurons <- do.call(c, lapply(c("Hippocampus", "Entorhinal_Cortex"), function(br) {
      ids <- df_brain_region$Neuron_Num[df_brain_region$ID2 == ID & df_brain_region$Brain_Region == br]
      tbl <- sort(table(LGCP_data[[1]]$feature_id[LGCP_data[[1]]$feature_id %in% ids]), decreasing = TRUE)
      as.integer(names(head(tbl, 75)))
    })) %>% sort()
    
  } else if (region == 'BOTH_100_NORMALIZED') {
    
    # -----------------------------------------------------------------------
    # BOTH_100_NORMALIZED: Joint clique across m0vr1 and m1vr1, then trim to
    # top 50 per region by VR-on spike counts.
    #
    # Strategy:
    #   1. Pre-filter to top 75 neurons per region by total spike count across
    #      all replicates (mirroring BOTH_150 selection).
    #   2. Pool all replicates where VR == 1 (union of m0vr1 and m1vr1) to
    #      rank neurons by total spike activity under VR-on conditions.
    #   3. Separately run the iterative pruning (min_events + max clique) on
    #      the m0vr1 and m1vr1 strata independently, restricted to the top-75
    #      candidate set, obtaining neurons non-degenerate within each stratum.
    #   4. Take the intersection of surviving neuron sets across both strata —
    #      this is the joint feasible set where all pairs are non-degenerate
    #      in BOTH conditions simultaneously.
    #   5. From that joint feasible set, select the top 50 neurons per region
    #      ranked by VR-on spike counts (step 2).
    # -----------------------------------------------------------------------
    
    # Step 1: Top 75 per region by total spike count across all replicates
    top75_neurons <- do.call(c, lapply(c("Hippocampus", "Entorhinal_Cortex"), function(br) {
      ids <- df_brain_region$Neuron_Num[df_brain_region$ID2 == ID & df_brain_region$Brain_Region == br]
      tbl <- sort(table(LGCP_data[[1]]$feature_id[LGCP_data[[1]]$feature_id %in% ids]), decreasing = TRUE)
      as.integer(names(head(tbl, 75)))
    }))
    message(sprintf("  [BOTH_100_NORMALIZED] Top-75-per-region candidate set size: %d", length(top75_neurons)))
    
    # Step 2: Identify all VR-on replicates (movement = 0 or 1, VR = 1)
    vr_on_subjects <- LGCP_data[[2]] %>%
      filter(VR == 1) %>%
      pull(subject_num)
    
    # Rank all neurons by spike count within VR-on replicates
    dt_vr_on <- as.data.table(LGCP_data[[1]])[subject_num %in% vr_on_subjects]
    vr_on_spike_counts <- sort(table(dt_vr_on$feature_id), decreasing = TRUE)
    
    # Step 3: Run pruning independently for m0vr1 and m1vr1, restricted to top75
    prune_stratum <- function(mov, vr) {
      subj <- LGCP_data[[2]] %>% filter(movement == mov, VR == vr) %>% pull(subject_num)
      dt_s <- as.data.table(LGCP_data[[1]])[subject_num %in% subj & feature_id %in% top75_neurons]
      message(sprintf("  [BOTH_100_NORMALIZED] Pruning stratum m%dvr%d: %d subjects, %d neurons before pruning",
                      mov, vr, length(subj), length(unique(dt_s$feature_id))))
      run_pruning(dt_s)
    }
    
    dt_m0vr1 <- prune_stratum(0, 1)
    dt_m1vr1 <- prune_stratum(1, 1)
    
    neurons_m0vr1 <- unique(dt_m0vr1$feature_id)
    neurons_m1vr1 <- unique(dt_m1vr1$feature_id)
    
    # Step 4: Joint feasible set — neurons that survive pruning in BOTH strata
    joint_feasible <- intersect(neurons_m0vr1, neurons_m1vr1)
    message(sprintf("  [BOTH_100_NORMALIZED] Joint feasible neuron set size: %d", length(joint_feasible)))
    
    if (length(joint_feasible) == 0) {
      stop("BOTH_100_NORMALIZED: No neurons survived joint pruning across m0vr1 and m1vr1.")
    }
    
    # Step 5: From the joint feasible set, pick top 50 per region by VR-on spikes
    relevant_neurons <- do.call(c, lapply(c("Hippocampus", "Entorhinal_Cortex"), function(br) {
      
      # Neurons in this region that are in the joint feasible set
      region_ids <- df_brain_region$Neuron_Num[
        df_brain_region$ID2 == ID & df_brain_region$Brain_Region == br]
      candidates <- intersect(joint_feasible, region_ids)
      message(sprintf("  [BOTH_100_NORMALIZED] %s: %d joint-feasible candidates before top-50 trim.", br, length(candidates)))
      
      # Rank by VR-on spike counts; neurons absent from vr_on_spike_counts get 0
      counts     <- as.integer(vr_on_spike_counts[as.character(candidates)])
      counts[is.na(counts)] <- 0L
      ranked     <- candidates[order(counts, decreasing = TRUE)]
      
      as.integer(head(ranked, 50))
    })) %>% sort()
    
    message(sprintf("  [BOTH_100_NORMALIZED] Selected %d neurons after trimming to top 50 per region.",
                    length(relevant_neurons)))
    
  } else {
    
    # Fallback: all neurons for this mouse
    relevant_neurons <- df_brain_region$Neuron_Num[df_brain_region$ID2 == ID]
  }
  
  print(paste0('Num relevant neurons: ', length(relevant_neurons)))
  
  # Apply max_processes cap (by sorted neuron index)
  n_temp           <- min(length(relevant_neurons), max_processes)
  relevant_neurons <- sort(relevant_neurons)[1:n_temp]
  
  print(paste0('Neurons after max_processes: ', length(relevant_neurons)))
  
  # --- 2) Filter to target stratum and relevant neurons ---
  valid_subjects <- LGCP_data[[2]] %>%
    filter(movement == movement_num, VR == vr_num) %>%
    pull(subject_num)
  
  print(paste0('Num subjects: ', length(valid_subjects)))
  
  dt <- as.data.table(LGCP_data[[1]])
  dt <- dt[subject_num %in% valid_subjects & feature_id %in% relevant_neurons]
  
  # --- 3) Iterative Pruning for the target stratum ---
  # For BOTH_100_NORMALIZED the relevant_neurons are already guaranteed to
  # survive pruning in BOTH strata, but we still run pruning here to enforce
  # min_events within the specific target stratum and to re-check the clique
  # after the max_processes cap may have reduced the neuron set.
  dt <- run_pruning(dt)
  
  if (nrow(dt) == 0) stop("No data left after filtering constraints.")
  
  # --- 4) Remapping and Formatting ---
  
  # Final IDs for recovery
  final_subjects     <- sort(unique(dt$subject_num))
  final_features     <- sort(unique(dt$feature_id))
  final_brain_region <- df_brain_region$Brain_Region[
    df_brain_region$ID2 == ID & df_brain_region$Neuron_Num %in% final_features]
  
  # Remap to continuous integers (1...n, 1...p)
  dt[, subject_num_map := match(subject_num, final_subjects)]
  dt[, feature_id_map  := match(feature_id,  final_features)]
  
  # Create event_times list: each element named "k_i" (Subject_Process)
  event_times <- split(dt$time, paste0(dt$subject_num_map, "_", dt$feature_id_map))
  
  # --- 5) Y_continuous processing ---
  y_cont_raw <- as.data.table(LGCP_data[[3]])[subject_num %in% final_subjects]
  # Ensure ordering matches the mapped subject_num_map
  y_cont_raw <- y_cont_raw[order(match(subject_num, final_subjects))]
  
  min_age <- min(y_cont_raw$age)
  max_age <- max(y_cont_raw$age)
  
  # Build the week query grid: either every integer in range or evenly spaced
  if (n_weeks >= (max_age - min_age + 1)) {
    vals <- min_age:max_age
  } else {
    vals <- round(seq(from = min_age, to = max_age, length.out = n_weeks))
  }
  
  if (y_c_structure == 'week_only') {
    Y_continuous <- matrix(y_cont_raw$age)
    y_c_query    <- matrix(vals, ncol = 1)
  } else {
    Y_continuous      <- as.matrix(y_cont_raw[, .(age, timestamp)])
    max_time          <- max(LGCP_data[[3]]$timestamp)
    y_c_query_time    <- seq(0, max_time, by = 120)   # one query every 2 minutes
    y_c_query         <- expand.grid(v1 = vals, v2 = y_c_query_time)
  }
  
  # --- 6) Output ---
  list(
    event_times = event_times,
    Y_continuous = Y_continuous,
    simulation_params = list(
      ID         = ID,
      movement   = movement_num,
      vr         = vr_num,
      time_scale = time_scale,
      min_events = min_events,
      n          = nrow(Y_continuous),
      p          = length(final_features),
      Tmax       = 1,
      query_y_cs = y_c_query,
      n_query    = nrow(y_c_query),
      time_grid_est = time_grid_est,
      seed       = seed
    ),
    recovery_params = list(
      p_og                = p_og,
      kept_neurons        = final_features,
      kept_neuron_regions = final_brain_region,
      kept_subjects       = final_subjects
    )
  )
}



