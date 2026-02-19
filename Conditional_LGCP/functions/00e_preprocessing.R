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

convert_data_for_storage <- function(LGCP_data, ID, y_c_structure, movement_num, vr_num, time_scale,
                                     time_grid_est, min_events, n_weeks, max_processes = Inf, seed = NULL){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: convert data from the mice pipeline and wraps it in a format ready for estimation
  #
  #   - re-number replicates if they are discarded due to movement and VR filtering
  #
  # 
  # input:
  #
  # - LGCP_data   (list of 3 items)
  #
  # 
  #   - [[1]] (data.frame with 'feature_id', 'time', and 'subject_num')
  #     - feature_id
  #     - time
  #     - subject_num
  #
  #   - [[2]] (nx3 data.frame with 'movement', 'VR', and 'subject_num')
  #   - [[3]] (nx3 data.frame with 'subject_num', 'age', and 'timestamp')
  #
  # - ID                (string)      mouse name like "Tau1"
  # - y_c_structure     (string)      "week_only" or "time_and_week"
  # - movement_num      (0 or 1)
  # - vr_num            (0 or 1)
  # - time_scale        (integer)     how many seconds per replicate?
  # - time_grid_est
  # - min_events        (integer)     minimum number of spikes for a replicate-process to be included
  # - n_weeks           (integer)     how many weeks do we want?
  # - max_processes     (integer)     how many neurons to look at 
  # - seed              (integer) 
  #
  # 
  # Output:
  # 
  # - output_list (list to replicate dataset)
  #
  #   - event_times (n*p-dim list)  each item is named 'k_i' is a vector of timestamps for the i-th process and k-th subject
  #   - Y_continuous
  #   - simulation_params
  #
  # ----------------------------------------------------------------------------
  
  # --- 1) Initial Extraction and Filtering ---
  p_og <- max(LGCP_data[[1]]$feature_id)
  
  # Get valid subject pool based on Movement and VR
  valid_subjects <- LGCP_data[[2]] %>% 
    filter(movement == movement_num, VR == vr_num) %>% 
    pull(subject_num)
  
  # Initial subset of the main data
  dt <- as.data.table(LGCP_data[[1]])
  dt <- dt[subject_num %in% valid_subjects & feature_id <= max_processes]
  
  # --- 2) Iterative Pruning (The While Loop) ---
  # We loop until the set of subjects and processes stabilizes
  converged <- FALSE
  
  while (!converged) {
    n_start <- nrow(dt)
    
    # A) Constraint: Minimum events per (Process x Subject)
    dt <- dt[, n_spikes := .N, by = .(feature_id, subject_num)][n_spikes >= min_events]
    dt[, n_spikes := NULL]
    
    # B) Constraint: Subject must have events on at least one process
    # (Automatically handled by data.table row removal, but we ensure subject pool is fresh)
    current_subjects <- unique(dt$subject_num)
    
    # C) Constraint: Pairwise Connectivity (Maximum Clique)
    # Build adjacency: Processes connected by shared subjects
    if (nrow(dt) > 0) {
      incidence <- table(dt$feature_id, dt$subject_num)
      incidence[incidence > 1] <- 1
      adj_matrix <- incidence %*% t(incidence)
      diag(adj_matrix) <- 0
      adj_matrix[adj_matrix > 0] <- 1
      
      # Find the largest subset of processes that are all mutually connected
      g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
      cliques <- largest_cliques(g)
      
      if (length(cliques) > 0) {
        # Keep the first largest clique found
        kept_features <- as.numeric(names(V(g)[cliques[[1]]]))
        dt <- dt[feature_id %in% kept_features]
      } else {
        dt <- dt[0] # Empty if no cliques
      }
    }
    
    # Check if any rows were removed in this iteration
    if (nrow(dt) == n_start) {
      converged <- TRUE
    }
    
    if (nrow(dt) == 0) break
  }
  
  if (nrow(dt) == 0) stop("No data left after filtering constraints.")
  
  # --- 3) Remapping and Formatting ---
  
  # Final IDs for recovery
  final_subjects <- sort(unique(dt$subject_num))
  final_features <- sort(unique(dt$feature_id))
  
  # Remap to continuous integers (1...n, 1...p)
  dt[, subject_num_map := match(subject_num, final_subjects)]
  dt[, feature_id_map := match(feature_id, final_features)]
  
  # Create event_times list: named "k_i" (Subject_Process)
  event_times <- split(dt$time, paste0(dt$subject_num_map, "_", dt$feature_id_map))
  
  # --- 4) Y_continuous processing ---
  y_cont_raw <- as.data.table(LGCP_data[[3]])[subject_num %in% final_subjects]
  # Ensure Y_continuous order matches the mapped subject_num_map
  y_cont_raw <- y_cont_raw[order(match(subject_num, final_subjects))]
  
  
    
  min_age <- min(y_cont_raw$age)
  max_age <- max(y_cont_raw$age)
  
  # 2. Check the condition
  # The number of integers in the range [min, max] is (max - min + 1)
  if (n_weeks >= (max_age - min_age + 1)) {
    # Case: Count every integer
    vals <- min_age:max_age
  } else {
    # Case: Space out the selection to get exactly n_weeks integers
    # We use round() to ensure they remain integers
    vals <- round(seq(from = min_age, to = max_age, length.out = n_weeks))
  }
  
  if (y_c_structure == 'week_only') {
    Y_continuous <- matrix(y_cont_raw$age)
    y_c_query <- matrix(vals, ncol = 1)
  } else {
    Y_continuous <- as.matrix(y_cont_raw[, .(age, timestamp)])
    max_time <- max(LGCP_data[[3]]$timestamp)
    y_c_query_time <- seq(0, max_time, by = 120)              # time is every 2 minutes
    y_c_query <- expand.grid(v1 = vals, v2 = y_c_query_time)
  }
  
  
  # --- 5) Output ---
  list(
    event_times = event_times,
    Y_continuous = Y_continuous,
    simulation_params = list(
      ID = ID,
      movement = movement_num,
      vr = vr_num,
      time_scale = time_scale,
      min_events = min_events,
      n = nrow(Y_continuous),
      p = length(final_features),
      Tmax = 1,
      query_y_cs = y_c_query,
      n_query = nrow(y_c_query),
      time_grid_est = time_grid_est,
      seed = seed
    ),
    recovery_params = list(
      p_og = p_og,
      kept_neurons = final_features,
      kept_subjects = final_subjects
    )
  )
  
}



