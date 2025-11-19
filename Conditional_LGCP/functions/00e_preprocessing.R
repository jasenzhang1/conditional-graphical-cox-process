make_time_grid <- function(m){
  
  # make m bins from 0 to 1, and create each timepoint to be in the middle of each range.
  # 
  # ex: m = 10 --> (0.05, 0.15, ..., 0.95)
  
  return((2*(1:m) - 1) / (2 * m))
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


convert_data_for_storage <- function(LGCP_data, y_c_structure, movement_num, vr_num, time_grid_est, min_events = 0, max_processes = Inf, seed = NULL){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: convert data from the mice pipeline and wraps it in a format ready for estimation
  #
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
  # - y_c_structure
  # - movement_num      (0 or 1)
  # - vr_num            (0 or 1)
  # - time_grid_est
  # - min_events        (integer)     minimum number of spikes for a replicate-process to be included
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
  

  # 0) subject_nums of this discrete strata
  
  y_d <- LGCP_data[[2]] %>% filter(movement == movement_num) %>% 
    filter(VR == vr_num) %>% 
    dplyr::pull(subject_num) %>% 
    sort()
  
  
  # 1) event_times
  
  # Ensure input is a data.table
  df <- as.data.table(LGCP_data[[1]]) %>% 
    filter(feature_id <= max_processes) %>%      # only first p-processes 
    filter(subject_num %in% y_d) %>%             # only subjects in discrete layer
    group_by(feature_id, subject_num) %>%        # only include items with >= min_events
    filter(n() >= min_events) %>%
    ungroup()
  
  y_d2 <- df$subject_num %>% unique() %>% sort()
    
  # remap the subject numbers 
  setDT(df) 
  df[, subject_num := match(subject_num, sort(unique(y_d2)))]
  
  # Split data by subject-feature combination
  event_times <- split(df$time, paste0(df$subject_num, "_", df$feature_id))

  # 2) Y_c_k and y_c_query
  
  if(y_c_structure == 'week_only'){
    Y_continuous <- LGCP_data[[3]] %>% filter(subject_num %in% y_d) %>% dplyr::pull(age) %>% matrix()

    
    y_c_query <- seq(min(LGCP_data[[3]]$age), max(LGCP_data[[3]]$age), 4) %>% matrix()
    
  } else if(y_c_structure == 'time_and_week'){
    Y_continuous <- as.matrix(LGCP_data[[3]] %>% filter(subject_num %in% y_d) %>% select(age, timestamp), ncol = 2)
    
    # y_c_query - every combination of week and time
    y_c_query_week <- seq(min(LGCP_data[[3]]$age), max(LGCP_data[[3]]$age), 4)
    
    max_time <- max(LGCP_data[[3]]$timestamp)
    last <- 120 + floor((max_time - 120) / 240) * 240
    y_c_query_time <- seq(120, last, by = 240)
    
    y_c_query <- expand.grid(v1 = y_c_query_week, v2 = y_c_query_time)
    
  }
  
  

  
  
  output_list <- list(event_times = event_times,
                      Y_continuous = Y_continuous,
                      simulation_params = list(n = nrow(Y_continuous),
                                               p = max(df$feature_id),
                                               Tmax = 1,
                                               query_y_cs = y_c_query,
                                               time_grid_est = time_grid_est,
                                               seed = seed))
  
  return(output_list)
  
}



