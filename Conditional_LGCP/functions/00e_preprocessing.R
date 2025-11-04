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
