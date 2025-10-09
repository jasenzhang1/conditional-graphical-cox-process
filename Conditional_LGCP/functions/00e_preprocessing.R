convert_data_for_estimation <- function(subject_list, Tmax){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: convert data that was generated in simulation to a format ready for estimation
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