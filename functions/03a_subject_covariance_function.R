# 8/28/2025
# - for the new estimation procedure, we need to get G_{i, j}(s,t)^k 
# - for each subject


subject_specific_log_intensity <- function(data_df4, Tseq_est){
  
  # ---------------------------------------------------------------------------- 
  #
  #
  # GOAL: calcualte subject-specific log intensities for all (n) subjects and (p) processes
  # 
  #
  # input:
  # 
  # - data_df4   ('time', 'feature_id', 'subject_num' dataframe)
  # - Tseq_est   (m-dim vector)   time discretizations
  # 
  # output:
  # 
  # - X_k_est  (p x m x n)  matrix of estimated log-intensities for each subject (n) and each process (p)
  #
  # 
  # ----------------------------------------------------------------------------
  
  m_est <- length(Tseq_est)
  n <- length(unique(data_df4$subject_num))
  p <- length(unique(data_df4$feature_id))
  
  # result_list_og <- data_df4 %>%
  #   group_by(feature_id, subject_num) %>%
  #   summarise(
  #     vec = list(estimate_log_intensity_function(time, Tseq_est)),  
  #     .groups = "drop"
  #   )
  
  
  # for missing groups, make their log-intensities N/A
  
  empty_vec <- rep(NA, m_est) 
  
  result_list <- data_df4 %>%
    group_by(feature_id, subject_num) %>%
    summarise(
      vec = list(estimate_log_intensity_function(time, Tseq_est)),  
      .groups = "drop"
    ) %>%
    # 2. Force all combinations to exist
    complete(
      feature_id = 1:p,          # Ensure all p features are present
      subject_num = 1:n,         # Ensure all n subjects are present
      fill = list(vec = list(empty_vec)) # Fill missing pairs with the empty vector
    ) %>%
    # 3. Sort to ensure the array alignment matches your dim() logic
    arrange(subject_num, feature_id) 
  
  # 4. Now unlist and array-ify safely
  X_k_est <- array(
    data = unlist(result_list$vec),
    dim = c(m_est, p, n) # Note: matching the arrange order (p changes fastest)
  ) %>% aperm(c(2, 1, 3)) # Reorder to [p, m_est, n]
  
  
  return(X_k_est)
  
}


