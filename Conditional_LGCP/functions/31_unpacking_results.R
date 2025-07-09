# once I have a set of saved files, how do I unpack and read them?

extract_pieces_cond <- function(files) {
  
  # 
  # 
  #
  # if our string is 
  #
  #           "WT1_m0vr0_w17.rda"
  #
  #
  # I want to keep
  # - 'WT1' to denote the mouse
  # - 0 to denote the movement
  # - 0 to denote the VR
  # - 17 to denote week
  #
  #
  # input:
  # 
  # - files (vector of strings)    .rda files
  
  
  rda_files <- sub("\\.rda$", "", files)
  # 1) WT1
  mouse_ID <- sub("_.*", "", rda_files)
  
  # 2) middle substring - m0vr0
  second_substring <- sub("^[^_]*_([^_]*)_.*", "\\1", rda_files)
  
  # 2a) movement, VR
  movement <- sub(".*m(.)?.*", "\\1", second_substring)
  VR <- sub(".*vr(.)?.*", "\\1", second_substring)
  
  # 3) third substring
  third_substring <- sub("^[^_]*_[^_]*_(.*)", "\\1", rda_files)
  week <- substr(third_substring, 2, nchar(third_substring))
  
  
  # Return all pieces as a named list
  df <- data.frame(
    mouse_ID = mouse_ID,
    movement = movement,
    VR = VR,
    ew_num = week
  )
  
  
  return(df)
}

extract_pieces_cond_v2 <- function(files) {
  
  # ----------------------------------------------------------------------------
  # 
  #
  # if our string is 
  #
  #           "WT1_m0vr0_age_0.5_ts_0.5.rda"
  #
  #
  # I want to keep
  # - 'WT1' to denote the mouse
  # - 0 to denote the movement
  # - 0 to denote the VR
  # - 0.5 to denote normalized age
  # - 0.5 to denote normalized timestamp
  #
  #
  # input:
  # 
  # - files (vector of strings)    .rda files
  #
  # output:
  #
  # -df (data.frame)
  #
  # ----------------------------------------------------------------------------
  
  
  rda_files <- sub("\\.rda$", "", files)
  # 1) WT1
  mouse_ID <- sub("_.*", "", rda_files)
  
  # 2) middle substring - m0vr0
  second_substring <- sub("^[^_]*_([^_]*)_.*", "\\1", rda_files)
  
  # 2a) movement, VR
  movement <- sub(".*m(.)?.*", "\\1", second_substring)
  VR <- sub(".*vr(.)?.*", "\\1", second_substring)
  
  # 3) third substring
  third_substring <- sub("^[^_]*_[^_]*_(.*)", "\\1", rda_files)
  third_split <- strsplit(third_substring, "_")
  
  age_norm <- sapply(third_split, `[`, 2) # in a list of vectors, extract each 2nd element
  ts_norm <- sapply(third_split, `[`, 4)  # in a list of vectors, extract each 4th element
  
  # Return all pieces as a named list
  df <- data.frame(
    mouse_ID = mouse_ID,
    movement = movement,
    VR = VR,
    age_norm = age_norm,
    ts_norm = ts_norm
  )
  
  
  return(df)
}

get_mvr_interaction <- function(summary_df){
  
  # assume summary_df has movement (0, 1) and VR (0, 1) as columns
  #
  # create m_vr_factor
  #
  # returns summary_df
  
  movement_factor <- rep('Resting', nrow(summary_df))
  movement_factor[summary_df$movement == 1] <- 'Running'
  movement_factor <- factor(movement_factor, levels = c('Resting', 'Running'))
  summary_df$movement_factor <- movement_factor
  
  # 2.2) same for VR: 0, 1,  are 'off', 'on',  respectively
  
  VR_factor <- rep('VR_Off', nrow(summary_df))
  VR_factor[summary_df$VR == 1] <- 'VR_On'
  VR_factor <- factor(VR_factor, levels = c('VR_Off', 'VR_On'))
  summary_df$VR_factor <- VR_factor
  
  # 2.3) and create a combined factor
  summary_df$m_vr_factor <- interaction(summary_df$movement_factor, summary_df$VR_factor)  
  
  return(summary_df)
}

graph_stat_trends <- function(df_stats){
  
  # once we have a dataframe of all the graph stats, we are ready to plot
  #
  # 
  
}
