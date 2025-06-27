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

graph_stat_trends <- function(df_stats){
  
  # once we have a dataframe of all the graph stats, we are ready to plot
  #
  # 
  
}
