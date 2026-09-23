make_edge_weight_graph <- function(silence_df, ID, ID_num, neuron_i, neuron_j, movements, VRs, weeks, results_dir){
  
  # silence_df () object that is saved
  # ID         (string) 'Tau1'
  # ID_num     (string) '346'
  # neuron_i   (integer) 
  # neuron_j   (integer)
  # movements  (vector)
  # VRs        (vector)
  # weeks      (vector)

  # start ------------------------------------------------------------------------
  
  # obtain silence status
  
  silence_i <- silence_df[[ID_num]][, neuron_i]
  silence_j <- silence_df[[ID_num]][, neuron_j]  
  
  silence_graph_df <- data.frame(silence_i, silence_j, rownames(silence_df[[ID_num]]))
  colnames(silence_graph_df) <- c('silence_i', 'silence_j', 'week')
  silence_graph_df$week <- as.numeric(silence_graph_df$week)
  
  weight_df <- data.frame()
  
  for(movement in movements){
    for(VR in VRs){
      
      setting_name <- paste0('m', movement, 'vr', VR)
      
      
      for(week in weeks){
        
        # 1.1) craft file name and load it
        
        f <- paste0(ID, '_m', movement, 'vr', VR, '_w', week, '.rda')
        f2 <-  sub("\\.rda$", "", f)
        load(paste0(results_dir, '/', f)) # called final_graph_estimates
        
        
        index_name <- paste0(min(neuron_i, neuron_j), '_', max(neuron_i, neuron_j))
        weight_ij <- final_graph_estimates$edge_strengths[[index_name]]
        
        
        weight_df <- rbind(weight_df, data.frame(weight_ij, week, setting_name))
        
      }
    }
  }
  
  colnames(weight_df) <- c('weight', 'week', 'setting')
  weight_df$setting <- factor(weight_df$setting)
  
  # plot
  
  df_shade_i <- silence_graph_df %>%
    filter(silence_i == 1) %>%
    mutate(xmin = week - 0.5,
           xmax = week + 0.5,
           ymin = -Inf,
           ymax = Inf)
  
  df_shade_j <- silence_graph_df %>%
    filter(silence_j == 1) %>%
    mutate(xmin = week - 0.5,
           xmax = week + 0.5,
           ymin = -Inf,
           ymax = Inf)
  
  figure_name <- paste0(ID, '_neurons_', neuron_i, '_', neuron_j)
  
  g_edge_silence <- ggplot() + geom_point(data = weight_df, aes(x = week, y = log(weight), group = setting, color = setting)) +
    geom_line(data = weight_df, aes(x = week, y = log(weight), group = setting, color = setting)) + 
    theme_bw() + 
    ylab('Log(Edge Weight)') + 
    xlab('Week') + 
    ggtitle(figure_name) + 
    geom_rect(data = df_shade_i,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "gray", alpha = 0.3, inherit.aes = FALSE) +
    geom_rect(data = df_shade_j,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "red", alpha = 0.3, inherit.aes = FALSE)  
  
  return(g_edge_silence)
}
