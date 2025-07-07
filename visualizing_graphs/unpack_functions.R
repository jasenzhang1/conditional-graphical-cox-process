# - After fitting LGCP models, I have several .rda results
# - I want to calculate and plot graph statistics and visualize them over time
#
# - good for unconditional and conditional models!

library(ggplot2)
library(igraph) # visualize graphs

adjacency_heatmap <- function(adj_mat, ID, region_border=NULL){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: plot an adjacency matrix graph 
  #
  # 
  # input:
  #
  # - adj_mat          (p x p matrix of 0's and 1's)
  # - ID               (string, name of mouse for the title)
  # - region_border    (integer)    neuron num cutoff between HIP and EHC
  #
  # output:
  #
  # - g3 (ggplot2)  adjacency matrix heatmap
  #
  # ----------------------------------------------------------------------------
  
  n <- dim(adj_mat)[1]
  m2_melt <- reshape2::melt(adj_mat)

  
  border_id3 <- 0:n + 0.5
  
  g3 <- ggplot(m2_melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Column", y = "Row", fill = "Value") +
    scale_y_reverse() +  # reverse y-axis scale
    geom_vline(xintercept = border_id3, alpha = 0.3) +  # neuron grid
    geom_hline(yintercept = border_id3, alpha = 0.3) +  # neuron grid
    
    # titles and theme
    xlab('Neuron') + 
    ylab('Neuron') + 
    ggtitle(ID) +
    theme_minimal() +
    theme(legend.position = "none")
    
  
  if(!is.null(region_border)){
    g3 <- g3 + 
      geom_vline(xintercept = region_border, color = 'dodgerblue', alpha = 0.6) + 
      geom_hline(yintercept = region_border, color = 'dodgerblue', alpha = 0.6)     
  }
  
  g3
}

insert_na_symmetric <- function(mat, indices) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: 
  #
  # - recall that in the unconditional model, some neurons were discarded and not included in the fitting procedure
  # - thus, our adjancency matrix is smaller
  # - we want to reconstruct the adjacency matrix with all neurons, since we removed some originally
  # - so we insert NA rows and columns in ascending order of neuron ID number
  #
  # 
  # input:
  # 
  # - mat     (matrix)                     truncated adjacency matrix
  # - indices (sorted array of integers)   sorted array of neuron numbers that were removed during estimation
  #
  # output:
  #
  # - mat    (p x p matrix)   recovered full adjacency matrix
  #
  # ----------------------------------------------------------------------------
  
  # print('starting conditions')
  # print(dim(mat)[1])
  # print(indices)
  # print(dim(mat)[1] + length(indices))
  
  if(length(indices) == 0){
    return(mat)
  }
  
  for(index in indices){
    
    d_i <- dim(mat)[1]
    
    
    if(index == 1){          # case 1: missing neuron is neuron 1
      mat <- rbind(NA, mat)
      mat <- cbind(NA, mat)
    } else if(index <= d_i){ # case 2: missing neuron is inserted inside the matrix
    
      mat_1 <- rbind( mat[1:(index-1), ],   # first rows
                      NA,                   # pad NA
                      mat[index:d_i, ])     # last rows
      
      mat_2 <- cbind(mat_1[, 1:(index-1)],  # first columns
                     NA,                    # pad NA
                     mat_1[,index:d_i])     # last columns
      
      mat <- mat_2
      
    } else if(index == d_i + 1){ # case 3: missing neuron is the very next one, so we append
      
      mat_1 <- rbind(mat, NA)    # pad a row
      mat_2 <- cbind(mat_1, NA)  # pad a column
      mat <- mat_2
      
    } else{ # case 4: the index is way beyond d_i
      stop('ERROR IN insert_na_symmetric')
    }
    
  }
  
  return(mat)
}

get_network <- function(adj_mat, ID){
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL:
  #
  # - with an unweighted adjacency matrix, we wish to visualize the graph
  # - visualize the graph with nodes and edges
  # - omit singletons
  #
  #
  # input:
  #
  # - adj_mat  (p x p matrix of 0's and 1's) adjacency matrix
  # - ID = 
  #
  # output:
  #
  # ----------------------------------------------------------------------------
  
  diag(adj_mat) <- 0
  verts <- which(apply(adj_mat,2,function(x){sum(x, na.rm = T)}) > 0)
  
  adj_mat2 <- adj_mat[verts, verts]
  colnames(adj_mat2) <- verts
  rownames(adj_mat2) <- verts
  
  
  g <- graph_from_adjacency_matrix(adj_mat2, mode = "undirected")  # or "directed"
  
  return(g)
  
  # plot(g, 
  #      vertex.label = V(g)$name,  # Optional: show vertex labels
  #      vertex.color = "lightblue",
  #      edge.color = "gray",
  #      edge.arrow.size = 0.5,
  #      main = ID)  
}

get_summary_statistics <- function(adj_mat_fill, weight_mat, brain_region_cutoff){
  
  # ---------------------------------------------------------------------------- 
  #
  # graphical summary statistics such as:
  # - average degree
  # - number of singletons
  #
  # input:
  #
  # - adj_mat_fill           (matrix):    the full nxn matrix of neurons with NA's denoting discarded neurons
  # - weight_mat             (matrix):    symmetric weight matrix that discards silent neurons
  #   - BOTH OF THESE MATRICES NEED 1'S ON THE DIAGONAL
  #   - diagonals of missing neurons are NA
  # 
  # 
  # - brain_region_cutoff (number): y + 0.5, which divides the y-th neuron and y+1-th neuron which are in different brain regions
  #
  # output:
  # 
  # -
  #
  # ---------------------------------------------------------------------------- 
  
  adj_mat_og <- adj_mat_fill[!apply(is.na(adj_mat_fill), 1, all),  # remove rows where all are NA
                             !apply(is.na(adj_mat_fill), 2, all)]  # still has 1's on the diagonal
  
  adj_mat2 <- adj_mat_og
  diag(adj_mat2) <- 0    # remove 1's on the diagonal
  
  # vertex statistics
  num_neurons <- dim(adj_mat_fill)[1]                           # number of neurons in total (169)
  num_NA <- sum(is.na(diag(adj_mat_fill)))                      # number of NA neurons       (42)
  num_candidates <- num_neurons - num_NA                        # number of non-NA neurons   (127)
  
  
  num_islands <- length(which(apply(adj_mat_og,2,sum) == 1))    # number of island neurons   (86)
  

  
  verts <- which(apply(adj_mat2,2,sum) > 0)
  num_con_verts <- length(verts)                                # number of connected neurons (41)
  
  
  # edge
  num_edges <- 0.5 * sum(adj_mat2)                              # number of edges            (61)
  
  off_diag_adj_mat <- adj_mat_fill[1:(brain_region_cutoff - 0.5), (brain_region_cutoff + 0.5):num_neurons]
  HH_adj_mat <- adj_mat_fill[1:(brain_region_cutoff - 0.5), 1:(brain_region_cutoff - 0.5)]
  EE_adj_mat <- adj_mat_fill[(brain_region_cutoff + 0.5):num_neurons, (brain_region_cutoff + 0.5):num_neurons]
  
  num_HE_edges <- sum(off_diag_adj_mat, na.rm = T)                                        # Hippocampus-EHC edges (17)
  num_HH_edges <- 0.5 * (sum(HH_adj_mat, na.rm = T) - sum(diag(HH_adj_mat), na.rm = T))   # Hip-Hip edges (28)
  num_EE_edges <- 0.5 * (sum(EE_adj_mat, na.rm = T) - sum(diag(EE_adj_mat), na.rm = T))   # EHC-EHC edges (16)
  
  
  # connectivity statistics
  avg_deg <- round(2 * num_edges/num_candidates, 2)                                       # average degree of candidate neurons (1.59) 
  avg_deg_normalized <- avg_deg/num_candidates

  # distribution of degrees
  degs <- table(apply(adj_mat2, 2, sum))

  
  # number of components
  adj_mat3 <- adj_mat2[verts, verts] # only keep the vertices that are connected
  g <- graph_from_adjacency_matrix(adj_mat3, mode = "undirected")
  num_comps <- components(g)$no                                                           # number of components (2)
  
  # subgraph statistics
  
  num_neurons_HIP <- dim(HH_adj_mat)[1]                           # number of neurons in total 
  num_NA_HIP <- sum(is.na(diag(HH_adj_mat)))                      # number of NA neurons       
  num_candidates_HIP <- num_neurons_HIP - num_NA_HIP              # number of non-NA neurons   
  avg_deg_HIP <- round(2 * num_HH_edges/num_candidates_HIP, 2)
  avg_deg_HIP_normalized <- avg_deg_HIP/num_candidates_HIP
    
  num_neurons_EHC <- dim(EE_adj_mat)[1]                           # number of neurons in total 
  num_NA_EHC <- sum(is.na(diag(EE_adj_mat)))                      # number of NA neurons       
  num_candidates_EHC <- num_neurons_EHC - num_NA_EHC              # number of non-NA neurons   
  avg_deg_EHC <- round(2 * num_EE_edges/num_candidates_EHC, 2)
  avg_deg_EHC_normalized <- avg_deg_EHC/num_candidates_EHC
  
  # weighted graph statistics
  
  g <- graph_from_adjacency_matrix(weight_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  avg_strength <- mean(strength(g, weights = E(g)$weight)) # average weight per node
  avg_edge_weight <- mean(E(g)$weight)                     # average weight per edge
  
  E(g)$weight <- 1 / E(g)$weight # now weight is inverted. The lower, the better.
  

  dist_mat <- distances(g, weights = E(g)$weight)  # path lengths between all pairs of nodes
  avg_dist <- mean(dist_mat[lower.tri(dist_mat)])  # average path length 
  
  graph_diameter <- diameter(g)        # diameter (worst case scenario for shortest path length)
  
  

  
  sum_stats <- list()
  detailed_stats <- list()
  total_stats <- list()
  
  # number of nodes
  sum_stats[['num_neurons']] <- num_neurons
  sum_stats[['num_NA']] <- num_NA
  sum_stats[['num_candidates']] <- num_candidates
  sum_stats[['num_islands']] <- num_islands
  sum_stats[['num_con_verts']] <- num_con_verts
  
  # number of edges
  
  sum_stats[['num_edges']] <- num_edges
  sum_stats[['num_HE_edges']] <- num_HE_edges
  sum_stats[['num_HH_edges']] <- num_HH_edges
  sum_stats[['num_EE_edges']] <- num_EE_edges
  
  # misc stats
  sum_stats[['avg_deg']] <- avg_deg
  sum_stats[['avg_deg_HIP']] <- avg_deg_HIP
  sum_stats[['avg_deg_EHC']] <- avg_deg_EHC
  sum_stats[['avg_deg_normalized']] <- avg_deg_normalized
  sum_stats[['avg_deg_HIP_normalized']] <- avg_deg_HIP_normalized
  sum_stats[['avg_deg_EHC_normalized']] <- avg_deg_EHC_normalized
  
  sum_stats[['num_comps']] <- num_comps
  
  # weighted graph stats
  sum_stats[['avg_node_strength']] <- avg_strength
  sum_stats[['avg_edge_strength']] <- avg_edge_weight
  sum_stats[['avg_dist']] <- avg_dist
  sum_stats[['diameter']] <- graph_diameter
  
  # detailed stats
  detailed_stats[['degs']] <- degs
  detailed_stats[['verts']] <- verts

  
  # merge
  total_stats[['sum_stats']] <- sum_stats
  total_stats[['detailed_stats']] <- detailed_stats
  
  
  return(total_stats)
  
}

save_summary_statistics <- function(figures_root, figure_name, settings_kept, sum_stats){
  
  txt_name <- paste(figure_name, '.txt', sep = '_summary')
  txt_name <- paste(figures_root, txt_name, sep = '')
  
  sink(txt_name)
  
  for(i in 1:length(sum_stats)){
    
    cat('-----------------------------------------\n')
    
    cat('Settings:\n')
    
    cat('Mouse: ')
    cat(settings_kept[[i]][['mouse']])
    cat('\n')
    cat('Week/epoch num: ')
    cat(settings_kept[[i]][['time']])
    cat('\n')
    cat('Movement: ')
    cat(settings_kept[[i]][['movement']])
    cat('\n') 
    cat('VR: ')
    cat(settings_kept[[i]][['VR']])
    cat('\n')
    
    cat('-----------------------------------------\n')
    
    cat('Total Neurons:\t')
    cat(sum_stats[[i]][['num_neurons']])
    cat('\n')
    cat('Discarded Neurons:\t')
    cat(sum_stats[[i]][['num_NA']])
    cat('\n')
    cat('Island Neurons:\t')
    cat(sum_stats[[i]][['num_islands']])
    cat('\n')
    cat('Connected Neurons:\t')
    cat(sum_stats[[i]][['num_con_verts']])
    cat('\n')   
    cat('Number of Edges:\t')
    cat(sum_stats[[i]][['num_edges']])
    cat('\n')  
    cat('Average Degree:\t')
    cat(sum_stats[[i]][['avg_deg']])
    cat('\n')      
    cat('Number of Components:\t')
    cat(sum_stats[[i]][['num_comps']])
    cat('\n\n')  
    cat('Degree Distribution: \n')
    cat('Degree:\t')
    cat(names(sum_stats[[i]][['degs']])) 
    cat('\n')
    cat('Count:\t')
    cat(unname(sum_stats[[i]][['degs']]))  
    
    
    cat('\n\n--------------\n\n')  
  }
  
  sink()  
}

get_graph_laplacian_stats <- function(adj_mat, zero_tol = 1e-8){
  
  # 
  # GOAL: get graph laplacian statistics
  #
  # 
  # input:
  # 
  # - adj_mat (p x p adjacency matrix): symmetric 
  # - zero_tol (number): if eigenvalues are below this, set them equal to 0
  #
  #
  # output:
  #
  # - list of graph laplacian statistics
  
  # if we feed in the adj_mat with NA's, set them to 0 except the diagonal (to make them singletons)
  adj_mat[is.na(adj_mat)] <- 0
  diag(adj_mat) <- 1
  
  p <- dim(adj_mat)[1]
  lap_mat <- diag(rowSums(adj_mat)) - adj_mat
  
  
  # eigendecomposition
  
  lap_eigen <- eigen(lap_mat, symmetric = T)
  lap_eval <- lap_eigen$values
  lap_eval[abs(lap_eval) < zero_tol] <- 0
  
  fiedler_value <- lap_eval[length(lap_eval) - 1] # lambda_2
  lambda_max <- max(lap_eval)     # largest eigenvalue
  num_zeros <- sum(lap_eval == 0) # number of zeros
  lambda_median <- median(lap_eval)
  lambda_lower <- quantile(lap_eval, 0.25) %>% unname()
  lambda_upper <- quantile(lap_eval, 0.75) %>% unname()
  
  # weighted laplacian
  
  deg_inv_sqrt <- diag(1 / sqrt(rowSums(adj_mat)))
  lap_mat_sym <- diag(rep(1, p)) - deg_inv_sqrt %*% adj_mat %*% deg_inv_sqrt
  
  lap_sym_eval <- eigen(lap_mat_sym, symmetric = T)$values
  lap_sym_eval[abs(lap_sym_eval) < zero_tol] <- 0
  
  lambda_max_sym <- max(lap_sym_eval)  # largest eigenvalue of symmetric laplacian
  lambda_median_sym <- median(lap_sym_eval)  
  lambda_lower_sym <- quantile(lap_sym_eval, 0.25) %>% unname() 
  lambda_upper_sym <- quantile(lap_sym_eval, 0.75) %>% unname()
  # compile and print results
  
  results <- list(num_zeros=num_zeros,            # number of zero eigenvalues
                  fiedler_value=fiedler_value,    # fiedler value
                  lambda_max=lambda_max,          # lambda max
                  lambda_median=lambda_median,
                  lambda_25=lambda_lower,
                  lambda_75=lambda_upper,
                  lambda_max_sym=lambda_max_sym,  # lambda max of symmetric matrix
                  lambda_median_sym=lambda_median_sym,
                  lambda_25_sym=lambda_lower_sym,
                  lambda_75_sym=lambda_upper_sym)
  
  return(unlist(results))
  
}

get_node_stats_over_time <- function(weighted_edge_lists, week_nums, node_num, threshold = 0){
  
  #
  # GOAL: 
  #
  # - We have a longitudinal collection of lists of weighted edges
  # - For a single neuron, how does its stats change over time?
  #
  #
  # input:
  #
  # - weighted_edge_lists (list of lists)
  #   - each inner list contains entries of the form [['i_j']], which stores scalars (weight value)
  # - week_nums           (vector of week numbers)
  # - node_num            (number)
  # - threshold           (number) anything below threshold will be set to 0
  #
  # output:
  #
  # - df_total (data.frame) 'weight', 'week', 'connected_node'
  # 
  
  list_names <- names(weighted_edge_lists)
  n_lists <- length(weighted_edge_lists)
  
  df_total <- data.frame()
  
  for(i in 1:n_lists){
    list_name_i <- list_names[i]
    weighted_edge_list <- weighted_edge_lists[[list_name_i]]
    
    # unpack
    
    v_i <- sub("_.*", "", names(weighted_edge_list)) %>% as.numeric()
    v_j <- sub(".*_", "", names(weighted_edge_list)) %>% as.numeric()
    
    df_adj <- data.frame(v_i, v_j, unlist(weighted_edge_list))   
    colnames(df_adj) <- c('node_i', 'node_j', 'weight')
    
    df_node <- df_adj %>% dplyr::filter(node_i == node_num | node_j == node_num) %>% 
      mutate(weight = weight * (weight > threshold)) %>% 
      mutate(week = week_nums[i]) %>% 
      mutate(connected_node = ifelse(node_i == node_num, node_j, node_i)) %>% 
      select(-c('node_i', 'node_j'))
    
    # store
    df_total <- rbind(df_total, df_node)
  }
  
  return(df_total)
  
}

extract_pieces <- function(x) {
  
  # if our string is "week_move_2/w33_m2vr0_n50_me1"
  # I want to keep
  # - 'w' to denote that we are working in weeks
  # - 33 to denote the week
  # - 2 to denote the movement
  # - 0 to denote the VR
  # - 50 to denote the replicate count
  # - 1 to denote the minimum edges
  

  epoch_or_week <- sub(".*/(.)?.*", "\\1", x) #w
  
  # 1️⃣ Between last '/' and next '_', then delete first character
  first_piece <- sub("^.*/([^_]+)_.*$", "\\1", x) %>% substring(2) %>% as.numeric() #33
  
  # 2️⃣ After 'm' after the 3rd-to-last underscore
  # Find positions of underscores
  us_pos <- gregexpr("_", x)[[1]]
  third_last_us <- us_pos[length(us_pos) - 2]
  char_pos <- third_last_us + 2
  second_piece <- substring(x, char_pos, char_pos) %>% as.numeric() #2
  
  # 3️⃣ After 'vr' before second-to-last underscore
  second_last_us <- us_pos[length(us_pos) - 1]
  char_pos <- second_last_us - 1
  third_piece <- substring(x, char_pos, char_pos) %>% as.numeric() #0

  
  char_pos <- second_last_us + 1
  rep_or_time_scale = substring(x, char_pos, char_pos) #n
  
  # 4️⃣ Before last underscore
  last_us <- us_pos[length(us_pos)]
  start_pos <- second_last_us + 2
  end_pos <- last_us - 1
  fourth_piece <- substring(x, start_pos, end_pos) %>% as.numeric() #50
  
  # 5️⃣ Last character
  fifth_piece <- substring(x, nchar(x)) %>% as.numeric() #1
  
  # Return all pieces as a named list
  list(
    epoch_or_week = epoch_or_week,
    ew_num = first_piece,
    movement = second_piece,
    VR = third_piece,
    rep_or_ts = rep_or_time_scale,
    rep_ts_num = fourth_piece,
    min_edge = fifth_piece
  )
}

extract_pieces2 <- function(x) {
  
  # character between last / and next _
  # /w17_m0vr0_n50_me1/Tau1_w17_m0vr0_n50_me1.rda
  # return 'Tau1'
  
  result <- sapply(x, function(s) {
    last_slash_pos <- max(gregexpr("/", s)[[1]])
    next_us_pos <- regexpr("_", substring(s, last_slash_pos + 1))
    if (next_us_pos == -1) return(NA_character_)
    substring(s, last_slash_pos + 1, last_slash_pos + next_us_pos - 1)
  }) 
  
  return(unname(result))
}

unpack <- function(n, movement, VR, epoch_or_week, ew_num, min_edges, df_brain_region, first_root){
  
  # 
  # inputs:
  # 
  # - n (number): replicate count
  # - movement(number): 2 = running, 1 = resting, 0 = both
  # - VR (number):      2 = on, 1 = off, 0 = do not filter
  # - epoch_or_week (string): 'week' or 'epoch' to describe which unit of time we're working with
  # - ew_num (integer): week or epoch number 
  # - min_edges (integer): mininum number of edges in our network, usually 0 or 1
  # - df_brain_region (dataframe)
  # - first_root (string): path to the set of folders
  #
  #
  # 
  
  # 1) obtain file directory
  if(epoch_or_week == 'week'){
    setting_ID <- paste('w', '_m', sep = as.character(ew_num))          # w1_m
  } else{
    setting_ID <- paste('e', '_m', sep = as.character(ew_num))          # e1_m
  }

  setting_ID <- paste(setting_ID, 'vr', sep = as.character(movement))    # e1_m2vr
  setting_ID <- paste(setting_ID, '_n', sep = as.character(VR))          # e1_m2vr0_n
  setting_ID <- paste(setting_ID, as.character(n), sep = '')             # e1_m2vr0_n400
  setting_ID <- paste(setting_ID, as.character(min_edges), sep = '_me')  # e1_m2vr0_n400_me1
  
  results_root <- paste(first_root, #"/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/"
                        '/',
                        sep = setting_ID)
  
  figures_root <- paste(first_root, #"/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/"
                        '/',
                        sep = 'figures')  
  
  # if the folder doesn't exist, create it 
  if (!dir.exists(figures_root)) {
    dir.create(figures_root)
  }     
  
  # 2) prepare to store results
  graphs <- list()
  g_vis <- list()
  edge_vec <- c()
  neuron_vec <- c()
  sum_stats <- list()
  
  file_names <- list.files(path = results_root, full.names = TRUE)
  
  print(file_names)
  
  m <- length(file_names)
  no_ext <- sub("\\.[^.]*$", "", file_names)
  rda_names  <- sub(".*/", "", no_ext)                     # rda names
  mice_names <-  sub(".*/([^/_]+)_.*", "\\1", file_names)  # get mice names, between last '/' and very next '_'
  
  for(i in 1:m){
    
    
    # 2) loading
    print(file_names[i])
    load(file_names[i])
    
    # 3) fill in NA's
    adj_mat_og <- graph_all[['GPP_BIC']]
    missing_neurons <- sort(graph_all[['missing_neurons']])
    adj_mat_fill <- insert_na_symmetric(adj_mat_og, missing_neurons)
    
    # 4) graph and store
    region_border <- max(df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == mice_names[i]]) + 0.5
    
    g_i <- adjacency_heatmap(adj_mat_fill, mice_names[i], region_border=region_border)
    
    graphs[[i]] <- g_i
    
    # 5) visualization of network
    g_vis[[i]] <- get_network(adj_mat_fill, rda_names[i])
    
    # 5) summary statistics
    
    sum_stats[[i]] <- get_summary_statistics(adj_mat_fill, adj_mat_og)
    
  }
  
  # save heatmap ---------------------------------------------------
  
  pdf_name <- paste(setting_ID, '.pdf', sep = '_heatmap')
  pdf_name <- paste(figures_root, pdf_name, sep = '')
  pdf(pdf_name, width = 10, height = 10)
  
  for(i in 1:m){
    print(graphs[[i]])
  }
  
  dev.off()
  
  # save graph visualization --------------------------------------
  
  pdf_name <- paste(setting_ID, '.pdf', sep = '_graph_vis')
  pdf_name <- paste(figures_root, pdf_name, sep = '')
  pdf(pdf_name, width = 10, height = 10)
  
  for(i in 1:m){
    
    HIP_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == mice_names[i]] 
    EHC_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Entorhinal_Cortex" & df_brain_region$ID2 == mice_names[i]] 
    
    plot_pts <- sum_stats[[i]][['verts']]
    
    plot_pts_color <- rep('gold', length(plot_pts))
    plot_pts_color[plot_pts %in% EHC_ID] <- 'dodgerblue'
    
    plot(g_vis[[i]], 
         vertex.size = 8,
         vertex.label = V(g_vis[[i]])$name,
         vertex.color = plot_pts_color,
         edge.color = "black",
         edge.arrow.size = 0.5,
         main = mice_names[i])
  }
  
  dev.off()
  
  # save summary -------------------------------------------------
  
  txt_name <- paste(setting_ID, '.txt', sep = '_summary')
  txt_name <- paste(figures_root, txt_name, sep = '')
  
  sink(txt_name)
  
  for(i in 1:m){
    
    cat('-----------------------------------------\n')
    
    cat('Mouse ID: ')
    cat(mice_names[i])
    cat('\n\n')
    
    cat('Total Neurons:\t')
    cat(sum_stats[[i]][['num_neurons']])
    cat('\n')
    cat('Discarded Neurons:\t')
    cat(sum_stats[[i]][['num_NA']])
    cat('\n')
    cat('Island Neurons:\t')
    cat(sum_stats[[i]][['num_islands']])
    cat('\n')
    cat('Connected Neurons:\t')
    cat(sum_stats[[i]][['num_con_verts']])
    cat('\n')   
    cat('Number of Edges:\t')
    cat(sum_stats[[i]][['num_edges']])
    cat('\n')      
    cat('Number of Components:\t')
    cat(sum_stats[[i]][['num_comps']])
    cat('\n\n')  
    cat('Degree Distribution: \n')
    cat('Degree:\t')
    cat(names(sum_stats[[i]][['degs']])) 
    cat('\n')
    cat('Count:\t')
    cat(unname(sum_stats[[i]][['degs']]))  
    
    
    cat('\n\n--------------\n\n')  
  }
  
  sink()
  
  
  return(NULL)
}

unpack_tabular_summary <- function(n, movement, VR, epoch_or_week, ew_num, rep_or_ts, rep_ts_num, mouse_ID, min_edges, df_brain_region, file_name){
  
  # 
  # goal: 
  #
  # for each fitted model, we want to extract key settings, fitting parameters, and graph statistics
  # 
  # 
  # 
  # inputs:
  # 
  # - n (number): replicate count
  # - movement(number): 2 = running, 1 = resting, 0 = both
  # - VR (number):      2 = on, 1 = off, 0 = do not filter
  # - epoch_or_week   (string): 'week' or 'epoch' to describe which unit of time we're working with
  # - ew_num          (integer): week or epoch number 
  # - rep_or_ts        (string): 't' or 'n' denoting whether we are partitioning by timescale or number of replicates
  # - rep_ts_num       (integer): length of time scale, or number of replicates
  #
  # - mouse_ID (string): 'Tau1'
  # - min_edges (integer): mininum number of edges in our network, usually 0 or 1
  # - df_brain_region (dataframe)
  # - file_name (string): the file of interest "/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/week_move_2/w17_m0vr0_n50_me1/Tau1_w17_m0vr0_n50_me1.rda"
  #
  #
  # outputs:
  # 
  # a single vector with the following attributes:
  # 
  # 1) mouse_ID
  # 2) ew_num
  # 3) movement
  # 4) VR
  # 5) min_edges
  # 6) num_neurons
  # 7) num_NA
  # 8) num_candidates
  # 9) num_islands
  # 10) num_con_verts
  # 11) num_edges
  # 12) num_HE_edges
  # 13) num_HH_edges
  # 14) num_EE_edges
  # 15) avg_deg
  # 16) num_comps
  # 
  
  
  # 1) prepare to store results
  result_df <- data.frame()
  
  # 2) loading
  load(file_name)
  
  # 3) fill in NA's
  adj_mat_og <- graph_all[['GPP_BIC']]
  missing_neurons <- sort(graph_all[['missing_neurons']])
  adj_mat_fill <- insert_na_symmetric(adj_mat_og, missing_neurons)
  
  adj_mat_w_og <- graph_all[['weighted_GPP_BIC']]
  adj_mat_w_fill <- insert_na_symmetric(adj_mat_w_og, missing_neurons)
  
  # 3.1) summary stats
  brain_region_cutoff <- max(which(df_brain_region$Brain_Region[df_brain_region$ID2 == mouse_ID] == 'Hippocampus')) + 0.5
  
  all_stats <- get_summary_statistics(adj_mat_fill, adj_mat_og, brain_region_cutoff)

  
  # 4) package results
  results <- c(mouse_ID,
               ew_num,
               movement,
               VR,
               graph_all[['time_scale']],
               graph_all[['num_replicates']],
               min_edges,
               graph_all[['tuning_parameters']],
               graph_all[['tuning_parameter_indices']],
               graph_all[['tuning_parameter_max_indices']],
               unlist(all_stats[['sum_stats']])

               )
  
  names(results) <- c('mouse_ID', 'ew_num', 'movement', 'VR', 'time_scale', 'num_replicates', 'min_edges',   # settings for the simulation
                      'tau_c', 'tau_p',                                                                      # parameters of the result of the simulation
                      'tau_c_i', 'tau_p_i',
                      'tau_c_max', 'tau_p_max',
                      names(unlist(all_stats[['sum_stats']])))                                                # graph statistics from the simulation
    

  results_df <- t(as.data.frame(results))
  rownames(results_df) <- NULL
  
  return(results_df)
}

unpack_over_time <- function(n, movement, VR, epoch_or_week, ew_nums, mouse_ID, min_edges, df_brain_region, first_root){
  
  # 1) prepare to store results
  graphs <- list()
  g_vis <- list()
  edge_vec <- c()
  neuron_vec <- c()
  sum_stats <- list()  
  
  graph_index <- 0
  ew_kept <- c()      # times kept
  settings_kept <- list()
  
  for(ew_num in ew_nums){ # 2) FOR EACH TIMEPOINT

    # 2.1) obtain file directory
    if(epoch_or_week == 'week'){
      setting_ID <- paste('w', '_m', sep = as.character(ew_num))          # w1_m
    } else{
      setting_ID <- paste('e', '_m', sep = as.character(ew_num))          # e1_m
    }
    
    setting_ID <- paste(setting_ID, 'vr', sep = as.character(movement))    # e1_m2vr
    setting_ID <- paste(setting_ID, '_n', sep = as.character(VR))          # e1_m2vr0_n
    setting_ID <- paste(setting_ID, as.character(n), sep = '')             # e1_m2vr0_n400
    setting_ID <- paste(setting_ID, as.character(min_edges), sep = '_me')  # e1_m2vr0_n400_me1
    
    results_root <- paste(first_root, #"/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/"
                          '/',
                          sep = setting_ID)
    
    figures_root <- paste(first_root, #"/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/"
                          '/',
                          sep = 'figures')  
    
    # if the folder doesn't exist, create it 
    if (!dir.exists(figures_root)) {
      dir.create(figures_root)
    }     
    
    # 2.2) Obtain file names
    
    file_names <- list.files(path = results_root, full.names = TRUE)
    
    m <- length(file_names)
    no_ext <- sub("\\.[^.]*$", "", file_names)
    rda_names  <- sub(".*/", "", no_ext)                     # rda names
    mice_names <-  sub(".*/([^/_]+)_.*", "\\1", file_names)  # get mice names, between last '/' and very next '_'  
    
    if(mouse_ID %in% mice_names){ # 3) if this mouse is in this time, examine it
      
      ew_kept <- c(ew_kept, ew_num)
      mouse_index <- which(mice_names == mouse_ID)
      
      # 2) loading
      load(file_names[mouse_index])
      
      # 3) fill in NA's
      adj_mat_og <- graph_all[['GPP_BIC']]
      missing_neurons <- sort(graph_all[['missing_neurons']])
      adj_mat_fill <- insert_na_symmetric(adj_mat_og, missing_neurons)
      
      # 4) graph and store
      graph_index <- graph_index + 1
      region_border <- max(df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == mouse_ID]) + 0.5
      
      g_i <- adjacency_heatmap(adj_mat_fill, mouse_ID, region_border=region_border)
      
      graphs[[graph_index]] <- g_i
      
      # 5) visualization of network
      g_vis[[graph_index]] <- get_network(adj_mat_fill, rda_names[mouse_index])
      
      # 5) summary statistics
      sum_stats[[graph_index]] <- get_summary_statistics(adj_mat_fill, adj_mat_og)
      
      # 6) settings
      settings_kept[[graph_index]] <- list(mouse = mouse_ID,
                                           time = ew_num,
                                           movement = movement,
                                           VR = VR)
      
    }
        
  }

  
  # save heatmap ---------------------------------------------------
  
  figure_name <- paste(mouse_ID, movement, sep = '_m') #Tau1_m2
  figure_name <- paste(figure_name, VR, sep = 'vr') #Tau1_m2vr0
  
  pdf_name <- paste(figure_name, '.pdf', sep = '_heatmap')
  pdf_name <- paste(figures_root, pdf_name, sep = '')
  pdf(pdf_name, width = 10, height = 10)
  
  for(i in 1:graph_index){
    print(graphs[[i]])
  }
  
  dev.off()
  
  # save graph visualization --------------------------------------
  

  
  pdf_name <- paste(figure_name, '.pdf', sep = '_graph_vis')
  pdf_name <- paste(figures_root, pdf_name, sep = '')
  pdf(pdf_name, width = 20, height = 10)
  
  n <- length(g_vis) # start plotting
  plots_per_page <- 3
  
  for (i in seq(1, n, by=plots_per_page)) {
    # Start a new page
    par(mfrow = c(1, 3))
    
    # Plot up to 3 plots in this page
    for (j in i:min(i + plots_per_page - 1, n)) {
      
      HIP_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Hippocampus" & df_brain_region$ID2 == mice_names[j]] 
      EHC_ID <- df_brain_region$Neuron_Num[df_brain_region$Brain_Region == "Entorhinal_Cortex" & df_brain_region$ID2 == mice_names[j]] 
      
      plot_pts <- sum_stats[[j]][['verts']]
      
      plot_pts_color <- rep('gold', length(plot_pts))
      plot_pts_color[plot_pts %in% EHC_ID] <- 'dodgerblue'      
      
      
      plot(g_vis[[j]], 
           vertex.size = 8,
           vertex.label = V(g_vis[[j]])$name,
           vertex.color = plot_pts_color,
           edge.color = "black",
           edge.arrow.size = 0.5,
           main = ew_kept[j])
    }
    
    # Fill empty slots if the last page has fewer than 3 plots
    for (k in (j+1):(i + plots_per_page - 1)) {
      plot.new()  # blank plot
    }
  }
  
  dev.off()
  
  # save summary -------------------------------------------------
  
  save_summary_statistics(figures_root, figure_name, settings_kept, sum_stats)  
  
  
  return(NULL)
}

type_1_graphs <- function(summary_df, value, color_var, group_var, mice_strains2){
  
  #
  # 
  
  g <- ggplot(summary_df, aes(x = ew_num, y = .data[[value]], color = .data[[color_var]])) + 
    geom_point() + 
    geom_line() + 
    facet_wrap(vars(.data[[group_var]])) + 
    ylab(value) + 
    xlab('Week') + 
    theme_bw() + 
    scale_color_manual(values = mice_strains2)
  
  return(g)
}

plot_graph_stats <- function(summary_df, final_results_dir, thresh_value){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: plot graph statistics over time
  #
  # - average degree
  # - percent of HH edges
  # - percent of EE edges
  # 
  # 
  # Input:
  # 
  # - thresh_value   (number)
  # - summary_df     (data.frame)
  #
  #   - mouse_ID (string)
  #   - movement (0 = rest, 1 = run, 2 = both)
  #   - VR       (0 = off, 1 = on, 2 = both)
  #   - ew_num
  #   - num_neurons
  #   - num_NA
  #   - num_candidates
  #   - num_islands
  #   - num_con_verts
  #   - num_edges
  #   - num_HE_edges
  #   - num_EE_edges
  #   - avg_deg
  #   - avg_deg_HIP
  #   - avg_deg_EHC
  #   - num_comps
  #
  # Output:
  # 
  # - saved pdf of all the graphs
  # - save summary_df as well
  #
  # ----------------------------------------------------------------------------
  
  
  
  # 1) make everything but mouse_ID numeric, make mouse_ID factor
  summary_df[-which(names(summary_df) == "mouse_ID")] <- lapply(summary_df[-which(names(summary_df) == "mouse_ID")], as.numeric)
  summary_df$mouse_ID <- factor(summary_df$mouse_ID)
  
  
  # 2) factors
  
  # change movement from 0, 1, 2 to 'resting', 'running', 'both' respectively
  
  movement_factor <- rep('Resting', nrow(summary_df))
  # movement_factor[summary_df$movement == 0] <- 'Resting'
  movement_factor[summary_df$movement == 1] <- 'Running'
  movement_factor <- factor(movement_factor, levels = c('Resting', 'Running'))
  summary_df$movement_factor <- movement_factor
  
  # same for VR: 0, 1, 2 are 'off', 'on', 'both' respectively
  
  VR_factor <- rep('VR_Off', nrow(summary_df))
  # VR_factor[summary_df$VR == 0] <- 'Off'
  VR_factor[summary_df$VR == 1] <- 'VR_On'
  VR_factor <- factor(VR_factor, levels = c('VR_Off', 'VR_On'))
  summary_df$VR_factor <- VR_factor
  
  # and create a combined factor
  summary_df$m_vr_factor <- interaction(summary_df$movement_factor, summary_df$VR_factor)
  
  # more statistics
  summary_df$pct_non_singletons <- 100 * summary_df$num_con_verts/summary_df$num_candidates
  
  
  summary_df$pct_inter_region_edge <- 100 * summary_df$num_HE_edges/summary_df$num_edges
  summary_df$pct_HH_edge <- 100 * summary_df$num_HH_edges/summary_df$num_edges
  summary_df$pct_EE_edge <- 100 * summary_df$num_EE_edges/summary_df$num_edges  
  
  summary_df$pct_candidate_neurons <- 100 * summary_df$num_candidates/summary_df$num_neurons
  
  ## 3) graphing the summary statistics over time
  
  mice_strains2 <- c(
    "Tau1" = "lightcoral",
    "Tau2" = "red",
    "Tau3" = "darkred",
    "WT1" = "lightgreen",
    "WT2" = "green",
    "WT3" = "darkgreen"
  )
  
  x_graph_title <- ('Week')
  graphs <- list()
  graphs_group_by_mouse <- list()
  
  # 3.1) average degree over time
  # - ew_num
  # - avg_deg
  # - mouse_ID
  
  value_names <- c('avg_deg', 'normalized_avg_deg', 'avg_deg_HIP', 'avg_deg_HIP_normalized',
                   'avg_deg_EHC', 'avg_deg_EHC_normalized',
                   'pct_non_singleton', 'pct_HE', 'pct_HH', 'pct_EE', 'pct_candidates',
                   'num_components',
                   'lap_eval', 'median_lap_eval', 'max_lap_eval',
                   'lap_eval_sym', 'median_lap_eval_sym', 'max_lap_eval_sym',
                   'avg_node_strength', 'avg_edge_strength', 'avg_dist', 'diameter')
  color_var <- 'mouse_ID'
  group_var <- 'm_vr_factor'
  
  for(k in 1:length(value_names)){
    graphs[[value_names[k]]] <-                type_1_graphs(summary_df, value_names[k], color_var, group_var, mice_strains2)
    graphs_group_by_mouse[[value_names[k]]] <- type_1_graphs(summary_df, value_names[k], group_var, color_var, mice_strains2)
  
  }
  

  
  # 6) number of replicates  - not relevant for conditional model
  
  # graphs[['num_replicates']] <- ggplot() + 
  #   geom_point(data = summary_df, aes(x = ew_num, y = num_replicates, color = mouse_ID)) + 
  #   geom_line(data = summary_df, aes(x = ew_num, y = num_replicates, color = mouse_ID)) + 
  #   facet_wrap(~ m_vr_factor) + 
  #   ylab('Number of Replicates') + 
  #   xlab(x_graph_title) +  
  #   theme_bw() + 
  #   scale_color_manual(values = mice_strains2)
  
  # save data
  
  
  
  if (!dir.exists(final_results_dir)) {
    dir.create(final_results_dir)
  }
  
  # store summary statistics of each fit
  
  write.csv(summary_df, file = paste0(final_results_dir, 'graph_statistics_thresh_', thresh_value, '.csv'), row.names = FALSE)
  
  # graphs
  
  pdf(paste0(final_results_dir, 'graph_statistics_group_by_strata_thresh_', thresh_value, '.pdf'), width = 8, height = 3)
  for(name_i in names(graphs)){
    print(graphs[[name_i]])
  }
  dev.off()  
  
  pdf(paste0(final_results_dir, 'graph_statistics_group_by_mouse_thresh_', thresh_value, '.pdf'), width = 8, height = 3)
  for(name_i in names(graphs)){
    print(graphs[[name_i]])
  }
  dev.off()    
  
  # no return
}
