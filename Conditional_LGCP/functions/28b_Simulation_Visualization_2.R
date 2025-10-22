source('functions/28z_Visualization_helpers.R')
source('functions/28y_Visualization_Blocks.R')
source('functions/28_Simulation_Visualization.R')

library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)

visualize_points <- function(step_0_events){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: plot the events from t = 0 to t = 1 for various processes
  #
  # 
  # input: 
  #
  # - step_0_events   (list)  each item denotes a process and is a vector of points
  #
  # 
  # output:
  #
  # - g     (graph)
  #
  # ----------------------------------------------------------------------------
  
  
  
  df <- do.call(rbind, lapply(seq_along(step_0_events), function(i) {
    data.frame(process = i, time = events[[i]])
  }))  
  
  g <- ggplot(df, aes(x = time, y = process)) +
    geom_point(size = 0.8, alpha = 0.5) +
    labs(x = "Time", y = "Process", title = "Point Process Events of First Subject") +
    theme_bw() + 
    theme(panel.grid = element_blank())

  
}

visualize_points_on_intensity <- function(step_0_events, step_1, time_grid){
  
  
  # ----------------------------------------------------------------------------
  #
  # 
  # GOAL: map the events onto the estimated and true intensities from step_1
  #
  # inputs:
  #
  # - step_0_events   (list)    each item denotes a process and is a vector of points
  # - step_1          (list)    
  #   - X_k_est                   (p x m_est   x n)
  #   - X_k_truth                 (p x m_truth x n)
  #   - X_k_coarse_truth          (p x m_est   x n)
  #   - X_k_both_truth            (p x m_both  x n)
  #   - Lambda_k_truth            (p x m_truth x n)
  #   - Lambda_k_coarse_truth     (p x m_est   x n)
  #   - Lambda_k_est              (p x m_est   x n)
  #
  # 
  # - time_grid_est  (m_est-dim vec)    vector of timepoints
  #
  #
  # outputs:
  #
  # - graphs   (list of graphs)       graph of all processes mapped on their Lambda_k_truth and Lambda_k_est
  #
  #
  # ----------------------------------------------------------------------------
  
  process_id <- 1:length(step_0_events)
  
  graphs <- list()
  
  nbins <- 30

  
  for(i in process_id){
    x <- step_0_events[[i]]
    
    df_est <- data.frame(x = time_grid_est, y = step_1$Lambda_k_est[i, , 1] / nbins)
    df_truth <- data.frame(x = time_grid, y = step_1$Lambda_k_truth[i, , 1] / nbins)
    
    graphs[[i]] <- ggplot() + geom_histogram(data = data.frame(x), aes(x = x, fill = 'Events'), bins = nbins) + 
      geom_line(data = df_est, aes(x = x, y = y, color = 'Estimate')) + 
      geom_line(data = df_truth, aes(x = x, y = y, color = 'Truth')) + 
      labs(
        title = paste0("Process ", i),
        x = "Time",
        y = "Count or Intensity / Bin Width",
        fill = "Event",
        color = "Intensity",
      ) +
      scale_fill_manual(values = c("Events" = "lightblue")) +
      scale_color_manual(values = c("Estimate" = "red", "Truth" = "green")) +
      theme_minimal()
  }
  
  return(graphs)
  
  
  
}

# visualize ||V_cond - V_cond_est||_HS convergence

visualize_V_cond_convergence <- function(folder_name, mat_name, i, j){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: visualize V_cond_convergence at block matrix (i, j)
  #
  #       error =  || V_{X_i, X_j}^(y_c) - \hat{V}_{X_i, X_j}^(y_c) ||_HS
  #
  #
  # input:
  #
  # - folder_name (string)  'simu_results_banded_c1_3'
  # - mat_name    (string)  'V', 'C', or 'P'
  # - i, j        (scalar)  block matrix numbers
  #
  #
  # output:
  #
  # - graph of V_12 error as we change sample size 
  #
  #
  # ----------------------------------------------------------------------------
  
  
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <-  as.numeric(sub(".*_(.*)\\.RData$", "\\1", files)) # retrieve numbers
  
  results_list <- lapply(files, function(f) {
    e <- new.env()          # create an isolated environment
    load(f, envir = e)      # load into that environment
    as.list(e)              # convert to a list (in case multiple objects)
  })
  
  names(results_list) <- ns
  
  
  # retrieve the 1_2 entry from metrics of V_cond 9written by chatgpt)
  
  results_df <- do.call(rbind, lapply(names(results_list), function(top_name) {
    
    top_entry <- results_list[[top_name]]$graph_results_i$estimated_graphs_part_2
    
    do.call(rbind, lapply(names(top_entry), function(low_name) {
      low_entry <- top_entry[[low_name]]
      
      # Dynamically select the matrix
      mat <- if (mat_name == "V") {
        low_entry$metrics$V_HS
      } else if (mat_name == "C") {
        low_entry$metrics$C_HS
      } else if (mat_name == "P") {
        low_entry$metrics$P_HS
      } else {
        stop("Unknown mat_name: must be 'V', 'C', or 'P'")
      }
      
      value <- mat[i, j]
      
      data.frame(
        top_level = top_name,
        low_level = low_name,
        value = value,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  colnames(results_df) <- c('n', 'y_c_query', 'V_12_error')
  
  results_df$n <- as.numeric(results_df$n)
  results_df$y_c_query <- as.numeric(results_df$y_c_query)
  
  if(mat_name == 'V'){
    title_name <- 'Covariance'
  } else if(mat_name == 'C'){
    title_name <- 'Correlation'
  } else if (mat_name == 'P'){
    title_name <- 'Precision'
  } else{
    stop("Unknown mat_name: must be 'V', 'C', or 'P'")
  }
  
  g <- ggplot() + geom_line(data = results_df, aes(x = y_c_query, y = V_12_error, group = n, color = n)) + 
    ylab(paste0(mat_name, '_' , i, '_', j, ' Error')) + 
    xlab('Y_c Query') + 
    ggtitle(paste0('Conditional ', title_name, ' Operator Convergence')) + 
    ylim(0, NA) + 
    theme_bw() 
  
  return(g)
  
  
}

# across all n's, plot Y_c_query vs AUC

visualize_AUC_across_n <- function(folder_name){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: visualize AUC metric across n and y_c_query
  #
  #
  #
  # input:
  #
  # - folder_name (string)  'simu_results_banded_c1_3'
  #
  #
  # output:
  #
  # - graph of AUC vs y_c_query (x-axis) and n (color)
  #
  #
  # ----------------------------------------------------------------------------
  
  
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <-  as.numeric(sub(".*_(.*)\\.RData$", "\\1", files)) # retrieve numbers
  
  results_list <- lapply(files, function(f) {
    e <- new.env()          # create an isolated environment
    load(f, envir = e)      # load into that environment
    as.list(e)              # convert to a list (in case multiple objects)
  })
  
  names(results_list) <- ns
  
  
  # retrieve the 1_2 entry from metrics of V_cond 9written by chatgpt)
  
  results_df <- do.call(rbind, lapply(names(results_list), function(top_name) {
    
    top_entry <- results_list[[top_name]]$graph_results_i$estimated_graphs_part_2
    
    do.call(rbind, lapply(names(top_entry), function(low_name) {
      low_entry <- top_entry[[low_name]]
      
      
      
      AUC <- low_entry$metrics$auc
      
      data.frame(
        top_level = top_name,
        low_level = low_name,
        value = AUC,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  colnames(results_df) <- c('n', 'y_c_query', 'AUC')
  
  results_df$n <- as.numeric(results_df$n)
  results_df$y_c_query <- as.numeric(results_df$y_c_query)
  
  title_name <- 'AUC'
  
  
  
  g <- ggplot() + geom_line(data = results_df, aes(x = y_c_query, y = AUC, group = n, color = n)) + 
    ylab('AUC') + 
    xlab('Y_c Query') + 
    ggtitle('AUC versus n and y_c_query') + 
    ylim(0, 1) + 
    theme_bw() 
  
  return(g)
  
  
}


visualize_metrics <- function(folder_name, metrics, i = NULL, j = NULL){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: visualize convergences of various metrics:
  #
  # - metrics
  #   - rho_i_dist (scalar)
  #   - rho_ij_dist (pxp matrix, each value is HS norm of m_est x m_est rho_ij)
  #   - g_ij_dist   (pxp matrix)
  #   - P_HS        (pxp matrix, each value is HS norm of difference of P_hat - P)
  #   - C_HS        (pxp matrix)
  #   - C_HS_v2     (pxp matrix)
  #   - V_HS        (pxp matrix)
  #   - sens        (scalar)
  #   - spec        (scalar)
  #   - auc         (scalar)
  #   - accuracy    (scalar)
  #
  #
  #
  # input:
  #
  # - folder_name (string)  'simu_results_banded_c1_3'
  # - metric      (string)  
  #   - 'rho_i_dist'
  #   - 'rho_ij_dist'
  #   - 'g_ij_dist'
  #   - 'P_HS', 'C_HS', 'V_HS'
  # - i and j    (integers)  indices for matrix metrics
  #
  # output:
  #
  # - table and graph of intermediate convergence metrics:
  # 
  #   - ||rho_i(t) - rho_i_est(t)||
  #   - ||rho_ij(s,t) - rho_ij_est(s,t)||
  #   - ||g_ij(s,t) - g_ij_est(s,t)||
  #   - ||C_ij - C_ij_est||
  #   - ||P_ij - P_ij_est||
  #   - AUC
  #
  # ----------------------------------------------------------------------------
  
  
  # 1) loading all datasets in a folder
  
  results_list <- load_all_results(folder_name)
  
  # 2) getting all metrics
  
  # top_name = '100', '200' etc n
  # low_name = '0', '0.125', '0.25', etc y_c_query  
  results_df <- do.call(rbind, lapply(names(results_list), function(top_name) {
    
    top_entry <- results_list[[top_name]]$graph_results_i$metrics
    
    low_names <- results_list[[top_name]]$graph_results_i$y_c_query %>% as.character()
    
    names(top_entry) <- low_names
    
    do.call(rbind, lapply(names(top_entry), function(low_name) {
      low_entry <- top_entry[[low_name]]
      
      values <- c()
      for(metric in metrics){
        
        metric_value <- low_entry[[metric]]
        
        if(metric %in% c('rho_ij_dist', 'g_ij_dist', 'P_HS', 'C_HS', 'C_HS_v2', 'V_HS')){
          metric_value = metric_value[i, j]
        }
        
        values <- c(values, metric_value)
        
      }
      
      
      c(top_name, low_name, values)
    }))
  }))
  
  # names 
  metric_names <- c()
  for(metric in metrics){
    
    
    if(metric %in% c('rho_ij_dist', 'g_ij_dist', 'P_HS', 'C_HS', 'V_HS')){
      metric_name <- paste0(metric, '_', i, '_', j)
    } else{
      metric_name <- metric  
    }
    
    metric_names <- c(metric_names, metric_name)
  }
  
  results_df <- data.frame(results_df) %>% mutate_all(as.numeric)
  colnames(results_df) <- c('n', 'y_c_query', metric_names)
  
  
  results_df$n <- as.factor(results_df$n)
  
  # graph
  graphs <- list()
  
  for(metric_name in metric_names){
    
    title_name <- paste0(metric_name, ' versus n and y_c_query')
    
    g <- ggplot() + geom_line(data = results_df, aes(x = y_c_query, y = .data[[metric_name]], group = n, color = n)) + 
      geom_point(data = results_df, aes(x = y_c_query, y = .data[[metric_name]], group = n, color = n)) + 
      ylab(metric_name) + 
      xlab('Y_c Query') + 
      # ggtitle(title_name) + 
      theme_bw() 
    
    if(metric_name %in% c('auc')){
      g <- g + ylim(0, 1)
    } else{
      g <- g + scale_y_log10(limits = c(NA, NA))
    }
    graphs[[metric_name]] <- g
  }
  
  # grid arrange
  
  arranged_plots <- do.call(arrangeGrob, c(graphs, ncol = 3))
  
  # Display it
  return(list(metric_graph = arranged_plots,
              metric_table = results_df))
  
  
}


visualize_truths_from_est <- function(step_list, graph_ids, time_grid_est, time_grid, time_grid_both, p){
  
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: visualize intermediate metrics for simulations
  #
  # inputs:
  #
  # - step_list       (list)
  # - graph_ids       (vector)             numerical string ID's of graphs we want
  # - time_grid_est   (m_est-dim vector)
  # - time_grid       (m-dim vector)
  # - time_grid_both  (m_both-dim vector)
  # - p               (integer)
  #
  # 
  # outputs:
  #
  # - graphs (list)   list of grpahs for each correspondign graph_id
  # 
  #
  # ----------------------------------------------------------------------------
  
  
  # load step_list
  
  y_c_id <- 1 # first y_c_query
  
  step_0 <- step_list$step_0 
  step_1 <- step_list$step_1 
  step_2 <- step_list$step_2[[y_c_id]] 
  step_3 <- step_list$step_3[[y_c_id]] 
  step_4 <- step_list$step_4[[y_c_id]]  
  step_5 <- step_list$step_5[[y_c_id]]  
  step_6 <- step_list$step_6[[y_c_id]]  
  step_7 <- step_list$step_7[[y_c_id]]  
  step_8 <- step_list$step_8[[y_c_id]]  
  step_9 <- step_list$step_9[[y_c_id]]  
  step_10 <- step_list$step_10[[y_c_id]]  
  step_11 <- step_list$step_11[[y_c_id]]  
  step_12 <- step_list$step_12[[y_c_id]]  
  step_2b <- step_list$step_2b[[y_c_id]]  
  

  

  
  graphs = list()
  
  arr_mat_3 <- matrix(1:3, nrow = 1, byrow = F)
  arr_mat_4 <- matrix(1:4, nrow = 2, byrow = F)
  arr_mat_6 <- matrix(1:6, nrow = 2, byrow = F)
  arr_mat_8 <- matrix(1:8, nrow = 2, byrow = F)
  arr_mat_10 <- matrix(1:10, nrow = 2, byrow = F)  
  
  m <- length(time_grid)
  m_est <- length(time_grid_est)
  
  
  if('01' %in% graph_ids){ graphs[['g_01']] <- result_01(step_0) }
  
  if('11' %in% graph_ids){ 

    g_list <- result_11_prep(step_1, time_grid, time_grid_est)
    grob_caption <- "0. Log Intensity\n of first 5 processes\nof subject 1"
    graphs[['g_11']] <- result_arr_mat(g_list, grob_caption, arr_mat_6)

    }
  
  
  # rho_i for first 5 processes 
  if('22' %in% graph_ids){ 

    
    g_list <- result_22_prep(step_2, time_grid, time_grid_est)
    grob_caption <- "2. Rho_i\nEstimation for \n first 5 processes"
    graphs[['g_22']] <- result_arr_mat(g_list, grob_caption, arr_mat_6)
    

    }
  
  
  # rho_ij(s,t) for process pair 1_1
  if('24' %in% graph_ids){ 
    
    g_list <- result_20s_prep(step_2b, i = 1, j = 1)
    grob_caption <- "2. Rho ii\nEstimation \n for processes 1_1"
    graphs[['g_24']] <- result_arr_mat(g_list, grob_caption, arr_mat_6)    
    
  }
  
  
  # rho_ij(s,t) for process pair 1_2
  if('25' %in% graph_ids){ 
    
    g_list <- result_20s_prep(step_2b, i = 1, j = 2)
    grob_caption <- "2. Rho ii\nEstimation \n for processes 1_2"
    graphs[['g_25']] <- result_arr_mat(g_list, grob_caption, arr_mat_6)  
    
  }
  
  # rho_ij (s,t) for process pair 3_5
  if('26' %in% graph_ids){ 
    g_list <- result_20s_prep(step_2b, i = 3, j = 5)
    grob_caption <- "2. Rho ii\nEstimation \n for processes 3_5"
    graphs[['g_26']] <- result_arr_mat(g_list, grob_caption, arr_mat_6)  
  
    }
  
  # distribution of weights due to continuous covariate

  if('29' %in% graph_ids){ 
    graphs[['g_29']] <- result_29(step_2, y_c_id) 
  }
  
  # g_ij(s,t) at 1_1
  if('31' %in% graph_ids){ 
    
    g_list <- result_30s_prep(step_3, i = 1, j = 1)
    grob_caption <- "3. Covariance Function\nEstimation (G_ii)\nfor 1_1"
    graphs[['g_31']] <- result_arr_mat(g_list, grob_caption, arr_mat_8) 
    
  }
  
  # g_ij(s,t) at 1_2
  if('32' %in% graph_ids){ 
    g_list <- result_30s_prep(step_3, i = 1, j = 2)
    grob_caption <- "3. Covariance Function\nEstimation (G_ii)\nfor 1_2"
    graphs[['g_32']] <- result_arr_mat(g_list, grob_caption, arr_mat_8) 
  }
  

  # g_ij(s,t) at 3_5
  if('33' %in% graph_ids){ 
    g_list <- result_30s_prep(step_3, i = 3, j = 5)
    grob_caption <- "3. Covariance Function\nEstimation (G_ii)\nfor 3_5"
    graphs[['g_33']] <- result_arr_mat(g_list, grob_caption, arr_mat_8) 
  }    
  
  # 4.1) plot eigenfunctions
  
  if('41' %in% graph_ids){ 
    
    g_list <- result_41_prep(step_4, time_grid, time_grid_est)
    grob_caption <- "4. Eigenfunction\nBasis of\nProcess 1"
    graphs[['g_41']] <- result_arr_mat(g_list, grob_caption, arr_mat_6)
    
  }

  # visualize g_ij vs its eigendecomposition reconstruction
  
  if('42' %in% graph_ids){ 
    
    g_list <- result_42_prep(step_3, step_4, p)
    grob_caption <- "4. Eigenfunction\nReconstruction"
    graphs[['g_42']] <- result_arr_mat(g_list, grob_caption, arr_mat_8)
    
  }
  
  #  are all eigenfunctions orthogonal? Orthonormal? 
  if('43' %in% graph_ids){ 
    
    g_list <- result_43_prep(step_4)
    grob_caption <- "4. Eigenfunction\nOrthogonality"
    graphs[['g_43']] <- result_arr_mat(g_list, grob_caption, arr_mat_6)
    
  }
  
  # reconstruction error histogram
  if('44' %in% graph_ids){ 
    
    g_list <- result_44_prep(step_3, step_4, p)
    grob_caption <- "4. Eigendecomposition\nReconstruction Error"
    graphs[['g_44']] <- result_arr_mat(g_list, grob_caption, arr_mat_6)
    
  }
  
  # eigenvalue decay - check if G_ii / m makes the eigenvalues similar between coarse/fine settings
  
  if('45' %in% graph_ids){
    
    # keep eigenvalues of at most 5 processes
    step_4_evals <- lapply(step_4, function(x) x$eigenvalues[1:min(5, length(x$eigenvalues))])
    
    
    # reorganize dataframe of eigenvalues 
    df <- data.frame()
    
    for (outer_name in names(step_4_evals)) {
      inner_list <- step_4_evals[[outer_name]]
      
      for (j in seq_along(inner_list)) {
        vals <- inner_list[[j]]
        temp <- data.frame(
          description = outer_name,
          process = j,
          eigen_index = seq_along(vals),
          value = vals
        )
        df <- rbind(df, temp)
      }
    }
    
    # Plot: one plot per component, colored by process
    graphs[['g_45']] <- ggplot(df, aes(x = eigen_index, y = value, color = description)) +
      geom_line() + 
      geom_point() + 
      facet_wrap(~process, scales = "free_y") +
      labs(x = "Eigenvalue index", y = "Eigenvalue", color = "Process") +
      theme_minimal()    
    
  }
  
  if('46' %in% graph_ids){
    
    # keep eigenvalues of at most 5 processes
    step_4_evecs <- lapply(step_4, function(x) x$eigenfunctions)    
    
    Y1 <- data.frame(step_4_evecs[[1]][[1]])
    Y2 <- data.frame(step_4_evecs[[2]][[1]])
    
    
    df1 <- as.data.frame(Y1) %>%
      mutate(x = time_grid) %>%
      pivot_longer(
        cols = -x,
        names_to = "series",
        values_to = "y"
      ) %>%
      mutate(group = "50-dim")
    
    df2 <- as.data.frame(Y2) %>%
      mutate(x = time_grid_est) %>%
      pivot_longer(
        cols = -x,
        names_to = "series",
        values_to = "y"
      ) %>%
      mutate(group = "19-dim")
    
    # Combine both
    df_all <- bind_rows(df1, df2)
    
    # Plot
    graphs[['g_45']] <- ggplot(df_all, aes(x = x, y = y, group = interaction(series, group), color = group)) +
      geom_line(size = 0.8, alpha = 0.9) +
      scale_color_manual(values = c("50-dim" = "red", "19-dim" = "blue")) +
      theme_minimal() +
      labs(
        x = "Time",
        y = "Value",
        color = "Series type",
        title = "Eigenvectors of Fine and Coarse Grids"
      )
    
    
  }
  
  
  # 5.3) Graph X_k (true, true_reconstruct, estimate, estimate_reconstruct)  for the first subject
  
  if (any(c('51', '52', '53') %in% graph_ids)) {
    #   - X_k_est                   (p x m_est   x n)
    #   - X_k_truth                 (p x m_truth x n)
    #   - X_k_coarse_truth          (p x m_est   x n)
    #   - X_k_both_truth            (p x m_both  x n)
    #   - Lambda_k_truth            (p x m_truth x n)
    #   - Lambda_k_coarse_truth     (p x m_est   x n)    
    


    df1 <- validate_kl_full(step_1[[1]], step_4[[5]]$eigenfunctions, time_grid_est,  'Estimate',            ncores)
    df2 <- validate_kl_full(step_1[[3]], step_4[[4]]$eigenfunctions, time_grid_est,  'Coarse Truth X',      ncores)
    df3 <- validate_kl_full(step_1[[2]], step_4[[3]]$eigenfunctions, time_grid,      'Truth X',             ncores)
    df4 <- validate_kl_full(step_1[[3]], step_4[[2]]$eigenfunctions, time_grid_est,  'Coarse Truth Theory', ncores)
    df5 <- validate_kl_full(step_1[[2]], step_4[[1]]$eigenfunctions, time_grid,      'Truth Theory',        ncores)

    df_both <- rbind(df1, df5)
    df_both$cat <- factor(df_both$cat)

    df_all <- rbind(df1, df2, df3, df4, df5)
    df_all$cat <- factor(df_all$cat)

    df_no_reconstruct <- df_all %>% filter(cat %in% c('Estimate', 'Coarse Truth X', 'Truth X', 'Coarse Truth Theory', 'Truth Theory'))


    # estimate and truth theory + reconstruct
    graphs[['g_51']] <- ggplot() +
      geom_line(data = df_both, aes(x = Var2, y = value, color = cat)) +
      facet_wrap(~ Var1) +
      ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
      ylab('Log Intensity') +
      xlab('Time') +
      labs(color = "Truth or Estimate") +
      theme_bw()

    # estimate, truth x and truth theory + reconstruct
    graphs[['g_52']] <- ggplot() +
      geom_line(data = df_all, aes(x = Var2, y = value, color = cat)) +
      facet_wrap(~ Var1) +
      ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
      ylab('Log Intensity') +
      xlab('Time') +
      labs(color = "Truth or Estimate") +
      theme_bw()

    # estimate, truth x and truth theory + no reconstruct
    graphs[['g_53']] <- ggplot() +
      geom_line(data = df_no_reconstruct, aes(x = Var2, y = value, color = cat)) +
      facet_wrap(~ Var1) +
      ggtitle('KL Reconstruction of Log Intensities for Subject 1')  +
      ylab('Log Intensity') +
      xlab('Time') +
      labs(color = "Truth or Estimate") +
      theme_bw()
  }
  
  # 54 = how do the eigenvalues in step 5, var(alpha) = mu, differ from those in step 4?
  
  if('54' %in% graph_ids){
    

    p_min <- min(p, 5)
    
    # extract diagonal entries from i_i
    step_5_diags <- lapply(step_5, function(inner_list) {

      ii_names <- paste0(1:p_min, "_", 1:p_min)
      
      # keep only those entries that exist in the inner list
      ii_mats <- inner_list[ii_names]
      
      # extract diagonals from each
      lapply(ii_mats, diag)
    })
    
    # eigenvalues from step 4
    step_4_evals <- lapply(step_4, function(x) x$eigenvalues[1:p_min])
    
    # 3) for each process, est/truth combo, plot the eigenvalues to make sure they are equal
    
    df <- data.frame()
    
    names_truth <- factor(c('Truth', 'Coarse_Truth', 'X_Truth', 'X_Coarse_Truth', 'Est'), 
                          levels = c('Truth', 'Coarse_Truth', 'X_Truth', 'X_Coarse_Truth', 'Est'))
    
    for (k in seq_along(step_5_diags)) {  # over a,b,c,d,e
      name_k <- names_truth[k]
      
      for (p in seq_along(step_5_diags[[k]])) {  # over 1:p eigenvalue vectors
        v1 <- step_5_diags[[k]][[p]]
        v2 <- step_4_evals[[k]][[p]]
        
        n <- min(length(v1), length(v2))  # ensure same length
        
        temp <- data.frame(
          index = seq_len(n),
          eigen1 = v1[1:n],
          eigen2 = v2[1:n],
          process = name_k,
          component = p
        )
        
        df <- rbind(df, temp)
      }
    }    
    
    graphs[['g_54']] <- ggplot(df, aes(x = index)) +
      geom_line(aes(y = eigen1, color = "Var(KL)")) +
      geom_line(aes(y = eigen2, color = "Eigenvalue"), linetype = "dashed") +
      geom_point(aes(y = eigen1, color = "Var(KL)")) +
      geom_point(aes(y = eigen2, color = "Eigenvalue")) +      
      facet_grid(process ~ component, scales = "free_y") +
      labs(
        title = "Eigenvalue Comparison between Two Lists",
        x = "Eigenvalue index",
        y = "Eigenvalue",
        color = "Source"
      ) +
      theme_minimal()
    
  }
  
  
  # for pair 1_2, plot estimate and truths - this one should look correlated
  
  if('81' %in% graph_ids){
    
    g_list <- result_80s_prep(step_8, m, m_est, i = 1, j = 2)
    grob_caption <- "8. Conditional\nCovariance Operator\nProcess 1 with 2"
    graphs[['g_81']] <- result_arr_mat(g_list, grob_caption, arr_mat_8)
    
  }
  
  # visualize again for process 3 with 6 - they should not look correlated  
  if('82' %in% graph_ids){
    
    g_list <- result_80s_prep(step_8, m, m_est, i = 3, j = 6)
    grob_caption <- "8. Conditional\nCovariance Operator\nProcess 3 with 6"
    graphs[['g_82']] <- result_arr_mat(g_list, grob_caption, arr_mat_8)    
  }
  

  # est vs truth vs ground truth
  if('83' %in% graph_ids){
    graphs[['g_83']] <- grid.arrange(visualize_pm_block_matrix_heatmap(step_8[[1]], 'Ground Truth'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[3]], 'Truth Theory'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[7]], 'Estimate'),
                                     nrow = 1)    
  }

  
  # all 5 + ground truth
  if('84' %in% graph_ids){
    graphs[['g_84']] <- grid.arrange(visualize_pm_block_matrix_heatmap(step_8[[1]], 'Ground Truth'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[2]], 'Coarse Ground Truth'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[3]], 'Truth Theory'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[4]], 'Coarse Truth Theory'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[5]], 'Truth X'), 
                                     visualize_pm_block_matrix_heatmap(step_8[[6]], 'Coarse Truth X'), 
                                     textGrob("8. Conditional\nCovariance Operator\nAll Processes", gp = gpar(fontsize = 14)),
                                     visualize_pm_block_matrix_heatmap(step_8[[7]], 'Estimate'),
                                     layout_matrix = arr_mat_8)      
  }
  

  

  

  if('91' %in% graph_ids){
    i <- 1
    j <- 2
    g_list <- result_90s_prep_ij(step_9, m, m_est, i, j)
    caption <- "9. Conditional\nCorrelation Operator\nProcess 1 with 2"
    graphs[['g_91']] <- result_arr_mat(g_list, caption, arr_mat_10)    
  }

  if('92' %in% graph_ids){
    i <- 3
    j <- 5
    g_list <- result_90s_prep_ij(step_9, m, m_est, i, j)
    caption <- "9. Conditional\nCorrelation Operator\nProcess 3 with 6"
    graphs[['g_92']] <- result_arr_mat(g_list, caption, arr_mat_10)
  }
  
  if('93' %in% graph_ids){

    g_list <- result_90s_prep_pm(step_9, m, m_est)
    grob_caption <- "9. Conditional\nCorrelation Operator\nAll Processes"
    g_list <- g_list[c(1,3,7)]
    graphs[['g_93']] <- result_arr_mat(g_list, grob_caption, arr_mat_4)  
    
  }
  
  
  
  # all 5 + ground truth
  
  if('94' %in% graph_ids){
    
    g_list <- result_90s_prep_pm(step_9, m, m_est)
    grob_caption <- "9. Conditional\nCorrelation Operator\nAll Processes"
    graphs[['g_94']] <- result_arr_mat(g_list, grob_caption, arr_mat_8)     
    
  }

  
  # sanity check - hilbert schmidt norms of the C_Cond
  if('95' %in% graph_ids){
    
    g_list <- result_95_prep(step_9, m, m_est, p)
    grob_caption <- "9. Conditional\nCorrelation Operator\nHS Norm"
    graphs[['g_95']] <- result_arr_mat(g_list, grob_caption, arr_mat_8)    

  }
 
  
  if('101' %in% graph_ids){
    
    i <- 1
    j <- 2
    g_list <- result_100s_prep_ij(step_10, m, m_est, i, j)
    caption <- "10. Conditional\nPrecision Operator\n Process 1 and 2"
    graphs[['g_101']] <- result_arr_mat(g_list, caption, arr_mat_8)       
    
  }

  # P_Xi_Xj for processes 3 and 5
  
  if('102' %in% graph_ids){
    i <- 3
    j <- 5
    g_list <- result_100s_prep_ij(step_10, m, m_est, i, j)
    caption <- "10. Conditional\nPrecision Operator\n Process 3 and 5"
    graphs[['g_102']] <- result_arr_mat(g_list, caption, arr_mat_8)       
  }
  

  # 10.2) visualize the entire pm x pm block
  
  if('103' %in% graph_ids){
    
    g_list <- result_100s_prep_pm(step_10, m, m_est)
    grob_caption <- "10. Conditional\nPrecision Operator\nAll Processes"
    g_list <- g_list[c(1,3,7)]
    graphs[['g_103']] <- result_arr_mat(g_list, grob_caption, arr_mat_4)      
    
  }
      
  
  
  # all 5 + ground truth
  if('104' %in% graph_ids){
    
    
    g_list <- result_100s_prep_pm(step_10, m, m_est)
    grob_caption <- "10. Conditional\nPrecision Operator\n All Processes"
    graphs[['g_104']] <- result_arr_mat(g_list, grob_caption, arr_mat_8)
    
  }
 

  # visualize HS norms with adj_mat (0 or 1)
  # g_111 <- grid.arrange(visualize_pm_block_matrix_heatmap(adj_mat_i, 'Adjacency Truth'), 
  #                       visualize_pm_block_matrix_heatmap(adj_mat_i, 'Coarse Adjacency Truth'), 
  #                       visualize_pm_block_matrix_heatmap(w_mat_truth, 'Truth Theory'), 
  #                       visualize_pm_block_matrix_heatmap(w_mat_coarse_truth, 'Coarse Truth Theory'), 
  #                       visualize_pm_block_matrix_heatmap(w_mat_X_truth, 'Truth X'), 
  #                       visualize_pm_block_matrix_heatmap(w_mat_X_coarse_truth, 'Coarse Truth X'), 
  #                       textGrob("11. Final Estimates\nvs Adj Truth", gp = gpar(fontsize = 14)),
  #                       visualize_pm_block_matrix_heatmap(w_mat_est, 'Estimate'),
  #                       layout_matrix = arr_mat_8)  
  
  
  
  # g_111b <- grid.arrange(visualize_pm_block_matrix_heatmap(adj_mat_i, 'Coarse Ground Truth'),
  #                        visualize_pm_block_matrix_heatmap(adj_mat_i, 'Ground Truth'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_coarse_truth, 'Coarse Truth Theory'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_truth, 'Truth Theory'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_X_coarse_truth, 'Coarse Truth X'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_X_truth, 'Truth X'), 
  #                        visualize_pm_block_matrix_heatmap(w_mat_normalized_est, 'Estimate'),
  #                        textGrob("11b. Normalized\nFinal Estimates\nvs Adj Truth", gp = gpar(fontsize = 14)),
  #                        layout_matrix = arr_mat_8)      
  
  
  
  
  # visualize HS norms with w_mat ground truth
  
  if('112' %in% graph_ids){

    g_list <- result_112_prep(step_11, remove_diag = F)
    grob_caption <- "11. Hilbert Schmidt\n Norm"
    graphs[['g_112']] <- result_arr_mat(g_list, grob_caption, arr_mat_8)
    
  }
  
  # remove the diagonals
  if('112b' %in% graph_ids){
    
    g_list <- result_112_prep(step_11, remove_diag = T)
    grob_caption <- "11. Hilbert Schmidt\n Norm"
    graphs[['g_112']] <- result_arr_mat(g_list, grob_caption, arr_mat_8)
    
  }
  
  if('113' %in% graph_ids){
    
    g_list <- result_113_prep(step_12)
    grob_caption <- "12. ROC Curve"
    graphs[['g_113']] <- result_arr_mat(g_list, grob_caption, arr_mat_8)
    
    
  }
  
  return(graphs)
 
  
}

unpack_step_list <- function(folder_name, n, query_id){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: get step_list of a specific dataset in a folder, to prepare to visualize results
  #
  # - metrics
  #   - rho_i_dist (scalar)
  #   - rho_ij_dist (pxp matrix, each value is HS norm of m_est x m_est rho_ij)
  #   - g_ij_dist   (pxp matrix)
  #   - P_HS        (pxp matrix, each value is HS norm of difference of P_hat - P)
  #   - C_HS        (pxp matrix)
  #   - V_HS        (pxp matrix)
  #   - sens        (scalar)
  #   - spec        (scalar)
  #   - auc         (scalar)
  #   - accuracy    (scalar)
  #
  #
  #
  # input:
  #
  # - folder_name    (string)         'simu_results_banded_c1_3'
  # - n              (integer)        sample size result in the folder that we want
  # - query_id       (integer)        integer denoting which query to look at
  #
  # output:
  #
  # - step_list (list with all relevant estimates)
  #
  # ----------------------------------------------------------------------------
  
  
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <- suppressWarnings({
    as.numeric(sub(".*_(.*)\\.RData$", "\\1", files))   # retrieve numbers
  })
  
  
  non_na_idx <- which(!is.na(ns))
  
  method <- sub("_.*", "", sub(".*/", "", files[non_na_idx[1]]))
  
  ns <- ns[!is.na(ns)]
  
  if(n %in% ns){
    idx <- which(ns == n)
  } else{
    stop('n not available')
  }
  
  load(files[non_na_idx[idx]])
  
  if(method %in% c('OG', 'JASA')){
    step_list <- list(step_0 = graph_results_i$estimated_graphs_part_1$step_0,
                      step_1 = graph_results_i$estimated_graphs_part_1$step_1,
                      step_2 = graph_results_i$estimated_graphs_part_1$step_2,
                      step_3 = graph_results_i$estimated_graphs_part_1$step_3,
                      step_4 = graph_results_i$estimated_graphs_part_1$step_4,
                      step_5 = graph_results_i$estimated_graphs_part_1$step_5,
                      step_6 = NA,
                      step_7 = NA,
                      step_8 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_8,
                      step_9 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_9,
                      step_10 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_10,
                      step_11 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_11,
                      step_12 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_12,
                      step_2b = graph_results_i$estimated_graphs_part_1$step_2b
    )    
  } else if(method == 'CPGM'){
    step_list <- list(step_0 = graph_results_i$estimated_graphs_part_1$step_0,
                      step_1 = graph_results_i$estimated_graphs_part_1$step_1,
                      step_2 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_2,
                      step_3 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_3,
                      step_4 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_4,
                      step_5 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_5,
                      step_6 = NA,
                      step_7 = NA,
                      step_8 = NA,
                      step_9 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_9,
                      step_10 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_10,
                      step_11 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_11,
                      step_12 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_12,
                      step_2b = graph_results_i$estimated_graphs_part_2[[query_id]]$step_2b
    )    
  } else{
    stop('unknown method')
  }
  
  return(step_list)
  
}

unpacking_pipeline <- function(folder_name, graph_id, time_grid_est, time_grid, n, p, query_id){
  
  # ----------------------------------------------------------------------------
  #
  #
  # GOAL: for a specific query_id and specific n, visualize a specific result
  #
  # - metrics
  #   - rho_i_dist (scalar)
  #   - rho_ij_dist (pxp matrix, each value is HS norm of m_est x m_est rho_ij)
  #   - g_ij_dist   (pxp matrix)
  #   - P_HS        (pxp matrix, each value is HS norm of difference of P_hat - P)
  #   - C_HS        (pxp matrix)
  #   - V_HS        (pxp matrix)
  #   - sens        (scalar)
  #   - spec        (scalar)
  #   - auc         (scalar)
  #   - accuracy    (scalar)
  #
  #
  #
  # input:
  #
  # - folder_name    (string)         'simu_results_banded_c1_3'
  # - graph_id       (string)         '113' for ROC curve 
  # - time_grid_est  (m_est-dim vec)  time grid discretization
  # - time_grid      (m-dim vec)      time grid discretization
  # - n              (integer)
  # - p              (integer)
  # - query_id       (integer)        integer denoting which query to look at
  #
  # output:
  #
  # - visualize_truths_from_est() output
  #
  # ----------------------------------------------------------------------------
  
  # 1) load file and method 
  graph_results_i <- get_file_name_and_load(folder_name, n)  
  
  method <- get_method(folder_name)
  
  if(method %in% c('OG', 'JASA')){
    step_list <- list(step_0 = graph_results_i$estimated_graphs_part_1$step_0,
                      step_1 = graph_results_i$estimated_graphs_part_1$step_1,
                         step_2 = graph_results_i$estimated_graphs_part_1$step_2,
                         step_3 = graph_results_i$estimated_graphs_part_1$step_3,
                         step_4 = graph_results_i$estimated_graphs_part_1$step_4,
                         step_5 = graph_results_i$estimated_graphs_part_1$step_5,
                         step_6 = NA,
                         step_7 = NA,
                         step_8 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_8,
                         step_9 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_9,
                         step_10 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_10,
                         step_11 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_11,
                         step_12 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_12,
                         step_2b = graph_results_i$estimated_graphs_part_1$step_2b
    )    
  } else if(method == 'CPGM'){
    step_list <- list(step_0 = graph_results_i$estimated_graphs_part_1$step_0,
                      step_1 = graph_results_i$estimated_graphs_part_1$step_1,
                           step_2 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_2,
                           step_3 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_3,
                           step_4 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_4,
                           step_5 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_5,
                           step_6 = NA,
                           step_7 = NA,
                           step_8 = NA,
                           step_9 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_9,
                           step_10 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_10,
                           step_11 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_11,
                           step_12 = graph_results_i$estimated_graphs_part_2[[query_id]]$step_12,
                           step_2b = graph_results_i$estimated_graphs_part_2[[query_id]]$step_2b
    )    
  } else{
    stop('unknown method')
  }
  
  # perform visualization
  return(visualize_truths_from_est(step_list, graph_id, time_grid_est, time_grid, time_grid_both, p))
  
}







