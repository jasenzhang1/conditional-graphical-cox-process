library(dplyr)

# load data that has been pre-made

load_realigned_dataset <- function(data_dir){
  
  # data_dir = /u/home/j/jasenzz/Alzheimers/Data/
  
  
  file_name <- paste(data_dir, 'Time_Realigned_Datasets.RData', sep = '')
  
  env <- new.env()
  
  load(file_name, envir = env) 
  
  return(as.list(env))
}

load_brain_region_dataset <- function(data_dir){
  
  # data_dir = /u/home/j/jasenzz/Alzheimers/Data/
  
  file_name <- paste(data_dir, 'Brain_Region.RData', sep = '')
  
  env <- new.env()
  
  load(file_name, envir = env) 
  
  return(as.list(env))
}

get_spiketrain_dataset_conditional_LGCP <- function(ID, time_scale, discrete_covariates, continuous_covariates, weekly_dataset, neuron_df){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: obtain spiketrain data from a specific mouse across multiple weeks 
  # 
  # - we partition the dataset into replicates of duration = time_scale
  # - and we keep replicates where the discrete settings don't change within each replicate
  # - so if VR switches during a replicate, we discard it
  #
  # our discrete settings:
  # 
  # 1) running/resting
  # 2) VR on/off
  #
  # 7/7/2025:
  #
  # - required to pass in discrete and continuous covariate names
  #
  # 10/13/2025: 
  # 
  # - updating 
  #
  # ------------------------------------------------------------------------
  # 
  # inputs:
  # - ID                (string):             3-digit mouse ID number in string form
  # - time_scale        (integer):            Length of time (secs) for each replicate
  # - covariates        (vector of strings):  Which covariates are we stratifying by? Options include 'movement', 'VR'
  # - weekly_dataset    (enormous dataset)
  # - neuron_df         (dataset of neuron info)
  #
  # output:
  # - list of (3) items:
  # 
  # - df_total     (data.frame):  Nx3 dataframe with the following column names
  #                                 subject_num (replicate num)
  #                                 feature_id  (neuron num)
  #                                 time        (timestamp)
  # 
  # - y_d_df (data.frame): dataframe with   discrete covariate values (movement and VR) and their subject_num
  # - y_c_df (data.frame): dataframe with continuous covariate values (possibly age and timestamp) and their subject_num
  #
  #
  # ------------------------------------------------------------------------  
  
  # 0) prepare datasets
  
  dataset <- weekly_dataset[[ID]]
  n_neurons <- neuron_info(ID, neuron_df) %>% nrow()
  elapsed_time <- 0 # how much time has passed, to pad
  first <- T
  weeks2 <- names(dataset)
  
  kept_iter <- 1
  
  for(j in weeks2){ # 1) for each week
    
    
    
    border_times_j <- c()
    
    # 1a) load data
    raw_j <- dataset[[j]][[1]]
    dataset_j <- dataset[[j]][[2]]
    tmax_j <- dataset[[j]][[3]]
    
    # long format, dataframe of ('neuron', 'spike_times')
    spikes_j <- enframe(dataset_j, name = "neuron", value = "spike_times") %>%
      unnest(spike_times)
    
    # 1b) create binary covariate labels
    
    df_cov <- data.frame(time = raw_j$time)
    for(cov_i in discrete_covariates){
      if(cov_i == 'movement'){
        running_mode <- label_activity_mode(raw_j$velocity)
        binary_label <- rep(0, length(running_mode))
        binary_label[running_mode == 'running'] <- 1      # 1 = running, 0 = resting          
      }
      
      if(cov_i == 'VR'){
        binary_label <- raw_j$isVisible                         # 1 = visible, 0 = invisible
      }
      
      
      df_cov <- df_cov %>% mutate(!!cov_i := binary_label)
      
    }
    
    
    
    # 2) now, create replicates of length time_scale 
    
    n2 <- floor(tmax_j / time_scale)
    tmax_trunc <- n2 * time_scale
    
    
    for(iter in 1:n2){ # for each iter
      t_start <- (iter - 1) * time_scale
      t_end <- iter * time_scale
      
      iter_index <- which(raw_j$time > t_start & raw_j$time < t_end)
      
      
      # filter raw_j, df_cov, spike_times to be in [t_start, t_end]
      raw_j_iter <- raw_j[iter_index, ]
      df_cov_iter <- df_cov[iter_index, ]
      spikes_iter <- spikes_j %>% filter(spike_times > t_start & spike_times < t_end)
      spikes_iter$spike_times <- (spikes_iter$spike_times - t_start) / time_scale      # time is now between 0 and 1
      
      # are there multiple levels (discrete strata) here?
      strata_iter <- df_cov_iter %>% 
        dplyr::select(-time) %>%
        unique()
      
      if(nrow(strata_iter) == 1){  # only continue if there is only 1 discrete strata level
        y_d <- strata_iter[1,]
        
        # store data in a df
        
        spikes_iter$subject_num <- kept_iter
        y_d$subject_num <- kept_iter
        
        y_c <- data.frame(kept_iter)
        colnames(y_c) <- 'subject_num'
        
        if('week' %in% continuous_covariates){
          y_c$age <- as.numeric(j)
        }
        
        if('timestamp' %in% continuous_covariates){
          y_c$timestamp  <- t_start
        }
        
        
        # update kept_iter and update df's 
        kept_iter <- kept_iter + 1
        
        if(first){
          df_total <- spikes_iter
          y_d_df <- y_d
          y_c_df <- y_c
          first <- F
        } else{
          df_total <- rbind(df_total, spikes_iter)
          y_d_df <- rbind(y_d_df, y_d)
          y_c_df <- rbind(y_c_df, y_c)
        }
        
      }
      
    } # end of each iterate
  } # end of each week
  
  # modify and return 
  
  colnames(df_total) <- c('feature_id', 'time', 'subject_num')
  df_total <- data.table(df_total)
  
  # # normalize the continuous covariates
  # 
  # if('week' %in% continuous_covariates){
  #   y_c_df$age <- (y_c_df$age - 17) / (38 - 17)
  # }
  # 
  # if('timestamp' %in% continuous_covariates){
  #   y_c_df$timestamp  <- (y_c_df$timestamp - 0) / (1500 - 0)
  # }  
  
  return(list(df_total, y_d_df, y_c_df))
  
}


label_activity_mode <- function(velocity, threshold = 4, rest_window = 300) {
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: label a velocity vector as "running" or "resting"
  #
  #
  # inputs:
  # 
  # - velocity        (N-dim vector):       vector of velocity values
  # - threshold             (number):       what speed do we differentiate between running and resting
  # - rest_window          (integer):       after how many consecutive events of "resting" do we see before we conclude that we are actually resting.
  #
  # outputs:
  #
  # - mode            (N-dim vector):       vector of "running" or "resting" elements
  #
  # ----------------------------------------------------------------------------
  
  # each dt is 1/60 of a second, so rest window is 5 seconds?
  
  n <- length(velocity)
  mode <- rep("resting", n)
  
  # Binary indicator for running
  is_running <- velocity >= threshold
  
  # Initialize state
  state <- "resting"
  rest_counter <- 0
  last_running_index <- 0
  
  for (i in seq_len(n)) { # for each timepoint
    
    if (state == "resting") { # if the state is resting but the vector says im running, turn the state to running 
      if (is_running[i]) {  
        
        state <- "running"
        mode[i] <- 'running'
        last_running_index <- i
        rest_counter <- 0
      }
    } else if (state == "running") { # if the state is running...
      if (!is_running[i]) { # ... but the vector says im not running, start the rest counter
        
        rest_counter <- rest_counter + 1
        # Only switch to resting after rest_window consecutive low values
        if (rest_counter >= rest_window) {
          
          mode[(last_running_index):i] <- "resting"
          state <- "resting"
          rest_counter <- 0          
          
        } else{
          mode[i] <- "running"
        }
      } else { # ... and the vector says i'm also running, do this
        mode[i] <- "running"
        rest_counter <- 0
        last_running_index <- i
      }
    }
  }
  
  return(mode)
}  
