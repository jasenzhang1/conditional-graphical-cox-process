# 1) obtain mice data

# 2) for each mouse, fit the model


t0 <- Sys.time()

source('functions/00_function_wrapper.R')

# 1) system parameters
seed = 1
ncores = parallel::detectCores() - 1

# 2) data/output parameters
terse = TRUE

data_folder_name <- 'data/'
results_folder_name <- 'mouse_results_1'
if (!dir.exists(results_folder_name)) dir.create(results_folder_name)

# 3) continuous covariate parameters
query_y_cs <- c(0.25, 0.5, 0.75) #normalized age
q_c = 1                         # Continuous conditioning dimension (q_c)



# 4) time discretization
time_grid <- 1:19/20
m <- length(time_grid)


# 7) sample size and # of processes
ns <- c(100, 300, 1000, 3000, 10000)     # Sample size (n)
n_large <- max(ns)
p = 10                                   # Number of processes (p)  

# 8) mouse parameters
IDs <- c('346', '351', '366', '361', '362', '368')  # mouse ID
ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3') # our name 
num_neurons <- c(169, 250, 240, 249, 235, 294)

# 8.1) truncated mouse parameters
IDs <- c('351',  '368')  
ID2 <- c('Tau2', 'WT3') 
num_neurons <- c(250,  294)




# starting ---------------------------------------------------------------------

print('Starting ----------------------------')
print(paste0('number of cores: ', ncores))

for(i in 1:length(IDs)){ 
  
  # 1) load data of mouse_i
  
  
  
  load(paste0(data_folder_name, ID2[i], '_data.rda'))
  p <- num_neurons[i]
  
  # recall that this is called LGCP_data
  # [[1]] = main dataframe with 'feature_id', 'time', and 'subject_num'
  # [[2]] = dataframe of discrete covariates
  # [[3]] = dataframe of continuous covariates
  
  data_df <- LGCP_data[[1]] 
  y_d_df <- LGCP_data[[2]]
  y_c_df <- LGCP_data[[3]] %>% select(c('subject_num', 'age'))
  
  
  
  # 1.1) obtain discrete strata from y_d_df
  
  cov_df <- y_d_df %>% dplyr::select(- subject_num)
  discrete_strata <- cov_df %>% unique()
  strata_stats <- data.frame()
  
    
  for(y_ind in 1:nrow(discrete_strata)){ # 2) for each discrete variable level 
    
    t0 <- Sys.time()
    
    # 4) filtering by strata + included neurons to get data_df4
    
    yd_i <- discrete_strata[y_ind, ]
    discrete_strata_name <- paste0('m', yd_i$movement, 'vr', yd_i$VR)
    
    s_yd <- y_d_df$subject_num[apply(cov_df, 1, function(x) all(x == yd_i))]  # s_yd = indices
    data_df4 <- data_df %>% filter(subject_num %in% s_yd)
    
    
    # 5) counting neurons
    
    y_c_strata <- y_c_df %>% filter(subject_num %in% s_yd) %>% dplyr::select(-subject_num) %>% as.matrix()
    
    silent_neurons <- setdiff(1:neuron_count,
                              unique(data_df4$feature_id))
    
    
    # size of dataset
    print(paste0('number of subjects: ', length(s_yd)))
    print(paste0('number of spikes: ', nrow(data_df4)))
    
    # print number of neurons
    print(paste('total neurons: ', unname(neuron_count)))                     # total neurons
    print(paste('active neurons: ', length(unique(data_df4$feature_id))))     # active neurons
    print(paste('silent neurons: ', length(silent_neurons)))                  # silent neurons
    print(paste('do they add up? ', length(unique(data_df4$feature_id)) + length(silent_neurons) == unname(neuron_count)))
    print('-------------------------------------')
    
    # 6) retroactively get parameters for subjects (ntrain, patient_sel) and neurons (p, feature_sel)
    
    ntrain <- length(unique(data_df4$subject_num))
    p <- length(unique(data_df4$feature_id))
    
    
    patient_sel = unique(data_df4$subject_num) %>% sort()
    feature_sel = unique(data_df4$feature_id) %>% sort()
    
    # 7) estimation 
    
  }

}



all_results <- list()

for(n in ns){
  
  t_n_start <- Sys.time()
  
  dataset_i <- dataset
  dataset_i$subject_data <- dataset$subject_data[1:n]
  dataset_i$X_k_truth <- dataset$X_k_truth[,,1:n]
  dataset_i$X_k_coarse_truth <- dataset$X_k_coarse_truth[,,1:n]
  dataset_i$X_k_both_truth <- dataset$X_k_both_truth[,,1:n]
  dataset_i$Y_continuous <- matrix(dataset$Y_continuous[1:n,], nrow = n)
  dataset_i$simulation_params$n <- n
  
  graph_results_i <- full_conditional_estimation_with_truths(dataset_i, terse, ncores)
  
  file_dir <- paste0(results_folder_name, '/n_', n, '.RData')
  save(graph_results_i, file = file_dir)
  
  t_n_end <- Sys.time()
  
  rm(graph_results_i)
  rm(dataset_i)
  
  unlink("~/.RData")
  unlink("~/.Rhistory")
  unlink("~/.local/share/rstudio/sessions", recursive = TRUE)  
  
  print('Done --------------------------------')
  print(paste0('All Loops: ', paste(ns, collapse = ' ')))
  print(paste0('Current Loop: ', n))
  print(paste0('Time to finish: ', round(as.numeric(t_n_end - t_n_start, units = "mins"), 2), ' minutes'))  
  print(strrep("-", 50))
}


print(paste0('Grand total time: ', round(as.numeric(t_n_end - t0, units = "hours"), 2), ' hours'))  