
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP


source('functions/00_function_wrapper.R')

# args 
args <- commandArgs(trailingOnly = TRUE)

model_type <- args[1]                 # model_type = 'simu', or 'mice'

if(model_type == 'mice'){
  ID <- args[2]                       # ID <- 'Tau1'
  y_c_structure <- args[3]            # y_c_structure <- "week_only" or "time_and_week"
  time_scale <- as.numeric(args[4])   # time_scale <- 10   (each replicate is 5 seconds)
  method <- args[5]                   # method <- 'CPGM'
  movement <- as.numeric(args[6])     # movement <- 0
  VR <- as.numeric(args[7])           # VR <- 0
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level)
  
  data_folder <- 'mice_data'
  if (!dir.exists(data_folder)) dir.create(data_folder)  # /mice_data
  
  temp_file_dir <- 'temp_data'
  if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # /temp_data
  
  temp_file_dir <- paste0(temp_file_dir, '/mice')
  if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data/mice
  mouse <- T
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  adj_type <- args[4]
  method <- args[5]
  
  setting_info_list <- list(n_large = n_large,
                            n = n,
                            adj_type = adj_type,
                            method = method)
  
  temp_file_dir <- 'temp_data'
  if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data
  
  temp_file_dir <- paste0(temp_file_dir, '/simu')
  if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data/simu
  mouse <- F
  
} else{
  stop('model_type not supported')
}

ncores <- 1


# ---------------------------
# estimation
# ---------------------------

if(mouse){
  load(paste0('mice_data/', ID, '_', discrete_level, '_t', time_scale, '.RData')) # dataset_k
} else{
  load(paste0('simu_data/', adj_type, '_n_', n_large, '.RData')) #dataset --> dataset_k
  
  idx <- round(seq(1, n_large, length.out = n))
  
  dataset_k <- dataset
  
  # Split strings by "_"
  split_list <- strsplit(names(dataset_k$event_times), "_")
  
  # Extract the first part as numeric
  i_values <- sapply(split_list, function(x) as.numeric(x[1]))
  selected_event_times <- i_values %in% idx
  
  dataset_k$event_times <- dataset$event_times[selected_event_times]
  dataset_k$X_k_truth <- dataset$X_k_truth[,,idx]
  dataset_k$X_k_coarse_truth <- dataset$X_k_coarse_truth[,,idx]
  dataset_k$X_k_both_truth <- dataset$X_k_both_truth[,,idx]
  #dataset_k$Y_continuous <- matrix(dataset$Y_continuous[idx,], nrow = length(idx))  # filtering y's as well
  dataset_k$simulation_params$n <- length(idx)  
}

  

if(method == 'CPGM'){
  full_conditional_estimation_with_no_truth_part1(dataset_k, setting_info_list, ncores, temp_file_dir, mouse)
} else{
  stop('Invalid method. Must be CPGM')
}



