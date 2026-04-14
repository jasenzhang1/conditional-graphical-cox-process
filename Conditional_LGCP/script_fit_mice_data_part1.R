# ------------------------------------------------------------------------------
# 
# GOAL: prepare for full_conditional_estimation_with_no_truth_part1
#
#
# inputs:
#
# - ID               (string)    mouse name such as 'Tau1'
# - y_c_structure    (string)    "week_only" or "time_and_week"
# - time_scale       (integer)   how many seconds is each replicate? Values may be 1, 2, 5, 10
# - method           (string)    estimation method, "CPGM"
# - m                (integer)   time_grid spacing
# - movement         (integer)   0 (resting) or 1 (moving)
# - VR               (integer)   0 (off) or 1 (on)
# - max_processes    (integer)   how many processes should we truncate? 
#
# 
# outputs:
#
# 
# ------------------------------------------------------------------------------





source('functions/00_function_wrapper.R')

# args 
args <- commandArgs(trailingOnly = TRUE)

model_type <- args[1]                 # model_type = 'simu', or 'mice'

if(model_type == 'mice'){
  ID <- args[2]                       # ID <- 'Tau3'
  y_c_structure <- args[3]            # y_c_structure <- "week_only" or "time_and_week"
  time_scale <- as.numeric(args[4])   # time_scale <- 10   (each replicate is 5 seconds)
  method <- args[5]                   # method <- 'CPGM'
  movement <- as.numeric(args[6])     # movement <- 0
  VR <- as.numeric(args[7])           # VR <- 0
  cluster <- args[8]
  
  if(cluster == 'andrew'){
    library(RhpcBLASctl)
    
    # limit threads in BLAS/LAPACK
    blas_set_num_threads(1)   # limit BLAS
    omp_set_num_threads(1)    # limit OpenMP
  }
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level)
  
  data_folder <- 'mice_data'
  if (!dir.exists(data_folder)) dir.create(data_folder)  # /mice_data
  
  data_folder <- paste0(data_folder, '/', y_c_structure) 
  if (!dir.exists(data_folder)) dir.create(data_folder)  # /mice_data/week_only
  
  if(cluster == 'andrew'){
    temp_file_dir <- 'temp_data'
  } else{
    temp_file_dir <- '../../../temp_data'
  }
  
  if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # /temp_data
  
  temp_file_dir <- paste0(temp_file_dir, '/mice')
  if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data/mice
  
  temp_file_dir <- paste0(temp_file_dir, '/', ID, '_', discrete_level, '_t', time_scale)
  if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data/mice/Tau3_m0vr0_t10/
  
  mouse <- T
  X_truth <- F
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  rep_i <- as.numeric(args[4])
  adj_type <- args[5]
  method <- args[6]
  X_truth <- as.logical(args[7])
  
  setting_info_list <- list(n_large = n_large,
                            n = n,
                            rep_i = rep_i,
                            adj_type = adj_type,
                            method = method,
                            X_truth = X_truth)
  
  temp_file_dir <- 'temp_data'
  if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data
  
  temp_file_dir <- paste0('temp_data/simu_', adj_type, '_n_', n_large, '_rep_', rep_i)
  if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # temp_data/simu_adj_type
  mouse <- F
  
} else{
  stop('model_type not supported')
}

ncores <- 1


# ---------------------------
# estimation - read in data from simu_data/ or mice_data/
# ---------------------------

if(mouse){
  load(paste0(data_folder, '/', ID, '_', discrete_level, '_t', time_scale, '.RData')) # dataset_k
  dataset_k$Y_continuous_k <- dataset_k$Y_continuous  
} else{
  simu_data_dir <- paste0('simu_data/', adj_type, '_n_', n_large, '_rep_', rep_i)
  load(paste0(simu_data_dir, '/', adj_type, '_n_', n_large, '_rep_', rep_i, '.RData')) #dataset --> dataset_k
  
  # 1) choose indices that thin the sample size
  idx <- round(seq(1, n_large, length.out = n))
  
  dataset_k <- dataset
  
  # 2) Split strings and extract k and i from event_times
  split_list <- strsplit(names(dataset_k$event_times), "_")

  k_values <- sapply(split_list, function(x) as.numeric(x[1]))
  i_values <- sapply(split_list, function(x) as.numeric(x[2]))
  
  # 3) keep only entries where i is in idx
  keep <- k_values %in% idx
  dataset_k$event_times <- dataset_k$event_times[keep]
  k_values <- k_values[keep]
  i_values <- i_values[keep]
  
  # 4) Compress k to 1:n (preserving order of appearance)
  k_map <- match(k_values, unique(k_values))
  
  # 5) Rename entries as "newk_i"
  names(dataset_k$event_times) <- paste0(k_map, "_", i_values)
  
  
  # 6) update the rest of the dataset
  dataset_k$X_k_truth <- dataset$X_k_truth[,,idx]
  dataset_k$X_k_coarse_truth <- dataset$X_k_coarse_truth[,,idx]
  dataset_k$X_k_both_truth <- dataset$X_k_both_truth[,,idx]
  dataset_k$Y_continuous_k <- dataset_k$Y_continuous  
  dataset_k$Y_continuous_k <- matrix(dataset_k$Y_continuous_k[idx,], nrow = length(idx))  # creating filtered Y_c and unfiltered Y_c
  dataset_k$beta_coeffs <- dataset$beta_coeffs[,,idx]
  dataset_k$simulation_params$n <- length(idx)  
}

  

if(method == 'CPGM'){
  full_conditional_estimation_with_no_truth_part1(dataset_k, setting_info_list, ncores, temp_file_dir, mouse, X_truth)
} else{
  stop('Invalid method. Must be CPGM')
}



