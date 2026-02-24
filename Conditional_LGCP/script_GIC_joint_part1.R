suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(1); omp_set_num_threads(1)
source('functions/00_function_wrapper.R')

args <- commandArgs(trailingOnly = TRUE)
# Reconstruct GIC_folder from args as you did in your previous scripts
model_type <- args[1] 

if(model_type == 'mice'){
  ID <- args[2]                       # ID <- 'Tau1'
  y_c_structure <- args[3]            # y_c_structure <- "week_only" or "time_and_week"
  time_scale <- as.numeric(args[4])   # time_scale <- 10   (each replicate is 5 seconds)
  method <- args[5]                   # method <- 'CPGM'
  movement <- as.numeric(args[6])     # movement <- 0
  VR <- as.numeric(args[7])           # VR <- 0
  cont_ind <- as.numeric(args[8])     # cont_ind <- 1  
  eigen_setting <- args[9]            # eigen_seting <- 'only_joint'
  
  X_truth <- F
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  
  temp_file_dir <- 'temp_data/mice'
  GIC_folder <- paste0('temp_data/mice/GIC_global_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind)
  mouse <- T
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  rep_i <- as.numeric(args[4])
  adj_type <- args[5]
  method <- args[6]
  X_truth <- as.logical(args[7])
  cont_inds <- as.numeric(args[8])
  eigen_setting <- args[9]
  

  
  temp_file_dir <- 'temp_data/simu'
  GIC_folder <- paste0('temp_data/simu/GIC_global_', adj_type, '_n_', n, '_nquery', cont_ind, '_rep', rep_i)
  mouse <- F
  
} else{
  stop('model_type not supported')
}

# 0) load 


# load everything from step_1
if(mouse){
  step_3_info_list <- paste0(temp_file_dir, '/part3_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
  step_2_list_names <- paste0(temp_file_dir, '/part2_', ID, '_', discrete_level, '_t', time_scale, '_nquery', 1:cont_inds, '.rds')
} else{
  step_3_info_list <- paste0(temp_file_dir, '/part3_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_rep_', rep_i, '.rds')
  step_2_list_names <- paste0(temp_file_dir, '/part2_', adj_type, '_n_', n, '_nquery', 1:cont_inds, '_rep_', rep_i, '.rds')
}

step_2_all_data <- lapply(step_2_list_names, readRDS)
W_y_list <- lapply(step_2_all_data, `[[`, "W_y")
p <- step_2_all_data[[1]]$p

results <- lapply(step_3_info_list, readRDS)

# KL_cor_est_eig1, etc...
est_names <- names(results[[1]]$step_5b)

# 2. Run Setup
max_k <- GIC_joint_part1_setup(C_cond_list, p, W_y_list, GIC_folder)

# 3. Output for Bash
cat(max_k)