library(RhpcBLASctl)
blas_set_num_threads(1); omp_set_num_threads(1)
source('functions/00_function_wrapper.R')


args <- commandArgs(trailingOnly = TRUE)
model_type <- args[1] 


if(model_type == 'mice'){
  ID <- args[2]                       # ID <- 'Tau1'
  y_c_structure <- args[3]            # y_c_structure <- "week_only" or "time_and_week"
  time_scale <- as.numeric(args[4])   # time_scale <- 10   (each replicate is 5 seconds)
  method <- args[5]                   # method <- 'CPGM'
  movement <- as.numeric(args[6])     # movement <- 0
  VR <- as.numeric(args[7])           # VR <- 0
  cont_inds <- as.numeric(args[8])
  
  X_truth <- F
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level)
  
  temp_file_dir <- 'temp_data/simu'
  GIC_folder <- paste0('temp_data/mice/GIC_joint_', ID, '_', discrete_level, '_t', time_scale)
  mouse <- T
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  rep_i <- as.numeric(args[4])
  adj_type <- args[5]
  method <- args[6]
  cont_inds <- as.numeric(args[7])
  
  setting_info_list <- list(n_large = n_large,
                            n = n,
                            rep_i = rep_i,
                            adj_type = adj_type,
                            method = method)
  
  temp_file_dir <- paste0('temp_data/simu_', adj_type, '_n_', n_large, '_rep_', rep_i)
  GIC_folder <- paste0(temp_file_dir, '/GIC_joint_', adj_type, '_n_', n, '_rep', rep_i)
  mouse <- F
  
} else{
  stop('model_type not supported')
}


GIC_joint_part4_finalize(GIC_folder, temp_file_dir, setting_info_list, cont_inds, mouse)
