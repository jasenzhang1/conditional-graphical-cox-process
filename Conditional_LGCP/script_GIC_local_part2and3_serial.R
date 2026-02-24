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
  cont_ind <- as.numeric(args[8])     # cont_ind <- 1  
  eigen_setting <- args[9]            # eigen_seting <- 'only_joint'
  id_suffix <- args[10]               # est or X_truth
  k         <- args[11]               # tau_c index
  
  
  X_truth <- F
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level)
  
  temp_file_dir <- paste0('temp_data/mice/GIC_local_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind)
  mouse <- T
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  rep_i <- as.numeric(args[4])
  adj_type <- args[5]
  method <- args[6]
  cont_ind <- as.numeric(args[7])
  eigen_setting <- args[8]
  id_suffix <- as.numeric(args[9])   # est or X_truth
  k         <- as.numeric(args[10])   # tau_c index

  
  temp_file_dir <- paste0('temp_data/simu/GIC_local_', adj_type, '_n_', n, '_nquery', cont_ind, '_rep', rep_i)
  mouse <- F
  
} else{
  stop('model_type not supported')
}

# retrieve the suffix_name

task_csv <- read.csv(paste0(temp_file_dir, '/task_map.csv'))

suffix_name <- task_csv[1, id_suffix]
GIC_step2and3_serial_tau_c(temp_file_dir, suffix_name, k)


