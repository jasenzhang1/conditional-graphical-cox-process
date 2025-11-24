
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
  cont_inds <- as.numeric(args[8])    # cont_inds <- 1  
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level)
  
  
  mouse <- T
  temp_file_dir <- 'temp_data/mice'
  
  folder_1_name <- 'mice_results'
  if (!dir.exists(folder_1_name)) dir.create(folder_1_name)  # /mice_results
  
  folder_2_name <- paste0(folder_1_name, "/", y_c_structure)
  if (!dir.exists(folder_2_name)) dir.create(folder_2_name)   # /mice_results/week_only
  
  results_folder_name <- paste0(folder_2_name, "/", method) 
  if (!dir.exists(results_folder_name)) dir.create(results_folder_name)  # /mice_results/week_only/CPGM
  
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  adj_type <- args[4]
  method <- args[5]
  cont_inds <- as.numeric(args[6])
  
  setting_info_list <- list(n_large = n_large,
                            n = n,
                            adj_type = adj_type,
                            method = method)
  
  mouse <- F
  temp_file_dir <- 'temp_data/simu'
  
  
  folder_1_name <- 'simu_results'
  if (!dir.exists(folder_1_name)) dir.create(folder_1_name)  # simu_results
  
  folder_2_name <- paste0(folder_1_name, "/", adj_type)
  if (!dir.exists(folder_2_name)) dir.create(folder_2_name)  # simu_results/single_c2
  
  results_folder_name <- paste0(folder_2_name, "/", method)
  if (!dir.exists(results_folder_name)) dir.create(results_folder_name)  # simu_results/single_c2/CPGM
  
} else{
  stop('model_type not supported')
}

# ---------------------------
# estimation
# ---------------------------

if(method == 'CPGM'){
  graph_results_i <- full_conditional_estimation_with_no_truth_part3(temp_file_dir, setting_info_list, cont_inds, mouse)
} else{
  stop('Invalid method. Must be CPGM')
}

if(mouse){
  file_dir <- paste0(results_folder_name, '/', ID, '_', discrete_level, '_t', time_scale, '.RData')
} else{
  file_dir <- paste0(results_folder_name, '/', adj_type, '_n_', n, '.RData')
}


save(graph_results_i, file = file_dir)






