


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
  n_queries <- as.numeric(args[8])    # n_queries <- 22
  
  
  ncores <- parallel::detectCores() - 1
  
  X_truth <- F
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level)
  
  temp_file_dirs <- c(paste0('temp_data/mice/', ID, '_', discrete_level, '_t', time_scale), 'temp_data/mice_data')
  mouse <- T
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  rep_i <- as.numeric(args[4])
  adj_type <- args[5]
  method <- args[6]
  X_truth <- as.logical(args[7])
  ij <- as.numeric(args[8])
  
  setting_info_list <- list(n_large = n_large,
                            n = n,
                            rep_i = rep_i,
                            adj_type = adj_type,
                            method = method)
  
  simu_dir <- paste0('temp_data/simu_', adj_type, '_n_', n_large, '_rep_', rep_i)
  simu_data_dir <- paste0('temp_data/simu_data_', adj_type, '_n_', n_large, '_rep_', rep_i)
  temp_file_dirs <- c(simu_dir, simu_data_dir)
  
  mouse <- F
  
} else{
  stop('model_type not supported')
}



# ---------------------------
# estimation
# ---------------------------

if(method == 'CPGM'){
  run_pipeline_all_queries_hoffman(temp_file_dirs, setting_info_list,
                                   n_queries, ncores, mouse, X_truth, eigen_setting)
} else{
  stop('Invalid method. Must be CPGM')
}







