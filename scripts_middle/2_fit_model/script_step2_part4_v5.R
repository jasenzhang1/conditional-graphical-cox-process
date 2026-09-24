library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP


source('functions/00_function_wrapper.R')

# args 
args <- commandArgs(trailingOnly = TRUE)

model_type <- args[1]                 # model_type = 'simu', or 'mice'

if(model_type == 'mice'){
  ID <- args[2]                              # ID <- 'Tau1'
  y_c_structure <- args[3]                   # y_c_structure <- "week_only" or "time_and_week"
  time_scale <- as.numeric(args[4])          # time_scale <- 10   (each replicate is 5 seconds)
  method <- args[5]                          # method <- 'CPGM'
  movement <- as.numeric(args[6])            # movement <- 0
  VR <- as.numeric(args[7])                  # VR <- 0
  cont_ind <- as.numeric(args[8])    
  
  use_default <- grepl('_default$', args[9])
  if (use_default) {
    y_c_bandwidth <- as.numeric(sub('_default$', '', args[9]))
  } else {
    y_c_bandwidth <- as.numeric(args[9])
  }
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level,
                            y_c_bandwidth = y_c_bandwidth)
  
  temp_file_dirs <- c(paste0('temp_data/mice/', ID, '_', discrete_level, '_t', time_scale), 'temp_data/mice_data')
  mouse <- T
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])      # which sub sample size  
  rep_i <- as.numeric(args[4])
  adj_type <- args[5]
  method <- args[6]
  cont_ind <- as.numeric(args[7])       # which y_c query
  
  use_default <- grepl('_default$', args[8])
  if (use_default) {
    y_c_bandwidth <- as.numeric(sub('_default$', '', args[8]))
  } else if (is.na(args[8]) || args[8] == '') {
    y_c_bandwidth <- NULL   # use the default gamma_c
  } else {
    y_c_bandwidth <- as.numeric(args[8])
  }

  
  setting_info_list <- list(n_large = n_large,
                            n = n,
                            rep_i = rep_i,
                            adj_type = adj_type,
                            method = method,
                            y_c_bandwidth = y_c_bandwidth,
                            use_default = use_default)
  
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
  estimate_intensities_stratum_parallel_with_yc_part4_v5(temp_file_dirs, setting_info_list, cont_ind, mouse)
} else{
  stop('Invalid method. Must be CPGM')
}







