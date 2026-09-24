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
  min_connect_pcts <- as.numeric(args[9:length(args)])
  
  X_truth <- F
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level)
  
  temp_file_dir <- paste0('temp_data/mice/', ID, '_', discrete_level, '_t', time_scale) # mice/WT2_m0vr1_t10
  temp_file_dir <- paste0(temp_file_dir, '/GIC_local_', ID, '_', discrete_level, '_t', time_scale, '_nquery', cont_ind) # mice/WT2_m0vr1_t10/GIC_local_WT2_m0vr1_t10_nquery1
  
  
  min_connect_pcts_string <- sapply(min_connect_pcts, function(min_connect_pct) {
    # min_pct string
    if (min_connect_pct == 0) {
      min_connect_string <- 'min_00'
    } else if (min_connect_pct == 1) {
      min_connect_string <- 'min_100'
    } else {
      decimal_digits <- sub(".*\\.", "", formatC(min_connect_pct, digits = 2, format = "f"))
      min_connect_string <- paste0('min_', decimal_digits)
    }
  })
  

  
  temp_file_min_dirs <- paste0(temp_file_dir, '/', min_connect_pcts_string) # mice/WT2_m0vr1_t10/GIC_local_WT2_m0vr1_t10_nquery1/min_xx
  
  mouse <- T
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  rep_i <- as.numeric(args[4])
  adj_type <- args[5]
  method <- args[6]
  cont_ind <- as.numeric(args[7])
  min_connect_pcts <- as.numeric(args[8:length(args)])
  
  setting_info_list <- list(n_large = n_large,
                            n = n,
                            rep_i = rep_i,
                            adj_type = adj_type,
                            method = method)
  
  
  temp_file_dir <- paste0('temp_data/simu_', adj_type, '_n_', n_large, '_rep_', rep_i)
  temp_file_dir <- paste0(temp_file_dir, '/GIC_local_', adj_type, '_n_', n, '_nquery', cont_ind, '_rep', rep_i)
  temp_file_min_dirs <- paste0(temp_file_dir, '/', min_connect_pct_strings(min_connect_pcts))
  mouse <- F
  
} else{
  stop('model_type not supported')
}


GIC_step4_finalize(temp_file_dir, temp_file_min_dirs)
