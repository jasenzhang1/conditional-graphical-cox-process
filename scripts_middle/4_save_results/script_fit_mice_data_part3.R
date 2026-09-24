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
  
  
  use_default <- grepl('_default$', args[9])
  if (use_default) {
    y_c_bandwidth <- as.numeric(sub('_default$', '', args[9]))
  } else {
    y_c_bandwidth <- as.numeric(args[9])
  }
  
  region            <- args[10]                  # region <- 'HIP'
  experiment_folder <- args[11]                  # experiment_folder <- 'experiment_11'
  min_connect_pcts <- as.numeric(args[12:length(args)])
  
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
  
  discrete_level <- paste0('m', movement, 'vr', VR)
  
  setting_info_list <- list(ID = ID,
                            y_c_structure = y_c_structure,
                            time_scale = time_scale,
                            method = method,
                            discrete_level = discrete_level)
  
  
  mouse <- T
  
  library(RhpcBLASctl)
  
  # limit threads in BLAS/LAPACK
  blas_set_num_threads(1)   # limit BLAS
  omp_set_num_threads(1)    # limit OpenMP
  
  temp_file_dir <- paste0('temp_data/mice/', ID, '_', discrete_level, '_t', time_scale) # temp_data/mice/WT1_m1vr1_t10/
  temp_file_dir2 <- NA
  
  folder_1_name <- 'mice_results'
  if (!dir.exists(folder_1_name)) dir.create(folder_1_name)  # /mice_results
  
  # bw string
  # "1.0_default" or "001" given the bandwidth
  if (use_default) {
    bw_string <- paste0(y_c_bandwidth, '_default')
  } else {
    bw_string <- sub(".*\\.", "", format(y_c_bandwidth, scientific = FALSE))
  }
  
  
  
} else if(model_type == 'simu'){
  
  n_large <- as.numeric(args[2])
  n <- as.numeric(args[3])
  rep_i <- as.numeric(args[4])
  adj_type <- args[5]
  method <- args[6]
  cont_inds <- as.numeric(args[7])
  min_connect_pcts <- as.numeric(args[8:length(args)])
  min_connect_pcts_string <- min_connect_pct_strings(min_connect_pcts)
  
  setting_info_list <- list(n_large = n_large,
                            n = n,
                            rep_i = rep_i,
                            adj_type = adj_type,
                            method = method,
                            min_connect_pcts_string = min_connect_pcts_string)
  
  mouse <- F
  temp_file_dir  <- paste0('temp_data/simu_', adj_type, '_n_', n_large, '_rep_', rep_i)
  temp_file_dir2 <- paste0('temp_data/simu_data_', adj_type, '_n_', n_large, '_rep_', rep_i)
  
  folder_1_name <- 'simu_results'
  if (!dir.exists(folder_1_name)) dir.create(folder_1_name)  # simu_results
  
  folder_2_name <- paste0(folder_1_name, "/", adj_type)
  if (!dir.exists(folder_2_name)) dir.create(folder_2_name)  # simu_results/single_c2
  
  results_folder_name <- paste0(folder_2_name, "/", method)
  if (!dir.exists(results_folder_name)) dir.create(results_folder_name)  # simu_results/single_c2/CPGM
  
  results_folder_name <- paste0(results_folder_name, '/rep', rep_i)
  if (!dir.exists(results_folder_name)) dir.create(results_folder_name)  # simu_results/single_c2/CPGM/rep1
  
} else{
  stop('model_type not supported')
}

# ---------------------------
# estimation - read everything from temp_data/simu (part_1, part_2, part_3)
#
#              save in simu_results/adj_type/CPGM/adj_type_n.RData
# ---------------------------

for(min_pct_string in min_connect_pcts_string){
  
  setting_info_list[[min_pct_string]] <- min_pct_string
  
  if(method == 'CPGM'){
    graph_results_i <- full_conditional_estimation_with_no_truth_part3(temp_file_dir, temp_file_dir2, setting_info_list, cont_inds, mouse)
  } else{
    stop('Invalid method. Must be CPGM')
  }
  
  if(mouse){
    results_folder_name <- paste0(folder_1_name, "/", experiment_folder, '/', y_c_structure, '_bw_', bw_string, '_', min_pct_string, "/", method, '_', region)
    if (!dir.exists(results_folder_name)) dir.create(results_folder_name, recursive = TRUE) # /mice_results/experiment_11/week_only_bw_0001_min_01/CPGM_HIP
    
    file_dir <- paste0(results_folder_name, '/', ID, '_', discrete_level, '_t', time_scale, '.RData')
  } else{
    file_dir <- paste0(results_folder_name, '/', adj_type, '_n_', n, '_rep_', rep_i, '.RData')
  }
  
  save(graph_results_i, file = file_dir)
}


