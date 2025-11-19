t0 <- Sys.time()
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP


source('functions/00_function_wrapper.R')

# args 
args <- commandArgs(trailingOnly = TRUE)

ID <- args[1]                       # ID <- 'Tau1'
y_c_structure <- args[2]            # y_c_structure <- "week_only" or "time_and_week"
time_scale <- as.numeric(args[3])   # time_scale <- 10   (each replicate is 5 seconds)
method <- args[4]                   # method <- 'CPGM'
movement <- as.numeric(args[5])     # movement <- 0
VR <- as.numeric(args[6])           # VR <- 0



discrete_level <- paste0('m', movement, 'vr', VR)

setting_info_list <- list(ID = ID,
                          y_c_structure = y_c_structure,
                          time_scale = time_scale,
                          method = method,
                          discrete_level = discrete_level)

ncores <- 1


data_folder <- 'mice_data'
if (!dir.exists(data_folder)) dir.create(data_folder)  # /mice_data

temp_file_dir <- 'temp_data'
if (!dir.exists(temp_file_dir)) dir.create(temp_file_dir)  # /temp_data

print(paste0("mouse: ", ID))
print(paste0("y_c_structure: ", y_c_structure))
print(paste0("time_scale: ", time_scale))
print(paste0("estimation method: ", method))
print(paste0("discrete level: ", discrete_level))

# ---------------------------
# estimation
# ---------------------------

load(paste0(data_folder, '/', ID, '_', discrete_level, '_t', time_scale, '.RData')) # dataset_k
  

if(method == 'CPGM'){
  full_conditional_estimation_with_no_truth_part1(dataset_k, setting_info_list, method, ncores, temp_file_dir, mouse = T)
  
  n_queries <- nrow(dataset_k$simulation_params$query_y_cs)
  cat("n_queries=", n_queries, "\n")
  
} else{
  stop('Invalid method. Must be CPGM')
}


t1 <- Sys.time()

print(paste0('Time to finish: ', round(as.numeric(t1 - t0, units = "mins"), 2), ' minutes'))  
print(strrep("-", 50))


# ALL WARNINGS
print('all warnings below:')
print(warnings())

