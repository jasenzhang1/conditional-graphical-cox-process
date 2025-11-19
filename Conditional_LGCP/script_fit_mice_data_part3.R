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
movement <- args[5]                 # movement <- 0
VR <- args[6]                       # VR <- 0
cont_inds <- args[6]                # all of the cont_inds

discrete_level <- paste0('m', movement, 'vr', VR)


temp_file_dir <- 'temp_data'

folder_1_name <- 'mice_results'
if (!dir.exists(folder_1_name)) dir.create(folder_1_name)  # /mice_results

folder_2_name <- paste0(folder_1_name, "/", y_c_structure)
if (!dir.exists(folder_2_name)) dir.create(folder_2_name)   # /mice_results/week_only

results_folder_name <- paste0(folder_2_name, "/", method) 
if (!dir.exists(results_folder_name)) dir.create(results_folder_name)  # /mice_results/week_only/CPGM

# ---------------------------
# estimation
# ---------------------------

if(method == 'CPGM'){
  graph_results_i <- full_conditional_estimation_with_no_truth_part3(temp_file_dir, cont_inds)
} else{
  stop('Invalid method. Must be CPGM')
}

print('obtained estimate')

file_dir <- paste0(results_folder_name, '/', ID, '_', discrete_level, '_t', time_scale, '.RData')
save(graph_results_i, file = file_dir)

print('saved estimate')



t1 <- Sys.time()

print(paste0('Time to finish: ', round(as.numeric(t1 - t0, units = "mins"), 2), ' minutes'))  
print(strrep("-", 50))


# ALL WARNINGS
print('all warnings below:')
print(warnings())

