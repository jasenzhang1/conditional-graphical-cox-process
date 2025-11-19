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
m <- args[5]                        # m <- 20
movement <- args[6]
VR <- args[7]

ncores <- 1


data_folder <- 'mice_data'
if (!dir.exists(data_folder)) dir.create(data_folder)  # /mice_data

folder_1_name <- 'mice_results'
if (!dir.exists(folder_1_name)) dir.create(folder_1_name)  # /mice_results

folder_2_name <- paste0(folder_1_name, "/", y_c_structure)
if (!dir.exists(folder_2_name)) dir.create(folder_2_name)   # /mice_results/week_only

results_folder_name <- paste0(folder_2_name, "/", method) 
if (!dir.exists(results_folder_name)) dir.create(results_folder_name)  # /mice_results/week_only/CPGM

print(paste0("mouse: ", ID))
print(paste0("y_c_structure: ", y_c_structure))
print(paste0("time_scale: ", time_scale))
print(paste0("estimation method: ", method))

# load `dataset`
load(paste0('data/with_ts/', ID, '_t', time_scale, '_data.rda'))

# ------------------------------------------------------------
# pre-processing to make it look just like "dataset"
# ------------------------------------------------------------


time_grid_est <- make_time_grid(m)
discrete_covariates <- c('movement', 'VR')
discrete_covariate_options <- expand.grid('movement' = c(0,1), 'VR' = c(0,1))


  
#m0vr0
discrete_levels = discrete_covariate_options[k,]
discrete_name <- paste0('m', movement, 'vr', VR)

dataset_k <- convert_data_for_storage(LGCP_data, y_c_structure, discrete_levels, time_grid_est, min_events = 5, max_processes = 20, seed = NULL) # 00e

# store

file_dir <- paste0('mice_data', '/', ID, '_', discrete_name, '_t', time_scale, '.RData')
save(dataset_k, file = file_dir)



print('saved dataset')


t1 <- Sys.time()

print(paste0('Time to finish: ', round(as.numeric(t1 - t0, units = "mins"), 2), ' minutes'))  
print(strrep("-", 50))


# ALL WARNINGS
print('all warnings below:')
print(warnings())

