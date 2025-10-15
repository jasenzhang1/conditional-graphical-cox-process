t0 <- Sys.time()
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP


source('functions/00_function_wrapper.R')

# args 
args <- commandArgs(trailingOnly = TRUE)

n_large <- as.numeric(args[1])
n <- as.numeric(args[2])
adj_type <- args[3]
method <- args[4]

terse <- T
# ncores <- parallel::detectCores() - 1
ncores <- 1

#simu_results_a_sparse_v2_OG
folder_1_name <- 'simu_results'
if (!dir.exists(folder_1_name)) dir.create(folder_1_name)

folder_2_name <- paste0(folder_1_name, "/", adj_type)
if (!dir.exists(folder_2_name)) dir.create(folder_2_name)

results_folder_name <- paste0(folder_2_name, "/", method)
if (!dir.exists(results_folder_name)) dir.create(results_folder_name)

print(paste0("n_large: ", n_large))
print(paste0("n: ", n))
print(paste0("adj_type: ", adj_type))
print(paste0("method: ", method))

# load `dataset`
load(paste0('simu_data/', adj_type, '_n_', n_large, '.RData'))

idx <- round(seq(1, n_large, length.out = n))

dataset_i <- dataset
dataset_i$event_times <- NULL
dataset_i$subject_data <- dataset$subject_data[idx]
dataset_i$X_k_truth <- dataset$X_k_truth[,,idx]
dataset_i$X_k_coarse_truth <- dataset$X_k_coarse_truth[,,idx]
dataset_i$X_k_both_truth <- dataset$X_k_both_truth[,,idx]
dataset_i$Y_continuous <- matrix(dataset$Y_continuous[idx,], nrow = length(idx))
dataset_i$simulation_params$n <- length(idx)


# estimate
graph_results_i <- full_conditional_estimation_with_truths_v4(dataset_i, method, terse, ncores, results_folder_name)


print('obtained estimate')

file_dir <- paste0(results_folder_name, '/', method, '_', 'n_', n, '.RData')
save(graph_results_i, file = file_dir)

print('saved estimate')


t1 <- Sys.time()

print(paste0('Time to finish: ', round(as.numeric(t1 - t0, units = "mins"), 2), ' minutes'))  
print(strrep("-", 50))


# ALL WARNINGS
print('all warnings below:')
print(warnings())

