# fitting finite basis data

t0 <- Sys.time()
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP


source('functions/00_function_wrapper.R')

# args 
args <- commandArgs(trailingOnly = TRUE)

n_large <- as.numeric(args[1])  # n_large <- 200
n <- as.numeric(args[2])        # n <- 200
adj_type <- args[3]             # adj_type <- 'block_banded_v2'
method <- args[4]               # method <- 'CPGM'


# ncores <- parallel::detectCores() - 1
ncores <- 1


folder_1_name <- 'simu_results'
if (!dir.exists(folder_1_name)) dir.create(folder_1_name)  # simu_results

folder_2_name <- paste0(folder_1_name, "/", adj_type)
if (!dir.exists(folder_2_name)) dir.create(folder_2_name)  # simu_results/single_c2

results_folder_name <- paste0(folder_2_name, "/", method)
if (!dir.exists(results_folder_name)) dir.create(results_folder_name)  # simu_results/single_c2/CPGM

print(paste0("n_large: ", n_large))
print(paste0("n: ", n))
print(paste0("adj_type: ", adj_type))
print(paste0("method: ", method))

# load `dataset`
load(paste0('simu_data/', adj_type, '_n_', n_large, '.RData'))

idx <- round(seq(1, n_large, length.out = n))

dataset_i <- dataset

# Split strings by "_"
split_list <- strsplit(names(dataset_i$event_times), "_")

# Extract the first part as numeric
i_values <- sapply(split_list, function(x) as.numeric(x[1]))
selected_event_times <- i_values %in% idx

dataset_i$event_times <- dataset$event_times[selected_event_times]
dataset_i$X_k_truth <- dataset$X_k_truth[,,idx]
dataset_i$X_k_coarse_truth <- dataset$X_k_coarse_truth[,,idx]
dataset_i$X_k_both_truth <- dataset$X_k_both_truth[,,idx]
dataset_i$Y_continuous <- matrix(dataset$Y_continuous[idx,], nrow = length(idx))
dataset_i$simulation_params$n <- length(idx)


# estimate
graph_results_i <- full_conditional_estimation_with_no_truth(dataset_i, method, ncores, results_folder_name)

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

