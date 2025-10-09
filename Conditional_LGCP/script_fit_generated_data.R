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
ncores <- 1

#simu_results_a_sparse_v2_OG
folder_1_name <- 'simu_results'
if (!dir.exists(folder_1_name)) dir.create(folder_1_name)
results_folder_name <- paste0(folder_1_name, "/", adj_type)
if (!dir.exists(results_folder_name)) dir.create(results_folder_name)

cat("n_large: ", n_large, "\n")
cat("n: ", n, "\n")
cat("adj_type: ", adj_type, "\n")


# load `dataset`
load(paste0('simu_data/', adj_type, '_n_', n, '.RData'))

idx <- round(seq(1, n_large, length.out = n))

dataset_i <- dataset
dataset_i$event_times <- NULL
dataset_i$subject_data <- dataset$subject_data[idx]
dataset_i$X_k_truth <- dataset$X_k_truth[,,idx]
dataset_i$X_k_coarse_truth <- dataset$X_k_coarse_truth[,,idx]
dataset_i$X_k_both_truth <- dataset$X_k_both_truth[,,idx]
dataset_i$Y_continuous <- matrix(dataset$Y_continuous[idx,], nrow = length(idx))
dataset_i$simulation_params$n <- length(idx)


if(method == 'CPGM'){
  graph_results_i <- full_conditional_estimation_with_truths_v3(dataset_i, method, terse, ncores, results_folder_name)
} else if(method %in% c('OG', 'JASA')){
  graph_results_i <- full_conditional_estimation_with_truths_v2(dataset_i, method, terse, ncores, results_folder_name)
} else{
  stop('Invalid method. Must be CPGM, OG, or JASA')
}

cat('obtained estimate')

file_dir <- paste0(results_folder_name, '/', method, '_', 'n_', n, '.RData')
save(graph_results_i, file = file_dir)

cat('saved estimate')


# keep warnings
sink(paste0(results_folder_name, "/warnings.txt"))
print(warnings())
sink()

t_n_end <- Sys.time()

print(paste0('Time to finish: ', round(as.numeric(t_n_end - t_n_start, units = "mins"), 2), ' minutes'))  
print(strrep("-", 50))
