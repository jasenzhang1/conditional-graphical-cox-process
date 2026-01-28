
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP

args <- commandArgs(trailingOnly = TRUE)

n_large <- as.numeric(args[1])
adj_type <- args[2]
group_idx <- as.numeric(args[3])
n_group <- as.numeric(args[4])
min_events <- as.numeric(args[5])
max_events <- as.numeric(args[6])

setting_info_list <- list(n = n_large,
                          adj_type = adj_type)

source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')


# 1) generate dataset ----------------------------------------------------------

temp_file_dir <- 'temp_data/simu_data'

simulate_finite_basis_cox_data_parts1_and_2(temp_file_dir, setting_info_list, group_idx, n_group, min_events, max_events)
                                          





