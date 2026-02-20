
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP

args <- commandArgs(trailingOnly = TRUE)

n_large <- as.numeric(args[1])
rep_i <- as.numeric(args[2])
adj_type <- args[3]
cont_ind <- as.numeric(args[4])
beta_truth <- as.logical(args[5])
adj_params <- as.numeric(args[6:length(args)])

setting_info_list <- list(n = n_large,
                          rep_i = rep_i,
                          beta_truth = beta_truth,
                          adj_type = adj_type,
                          adj_params = adj_params)

source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')


# 1) generate dataset ----------------------------------------------------------

temp_file_dir <- 'temp_data/simu_data'

simulate_finite_basis_cox_data_part4(temp_file_dir, setting_info_list, cont_ind)




