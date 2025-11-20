
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP

args <- commandArgs(trailingOnly = TRUE)

n <- as.numeric(args[1])
adj_type <- args[2]
cont_ind <- as.numeric(args[3])

source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')


# 1) generate dataset ----------------------------------------------------------

temp_file_dir <- 'temp_data/simu_data'

simulate_finite_basis_cox_data_part4(temp_file_dir, setting_info_list, cont_ind)




