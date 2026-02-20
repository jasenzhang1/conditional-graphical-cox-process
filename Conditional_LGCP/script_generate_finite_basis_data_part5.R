
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP

args <- commandArgs(trailingOnly = TRUE)

n_large <- as.numeric(args[1])
rep_i <- as.numeric(args[2])
adj_type <- args[3]
cont_inds <- as.numeric(args[4])
groups <- as.numeric(args[5])

setting_info_list <- list(n = n_large,
                          rep_i = rep_i,
                          adj_type = adj_type)

source('functions/00_function_wrapper.R')
source('functions/20_simulation_function_wrapper.R')


# 1) generate dataset ----------------------------------------------------------

temp_file_dir <- 'temp_data/simu_data'

result <- simulate_finite_basis_cox_data_part5(temp_file_dir, setting_info_list, cont_inds, groups)
                                          

dataset <- result$dataset                                           
truths <- result$all_truths
rm(result)

if (!dir.exists('simu_data')) dir.create('simu_data')

save(dataset, file = paste0('simu_data/', adj_type, '_n_', n_large, '_rep_', rep_i, '.RData'))
save(truths, file = paste0('simu_data/', adj_type, '_n_', n_large, '_rep_', rep_i, '_truths.RData'))



