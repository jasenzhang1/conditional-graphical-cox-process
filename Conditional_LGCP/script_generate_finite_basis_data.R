t0 <- Sys.time()
library(RhpcBLASctl)

# limit threads in BLAS/LAPACK
blas_set_num_threads(1)   # limit BLAS
omp_set_num_threads(1)    # limit OpenMP

args <- commandArgs(trailingOnly = TRUE)

n <- as.numeric(args[1])
adj_type <- args[2]
adj_params <- as.numeric(args[3:length(args)])

# n <- 200
# adj_type <- 'block_banded_v2'
# adj_params <- c(0, 1, 0, 0.5, 3)
d <- 2
p <- 12
m <- 50
T_max <- 1


beta_0 <- 5

time_grid <- seq(0, T_max, length.out = m)
time_grid_est <- 1:19/20
time_grid_both <- sort(union(time_grid, time_grid_est))

seed <- NULL

y_c_query <- matrix(0:2/2, nrow = 3)



cat("n: ", n, "\n")
cat("adj_type: ", adj_type, "\n")
cat("adj_params: ", paste(adj_params, collapse = ", "), "\n")

source('functions/20_simulation_function_wrapper.R')
source('functions/00_function_wrapper.R')

# 1) system parameters
seed = NULL
# ncores <- parallel::detectCores() - 1
ncores = 1





# 1) generate dataset ----------------------------------------------------------

result <- simulate_finite_basis_cox_data(n, d, p, adj_type, adj_params, beta_0, 
                                         time_grid, time_grid_est, time_grid_both,
                                         T_max, y_c_query, seed)
                                          
dataset <- result$dataset                                           
truths <- result$all_truths
rm(result)

if (!dir.exists('simu_data')) dir.create('simu_data')

save(dataset, file = paste0('simu_data/', adj_type, '_n_', n, '.RData'))
save(truths, file = paste0('simu_data/', adj_type, '_n_', n, '_truths.RData'))

# time taken
t1 <- Sys.time()
elapsed_time <- as.numeric(difftime(t1, t0, units = "mins"))




cat('Time to generate ', n, ' subjects: ', elapsed_time, ' mins\n')


