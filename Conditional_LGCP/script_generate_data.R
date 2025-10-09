t0 <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)



# adj_type = "sparse_v2"
# adj_params <- c(0, 1, 2, 0.3, -1, 2, 0.01)

n_large <- as.numeric(args[1])
adj_type <- args[2]
adj_params <- as.numeric(args[3:length(args)])


cat("n: ", n_large, "\n")
cat("adj_type: ", adj_type, "\n")
cat("adj_params: ", paste(adj_params, collapse = ", "), "\n")

source('functions/20_simulation_function_wrapper.R')
source('functions/00_function_wrapper.R')

# 1) system parameters
seed = 1
ncores = 6

# 2) output parameters
terse = TRUE



# 3) time discretization
T_max = 1                                           # Time horizon (T)
time_grid_size = 50                                 # Time discretization (m)

time_grid <- seq(0, T_max, length.out = time_grid_size)
time_grid_est <- 1:19/20
time_grid_both <- sort(union(time_grid, time_grid_est))

# 4) base covariance
# kernels = 'exponential', 'rbf', 'rbf_pd', 'polynomial'
base_kernel_params <- list(base_gamma = 20,      
                           base_kernel = 'rbf',   
                           base_variance = 1,
                           base_GP_mean = 5)

if(base_kernel_params$base_kernel == 'rbf_pd'){
  base_kernel_params$base_gamma = 2 * time_grid_size * log(2 * time_grid_size)
}


# 6) sample size and # of processes

query_y_cs <- matrix(0:9/9)               # query y_values

p = 10                                   # Number of processes (p)  


# 1) generate dataset ----------------------------------------------------------

dataset <- simulate_conditional_cox_data_v4(n_large, p, T_max, query_y_cs,
                                            adj_type,
                                            adj_params,
                                            time_grid,
                                            time_grid_est,
                                            base_kernel_params,
                                            ncores,
                                            seed = NULL,
                                            verbose = FALSE)

if (!dir.exists('simu_data')) dir.create('simu_data')

save(dataset, file = paste0('simu_data/', adj_type, '_n_', n_large, '.RData'))

# time taken
t1 <- Sys.time()
elapsed_time <- as.numeric(difftime(t1, t0, units = "mins"))




cat('Time to generate ', n, ' subjects: ', elapsed_time, ' mins\n')


