# 1) generate simulated data
# 2) estimate ground truth
# 3) calculate accuracy metrics


t0 <- Sys.time()

source('functions/20_simulation_function_wrapper.R')
source('functions/00_function_wrapper.R')

# 1) system parameters
seed = 1
ncores = parallel::detectCores() - 1

# 2) output parameters
terse = TRUE
results_folder_name <- "simu_results_banded_v2"
if (!dir.exists(results_folder_name)) dir.create(results_folder_name)

# 3) continuous covariate parameters
q_c = 1                       # Continuous conditioning dimension (q_c)
y_c_borders = list(1:4)       # Border values
K = 4                         # Discrete combinations (K)
dependence_type = "constant"  # Conditional dependence type

# 4) time discretization
T_max = 1                     # Time horizon (T)
time_grid_size = 50           # Time discretization (m)
m <- time_grid_size

time_grid <- seq(0, T_max, length.out = time_grid_size)
time_grid_est <- 1:19/20
time_grid_both <- sort(union(time_grid, time_grid_est))

# 5) base covariance
base_kernel_params <- list(base_gamma = 20,      # won't be pd, but we can massage it 
                           base_kernel = 'rbf',
                           base_variance = 1,
                           base_GP_mean = 5)

# 6) adj matrix params
sparsity = 0.2                           # Graph sparsity (s)
theta = 1.0                              # Signal strength (theta)
adj_type = "banded_v2"                   # Graph topology
adj_params <- c(3,  0.4, 0.05)           # associated parameters 

# 7) sample size and # of processes
ns <- c(100, 300, 1000, 3000, 10000)     # Sample size (n)
n_large <- max(ns)
p = 10                                   # Number of processes (p)  




# 1) generate dataset ----------------------------------------------------------

dataset <- simulate_conditional_cox_data_v4(n_large, p, T_max, q_c, y_c_borders,
                                            sparsity,
                                            theta,
                                            dependence_type,
                                            adj_type,
                                            adj_params,
                                            time_grid,
                                            time_grid_est,
                                            base_kernel_params,
                                            seed,
                                            ncores)




t1 <- Sys.time()

print(paste0('Time to generate data: ', round(as.numeric(t1 - t0, units = "mins"), 2), ' minutes'))

print(strrep("-", 50))



all_results <- list()

for(n in ns){
  
  t_n_start <- Sys.time()
  
  dataset_i <- dataset
  dataset_i$subject_data <- dataset$subject_data[1:n]
  dataset_i$X_k_truth <- dataset$X_k_truth[,,1:n]
  dataset_i$X_k_coarse_truth <- dataset$X_k_coarse_truth[,,1:n]
  dataset_i$X_k_both_truth <- dataset$X_k_both_truth[,,1:n]
  dataset_i$Y_continuous <- matrix(dataset$Y_continuous[1:n,], nrow = n)
  dataset_i$simulation_params$n <- n
  
  graph_results_i <- full_conditional_estimation_with_truths(dataset_i, terse, ncores)
  
  file_dir <- paste0(results_folder_name, '/n_', n, '.RData')
  save(graph_results_i, file = file_dir)
  
  t_n_end <- Sys.time()
  
  rm(graph_results_i)
  rm(dataset_i)
  
  unlink("~/.RData")
  unlink("~/.Rhistory")
  unlink("~/.local/share/rstudio/sessions", recursive = TRUE)  
  
  print('Done --------------------------------')
  print(paste0('All Loops: ', paste(ns, collapse = ' ')))
  print(paste0('Current Loop: ', n))
  print(paste0('Time to finish: ', round(as.numeric(t_n_end - t_n_start, units = "mins"), 2), ' minutes'))  
  print(strrep("-", 50))
}


print(paste0('Grand total time: ', round(as.numeric(t_n_end - t0, units = "hours"), 2), ' hours'))  



